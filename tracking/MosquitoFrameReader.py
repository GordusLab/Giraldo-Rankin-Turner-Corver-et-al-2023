import gc, numpy as np, imageio, glob, os, regex as re, datetime, \
    traceback, numba as nb, time, psutil, imageio.v3 as iio
from tqdm import tqdm
from motmot.ufmf.ufmf import UfmfSaverV3
from dask.distributed import Client

def lowPriority():
    # Before starting processing, set the priority to the lowest class to avoid freezing the system
    p = psutil.Process(os.getpid())
    p.nice(psutil.IDLE_PRIORITY_CLASS)

class ImageIOReader():
    def __init__(self, fname):
        try:
            self.file = iio.imopen(fname, "r", plugin="pyav")
        except:
            self.file = None
            self.gen = None
            raise Exception()
        try:
            self.gen = self.file.iter(format='rgb24', thread_type='FRAME', thread_count=16)
        except:
            self.file.close()
            self.file = None
            self.gen = None
            raise Exception()

    def get_next_data(self):
        try:
            return self.gen.__next__()
        except:
            raise StopIteration()

    def close(self):
        try:
            self.file.close()
        except:
            pass

class MosquitoFrameReader():
    def __init__(self, directory, maxFrames=-1, subtractbg=True):
        self.directory = directory
        self.maxFrames = maxFrames
        self.subtractbg = subtractbg

        # Determine file chunks in directory and sort by timestamp
        fnames = glob.glob(os.path.join(directory, '*.h264'))
        fileTimes = []
        for fname in fnames:
            try:
                month, day, year, hour, minute, second = [int(x) for x in re.search(
                    '([0-9]*)-([0-9]*)-([0-9]*)_([0-9]*)-([0-9]*)-([0-9]*).h264', fname).groups()]
                fileTimes.append(datetime.datetime(
                    year=year, month=month, day=day, hour=hour, minute=minute, second=second))
            except:
                fileTimes.append(None)
        fnames = [fn for fn, ft in zip(fnames, fileTimes) if ft is not None]
        fileTimes = [x for x in fileTimes if x is not None]
        self.fnames = [fnames[i] for i in np.argsort(fileTimes)]

        # Variables for the rolling background computation
        self.rollingBgIdx = -1000000
        self.rollingBg = None
        self.rollingBgFrameCacheSize = 300 # Keep 3,000 frame duration (300 actual frames) in memory
        self.rollingBgFrameCache = []
        self.rollingBgFrameCacheSubsample = 10

        # Remember we're currently reading the first file
        self.frameIdx = 0
        self.reader = None
        self.fileIdx = -1
        self.openChunk(0)

        self.frameIdxInChunk = 0

    def openChunk(self, idx):
        self.fileIdx = idx
        self.frameIdxInChunk = 0
        if self.reader is not None:
            self.reader.close()
            del self.reader
            gc.collect()
        if self.fileIdx < len(self.fnames):
            for tries in range(10):
                try:
                    # Option 1
                    #self.reader = imageio.get_reader(self.fnames[self.fileIdx], format='mp4')
                    # Option 2
                    self.reader = ImageIOReader(self.fnames[self.fileIdx])
                    break
                except:
                    gc.collect()
                    time.sleep(2)
                raise StopIteration
        else:
            raise StopIteration

    def __iter__(self):
        return self

    def __getFrameFromReader(self):
        try:
            fr = self.reader.get_next_data()[:,:,0].copy()
            self.frameIdxInChunk += 1
            return fr
        except:
            raise StopIteration()

    def getFrame(self):
        # If this function is run for the first time, populate the cache array by reading-ahead,
        # then resetting the reader
        if len(self.rollingBgFrameCache) == 0:
            self.rollingBgFrameCache = []
            if self.subtractbg:
                try:
                    for i in tqdm(range(self.rollingBgFrameCacheSubsample *
                            self.rollingBgFrameCacheSize), desc='Pre-populating frame cache'):
                        if (i % self.rollingBgFrameCacheSubsample) == 0:
                            self.rollingBgFrameCache.append(self.__getFrameFromReader())
                except:
                    # Resample to ensure correct size
                    while len(self.rollingBgFrameCache) < self.rollingBgFrameCacheSize:
                        w = int(np.random.randint(0, len(self.rollingBgFrameCache)))
                        self.rollingBgFrameCache.append(self.rollingBgFrameCache[w].copy())
                # Reset reader
                self.openChunk(0)
        # Decode frame data and cast to monochrome
        fr = self.__getFrameFromReader()
        # Update rolling cache
        if self.subtractbg:
            if self.frameIdx % self.rollingBgFrameCacheSubsample == 0:
                self.rollingBgFrameCache[(self.frameIdx//
                    self.rollingBgFrameCacheSubsample) % self.rollingBgFrameCacheSize] = fr.copy()
            # Update the rolling background if necessary
            if self.frameIdx - self.rollingBgIdx > 100:
                self.rollingBg = np.percentile(self.rollingBgFrameCache, 15, axis=0)
                self.rollingBgIdx = self.frameIdx
        # Compute the background-subtracted frame
        if self.subtractbg:
            frBgSubtracted = np.clip(fr.astype(int) - self.rollingBg, 0, 255).astype(np.uint8)
        else:
            frBgSubtracted = fr
        # Done! Increment frame index
        self.frameIdx += 1
        return frBgSubtracted, fr, (self.rollingBg.astype(np.uint8) if self.rollingBg is not None else None)

    def __next__(self):
        if self.frameIdx >= self.maxFrames and self.maxFrames > 0:
            raise StopIteration()
        try:
            return self.getFrame()
        except Exception as e:
            if self.fileIdx == len(self.fnames):
                raise StopIteration()
            else:
                self.openChunk(self.fileIdx + 1)
                return self.getFrame()

@nb.njit(nogil=True)
def label(d):
    localThreshold = 1
    labels = np.zeros(d.shape, dtype=np.uint16)
    numclusters = 0
    # First pass
    for y in range(labels.shape[0]):
        for x in range(labels.shape[1]):
            if d[y, x] > 0:
                # If a neighbor (up, down, left or right) is present with a label, copy that label number
                if y > 0 and labels[y - 1, x] > 0:
                    labels[y, x] = labels[y - 1, x]
                elif x > 0 and labels[y, x - 1] > 0:
                    labels[y, x] = labels[y, x - 1]
                elif y < labels.shape[0] - 1 and labels[y + 1, x] > 0:
                    labels[y, x] = labels[y + 1, x]
                elif x < labels.shape[1] - 1 and labels[y, x + 1] > 0:
                    labels[y, x] = labels[y, x + 1]
                else:
                    # No clusters nearby?
                    # -- Ensure this cluster is big enough before proceeding
                    if np.sum(d[(y - 1):(y + 2), (x - 1):(x + 2)]) > localThreshold:
                        numclusters += 1
                        labels[y, x] = numclusters

    # Second pass
    for y in range(labels.shape[0]):
        for x in range(labels.shape[1]):
            if labels[y, x] > 0:
                # Ensure that any neighboring pixels have the same label identity
                nv = labels[y, x]
                if y > 0 and d[y - 1, x] > 0:
                    v = labels[y - 1, x]
                    # Search-and-replace label
                    if nv != v:
                        for a in range(max(0, x-100), min(labels.shape[1]-1, x+100)):
                            for b in range(max(0, y-100), min(labels.shape[0]-1, y+100)):
                                if labels[b,a] == v or labels[b,a] == nv:
                                    labels[b,a] = min(v, nv)
                if x > 0 and d[y, x - 1] > 0:
                    v = labels[y, x - 1]
                    # Search-and-replace label
                    if nv != v:
                        for a in range(max(0, x-100), min(labels.shape[1]-1, x+100)):
                            for b in range(max(0, y-100), min(labels.shape[0]-1, y+100)):
                                if labels[b,a] == v or labels[b,a] == nv:
                                    labels[b,a] = min(v, nv)

    return labels

@nb.jit(nogil=True)
def dilation(img, k = 1):
    imgout = img.copy()
    h, w = img.shape

    for _k in range(k):
        for y in range(h):
            for x in range(w):
                imgout[y, x] = img[max(0, min(h-1, y-1)), x] | img[
                    max(0, min(h-1, y+1)), x] | img[y, max(0, min(w-1, x-1))] | img[y, max(0, min(w-1, x+1))]
        img = imgout.copy()
    return imgout

# Custom function to erode binary array ('k' specifies number of erosions to perform)
# Note: Wrote custom function to make this Numba (optimization) compatible
@nb.njit(fastmath=True, nogil=True, cache=True)
def erosion(img, k=1):
    imgout = img.copy()
    h, w = img.shape

    for _k in range(k):
        for y in range(h):
            for x in range(w):
                if imgout[y, x] > 0:
                    if np.sum(img[(y - 1):(y + 2), (x - 1):(x + 2)]) < 9:
                        imgout[y, x] = 0
        img = imgout.copy()
    return imgout

@nb.njit(nogil=True)
def getChangedROIs(diff, org):
    # Get thresholded image
    diffFr = diff.astype(np.float32) / org
    thresholded = dilation(erosion((diffFr >= 0.25)&(diff >= 25), 2), 24)
    # Get labels
    labels = label(thresholded)
    # Get coords
    idxs = {}
    bglabel = labels[0, 0]
    for x in range(labels.shape[0]):
        for y in range(labels.shape[1]):
            if labels[x, y] != bglabel:
                if labels[x, y] not in idxs:
                    lst = {}
                    lst[(x, y)] = True
                    idxs[labels[x, y]] = lst
                idxs[labels[x, y]][(x, y)] = True
    # Done
    ROIs = []
    irow = 0
    for z, xyKeys in enumerate(idxs.values()):
        if irow >= 1000:
            break
        xmin, xmax, ymin, ymax = 9999, -9999, 9999, -9999
        for x, y in xyKeys.keys():
            xmin = min(xmin, x)
            xmax = max(xmax, x)
            ymin = min(ymin, y)
            ymax = max(ymax, y)
        xcenter = 0.5 * xmin + 0.5 * xmax
        ycenter = 0.5 * ymin + 0.5 * ymax
        xwidth = xmax - xmin
        ywidth = ymax - ymin
        xcenter = int(0.5 + xcenter)
        ycenter = int(0.5 + ycenter)
        xwidth = int(0.5 + xwidth)
        ywidth = int(0.5 + ywidth)
        if xwidth >= 3 and ywidth >= 3:
            ROIs.append((ycenter, xcenter, ywidth, xwidth))
            irow += 1

    return ROIs

def run(directory, maxFrames=-1):
    lowPriority()

    reader = MosquitoFrameReader(directory, maxFrames=maxFrames)

    chunkID = 0

    writerBgSubtracted = None
    writerRaw = None

    for iframe, (frBgSubtracted, fr, bg) in tqdm(enumerate(reader)):
        # Open writers if necessary
        if writerRaw is None:
            fnameBgSubtracted = os.path.join(directory, 'bgsubtracted_{}c.ufmf'.format(chunkID))
            fnameRaw = os.path.join(directory, 'raw_{}c.ufmf'.format(chunkID))
            writerBgSubtracted = UfmfSaverV3(fnameBgSubtracted, max_width=fr.shape[0], max_height=fr.shape[1])
            writerRaw = UfmfSaverV3(fnameRaw, max_width=fr.shape[0], max_height=fr.shape[1])
        # Add keyframes?
        if iframe%1000 == 0:
            writerBgSubtracted.add_keyframe('mean', np.zeros(fr.shape, dtype=fr.dtype), iframe)
        if iframe%1000 == 0:
            writerRaw.add_keyframe('mean', bg, iframe)
            roisRaw = [(int(fr.shape[0]/2), int(fr.shape[1]/2), fr.shape[0], fr.shape[1]), ]
        else:
            # Compute and store changed ROIs
            roisRaw = getChangedROIs(frBgSubtracted, bg)
        # When viewed through certain GUIs, an error is thrown if the first image is all-black.
        # Change the intensity of the top pixel if necessary to fix this for convenience.
        if iframe == 0:
            if fr.max() == 0:
                fr[0, 0] = 1
            if frBgSubtracted.max() == 0:
                frBgSubtracted[0, 0] = 1
            roisRaw.append((1, 1, 2, 2))
        writerRaw.add_frame(fr, iframe, roisRaw)
        writerBgSubtracted.add_frame(frBgSubtracted, iframe, roisRaw)
        # Close chunk?
        if (iframe+1)%10000 == 0:
            writerRaw.close()
            writerBgSubtracted.close()
            writerRaw = None
            writerBgSubtracted = None
            chunkID += 1
    # Close I/O
    writerRaw.close()
    writerBgSubtracted.close()

def runSafe(directory):
    try:
        run(directory)
    except Exception as e:
        print(traceback.format_exc())

def getRecordingDirectories():
    fnamesAll = glob.glob('Z:\\Abel\\MOSQUITOS\\*\\*\\*\\*.h264') + glob.glob(
        'Z:\\Abel\\MOSQUITOS\\*\\*\\*\\*\\*.h264')
    directories = [x.replace('/', '\\') for x in set([os.path.dirname(x) for x in fnamesAll])]
    return directories

if __name__ == "__main__":
    # Get all directories
    directories = getRecordingDirectories()
    # List directories
    for d in directories:
        print('Found directory: {}'.format(d))
    time.sleep(1)
    # Process in parallel
    client = Client('10.99.66.244:8786')
    client.gather(client.map(run, directories))
