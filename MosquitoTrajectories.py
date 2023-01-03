#
# This code extracts landing trajectory events from pre-processed landing videos.
#

import numpy as np, time, os, joblib as jl, numba as nb, skimage.morphology, skimage.measure, \
    skimage.feature, glob, skimage.draw, scipy.ndimage
from tqdm import tqdm
from motmot.ufmf.ufmf import UfmfV3, FlyMovieEmulator

def getBlobsManualForFrame(tracesIdx, idx):
    if idx in tracesIdx:
        centers = []
        datarow = tracesIdx[idx]
        for k in range(10):
            if k * 5 + 3 <= len(datarow):
                # np.flip here corrects for the flipped X/Y meaning relative to the
                # predicted blob processing pipeline
                centers.append(np.flip(datarow[
                    (1 + k * 5):(3 + k * 5)]).tolist() + [0, 0, 200])
        centers = np.array(centers, dtype=float)
        return centers
    else:
        return np.zeros((0, 5), dtype=float)

def getBlobsManual(directory):
    # Load manual annotations
    fnameTraces = os.path.join(directory.replace('/', '\\').replace(
        '\\2020', '\\traces\\2020'), 'list_script_roi.txt')
    if os.path.exists(fnameTraces):
        print('Found manual traces for: {}'.format(directory))
        tracesIdx = None
        txtTraces = ''
        with open(fnameTraces, 'r') as f:
            txtTraces = f.read()
        arrTraces = [[float(y) for y in x.split(' ') if len(y.strip()) > 0] for x in txtTraces.split('\n')]
        tracesIdx = {}
        for x in arrTraces:
            if len(x) > 0:
                tracesIdx[int(x[0])] = x

        # Get the frame range
        numFrames = np.max(list(tracesIdx.keys())) + 1

        # Convert to blob format
        blobs = [getBlobsManualForFrame(tracesIdx, i) for i in range(numFrames)]

        # Done!
        return blobs
    else:
        print('Could not find manual traces for: {}'.format(directory))
        return None

def getBlobs(directory, maxNumChunks=-1):
    # Get blobs for each frame
    blobs = []
    for chunkID in range(1000 if maxNumChunks < 0 else maxNumChunks):
        typeSuccessful = False
        for type in ['c',]: #['', 'b']:
            try:
                mov = UfmfV3(os.path.join(directory, 'bgsubtracted_{}{}.ufmf'.format(chunkID, type)))
            except:
                break
            img = np.zeros((1080, 1080), dtype=np.uint8)
            for i, (ts, rois) in tqdm(enumerate(mov.readframes()), total=mov.get_number_of_frames(), leave=False):
                img[:, :] = 0
                blobsForFrame = []
                for ox, oy, roi in rois:
                    rps = skimage.measure.regionprops(
                        skimage.measure.label(roi >= 40))
                    prop = [(x.centroid[0] + oy, x.centroid[1] + ox) + (x.orientation, x.eccentricity, x.area) for x in rps]
                    for p in prop:
                        if p[4] >= 50:
                            blobsForFrame.append(p)
                blobsForFrame = np.array(blobsForFrame)
                blobs.append(blobsForFrame)
            mov.close()
            typeSuccessful = True
            break
        if not typeSuccessful:
            break
    # Done!
    return blobs

def linkDistance(a, b):
    posdiff = np.linalg.norm(a[0:2] - b[0:2])
    anglediff = np.mod(abs(a[2] - b[2]), 2 * np.pi)
    return posdiff + anglediff * 0

def loadManualAnnotations(directory):
    tracesIdx = None
    fnameTraces = directory.replace('\\2020', '\\traces\\2020') + 'list_script_roi.txt'
    if os.path.exists(fnameTraces):
        txtTraces = ''
        with open(fnameTraces, 'r') as f:
            txtTraces = f.read()
        arrTraces = [[float(y) for y in x.split(' ') if len(y.strip()) > 0] for x in txtTraces.split('\n')]
        tracesIdx = {}
        for x in arrTraces:
            if len(x) > 0:
                tracesIdx[int(x[0])] = x

    return tracesIdx

def trajSpan(x):
    return np.linalg.norm(np.max(x, axis=0) - np.min(x, axis=0))

# Filter trajectories by speed
def filterTrajectory(tr, maxStableVel, trajDilation, minTrajSpan):
    # Determine trajectory stability
    v = np.hstack((0, np.linalg.norm(tr[:-1, 1:3] - tr[1:, 1:3], axis=1)))
    isSlow = v <= maxStableVel

    # Allow high-speed movement within two stable segments
    bridgeInstability = 50
    isSlow = np.hstack((np.zeros(bridgeInstability, dtype=bool),
                        isSlow, np.zeros(bridgeInstability, dtype=bool)))
    isSlow = scipy.ndimage.binary_erosion(scipy.ndimage.binary_dilation(
        isSlow, iterations=bridgeInstability), iterations=bridgeInstability)
    isSlow = isSlow[bridgeInstability:(isSlow.size - bridgeInstability)]
    if trajDilation > 0:
        isSlow = scipy.ndimage.binary_dilation(isSlow, iterations=trajDilation)
    elif trajDilation < 0:
        isSlow = scipy.ndimage.binary_erosion(isSlow, iterations=-trajDilation)

    # Is too slow?
    if trajSpan(tr[:, 1:3]) < minTrajSpan:
        isSlow[:] = False

    return tr[isSlow]

def buildTrajectories(blobs, linkThreshold=300, minStableTrajSize=5,
                      maxStableVel=20, trajDilation=1, minTrajSpan=25, maxNum=100):
    trajActive = []
    resultCount = 0

    for fid in tqdm(range(len(blobs) - 6)):
        linksForFrame = {}
        linksForFrameDt = {}

        for dt in range(1, 6):
            # Compute linkage scores for this time delta
            linkages = np.zeros((len(blobs[fid]), len(blobs[fid + dt])), dtype=np.float32)
            for iA, blobA in enumerate(blobs[fid]):
                for iB, blobB in enumerate(blobs[fid + dt]):
                    linkages[iA, iB] = linkDistance(blobA, blobB)

            # Link best-matching blobs in turn
            if linkages.size > 0:
                while np.min(linkages) < linkThreshold:
                    linkA, linkB = np.unravel_index(np.argmin(linkages), linkages.shape)
                    if linkA not in linksForFrame:
                        linksForFrame[linkA] = linkB
                        linksForFrameDt[linkA] = dt
                    linkages[linkA, :] = 99999
                    linkages[:, linkB] = 99999

        # Append new links to existing trajectories if they exist
        for linkA, linkB in linksForFrame.items():
            existingTrajectory = -1
            for i in range(len(trajActive)):
                if trajActive[i][-1][1] == linkA and trajActive[i][-1][0] == fid:
                    existingTrajectory = i
                    break
            dt = linksForFrameDt[linkA]
            if existingTrajectory >= 0:
                trajActive[existingTrajectory].append((fid + dt, linkB))
            else:
                trajActive.append([(fid, linkA), (fid + dt, linkB), ])

        # Now remove outdated trajectories
        while True:
            anythingRemoved = False
            for itraj, tr in enumerate(trajActive):
                if fid - tr[-1][0] >= 10:
                    # Before we add this trajectory to the results,
                    # check that it is valid
                    trXY = np.array([[iframe, ] + blobs[iframe][blobid, :].tolist() for iframe, blobid in tr])
                    trStable = filterTrajectory(trXY, maxStableVel, trajDilation, minTrajSpan)
                    if trStable.shape[0] >= minStableTrajSize:
                        yield trStable
                        resultCount += 1
                        if resultCount >= maxNum and maxNum >= 0:
                            return
                    del trajActive[itraj]
                    anythingRemoved = True
                    break
            if not anythingRemoved:
                break

    # Return remaining trajectories
    for itraj, tr in enumerate(trajActive):
        # Before we add this trajectory to the results,
        # check that it is valid
        trXY = np.array([[iframe, ] + blobs[iframe][blobid, 0:2].tolist() for iframe, blobid in tr])
        trStable = filterTrajectory(trXY, maxStableVel, trajDilation, minTrajSpan)
        if trStable.shape[0] >= minStableTrajSize:
            yield trStable
            resultCount += 1
            if resultCount >= maxNum and maxNum >= 0:
                return
        del trajActive[itraj]
        anythingRemoved = True
        break

def indexTrajectoryPositions(trajectories):
    blobs = {}
    for tri, tr in enumerate(trajectories):
        for k, _fid in enumerate(tr[:, 0]):
            fid = int(_fid)
            if fid in blobs:
                blobs[fid].append(tr[k, 1:3])
            else:
                blobs[fid] = [tr[k, 1:3], ]
    # Convert to NumPy arrays
    for fid in blobs:
        blobs[fid] = np.array(blobs[fid], dtype=float)
    # Done!
    return blobs

def safemax(x):
    return np.max(x) if len(x) > 0 else 0

def computeClassificationMetrics(trajectoriesManual, trajectoriesPred):
    # MATCH_CRITERIA determines how close two markers have to be to be considered
    # a manual-to-predicted match
    MATCH_CRITERIA = 20

    # Get the maximum number of frames present in these trajectories
    numFrames = 1 + int(max(
        safemax([safemax(x[:, 0]) for x in trajectoriesManual]),
        safemax([safemax(x[:, 0]) for x in trajectoriesPred])
    ))

    # Extract blob positions for all trajectories
    blobsManual = indexTrajectoryPositions(trajectoriesManual)
    blobsPred = indexTrajectoryPositions(trajectoriesPred)

    # Count performance metric
    truePos = 0
    falsePos = 0
    falseNeg = 0
    truePosFr = []
    falsePosFr = []
    falseNegFr = []

    # Now determine overlap of manual and predicted positions
    for frameID in tqdm(range(numFrames)):
        markersManual = blobsManual[frameID] if frameID in blobsManual \
            else np.zeros((0, 2), dtype=float)
        markersPred = blobsPred[frameID] if frameID in blobsPred \
            else np.zeros((0, 2), dtype=float)

        # Count true(/false) positives (predicted marker with(/without) nearby manual marker)
        for xy in markersPred:
            if markersManual.shape[0] == 0:
                falsePos += 1
                falsePosFr.append(frameID)
            elif np.min(np.linalg.norm(markersManual - xy, axis=1)) < MATCH_CRITERIA:
                truePos += 1
                truePosFr.append(frameID)
            else:
                falsePos += 1
                falsePosFr.append(frameID)

        # Count false negatives (manual marker without nearby predicted marker)
        for xy in markersManual:
            if markersPred.shape[0] == 0:
                falseNeg += 1
                falseNegFr.append(frameID)
            elif np.min(np.linalg.norm(markersPred - xy, axis=1)) > MATCH_CRITERIA:
                falseNeg += 1
                falseNegFr.append(frameID)

    # Compute summary metrics
    precision = -1
    if truePos + falsePos > 0:
        precision = truePos / (truePos + falsePos)
    recall = -1
    if truePos + falseNeg > 0:
        recall = truePos / (truePos + falseNeg)
    F1 = -1
    if precision + recall > 0:
        F1 = 2 * precision * recall / (precision + recall)

    # Done!
    return F1, precision, recall

@nb.njit(boundscheck=True, nogil=True)
def trajectoryTooCloseToLights(lightsMask, tr, tooClose):
    # Don't run tooClose search parameters that are too large
    if tooClose > 100:
        raise Exception()
    isTooClose = np.zeros(tr.shape[0], dtype=np.bool_)
    for i in range(tr.shape[0]):
        for dx in range(-tooClose, tooClose + 1):
            for dy in range(-tooClose, tooClose + 1):
                if dx ** 2 + dy ** 2 < tooClose ** 2:
                    _x = int(0.5 + tr[i, 1] + dx)
                    _y = int(0.5 + tr[i, 2] + dy)
                    if _x >= 0 and _y >= 0 and _x < lightsMask.shape[0] and _y < lightsMask.shape[1]:
                        if lightsMask[_x, _y]:
                            isTooClose[i] = True
                            break
            if isTooClose[i]:
                break

    return isTooClose

def filterTrajectories(directory, trajectories, tooClose):
    # Compute mask for lights
    keyframes = []
    for i in range(100):
        try:
            mov = UfmfV3(os.path.join(directory, 'raw_{}c.ufmf'.format(i)))
            for k in range(0, 10000, 100):
                keyframes.append(mov.get_keyframe_for_timestamp('mean', 10000 * i + k)[0])
        except Exception as e:
            print(i)
            break
    if len(keyframes) == 0:
        return [x for x in trajectories]
    else:
        keyframeAvg = np.percentile(keyframes, 95, axis=0)
        lightsMask = keyframeAvg > np.percentile(keyframes, 90)

        # Filter trajectories
        newTraj = []
        for tr in trajectories:
            if not np.all(trajectoryTooCloseToLights(lightsMask, tr, tooClose)):
                newTraj.append(tr)
        return newTraj

def run(directory, maxNumChunks = -1, settings = {}):
    # Does this recording have corresponding manual annotations?
    blobsManual = getBlobsManual(directory)
    hasManualAnnotations = (blobsManual is not None)

    # Set missing settings to defaults
    if 'LINK_THRESHOLD' not in settings:
        settings['LINK_THRESHOLD'] = 50
    if 'MIN_TRAJ_DUR' not in settings:
        settings['MIN_TRAJ_DUR'] = 10
    if 'MAX_LANDED_VEL' not in settings:
        settings['MAX_LANDED_VEL'] = 50
    if 'TRAJ_DILATION' not in settings:
        settings['TRAJ_DILATION'] = 0
    if 'MIN_TRAJ_SPAN' not in settings:
        settings['MIN_TRAJ_SPAN'] = 50

    # Build trajectories for manual annotations, if present
    trajectoriesManual = []
    if blobsManual is not None:
        for tr in buildTrajectories(blobsManual,
                linkThreshold=settings['LINK_THRESHOLD'],
                minStableTrajSize=settings['MIN_TRAJ_DUR'],
                maxStableVel=settings['MAX_LANDED_VEL'],
                trajDilation=settings['TRAJ_DILATION'],
                minTrajSpan=settings['MIN_TRAJ_SPAN'],
                maxNum=-1):
            trajectoriesManual.append(tr)

    # Get blobs for each frame
    blobsPred = getBlobs(directory, maxNumChunks = maxNumChunks)

    # Get Trajectories
    trajectoriesPred = []
    for tr in buildTrajectories(blobsPred,
                linkThreshold=settings['LINK_THRESHOLD'],
                minStableTrajSize=settings['MIN_TRAJ_DUR'],
                maxStableVel=settings['MAX_LANDED_VEL'],
                trajDilation=settings['TRAJ_DILATION'],
                minTrajSpan=settings['MIN_TRAJ_SPAN'],
                maxNum=-1):
        trajectoriesPred.append(tr)

    # Now filter out trajectories too close to the lights
    trajectoriesManualFilt = filterTrajectories(directory, trajectoriesManual, 20)
    trajectoriesPredFilt = filterTrajectories(directory, trajectoriesPred, 20)

    # Compute classification performance, if manual annotations were present
    F1, precision, recall = computeClassificationMetrics(
        trajectoriesManual, trajectoriesPred)
    F1_filt, precision_filt, recall_filt = computeClassificationMetrics(
        trajectoriesManualFilt, trajectoriesPredFilt)

    # Store results
    fnameResults = os.path.join(directory, 'trajectories_{}_{}_{}_{}_{}.pickle'.format(
        settings['LINK_THRESHOLD'],
        settings['MIN_TRAJ_DUR'],
        settings['MAX_LANDED_VEL'],
        settings['TRAJ_DILATION'],
        settings['MIN_TRAJ_SPAN']))
    jl.dump({
        'F1': F1,
        'precision': precision,
        'recall': recall,
        'F1_filt': F1_filt,
        'precision_filt': precision_filt,
        'recall_filt': recall_filt,
        'trajectoriesManual': trajectoriesManual,
        'trajectoriesPred': trajectoriesPred,
        'trajectoriesManualFilt': trajectoriesManualFilt,
        'trajectoriesPredFilt': trajectoriesPredFilt,
        'hasManualAnnotations': hasManualAnnotations
    }, fnameResults)

def run():
    # Get recording directories
    from MosquitoFrameReader import getRecordingDirectories
    directories = getRecordingDirectories()
    for d in directories:
        print('Found directory: {}'.format(d))
    time.sleep(1)
    # Configurations to run
    settings = []
    for LINK_THRESHOLD, MIN_TRAJ_SPAN, MAX_LANDED_VEL in [(50, 20, 20), ]:
        for d in directories:
            settings.append((
                d, {
                    'LINK_THRESHOLD': LINK_THRESHOLD,
                    'MAX_LANDED_VEL': MAX_LANDED_VEL,
                    'MIN_TRAJ_SPAN': MIN_TRAJ_SPAN,
                }))
    # Process in parallel
    jl.Parallel(n_jobs=40)(jl.delayed(run)(s[0], -1, s[1]) for s in settings)
    
if __name__ == "__main__":
    run()