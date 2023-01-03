import os, glob, re, pandas as pd, joblib as jl, datetime
from tqdm import tqdm

def trysearch(p, s):
    try:
        return re.search(p, s).groups()
    except Exception as e:
        return []

def run(LINK_THRESHOLD, MIN_TRAJ_SPAN, MAX_LANDED_VEL):
    # Find data directories
    trajFname = 'trajectories_{}_{}_{}_{}_{}.pickle'.format(
        LINK_THRESHOLD,
        10,
        MAX_LANDED_VEL,
        0,
        MIN_TRAJ_SPAN
    )
    dirs = glob.glob('Z:\\Abel\\MOSQUITOS\\*\\*\\*\\{}'.format(trajFname)) + \
           glob.glob('Z:\\Abel\\MOSQUITOS\\*\\*\\*\\*\\{}'.format(trajFname))

    # Save multiple versions
    for filt in [True, False]:
        for startendonly in [False, True]:
            tblAll = []
            for d in tqdm(dirs):
                data = jl.load(d)['trajectoriesPredFilt' if filt else 'trajectoriesPred']
                data = [x[:, 0:3] for x in data]

                # Process timestamps
                _dates = [[int(x) for x in trysearch(
                    '([0-9]{2})-([0-9]{2})-([0-9]{4})_([0-9]{2})-([0-9]{2})-([0-9]{2})',
                        y)] for y in glob.glob(os.path.join(os.path.dirname(d), '*.h264'))]
                _dates = [a for a in _dates if len(a) == 6]
                dates = [datetime.datetime(
                    day=d, month=mo, year=y, hour=h, minute=m, second=s) for mo, d, y, h, m, s in _dates]
                startTime = min(dates)

                # Process metadata
                experimentYear, experimentType, experimentDate, experimentPos, _ = [
                    x.replace('\\', '') if x is not None else '' for x in re.search(
                        '(\\\\[0-9]{4})(\\\\[a-zA-Z0-9 -]*)?(\\\\[0-9-]*)(\\\\[a-zA-Z0-9 -_]*)(\\\\[a-zA-Z0-9 -_]*)',
                            d).groups()]

                # Store
                for landingid, tr in enumerate(data):
                    for relfr, (fr, x, y) in enumerate(tr):
                        if (not startendonly) or (relfr == 0 or relfr == len(tr) - 1):
                            dt = (startTime + datetime.timedelta(seconds=fr)).strftime('%Y-%m-%d %H-%M-%S')
                            tblAll.append((experimentYear, experimentType, experimentDate,
                                           experimentPos, landingid, int(relfr), int(fr),
                                           dt, int(0.5 + x), int(0.5 + y)))

            runName = os.path.basename(dirs[0]).replace('.pickle', '').replace('trajectories_', '')

            tbl = pd.DataFrame(tblAll, columns=['experiment_year', 'experiment_type', 'experiment_date',
                'participant_position', 'landing_id', 'frame_rel', 'frame', 'timestamp', 'x', 'y'])
            tbl.to_csv('tracking-predictions-aggregated{}{}-{}.csv'.format(
                '-startendonly' if startendonly else '', '-filtered' if filt else '', runName))

if __name__ == "__main__":
    run(50, 20, 20)
