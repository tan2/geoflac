#!/usr/bin/python3

import pandas as pd
import flac


def read_all_markers(start=0, end=-1):
    """Read markers for all frames and return a single pandas DataFrame.

    The returned DataFrame contains markers from all frames concatenated
    together. A column named 'frame' identifies which frame each marker
    belongs to. The columns are:
      ['frame', 'x', 'z', 'age', 'phase', 'ID', 'a1', 'a2', 'ntriag']
    """
    fl = flac.Flac()
    frames = fl.frames.astype(int)
    dfs = []
    for i in frames[start:end]:
        x, z, age, phase, ID, a1, a2, ntriag = fl.read_markers(i)
        if len(x) == 0:
            raise RuntimeError('No markers found in frame %d' % i)
        df = pd.DataFrame({
            'frame': int(i),
            'x': x,
            'z': z,
            'age': age,
            'phase': phase,
            'ID': ID,
            'a1': a1,
            'a2': a2,
            'ntriag': ntriag,
        })
        dfs.append(df)

    if len(dfs) == 0:
        # no markers at all
        raise RuntimeError('No frames')

    # Concatenate and reset index so we have tidy DataFrame
    result = pd.concat(dfs, ignore_index=True)
    # ensure column ordering: put 'frame' first
    cols = ['frame', 'x', 'z', 'age', 'phase', 'ID', 'a1', 'a2', 'ntriag']
    result = result[cols]
    return result


if __name__ == '__main__':
    startframe = 0
    lastframe = -1
    mks = read_all_markers(startframe, lastframe)