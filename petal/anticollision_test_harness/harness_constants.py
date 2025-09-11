import os
import cProfile
import pstats
import sys
sys.path.append(os.path.abspath('../../petal/'))
import posconstants as pc
import time
import numpy as np
import math

device_locations_path = '../positioner_locations_0530v14.csv'

pos_dir = os.path.join(pc.FP_SETTINGS_PATH, 'harness_settings', 'pos_parameter_sets')
req_dir = os.path.join(pc.FP_SETTINGS_PATH, 'harness_settings', 'move_request_sets')
pos_prefix = 'params'
req_prefix = 'requests'

def bool(s):
    if not s or s in {'False','0','None','FALSE','NONE','false','none'}:
        return False
    return True

data_types = {'POS_ID':str, 'DEVICE_LOC':int, 'CTRL_ENABLED':bool,
              'LENGTH_R1':float, 'LENGTH_R2':float, 'OFFSET_T':float,
              'OFFSET_P':float, 'OFFSET_X':float, 'OFFSET_Y':float,
              'PHYSICAL_RANGE_T':float, 'PHYSICAL_RANGE_P':float,
              'command':str, 'u':float, 'v':float,
              'KEEPOUT_EXPANSION_PHI_RADIAL':float,
              'KEEPOUT_EXPANSION_PHI_ANGULAR':float,
              'KEEPOUT_EXPANSION_THETA_RADIAL':float,
              'KEEPOUT_EXPANSION_THETA_ANGULAR':float,
              'CLASSIFIED_AS_RETRACTED':bool,
              }

def filenumber(filename, prefix):
    base = os.path.basename(filename)
    num_str = os.path.splitext(base)[0].split(prefix)[1]
    return int(num_str)

def filename(prefix, suffix=None):
    if isinstance(suffix, int):
        suffix = f'{suffix:05}'
    suffix = f'_{suffix}' if suffix != None else ''
    return f'{prefix}{suffix}.csv'

def filepath(directory, prefix, integer=None):
    return os.path.join(directory, filename(prefix, integer))

def _make_prefix(lead, ptlnum, setnum=None):
    assert isinstance(ptlnum, int)
    prefix = f'{lead}_ptl{ptlnum:02}'
    if setnum != None:
        assert isinstance(setnum, int)
        prefix += f'_set{setnum:02}'
    return prefix

def make_request_prefix(ptlnum, setnum):
    return _make_prefix(req_prefix, ptlnum, setnum)

def make_params_prefix(ptlnum):
    return _make_prefix(pos_prefix, ptlnum)

def compact_timestamp(nowtime=None, basetime=1582915648):
    '''Compact, readable time code. Default return string is six characters
    in length; will exceed this length at basetime + 69 years. Precision is
    rounded to seconds. Default argument baselines it at a recent time on
    Feb 28, 2020, 10:47 AM PST. The argument nowtime is just there for testing.
    '''
    maxchar = 6
    nowtime = time.time() if not nowtime else nowtime
    relative_now = math.floor(nowtime - basetime)
    converted = np.base_repr(relative_now, base=36)
    padded = converted.rjust(maxchar,'0') if len(converted) < maxchar else converted
    return padded

# Timing profiler wrapper function
statsfile = os.path.join(pc.dirs['temp_files'],'stats_harness')
def profile(evaluatable_string, sorting='cumtime', n_lines=20):
    assert sorting in {'ncalls', 'tottime', 'percall', 'cumtime', 'percall'}
    print(evaluatable_string)
    cProfile.run(evaluatable_string,statsfile)
    p = pstats.Stats(statsfile)
    p.strip_dirs()
    p.sort_stats(sorting)
    p.print_stats(n_lines)
    return p