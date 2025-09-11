'''Reads in a table of measured (x,y) values for positioners. Triggers the
internal function test_and_update_TP. (That function calculates new posintT,
posintP --- same as POS_T, POS_P --- using the current calibration parameters
in the online database. Then it stores the new theta and phi values to the db.)
This is a powerful function with risk of major operational errors if used
incorrectly. Only for focal plane experts, and only to be used with consensus
of the focal plane team.
'''

import os
script_name = os.path.basename(__file__)

# input file format
required_keys = ['DEVICE_ID', 'X_FP', 'Y_FP']
xy_keys = ['X_FP', 'Y_FP']
xy_remap = {'X_FP': 'FPA_X',
            'Y_FP': 'FPA_Y'}

# command line argument parsing
import argparse
doc = f'{__doc__} Input CSV file should contain the following columns: {required_keys}'
parser = argparse.ArgumentParser(description=doc)
parser.add_argument('-i', '--infile', type=str, required=True, help='path to input csv file')
args = parser.parse_args()

# read input data
import pandas as pd
frame = pd.read_csv(args.infile)

# set up a logger
import logging
import time
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
[logger.removeHandler(h) for h in logger.handlers]
formatter = logging.Formatter(fmt='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S %z')
formatter.converter = time.gmtime
sh = logging.StreamHandler()
sh.setFormatter(formatter)
logger.addHandler(sh)
logger.info(f'Running {script_name} to set positioner POS_T and POS_P.')
logger.info(f'Input file is: {args.infile}')
logger.info(f'Input table contains {len(frame)} rows')
    
# validate the table format
for key in required_keys:
    assert key in frame.columns, f'Error no {key} column found'

# some basic data validation, not a guarantee of quality of parameters
import numpy as np
for key in xy_keys:
    assert all(np.isfinite(frame[key])), 'Error non-finite input data'
posids = set(frame['DEVICE_ID'])
assert len(posids) == len(frame), 'Error multiple rows for same device'

# delete extraneous columns. I do this to clarify any debugging situations
for key in frame.columns:
    if key not in required_keys:
        del frame[key]

# set up pecs (access to online system)
from pecs import PECS
pecs = PECS(interactive=False)
logger.info(f'PECS initialized, discovered PC ids {pecs.pcids}')
pecs.ptl_setup(pecs.pcids, posids=posids)

# gather some interactive information from the user
user = ''
while not user:
    user = input('For the log, please enter your NAME or INITIALS:\n >> ')
comment = ''
while not comment:
    comment = input('For the log, enter any additional COMMENT, giving context/rationale for posting these new values:\n >> ')

# posconstants import
try:
    import posconstants as pc
except:
    import sys
    path_to_petal = '../petal'
    sys.path.append(os.path.abspath(path_to_petal))
    logger.info('Couldn\'t find posconstants the usual way, resorting to sys.path.append')
    import posconstants as pc

# perform the update
log_note = pc.join_notes(script_name, f'user {user}', f'comment {comment}', f'input_file {args.infile}')
frame.rename(columns=xy_remap, inplace=True)
kwargs = {'mode': 'posTP',
          'tp_updates_tol': 0,
          'tp_updates_fraction': 1,
          'auto_update': True,
          'verbose': True,
          'log_note': log_note,
          'commit': True,
          }
logger.info(f'Now performing test_and_update_tp on positioners {posids}.\nInput values:\n{frame}\n{kwargs}')
updates = pecs.ptlm.test_and_update_TP(measured_data=frame, **kwargs)
logger.info(f'Result data:\n{updates}')
#logger.info(f'DEBUG Result data to_string():\n{updates.to_string()}')
#logger.info(f'DEBUG Result data columns:\n{updates.columns}')
#updates.to_csv('debug_test.csv')
logger.info('Complete. Please double-check online db to confirm results are what you expected.')
