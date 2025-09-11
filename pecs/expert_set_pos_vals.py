'''
expert_set_pos_vals.py

This script takes in a csv file and sets positioner values contained in the file.
WARNING: This script starts up its own posstates and can commit to the DB. This
can cause undesirable interactions with an active instance. DO NOT use if an instance
is up. DO NOT use if you can you set_calibrations.py to do the same task.
'''

'''
2020/11/12 Discussion between KF and JHS led to creation of tool to dig us out
           of a hole where petals won't start because of bad OFFSET_X/Y
'''
allowed_keys = {'OFFSET_X', 'OFFSET_Y'}
allowed_date = '2020/11/12' # Change this when changing above

import posstate
from DOSlib.positioner_index import PositionerIndex
from DBSingleton import DBSingleton
import os
import pandas as pd
import argparse
import posconstants as pc

index = PositionerIndex()
altered_states = {} #contains sets of states keyed by ptlid
altered_calib_states = {} #contains sets of states keyed by ptlid

def add_to_altered_states(ptlid):
    if ptlid not in set(altered_states.keys()):
        altered_states[ptlid] = set()
    return lambda state: altered_states[ptlid].add(state)
def add_to_altered_calib_states(ptlid):
    if ptlid not in set(altered_calib_states.keys()):
        altered_calib_states[ptlid] = set()
    return lambda state: altered_calib_states[ptlid].add(state)

def get_posstates(posids):
    posstates = {}
    for posid in posids:
        ids = index.find_by_device_id(posid)
        assert ids['DEVICE_TYPE'] in ['POS','ETC'], 'POS_ID must be a positioner!'
        ptlid = ids['PETAL_ID']
        posstates[posid] = posstate.PosState(posid, logging=False,
                                             device_type='pos', petal_id=ptlid,
                                             alt_move_adder=add_to_altered_states(ptlid),
                                             alt_calib_adder=add_to_altered_calib_states(ptlid))
    return posstates

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-i', '--infile', type=str, required=True, help='path to input csv file')
uargs = parser.parse_args()

data = pd.read_csv(uargs.infile)
assert 'POS_ID' in data.columns, 'POS_ID must be present in the input file'

input('WARNING: this script may have undesirable interactions with an active instance. DO NOT use with an active instance. DO NOT use in place of set_calibrations.py. Hit enter to agree.')
initials = input('Please enter your INITIALS for logging purposes: ')
comment = input('Please give a short comment why this script is being used: ')
input(f'Initials are {initials} and comment is {comment}. Hit enter to confirm.')

columns = set()
for key in data.columns:
    if not(pc.is_constants_key(key)) and key != 'POS_ID': #Try to filter keys, not a complete filter
        columns.add(key)

disallowed = columns - allowed_keys
columns_to_update = columns & allowed_keys

assert len(columns_to_update) != 0, 'No columns to update'
print(f'Accepted columns: {columns_to_update}')
if disallowed:
    print(f'Columns {disallowed} will *not* be changed. Only {allowed_keys} are allowed.'
          f'This is our current policy (as of {allowed_date}) in order to restrict ' 
          f'the power/risk of this script. The more general, "proper" method of ' 
          f'updating calibration settings is via set_calibrations.py. Procedures ' 
          f'are provided in DESI-5732. Additionally, certain limited parameters ' 
          f'may be safely set at the PETAL console. Enter command "readme" there ' 
          f'for more information.')

posids = set(data['POS_ID'])
posstates = get_posstates(posids)

log_note = f'expert_set_pos_vals; Initials {initials}; reason {comment}'
print(f'Log note will be {log_note}')

# Store values in states
for i, row in data.iterrows():
    altered_calib = False
    altered_move = False
    posid = row['POS_ID']
    for key in columns_to_update:
        if not pd.isnull(row[key]):
            success = posstates[posid].store(key, row[key])
            if success:
                if pc.is_calib_key(key):
                    altered_calib = True
                else:
                    altered_move = True
    if altered_move:
        posstates[posid].store('LOG_NOTE', log_note)
    if altered_calib:
        posstates[posid].store('CALIB_NOTE', log_note)

# Commit values to DB
ptlids = set(altered_states.keys()) | set(altered_calib_states.keys())
os.environ['DOS_POSMOVE_WRITE_TO_DB'] = 'True'
for ptlid in ptlids:
    posmoveDB = DBSingleton(petal_id=ptlid)
    if ptlid in set(altered_states.keys()):
        if altered_states[ptlid]: #make sure its not an empty set
            posmoveDB.WriteToDB(altered_states[ptlid], ptlid, 'pos_move')
    if ptlid in set(altered_calib_states.keys()):
        if altered_calib_states[ptlid]: #make sure its not an empty set
            posmoveDB.WriteToDB(altered_calib_states[ptlid], ptlid, 'pos_calib')
