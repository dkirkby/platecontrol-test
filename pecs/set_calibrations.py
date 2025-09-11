'''Store fiber positioner calibration values to the online database. This is a
powerful function with risk of major operational errors if used incorrectly.
Only for focal plane experts, and only to be used with consensus of the focal
plane team.
'''

import os
script_name = os.path.basename(__file__)

# input file format
format_info = 'For data model and procedures to generate these values, see DESI-5732.'
valid_keys = {'LENGTH_R1', 'LENGTH_R2', 'OFFSET_T', 'OFFSET_P', 'OFFSET_X',
              'OFFSET_Y', 'PHYSICAL_RANGE_T', 'PHYSICAL_RANGE_P',
              'GEAR_CALIB_T', 'GEAR_CALIB_P', 'SCALE_T', 'SCALE_P',
              'DEVICE_CLASSIFIED_NONFUNCTIONAL','CLASSIFIED_AS_RETRACTED',
              'POS_T', 'POS_P', 'LOC_T', 'LOC_P',
              'KEEPOUT_EXPANSION_THETA_RADIAL', 'KEEPOUT_EXPANSION_PHI_RADIAL',
              'KEEPOUT_EXPANSION_THETA_ANGULAR', 'KEEPOUT_EXPANSION_PHI_ANGULAR',
              'ZENO_MOTOR_P', 'ZENO_MOTOR_T', 'SZ_CW_P', 'SZ_CCW_P',  'SZ_CW_T', 'SZ_CCW_T'}
fit_err_keys = {'FIT_ERROR_STATIC', 'FIT_ERROR_DYNAMIC', 'FIT_ERROR',
                'NUM_POINTS_IN_FIT_STATIC', 'NUM_POINTS_IN_FIT_DYNAMIC',
                'NUM_OUTLIERS_EXCLUDED_STATIC', 'NUM_OUTLIERS_EXCLUDED_DYNAMIC'}
commit_prefix = 'COMMIT_'
commit_keys = {key: commit_prefix + key for key in valid_keys}
boolean_keys = set(commit_keys.values()) | {'DEVICE_CLASSIFIED_NONFUNCTIONAL', 'CLASSIFIED_AS_RETRACTED',}
float_keys = (valid_keys | fit_err_keys) - boolean_keys
no_nominal_val = {'DEVICE_CLASSIFIED_NONFUNCTIONAL', 'CLASSIFIED_AS_RETRACTED', 'POS_P', 'POS_T',
                  'KEEPOUT_EXPANSION_THETA_RADIAL', 'KEEPOUT_EXPANSION_PHI_RADIAL',
                  'KEEPOUT_EXPANSION_THETA_ANGULAR', 'KEEPOUT_EXPANSION_PHI_ANGULAR',
                  'ZENO_MOTOR_P', 'ZENO_MOTOR_T', 'SZ_CW_P', 'SZ_CCW_P',  'SZ_CW_T', 'SZ_CCW_T'}
def dbkey_for_key(key):
    '''Maps special cases of keys that may have different terminology in input file
    online database.'''
    remap = {'SCALE_T': 'GEAR_CALIB_T',
             'SCALE_P': 'GEAR_CALIB_P',
             'COMMIT_SCALE_T': 'COMMIT_GEAR_CALIB_T',
             'COMMIT_SCALE_P': 'COMMIT_GEAR_CALIB_P',
             'LOC_T': 'POS_T',
             'LOC_P': 'POS_P',
             'COMMIT_LOC_T': 'COMMIT_POS_T',
             'COMMIT_LOC_P': 'COMMIT_POS_P',}
    if key in remap:
        return remap[key]
    return key

# command line argument parsing
import argparse
doc = f'{__doc__} Input CSV file must contain a POS_ID column. Input should'
doc += f' contain all or any subset of the following parameter columns: {valid_keys}.'
doc +=  ' For every parameter column there must be a corresponding boolean column'
doc += f' prefixed with "{commit_prefix}", stating in each row whether that value is'
doc += f' intended to be saved to the online db. {format_info}'
parser = argparse.ArgumentParser(description=doc)
parser.add_argument('-i', '--infile', type=str, required=True, help='path to input csv file')
parser.add_argument('-s', '--simulate', action='store_true', help='perform the script "offline" with no storing to memory or database')
parser.add_argument('-t', '--talkative', action='store_true', help='be talkative')
args = parser.parse_args()

# read input data
from astropy.table import Table
table = Table.read(args.infile)

# set up a log file
import simple_logger
try:
    import posconstants as pc
except:
    import os, sys
    path_to_petal = '../petal'
    sys.path.append(os.path.abspath(path_to_petal))
    print('Couldn\'t find posconstants the usual way, resorting to sys.path.append')
    import posconstants as pc

# yes, as of 2020-06-17, I'm saving it in two places
# to be absolutely sure we have a record
# can cut back to one, once we are confident of logs getting properly saved at KPNO
log_dirs = [os.path.dirname(args.infile), pc.dirs['calib_logs']]
log_name = pc.filename_timestamp_str() + '_set_calibrations.log'
log_paths = [os.path.join(d, log_name) for d in log_dirs]
logger, _, _ = simple_logger.start_logger(log_paths)
logger.info(f'Running {script_name} to set positioner calibration parameters.')
logger.info(f'Input file is: {args.infile}')
logger.info(f'Table contains {len(table)} rows')
if args.simulate:
    logger.info('Running in simulation mode. No data will stored to petal memory nor database.')
assert2 = simple_logger.assert2
input2 = simple_logger.input2

# deal with astropy's idiotic handling of booleans as strings
for key in boolean_keys & set(table.columns):
    table[key] = [pc.boolean(x) for x in table[key]]

# deal with astropy's annoying restrictions on integer values
truefalse = {'TRUE': True, 'FALSE': False}

for key in float_keys & set(table.columns):
    table[key] = [float(x) if x not in truefalse else truefalse[x] for x in table[key]]

# validate the table format
if args.talkative:
    logger.info(f'Columns: "{table.columns}"')
assert2('POS_ID' in table.columns, 'No POS_ID column found in input table')
for key in set(table.columns):  # set op provides a copy
    remapped_key = dbkey_for_key(key)
    if key != remapped_key:
        remap_conflict = key in table.columns and remapped_key in table.columns
        assert2(not(remap_conflict), f'Columns {key} and {remapped_key} conflict in the input file (two columns with same meaning, so cannot disambiguate)')
        table.rename_column(key, dbkey_for_key(key))
        logger.warning(f'For compatibility with online DB, renamed column {key} to {remapped_key}')
keys = valid_keys & set(table.columns)
if args.talkative:
    logger.info(f'keys: "{keys}"')
assert2(len(keys) > 0, 'No valid parameter columns found in input table')
no_commit_key = {key for key in keys if commit_keys[key] not in table.columns}
assert2(len(no_commit_key) == 0, f'Input table params {no_commit_key} lack {commit_prefix}-prefixed fields. {format_info}')
logger.info('Checked formatting of input table')

# some basic data validation (type and bounds checking)
# not a guarantee of quality of parameters
import numpy as np
requested_posids = set()
for key in keys:
    column = table[key]
    commit_key = commit_keys[key]
    commit_type_ok = table[commit_key].dtype in [int, bool]     # [np.int, np.bool]
    assert2(commit_type_ok, f'{commit_key} data type must be boolean or integer representing boolean')
    data_type_ok = column.dtype in [int, float, bool]    # [np.int, np.float, np.bool]
    assert2(data_type_ok, f'{key} data type must be numeric or boolean')
    commit_requested = table[commit_key]
    def assert_ok(is_valid_array, err_desc):
        rng = range(len(table))
        cannot_commit = [table['POS_ID'][i] for i in rng if not is_valid_array[i] and commit_requested[i]]
        assert2(not(any(cannot_commit)), f'Input table rejected: {key} contains {err_desc} for posid(s) {cannot_commit}')
    isfinite = np.isfinite(column)
    assert_ok(isfinite, 'non-finite value(s)')
    if key not in no_nominal_val:
        nom = pc.nominals[key]
        nom_min = nom['value'] - nom['tol']
        nom_max = nom['value'] + nom['tol']
        column[isfinite == False] = np.inf  # dummy value to suppress astropy warnings when performing >= or <= ops below, on values that we don't plan to commit anyway. ok to do here because already checked commit values for finitude above
        assert_ok(column >= nom_min, f'value(s) below limit {nom_min}')
        assert_ok(column <= nom_max, f'value(s) above limit {nom_max}')
    requested_posids |= set(table['POS_ID'][commit_requested])
logger.info('Checked data types and bounds')

# set up pecs (access to online system)
try:
    from pecs import PECS
    pecs = PECS(interactive=False, fvc=False, no_expid=True)
    logger.info(f'PECS initialized, discovered PC ids {pecs.pcids}')
    posids = set(table['POS_ID'])
    pecs.ptl_setup(pecs.pcids, posids=posids)
    pecs_on = True
except Exception as e:
    logger.warning(f'PECS initialization failed: {e}')
    pecs_on = False

# gather some interactive information from the user
user = ''
while not user:
    user = input2('For the log, please enter your NAME or INITIALS:')
archive_ref = ''
while not archive_ref:
    archive_ref = input2('For the log, enter DESI-XXXX DOCUMENT where this input csv table is archived:')
comment = ''
while not comment:
    comment = input2('For the log, enter any additional COMMENT, giving context/rationale for posting these new values:')

# interactive human checks, since it is important to get everything right
if pecs_on:
    # Most checks should be done ahead of time with preparation script in desimeter.
    # Here we would only check about stuff that is known to the online system.
    # Like maybe validating the posids, etc. Might not be essential to do so here.
    pass
else:
    #logger.warning('Skipping interactive checks, since PECS not initialized')
    pass

# store the data
logger.info(f'Storing data to memory (not yet to database) for {len(table)} positioners.')
any_stored = False
for row in table:
    posid = row['POS_ID']
    role = -1 if args.simulate else pecs.ptl_role_lookup(posid)
    stored = {}
    for key in keys:
        if row[commit_keys[key]]:
            value = row[key]
            kwargs = {'device_id': posid,
                      'key': key,
                      'value': value,
                      'participating_petals': role
                      }
            updates = [kwargs]
            if key == 'DEVICE_CLASSIFIED_NONFUNCTIONAL':
                kwargs['comment'] = row['ENABLE_DISABLE_RATIONALE']
                fiber_intact = True if args.simulate else pecs.ptlm.get_posfid_val(posid, 'FIBER_INTACT', participating_petals=role)
                ctrl_enabled = True if args.simulate else pecs.ptlm.get_posfid_val(posid, 'CTRL_ENABLED', participating_petals=role)
                new_ctrl_enabled = not(value) & fiber_intact
                if ctrl_enabled != new_ctrl_enabled:
                    kwargs2 = kwargs.copy()
                    kwargs2['key'] = 'CTRL_ENABLED'
                    kwargs2['value'] = new_ctrl_enabled
                    kwargs2['comment'] = f'auto-{"enabled" if new_ctrl_enabled else "disabled"} by set_calibrations.py upon setting {key}={value}'
                    updates = [kwargs, kwargs2]
            elif key == 'CLASSIFIED_AS_RETRACTED':
                kwargs['comment'] = row['ENABLE_DISABLE_RATIONALE']
            for kwargs in updates:
                val_accepted = True if args.simulate else pecs.ptlm.set_posfid_val(**kwargs)
                if val_accepted:
                    stored[kwargs['key']] = kwargs['value']
                else:
                    logger.error(f'set_posfid_val({kwargs})')
    if stored:
        # [JHS] Would be nice to include analysis metadata fields in the note, drawn
        # from the input table. Presumably when that table is in ecsv format.
        note = pc.join_notes(script_name, f'user {user}', f'comment {comment}',
                                 f'input_file {args.infile}', f'archive_ref {archive_ref}',
                                 f'params {stored}')
        for key in fit_err_keys:
            if key in row.columns:
                note = pc.join_notes(note, f'{key.lower()} {row[key]}')
        if not args.simulate:
            pecs.ptlm.set_posfid_val(posid, 'CALIB_NOTE', note, participating_petals=role)
        logger.info(f'{posid}: {stored}')
        any_stored = True

# commit to online database
if any_stored:
    logger.info('Committing the data set to online database.')
    if not args.simulate:
        pecs.ptlm.commit(mode='both')
    logger.info('Commit complete.')
else:
    logger.warning('No data found to commit. Nothing will be changed in the online db.')
logger.info('Please remember to archive the input CSV file, and the .log file' +
            ' (generated in the same folder) to docdb.')
