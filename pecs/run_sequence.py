# -*- coding: utf-8 -*-
"""
Perform a sequence of moves + FVC measurements on hardware. (Note that any positioner
calibration parameters or motor settings altered during the sequence are stored only
in memory, so a re-initialization of the petal software should restore you to the
original state, as defined by the posmove and constants databases. Also, parking moves
will include final antibacklash rotations, regardless of the antibacklash settings
for the main body of the sequence.)
"""

import os
script_name = os.path.basename(__file__)

# command line argument parsing
import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-i', '--infile', type=str, required=True, help='Path to sequence file, readable by sequence.py. For some common pre-cooked sequences, run sequence_generator.py.')
parser.add_argument('-a', '--anticollision', type=str, default='adjust', help='anticollision mode, can be "adjust", "adjust_requested_only", "freeze" or None. Default is "adjust"')
parser.add_argument('-p', '--enable_phi_limit', action='store_true', help='turns on minimum phi limit for move targets, default is False')
parser.add_argument('-r', '--match_radius', type=int, default=None, help='int, specify a particular match radius, other than default')
parser.add_argument('-u', '--check_unmatched', action='store_true', help='turns on auto-disabling of unmatched positioners, default is False')
parser.add_argument('-t', '--test_tp', action='store_true', help='turns on auto-updating of POS_T, POS_P based on measurements, default is False')
parser.add_argument('-x', '--no_movement', action='store_true', help='for debugging purposes, this option suppresses sending move tables to positioners, so they will not physically move')
default_cycle_time = 60.0 # sec
parser.add_argument('-c', '--cycle_time', type=float, default=default_cycle_time, help=f'min period of time (seconds) for successive moves (default={default_cycle_time})')
max_fvc_iter = 10
parser.add_argument('-nm', '--num_meas', type=int, default=1, help=f'int, number of measurements by the FVC per move (default is 1, max is {max_fvc_iter})')
max_corr = 5
parser.add_argument('-nc', '--num_corr', type=int, default=0, help=f'int, number of correction moves, for those rows where allowed in sequence definition (default is 0, max is {max_corr})')
default_n_best = 5
default_n_worst = 5
parser.add_argument('-nb', '--num_best', type=int, default=default_n_best, help=f'int, number of best performers to display in log messages for each measurement (default is {default_n_best})')
parser.add_argument('-nw', '--num_worst', type=int, default=default_n_worst, help=f'int, number of worst performers to display in log messages for each measurement (default is {default_n_worst})')
park_options = ['posintTP', 'poslocTP', 'None', 'False']
default_park = park_options[0]
parser.add_argument('-prep', '--prepark', type=str, default=default_park, help=f'str, controls initial parking move, prior to running the sequence. Valid options are: {park_options}, default is {default_park}')
parser.add_argument('-post', '--postpark', type=str, default=default_park, help=f'str, controls final parking move, after running the sequence. Valid options are: {park_options}, default is {default_park}')
parser.add_argument('-ms', '--start_move', type=int, default=0, help='start the test at this move index (defaults to move 0)')
parser.add_argument('-mf', '--final_move', type=int, default=-1, help='finish the test at this move index (or defaults to last row)')
parser.add_argument('-curr', '--motor_current', type=int, default=None, help='set motor currents (duty cycle) for the duration of the test. Overrides any values in the online system or sequence file. Must be an integer between 1 and 100')
parser.add_argument('-v', '--verbose', action='store_true', help='turn on verbosity at terminal window (note that the log file on disk will always be verbose)')
parser.add_argument('-cf', '--cycle_fiducials', action='store_true', help='Set flag to cycle fiducials between moves like in an OCS sequence (they are left on between submoves - it is ok to turn fiducials on or off with the pecs prompt in the beginning)')
parser.add_argument('-ncp','--no_correction_pause', action='store_true', help='Set this flag to skip pauses for correction moves (example use is to act more like on-mountain moves)')

uargs = parser.parse_args()
if uargs.anticollision == 'None':
    uargs.anticollision = None
assert uargs.anticollision in {'adjust', 'adjust_requested_only', 'freeze', None}, f'bad argument {uargs.anticollision} for anticollision parameter'
assert 1 <= uargs.num_meas <= max_fvc_iter, f'out of range argument {uargs.num_meas} for num_meas parameter'
assert 0 <= uargs.num_corr <= max_corr, f'out of range argument {uargs.num_corr} for num_corr parameter'
assert uargs.prepark in park_options, f'invalid park option, must be one of {park_options}'
assert uargs.postpark in park_options, f'invalid park option, must be one of {park_options}'
uargs.prepark = None if uargs.prepark in ['None', 'False'] else uargs.prepark
uargs.postpark = None if uargs.postpark in ['None', 'False'] else uargs.postpark
assert uargs.motor_current == None or 0 < uargs.motor_current <= 100, f'out of range argument {uargs.motor_current} for motor current parameter'
cycle_fids = uargs.cycle_fiducials
correction_pause = not(uargs.no_correction_pause)

# read sequence file
import sequence
seq = sequence.Sequence.read(uargs.infile)

# check start / finish moves
try:
    start_move = seq.index(seq[uargs.start_move])  # normalizes negative index cases
    final_move = seq.index(seq[uargs.final_move])  # normalizes negative index cases
except:
    assert False, f'start_move={uargs.start_move} or final_move={uargs.final_move} is out of range of the sequence (which has length={len(seq)})'
assert start_move <= final_move, f'start_move {start_move} > final_move {final_move}'
is_subsequence = start_move != 0 or final_move != len(seq) - 1

# set up a log file
import logging
import simple_logger
import traceback # Log excepted exceptions
try:
    import posconstants as pc
except:
    import sys
    path_to_petal = '../petal'
    sys.path.append(os.path.abspath(path_to_petal))
    print('Couldn\'t find posconstants the usual way, resorting to sys.path.append')
    import posconstants as pc
log_dir = pc.dirs['sequence_logs']
log_timestamp = pc.filename_timestamp_str()
log_name = log_timestamp + '_run_sequence.log'
log_path = os.path.join(log_dir, log_name)
logger, logger_fh, logger_sh = simple_logger.start_logger(log_path)
logger_fh.setLevel(logging.DEBUG)
if uargs.verbose:
    logger_sh.setLevel(logging.DEBUG)
else:
    logger_sh.setLevel(logging.INFO)
logger.info(f'Running {script_name} to perform a positioner move + measure sequence.')
logger.info(f'Input file: {uargs.infile}')
subseq_str = f'\n\nA subset of this sequence is to be performed:\n start move_idx = {start_move:3}\n final move_idx = {final_move:3}' if is_subsequence else ''
logger.debug(f'Complete contents:\n\n{seq}{subseq_str}')
logger.info(f'Contents:\n\n{seq.str(max_lines=50)}{subseq_str}')
assert2 = simple_logger.assert2
input2 = simple_logger.input2

# other imports
import sys
import time

def quit_query(question):
    response = input2(f'\n{question} (y/n) >> ')
    if 'n' in response.lower():
        logger.info('User rejected the sequence prior to running. Now quitting.')
        sys.exit(0)

quit_query('Does the sequence look correct?')
logger.info(f'Number of correction submoves: {uargs.num_corr}')
logger.info(f'Number of fvc images per measurement: {uargs.num_meas}')
logger.info(f'Minimum move cycle time: {uargs.cycle_time} sec')

# set up PECS (online control system)
try:
    from pecs import PECS
    kwargs = {'interactive': True, 'logger': logger, 'inputfunc': input2}
    specified_posids = seq.get_posids()
    if any(specified_posids):
        kwargs['posids'] = specified_posids
    pecs = PECS(**kwargs)
    logger.info(f'PECS initialized, discovered PC ids {pecs.pcids}')
    pecs_on = True
    _get_posids = lambda: list(pecs.get_enabled_posids('sub', include_posinfo=False))
    get_all_enabled_posids = lambda: list(pecs.get_enabled_posids('all', include_posinfo=False))
except:
    # still a useful case, for testing some portion of the script offline
    logger.warning('PECS initialization failed (hint: double-check whether you need to join instance in this terminal')
    pecs_on = False
    specified_posids = list(seq.get_posids())
    max_posids_for_debug = 10
    if len(specified_posids) == 0:
        specified_posids = [f'DUMMY{i:05d}' for i in range(max_posids_for_debug)]
    if len(specified_posids) > max_posids_for_debug:
        specified_posids = specified_posids[:max_posids_for_debug]
    _get_posids = lambda: specified_posids

class NoPosidsError(Exception):
    pass

def get_posids():
    '''Wrapper function to get list of currently enabled + selected posids, including
    raising an exception if that list is empty.'''
    posids = _get_posids()
    if len(posids) == 0:
        raise NoPosidsError
    return posids

# check that there are at least some valid posids selected
try:
    initial_selected_posids = get_posids()
    logger.info(f'selected posids: {get_posids()}')
except NoPosidsError:
    logger.error('No positioners were found matching selection(s) for this test.' +
                 ' Suggested things to check: CTRL_ENABLED flags, posfid power,' +
                 ' canbus status, petal ops_state. Now quitting.')
    sys.exit(0)

if pecs_on:
    these = len(get_posids())
    allofthem = len(pecs.get_enabled_posids('all'))
    if these == allofthem:
        quit_query(f'Are you sure you want to be running ALL {allofthem} positioners?')

# helpers for generating move requests
import pandas as pd
import numpy as np
import math

def is_number(x):
    '''Check type to see if it's a common number type.'''
    return isinstance(x, (int, float, np.integer, np.floating))

def is_boolean(x, include01=True):
    '''Check type to see if it's a common boolean type.'''
    if isinstance(x, (bool, np.bool_)):
        return True
    if include01 and is_number(x) and (x == 0 or x == 1):
        return True
    return False

def ptlcall(funcname, posid, *args, **kwargs):
    '''Passes function calls off to the petal instance. Works exclusively for
    functions whose first argument accepts a single positioner id.'''
    role = pecs.ptl_role_lookup(posid)
    extra_kwargs = {'participating_petals': role}
    if kwargs:
        kwargs.update(extra_kwargs)
    else:
        kwargs = extra_kwargs
    function = getattr(pecs.ptlm, funcname)
    result = function(posid, *args, **kwargs)
    if isinstance(result, dict) and len(result) == 1 and list(result.keys())[0] == role:
        # 2020-09-10 [JHS] I don't know how to predict why / when this happens, but some
        # calls through petalman return a single element dict that needs to be decomposed.
        result = list(result.values())[0]
    return result

def trans(posid, method, *args, **kwargs):
    '''Passes calls off to the petal instance to do postransforms conversions.
    The argument "method" is the name of the function in postransforms that
    you want to call.
    '''
    result = ptlcall('postrans', posid, method, *args, **kwargs) 
    return result

# general settings for the move measure function
if uargs.match_radius == None and pecs_on:
    uargs.match_radius = pecs.match_radius
move_meas_settings = {key: uargs.__dict__[key] for key in ['match_radius', 'check_unmatched', 'test_tp', 'anticollision']}
logger.info(f'move_measure() general settings: {move_meas_settings}')

# organize positioner settings
motor_current_keys = {'CURR_SPIN_UP_DOWN', 'CURR_CRUISE', 'CURR_CREEP'}
motor_settings_keys = motor_current_keys | {'CREEP_PERIOD', 'SPINUPDOWN_PERIOD'}
pos_settings = seq.pos_settings.copy()
if uargs.motor_current:
    new = {key: uargs.motor_current for key in motor_current_keys}
    pos_settings.update(new)
    logger.info(f'Motor currents: uniformly set at the run_sequence command line, value={uargs.motor_current}')
logger.info(f'Positioner settings for this test: {pos_settings}')

# caching / retrieval / application of positioner settings
def cache_current_pos_settings(posids):
    '''Gathers current positioner settings and caches them to disk. Only does so
    for those which will actually be changed by pos_settings. Returns a file path.
    '''
    cache_name = f'{log_timestamp}_pos_settings_cache.ecsv'
    cache_path = os.path.join(log_dir, cache_name)
    current = {}
    for key in pos_settings:
        current[key] = pecs.quick_query(key=key, posids=posids, mode='iterable')  # quick_query returns a dict with keys=posids
    reshaped = {key:[] for key in set(current) | {'POS_ID'}}
    for posid in posids:
        reshaped['POS_ID'].append(posid)
        for key in current:
            reshaped[key].append(current[key][posid])
    frame = pd.DataFrame(data=reshaped)
    frame.to_csv(cache_path, index=False)
    return cache_path

def retrieve_cached_pos_settings(path):
    '''Retrieves positioner settings from file cached to disk. Returns a dict with
    keys = posids and values = settings subdicts. This can be directly used in the
    apply_pos_settings function.'''
    frame = pd.read_csv(path)
    settings = {}
    keys = frame.columns.tolist()
    for idx, row in frame.iterrows():
        posid = row['POS_ID']
        these_settings = {key: row[key] for key in keys if key != 'POS_ID'}
        settings[posid] = these_settings
    return settings

def apply_pos_settings(settings):
    '''Push new settings values to the positioners in posids. The dict settings
    should have keys = posids, and values = sub-dictionaries. Each sub-dict
    should have keys / value types like some subset of pos_defaults in the sequence
    module.
    
    Nothing will be committed to the posmoves database (settings will be stored
    only in memory, and can be reset by reinitializing the petal software.
    '''
    if not settings:
        logger.info('apply_pos_settings: none to apply')
        return
    motor_update_petals = set()
    for posid, these_settings in settings.items():
        role = pecs.ptl_role_lookup(posid)
        for key, value in these_settings.items():
            assert2(key in sequence.pos_defaults, f'unexpected pos settings key {key} for posid {posid}')
            default = sequence.pos_defaults[key]
            if is_boolean(default, include01=False):
                test = is_boolean(value, include01=True)
            elif is_number(default):
                test = is_number(value)
            else:
                test = isinstance(value, type(default))
            assert2(test, f'unexpected type {type(value)} for value {value} for posid {posid}')
    accepted_by_petal = pecs.ptlm.batch_set_posfid_val(settings, check_existing=True)
    for accepted in accepted_by_petal.values():
        for posid, these_accepted in accepted.items():
            for key, val_accepted in these_accepted.items():
                value = settings[posid][key]
                if val_accepted == False:  # val_accepted == None is in fact is ok --- just means no change needed
                    assert2(False, f'unable to set {key}={value} for {posid}')
                elif val_accepted == True:  # again, distinct from the None case
                    if key in motor_settings_keys:
                        motor_update_petals.add(role)
    logger.info('apply_pos_settings: Positioner settings updated in memory')
    if motor_update_petals:
        logger.info('apply_pos_settings: Positioner settings include change(s) to motor' +
                    ' parameters. Now pushing new settings to petal controllers on petals' +
                    f' {motor_update_petals}')
        pecs.ptlm.set_motor_parameters(participating_petals=list(motor_update_petals))

def calc_poslocXY_errors(requests, results):
    '''Calculates positioning errors based on initial requested targets and
    FVC measured results. Error values are calculated and returned in poslocXY
    coordinates. Only works for initial requests that were given in absolute
    coordinates.
    
    INPUT:   requests ... dataframe as generated by make_requests()
             results ... dataframe as generated by pecs.move_measure()
                         assumed to contain 'obsX', 'obsY', 'posintT', 'posintP', and index 'DEVICE_ID'
             
    OUTPUT:  tuple of two dicts, both having keys = posids and values = [err_locX, err_locY]
                targ_errs ... calculated w.r.t. requested target
                trac_errs ... calculated w.r.t. posintTP tracking of positioner's location
    '''
    if results.index.name == 'DEVICE_ID':
        results = results.reset_index()
    useful_results = results[['DEVICE_ID', 'obsX', 'obsY', 'posintT', 'posintP']]
    combo = requests.merge(useful_results, on='DEVICE_ID')
    targ_errs = {}
    trac_errs = {}
    df = combo.set_index('DEVICE_ID')
    ret = pecs.ptlm.batch_transform(df, 'obsXY', 'poslocXY')
    df = pd.concat([x for x in ret.values()]).rename(columns={'poslocX': 'x_meas', 'poslocY': 'y_meas'})
    ret = pecs.ptlm.batch_transform(df, 'posintTP', 'poslocXY')
    df = pd.concat([x for x in ret.values()]).rename(columns={'poslocX': 'x_trac', 'poslocY': 'y_trac'})
    combo = df.reset_index()
    for row in combo.iterrows():
        data = row[1]
        posid = data['DEVICE_ID']
        xy_meas = (data['x_meas'], data['y_meas'])
        xy_trac = (data['x_trac'], data['y_trac'])
        command = data['COMMAND']
        uv_targ = [data['X1'], data['X2']]
        assert2(command in sequence.abs_commands, 'error calc when initial request was a delta ' +
                                                  f'move ({command}) is not currently supported')
        if command == 'poslocXY':
            xy_targ = uv_targ
        else:
            conversion = f'{command}_to_poslocXY'
            # Need to think how best to do in batch
            xy_targ = trans(posid, conversion, uv_targ)
        targ_err_x = xy_meas[0] - xy_targ[0]
        targ_err_y = xy_meas[1] - xy_targ[1]
        trac_err_x = xy_meas[0] - xy_trac[0]
        trac_err_y = xy_meas[1] - xy_trac[1]
        targ_errs[posid] = [targ_err_x, targ_err_y]
        trac_errs[posid] = [trac_err_x, trac_err_y]
    return targ_errs, trac_errs       

def summarize_errors(errs, prefix=''):
    '''Produces a string summarizing errors for a set of positioners.
    
    INPUT:   errs ... same format as output by calc_poslocXY_errors, a dict with
                      keys = posids and values = [err_locX, err_locY]
             prefix ... str, will be put at head of output
             
    OUTPUT:  string
    '''
    posids = sorted(errs.keys())
    vec_errs = [math.hypot(errs[posid][0], errs[posid][1]) for posid in posids]
    vec_errs_um = [err*1000 for err in vec_errs]
    first = True
    err_str = f'{prefix} ... '
    tab = ' ' * len(err_str)
    newline = f'\n{tab}'
    for name, func in {'max': max, 'min': min, 'mean': np.mean,
                       'median': np.median, 'std': np.std,
                       'rms': lambda X: (sum(x**2 for x in X)/len(X))**0.5,
                       }.items():
        if not first:
            err_str += ', '
        err_str += f'{name}: {func(vec_errs_um):5.1f}'
        first = False
    sorted_err_idxs = np.argsort(vec_errs_um)
    for desc, count in {'best': uargs.num_best, 'worst': uargs.num_worst}.items():
        n = len(vec_errs_um)
        if n == 0:
            break
        count_mod = min(count, int(n/2)) if n > 1 else 1
        err_str += f'{newline}{desc:<5} {count_mod}:'
        this_sort = list(sorted_err_idxs) if desc == 'best' else list(reversed(sorted_err_idxs))
        these_idxs = [this_sort[i] for i in range(count_mod)]
        for i in these_idxs:
            err_str += f' {posids[i]}: {vec_errs_um[i]:5.1f},'
        err_str = err_str[:-1]  # remove dangling comma
    return err_str

def pause_between_moves(last_move_time, cycle_time):
    '''Pauses the sequence between moves. User is given the option to safely
    CTRL-C exit the sequence during this period.
    
    INPUT:  last_move_time ... seconds since epoch, defining when the previous
                               move was done. last_move_time=None will wait the
                               full cycle_time
            cycle_time ... seconds to wait after last_move_time
    
    OUTPUT: Returns the time now, if the pause was succesfully completed.
            Returns KeyboardInterrupt if user did CTRL-C.
    '''
    if last_move_time == None:
        last_move_time = time.time()
    sec_since_last_move = time.time() - last_move_time
    need_to_wait = cycle_time - sec_since_last_move
    if need_to_wait > 0:
        logger.info(f'Pausing {need_to_wait:.1f} sec for positioner cool down. ' +
                    'It is safe to CTRL-C abort the test during this wait period.\n')
        try:
            dt = 0.2  # sec
            for i in range(int(np.ceil(need_to_wait/dt))):
                time.sleep(dt)
        except KeyboardInterrupt:
            return KeyboardInterrupt
    return time.time()

def get_parkable_neighbors(posids):
    '''Return set of neighbors of posids which can be "parked".
    '''
    parkable_neighbors = set()
    all_enabled = get_all_enabled_posids()
    ret = pecs.ptlm.batch_get_neighbors(posids)
    all_neighbors = {}
    for neigh in ret.values():
        all_neighbors.update(neigh)
    for posid in posids:
        neighbors = all_neighbors[posid]
        enabled_neighbors = set(neighbors) & set(all_enabled)
        parkable_neighbors |= enabled_neighbors
    
    # 2020-09-10 [JHS] Unfortunately, PECS seems not quite set up to move neighbors
    # unless they were also in the selected set of commanded positioners at the
    # beginning. So for now I am simply removing such neighbors --- rendering this
    # function right now a functionally useless placeholder.
    all_known_as_commandable_to_pecs = set(pecs.posids)
    parkable_neighbors &= all_known_as_commandable_to_pecs
    return parkable_neighbors

# setup prior to running sequence
if pecs_on:
    # cache the existing pos settings
    cache_posids = get_posids()
    cache_path = cache_current_pos_settings(cache_posids)
    logger.info(f'Initial settings of positioner(s) cached to: {cache_path}')
    
    # apply new pos settings
    new_settings = {posid: pos_settings.copy() for posid in cache_posids}
    apply_pos_settings(new_settings)
    
    # set phi limit angle for the test
    old_phi_limits = pecs.ptlm.get_phi_limit_angle()
    if uargs.enable_phi_limit:
        some_petal = pecs.ptl_roles[0]
        typical_phi_limit = pecs.ptlm.app_get('typical_phi_limit_angle', participating_petals=some_petal)
        uniform_phi_limit = typical_phi_limit
    else:
        uniform_phi_limit = None
    pecs.ptlm.set_phi_limit_angle(uniform_phi_limit)
    new_phi_limits = pecs.ptlm.get_phi_limit_angle()
    if new_phi_limits == old_phi_limits:
        logger.info(f'Phi limits unchanged: {old_phi_limits}')
    else:
        logger.info(f'Phi limits changed. Old phi limits: {old_phi_limits}, ' +
                    f'new phi limits: {new_phi_limits}')
real_moves = pecs_on and not uargs.no_movement
cooldown_margin = 1.0  # seconds
cycle_time = uargs.cycle_time + cooldown_margin if real_moves else 2.0
last_move_time = 'no moves done yet'

# settings info to store in LOG_NOTE fields
posids = get_posids()
lognote_settings_noms = sequence.pos_defaults.copy()
lognote_settings_noms.update({'GEAR_CALIB_T': 1.0, 'GEAR_CALIB_P': 1.0})
lognote_settings_strs = {posid: '' for posid in posids}
lognote_settings_strs_exist = False
for key, nominal in lognote_settings_noms.items():
    if pecs_on:
        current = pecs.quick_query(key=key, posids=posids, mode='iterable')  # quick_query returns a dict with keys=posids
    else:
        current = {posid: pos_settings[key] if key in pos_settings else nominal for posid in posids}        
        current[posids[0]] = nominal / 2  # just a dummy non-nominal value
    for posid, value in current.items():
        if value != nominal:
            lognote_settings_strs[posid] = pc.join_notes(lognote_settings_strs[posid], f'{key}={value}')
            lognote_settings_strs_exist = True
if lognote_settings_strs_exist:
    logger.info('The following positioners have \'special\' settings during this test:')
    for posid, note in lognote_settings_strs.items():
        if note:
            logger.info(f'{posid}: {note}')

def do_pause():
    '''Convenience wrapper for pause_between_moves(), utilizing a global to
    track the timing.'''
    global last_move_time
    if last_move_time == 'no moves done yet':
        last_move_time = time.time()
    else:
        last_move_time = pause_between_moves(last_move_time, cycle_time)
    if last_move_time == KeyboardInterrupt:
        logger.info('Keyboard interrupt detected.')
        raise StopIteration

parking_overrides = {'FINAL_CREEP_ON': True, 'ANTIBACKLASH_ON': True}      
  
def park(park_option, is_prepark=True):
    '''Perform (or skip if appropriate) parking based on argued park_option.'''
    global last_move_time
    assert park_option in park_options
    if not park_option:
        return
    posids = get_posids()
    if pecs_on:
        neighbors = get_parkable_neighbors(posids)
    else:
        neighbors = set()
    extras = set(neighbors) - set(posids)
    neighbor_text = '' if not extras else f' specified for this test, as well as {len(extras)} of their unspecificed neighbors: {sorted(extras)}'
    logger.info(f'Doing {"initial" if is_prepark else "final"} parking move (coords={park_option}) on {len(posids)} positioners' + neighbor_text)
    all_to_park = sorted(set(posids) | set(neighbors))
    extra_note = sequence.sequence_note_prefix + seq.normalized_short_name
    if is_prepark:
        last_move_time = time.time()
    logger.debug(f'Requested positioners: {sorted(set(all_to_park))}')
    original_settings = {}
    override_settings = {}
    for key in parking_overrides:
        overridden_keys = set()
        if pecs_on:
            originals = pecs.quick_query(key=key, posids=all_to_park, mode='iterable') # quick_query returns a dict with keys=posids
        else:
            originals = {posid: pos_settings[key] if key in pos_settings else sequence.pos_defaults[key] for posid in all_to_park}
        for posid, original in originals.items():
            if original != parking_overrides[key]:
                if posid not in override_settings:
                    override_settings[posid] = {}
                    original_settings[posid] = {}
                original_settings[posid][key] = original
                override_settings[posid][key] = parking_overrides[key]
                overridden_keys.add(key)
    for key in overridden_keys:
        logger.info(f'During parking, {key} will be temporarily overridden with value: {parking_overrides[key]}')
    if real_moves:
        kwargs = {'posids': all_to_park,
                  'mode': 'normal',
                  'coords': park_option,
                  'log_note': extra_note,
                  }
        kwargs.update(move_meas_settings)
        if 'anticollision' in kwargs:
            del kwargs['anticollision']  # not an arg to park_and_measure
        if uargs.test_tp:
            orig_tp_frac = pecs.tp_frac
            pecs.tp_frac = 1.0
            if orig_tp_frac != pecs.tp_frac:
                logger.info(f'For pre-parking, temporarily adjusted tp update error fraction from {orig_tp_frac} to {pecs.tp_frac}')
        apply_pos_settings(override_settings)
        if cycle_fids:
            logger.info('Turning on fiducials')
            pecs.ptlm.set_fiducials('on', participating_petals=pecs.ptlm.get('fid_petals'))
        pecs.park_and_measure(**kwargs)
        if cycle_fids:
            logger.info('Turning off fiducials')
            pecs.ptlm.set_fiducials('off', participating_petals=pecs.ptlm.get('fid_petals'))
        apply_pos_settings(original_settings)
        if uargs.test_tp and orig_tp_frac != pecs.tp_frac:
            pecs.tp_frac = orig_tp_frac
            logger.info(f'Restored tp update error fraction to {pecs.tp_frac}')
    if any(seq) and is_prepark:
        do_pause()
    logger.info('Parking complete\n')

# do the sequence
logger.info('Beginning the move sequence\n')
exception_here = None  # storage for exceptions
try:
    if uargs.prepark:
        park(park_option=uargs.prepark, is_prepark=True)
    move_counter = 0
    num_moves = final_move - start_move + 1
    for m in range(start_move, final_move + 1):
        move = seq[m]
        posids = get_posids()  # dynamically retrieved, in case some positioner gets disabled mid-sequence
        move_counter += 1
        move_num_text = f'target {move_counter} of {num_moves} (sequence_move_idx = {m})'
        if not move.is_defined_for_all_positioners(posids):
            logger.warning(f'Skipping {move_num_text}, because targets not defined for some positioners.\n')
            continue
        logger.info(f'Preparing {move_num_text} on {len(posids)} positioner{"s" if len(posids) > 1 else ""}.')
        correction_not_defined = move.command not in sequence.general_commands
        correction_not_allowed = not move.allow_corr
        not_correctable = correction_not_allowed or correction_not_defined 
        if uargs.num_corr > 0 and not_correctable:
            reason = f'not defined for {move.command} moves' if correction_not_defined else ''
            reason += ', and ' if correction_not_defined and correction_not_allowed else ''
            reason += 'disabled for this row in the sequence file' if correction_not_allowed else ''
            logger.info(f'Correction move skipped ({reason})')
        n_corr = 0 if not_correctable else uargs.num_corr
        targ_errs = None
        calc_errors = move.command not in sequence.delta_commands | sequence.homing_commands
        initial_request = None
        if (n_corr > 0) and cycle_fids and pecs_on:
            logger.info('Turing on fiducials')
            pecs.ptlm.set_fiducials('on', participating_petals=pecs.ptlm.get('fid_petals'))
        for submove_num in range(1 + n_corr):
            extra_log_note = pc.join_notes(f'sequence_move_idx {m}', f'move {move_counter}')
            submove_txt = f'submove {submove_num}'
            move_with_submove_txt = f'{move_num_text}, {submove_txt}'
            if n_corr > 0:
                extra_log_note = pc.join_notes(extra_log_note, submove_txt)
            if move.command in sequence.general_commands:
                if pecs_on:
                    move_measure_func = pecs.move_measure
                if submove_num == 0:
                    submove = move
                else:
                    operable = get_posids()  # dynamically retrieved, in case some positioner gets disabled mid-sequence
                    posids = [p for p in operable if p in targ_errs]  # exclude any pos that had no err result in previous submove (e.g. unmatched case)
                    no_err_val = set(operable) - set(posids)
                    if any(no_err_val):
                        logger.info(f'{len(no_err_val)} positioners excluded from next submove, due to missing error results.')
                        logger.debug(f'Excluded posids: {no_err_val}')
                    # note below how order is preserved for target0 and target1
                    # lists, on the assumption that get_posids() returns a list
                    submove = sequence.Move(command='poslocdXdY',
                                            target0=[-targ_errs[posid][0] for posid in posids],
                                            target1=[-targ_errs[posid][1] for posid in posids],
                                            posids=posids,
                                            log_note=move.get_log_notes(posids),
                                            allow_corr=move.allow_corr)                    
                request = submove.make_request(posids=posids, log_note=extra_log_note)
                these_notes = request[['DEVICE_ID', 'LOG_NOTE']].to_dict(orient='list')
                for i, posid in enumerate(these_notes['DEVICE_ID']):
                    these_notes['LOG_NOTE'][i] = pc.join_notes(these_notes['LOG_NOTE'][i], lognote_settings_strs[posid])
                request.update(these_notes)
                if submove_num == 0:
                    initial_request = request
                if pecs_on:
                    request = request.merge(pecs.posinfo, on='DEVICE_ID')
                kwargs = {'request': request, 'num_meas': uargs.num_meas}
                descriptive_dict = submove.to_dict(sparse=True, posids=posids)
            elif move.command in sequence.homing_commands:
                if pecs_on:
                    move_measure_func = pecs.rehome_and_measure
                kwargs = move.make_homing_kwargs(posids=posids, log_note=extra_log_note)
                submove = move  # just for consistency of operations below
            else:
                logger.warning(f'Skipping {move_with_submove_txt} due to unexpected command {move.command}\n')
                continue
            kwargs.update(move_meas_settings)
            descriptive_dict = submove.to_dict(sparse=True, posids=posids)
            logger.debug(f'Move command: {descriptive_dict}')
            if submove.is_uniform:
                logger.info(f'Num requested positioners: {len(posids)}')
                logger.debug(f'Requested positioners: {posids}')  # non-uniform cases already print a bunch of posids when displaying move command
            logger.info(f'Going to {move_with_submove_txt}')
            if real_moves:
                results = move_measure_func(**kwargs)
                prev_posids = posids
                posids = get_posids()  # dynamically retrieved, in case some positioner gets disabled mid-sequence
                if len(prev_posids) != len(posids):
                    logger.info(f'Num posids changed from {len(prev_posids)} to {len(posids)}. ' +
                                'Check for auto-disabling after comm error. Lost positioners: ' +
                                f'{set(prev_posids) - set(posids)}')
            if calc_errors:
                if real_moves:
                    targ_errs, trac_errs = calc_poslocXY_errors(initial_request, results)
                else:
                    dummy_err_x1 = np.random.normal(loc=0, scale=0.1, size=len(posids))
                    dummy_err_y1 = np.random.normal(loc=0, scale=0.1, size=len(posids))
                    targ_errs = {posids[i]: [dummy_err_x1[i], dummy_err_y1[i]] for i in range(len(posids))}
                    dummy_err_x2 = np.random.normal(loc=0, scale=0.1, size=len(posids))
                    dummy_err_y2 = np.random.normal(loc=0, scale=0.1, size=len(posids))
                    trac_errs = {posids[i]: [dummy_err_x2[i], dummy_err_y2[i]] for i in range(len(posids))}
                err_str = f'Results for {move_with_submove_txt}, n_pos={len(posids)}, errors given in um:'
                err_str += '\n' + summarize_errors(targ_errs, prefix='TARGETING')
                err_str += '\n' + summarize_errors(trac_errs, prefix=' TRACKING')
                logger.info(err_str + '\n')
            more_corrections_to_do = submove_num < n_corr
            if not(more_corrections_to_do) and cycle_fids and pecs_on:
                logger.info('Turning off fiducials')
                pecs.ptlm.set_fiducials('off', participating_petals=pecs.ptlm.get('fid_petals'))
            more_moves_to_do = more_corrections_to_do or m < final_move or uargs.postpark
            if not(correction_pause):
                last_move_time = time.time()
            if more_moves_to_do and (not(more_corrections_to_do) or correction_pause):
                do_pause()
    if uargs.postpark:
        park(park_option=uargs.postpark, is_prepark=False)
except StopIteration:
    logger.info('Safely aborting the sequence.')
except NoPosidsError:
    logger.error('All positioners have become disabled. Aborting the sequence.')
except Exception as e:
    exception_here = e
    logger.error('The sequence crashed! See traceback below:')
    logger.critical(traceback.format_exc())
    logger.info('Attempting to preform cleanup before hard crashing. Configure the instance before trying again.')
if exception_here is None:
    logger.info(f'Sequence "{seq.short_name}" complete!')

# cleanup after running sequence
if pecs_on:
    # Trigger FVCCollector to backup FVC images, should not wait for collection to finish
    pecs.fvc_collect()
    # restore the original pos settings
    orig_settings = retrieve_cached_pos_settings(cache_path)
    logger.info(f'Retrieved original positioner settings from {cache_path}')
    new_cache_path = cache_current_pos_settings(cache_posids)
    new_settings = retrieve_cached_pos_settings(new_cache_path)
    if orig_settings == new_settings:
        logger.info('No net change of positioner settings detected between start and finish' +
                    ' of sequence. No further modifications to postioner state required.')
    else:
        logger.info('Some net change(s) of positioner settings detected between start and' +
                    ' finish of sequence. Now restoring system to original state.')
        apply_pos_settings(orig_settings)
        logger.info('Restoration of original positioner settings complete.')
    
    # restore old phi limit angles
    if new_phi_limits != old_phi_limits:
        if isinstance(old_phi_limits, dict):
            # 2020-10-27 [KF] old_phi_limits will be a dict with rolename keys (PETALx)
            # if the petals had different limits. Otherwise will be a single value.
            roles_angles = [(role, angle) for role, angle in old_phi_limits.items()]
        else:
            roles_angles = [(None, old_phi_limits)]
        for ra in roles_angles:
            pecs.ptlm.set_phi_limit_angle(ra[1], participating_petals=ra[0])
        restored_phi_limits = pecs.ptlm.get_phi_limit_angle()
        if restored_phi_limits == old_phi_limits:
            logger.info(f'Phi limits restored to: {restored_phi_limits}')
        else:
            logger.warning('Some error when restoring phi limits. Old limits were' +
                           f' {old_phi_limits} but restored values are different:' +
                           f' {restored_phi_limits}')

# final thoughts...
logger.info(f'Log file: {log_path}')
simple_logger.clear_logger()

# raise exception from loop if we have one
if exception_here is not None:
    raise(exception_here)
