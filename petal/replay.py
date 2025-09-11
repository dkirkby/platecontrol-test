#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Script to replay in simulation, a sequence of moves done on hardware. Note:
Get data from online DB using desimeter/get_posmoves tool, including options
--with-calib and --tp-updates.
'''

import os
get_move_id = lambda expid, expiter: f'{expid:05}.{expiter:03}'
null_move_id = get_move_id(0,0)
null_move_id_large = get_move_id(99999,999)

# command line argument parsing
import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-i', '--infiles', type=str, required=True, nargs='*',
                    help='Path to input csv file(s), containing results from "get_posmoves -c -t" ' +
                         '(see that function''s help for full syntax). Regex is ok (like M*.csv). Multiple ' +
                         'file args are also ok (like M00001.csv M00002.csv M01*.csv), as is a directory ' +
                         'that contains the files.')
parser.add_argument('-a', '--anticollision', type=str, default='adjust', help='anticollision mode, can be "adjust", "adjust_requested_only", "freeze" or None. Default is "adjust"')
parser.add_argument('-p', '--enable_phi_limit', action='store_true', help='turns on minimum phi limit for move targets, default is False')
parser.add_argument('-ms', '--start_move', type=str, default=null_move_id, help='syntax: <exposure_id>.<exposure_iter>, start the simulation at this move (or defaults to first move in data)')
parser.add_argument('-mf', '--final_move', type=str, default=null_move_id_large, help='syntax: <exposure_id>.<exposure_iter>, finish the simulation at this move (or defaults to last move in data)')
parser.add_argument('-anim', '--animate', action='store_true', help='plot an animation of simulated moves (can be slow), defaults to False')
parser.add_argument('-f', '--focus', type=str, default=None, help='focus the animation in on a particular positioner and its neighbors. Identify it either by device location integer or POS_ID')
parser.add_argument('-x', '--focus_expand', action='store_true', help='when focus option is specified, this causes not just neighbors, but also neighbors-of-neighbors to be animated')
parser.add_argument('-v', '--verbose', action='store_true', help='turn on additional verbosity at console, may be helpful for debugging')
parser.add_argument('-t', '--tp_update', action='store_true', help='for each move, force POS_T and POS_P to adopt values from data file (rather than using the simulated values from the previous move), thus any tp_updates done by the online system will be incorporated in sim')
parser.add_argument('-d', '--display_posids', type=str, default='', help='one or more posids (comma-separated) for which you want extra information printed')
parser.add_argument('-c', '--reuse_cached_data', action='store_true', help='re-use cached data from previous run (faster to read from disk for repeat uses)')
parser.add_argument('-r', '--reconstruct_corrections', action='store_true', help='attempts to reconstruct correction delta targets (useful for debugging a crash during a correction move, before data got into DB')

uargs = parser.parse_args()
if uargs.anticollision == 'None':
    uargs.anticollision = None
assert uargs.anticollision in {'adjust', 'adjust_requested_only', 'freeze', None}, f'bad argument {uargs.anticollision} for anticollision parameter'

import time
import glob
import petal
import posconstants as pc
import posstate
from astropy.table import Table, vstack, MaskedColumn
import numpy as np

# data types / definitions
nulls = {'--', None, 'none', 'None', '', 0, False, 'False', 'FALSE', 'false'}
boolean = lambda x: x not in nulls
required = {'DATE': str, 'EXPOSURE_ID': int, 'EXPOSURE_ITER': int, 'POS_MOVE_INDEX': int,
            'POS_ID': str, 'PETAL_ID': int, 'DEVICE_LOC': int,
            'MOVE_CMD': str, 'MOVE_VAL1': str, 'MOVE_VAL2': str, 'LOG_NOTE': str,
            'POS_T': float, 'POS_P': float, 'PTL_X': float, 'PTL_Y': float,
            }
param_keys = {'LENGTH_R1': float, 'LENGTH_R2': float, 'OFFSET_X': float, 'OFFSET_Y': float,
              'OFFSET_T': float, 'OFFSET_P': float, 'PHYSICAL_RANGE_T': float, 'PHYSICAL_RANGE_P': float,
              'KEEPOUT_EXPANSION_PHI_RADIAL': float, 'KEEPOUT_EXPANSION_PHI_ANGULAR': float,
              'KEEPOUT_EXPANSION_THETA_RADIAL': float, 'KEEPOUT_EXPANSION_THETA_ANGULAR': float,
              'CLASSIFIED_AS_RETRACTED': boolean, 'CTRL_ENABLED': boolean,
             }
required.update(param_keys)
mask_fills = {}
for key, func in required.items():
    if func == float:
        mask_fills[key] = 0.0
    elif func == int:
        mask_fills[key] = 0
    elif func == boolean:
        mask_fills[key] = False
    elif func == str:
        mask_fills[key] = ''

# read in the move data
infiles = []
cache_basename = 'replay_cache.ecsv'
cache_path = os.path.join(pc.dirs['temp_files'], cache_basename)
timer_start = time.perf_counter()
print_read_complete = lambda: print(f'Read complete in {time.perf_counter() - timer_start:.2f} seconds.')
if uargs.reuse_cached_data:
    print(f'Reading data from cache at {cache_path}...')
    t = Table.read(cache_path)
    print_read_complete()
    
    # handling astropy masked columns
    fill_values = {key: '' for key in ['MOVE_CMD', 'MOVE_VAL1', 'MOVE_VAL2']}
    for key, fill in fill_values.items():
        t[key].fill_value = fill
        t[key] = t[key].filled()
else:
    print('Reading and parsing data from argued infiles...')
    for s in uargs.infiles:
        these = glob.glob(s)
        for this in these:
            if os.path.isdir(this):
                contents = os.listdir(this)
                contents = [os.path.join(this, file) for file in contents]
                infiles.extend(contents)
            else:
                infiles.append(this)
    infiles = [os.path.realpath(p) for p in infiles]
    assert any(infiles), 'check path args, no valid input files detected'
    print(f'{len(infiles)} infiles detected. Now reading data...')
    input_tables = []
    for path in infiles:
        if 'csv' not in path:
            continue  # i.e. skip non-data files
        table = Table.read(path)
        missing = set(required) - set(table.columns)
        if any(missing):
            print(f'Skipped data at {path} due to missing columns {missing}')
            continue
        new = {}
        for key, typefunc in required.items():
            if isinstance(table[key], MaskedColumn):
                vec = []
                for i in range(len(table)):
                    if table[key].mask[i]:
                        x = mask_fills[key]
                    else:
                        x = typefunc(table[key][i])
                    vec.append(x)
            else:
                vec = [typefunc(x) for x in table[key]]
            new[key] = vec
        new_table = Table(new)
        if len(new_table) > 0:
            input_tables.append(Table(new))
    t = vstack(input_tables)
    t.sort(keys=['POS_ID', 'DATE'])
    print_read_complete()
    t.write(cache_path, overwrite=True)
    print(f'Wrote data to cache at {cache_path}. (For repeat runs with same data, you can use -c option to read this data from cache, which is much faster.)')

print('Setting up petal data...')

# confirm only one petal
petal_ids = set(t['PETAL_ID'])
assert len(petal_ids) == 1, f'replay script does not support {len(petal_ids)} petals'
petal_id = petal_ids.pop()

# identify moving positioners and their neighbors
no_cmd = [x in nulls for x in t['MOVE_CMD']]
has_cmd = [not(x) for x in no_cmd]
t['HAS_MOVE_CMD'] = has_cmd
moving = set(t['POS_ID'][has_cmd])
nonmoving = set(t['POS_ID']) - moving
assert any(moving), 'no moving positioners found in data'
nonmoving_neighbors = set()
known_device_locs = set(t['DEVICE_LOC'])
for posid in moving:
    this_loc = t['DEVICE_LOC'][t['POS_ID'] == posid][0]
    neighbor_locs = pc.generic_pos_neighbor_locs[this_loc]
    for loc in neighbor_locs:
        if loc not in known_device_locs:
            continue
        neighbor_posid = t['POS_ID'][t['DEVICE_LOC'] == loc][0]
        if neighbor_posid in nonmoving:
            nonmoving_neighbors.add(neighbor_posid)
all_posids = moving | nonmoving_neighbors

# organize the display posids
display_posids = set(uargs.display_posids.split(','))
display_posids -= {''}
invalid_display_posids = display_posids - all_posids
assert len(invalid_display_posids) == 0, f'requested display posid(s) were not found in data: {invalid_display_posids}'

# identify / define which rows should be run in simulation
t['MOVE_ID'] = [get_move_id(t["EXPOSURE_ID"][i], t["EXPOSURE_ITER"][i]) for i in range(len(t))]
has_cmd_idxs = np.where(t['HAS_MOVE_CMD'])[0].tolist()
first_move_found = t[has_cmd_idxs[0]]['MOVE_ID']
last_move_found = t[has_cmd_idxs[-1]]['MOVE_ID']
def check_move_id(val):
    try:
        split = val.split('.')
        assert len(split) == 2
        move_id = get_move_id(int(split[0]), int(split[1]))
        return move_id
    except:
        assert False, f'invalid move id argument {val}'
start_move = check_move_id(uargs.start_move)
final_move = check_move_id(uargs.final_move)
start_move = max(start_move, first_move_found)
final_move = min(final_move, last_move_found)
assert start_move <= final_move, f'final move index {final_move} must be >= than start {start_move}'
should_run = np.logical_and(t['MOVE_ID'] >= start_move, t['MOVE_ID'] <= final_move)
should_run = np.logical_and(should_run, t['HAS_MOVE_CMD'])
t['SHOULD_RUN'] = should_run

# set up columns for reconstructed correction calcs
t['RECORR_VAL1'] = [0.0]*len(t)
t['RECORR_VAL2'] = [0.0]*len(t)

# helper functions
def get_row(posid, move_id=None, date=None, prior=False):
    '''Get matching row in table for a given posid and either move_id or date.
    Argue prior=True to get the last row just before the match. If no prior rows
    exist, then it will return the matched one.
    '''
    assert (move_id and not date) or (date and not move_id), 'must only argue either move_id or date'
    match_field = 'MOVE_ID' if move_id else 'DATE'
    match_value = move_id if move_id else date
    if prior:
        matches = t[match_field] < match_value
    else:
        matches = t[match_field] == match_value
    matches = np.logical_and(t['POS_ID'] == posid, matches)
    match_idxs = np.where(matches)[0].tolist()
    if prior and not any(matches):
        return get_row(posid, move_id=move_id, date=date, prior=False)
    row = t[match_idxs[-1]]
    return row
        
def set_params(move_id, ptl=None, tp_update_mode=None):
    '''Set state parameters for all positioners at a given move index.
    INPUTS:  move_id ... in the general astropy table "t"
             ptl ... Petal instance
             tp_update_mode ... None --> no update
                                'data' --> set POS_T, POS_P using values from CSV data file
                                'parked' --> use a typical "parked" value for POS_T, POS_P
    '''
    posids = ptl.posids if ptl else all_posids
    keys = set(param_keys) | {'DEVICE_LOC'}
    for posid in posids:
        if ptl:
            state = ptl.states[posid]
        else:
            state = posstate.PosState(unit_id=posid, device_type='pos', petal_id=petal_id)
        row = get_row(posid, move_id=move_id)
        row_dict = {key: row[key] for key in keys}
        if tp_update_mode == 'data':
            position_row = get_row(posid, date=row['DATE'], prior=True)
            row_dict['POS_T'] = position_row['POS_T']
            row_dict['POS_P'] = position_row['POS_P']
        elif tp_update_mode == 'parked':
            row_dict['POS_T'] = 0 
            row_dict['POS_P'] = 150
        for key, value in row_dict.items():
            state.store(key, value, register_if_altered=False)
        state.write()
        
def make_request(move_cmd_str):
    '''Generate a move request, by parsing a MOVE_CMD string.'''
    a = move_cmd_str.split('=')
    command = a[0]
    b = a[1].split(';')
    targ_strs = b[0].strip('[').strip(']').split(',')
    target = [float(x) for x in targ_strs]
    request = {'command':command, 'target':target, 'log_note':''}
    return request

def print_positions(prefix, move_id, posids=display_posids):
    '''Print out current positions for positioners.'''
    if posids:
        print(f'\n{prefix} MOVE_ID {move_id}...')
        print(ptl.quick_table(posids))
        print('\n')

# initialize petal
move_idxs = np.where(t['SHOULD_RUN'])[0]
first_move_id = t['MOVE_ID'][move_idxs[0]]
set_params(first_move_id, ptl=None, tp_update_mode='parked')
print('Initializing petal instance...')
ptl = petal.Petal(petal_id        = petal_id,
                  petal_loc       = 3,
                  posids          = all_posids,
                  fidids          = {},
                  simulator_on    = True,
                  db_commit_on    = False,
                  local_commit_on = False,
                  local_log_on    = False,
                  collider_file   = None,
                  sched_stats_on  = True,
                  anticollision   = uargs.anticollision,
                  verbose         = False,  # typically distinct from uargs.verbose
                  phi_limit_on    = uargs.enable_phi_limit,
                  )

# initialize animator
anim_label_size = 'medium' # size in points, 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'
anim_cropping_on = True # crops the plot window to just contain the animation
if uargs.animate:
    ptl.animator.cropping_on = True
    ptl.animator.label_size = anim_label_size
    if uargs.focus != None:
        try:
            focus_loc = int(uargs.focus)
            focus_posid = ptl.devices[focus_loc]
        except:
            focus_posid = uargs.focus
            focus_loc = ptl.posmodels[focus_posid].deviceloc
        assert focus_posid in all_posids, f'cannot focus video on missing posid {focus_posid}'
        neighbor_locs = pc.generic_pos_neighbor_locs[focus_loc]
        if uargs.focus_expand:  # include neighbors of neighbors
            for n in neighbor_locs.copy():
                neighbor_locs |= pc.generic_pos_neighbor_locs[n]
        neighbor_posids = {ptl.devices[loc] for loc in neighbor_locs if loc in ptl.devices}
        all_foci = set([focus_posid]) | neighbor_posids
        ptl.collider.posids_to_animate = all_foci
        ptl.collider.fixed_items_to_animate = set()
    ptl.start_gathering_frames()

# run the sequence
move_ids = sorted(set(t['MOVE_ID']) - {null_move_id})
for move_id in move_ids:
    run_now = np.logical_and(t['SHOULD_RUN'], move_id == t['MOVE_ID'])
    run_now_idxs = np.where(run_now)[0].tolist()
    if len(run_now_idxs) == 0:
        continue
    rows = t[run_now]
    move_id_str = pc.join_notes(*[f'{key}: {rows[key][0]}' for key in ['MOVE_ID', 'DATE', 'EXPOSURE_ID', 'EXPOSURE_ITER']])
    print('\n\n', move_id_str, '\n')
    tp_update_mode = None
    if move_id == first_move_id or uargs.tp_update:
        tp_update_mode = 'data'
    set_params(move_id, ptl=ptl, tp_update_mode=tp_update_mode)
    if ptl.schedule_stats.is_enabled():
        ptl.schedule_stats.add_note(move_id_str)
    requests = {}
    for row in rows:
        if row['HAS_MOVE_CMD']:
            posid = row['POS_ID']
            requests[posid] = make_request(row['MOVE_CMD'])
    print_positions('BEFORE', move_id)
    ptl.request_targets(requests)
    ptl.schedule_moves()
    for posid in display_posids:
        print(f'\n---- {posid} @ MOVE_ID {move_id} ----')
        row = get_row(posid, move_id=move_id)
        for key in ['MOVE_CMD', 'MOVE_VAL1', 'MOVE_VAL2', 'LOG_NOTE']:
            print(f' {key}: {row[key]}')
        if posid in ptl.schedule.move_tables:
            ptl.schedule.move_tables[posid].for_hardware() # dummy call, to ensure inclusion of final creep rows in next display line
            ptl.schedule.move_tables[posid].display()
            ptl.schedule.move_tables[posid].display_for('schedule')
            ptl.schedule.move_tables[posid].display_for('hardware')
        print('\n')
    failed_posids = ptl.send_and_execute_moves()
    print_positions('AFTER', move_id)
    if uargs.verbose:
        tab = '   '
        print(f'{"POSID":7}{tab}{"COORD":8}{tab}{"SIMULATED":20}{tab}{"FROM_FILE":20}{tab}{"ERROR":>6}')
        coords = {'posintTP': ['POS_T', 'POS_P'], 'ptlXY': ['PTL_X', 'PTL_Y']}
        for posid in set(requests):
            expected = ptl.posmodels[posid].expected_current_position
            posid_str = f'{posid:7}'
            for ptl_coord, data_coord in coords.items():
                exp = expected[ptl_coord]
                dat = list(get_row(posid, move_id=move_id)[data_coord])
                err = np.hypot(exp[0]-dat[0], exp[1]-dat[1])
                exp_str = f'({exp[0]:8.3f}, {exp[1]:8.3f})'
                dat_str = f'({dat[0]:8.3f}, {dat[1]:8.3f})'
                print(f'{posid_str}{tab}{ptl_coord:8}{tab}{exp_str:20}{tab}{dat_str:20}{tab}{err:>6.3f}')
                if posid_str:
                    posid_str = ' ' * len(posid_str)
    if uargs.reconstruct_corrections:
        # 2020-11-16 [JHS] currently only supporting poslocXY, and only the first correction move
        # that is what I need "today" for some debugging, and more generality will take too much time
        recorr_requests = {}
        for i in run_now_idxs:
            posid = t[i]['POS_ID']
            req = requests[posid]
            if req['command'] == 'poslocXY':
                trans = ptl.posmodels[posid].trans
                meas_ptl_xy = [t[i][key] for key in ['PTL_X', 'PTL_Y']]
                meas_loc_xy = trans.ptlXY_to_poslocXY(meas_ptl_xy)
                delta_xy = trans.delta_XY(meas_loc_xy, req['target'])
                t[i]['RECORR_VAL1'] = delta_xy[0]
                t[i]['RECORR_VAL2'] = delta_xy[1]
                recorr_requests[posid] = {'command': 'poslocdXdY', 'target':delta_xy}
        sel = [i for i in run_now_idxs if t[i]['POS_ID'] in display_posids]
        if sel:
            print(f'\nDoing a reconstructed correction move on {len(recorr_requests)} positioner(s)...')
            print('\nCorrection amounts for display posids are:')
            print(t[sel]['POS_ID', 'RECORR_VAL1', 'RECORR_VAL2'])
            print('')
            ptl.request_targets(recorr_requests)
            ptl.schedule_moves()
            ptl._cancel_move()
        
if ptl.schedule_stats.is_enabled():
    ptl.schedule_stats.save(path=ptl.sched_stats_path, footers=True)
    print(f'Stats saved to {ptl.sched_stats_path}')
if uargs.animate and not ptl.animator.is_empty():
    ptl.stop_gathering_frames()
    print('Generating animation (this can be quite slow)...')
    ptl.animator.filename_suffix = pc.filename_timestamp_str()
    ptl.animator.add_timestamp_prefix_to_filename = False
    ptl.generate_animation() 
