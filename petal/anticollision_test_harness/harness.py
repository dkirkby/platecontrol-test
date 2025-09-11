"""Script to run the anticollision test harness.
"""

import os
import sys
sys.path.append(os.path.abspath('../../petal/'))
import petal
import posconstants as pc
import harness_constants as hc
import sequences as harness_sequences
import posstate
import random
from astropy.table import Table
import numpy as np
sys.path.append(os.path.abspath('../../pecs/'))
import sequence as move_sequence

# uniform starting position
start_posintTP = [0.0, 150.0]

# Some pre-cooked device location selections (can be used below)
locations_all = 'all'
locations_near_gfa = {309,328,348,368,389,410,432,454,474,492,508,521,329,349,369,390,411,433,455,475,493,509,522,532,330,350,370,391,412,434,456,476,494,510,523,533,351,371,392,413,435,457,477,495,511,524,372,393,414,436,458,512,525,415,437,459,526}
subset7a = {78, 79, 87, 88, 89, 98, 99}
subset7b = {21, 22, 26, 27, 28, 33, 34}

# Selection of which device location ids to send move requests to
# (i.e. which positioners on petal to directly command)
# either a set of device locations, or the keyword 'all' or 'near_gfa'
device_loc_to_command = subset7b # note pre-cooked options above

# Select devices to CLASSIFY_AS_RETRACTED and disable
retract_and_disable_sets = {0: {},
                            10: {231, 201, 170, 425, 396, 45, 272, 179, 84, 220},
                            25: {258, 387, 520, 137, 269, 24, 28, 286, 288, 165, 170, 305, 307, 57, 443, 319, 193, 322, 67, 228, 240, 113, 243, 116, 123},
                            50: {6, 263, 520, 266, 523, 524, 142, 271, 16, 145, 18, 400, 149, 157, 288, 416, 165, 297, 425, 45, 431, 176, 304, 59, 446, 64, 450, 69, 72, 202, 81, 83, 339, 213, 86, 354, 228, 104, 365, 494, 372, 380, 250, 124, 126},
                            100: {2, 8, 520, 15, 18, 26, 29, 31, 42, 44, 50, 74, 81, 87, 90, 91, 92, 95, 98, 104, 108, 111, 113, 119, 122, 127, 130, 137, 139, 152, 165, 168, 175, 190, 192, 196, 203, 207, 210, 215, 218, 223, 225, 226, 235, 238, 258, 262, 269, 273, 284, 292, 293, 298, 301, 307, 314, 315, 322, 324, 326, 334, 337, 342, 347, 348, 349, 354, 357, 360, 371, 372, 378, 381, 382, 386, 392, 428, 431, 437, 452, 457, 461, 462, 464, 465, 467, 469, 475, 483, 486, 487, 493, 495, 508},
                            }
retract_and_disable = retract_and_disable_sets[0] #{87,88} # enter device locations to simulate those positioners as retracted and disabled
retracted_TP = [0, 110]

# Set any non 1.0 output ratio scales
# format of dict: keys = posids, values = subdicts like {'T': 0.7} or {'P': 0.2, 'T': 0.4} etc
# scale == 0 --> disable that axis entirely
scale_changes = {'M02182': {'T': 0.5},
                 'M01981': {'P': 0.3},
                 'M06389': {'T': 0.0, 'P': 0.4},
                 }
set_scale_changes_as_retracted = True

# Whether to include any untargeted neighbors in the calculations
include_neighbors = True

# Whether to test some "expert" mode commands
test_direct_dTdP = False
test_homing = False  # note that this one will look a bit weird, since there are no hardstops in simulation. So the results take a bit of extra inspection, but still quite useful esp. to check syntax / basic function

# Override for petal simulated hardware failure rates
sim_fail_freq = {'send_tables': 0.0} 

# Selection of which pre-cooked sequences to run. See "sequences.py" for more detail.
runstamp = hc.compact_timestamp()
pos_param_sequence_id = 'ptl01_sept2020_nominal' # 'cmds_unit_test'
move_request_sequence_id = 'ptl01_set00_double' # 'cmds_unit_test'
ignore_params_ctrl_enabled = False # turn on posids regardless of the CTRL_ENABLED column in params file
new_stats_per_loop = True # save a new stats file for each loop of this script

# Other ids and notes
fidids = {}
petal_id = 0
note = ''
filename_suffix = str(runstamp) + '_' + str(move_request_sequence_id) + ('_' + str(note) if note else '')

# Animation on/off options
should_animate = False
anim_label_size = 'medium' # size in points, 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'
anim_cropping_on = True # crops the plot window to just contain the animation

# Animation foci
#   {} or 'all' to animate everything, including PTL and GFA.
#   'colliding' to only animate colliding bots.
#   'commanded' to only animate the robots that receive move requests
#   Otherwise, give a specific set of posids to animate. Can include 'GFA' or 'PTL' as desired.
animation_foci = 'all'

# other options
n_corrections = 0 # number of correction moves to simulate after each target
max_correction_move = 0.1/1.414 # mm
should_profile = False
should_inspect_some_TP = False # some *very* verbose printouts of POS_T, OFFSET_T, etc, sometimes helpful for debugging
expand_keepouts = False

# randomizer for correction moves
randomizer_seed = 0
random.seed(randomizer_seed)

# resolve which device locations to operate on
if device_loc_to_command == 'all':
    device_loc_to_command = set(pc.generic_pos_neighbor_locs.keys())
all_device_loc = set(device_loc_to_command)
if include_neighbors:
    for loc in device_loc_to_command:
        neighbors = pc.generic_pos_neighbor_locs[loc]
        all_device_loc = all_device_loc.union(neighbors)

# saving of target sets for later use on hardware
# formatted for direct input to pecs/run_sequence.py
should_export_targets = True
encode_parking_in_export_targets = False  # 2020-09-08 [JHS] for now, plan to do this with command line arg in run_sequence.py
if should_export_targets:
    hw_seq = move_sequence.Sequence(short_name=f'xytest_{filename_suffix}',
                                    long_name=f'xy test sequence {filename_suffix}, generated by harness.py',
                                    details='')
    if encode_parking_in_export_targets:
        initial_move = move_sequence.Move(command='posintTP',
                                          target0=[start_posintTP[0]] * len(device_loc_to_command),
                                          target1=[start_posintTP[1]] * len(device_loc_to_command),
                                          device_loc=sorted(device_loc_to_command),
                                          log_note=f'{runstamp} initial position',
                                          pos_settings={},
                                          allow_corr=False)
        hw_seq.append(initial_move)

def update_command_format(command_str):
    '''handles older file format'''
    if command_str == 'posXY': 
        return 'poslocXY'
    return command_str

# Run the sequences.
pos_param_sequence = harness_sequences.get_positioner_param_sequence(pos_param_sequence_id, all_device_loc)
move_request_sequence = harness_sequences.get_move_request_sequence(move_request_sequence_id, device_loc_to_command)
expected_result_keys = {('target_posintT', 'target_posintP'): 'posintTP',
                        ('target_poslocT', 'target_poslocP'): 'poslocTP',
                        ('target_poslocX', 'target_poslocY'): 'poslocXY',
                        ('target_ptlX', 'target_ptlY'): 'ptlXY',
                        ('expected_posintT', 'expected_posintP'): 'posintTP',
                        ('expected_poslocT', 'expected_poslocP'): 'poslocTP',
                        ('expected_poslocX', 'expected_poslocY'): 'poslocXY',
                        ('expected_ptlX', 'expected_ptlY'): 'ptlXY',
                        }
exp_colnames = ['POS_ID', 'DEVICE_LOC', 'move_idx']
exp_dtypes = [str, int, int]
for key in expected_result_keys.keys():
    exp_colnames += list(key)
    exp_dtypes += [float, float]
expected_results = Table(names=exp_colnames, dtype=exp_dtypes)
orig_params = {}
loop = 0
for pos_param_id, pos_params in pos_param_sequence.items():
    for posid, params in pos_params.items():
        state = posstate.PosState(unit_id=posid, device_type='pos', petal_id=petal_id)
        state.store('POS_T', start_posintTP[0], register_if_altered=False)
        state.store('POS_P', start_posintTP[1], register_if_altered=False)
        for key,val in params.items():
            state.store(key, val, register_if_altered=False)
        if ignore_params_ctrl_enabled:
            state.store('CTRL_ENABLED', True, register_if_altered=False)
        if params['DEVICE_LOC'] in retract_and_disable:
            orig_params[posid] = {key:state._val[key] for key in ['CTRL_ENABLED', 'CLASSIFIED_AS_RETRACTED']}
            state.store('POS_T', retracted_TP[0], register_if_altered=False)
            state.store('POS_P', retracted_TP[1], register_if_altered=False)
            state.store('CTRL_ENABLED', False, register_if_altered=False)
            state.store('CLASSIFIED_AS_RETRACTED', True, register_if_altered=False)
        state.write()
    ptl = petal.Petal(petal_id        = petal_id,
                      petal_loc       = 3,
                      posids          = pos_params.keys(),
                      fidids          = fidids,
                      simulator_on    = True,
                      db_commit_on    = False,
                      local_commit_on = False,
                      local_log_on    = False,
                      collider_file   = None,
                      sched_stats_on  = True, # minor speed-up if turn off
                      anticollision   = 'adjust',
                      verbose         = False,
                      phi_limit_on    = False,
                      save_debug      = True,
                      anneal_mode     = 'ramped',
#                      auto_disabling_on = True,
                      )
    if expand_keepouts:
        ptl.set_keepouts(posids='all', radT=0.05, radP=0.05, angT=5, angP=10)
    for key, val in sim_fail_freq.items():
        ptl.sim_fail_freq[key] = val
    ptl.limit_radius = None
    if ptl.schedule_stats.is_enabled():
        if new_stats_per_loop:
            loop_suffix = f'_loop{loop}' if len(pos_param_sequence) > 1 else ''
            stats_path = os.path.join(pc.dirs['temp_files'], f'sched_stats_{filename_suffix}{loop_suffix}.csv')
            ptl.sched_stats_path = stats_path
        ptl.schedule_stats.clear_cache_after_save = False
        ptl.schedule_stats.add_note('POS_PARAMS_ID: ' + str(pos_param_id))
    for posid in scale_changes:
        if posid not in ptl.posids:
            continue
        for axis, scale in scale_changes[posid].items():
            ptl.set_posfid_val(posid, key=f'GEAR_CALIB_{axis}', value=scale)
        if set_scale_changes_as_retracted:
            ptl.set_posfid_val(posid, key='CLASSIFIED_AS_RETRACTED', value=True)
    ptl.collider.refresh_calibrations()
    if should_animate:
        ptl.animator.cropping_on = True
        ptl.animator.label_size = anim_label_size
        if animation_foci == 'colliding':
            ptl.collider.animate_colliding_only = True  
        elif animation_foci == 'commanded':
            ptl.collider.posids_to_animate = {ptl.devices[loc] for loc in device_loc_to_command}
            ptl.collider.fixed_items_to_animate = set()
        elif animation_foci != 'all' and len(animation_foci) > 0:
            posids_to_animate = set()
            fixed_items_to_animate = set()
            for identifier in animation_foci:
                if identifier in ptl.posids:
                    posids_to_animate.add(identifier)
                    posids_to_animate.update(ptl.collider.pos_neighbors[identifier])
                elif identifier in ['PTL','GFA']:
                    fixed_items_to_animate.add(identifier)
                else:
                    print('Warning: ' + str(identifier) + ', requested as an animation focus, is not a known item on this petal, therefore was ignored.')
                ptl.collider.posids_to_animate = posids_to_animate
                ptl.collider.fixed_items_to_animate = fixed_items_to_animate
        ptl.start_gathering_frames()
    m = 0
    mtot = len(move_request_sequence)
    for move_requests_id, move_request_data in move_request_sequence.items():
        m += 1
        if ptl.schedule_stats.is_enabled():
            ptl.schedule_stats.add_note('MOVE_REQUESTS_ID: ' + str(move_requests_id))
        for n in range(n_corrections + 1):
            print(' move: ' + str(m) + ' of ' + str(mtot) + ', submove: ' + str(n))
            log_note = f'harness move {m} submove {n}'
            requests = {}
            results_rows = {}
            export_targets_this_submove = should_export_targets and n == 0
            if export_targets_this_submove:
                commands = {d['command'] for d in move_request_data.values()}
                assert len(commands) == 1
                command = update_command_format(commands.pop())
                hw_move_args = {'command': command,
                                'target0': [],
                                'target1': [],
                                'posids': [],
                                'log_note': runstamp,
                                'allow_corr': True}
            for loc_id, data in move_request_data.items():
                data['command'] = update_command_format(data['command'])
                if loc_id in ptl.devices:
                    posid = ptl.devices[loc_id]
                    command = data['command']
                    target = [data['u'], data['v']]
                    
                    # hack to dynamically contract targets when testing retracted + scale change
                    if posid in scale_changes and ptl.get_posfid_val(posid, 'CLASSIFIED_AS_RETRACTED') and command=='poslocXY':
                        radius = np.hypot(*target)
                        if radius > 3.2:
                            new_radius = random.random() * 3.2
                            angle = np.arctan2(target[1], target[0])
                            target = [new_radius * np.cos(angle), new_radius* np.sin(angle)]
                            
                    requests[posid] = {'command':command, 'target':target, 'log_note':log_note}
                    if export_targets_this_submove:
                        hw_move_args['target0'].append(target[0])
                        hw_move_args['target1'].append(target[1])
                        hw_move_args['posids'].append(posid)
                    trans = ptl.posmodels[posid].trans
                    results_row = {key:None for key in expected_results.columns}
                    results_row['move_idx'] = m - 1
                    results_row['POS_ID'] = posid
                    results_row['DEVICE_LOC'] = loc_id
                    for keys, coord in expected_result_keys.items():
                        if 'target' in keys[0]:
                            if coord == command:
                                expected = target
                            else:
                                transform_func = trans.construct(coord_in=data['command'], coord_out=coord)
                                expected = transform_func(target)
                                if isinstance(expected[1], (bool, np.bool_)):
                                    expected = expected[0] # because some transform funcs return "unreachable" bool as a 2nd tuple element
                            results_row[keys[0]] = expected[0]
                            results_row[keys[1]] = expected[1] 
                    results_rows[posid] = results_row
            anticollision = 'adjust'
            if n > 0:
                for request in requests.values():
                    request['command'] = 'poslocdXdY'
                    request['target'][0] = random.uniform(-max_correction_move,max_correction_move)
                    request['target'][1] = random.uniform(-max_correction_move,max_correction_move)
                    request['log_note'] = log_note
                anticollision = 'freeze'
            if should_profile:
                hc.profile('ptl.request_targets(requests)')
            else:
                ptl.request_targets(requests)
            if should_profile:
                hc.profile('ptl.schedule_send_and_execute_moves(anticollision="'+anticollision+'")')
            else:
                ptl.schedule_send_and_execute_moves(anticollision=anticollision)
            if should_inspect_some_TP:
                for dev in device_loc_to_command:
                    posid = ptl.devices[dev]
                    print('---------------')
                    print(f'posid {posid}, device {dev}')
                    vals = ptl.states[posid]._val
                    print(f'POS_T = {vals["POS_T"]}')
                    print(f'POS_P = {vals["POS_P"]}')
                    print(f'OFFSET_T = {vals["OFFSET_T"]}')
                    print(f'OFFSET_P = {vals["OFFSET_P"]}')
                    print(f'POS_T + OFFSET_T = {vals["POS_T"] + vals["OFFSET_T"]}')
                    print(f'POS_P + OFFSET_P = {vals["POS_P"] + vals["OFFSET_P"]}')
                print('---------------')
        if should_export_targets:
            move = move_sequence.Move(**hw_move_args)
            hw_seq.append(move)
        for posid, row in results_rows.items():
            posmodel = ptl.posmodels[posid]
            current = posmodel.expected_current_position
            for keys, coord in expected_result_keys.items():
                if 'expected' in keys[0]:
                    row[keys[0]] = current[coord][0]
                    row[keys[1]] = current[coord][1]
            expected_results.add_row(row)
        if test_direct_dTdP:
            posids_to_test = list(requests.keys())
            for dtdp in [[30,0], [0,-30], [-30,30]]:
                print(f'direct_dTdP {dtdp}')
                direct_requests = {posid: {'target': dtdp, 'log_note':''} for posid in posids_to_test}
                ptl.request_direct_dtdp(direct_requests)
                ptl.schedule_send_and_execute_moves(anticollision='adjust') # 'adjust' here *should* internally be ignored in favor of 'freeze'
        if test_homing:
            posids_to_test = list(requests.keys())
            for axis in ['phi_only', 'theta_only', 'both']:
                print(f'homing {axis}')
                ptl.request_homing(posids_to_test, axis=axis)
                ptl.schedule_send_and_execute_moves(anticollision='adjust') # 'adjust' here *should* internally be ignored in favor of 'freeze'
    if ptl.schedule_stats.is_enabled():
        ptl.schedule_stats.save(path=ptl.sched_stats_path, footers=True)
        print(f'Stats saved to {ptl.sched_stats_path}')
    if should_export_targets:
        hw_seq.save(directory=pc.dirs['temp_files'])
        filename = 'expected_' + filename_suffix + '.csv'
        path = os.path.join(pc.dirs['temp_files'], filename)
        expected_results.write(path)
    if should_animate and not ptl.animator.is_empty():
        ptl.stop_gathering_frames()
        print('Generating animation (this can be quite slow)...')
        ptl.animator.filename_suffix = filename_suffix
        ptl.animator.add_timestamp_prefix_to_filename = False
        ptl.generate_animation()
    for posid, params in orig_params.items():
        for key,value in params.items():
            ptl.states[posid].store(key, value, register_if_altered=False)
        ptl.states[posid].write()
    loop += 1
