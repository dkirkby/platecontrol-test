#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Produces move request sequence. In particular for testing a specific group
of positioners, using their particular calibration parameters to do some intelligent
target selection.
"""

# command line argument parsing
import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-n', '--num_moves', type=int, required=True, help='integer, number of moves to generate')
parser.add_argument('-ptl', '--petal_ids', type=str, default='kpno', help='Comma-separated integers, specifying one or more PETAL_ID number(s) for which to retrieve data. Defaults to all known petals at kpno. Alternately argue "lbnl" for all known petals at lbnl.')
parser.add_argument('-pos', '--posids', type=str, default='all', help='optional, comma-separated POS_ID strings, saying which positioners to generate targets for (defaults to "all")')
parser.add_argument('-s', '--num_stress_select', type=int, default=1, help='optional, integer, for every selected move, code will internally test this many moves and pick the one with the most opportunities for collision.')
parser.add_argument('-ai', '--allow_interference', action='store_true', help='optional, allows targets to interfere with one another')
parser.add_argument('-lim', '--enable_phi_limit', action='store_true', help='optional, turns on minimum phi limit for move targets, default is False')
parser.add_argument('-nb', '--no_buffer', action='store_true', help='optional, turns off extra buffer space added to polygons for target selection')
parser.add_argument('-p', '--profile', action='store_true', help='optional, turns on timing profiler for move scheduling')
parser.add_argument('-i', '--infile', type=str, default=None, help='optional, either "cache" to use the most recently cached calib data, or else a path to an offline csv file, containing positioner calibration parameters. If not argued, will try getting current values from online db instead (for this you may have to be running from a machine at kpno or beyonce)')
parser.add_argument('-o', '--outdir', type=str, default='.', help='optional, path to directory where to save output file, defaults to current dir')
parser.add_argument('-m', '--comment', type=str, default='', help='optional, comment string (must enclose in "") that will be included in output file metadata')
parser.add_argument('-np', '--n_processes_max', type=int, default=None,  help='max number of processors to use')
parser.add_argument('-d', '--debug_mode', action='store_true',  help='restricts processors and any other debug options')
uargs = parser.parse_args()

import os
save_dir = os.path.realpath(uargs.outdir)
assert uargs.num_moves > 0
assert uargs.num_stress_select > 0
assert uargs.n_processes_max == None or uargs.n_processes_max > 0, f'unrecognized arg {uargs.n_processes_max} for n_processes_max'

# proceed with bulk of imports
from astropy.table import Table
import sequence
import random
import numpy as np
import multiprocessing

# desi imports
try:
    import petal
except:
    import sys
    path_to_petal = '../petal'
    sys.path.append(os.path.abspath(path_to_petal))
    import petal
import posconstants as pc

# required data fields
required = {'POS_ID', 'DEVICE_LOC', 'LENGTH_R1', 'LENGTH_R2', 'OFFSET_T',
            'OFFSET_P', 'OFFSET_X', 'OFFSET_Y', 'PHYSICAL_RANGE_T', 'PHYSICAL_RANGE_P',
            'KEEPOUT_EXPANSION_PHI_RADIAL', 'KEEPOUT_EXPANSION_PHI_ANGULAR',
            'KEEPOUT_EXPANSION_THETA_RADIAL', 'KEEPOUT_EXPANSION_THETA_ANGULAR',
            'CLASSIFIED_AS_RETRACTED', 'CTRL_ENABLED'}

# Nominal polygon buffers. These are temporary angular and radial keepout expansions
# applied during target. They add to any existing expansions --- not override.
buffers = {'KEEPOUT_EXPANSION_PHI_RADIAL': 0.05,  # mm
           'KEEPOUT_EXPANSION_PHI_ANGULAR': 6.0,  # deg
           'KEEPOUT_EXPANSION_THETA_RADIAL': 0.0,  # mm
           'KEEPOUT_EXPANSION_THETA_ANGULAR': 3.0,  # deg
           }
if uargs.no_buffer:
    buffers = {}
    print('No keepout expansion buffers during target generation.')
else:
    print(f'Keepout expansion buffers during target generation: {buffers}')

# generate calibrations file or resolve location
calibs_file_path = uargs.infile
if calibs_file_path == None:
    outdir = pc.dirs['temp_files']
    get_calibs_comment = 'auto-retrieval of calibrations by make_sequence.py'
    print('\nGENERATING NEW CALIBRATIONS FILE')
    cmd = f'python get_calibrations.py -m "{get_calibs_comment}" -o {outdir} -ptl {uargs.petal_ids}'
    err = os.system(cmd)
    assert not(err), f'error calling \'{cmd}\''
if calibs_file_path in {None, 'cache'}:
    with open(pc.fp_calibs_path_cache, 'r') as file:
        calibs_file_path = file.read()
assert os.path.exists(calibs_file_path), f'calibration data file not found at path {calibs_file_path}'

# read positioner parameter data from csv file
booleans = {'CLASSIFIED_AS_RETRACTED', 'CTRL_ENABLED'}
nulls = {'--', None, 'none', 'None', '', 0, False, 'False', 'FALSE', 'false'}
boolean = lambda x: x not in nulls
print(f'\nREADING CALIBRATIONS FILE: {calibs_file_path}\n')
input_table = Table.read(calibs_file_path)
missing = required - set(input_table.columns)
assert len(missing) == 0, f'input positioner parameters file is missing columns {missing}'
for col in booleans:
    input_table[col] = [boolean(x) for x in input_table[col].tolist()]
print('calibration parameters: read from file complete!')

# determine which petals to make sequences for, and split up the data
petal_ids = set(input_table['PETAL_ID'])
tables = {}
for petal_id in petal_ids:
    this = input_table.copy()
    undesired = this['PETAL_ID'] != petal_id
    this.remove_rows(undesired)
    tables[petal_id] = this

# define multi-processable function for generating sequences
def make_sequence_for_one_petal(table):
    petal_ids = set(table['PETAL_ID'])
    assert len(petal_ids) == 1, 'table should contain only one, uniform petal_id'
    petal_id = petal_ids.pop()

    # ensure single, unique set of parameters per positioner
    all_posids = set(table['POS_ID'])
    params_dict = {key:[] for key in required}
    if 'DATE' in table.columns:
        table.sort('DATE', reverse=True)
    for posid in all_posids:
        row = table[table['POS_ID'] == posid][0]
        for key in params_dict:
            params_dict[key].append(row[key])
    params = Table(params_dict)
    
    # ensure single, unique positioner in any given device locations per petal
    possible_locs = set(pc.generic_pos_neighbor_locs)
    device_locs = params['DEVICE_LOC']
    invalid_locs = set(params['DEVICE_LOC']) - possible_locs
    assert len(invalid_locs) == 0, f'invalid device locations in input file: {invalid_locs}'
    assert len(device_locs) == len(set(device_locs)), f'found repeated device location(s) in input file for petal {petal_id}'
    
    # select which positioners get move commands
    if uargs.posids == 'all':
        movers = all_posids
    else:
        movers = set(uargs.posids.split(','))
        missing = movers - all_posids
        assert len(missing) == 0, f'argued posids {missing} were not found in parameters file'
    disabled = set(params['POS_ID'][params['CTRL_ENABLED'] == False])
    movers -= disabled
    movers = sorted(movers)  # code lower down assumes consistent order of these
    
    # reduce positioners collection to just movers + their neighbors, for speed-up
    # in cases where not sequencing the whole petal, etc
    loc2id = {params['DEVICE_LOC'][i]: params['POS_ID'][i] for i in range(len(params))}
    id2loc = {val: key for key, val in loc2id.items()}
    reduced = set(movers)
    for posid in movers:
        loc = id2loc[posid]
        neighbor_locs = pc.generic_pos_neighbor_locs[loc]
        reduced |= {loc2id[loc] for loc in neighbor_locs if loc in loc2id}
    eliminate = all_posids - reduced
    elim_rows = [np.where(params['POS_ID'] == posid)[0][0] for posid in eliminate]
    params.remove_rows(elim_rows)
    params.sort('POS_ID')
    
    # initialize a simulated petal instance
    print(f'Now initializing petal {petal_id}, will take a few seconds...')
    ptl = petal.Petal(petal_id        = petal_id,
                      petal_loc       = 3,
                      posids          = params['POS_ID'].tolist(),
                      fidids          = set(),
                      simulator_on    = True,
                      db_commit_on    = False,
                      local_commit_on = False,
                      local_log_on    = False,
                      collider_file   = None,
                      sched_stats_on  = True,
                      anticollision   = 'adjust',
                      verbose         = False,
                      phi_limit_on    = uargs.enable_phi_limit,
                      petalbox_id     = -1, # [2021-01-14 JHS] horrible hack, actual intent is to skip trying to read non-existent config file when at KPNO
                      save_debug      = False,
                      )
    
    # set up calibration parameters
    assert ptl.collider.use_neighbor_loc_dict, 'must configure petal to initialize using neighbor locs from file, since calib params not yet set'
    keys_to_store = set(params.columns) - {'POS_ID'}
    for posid, state in ptl.states.items():
        row = params[params['POS_ID'] == posid][0]
        for key in keys_to_store:
            state.store(key, row[key])
    ptl.collider.refresh_calibrations()
    
    def set_keepout_buffers(posids, sign=1):
        '''Set keepout expansion buffers. Sign is +1 (add) or -1 (subtract) from existing.'''
        if not buffers:
            return
        for key, val in buffers.items():
            for posid in posids:
                adjusted = ptl.states[posid]._val[key] + sign * val
                ptl.states[posid].store(key, adjusted, register_if_altered=False)
        ptl.collider.refresh_calibrations()
        print(f'Keepout polygons on {len(posids)} pos adjusted by {sign:+2} * {buffers}.')
    
    # helper functions
    def generate_target_set(posids):
        '''Produces one target set, for the argued collection of posids. Returns a
        dict with keys = posids and values = target locations (poslocXY coordinates).
        '''
        targets_locTP = {}
        targets_posXY = {}
        set_keepout_buffers(posids, sign=+1)
        print(f'Generating target set for {len(posids)} positioners.')
        for posid in posids:
            attempts_remaining = 10000 # just to prevent infinite loop if there's a bug somewhere
            model = ptl.posmodels[posid]
            min_patrol = abs(ptl.collider.R1[posid] - ptl.collider.R2[posid])
            max_patrol = ptl.collider.R1[posid] + ptl.collider.R2[posid]
            if uargs.enable_phi_limit:
                limit_xy = model.trans.poslocTP_to_poslocXY([0, ptl.typical_phi_limit_angle])
                max_patrol = min(max_patrol, np.hypot(limit_xy[0], limit_xy[1]))
            while posid not in targets_locTP and attempts_remaining:
                bad_target = False
                rangeT = model.targetable_range_posintT
                rangeP = model.targetable_range_posintP
                x = random.uniform(-max_patrol, max_patrol)
                y = random.uniform(-max_patrol, max_patrol)
                this_radius = (x**2 + y**2)**0.5
                if min_patrol > this_radius or this_radius > max_patrol:
                    bad_target = True
                else:
                    this_posXY = [x,y]
                    this_posTP, unreachable = model.trans.poslocXY_to_posintTP(this_posXY)
                    if unreachable:
                        bad_target = True
                    this_locTP = model.trans.posintTP_to_poslocTP(this_posTP)
                    if not uargs.allow_interference and not unreachable:
                        target_interference = False
                        for neighbor in ptl.collider.pos_neighbors[posid]:
                            static_neighbor = not ptl.posmodels[neighbor].is_enabled or \
                                              ptl.posmodels[neighbor].classified_as_retracted
                            if neighbor in targets_locTP or static_neighbor:
                                if static_neighbor:
                                    neighbor_locTP = ptl.posmodels[neighbor].expected_current_poslocTP
                                else:
                                    neighbor_locTP = targets_locTP[neighbor]
                                if ptl.collider.spatial_collision_between_positioners(
                                        posid, neighbor, this_locTP, neighbor_locTP):
                                    target_interference = True
                                    break
                        out_of_bounds = ptl.collider.spatial_collision_with_fixed(posid, this_locTP)
                        out_of_range_T = min(rangeT) > this_posTP[0] or max(rangeT) < this_posTP[0]
                        out_of_range_P = min(rangeP) > this_posTP[1] or max(rangeP) < this_posTP[1]
                        if target_interference or out_of_bounds or out_of_range_T or out_of_range_P:
                            bad_target = True
                if bad_target:
                    attempts_remaining -= 1
                else:
                    targets_locTP[posid] = this_locTP
                    targets_posXY[posid] = this_posXY
            if attempts_remaining <= 0:
                v = model.state._val
                print(f'Warning: no valid target found for posid: {posid} at location {v["DEVICE_LOC"]}' +
                      f' (x0, y0) = ({v["OFFSET_X"]:.3f}, {v["OFFSET_Y"]:.3f})!')        
        set_keepout_buffers(posids, sign=-1)
        return targets_posXY
    
    def generate_request_set(targets):
        '''Produces a dict of requests, ready for petal, from a dict with keys = posids
        and values = target locations (poslocXY coordinates).
        '''
        requests = {}
        for posid, target in targets.items():
            requests[posid] = {'command': 'poslocXY', 'target':target, 'log_note':''}
        return requests
    
    def get_current_posTP():
        '''Returns dict with keys=posid, values=posintTP for all positioners.
        '''
        return {posid: model.expected_current_posintTP for posid, model in ptl.posmodels.items()}
    
    def set_posTP(tp):
        '''Set system with POS_T, POS_P values. Input is a dict with keys = posid,
        values=posintTP.
        '''
        for posid, values in tp.items():
            state = ptl.states[posid]
            state.store('POS_T', values[0])
            state.store('POS_P', values[1])
            
    n_stats_lines = 20
    statsfile = os.path.join(pc.dirs['temp_files'], 'stats_make_sequence')
    import cProfile, pstats
    def profile(evaluatable_string):
        print(f'\nProfiling {evaluatable_string}:')
        cProfile.run(evaluatable_string,statsfile)
        p = pstats.Stats(statsfile)
        p.strip_dirs()
        p.sort_stats('cumtime')
        p.print_stats(n_stats_lines)
    
    # generate sequence
    print('\nGENERATING MOVE SEQUENCE')
    print(f'n_moves = {uargs.num_moves}\n')
    seq = sequence.Sequence(short_name=f'temp_ptl{petal_id:02}', long_name='', details='')  # strings will be properly filled after merging all sequences
    set_posTP({posid: (0, 150) for posid in ptl.posids})  # inital values like typical "parked" position
    n_collisions_resolved = []
    for m in range(uargs.num_moves):
        initial_posTP = get_current_posTP()
        candidates = {key: [] for key in ['targets', 'requests', 'final_posTP', 'n_resolved']}
        for s in range(uargs.num_stress_select):
            print(f'Move {m}: Generating candidate target set {s}.')
            candidates['targets'] += [generate_target_set(movers)]
            candidates['requests'] += [generate_request_set(candidates['targets'][-1])]
            requests = ptl.request_targets(candidates['requests'][-1])
            schedule_command = "ptl.schedule_moves(anticollision='default', should_anneal=True)"
            if uargs.profile:
                profile(schedule_command)
            else:
                eval(schedule_command)
            candidates['n_resolved'] += [ptl.schedule_stats.total_resolved]
            ptl.send_and_execute_moves()
            candidates['final_posTP'] += [get_current_posTP()]
            set_posTP(initial_posTP)
        selection = candidates['n_resolved'].index(max(candidates['n_resolved']))
        sel = {k:v[selection] for k,v in candidates.items()}
        n_collisions_resolved += [sel['n_resolved']]
        set_posTP(sel['final_posTP'])
        posids_with_targ = sorted(set(sel['targets']) & set(movers))
        print(f'Move {m}: Targets selected for {len(posids_with_targ)} positioners. Num collisions avoided = {n_collisions_resolved[-1]}')
        move = sequence.Move(command='poslocXY',
                             target0=[sel['targets'][posid][0] for posid in posids_with_targ],
                             target1=[sel['targets'][posid][1] for posid in posids_with_targ],
                             posids=posids_with_targ,
                             log_note='',
                             allow_corr=True,
                             )
        seq.append(move)
    seq.n_collisions_resolved = n_collisions_resolved  # hack, sneaks this value into the sequence meta data
    path = seq.save(pc.dirs['temp_files'])
    return path 

# single/multiprocess switching
n_processes_max = 1 if uargs.debug_mode else uargs.n_processes_max
single_process = n_processes_max == 1
def imap(function, iterable_data):
    '''Common wrapper for pooled single or multiprocess imap'''
    if single_process:
        results = []
        for data in iterable_data:
            result = function(data)
            results.append(result)
    else:
        with multiprocessing.Pool(processes=n_processes_max) as pool:
            results = pool.imap(function, iterable_data)
            pool.close()
            pool.join()
    return results

if __name__ == '__main__':
    big_seq = None
    results = imap(make_sequence_for_one_petal, tables.values())
    for path in results:
        seq = sequence.Sequence.read(path)
        big_seq = seq if big_seq == None else big_seq.merge(seq)
    timestamp = pc.compact_timestamp()
    details = f'Generated with settings: {uargs}'
    details += f'\nCalibrations data source: {calibs_file_path}'
    if uargs.comment:
        details += f'\nComment: {uargs.comment}'
    details += f'\nKeepout polygon buffers: {buffers}'
    posids = big_seq.get_posids()
    big_seq.short_name=f'nptls_{len(tables):02}_npos_{len(posids):03}_ntarg_{uargs.num_moves:03}_{timestamp}'
    big_seq.long_name='Test move sequence generated by make_sequence.py'
    big_seq.details=details
    if not os.path.isdir(save_dir):
        os.path.os.makedirs(save_dir)
    path = big_seq.save(save_dir)
    print('Sequence generation complete!')
    # print(f'\n{big_seq}\n')
    print(f'Saved to file: {path}\n')
