# -*- coding: utf-8 -*-
"""
Produces positioner move sequences, for any of several pre-cooked types of test.
"""

import sequence
import numpy as np
import itertools
import math
import os

seqs = []
make_sparse_csv = []

# DEBUGGING / CODE SYNTAX seqs
# -----------------------------

seq = sequence.Sequence(short_name='debug dTdP',
                        long_name='test the code with single tiny move',
                        pos_settings={'ANTIBACKLASH_ON': False}
                        )
move = sequence.Move(command='dTdP', target0=0.01, target1=0.01, log_note='', allow_corr=False)
seq.append(move)
seqs.append(seq)

seq = sequence.Sequence(short_name='debug dTdP full current',
                        long_name='test the code with single tiny move and motor current set to 100',
                        pos_settings={'ANTIBACKLASH_ON': False,
                                      'CURR_SPIN_UP_DOWN': 100,
                                      'CURR_CRUISE': 100,
                                      'CURR_CREEP': 100,
                                      },
                        )
move = sequence.Move(command='dTdP', target0=0.01, target1=0.01, log_note='', allow_corr=False)
seq.append(move)
seqs.append(seq)

x = 1.0
y = 1.0
seq = sequence.Sequence(short_name='debug poslocXY with corr',
                        long_name=f'single move to poslocXY = ({x}, {y}) and a followup correction move',
                        )
move = sequence.Move(command='poslocXY', target0=x, target1=y, log_note='', allow_corr=True)
seq.append(move)
seqs.append(seq)

t = 0.0
p = 150.0
seq = sequence.Sequence(short_name='debug posintTP with corr',
                        long_name=f'single move to posintTP = ({t}, {p}) and a followup correction move',
                        )
move = sequence.Move(command='posintTP', target0=t, target1=p, log_note='', allow_corr=True)
seq.append(move)
seqs.append(seq)
make_sparse_csv.append(seq)

seq = sequence.Sequence(short_name='debug multitarget',
                        long_name='test the code with differing simultaneous targets on several positioners',
                        )
posids = [f'M{i:05}' for i in range(10000)]
targ_val_options = [i*.5-2.0 for i in range(9)]
targ_val_unique_combos = list(itertools.combinations(targ_val_options,2))
fill_repeats = math.ceil(len(posids) / len(targ_val_unique_combos))
targ_val_combos = []
for i in range(fill_repeats):
    targ_val_combos += targ_val_unique_combos
targets = [[], []]
for axis in [0,1]:
    targets[axis] = [targ_val_combos[k][axis] for k in range(len(posids))]
n_moves = 2
shift = lambda X: X[1:] + [X[0]]
for i in range(n_moves):
    for axis in range(len(targets)):
        targets[axis] = shift(targets[axis])
    move = sequence.Move(command='poslocXY',
                         target0=targets[0],
                         target1=targets[1],
                         posids=posids,
                         log_note='',
                         allow_corr=True,
                         )
    seq.append(move)
seqs.append(seq)


# GENERIC XY TESTS
# ----------------
import xy_targets_generator
n_targs = [4, 24, 100]
calib_seqs = {}
for n in n_targs:
    for limited in [True]:  # [JHS] as of 2020-08-27 I'm not yet releasing the unlimited version into the wild, until anticollision is well-tested
        seq = sequence.Sequence(short_name=f'xytest uniform{" limited" if limited else ""} {n}',
                                long_name=f'rectilinear grid of test points, same local xy for all pos{", limited patrol" if limited else ""}',
                                )
        targs = xy_targets_generator.filled_annulus(n_points=n,
                                                    r_min=0.0,
                                                    r_max=3.1 if limited else 6.0,
                                                    random=False)
        for targ in targs:
            move = sequence.Move(command='poslocXY', target0=targ[0], target1=targ[1], allow_corr=True)
            seq.append(move)
        seqs.append(seq)


# BASIC HOMING SEQUENCES
# ----------------------

seq = sequence.Sequence(short_name='home_and_debounce',
                        long_name='run a single rehome on both axes, followed by debounce moves',
                        )
move = sequence.Move(command='home_and_debounce', target0=1, target1=1, log_note='', allow_corr=False)
seq.append(move)
seqs.append(seq)

seq = sequence.Sequence(short_name='home_no_debounce', 
                        long_name='run a single rehome on both axes, with no debounce moves',
                        )
move = sequence.Move(command='home_no_debounce', target0=1, target1=1, log_note='', allow_corr=False)
seq.append(move)
seqs.append(seq)


# SIMPLE ARC SEQUENCES
# --------------------

cmd = 'posintTP'
settings = {'ALLOW_EXCEED_LIMITS': True}
name_note = ', travel limits OFF'

def arc_targets(axis, n_points, start, step):
    if axis == 'theta':
        thetas = [-170+i*20 for i in range(18)]
        phi = 130
        targets = [[theta, phi] for theta in thetas]
    else:
        theta = 0
        phis = [120+i*3 for i in range(18)]
        targets = [[theta, phi] for phi in phis]
    return targets
targets = {'theta': arc_targets('theta', n_points=18, start=-170, step=20),
           'phi': arc_targets('phi', n_points=18, start=120, step=3),
           }
for axis, targs in targets.items():
    seq = sequence.Sequence(short_name=f'arc {axis}',
                            long_name=f'rotate {axis} repeatedly, for use in circle fits{name_note}',
                            pos_settings=settings)
    for i in range(len(targs)):
        target = targs[i]
        move = sequence.Move(command=cmd, target0=target[0], target1=target[1], log_note='', allow_corr=False)
        seq.append(move)
    seqs.append(seq)
    
cmd = 'dTdP'
deltas = [1.0 for i in range(10)]
for axis in ['theta', 'phi']:
    seq = sequence.Sequence(short_name=f'shortdeltas {axis}',
                            long_name=f'rotate {axis} small delta amounts, over a short distance{name_note}',
                            pos_settings=settings,
                            )
    if axis == 'theta':
        targets = [[delta, 0] for delta in deltas]
    else:
        targets = [[0, delta] for delta in deltas]
    for i in range(len(targets)):
        target = targets[i]
        move = sequence.Move(command=cmd, target0=target[0], target1=target[1], log_note='', allow_corr=False)
        seq.append(move)
    seqs.append(seq)


# CALIBRATION SEQUENCE
# --------------------
seq = sequence.Sequence(short_name='RC calib',
                        long_name='rehome and calibrate',
                        details='typical sequence of hitting hard limits, then running arcs on theta and phi, then a regular cartesian grid',
                        )
move = sequence.Move(command='home_and_debounce', target0=1, target1=1, log_note='homing', allow_corr=False)
seq.append(move)
targets = {'xy': xy_targets_generator.filled_annulus(n_points=24, r_min=0.0, r_max=3.1, random=False),
           'theta': arc_targets('theta', n_points=18, start=-170, step=20),
           'phi': arc_targets('phi', n_points=18, start=120, step=3),
           }
for key in ['theta', 'phi', 'xy']:
    targs = targets[key]
    cmd = 'posintTP' if key in {'theta', 'phi'} else 'poslocXY'
    for i in range(len(targs)):
        target = targs[i]
        move = sequence.Move(command=cmd, target0=target[0], target1=target[1],
                             log_note=f'{key} {"arc" if key in ["theta", "phi"] else "grid"}',
                             allow_corr=False)
        seq.append(move)
seqs.append(seq)


# MOTOR TESTS
# -----------

# Helper functions
def wiggle(forward, case=0):
    '''Generates a wiggly sequence with lots of back and forth, using a starting list
    of forward direction deltas. Several pre-cooked cases are provided.'''
    backward = [-x for x in forward]
    forward_reversed = [x for x in reversed(forward)]
    backward_reversed = [x for x in reversed(backward)]
    backandforth = [x for t in zip(forward, backward) for x in t]
    if case == 0:
        deltas = forward + backward + backward_reversed + forward_reversed + backandforth
    if case == 1:
        deltas = backward + backandforth + forward + backward_reversed + forward_reversed
    return deltas
    
def typ_motortest_sequence(prefix, short_suffix, long_suffix, details, forward_deltas, settings):
    seq = sequence.Sequence(short_name=f'motortest {prefix} {short_suffix}',
                            long_name=f'{prefix} motor test, {long_suffix}',
                            details=details,
                            pos_settings=settings,
                            )
    i = 0 if 'THETA' in seq.normalized_short_name else 1
    case = 0  # at one point I had case=i, to do phi differently from theta
    deltas = wiggle(forward_deltas, case=case)
    for j in range(len(deltas)):
        target = [0,0]
        target[i] = deltas[j]
        move = sequence.Move(command='dTdP', target0=target[0], target1=target[1], log_note='', allow_corr=False)
        seq.append(move)
    return seq

option_groups = {'NOMINAL': {},
                 'CRUISEONLY':
                     {'FINAL_CREEP_ON': False,
                      'MIN_DIST_AT_CRUISE_SPEED': sequence.nominals['stepsize_cruise'], # smallest finite value
                      'SPINUPDOWN_PERIOD': 8, # to keep total accel+decel distance < backlash. i.e. 8*_spinupdown_dist_per_period*2 + 3.3)/337 = 2.67 < 3.0
                      'BACKLASH': 3.0,  # just making sure, though this is almost always already the default
                      },
                 'CRUISEONLY_NOANTIBACKLASH':
                     {'FINAL_CREEP_ON': False,
                      'ANTIBACKLASH_ON': False,
                      'MIN_DIST_AT_CRUISE_SPEED': sequence.nominals['stepsize_cruise'], # smallest finite value
                      },                     
                 'CRUISEONLY_NOSPINUPDOWN':
                     {'FINAL_CREEP_ON': False,
                      'MIN_DIST_AT_CRUISE_SPEED': sequence.nominals['stepsize_cruise'], # smallest finite value
                      'CURR_SPIN_UP_DOWN': 0,
                      'SPINUPDOWN_PERIOD': 1, # smallest finite value
                      },
                 'CREEPONLY':
                     {'ONLY_CREEP': True,
                      'FINAL_CREEP_ON': False,
                      },
                 'FASTCREEP':
                     {'ONLY_CREEP': True,
                      'CREEP_PERIOD': 1,
                      'FINAL_CREEP_ON': False,
                      },
                 }

# Common values for "wiggle" inputs
forward_deltas = {'cruise0': [1, 4, 16, 32],
                  'cruise1': [1, 4, 8, 12],
                  'creep': [0.5, 1.0, 1.5, 2.0],
                  }

# Theta and phi performance at nominal settings
details = '''Setup: default
Moves: Several moves at cruise speed, each followed by the usual final creep moves for precision positioning and antibacklash.'
Purpose: Baseline tests for typical proper function.'''
options = {}
prefix = 'nominal'
long_suffix = 'performance at nominal settings'
seqs.append(typ_motortest_sequence(prefix, 'Theta', long_suffix, details, forward_deltas['cruise0'], options))
seqs.append(typ_motortest_sequence(prefix,   'Phi', long_suffix, details, forward_deltas['cruise1'], options))

# Theta and phi cruise-only at otherwise nominal settings
details = '''Setup: Turn off parameters FINAL_CREEP_ON and ANTIBACKLASH_ON
Moves: Several moves at cruise speed. In each direction, at several step sizes.
Purpose: Measure the effective output ratio in typical cruise mode.'''
prefix = 'cruiseonly'
options = option_groups[prefix.upper()]
long_suffix = 'cruise-only, at otherwise nominal settings'
seqs.append(typ_motortest_sequence(prefix, 'Theta', long_suffix, details, forward_deltas['cruise0'], options))
seqs.append(typ_motortest_sequence(prefix,   'Phi', long_suffix, details, forward_deltas['cruise1'], options))

# Theta and phi cruise-only, with spinup/down power disabled
details = '''Setup: Turn off parameters FINAL_CREEP_ON and ANTIBACKLASH_ON. Set CURR_SPIN_UP_DOWN = 0.
Moves: Several moves at cruise speed. In each direction, at several step sizes.
Purpose: Measure the effective output ratio in typical cruise mode.'''
prefix = 'cruiseonly nospinupdown'
options = option_groups[prefix.upper().replace(' ','_')]
long_suffix = 'cruise-only, with spinup/down power disabled'
seqs.append(typ_motortest_sequence(prefix, 'Theta', long_suffix, details, forward_deltas['cruise0'], options))
seqs.append(typ_motortest_sequence(prefix,   'Phi', long_suffix, details, forward_deltas['cruise1'], options))

# Theta and phi creep-only at otherwise nominal settings
details = '''Setup: Turn on parameter ONLY_CREEP.
Moves: Several moves at creep speed. In each direction, at several step sizes.
Purpose: Measure the effective output ratio in creep mode.'''
prefix = 'creeponly'
options = option_groups[prefix.upper()]
long_suffix = 'creep-only, at otherwise nominal settings'
seqs.append(typ_motortest_sequence(prefix, 'Theta', long_suffix, details, forward_deltas['creep'], options))
seqs.append(typ_motortest_sequence(prefix,   'Phi', long_suffix, details, forward_deltas['creep'], options))

# Theta and phi fast creep
details = '''Setup: Halve the creep period thereby doubling the creep speed.
Moves: Several moves at creep speed. In each direction, at several step sizes.
Purpose: Determine whether creep performance can be improved under existing firmware (v5.0) constraints.'''
prefix = 'fastcreep'
options = option_groups[prefix.upper()]
long_suffix = 'with fastest available creep speed under firmware v5.0'
seqs.append(typ_motortest_sequence(prefix, 'Theta', long_suffix, details, forward_deltas['creep'], options))
seqs.append(typ_motortest_sequence(prefix,   'Phi', long_suffix, details, forward_deltas['creep'], options))

# Simplest-possible test for SCALE of moves
details = '''Setup: Turn off parameters FINAL_CREEP_ON and ANTIBACKLASH_ON.
Moves: Several moves at cruise speed "forward" then "backward", with no creep.
Purpose: Determine whether the positioner has scale error in its output ratio.'''
options = {'FINAL_CREEP_ON': False,
           'ANTIBACKLASH_ON': False,
           'MIN_DIST_AT_CRUISE_SPEED': sequence.nominals['stepsize_cruise'], # smallest finite value
          }
prefix = 'scale'
for suffix, long_suffix in {'': 'cruise-only, forward and back',
                            'multi': 'cruise-only, multiple step sizes',
                           }.items():
    for axis in {'THETA', 'PHI'}:
        spaced_suffix = ' ' + suffix if suffix else ''
        seq = sequence.Sequence(short_name=f'motortest {prefix} {axis[0]}{spaced_suffix}',
                                long_name=f'output {prefix} motor test, {long_suffix}',
                                details=details,
                                pos_settings=options,
                                )
        if axis == 'THETA':
            if suffix == '':
                forward_deltas = [10 for i in range(4)]
                reverse_deltas = [-x for x in forward_deltas]
                deltas = forward_deltas + reverse_deltas + reverse_deltas + forward_deltas
            elif suffix == 'multi':
                forward_deltas = [5, 10, 15, -20, 30, -80, 40]
                reverse_deltas = [-x for x in forward_deltas]
                deltas = forward_deltas + reverse_deltas + reverse_deltas + forward_deltas
            else:
                assert False, f'unrecognized suffix identifier {suffix}'
            i = 0
        else:
            if suffix == '':
                forward_deltas = [-8 for i in range(5)]
                reverse_deltas = [-x for x in forward_deltas]
                deltas = forward_deltas + reverse_deltas
            elif suffix == 'multi':
                forward_deltas = [5, 8, 15, -20, 32]
                reverse_deltas = [-x for x in forward_deltas]
                deltas = forward_deltas + [-40, 40] + reverse_deltas
                deltas += deltas
            else:
                assert False, f'unrecognized suffix identifier {suffix}'
            i = 1
        for j in range(len(deltas)):
            target = [0,0]
            target[i] = deltas[j]
            move = sequence.Move(command='dTdP', target0=target[0], target1=target[1], log_note='', allow_corr=False)
            seq.append(move)
        seqs.append(seq)

# simple positioning tests for robots with calibrated scale errors
for axis, targets in {'THETA': [[x, 130] for x in [0, +60, -60, +120, -120]],
                      'PHI': [[0, x] for x in [150, 160, 140, 170, 130]],
                      }.items():
    for key, options in option_groups.items():
        seq = sequence.Sequence(short_name=f'SIMPLESCALE_{key}_{axis}',
                                long_name=f'Simplified test for positioner(s) with {axis} output scale != 1.0',
                                details='',
                                pos_settings=options.copy(),
                                )
        for target in targets:
            move = sequence.Move(command='poslocTP',
                                 target0=target[0],
                                 target1=target[1],
                                 log_note='',
                                 allow_corr=True,
                                 )
            seq.append(move)
        seqs.append(seq)

# SAVE ALL TO DISK
# ----------------
paths = []
sparse_paths = []
for seq in seqs:
    path = seq.save()
    paths.append(path)
    if seq in make_sparse_csv:
        table = seq.to_table()
        all_cols = set(table.columns)
        keep_cols = {sequence.move_idx_key, 'command', 'target0', 'target1', 'posids'}
        for col in all_cols - keep_cols:
            del table[col]
        sparse_path = os.path.splitext(path)[0] + '_sparse.csv'
        table.write(sparse_path, overwrite=True, delimiter=',')
        sparse_paths.append(sparse_path)
    paths.extend(sparse_paths)
    
# READ FROM DISK AND PRINT TO STDOUT
# ----------------------------------
for path in paths:
    if 'XYTEST' in path:
        print(f'Now reading: {path}')  # debug breakpoint
    seq = sequence.Sequence.read(path)
    print(seq,'\n')
    if 'MOTORTEST' in path:
        axis = 0 if 'THETA' in seq.normalized_short_name or 'SCALE_T' in seq.normalized_short_name else 1
        deltas = [getattr(move, f'target{axis}') for move in seq]
        cumsum = np.cumsum(deltas).tolist()
        print('num test points: '  + str(len(deltas)))
        print(f'min and max excursions: [{min(cumsum)} deg, {max(cumsum)} deg]')
        print('\n')
        # print('delta sequence: ' + str(deltas))
        # print('running total: ' + str(cumsum))
    
    
