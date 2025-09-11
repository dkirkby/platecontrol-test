#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=fixme, line-too-long, C0103, W0703
"""
move_to_ranges moves positioners to ranges specified in the input csv file which must contain the header
POSID, Theta_Lower_Limit, Theta_Upper_Limit, Phi_Lower_Limit, Phi_Upper_Limit
"""

import sys
import csv
import argparse
from pecs import PECS

max_iter = 10
park_options = ['posintTP', 'poslocTP', 'None', 'False']
default_park = park_options[0]
modes = ['posintTP', 'poslocTP']
default_padding = 5.0

parser = argparse.ArgumentParser()
parser.add_argument('-rc', '--ranges_csv', type=str, default=None, help='Name of csv file containing posids, theta ranges, and phi ranges')
parser.add_argument('-iter', '--iterations', type=int, default=2, help=f'Number of iterations to try to get positioners within desired rang before giving up. Must be less than {max_iter}.')
parser.add_argument('-prep', '--prepark', type=str, default=default_park, help=f'str, controls initial parking move, prior to running does not count as an iteration. Valid options are: {park_options}, default is {default_park}. NOTE: all positioners may be parked, not just those in the CSV file!')
parser.add_argument('-m', '--mode', type=str, default=modes[0], help=f'str, controls which tp coordinate type checks are done. Options: {modes}. Defaults to {modes[0]}.')
parser.add_argument('-p', '--angle_padding', type=float, default=default_padding, help=f'Float, angle padding by which to exceed limits when only one limit is specified. Defaults to {default_padding}.')
parser.add_argument('-d', '--disable', action='store_true', help='Argue to disable positioners that remain out of range at the end of iterations.')
parser.add_argument('-a', '--anticollision', type=str, default='freeze', help='anticollision mode, can be "adjust", "adjust_requested_only", "freeze" or None. Default is "freeze"')
parser.add_argument('-v', '--verbose', action='store_true', help='Argue to print extra messages about the axes to be moved.')
uargs = parser.parse_args()

if uargs.anticollision == 'None':
    uargs.anticollision = None  # fix None string

# basic checks on arguments
assert uargs.anticollision in {'adjust', 'adjust_requested_only', 'freeze', None}, f'bad argument {uargs.anticollision} for anticollision parameter'
assert uargs.iterations < max_iter, f'cannot exceed {max_iter} iterations.'
command = uargs.mode
assert command in modes, f'Invalid mode type {command}. Must be in {modes}.'
assert uargs.prepark in park_options, f'invalid park option, must be one of {park_options}'
uargs.prepark = None if uargs.prepark in ['None', 'False'] else uargs.prepark
rc = uargs.ranges_csv
assert rc is not None, 'Must specify csv file with positioners and ranges'

d_pos_limits = {}

with open(rc, newline='', encoding="utf-8") as f:
    rcsv = csv.DictReader(f)
    try:
        for row in rcsv:
#           print(row)  # DEBUG
#           NOTE order of list: [upper_limit, lower_limit]
            for lmt in ['Theta_Upper_Limit','Theta_Lower_Limit','Phi_Upper_Limit','Phi_Lower_Limit']:
                if str(row[lmt]).upper() == 'NONE':
                    row[lmt] = None
                else:
                    row[lmt] = float(row[lmt])
            if row['Theta_Upper_Limit'] is not None and row['Theta_Lower_Limit'] is not None and row['Theta_Upper_Limit'] < row['Theta_Lower_Limit']:
                t_ul = row['Theta_Upper_Limit']
                row['Theta_Upper_Limit'] = row['Theta_Lower_Limit']
                row['Theta_Lower_Limit'] = t_ul
                if uargs.verbose:
                    print(f"Theta limits reversed for {row['POSID']}, corrected.")
            if row['Phi_Upper_Limit'] is not None and row['Phi_Lower_Limit'] is not None and row['Phi_Upper_Limit'] < row['Phi_Lower_Limit']:
                t_ul = row['Phi_Upper_Limit']
                row['Phi_Upper_Limit'] = row['Phi_Lower_Limit']
                row['Phi_Lower_Limit'] = t_ul
                if uargs.verbose:
                    print(f"Phi limits reversed for {row['POSID']}, corrected.")
            d_pos_limits[row['POSID']] = {'T': [row['Theta_Upper_Limit'], row['Theta_Lower_Limit']], 'P': [row['Phi_Upper_Limit'], row['Phi_Lower_Limit']]}
    except csv.Error as e:
        sys.exit(f'file {rc}, line {rcsv.line_num}: {e}')

posids_to_check = list(d_pos_limits)

# cs = PECS(interactive=False, posids=posids_to_check)
cs = PECS(interactive=True)

# axis_limits = {'T': [tu, tl], 'P': [pu, pl]}
columns = {'T': 'X1', 'P': 'X2'}    # translate to the column name used in the dataframe
zeno_key = {'T': 'ZENO_MOTOR_T', 'P': 'ZENO_MOTOR_P'}   # names of zeno flags for TP axes
zeno_tol = 0.1174 # from posconstants.py schedule_checking_angular_tol_zeno
normal_tol = 0.01 # from posconstants.py schedule_checking_numeric_angular_tol
# t_motors_are_zeno = cs.ptlm.batch_get_posfid_val(uniqueids=posids_to_check, keys=[zeno_key['T']])
# p_motors_are_zeno = cs.ptlm.batch_get_posfid_val(uniqueids=posids_to_check, keys=[zeno_key['P']])
motors_are_zeno = cs.ptlm.batch_get_posfid_val(uniqueids=posids_to_check, keys=[zeno_key['T'], zeno_key['P']])
# print(f't_motors_are_zeno = {t_motors_are_zeno}')
# print(f'p_motors_are_zeno = {p_motors_are_zeno}')
print(f'motors_are_zeno = {motors_are_zeno}')

def axis_is_zeno(aposid, aaxis):
    ''' True if the specified axis is Zeno '''
    retval = False
    for ptl in motors_are_zeno:
        if aposid in motors_are_zeno[ptl]:
            if zeno_key[aaxis] in motors_are_zeno[ptl][aposid]:
                retval = motors_are_zeno[ptl][aposid][zeno_key[aaxis]] is True
                break
    return retval

def axis_tolerance(aposid, aaxis):
    ''' Returns the angular tolerance to use for the specified axis '''
    atol = zeno_tol if axis_is_zeno(aposid, aaxis) else normal_tol
    return atol

def check_if_out_of_limits():
    ''' Check if any axes are out of limits for all positioners in the csv file '''
    all_pos = cs.ptlm.get_positions(return_coord=command, drop_devid=False)
    pos = all_pos[all_pos['DEVICE_ID'].isin(cs.posids)]
    violating_pos = set()
    for c_posid, c_axis_limits in d_pos_limits.items():
        c_cond_a = pos['DEVICE_ID'] == c_posid
        for c_axis, c_limits in c_axis_limits.items():
            c_tol = axis_tolerance(c_posid, c_axis)
            u_violation = False
            l_violation = False
            if c_limits[0] is not None:
                c_limit = c_limits[0] + c_tol
                c_cond_l = pos[columns[c_axis]] > c_limit
                c_df = pos[c_cond_a & c_cond_l]
                if not c_df.empty:
                    c_val = round((pos.loc[c_cond_a & c_cond_l, columns[c_axis]]).iloc[0], 4)
                    if c_val > c_limit:
                        u_violation = True
                    if u_violation and uargs.verbose:
                        print(f'ulimit violation posid={c_posid}, axis={c_axis}, position={c_val}, tolerance={c_tol}, upper_limits: {c_limits[0]}, {c_limit}')
            if c_limits[1] is not None:
                c_limit = c_limits[1] - c_tol
                c_cond_l = pos[columns[c_axis]] < c_limit
                c_df = pos[c_cond_a & c_cond_l]
                if not c_df.empty:
                    c_val = round((pos.loc[c_cond_a & c_cond_l, columns[c_axis]]).iloc[0], 4)
                    if c_val < c_limit:
                        l_violation = True
                    if l_violation and uargs.verbose:
                        print(f'llimit violation posid={c_posid}, axis={c_axis}, position={c_val}, tolerance={c_tol}, lower_limits: {c_limit}, {c_limits[1]}')
            if u_violation or l_violation:
                violating_pos |= set([c_posid])
    if not violating_pos:
        return False
    return violating_pos

if uargs.prepark:
    cs.park_and_measure('all', test_tp=True)
else:
    cs.fvc_measure(test_tp=True)

for i in range(uargs.iterations):
    out_of_limits = check_if_out_of_limits()
    if not out_of_limits:
        break
    start = cs.ptlm.get_positions(return_coord=command, drop_devid=False)
    selected = start[start['DEVICE_ID'].isin(out_of_limits)].drop(columns='FLAGS')
#   print('selected\n',selected.to_string())    # DEBUG
#   print('d_pos_limits\n', "\n".join(f"{k}\t{v}" for k, v in d_pos_limits.items()))  # DEBUG
    for posid, axis_limits in d_pos_limits.items():
        cond_a = selected['DEVICE_ID'] == posid
        for axis, limits in axis_limits.items():
            tol = axis_tolerance(posid, axis)
#           print(f'axis = {axis}, column = {columns[axis]}, limits = {str(limits)}')   # DEBUG
            if set(limits) != {None}:
                if limits[0] is not None:
                    u_limit = limits[0] + tol
                    cond_b = selected[columns[axis]] > u_limit
                    above_df = selected[cond_a & cond_b]
                    if limits[1] is not None:
                        l_limit = limits[1] - tol
                        cond_c = selected[columns[axis]] < l_limit
                        below_df = selected[cond_a & cond_c]
                        # Have limits on both ends - move to middle
                        target_for_those_above = (limits[0] + limits[1])/2
                        target_for_those_below = (limits[0] + limits[1])/2
                        # Note these change selected!
#                       print(f'above_df = {str(above_df)}\n, below_df = {str(below_df)}\n')    # DEBUG
                        if not above_df.empty:
                            if uargs.verbose:
                                print(f'posid={posid}, axis={axis}, tolerance={tol}, Limits low to high: {l_limit}, {limits[1]}, {limits[0]}, {u_limit}')
                                val = round((selected.loc[cond_a & cond_b, columns[axis]]).iloc[0], 4)
                                print(f'Will move {posid} {axis} from {val} to {target_for_those_above}')
                            selected.loc[cond_a & cond_b, columns[axis]] = target_for_those_above
                        if not below_df.empty:
                            if uargs.verbose:
                                val = round((selected.loc[cond_a & cond_c, columns[axis]]).iloc[0], 4)
                                print(f'Will move {posid} {axis} from {val} to {target_for_those_below}')
                            selected.loc[cond_a & cond_c, columns[axis]] = target_for_those_below
                    else:
                        # Only have upper limit - target less than limit by angle_padding
                        target_for_those_above = limits[0] - float(uargs.angle_padding)
                        if not above_df.empty:
                            if uargs.verbose:
                                val = round((selected.loc[cond_a & cond_b, columns[axis]]).iloc[0], 4)
                                print(f'Will move {posid} {axis} from {val} to {target_for_those_above}')
                            selected.loc[cond_a & cond_b, columns[axis]] = target_for_those_above
                else:
                    # Only have lower limit (since we know both aren't None); target above it by angle_padding
                    cond_c = selected[columns[axis]] < l_limit
                    below_df = selected[cond_a & cond_c]
                    target_for_those_below = limits[1] + float(uargs.angle_padding)
                    if not below_df.empty:
                        if uargs.verbose:
                            val = round((selected.loc[cond_a & cond_c, columns[axis]]).iloc[0], 4)
                            print(f'Will move {posid} {axis} from {val} to {target_for_those_below}')
                        selected.loc[cond_a & cond_c, columns[axis]] = target_for_those_below
    selected['COMMAND'] = command
    selected['LOG_NOTE'] = 'Moving to specified limits; move_to_ranges.py'
    # Now selected is a proper move request
    cs.move_measure(selected, test_tp=True, anticollision=uargs.anticollision)

out_of_limits = check_if_out_of_limits()
if not out_of_limits:
    print(f'SUCCESS! All {len(cs.posids)} selected positioners within desired limits!')
else:
    print(f'FAILED. {len(out_of_limits)} positioners remain beyond desired limits. POSIDS: {sorted(list(out_of_limits))}')
    if uargs.disable:
        print('Disabling posids that remain out of range.')
        cs.ptlm.disable_positioners(out_of_limits, comment='move_to_range: disabled because still out of range.')
