import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-tu', '--theta_upper_limit', type=float, default=None, help='Angle which no posid should exceed on theta.')
parser.add_argument('-tl', '--theta_lower_limit', type=float, default=None, help='Angle which no posid should be below on theta.')
parser.add_argument('-pu', '--phi_upper_limit', type=float, default=None, help='Angle which no posid should exceed on phi.')
parser.add_argument('-pl', '--phi_lower_limit', type=float, default=None, help='Angle which no posid should be below on phi.')
max_iter = 10
parser.add_argument('-iter', '--iterations', type=int, default=2, help=f'Number of iterations to try to get positioners within desired rang before giving up. Must be less than {max_iter}.')
park_options = ['posintTP', 'poslocTP', 'None', 'False']
default_park = park_options[0]
parser.add_argument('-prep', '--prepark', type=str, default=default_park, help=f'str, controls initial parking move, prior to running does not count as an iteration. Valid options are: {park_options}, default is {default_park}.')
modes = ['posintTP', 'poslocTP']
parser.add_argument('-m', '--mode', type=str, default=modes[0], help=f'str, controls which tp coordinate type checks are done in. Options: {modes}. Defaults to {modes[0]}.')
default_padding = 10.0
parser.add_argument('-p', '--angle_padding', type=float, default=default_padding, help=f'Float, angle padding by which to exceed limits when only one limit is specified. Defaults to {default_padding}.')
parser.add_argument('-d', '--disable', action='store_true', help='Argue to disable positioners that remain out of range at the end of iterations.')
parser.add_argument('-a', '--anticollision', type=str, default='freeze', help='anticollision mode, can be "adjust", "adjust_requested_only", "freeze" or None. Default is "freeze"')
uargs = parser.parse_args()
if uargs.anticollision == 'None':
    uargs.anticollision = None
assert uargs.anticollision in {'adjust', 'adjust_requested_only', 'freeze', None}, f'bad argument {uargs.anticollision} for anticollision parameter'
assert uargs.iterations < max_iter, f'cannot exceed {max_iter} iterations.'
command = uargs.mode
assert command in modes, f'Invalid mode type {command}. Must be in {modes}.'
assert uargs.prepark in park_options, f'invalid park option, must be one of {park_options}'
uargs.prepark = None if uargs.prepark in ['None', 'False'] else uargs.prepark
tu = uargs.theta_upper_limit
tl = uargs.theta_lower_limit
pu = uargs.phi_upper_limit
pl = uargs.phi_lower_limit
assert {tu, tl, pu, pl} != {None}, 'Must specify at least one limit.'
if tu is not None and tl is not None:
    assert tu > tl, 'Theta upper limit must be greater than lower limit'
if pu is not None and pl is not None:
    assert pu > pl, 'Phi upper limit must be greater than lower limit'

from pecs import PECS
cs = PECS(interactive=True)

axis_limits = {'T': [tu, tl], 'P': [pu, pl]}
columns = {'T': 'X1', 'P': 'X2'}

def check_if_out_of_limits():
    all_pos = cs.ptlm.get_positions(return_coord=command, drop_devid=False)
    pos = all_pos[all_pos['DEVICE_ID'].isin(cs.posids)]
    violating_pos = set()
    for axis, limits in axis_limits.items():
        if limits[0] is not None:
            violating_pos |= set(pos[pos[columns[axis]] > limits[0]]['DEVICE_ID'])
        if limits[1] is not None:
            violating_pos |= set(pos[pos[columns[axis]] < limits[1]]['DEVICE_ID'])
    if not violating_pos:
        return False
    else:
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
    for axis, limits in axis_limits.items():
        if set(limits) != {None}:
            if limits[0] is not None:
                above_mask = selected[columns[axis]] > limits[0]
                if limits[1] is not None:
                    below_mask = selected[columns[axis]] < limits[1]
                    # Have limits on both ends - move to middle
                    target_for_those_above = (limits[0] + limits[1])/2
                    target_for_those_below = (limits[0] + limits[1])/2
                    # Note these change selected!
                    selected.loc[above_mask, columns[axis]] = target_for_those_above
                    selected.loc[below_mask, columns[axis]] = target_for_those_below
                else:
                    # Only have upper limit - target past it
                    target_for_those_above = limits[0] - uargs.angle_padding
                    selected.loc[above_mask, columns[axis]] = target_for_those_above
            else:
                # Only have lower limit (since we know both aren't None)
                below_mask = selected[columns[axis]] < limits[1]
                target_for_those_below = limits[1] + uargs.angle_padding
                selected.loc[above_mask, columns[axis]] = target_for_those_above
    selected['COMMAND'] = command
    selected['LOG_NOTE'] = 'Moving to specified limits for test setup; move_to_range.py'
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
