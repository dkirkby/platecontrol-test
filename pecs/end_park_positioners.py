#!/usr/bin/env python

# set up a log file
import argparse
import sys
parser = argparse.ArgumentParser()
def_mirror = "open"
def_ops = "observing"
parser.add_argument('-m', '--mirrorcovers', type=str, default=def_mirror, help=f'Status of mirror covers; "open" or "closed" (defaults to {def_mirror})')
parser.add_argument('-o', '--opsstate', type=str, default=def_ops, help=f'Operation state of DESI. Default is {def_ops}.')
uargs = parser.parse_args()
mirrors = uargs.mirrorcovers.lower()
opsstate = uargs.opsstate.lower()

import logging
import simple_logger
import traceback
import os
import posconstants as pc
from pecs import PECS

status = 0
update_error_report_thresh = 1.0 #warn user that an update of more than value mm was found
name = 'END_PARK_POSITIONERS'
parktype = 'intTlocP' #'poslocTP'
num_tries = 1 if parktype in ['poslocTP','posintTP'] else 1
t_angle = 0

axis_limits = {'T': [t_angle+5, t_angle-5], 'P': [155, 145]}
columns = {'T': 'X1', 'P': 'X2'}
def check_if_out_of_limits():
    all_pos = cs.ptlm.get_positions(return_coord=parktype, drop_devid=False)
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

log_dir = pc.dirs['sequence_logs']
log_timestamp = pc.filename_timestamp_str()
log_name = f'{name}_' + log_timestamp + '.log'
log_path = os.path.join(log_dir, log_name)
logger, logger_fh, logger_sh = simple_logger.start_logger(log_path)
logger_fh.setLevel(logging.DEBUG)
logger_sh.setLevel(logging.INFO)
logger.info(f'{name}: script is starting. The logging is rather verbose, but please try to follow along. A summary of important notes are provided at the end.')
#simple_logger.input2(f'Alert: you are about to run the park positioners script to park positioners at the end of the night. This takes about 2 minutes to execute. Keep power on while executing. Hit enter to continue. ')

if opsstate != 'observing':
    logger.info('DESI not in OBSERVING. Nothing to do. Script exiting.')
    sys.exit(0)

cs = PECS(interactive=False, test_name=f'FP_setup', logger=logger, inputfunc=simple_logger.input2)

logger.info(f'{name}: starting as exposure id {cs.exp.id}')

err = None #no exception

### Initial setup for proc: enable positioners then turn onf illumination and fiducials ###
### Exit script if any fails, do not continue ###

logger.info(f'{name}: turning on back illumination...')
try:
    cs.turn_on_illuminator()
    #logger.info(f'FP_SETUP: turning on back illumination returned: {ret}')
except (Exception, KeyboardInterrupt) as e:
    logger.info(f'{name}: back illumination failed to turn off with exception: {e}')
    logger.error(f'{name}: could not turn on back illumination! Please investigate before continuing!')
    sys.exit(1)

logger.info('FP_SETUP: turning on fiducials...')
try:
    cs.turn_on_fids()
    #logger.info(f'FP_SETUP: turning on fiducials returned: {ret}')
except (Exception, KeyboardInterrupt) as e:
    logger.info(f'{name}: turning on fiducials failed with exception: {e}')
    logger.error(f'{name}: could not turn on fiducials! Please investigate before continuing!')
    sys.exit(1)

### Here's the bulk of the script: 1p ###
try:
    if mirrors != 'open':
        res = cs.ptlm.petals('park_positioners', 'all', mode='normal', coords=parktype,
                            log_note='end of night park_positioners observer script', theta=t_angle)
    else:
        for i in range(num_tries):
            violating_pos = check_if_out_of_limits()
            if violating_pos:
                # Park
                logger.info(f'{name}: Parking positioners...')
                res = cs.park_and_measure(violating_pos, mode='normal', coords=parktype, log_note='end of night park_positioners observer script',
                                          match_radius=None, check_unmatched=True, test_tp=True, theta=t_angle)
                if res.empty:
                    logger.info(f'{name}: No matches in FVC image. Mirror covers may be closed. Continuing...')
                    break
            else:
                break
        out_of_limits = check_if_out_of_limits()
        if not out_of_limits:
            logger.info(f'{name}: SUCCESS! All {len(cs.posids)} selected positioners within desired limits!')
        else:
            logger.info(f'{name}: FAILED. {len(out_of_limits)} positioners remain beyond desired limits. POSIDS: {sorted(list(out_of_limits))}')

except (Exception, KeyboardInterrupt) as e:
    err = e
    status = 1
    logger.error(f'{name} crashed! See traceback below:')
    logger.critical(traceback.format_exc())

### Cleanup: turn off illuminator and fiducials, trigger fvc_collect ###
### Allow turning off to fail, observers should investigate and resolve issue before going on-sky ###
logger.info(f'{name}: turning off back illumination...')
try:
    cs.turn_off_illuminator()
    #logger.info(f'FP_SETUP: turning off back illumination returned: {ret}')
except (Exception, KeyboardInterrupt) as e:
    logger.info(f'{name}: back illumination failed to turn off with exception: {e}')
    logger.error(f'{name}: Could not turn off back illumination! Please investigate before continuing!')

logger.info(f'{name}: turning off fiducials...')
try:
    cs.turn_off_fids()
    #logger.info(f'FP_SETUP: turning off fiducials returned: {ret}')
except (Exception, KeyboardInterrupt) as e:
    logger.info(f'{name}: turning off fiducials failed with exception: {e}')
    logger.error(f'{name}: could not turn off fiducials! Please investigate before continuing!')

cs.fvc_collect()

### Print out summary for night log if successful ###
if err is None:
    logger.info(f'{name}: Successfully completed!')
else:
    logger.error(f'{name}: failed to complete. Please wait a moment to try again or contact an FP expert.')
### Clean up logger ###
logger.info(f'Log file: {log_path}')
simple_logger.clear_logger()
sys.exit(status)
