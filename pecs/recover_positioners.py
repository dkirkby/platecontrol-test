#!/usr/bin/env python

# set up a log file
import logging
import simple_logger
import traceback
import os
import posconstants as pc
from pecs import PECS

update_error_report_thresh = 1.0 #warn user that an update of more than value mm was found
name = 'RECOVER_POSITIONERS'

log_dir = pc.dirs['sequence_logs']
log_timestamp = pc.filename_timestamp_str()
log_name = f'{name}_' + log_timestamp + '.log'
log_path = os.path.join(log_dir, log_name)
logger, logger_fh, logger_sh = simple_logger.start_logger(log_path)
logger_fh.setLevel(logging.DEBUG)
logger_sh.setLevel(logging.INFO)
logger.info(f'{name}: script is starting. The logging is rather verbose, but please try to follow along. A summary of important notes are provided at the end.')
simple_logger.input2(f'Alert: you are about to run the recover_positioners script to enable positioners disabled by communication issues. This takes about 2 minutes to execute. Hit enter to continue. ')

cs = PECS(interactive=False, test_name=f'FP_setup', logger=logger, inputfunc=simple_logger.input2)

logger.info(f'{name}: starting as exposure id {cs.exp.id}')

#from 1p_calib import onepoint # doesn't work
import importlib
onept = importlib.import_module('1p_calib')
onepoint = onept.onepoint
from disambiguate_theta import disambig_class

def get_pos_set(which='enabled'):
    pos_set = set()
    if which == 'enabled':
        pos_dict = cs.ptlm.all_enabled_posids()
    elif which == 'disabled':
        pos_dict = cs.ptlm.all_disabled_posids()
    for poslist in pos_dict.values():
        pos_set |= set(poslist)
    return pos_set

initial_disabled = get_pos_set('disabled')
initial_enabled = get_pos_set('enabled')

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
    import sys; sys.exit(1)

logger.info('FP_SETUP: turning on fiducials...')
try:
    cs.turn_on_fids()
    #logger.info(f'FP_SETUP: turning on fiducials returned: {ret}')
except (Exception, KeyboardInterrupt) as e:
    logger.info(f'{name}: turning on fiducials failed with exception: {e}')
    logger.error(f'{name}: could not turn on fiducials! Please investigate before continuing!')
    import sys; sys.exit(1)

# Do enable last since we want to re-disable if we fail
logger.info(f'{name}: enabling positioners...')
ret = cs.ptlm.careful_enable_positioners('all', 'COMM', comment=f'{name}: script initial enable')
logger.info(f'{name}: careful_enable_positioners returned: {ret}')
for i in ret.values():
    if not(isinstance(i, set)):
        logger.error(f'{name}: failed to enable positioners! Exiting and try again, if failed twice contact FP expert!')
        ret = cs.ptlm.disable_positioners(ids=initial_disabled, comment=f'{name}: crashed in execution, resetting disabled devices')
        logger.info(f'{name}: disable_positioners returned: {ret}')
        import sys; sys.exit(1)

# Re-setup petal/posids in pecs to reflect newly enabled positioners
cs.ptl_setup(cs.pcids)

### Here's the bulk of the script: 1p ###
try:
    logger.info(f'{name}: running 1p calibration...')
    enabled_before_1p = get_pos_set('enabled')
    updates = onepoint(cs, mode='posTP', move=False, commit=True, tp_tol=0.065,
                       tp_frac=1.0, num_meas=1, use_disabled=True, check_unmatched=True)
    enabled_after_1p = get_pos_set('enabled')
    updates = updates[['DEVICE_ID', 'DEVICE_LOC', 'PETAL_LOC', 'ERR_XY']]
    updates.sort_values('ERR_XY', inplace=True, ascending=False)
    logger.debug('Updates sorted by descending ERR_XY:')
    logger.debug(updates.to_string())
    #logger.info(updates.to_string(max_rows=100))
    if update_error_report_thresh:
        selection = updates.loc[updates['ERR_XY'] > update_error_report_thresh]
        if len(selection):
            logger.warning(f'{name}: found {len(selection)} positioners with updates greater than {update_error_report_thresh}!')
            logger.warning(f'{name}: update deatils: \n{selection}')
    disabled_during_1p = enabled_before_1p - enabled_after_1p
    logger.info(f'{name}: {len(disabled_during_1p)} disabled during initial 1p calib. Posids {disabled_during_1p}')
    final_enabled = get_pos_set('enabled')
except (Exception, KeyboardInterrupt) as e:
    err = e
    logger.error(f'{name} crashed! See traceback below:')
    logger.critical(traceback.format_exc())
    logger.info('Attempting to preform cleanup before hard crashing. Configure the instance before trying again.')
    try:
        logger.info(f'{name}: Re-disabling initially disabled positioners for safety...')
        ret = cs.ptlm.disable_positioners(ids=initial_disabled, comment='FP_SETUP - crashed in execution, resetting disabled devices')
        logger.info(f'{name}: disable_positioners returned: {ret}')
    except:
        logger.info(f'{name}: disable_positioners returned: {ret}')
        logger.critical(f'{name}: failed to return to initial state. DO NOT continue without consulting FP expert.')

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
    logger.info(f'{name}: Successfully completed! Please record any following notes in the night log in addition to any other observations:')
    newly_enabled = final_enabled - initial_enabled
    if newly_enabled:
        logger.info(f'{name}_NOTE: recovered {len(newly_enabled)} positioner robots.')
    non_ambig_disabled = disabled_during_1p
    if non_ambig_disabled:
        logger.info(f'{name}_NOTE: {len(non_ambig_disabled)} positioners were disabled during the course of this script for other reasons, likely due to match errors or communication errors.')
        logger.info(f'{name}_NOTE: other disabled posids: {non_ambig_disabled}')
else:
    logger.error(f'{name}: failed to complete. Please wait a moment to try again or contact an FP expert.')
### Clean up logger ###
logger.info(f'Log file: {log_path}')
simple_logger.clear_logger()
