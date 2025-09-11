#!/usr/bin/env python

# set up a log file
import logging
import simple_logger
import traceback
import os
import sys
from datetime import datetime
import posconstants as pc
from pecs import PECS
import DOSlib.flags as flags
from DOSlib.util import obs_day
from DOSlib.join_instance import join_instance
from argparse import ArgumentParser

parser = ArgumentParser(description="FP_SETUP: nightly setup for the DESI focal plane.")
parser.add_argument("-i", "--instance",type=str, help = 'Instance name (desi_<obsday> is default>')
args = parser.parse_args()

inst = args.instance
# reset Pyro variables
os.environ['PYRO_NS_HOST'] = ''
os.environ['PYRO_NS_PORT'] = ''
os.environ['PYRO_NS_BCPORT'] = ''
if not isinstance(inst, str):
    inst = f'desi_{obs_day()}'
print(f'connecting to instance {inst}')
try:
    join_instance(inst, must_be_running=True)
except Exception as e:
    print(f'FP_SETUP: Exception joining instance {inst}: {str(e)}')
    sys.exit(1) 

update_error_report_thresh = 1.0 #warn user that an update of more than value mm was found

log_dir = pc.dirs['sequence_logs']
log_timestamp = pc.filename_timestamp_str()
log_name = 'FP_setup_' + log_timestamp + '.log'
log_path = os.path.join(log_dir, log_name)
logger, logger_fh, logger_sh = simple_logger.start_logger(log_path)
logger_fh.setLevel(logging.DEBUG)
logger_sh.setLevel(logging.INFO)
logger.info('FP_SETUP: script is starting. The logging is rather verbose, but please try to follow along. A summary of important notes are provided at the end.')
simple_logger.input2(f'Alert: you are about to run the focalplane setup script to recover positioners or prepare the focalplane. This takes *up to half an hour* to execute. Hit enter to continue. ')

cs = PECS(interactive=False, test_name=f'FP_setup', logger=logger, inputfunc=simple_logger.input2)
cs.ptlm.record_script_usage(script='fp_setup', alarm_id=1801, message='FP_SETUP starting...')
logger.info(f'FP_SETUP: starting as exposure id {cs.exp.id}')

if datetime.today().weekday() == 6: #this is sunday
    cs.fvc_feedback_timeout = 120.0 #2 minutes!

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

logger.info('FP_SETUP: turning on back illumination...')
try:
    cs.turn_on_illuminator()
    #logger.info(f'FP_SETUP: turning on back illumination returned: {ret}')
except (Exception, KeyboardInterrupt) as e:
    logger.info(f'FP_SETUP: back illumination failed to turn off with exception: {e}')
    logger.error('FP_SETUP: could not turn on back illumination! Please investigate before continuing!')
    sys.exit(1)

logger.info('FP_SETUP: turning on fiducials...')
try:
    cs.turn_on_fids()
    #logger.info(f'FP_SETUP: turning on fiducials returned: {ret}')
except (Exception, KeyboardInterrupt) as e:
    logger.info(f'FP_SETUP: turning on fiducials failed with exception: {e}')
    logger.error('FP_SETUP: could not turn on fiducials! Please investigate before continuing!')
    sys.exit(1)

logger.info('FP_SETUP: caching keepouts...')
try:
    keepouts = cs.ptlm.cache_keepouts()
except (Exception, KeyboardInterrupt) as e:
    logger.info(f'FP_SETUP: caching keepouts failed with exception: {e}')
    logger.error('FP_SETUP: could not cache keepouts! Please investigate before continuing!')
    sys.exit(1)

# Do enable last since we want to re-disable if we fail
logger.info('FP_SETUP: enabling positioners...')
ret = cs.ptlm.enable_positioners(ids='all', comment='FP setup script initial enable')
logger.info(f'FP_SETUP: enable_positioners returned: {ret}')
if not(ret == 'SUCCESS' or ret is None):
    logger.error('FP_SETUP: failed to enable positioners! Exiting and try again, if failed twice contact FP expert!')
    ret = cs.ptlm.disable_positioners(ids=initial_disabled, comment='FP_SETUP - crashed in execution, resetting disabled devices')
    logger.info(f'FP_SETUP: disable_positioners returned: {ret}')
    sys.exit(1)

# Re-setup petal/posids in pecs to reflect newly enabled positioners
cs.ptl_setup(cs.pcids)

### Here's the bulk of the setup: 1p, disambiguate, park overlaps, 1p ###
try:
    logger.info('FP_SETUP: running 1p calibration...')
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
            logger.warning(f'FP_SETUP: found {len(selection)} positioners with updates greater than {update_error_report_thresh}!')
            logger.warning(f'FP_SETUP: update deatils: \n{selection}')
    disabled_during_1p = enabled_before_1p - enabled_after_1p
    logger.info(f'FP_SETUP: {len(disabled_during_1p)} disabled during initial 1p calib. Posids {disabled_during_1p}')

    logger.info('FP_SETUP: shrinking keepouts for overlapping positioners')
    overlapping = cs.ptlm.get_overlaps(as_dict=True)
    overlapping_posids = set()
    for petal, overlaps in overlapping.items():
        posids = set(overlaps.keys())
        overlapping_posids |= posids
        if posids:
            logger.info(f'FP_SETUP: {petal} has {len(posids)} overlapping positioners! List: {posids}')
            logger.info(f'FP_SETUP: shrinking keepouts for above positioners...')
            cs.ptlm.set_keepouts(posids=posids, angT=-7.0, radT=-0.5, angP=-10.0, radP=-0.5, participating_petals=petal)

    logger.info('FP_SETUP: running disambiguation loops...')
    enabled_before_disambig = get_pos_set('enabled')
    disambig_obj = disambig_class(pecs=cs, logger=logger, num_meas=1, check_unmatched=True, num_tries=4)
    ambig = disambig_obj.disambig()
    logger.info('FP_SETUP: Disambiguation loops complete.')
    if ambig:
        details = cs.ptlm.quick_table(posids=ambig, coords=['posintTP', 'poslocTP'], as_table=False, sort='POSID')
        if pc.is_string(details):
            details_str = details
        else:
            details_str = ''
            for petal_id, table_str in details.items():
                details_str += f'\n{petal_id}\n{table_str}\n'
        logger.warning(f'{len(ambig)} positioners remain *unresolved* and will be disabled. Details:\n{details_str}')
        ret = cs.ptlm.disable_positioners(ids=ambig, comment='FP setup - disabling positioners that remain in ambiguous theta range.',
                                          set_flag=flags.POSITIONER_FLAGS_MASKS['BADPERFORMANCE'])
        logger.info(f'FP_SETUP: disable_positioners returned: {ret}')
    else:
        logger.info('FP_SETUP: All selected ambiguous cases were resolved!')

    enabled_after_disambig = get_pos_set('enabled')
    disabled_during_disambig = enabled_before_disambig - enabled_after_disambig
    logger.info(f'FP_SETUP: {len(disabled_during_disambig)} disabled during disambiguation loops. Posids {disabled_during_disambig}')

    #Park overlapping positioners to guarentee we tried to move them with shrunk keepouts
    #Do AFTER disambiguation in case overlapping positioner in ambiguous zone!!!
    if overlapping_posids:
        logger.info(f'FP_SETUP: tucking {len(overlapping_posids)} overlapping positioners.')
        #easiest to tuck by telling them to park
        cs.park_and_measure(overlapping_posids, coords='intTlocP', log_note='FP_Setup tucking overlapping positioners')

    # 1p calibration to clean up and tell us how we did
    logger.info(f'FP_SETUP: running final 1p calibration...')
    enabled_before_1p2 = get_pos_set('enabled')
    updates = onepoint(cs, mode='posTP', move=False, commit=True, tp_tol=0.065,
                       tp_frac=1.0, num_meas=1, use_disabled=True, check_unmatched=True)
    enabled_after_1p2 = get_pos_set('enabled')
    updates = updates[['DEVICE_ID', 'DEVICE_LOC', 'PETAL_LOC', 'ERR_XY']]
    updates.sort_values('ERR_XY', inplace=True, ascending=False)
    logger.debug('Updates sorted by descending ERR_XY:')
    logger.debug(updates.to_string())
    #logger.info(updates.to_string(max_rows=100))
    if update_error_report_thresh:
        selection = updates.loc[updates['ERR_XY'] > update_error_report_thresh]
        if len(selection):
            logger.warning(f'FP_SETUP: found {len(selection)} positioners with updates greater than {update_error_report_thresh}!')
            logger.warning(f'FP_SETUP: update deatils: \n{selection}')
    disabled_during_1p2 = enabled_before_1p2 - enabled_after_1p2
    logger.info(f'FP_SETUP: {len(disabled_during_1p2)} disabled during initial 1p calib. Posids {disabled_during_1p2}')
    final_enabled = get_pos_set('enabled')
except (Exception, KeyboardInterrupt) as e:
    err = e
    logger.error('FP_SETUP crashed! See traceback below:')
    logger.critical(traceback.format_exc())
    logger.info('Attempting to preform cleanup before hard crashing. Configure the instance before trying again.')
    try:
        logger.info('FP_SETUP: Re-disabling initially disabled positioners for safety...')
        ret = cs.ptlm.disable_positioners(ids=initial_disabled, comment='FP_SETUP - crashed in execution, resetting disabled devices')
        logger.info(f'FP_SETUP: disable_positioners returned: {ret}')
    except:
        logger.info(f'FP_SETUP: disable_positioners returned: {ret}')
        logger.critical('FP_SETUP: failed to return to initial state. DO NOT continue without consulting FP expert.')

### Cleanup: restore keepouts, turn off illuminator and fiducials, trigger fvc_collect ###
restore_keepout_err = False
logger.info('FP_SETUP: restoring positioner keepout zones...')
try:
    for petal, path in keepouts.items():
        cs.ptlm.restore_keepouts(path=path, participating_petals=petal)
except (Exception, KeyboardInterrupt) as e:
    logger.info(f'FP_SETUP: keepouts cache to restore by hand {keepouts}')
    logger.error(f'FP_SETUP: failed to restore keepouts, exception: {e}')
### Allow turning off to fail, observers should investigate and resolve issue before going on-sky ###
logger.info('FP_SETUP: turning off back illumination...')
try:
    cs.turn_off_illuminator()
    #logger.info(f'FP_SETUP: turning off back illumination returned: {ret}')
except (Exception, KeyboardInterrupt) as e:
    logger.info(f'FP_SETUP: back illumination failed to turn off with exception: {e}')
    logger.error('FP_SETUP: Could not turn off back illumination! Please investigate before continuing!')

logger.info('FP_SETUP: turning off fiducials...')
try:
    cs.turn_off_fids()
    #logger.info(f'FP_SETUP: turning off fiducials returned: {ret}')
except (Exception, KeyboardInterrupt) as e:
    logger.info(f'FP_SETUP: turning off fiducials failed with exception: {e}')
    logger.error('FP_SETUP: could not turn off fiducials! Please investigate before continuing!')

cs.fvc_collect()

### Print out summary for night log if successful ###
if err is None:
    logger.info('FP_SETUP: Successfully completed! Please record any following notes in the night log in addition to any other observations:')
    newly_enabled = final_enabled - initial_enabled
    if newly_enabled:
        logger.info(f'FP_SETUP NOTE: recovered {len(newly_enabled)} positioner robots.')
    if ambig:
        logger.info(f'FP_SETUP NOTE: {len(ambig)} positioners remaining in ambiguous theta range and are not used tonight: posids {ambig}')
    non_ambig_disabled = disabled_during_1p | disabled_during_1p2 | (disabled_during_disambig-ambig)
    if non_ambig_disabled:
        logger.info(f'FP_SETUP_NOTE: {len(non_ambig_disabled)} positioners were disabled during the course of this script for other reasons, likely due to match errors or communication errors.')
        logger.info(f'FP_SETUP_NOTE: other disabled posids: {non_ambig_disabled}')
    cs.ptlm.raise_script_alarm(script='fp_setup', alarm_id=1801, level='event', message='FP_SETUP completed successfully!')
else:
    logger.error('FP_SETUP: focal plane setup FAILED to complete! Please wait a moment before trying again or contact an FP expert.')
    cs.ptlm.raise_script_alarm(script='fp_setup', alarm_id=1801, level='error', message='FP_SETUP failed to complete!')
if restore_keepout_err:
    logger.error('FP_SETUP: could not restore keepouts!!!! DO NOT continue until this is resolved. Contact an expert.')
### Clean up logger ###
logger.info(f'Log file: {log_path}')
simple_logger.clear_logger()

