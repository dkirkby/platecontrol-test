'''
1p_calib.py

Contains a function to perform a onepoint calibration on positioners.

If executed standalone, main function will create a PECS instance to
interactively run the onepoint calibration.
'''

import pandas as pd

def onepoint(pecs, mode='posintTP', move=False, commit=True, tp_tol=0.0, tp_frac=1.0,
             match_radius=None, num_meas=1, use_disabled=False, check_unmatched=False):
    '''
    Function to do a onepoint calibration with given PECS isntance.
    args:
        pecs - instance of PECS
    kwargs:
        mode - string, posTP or offsetTP
        move - bool, whether or not to move positioners to nominal park position.
        commit - bool, whether or not to automatically commit updates.
        tp_tol - float, Minimum error in mm over which to update TP.
        tp_frac - float, Percentatge of error to apply in the TP update.

    returns 
        updates - pandas.DataFrame with columns 'DEVICE_ID','DEVICE_LOC',
                  'PETAL_LOC','ERR_XY','dT','dP' and a pair of TP parameters
    '''
    pecs.ptlm.set_exposure_info(pecs.exp.id, pecs.iteration)
    note = pecs.decorate_note(log_note=f'1p_calib_{mode}')
    if move:
        pecs.ptlm.park_positioners(pecs.posids, mode='normal', log_note=f'Moving for 1p_calib_{mode}')
    _, meas, matched, _ = pecs.fvc_measure(test_tp=False, check_unmatched=check_unmatched,
                                           match_radius=match_radius, num_meas=num_meas)
    if use_disabled:
        calib_pos = meas.loc[sorted(matched)]
    else:
        calib_pos = meas.loc[sorted(matched & (set(pecs.posids)))]
    updates = pecs.ptlm.test_and_update_TP(calib_pos.reset_index(), mode=mode,
                                           auto_update=commit, tp_updates_tol=tp_tol,
                                           tp_updates_fraction=tp_frac,
                                           log_note=note, verbose=True)
    updates = pd.concat([up for up in updates.values()]).reset_index(drop=True)
    return updates

if __name__ == '__main__':
    # set up a log file
    import logging
    import simple_logger
    import traceback
    import os
    import posconstants as pc
    from pecs import PECS
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    def_mode = 'posTP'
    def_tol = 0.0
    def_frac = 1.0
    parser.add_argument('-m', '--mode', type=str, default=def_mode, help=f'Which pair or tp angles to adjust, posTP or offsetsTP. Default is {def_mode}.')
    parser.add_argument('-t', '--tp_tol', type=float, default=def_tol, help=f'Minimum error in mm over which to update TP. Default is {def_tol} mm.')
    parser.add_argument('-f', '--tp_frac', type=float, default=def_frac, help=f'Percentatge of error to apply in the TP update. Maximum 1.0. Default is {def_frac}.')
    parser.add_argument('-r', '--match_radius', type=int, default=None, help='int, specify a particular match radius, other than default in PECS config.')
    parser.add_argument('-n', '--no_update', action='store_true', help='suppress auto-updating of TP values. Defaults to False.')
    parser.add_argument('-prep', '--prepark', action='store_true', help='automatically do an initial parking move prior to the measurement. Defaults to False.')
    parser.add_argument('-ud', '--use_disabled', action='store_true', help='Use disabled positioners in calibration as well. (All measured positioners will be used.)')
    max_fvc_iter = 10
    parser.add_argument('-nm', '--num_meas', type=int, default=1, help=f'int, number of measurements by the FVC per move (default is 1, max is {max_fvc_iter})')
    parser.add_argument('-u', '--check_unmatched', action='store_true', help='Whether or not to check and disable unmatched positioners and their neighbors. Defaults to False.')
    uargs = parser.parse_args()
    assert uargs.mode in ['posTP', 'offsetsTP'], 'mode argument must be either posTP or offsetsTP!'
    assert uargs.tp_frac <= 1.0, 'Updates fraction cannot be greater than 1.0!'
    assert 1 <= uargs.num_meas <= max_fvc_iter, f'out of range argument {uargs.num_meas} for num_meas parameter'

    log_dir = pc.dirs['sequence_logs']
    log_timestamp = pc.filename_timestamp_str()
    log_name = log_timestamp + '_1p_calib.log'
    log_path = os.path.join(log_dir, log_name)
    logger, logger_fh, logger_sh = simple_logger.start_logger(log_path)
    logger_fh.setLevel(logging.DEBUG)
    logger_sh.setLevel(logging.INFO)
    simple_logger.input2(f'Running one point calibration for {uargs.mode}. Hit enter to continue. ')
    cs = PECS(interactive=True, test_name=f'1p_calib_{uargs.mode}', logger=logger, inputfunc=simple_logger.input2)
    try:
        updates = onepoint(cs, mode=uargs.mode, move=uargs.prepark, commit=not(uargs.no_update),
                           tp_tol=uargs.tp_tol, tp_frac=uargs.tp_frac, match_radius=uargs.match_radius,
                           num_meas=uargs.num_meas, use_disabled=uargs.use_disabled, check_unmatched=uargs.check_unmatched)
        updates = updates[['DEVICE_ID', 'DEVICE_LOC', 'PETAL_LOC', 'ERR_XY']]
        updates.sort_values('ERR_XY', inplace=True, ascending=False)
        logger.info('Updates sorted by descending ERR_XY:')
        logger.debug(updates.to_string())
        logger.info(updates.to_string(max_rows=100))
    except Exception as e:
        logger.error('1p_calib crashed! See traceback below:')
        logger.critical(traceback.format_exc())
        logger.info('Attempting to preform cleanup before hard crashing. Configure the instance before trying again.')
    cs.fvc_collect()
    logger.info(f'Log file: {log_path}')
    simple_logger.clear_logger()
