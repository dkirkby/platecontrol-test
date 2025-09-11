# -*- coding: utf-8 -*-
"""
Measure positioners, and for any in the ambiguous zone near the theta hardstop,
do successive movements to disambiguate which side of the hardstop they are on.
This script is capable of striking hardstops on occasion, in the process. C.f.
DESI-5911.
"""

import os
script_name = os.path.basename(__file__)

# set up logger
import simple_logger
try:
    import posconstants as pc
except:
    import sys
    path_to_petal = '../petal'
    sys.path.append(os.path.abspath(path_to_petal))
    print('Couldn\'t find posconstants the usual way, resorting to sys.path.append')
    import posconstants as pc
# other imports
import random
import pandas
    
# common definitions
pos_settings_keys = ['ONLY_CREEP', 'CREEP_PERIOD']

class disambig_class():
    '''
    KF - 20210212
        Not the best but putting this in a class so it can accept an outside PECS or logger instance.
    '''
    def __init__(self, pecs=None, logger=None, num_meas=1, match_radius=None, check_unmatched=False, num_tries=4, only_creep=True, creep_period=0):
        if logger is None:
            # set up logger
            import simple_logger
            log_dir = pc.dirs['calib_logs']
            log_timestamp = pc.filename_timestamp_str()
            log_name = log_timestamp + '_disambig_theta.log'
            log_path = os.path.join(log_dir, log_name)
            self.logger, _, _ = simple_logger.start_logger(log_path)
        else:
            self.logger = logger
        self.logger.info(f'Running {script_name} to ascertain true theta near hardstops.')

        if pecs is None:
            # set up pecs
            from pecs import PECS
            self.pecs = PECS(interactive=True, logger=self.logger)
            self.logger.info(f'PECS initialized, discovered PC ids {self.pecs.pcids}')
        else:
            self.pecs = pecs
        self._tp_tol_old = self.pecs.tp_tol
        self._tp_frac_old = self.pecs.tp_frac

        self.num_tries = num_tries
        self.match_radius = match_radius
        self.check_unmatched = check_unmatched
        self.num_meas = num_meas
        self.settings = {'ONLY_CREEP': only_creep, 'CREEP_PERIOD': creep_period}
        # some boilerplate
        self.common_move_meas_kwargs = {'match_radius': self.match_radius,
                                        'check_unmatched': self.check_unmatched,
                                        'test_tp': True,
                                        'num_meas': self.num_meas,
                                        }

        # global collections
        self.enabled_posids = set(self.pecs.get_enabled_posids('all', include_posinfo=False))
        self.allowed_to_fix = set(self.pecs.get_enabled_posids('sub', include_posinfo=False))
        self.unambig = set()  # stores all posids that have been resolved
        self.neighbor_data = self.pecs.ptlm.app_get("collider.pos_neighbors")
        self.neighbors = {}
        for these in self.neighbor_data.values():
            self.neighbors.update(these)
   
        # gather initial settings
        pos_settings_by_petal = self.pecs.ptlm.batch_get_posfid_val(uniqueids=self.allowed_to_fix, keys=pos_settings_keys)
        self.orig_pos_settings = {}
        for dictionary in pos_settings_by_petal.values():
            self.orig_pos_settings.update(dictionary)

    # algorithm
    def disambig(self, n_try=None):
        '''See DESI-5911 for diagram and plain-language explanation of the basic algorithm.
        
        INPUTS:  
            n_try ... Integer, indicates number of times to repeat. Repetition is done by
                      recursing this function. If n_try == 0 --> break. The direction
                      of trial moves is flipped with each repeat, so one would almost always
                      want the initial call to have n_try >= 2. Even values do the first move
                      away from the presumed hardstop, odd values go toward it. The idea of
                      allowing more than 2 repeats is to deal with cases where a robot needs
                      its neighbors to be disambiguated before it can be safely moved itself.
        
        OUTPUTS:
            Set of posids that are still in ambiguous zone. Disabled positioners are
            *excluded* from this set, regardless of what zone they are in.
        '''
        if n_try is None:
            n_try = self.num_tries
        self.logger.info(f'Starting disambiguation iteration {self.num_tries - n_try + 1} of {self.num_tries}')

        self.pecs.tp_tol = 0.0   # correct POS_T, POS_P for FVC measurements above this err value
        self.pecs.tp_frac = 1.0  # when correcting POS_T, POS_P, do so by this fraction of err distance
        # get ambiguous and unambiguous posids
        #global unambig #replaced by class variable
        ambig_dict = self.pecs.quick_query('in_theta_hardstop_ambiguous_zone', posids=self.enabled_posids)
        in_ambig_zone = {posid for posid, val in ambig_dict.items() if val == True}
        all_ambig = in_ambig_zone - self.unambig  # because previous parking moves may have put already-resolved pos into ambig theta territory
        self.unambig |= self.enabled_posids - all_ambig
        self.logger.info(f'{len(all_ambig)} enabled positioner(s) are in theta hardstop ambiguous zone: {all_ambig}')
        self.logger.info(f'{len(self.unambig)} enabled positioner(s) are unambiguous')
        do_not_fix = all_ambig - self.allowed_to_fix
        ambig = all_ambig & self.allowed_to_fix
        if do_not_fix:
            self.logger.info(f'{len(do_not_fix)} positioner(s) are not within the allowed-to-fix selection group. Excluded posids: {do_not_fix}')
        if not ambig or n_try == 0:
            self.pecs.tp_tol = self._tp_tol_old
            self.pecs.tp_frac = self._tp_frac_old
            return ambig
        self.logger.info(f'Disambiguation attempt {self.num_tries - n_try + 1} of {self.num_tries}')
        self.logger.info(f'Will attempt to resolve {len(ambig)} posid(s): {ambig}')
        
        # targets for unambiguous pos
        active_posids = ambig | self.unambig
        locT_current = self.pecs.quick_query(key='poslocT', posids=active_posids)
        locT_targets = {posid: locT_current[posid] for posid in self.unambig}
        locP_target = 150.0
        
        # where possible, target unambiguous posids "opposite" ambiguous neighbors,
        # to maximize clearance
        for posid in self.unambig:
            these_neighbors = self.neighbors[posid]
            these_ambig = {n for n in these_neighbors if n in ambig}
            if these_ambig:
                selected = random.choice(list(these_ambig))
                locT_targets[posid] = locT_current[selected]  # this is a good config to minimize collision opportunity
            
        # move unambiguous positioners to targets
        sorted_unambig = sorted(self.unambig)
        anticollision = 'adjust_requested_only'
        request_data = {'DEVICE_ID': sorted_unambig,
                        'COMMAND': 'poslocTP',
                        'X1': [locT_targets[posid] for posid in sorted_unambig],
                        'X2': locP_target,
                        'LOG_NOTE': pc.join_notes(script_name, 'parking move for unambiguous positioner'),
                        }
        request = pandas.DataFrame(request_data)
        if request.empty:
            self.logger.info('No unambiguous + enabled positioners detected.')
        else:
            self.logger.info(f'Doing parking move for {len(request)} unambiguous positioners. Anticollision mode: {anticollision}')
            self.pecs.move_measure(request=request, anticollision=anticollision, **self.common_move_meas_kwargs)

        # retract ambiguous positioners' phi axes, again to maximize clearance
        # (this is safely done by *not* allowing extra theta anticollision moves)
        anticollision = 'freeze'
        intP_current = self.pecs.quick_query(key='posintP', posids=ambig)
        intP_target = 150.0  # this also helps *extend* fully-retracted arms a little --> improves the theta measurement
        dP = {posid: intP_target - intP_current[posid] for posid in ambig}
        dtdp_requests = {posid: {'target': [0.0, dP[posid]],
                                 'log_note': pc.join_notes(script_name, 'retraction move for ambiguous positioner'),
                                 } for posid in ambig
                         }
        self.logger.info(f'Doing retraction move for {len(dtdp_requests)} ambiguous positioners. Anticollision mode: {anticollision}')
        self.pecs.ptlm.request_direct_dtdp(dtdp_requests, return_posids_only=True)
        self.pecs.ptlm.schedule_moves(anticollision=anticollision)
        self.pecs.move_measure(request=None, anticollision=anticollision, **self.common_move_meas_kwargs)
        
        # theta test moves on ambiguous positioners
        anticollision = 'freeze'
        old_settings = {}
        new_settings = {}
        settings_note = ''
        for key, skipval in {'ONLY_CREEP': None, 'CREEP_PERIOD': 0}.items():
            uarg = self.settings[key] #getattr(uargs, key.lower())
            if uarg != skipval:
                for posid in ambig:
                    if posid not in new_settings:  # and by implication, not in old_settings yet, either
                        new_settings[posid] = {}
                        old_settings[posid] = {}
                    new_settings[posid][key] = uarg  # don't worry about whether this is actually different from old value, the batch setter function a few lines below handles this more generally
                    old_settings[posid][key] = self.orig_pos_settings[posid][key]
                settings_note = pc.join_notes(settings_note, f'{key}={uarg}')
        intT_current = self.pecs.quick_query(key='posintT', posids=ambig)
        ambig_max = self.pecs.quick_query(key='max_theta_hardstop_ambiguous_zone', posids=ambig)
        ambig_min = self.pecs.quick_query(key='min_theta_hardstop_ambiguous_zone', posids=ambig)
        dT_abs = {posid: ambig_max[posid] - ambig_min[posid] + pc.theta_hardstop_ambig_exit_margin for posid in ambig}
        presumed_no_hardstop_dir = {posid: 1 if intT_current[posid] < 0 else -1 for posid in ambig}
        move_dir = {posid: presumed_no_hardstop_dir[posid] * (-1 if n_try % 2 else 1) for posid in ambig}
        dT = {posid: dT_abs[posid] * move_dir[posid] for posid in ambig}
        dir_note = {posid: f'{"away from" if presumed_no_hardstop_dir[posid] == move_dir[posid] else "toward"} currently-presumed closest hardstop' for posid in ambig}
        dtdp_requests = {posid: {'target': [dT[posid], 0.0],
                                 'log_note': pc.join_notes(script_name, 'theta test move on ambiguous positioner', dir_note[posid], settings_note),
                                 } for posid in ambig
                         }
        if any(new_settings):
            self.logger.info(f'Applying pos settings: {new_settings}')
            self.pecs.ptlm.batch_set_posfid_val(settings=new_settings, check_existing=True)
        self.logger.info(f'Doing theta test move for {len(dtdp_requests)} ambiguous positioners. Anticollision mode: {anticollision}')
        self.pecs.ptlm.request_direct_dtdp(dtdp_requests, return_posids_only=True)
        self.pecs.ptlm.schedule_moves(anticollision=anticollision)
        self.pecs.move_measure(request=None, anticollision=anticollision, **self.common_move_meas_kwargs)
        if any(new_settings):
            self.logger.info(f'Restoring pos settings: {old_settings}')
            self.pecs.ptlm.batch_set_posfid_val(settings=old_settings, check_existing=True)
        return self.disambig(n_try=n_try - 1)

if __name__ == '__main__':
    # command line argument parsing
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    max_tries = 10
    parser.add_argument('-nt', '--num_tries', type=int, default=4, help=f'int, max number of tries of the algorithm (max is {max_tries})')
    max_fvc_iter = 10
    parser.add_argument('-nm', '--num_meas', type=int, default=1, help=f'int, number of measurements by the FVC per move (max is {max_fvc_iter})')
    parser.add_argument('-r', '--match_radius', type=int, default=None, help='int, specify a particular match radius')
    parser.add_argument('-u', '--check_unmatched', action='store_true', help='turns on auto-disabling of unmatched positioners')
    parser.add_argument('-oc', '--only_creep', type=str, default='True', help='True --> move ambigous positioners slowly when going toward their possible hard limits, False --> cruise speed, None --> use existing individual positioner values')
    #20220604 - default creep period for this script was 1 (faster) compared to normal 2, this functionality didn't work
    parser.add_argument('-cp', '--creep_period', type=int, default=0, help='int, overrides positioners\' values for parameter CREEP_PERIOD (1 is fastest possible, as of 2021-09-05 a setting of 2 is typical during observations, enter 0 to use existing individual positioner values')
    uargs = parser.parse_args()

    # input validation
    assert 1 <= uargs.num_tries <= max_tries, f'out of range argument {uargs.num_tries} for num_tries parameter'
    assert 1 <= uargs.num_meas <= max_fvc_iter, f'out of range argument {uargs.num_meas} for num_meas parameter'
    assert (2 == uargs.creep_period) or (0 == uargs.creep_period), 'Creep period different from 2 (default) not supported!' #Note 0 is no change
    # remove the above assertion when/if functionality to infrom petalcontroller of creep period (and clean up if crashes) is added
    assert 0 <= uargs.creep_period <= 2, f'out of range argument {uargs.creep_period} for creep_period parameter'
    if pc.is_none(uargs.only_creep):
        uargs.only_creep = None
    elif pc.is_boolean(uargs.only_creep):
        uargs.only_creep = pc.boolean(uargs.only_creep)
    else:
        assert False, f'out of range argument {uargs.only_creep} for only_creep parameter'

    disambig_obj = disambig_class(pecs=None, logger=None, num_meas=uargs.num_meas, num_tries=uargs.num_tries,
                                  match_radius=uargs.match_radius, check_unmatched=uargs.check_unmatched,
                                  only_creep=uargs.only_creep, creep_period=uargs.creep_period)
    ambig = disambig_obj.disambig()
    disambig_obj.logger.info('Disambiguation loops complete.')
    if ambig:
        details = disambig_obj.pecs.ptlm.quick_table(posids=ambig, coords=['posintTP', 'poslocTP'], as_table=False, sort='POSID')
        if pc.is_string(details):
            details_str = details
        else:
            details_str = ''
            for petal_id, table_str in details.items():
                details_str += f'\n{petal_id}\n{table_str}\n'
        disambig_obj.logger.warning(f'{len(ambig)} positioners remain *unresolved*. Details:\n{details_str}')
    else:
        disambig_obj.logger.info('All selected ambiguous cases were resolved!')

