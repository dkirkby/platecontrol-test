'''
P.E.C.S. - Petal Engineer Control System
Kevin Fanning fanning.59@osu.edu

Simple wrapper of PETAL and FVC proxies to make mock daytools that call
commands in PetalApp in a similar manner to an OCS sequence.
FVC/Spotmatch/Platemaker handling are hidden behind the FVC proxy measure
function. See FVChandler.py for use.
Requires a running DOS Instance with one PetalApp
(multiple petals are not supported at the moment)
along with FVC, Spotmatch and Platemaker.

'''
import os
import sys
import time
import numpy as np
import pandas as pd
from configobj import ConfigObj
import posconstants as pc
from DOSlib.proxies import FVC, PetalMan, SimpleProxy  # , Illuminator
from DOSlib.exposure import Exposure
from fvc_sim import FVC_proxy_sim

# Required environment variable
PECS_CONFIG_FILE = os.environ.get('PECS_CONFIG_FILE') # corresponds to e.g. 'pecs_default.cfg' or 'pecs_lbnl.cfg'
assert PECS_CONFIG_FILE

class PECS:
    '''if there is already a PECS instance from which you want to re-use
       the proxies, simply pass in pecs.fvc and pecs.ptls
       
        All available petal role names:     self.ptlm.Petals.keys()
        All selected petals:                self.ptlm.participating_petals
        Selected PCIDs:                     self.pcids from pecs_*.cfg
        
        For different hardware setups, you need to intialize including the
        correct config file for that local setup. This will be stored in
        fp_settings/hwsetups/ with a name like pecs_default.cfg or pecs_lbnl.cfg.

        KF - 20210201: A Note on PCIDs, PCID means Petal Controller ID, mistakenly
                       used as an alternate name for petal loc which it actually refers to.
                       pcid is only used internal to the code (to be purged later), all
                       printouts now properly refer to it as petal loc.
    '''
    def __init__(self, fvc=None, ptlm=None, printfunc=print, interactive=None,
                 test_name='PECS', device_locs=None, no_expid=False, posids=None,
                logger=None, inputfunc=input):
        # Allow local config so scripts do not always have to collect roles
        # and names from the user. No check for illuminator at the moment
        # since it is not used in tests.
        self._pcid2role = lambda pcid: f'PETAL{pcid}'
        self._role2pcid = lambda role: int(role.replace('PETAL', ''))
        self.use_desimeter = False
        self.match_to_positioner_centers = False
        self.match_radius = 50
        self.exptime = 1.
        self.start_time = pc.now()
        self.test_name = test_name
        self.interactive = interactive
        self.logger = logger
        self.input = inputfunc

        self.fvc_feedback_timeout = 120.0
        self.execute_move_timeout = 120.0

        pecs_local = ConfigObj(PECS_CONFIG_FILE, unrepr=True, encoding='utf-8')
        for attr in pecs_local.keys():
            setattr(self, attr, pecs_local[attr])
        if self.exposure_dir is None:
            self.exposure_dir = '/exposures/desi/'
        self.printfunc = printfunc
        if fvc is None:  # instantiate FVC proxy, sim or real
            if 'SIM' in self.fvc_role.upper():
                self.fvc = FVC_proxy_sim(max_err=0.0001)
            else:
                from Pyro4.errors import NamingError
                try:
                    self.fvc = FVC(self.pm_instrument, fvc_role=self.fvc_role,
                                   constants_version=self.constants_version,
                                   use_desimeter=self.use_desimeter,
                                   match_to_positioner_centers=self.match_to_positioner_centers,
                                   use_arcprep=self.use_arcprep, logger=self.logger)
                except NamingError as e:
                    self.print('Naming error! You did not join the instance or you joined the incorrect one! Try running join_instance desi_yyyymmdd.')
                    raise e
            self.print(f"FVC proxy created for instrument: "
                       f"{self.fvc.get('instrument')}")
        elif not(fvc): # if FVC is False don't startup FVC stuff
            self.fvc = fvc
            self.print('Starting PECS with no FVC proxy.')
        else:
            self.fvc = fvc
            self.print(f"Reusing existing FVC proxy: "
                       f"{self.fvc.get('instrument')}")
        if ptlm is None:
            self.ptlm = PetalMan(logger=self.logger)
            pcids = [self._role2pcid(role)
                     for role in self.ptlm.participating_petals]
            # Only use petals that are availible in Petalman
            self.pcids = list(set(self.pcids) & set(pcids))
            self.ptlm.participating_petals = [self._pcid2role(p) for p in self.pcids]
            self.print(f'PetalMan proxy initialized with active petal '
                       f'role numbers (locs): {pcids}')
        else:
            self.ptlm = ptlm
            self.print(f'Reusing existing PetalMan proxy with active petals: '
                       f'{self.ptlm.participating_petals}')
        self._all_petals = list(self.ptlm.Petals.keys())
        if self.illuminated_pcids == 'all':
            self.illuminated_ptl_roles = self._all_petals
        else:
            self.illuminated_ptl_roles = [self._pcid2role(pcid)
                                          for pcid in self.illuminated_pcids
                                          if self._pcid2role(pcid)
                                          in self._all_petals]
        assert set(self.illuminated_ptl_roles) <= set(
            self.ptlm.Petals.keys()), (
            'Illuminated petals must be in availible petals!')
        if self.interactive or (self.pcids is None):
            self.interactive_ptl_setup(device_locs=device_locs, posids=posids)
            self.fid_petals = self.ptlm.participating_petals
        elif self.interactive is False:
            # In non-interactive mode, check against obs_petals in petalman
            obs_petals = self.ptlm.get('obs_petals')
            obs_locs = [self._role2pcid(role) for role in obs_petals]
            self.pcids = list(set(self.pcids) & set(obs_locs))
            self.ptl_setup(self.pcids)  # use PCIDs specified in cfg
            # Check illum_petals and fid_petals
            illum_petals = self.ptlm.get('illum_petals')
            fid_petals = self.ptlm.get('fid_petals')
            self.illuminated_ptl_roles = list(set(self.illuminated_ptl_roles) & set(illum_petals))
            self.fid_petals = list(set(fid_petals) & set([self._pcid2role(p) for p in self.pcids]))
        # Do this after interactive_ptl_setup
        if self.interactive:
            self.home_adc() #asks to home, not automatic
            self.turn_on_fids()
            self.turn_on_illuminator()
        #Setup exposure ID last incase aborted doing the above
        if not(no_expid):
            self._get_expid()

    def _get_expid(self):
        '''
        Setup entry in exposureDB to get and exposure ID. 
        '''
        self.exp = Exposure(readonly=False)
        self.exp.sequence = 'Focalplane'
        self.exp.program = self.test_name
        self.exp.exptime = self.exptime
        self.iteration = 0
        self.print(f'DESI exposure ID set up as: {self.exp.id}')

    def print(self, msg):
        if self.logger is not None:
            self.logger.info(msg)  # use broadcast logger to log to all pcids
        else:
            self.printfunc(msg)

    def _parse_yn(self, yn_str):
        #  Trying to accept varieties like y/n, yes/no
        if 'y' in yn_str.lower():
            return True
        elif 'n' in yn_str.lower():
            return False
        else:
            return self._parse_yn(
                self.input(f'Invalid input: {yn_str}, must be y/n:'))

    def ptl_setup(self, pcids, posids=None, illumination_check=False, device_locs=None):
        '''input pcids must be a list of integers'''
        self.print(f'Setting up petals and positioners for {len(pcids)} '
                   f'selected petals, locs: {pcids}')
        if illumination_check:
            for pcid in pcids:  # illumination check
                assert self._pcid2role(pcid) in self.illuminated_ptl_roles, (
                    f'PETAL{pcid} must be illuminated.')
        self.ptlm.participating_petals = [self._pcid2role(pcid) for pcid in pcids]
        all_posinfo_dicts = self.ptlm.get_positioners(enabled_only=False)
        self.all_posinfo = pd.concat(all_posinfo_dicts.values(), ignore_index=True)
        self.petal_locs = self.all_posinfo[['DEVICE_ID', 'PETAL_LOC']]
        if posids is None:
            posids0, posinfo = self.get_enabled_posids(posids='all', include_posinfo=True)
            if device_locs:
                self.print(f'Looking for enabled positioners at {len(device_locs)} specific' +
                           ' device locations on each petal')
                drop_idxs = set()
                for idx, row in posinfo.iterrows():
                    if row['DEVICE_LOC'] not in device_locs:
                        drop_idxs.add(idx)
                posinfo = posinfo.drop(drop_idxs, axis=0)
                if posinfo.index.name == 'DEVICE_ID':
                    posids0 = sorted(posinfo.index)
                else:
                    posids0 = sorted(posinfo['DEVICE_ID'])
                self.print(f'Found {len(posids0)} enabled posids matching specified device locations')
            else:
                self.print('Defaulting to all enabled positioners')
        else:
            ret = self.ptlm.get_positioners(posids=posids, enabled_only=False)
            posinfo = pd.concat(list(ret.values())).set_index('DEVICE_ID')
            posids0 = sorted(set(posids) & set(posinfo.index))  # double check
            self.print(f'Validated {len(posids0)} of {len(posids)} positioners specified')
        self.posids = posids0
        self.posinfo = posinfo
        self.ptl_roles = self.ptlm.participating_petals

    def interactive_ptl_setup(self, device_locs=None, posids=None):
        self.print('Running interactive setup for PECS')
        pcids = self._interactively_get_pcid()  # set selected ptlid
        if posids:
            device_locs = None  # i.e. posids override device_locs
        elif device_locs == None:
            posids = self._interactively_get_posids()  # set selected posids
        self.ptl_setup(pcids, posids=posids, device_locs=device_locs)

    def pcid_lookup(self, posid):
        if posid in self.posinfo.index:
            return self.posinfo.loc[posid, 'PETAL_LOC']
        else:
            raise ValueError(f'Invalid posid {posid} not found in PCIDs: '
                             f'{self.pcids}')

    def ptl_role_lookup(self, posid):
        return self._pcid2role(self.pcid_lookup(posid))

    def _interactively_get_pcid(self):
        pcids = self.input('Please enter integer petal locs seperated by spaces. '
                           'Leave blank to select petals specified in cfg: ')
        if pcids == '':
            pcids = self.pcids
        for pcid in pcids:  # validate pcids against petalman available roles
            assert f'PETAL{pcid}' in self.ptlm.Petals.keys(), (
                f'PETAL{pcid} unavailable, all available ICS petal roles: '
                f'{self.ptlm.Petals.keys()}')
        self.print(f'Selected {len(self.pcids)} petals, Petal Locs: {self.pcids}')
        self.ptlm.participating_petals = [self._pcid2role(p) for p in pcids]
        return pcids

    def _interactively_get_posids(self):
        user_text = self.input('Please list canids (can??) or posids, seperated by '
                               'spaces. Leave blank to select all positioners: ')
        enabled_only, kwarg = True, {}
        if user_text == '':
            self.print('Defaulting to all enabled positioners...')
            posids = None
        else:
            selection = user_text.split()
            selection = list(set(selection)) # Stop repeated entries.
            kw = 'busids' if 'can' in selection[0] else 'posids'
            kwarg.update({kw: selection})
            enabled_only = False
            self.print(f'{len(selection)} items specified, '
                       'allowing disabled positioners to be selected...')
            ret = self.ptlm.get_positioners(enabled_only=enabled_only, **kwarg)
            posinfo = pd.concat(list(ret.values()))
            posids = sorted(posinfo['DEVICE_ID'])
            self.print(f'Selected {len(posids)} positioners')
        return posids

    def fvc_measure(self, exppos=None, match_radius=None, matched_only=True,
                    check_unmatched=False, test_tp=False, num_meas=1):
        '''Measure positioner locations.
        
        INPUTS:  exppos ... expected positions in some format (?), though if None
                            the code will just use the current expected positions
                            (almost always what you want)
                 match_radius ... passed to fvc proxy
                 matched_only ... passed to fvc proxy
                 check_unmatched ... passed to PetalApp.handle_fvc_feedback()
                 test_tp ... passed to PetalApp.handle_fvc_feedback()
                 num_meas ... how many FVC images to take (the results will be
                                                           median-ed in XY space)
        
        OUTPUTS: exppos ... expected positions in some format (?)
                 meapos ... measured results. pandas dataframe with columns:
                            ['DQ', 'DS', 'FLAGS', 'FWHM', 'MAG', 'MEAS_ERR', 'Q', 'S']
                            Additionally, each sub-measurement (see num_meas arg) will
                            have a complete set of identically-named columns, except with
                            a suffix like '0', '1', '2' ... apended. So for example you
                            would see columns ['Q0', 'Q1', 'Q2', ...]. If num_meas=1, you
                            would only find 'Q0', and it would be identical to the values
                            in column 'Q'.
                 matched ... set of positioner ids matched to measured centroids
                 unmatched ... set of unmatched positioner ids
        '''
        assert self.fvc, 'fvc_measure called when FVC not availible in PECS instance!'
        assert num_meas > 0, f'argued num_meas={num_meas}, but must do at least one measurement'
        if match_radius is None :
            match_radius = self.match_radius
        if exppos is None:
            # participating_petals=None gets responses from all
            # not just selected petals
            exppos = (self.ptlm.get_positions(
                          return_coord='QS', drop_devid=False,
                          participating_petals=self.illuminated_ptl_roles)
                      .sort_values(by='DEVICE_ID').reset_index(drop=True))
        if np.any(['P' in device_id for device_id in exppos['DEVICE_ID']]):
            self.print('Expected positions of positioners by PetalApp '
                       'are contaminated by fiducials.')
        centers = self.ptlm.get_centers(return_coord='QS', drop_devid=False,
                                        participating_petals=self.illuminated_ptl_roles)
        seqid = None
        if hasattr(self, 'exp'):
            seqid = self.exp.id
        ok = True
        for i in range(num_meas):
            self.print(f'Calling FVC.measure with exptime = {self.exptime} s, '
                   f'expecting {len(exppos)} backlit positioners. Image {i+1} of {num_meas}.')
            try:
                positions =self.fvc.measure(expected_positions=exppos, seqid=seqid,
                                            exptime=self.exptime, match_radius=match_radius,
                                            matched_only=matched_only,
                                            all_fiducials=self.all_fiducials,centers=centers)
                self.print('Finished FVC.measure')
                
                # Positions is either a pandas dataframe or else a dictionary (from np.rec_array).
                # Is may empty when no posiitoners are present. Sometimes I've also seen a list...                
                assert len(positions) > 0, 'Return from fvc.measure is empty! Check that positioners are back illuminated!'
            except:
                self.print('FVC measure failed! Check if fibers are visable in the image!')
                ok = False
                break
            else:
                this_meapos = pd.DataFrame(positions).rename(columns=
                                {'id': 'DEVICE_ID'}).set_index('DEVICE_ID').sort_index()
                if np.any(['P' in device_id for device_id in this_meapos.index]):
                    self.print('Measured positions of positioners by FVC '
                               'are contaminated by fiducials.')
                this_meapos.columns = this_meapos.columns.str.upper()  # clean up header to save
                this_meapos = self._batch_transform(this_meapos, 'QS', 'obsXY',
                                                    participating_petals=self.illuminated_ptl_roles)
                if i == 0:
                    meapos = this_meapos
                for column in ['obsX', 'obsY', 'Q', 'S']:
                    meapos[f'{column}{i}'] = this_meapos[column]
                self.iteration += 1
        if ok:
            median_columns = ['obsX', 'obsY']
            for column in median_columns:
                these_columns = [f'{column}{i}' for i in range(num_meas)]
                medians = meapos[these_columns].median(axis=1)
                meapos[column] = medians
            meapos = self._batch_transform(meapos, 'obsXY', 'QS',
                                           participating_petals=self.illuminated_ptl_roles)
            exppos = (exppos.rename(columns={'id': 'DEVICE_ID'})
                      .set_index('DEVICE_ID').sort_index())
            exppos.columns = exppos.columns.str.upper()
            # recover flags
            meapos['FLAGS'] |= exppos.loc[list(meapos.index)]['FLAGS']
            # find the posids that are unmatched, missing from FVC return
            matched = set(exppos.index) & set(meapos.index)
            unmatched = set(exppos.index) - matched
            if len(unmatched) == 0:
                self.print(f'All {len(exppos.index)} back-illuminated '
                           f'positioners measured by FVC')
            else:
                self.print(f'Missing {len(unmatched)} of expected backlit fibers:'
                           f'\n{sorted(unmatched)}')
            fvc_data = pd.concat([meapos,
                        exppos[exppos.index.isin(meapos.index)][['PETAL_LOC','DEVICE_LOC']]],
                        axis=1)
            submeas = {}
            if num_meas > 1:
                submeas = self.summarize_submeasurements(meapos)

            # gather tracked angles *prior* to their possible updates by handle_fvc_feedback()
            posids = exppos.index.tolist()
            for key in ['posintT', 'posintP']:
                this_data = self.quick_query_df(key=key, posids=posids,
                                                participating_petals=self.illuminated_ptl_roles)
                exppos = exppos.join(this_data, on='DEVICE_ID')

            control = {'timeout': self.fvc_feedback_timeout}
            self.ptlm.handle_fvc_feedback(fvc_data, check_unmatched=check_unmatched,
                                          test_tp=test_tp, auto_update=True,
                                          err_thresh=self.max_err, up_tol=self.tp_tol,
                                          up_frac=self.tp_frac, postscript=submeas,
                                          control=control, participating_petals=self.illuminated_ptl_roles)
        else:
            exppos, meapos, matched, unmatched = pd.DataFrame(), pd.DataFrame(), set(), set()

        return exppos, meapos, matched, unmatched

    def move_measure(self, request=None, match_radius=None, check_unmatched=False,
                     test_tp=False, anticollision='default', num_meas=1):
        '''
        Wrapper for often repeated moving and measuring sequence.
        Returns data merged with request.
        
        INPUTS:  request        ... pandas dataframe that includes columns:
                                    ['DEVICE_ID', 'COMMAND', 'X1', 'X2', 'LOG_NOTE']
                                ... if request == None, the "prepare_move" call is suppressed.
                                    This is for special cases where move requests have already
                                    been independently registered and scheduled by direct calls
                                    to petal(s).
                                    
                 match_radius   ... passed to fvc_measure()
                 
                 check_umatched ... passed to PetalApp.handle_fvc_feedback()
                 
                 test_tp        ... passed to PetalApp.handle_fvc_feedback()
                 
                 anticollison   ... mode like 'default', 'freeze', 'adjust', 'adjust_requested_only' or None
                 
                 num_meas       ... how many FVC images to take (the results will be median-ed)
                         
        OUTPUTS: result ... pandas dataframe that includes columns:
                            ['DQ', 'DS', 'FLAGS', 'FWHM', 'MAG', 'MEAS_ERR',
                             'mea_Q', 'mea_S', 'COMMAND', 'MOVE_VAL1',
                             'MOVE_VAL2', 'LOG_NOTE', 'BUS_ID', 'DEVICE_LOC',
                             'PETAL_LOC', 'STATUS', 'posintT', 'posintP']
                            
                            And also, ['Q0', 'Q1', 'Q2', etc ...]
                            for the sub-measurements, depending on the value of
                            argument num_meas. See fvc_measure() for details on
                            this. These sub-meas columns I am not bothering to
                            rename with the 'mea_' prefix.
                            
                            The dataframe has index column 'DEVICE_ID'
        '''
        self.print(f'Moving positioners... Exposure {self.exp.id}, iteration {self.iteration}')
        self.ptlm.set_exposure_info(self.exp.id, self.iteration)
        should_prepare_move = True
        if request is None:
            should_prepare_move = False
            self.print('Skipping "prepare_moves" call, under assumption of independently requested / scheduled moves.')
            requested_posids = self.ptlm.get_requested_posids(kind='all')
            all_requested = set()
            for posids in requested_posids.values():
                all_requested |= posids
            dummy_req = {'DEVICE_ID': list(all_requested), 'COMMAND': 'dummy_cmd', 'X1': 0.0, 'X2': 0.0, 'LOG_NOTE': ''}
            if len(all_requested) == 0:
                self.print('WARNING: No independently requested / scheduled moves. Check for rejected requests, IE invalid or all devices disabled.')
                dummy_req = {'DEVICE_ID': list(self.posids), 'COMMAND': 'dummy_cmd', 'X1': np.nan, 'X2': np.nan, 'LOG_NOTE': ''}
            request = pd.DataFrame(dummy_req)
        if 'PETAL_LOC' not in request.columns:
            request = request.merge(self.petal_locs, on='DEVICE_ID')
        if should_prepare_move:
            self.ptlm.prepare_move(request, anticollision=anticollision)
        self.ptlm.execute_move(reset_flags=False, control={'timeout': self.execute_move_timeout})
        exppos, meapos, matched, _ = self.fvc_measure(
            exppos=None, matched_only=True, match_radius=match_radius, 
            check_unmatched=check_unmatched, test_tp=test_tp, num_meas=num_meas)
        result = self._merge_match_and_rename_fvc_data(request, meapos, matched, exppos)
        self.ptlm.clear_exposure_info()
        return result
    
    def rehome_and_measure(self, posids, axis='both', debounce=True, log_note='',
                           match_radius=None, check_unmatched=False, test_tp=False,
                           anticollision='freeze'):
        '''Wrapper for sending rehome command and then measuring result.
        Returns whatever fvc_measure returns.
        '''
        assert axis in {'both', 'phi', 'phi_only', 'theta', 'theta_only'}
        assert debounce in {True, False}
        if anticollision not in {'freeze', None}:
            anticollision = 'freeze'
            self.print(f'Anticollision method {anticollision} is not supported during rehome.' +
                       f' Reverting to {anticollision}')
        self.print(f'Rehoming positioners, axis={axis}, anticollision={anticollision}' +
                   f', debounce={debounce}, exposure={self.exp.id}, iteration={self.iteration}')
        return self._rehome_or_park_and_measure(ids=posids, axis=axis, debounce=debounce,
                                                log_note=log_note, match_radius=match_radius, 
                                                check_unmatched=check_unmatched,
                                                test_tp=test_tp, anticollision=anticollision)
    
    def park_and_measure(self, posids, mode='normal', coords='poslocTP', log_note='',
                         match_radius=None, check_unmatched=False, test_tp=False, theta=0):
        '''Wrapper for sending park_positioners command and then measuring result.
        Returns whatever fvc_measure returns.
        '''
        assert mode in {'normal', 'center'}
        assert coords in {'posintTP', 'poslocTP', 'intTlocP', 'locTintP'}
        self.print(f'Parking positioners, mode={mode}, coords={coords}' +
                   f', exposure={self.exp.id}, iteration={self.iteration}')
        return self._rehome_or_park_and_measure(move='park', ids=posids, mode=mode,
                                                coords=coords, log_note=log_note,
                                                match_radius=match_radius, 
                                                check_unmatched=check_unmatched,
                                                test_tp=test_tp, theta=0,
                                                expid=self.exp.id, iteration=self.iteration)
        
    def _rehome_or_park_and_measure(self, move='rehome', **kwargs):
        '''Common operations for both "rehome_and_measure" and "park_and_measure".

        KF - note RE expid/iteration, park_positioners in petalman handles the setting
        of expid/iteration (for use by OCS), rehome needs this done here.
        '''
        funcs = {'rehome': self.ptlm.rehome_pos,
                 'park': self.ptlm.park_positioners}
        move_args = {'rehome': {'ids', 'axis', 'anticollision', 'debounce', 'log_note'},
                     'park': {'ids', 'mode', 'coords', 'log_note', 'theta', 'expid', 'iteration'}}
        meas_args = {'match_radius', 'check_unmatched', 'test_tp'}
        missing_args = (move_args[move] | meas_args) - set(kwargs)
        assert move in funcs, f'unrecognized move type {move}'
        assert not(any(missing_args)), f'missing args {missing_args}'
        posids = kwargs['ids']
        if move == 'rehome':
            self.ptlm.set_exposure_info(self.exp.id, self.iteration)
        enabled = self.get_enabled_posids(posids)
        move_kwargs = {key: kwargs[key] for key in move_args[move] if key != 'ids'}
        move_kwargs['ids'] = enabled
        if 'control' not in list(move_kwargs.keys()):
            move_kwargs['control'] = {'timeout': self.execute_move_timeout}
        else:
            move_kwargs['control']['timeout'] = self.execute_move_timeout
        funcs[move](**move_kwargs)
        
        # 2020-07-21 [JHS] dissimilar results than move_measure func, since no "request" data structure here
        meas_kwargs = {key:kwargs[key] for key in meas_args}
        meas_kwargs.update({'exppos': None, 'matched_only': True})
        # 20200-10-21 [KF] only use the meapos as result, not the whole tuple return
        _, result, _, _ = self.fvc_measure(**meas_kwargs)
        self.ptlm.clear_exposure_info()
        return result

    def home_adc(self):
        if self.adc_available:
            do_home = self._parse_yn(self.input('Home ADC? (y/n): ')) if self.interactive else True
            if do_home:
                try:
                    adc = SimpleProxy('ADC')
                    self.print('Homing ADC...')
                    retcode = adc._send_command('home', controllers=[1, 2])
                    self.print(f'ADC.home returned code: {retcode}')
                except Exception as e:
                    print(f'Exception homing ADC, {e}')

    def turn_on_fids(self):
        do_on = self._parse_yn(self.input('Turn on fiducials (y/n): ')) if self.interactive else True
        if do_on:
            responses = self.ptlm.set_fiducials(setting='on',
                            participating_petals=self.fid_petals)
            # Set fiducials for all availible petals
        else:
            responses = self.ptlm.get_fiducials()
        self.print(f'Petals report fiducials in the following states: {responses}')

    def turn_on_illuminator(self):
        if self.specman_available:
            do_on = self._parse_yn(self.input('Turn on illuminator? (y/n): ')) if self.interactive else True
            if do_on:
                try:
                    specman = SimpleProxy('SPECMAN')
                    self.print('Turning on illuminator')
                    specs = [f'SP{p}' for p in self.illuminated_pcids]
                    retcode = specman._send_command('illuminate', action='on', participating_spectrographs=specs)
                    self.print(f'SPECMAN.illuminate returned code: {retcode}')
                except Exception as e:
                    print(f'Exception when trying to turn on illumination, {e}') 

    def turn_off_fids(self):
        do_off = self._parse_yn(self.input('Turn off fiducials (y/n): ')) if self.interactive else True
        if do_off:
            responses = self.ptlm.set_fiducials(setting='off',
                            participating_petals=self.fid_petals)
            # Set fiducials for all availible petals
        else:
            responses = self.ptlm.get_fiducials()
        self.print(f'Petals report fiducials in the following states: {responses}')

    def turn_off_illuminator(self):
        if self.specman_available:
            do_off = self._parse_yn(self.input('Turn off illuminator? (y/n): ')) if self.interactive else True
            if do_off:
                try:
                    specman = SimpleProxy('SPECMAN')
                    self.print('Turning off illuminator')
                    specs = [f'SP{p}' for p in self.illuminated_pcids]
                    retcode = specman._send_command('illuminate', action='off', participating_spectrographs=specs)
                    self.print(f'SPECMAN.illuminate returned code: {retcode}')
                except Exception as e:
                    print(f'Exception when trying to turn off illumination, {e}') 


    def fvc_collect(self):
        destination = os.path.join(
            self.exposure_dir, pc.dir_date_str(t=self.start_time),
            f'{self.exp.id:08}')
        # os.makedirs(destination, exist_ok=True)  # no permission anyway
        try:
            fvccoll = SimpleProxy('FVCCOLLECTOR')
            retcode = fvccoll._send_command('configure')
            self.print(f'FVCCollector.configure returned code: {retcode}')
            retcode = fvccoll._send_command(
                'collect', expid=self.exp.id, output_dir=destination,
                logbook=False, remove_after_transfer=self.fvc_cleanup)
            self.print(f'FVCCollector.collect returned code: {retcode}')
            if retcode == 'SUCCESS':
                wait_time = 60 #wait for data transfer, maybe make this calculated time*iterations?
                exit_cond = lambda : os.path.isfile(f'{destination}/fvc-{self.exp.id:08}.fits.fz')
                self.print(f'Sleeping for maximum {wait_time} to wait for FVC data transfer...')
                self.countdown_sec(wait_time, break_condition=exit_cond) #wait a minute for data transfer/break early if file exists
                self.print('FVC data associated with exposure ID '
                           f'{self.exp.id} collected to: {destination}')
            expserv = SimpleProxy('EXPOSURESERVER')
            retcode = expserv._send_command('distribute', source=destination, expid=self.exp.id)
            self.print(f'ExposureServer.distribute returned code: {retcode}')
            if retcode == 'SUCCESS':
                self.print(f'symlink created for: {destination}')
        except Exception as e:
            self.print(f'FVC collection failed with exception: {e}')

    @staticmethod
    def countdown_sec(t, break_condition=lambda : False):  # in seconds
        t = int(t)
        t_wait = 0
        for i in reversed(range(t)):
            sys.stderr.write(f'\rSleeping... ({i} s / {t} s)')
            time.sleep(1)
            t_wait += 1
            sys.stdout.flush()
            if break_condition():
                break
        if t_wait < t:
            print(f'\nSleep finished early ({t_wait} s)')
        else:
            print(f'\nSleep finished ({t} s)')

    def pause(self, press_enter=False):
        '''pause operation between positioner moves for heat monitoring'''
        if self.pause_interval is None or press_enter:
            self.input('Paused for heat load monitoring for unspecified interval. '
                       'Press enter to continue: ')
        elif self.pause_interval == 0:  # no puase needed, just continue
            self.print('pause_interval = 0, continuing without pause...')
        elif self.pause_interval > 0:
            self.print(f'Pausing for {self.pause_interval} s...')
            self.countdown_sec(self.pause_interval)

    def decorate_note(self, log_note=''):
        '''Adds additional information to a positioner log note. The intended
        usage is to include important facts known only to PECS, when constructing
        the LOG_NOTE field that gets saved to the posmovedb.'''
        log_note = pc.join_notes(log_note, f'expid={self.exp.id}')
        log_note = pc.join_notes(log_note, f'use_desimeter={self.use_desimeter}')
        return log_note
                
    def get_enabled_posids(self, posids='sub', include_posinfo=False):
        '''Returns list of the enabled posids.
        
        INPUTS:
            posids ... 'sub' --> return enabled subset of PECS' self.posids property
                   ... 'all' --> return all currently enabled posids
                   ... some iterable --> return just the enabled subset of iterable
        
            include_posinfo ... boolean, whether to return a second item with
                                additional "posinfo" data
        '''
        if posids == 'sub':
            selected_posids = self.posids
        elif posids == 'all':
            selected_posids = None
        else:
            try:
                iterator = iter(posids)
                selected_posids = [p for p in iterator]
            except:
                assert False, f'error, could not iterate arg posids={posids}'
        ret = self.ptlm.get_positioners(enabled_only=True, posids=selected_posids)
        posinfo = pd.concat(list(ret.values()))
        posinfo = posinfo.set_index('DEVICE_ID')
        posids = sorted(posinfo.index)
        if include_posinfo:
            return posids, posinfo
        return posids
    
    def summarize_submeasurements(self, meapos):
        '''Collects submeasurements (like from fvc_measure with num_meas > 1)
        by positioner ID into suitable strings for logging purposes.
        
        INPUTS:   meapos ... dataframe including columns ['Q0', 'S0', 'Q1', 'S1', ...]
                             and index 'DEVICE_ID'

        OUTPUTS:  Dictionary with keys = posids and values = strings.
                  Only includes entries for posids as defined by calling
                  get_enabled_posids(posids, False)
        
        In cases where only the 0-th terms are found (i.e. we had num_meas==1)
        the strings in the return dict will be empty, ''.
        '''
        posids = set(meapos.index)
        num_meas = 0
        crazy_number_of_iterations = 100
        for i in range(crazy_number_of_iterations):
            if f'Q{i}' in meapos.columns and f'S{i}' in meapos.columns:
                num_meas += 1
            else:
                break
        assert num_meas != 0, 'meapos is missing columns "Q0" and/or "S0"'
        if num_meas == 1:
            return {posid: '' for posid in posids}
        out = {}
        for posid in posids:
            data = meapos.loc[posid]
            this_dict = {'Q':[], 'S':[], 'obsX':[], 'obsY':[]}
            for key in this_dict:
                for i in range(num_meas):
                    this_dict[key] += [data[f'{key}{i}']]
            fmt1 = lambda X: str([f'{x:.4f}' for x in X]).replace("'", '')
            s = 'submeas={'
            for k,v in this_dict.items():
                s += f'\'{k}\': {fmt1(v)}, '
            s = s[:-2] + '}'
            out[posid] = s
        return out
    
    def quick_query(self, key=None, op='', value='', posids='all', mode='iterable', participating_petals=None):
        '''Returns a collection containing values for all petals of a quick_query()
        call. See documentation for petal.quick_query() for details.
        '''
        # This type of call to PetalMan (in the next line) returns *either* dict with keys = like
        # 'PETAL1', 'PETAL2', etc, and corresponding return data from those petals, OR just the return
        # data from that one petal. It's not perfectly clear when this does or doesn't happen, but
        # from discussion with Kevin and some trials at the CONSOLE this ought to be ok here.
        if participating_petals is None:
            participating_petals = self.ptlm.participating_petals
        data_by_petal = self.ptlm.quick_query(key=key, op=op, value=value, posids=posids,
                                              mode=mode, skip_unknowns=True, participating_petals=participating_petals)
        if isinstance(data_by_petal, (str, list, tuple)):
            return data_by_petal
        check_val = data_by_petal[list(data_by_petal.keys())[0]]
        if isinstance(check_val, dict):
            combined = {}
            for data in data_by_petal.values():
                combined.update(data)
        elif isinstance(check_val, set):
            combined = set()
            for data in data_by_petal.values():
                combined |= data
        elif isinstance(check_val, list):
            combined = []
            for data in data_by_petal.values():
                combined += data
        elif isinstance(check_val, str):
            combined = '\n'.join(data_by_petal.values())
        else:
            assert False, f'unrecognized return type {type(check_val)} from quick_query()'
        return combined
        
    def quick_query_df(self, key, posids='all', participating_petals=None):
        '''Wrapper for quick_query which returns a pandas DataFrame, whose
        index is 'DEVICE_ID' and data column is key.
        '''
        if participating_petals is None:
            participating_petals = self.ptlm.participating_petals
        data = self.quick_query(key=key, posids=posids, participating_petals=participating_petals)
        ordered_posids = list(data)
        ordered_data = [data[posid] for posid in posids]
        listed = {'DEVICE_ID': ordered_posids, key: ordered_data}
        df = pd.DataFrame(listed)
        df = df.set_index('DEVICE_ID')
        return df
    
    def _merge_match_and_rename_fvc_data(self, request, meapos, matched, exppos):
        '''Returns results of fvc measurement after checking for target matches
        and doing some pandas juggling, very specifc to the other data interchange
        formats in pecs etc.
        
        request ... pandas dataframe with index column DEVICE_ID and request data
        meapos ... pandas dataframe with index column DEVICE_ID and fvc data
        matched ... set of posids
        exppos ... pandas dataframe with index column DEVICE_ID and expected positions data
        '''
        # meapos may contain not only matched but all posids in expected pos
        posids_to_match = set(request['DEVICE_ID']) | set(self.posids)
        matched_df = meapos.loc[sorted(matched & posids_to_match)]
        merged = matched_df.merge(request, on='DEVICE_ID')
        if not(merged.index.name == 'DEVICE_ID'):
            merged = merged.set_index('DEVICE_ID')
        
        # columns get renamed
        merged.rename(columns={'X1': 'MOVE_VAL1', 'X2': 'MOVE_VAL2',
                               'Q': 'mea_Q', 'S': 'mea_S', 'FLAG': 'FLAGS'},
                      inplace=True)
        mask = merged['FLAGS'].notnull()
        merged.loc[mask, 'STATUS'] = pc.decipher_posflags(merged.loc[mask, 'FLAGS'])
        
        # get expected (tracked) posintTP angles --- these are now *after* any
        # updating by test_and_update_tp, so may not match the original tracked
        # values --- hence the special suffix
        suffix = '_after_updates'
        posids = exppos.index.tolist()
        for key in ['posintT', 'posintP']:
            this_data = self.quick_query_df(key=key, posids=posids, participating_petals=self.illuminated_ptl_roles)
            exppos = exppos.join(this_data, on='DEVICE_ID', rsuffix=suffix)
        posint_keys = [key for key in exppos.columns if 'posint' in key]
        result = merged.merge(exppos[posint_keys], on='DEVICE_ID')
        
        return result
    
    def _batch_transform(self, frame, cs1, cs2, participating_petals=None):
        '''Calculate values in cs2 for a frame that has columns for coordinate
        system cs1. E.g. cs1='QS', cs2='obsXY'. Index is 'DEVICE_ID'. Then return
        a new frame that now includes new columns (or replaces existing columns)
        for coordinate system 2. E.g. you could start with only 'Q' and 'S' in
        frame, and then 'obsX' and 'obsY' will be added.

        DOES NOT change frame
        '''
        if participating_petals is None:
            participating_petals = self.ptlm.participating_petals
        if not(frame.index.name == 'DEVICE_ID'):
            df = frame.set_index('DEVICE_ID')
            reset_index = True
        else:
            df = frame.copy()
            reset_index = False
        ret = self.ptlm.batch_transform(df, cs1, cs2, participating_petals=participating_petals)
        df = pd.concat([x for x in ret.values()])
        if reset_index:
            df.reset_index(inplace=True)
        return df
