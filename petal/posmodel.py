import posstate
import postransforms
import posconstants as pc

class PosModel(object):
    """Software model of the physical positioner hardware.
    Takes in local (x,y) or (th,phi) targets, move speeds, converts
    to degrees of motor shaft rotation, speed, and type of move (such
    as cruise / creep / backlash / hardstop approach).

    One instance of PosModel corresponds to one PosState to physical positioner.
    """
    def __init__(self, state=None, petal_alignment=None, printfunc=print):
        self.DEBUG = 1
        if not(state):
            self.state = posstate.PosState()
        else:
            self.state = state
        self.printfunc = printfunc
        self.state.set_posmodel_cache_refresher(self.refresh_cache)
        self.trans = postransforms.PosTransforms(this_posmodel=self, petal_alignment=petal_alignment)
        self.axis = [None, None]
        self.axis[pc.T] = Axis(self, pc.T, printfunc=self.printfunc)
        self.axis[pc.P] = Axis(self, pc.P, printfunc=self.printfunc)
        posid = self.posid
        self._timer_update_rate          = 18e3   # Hz
        self._stepsize_creep             = 0.1    # deg
        self._motor_speed_cruise         = {pc.T: 9900.0 * 360.0 / 60.0, pc.P: 9900.0 * 360.0 / 60.0} # deg/sec (= RPM *360/60)
        self._stepsize_cruise            = {pc.T: 3.3, pc.P: 3.3} # deg/step
        if self.state._val['ZENO_MOTOR_P'] is True:
            if self.DEBUG:
                self.printfunc(f'PosModel: new linphi posid = {posid}')  # DEBUG
#           self.linphi_params['LAST_P_DIR'] = 1    # 1 is CCW, -1 is CW
            self._stepsize_cruise[pc.P] = 0.1 * float(pc.P_zeno_speed)
            self._motor_speed_cruise[pc.P] = (18000 * 60/3600 * pc.P_zeno_speed) * 360.0/60.0 # RPM * 360/60 = deg/sed
        self._spinupdown_dist_per_period = {pc.T: sum(range(round(self._stepsize_cruise[pc.T]/self._stepsize_creep) + 1))*self._stepsize_creep,
                                            pc.P: sum(range(round(self._stepsize_cruise[pc.P]/self._stepsize_creep) + 1))*self._stepsize_creep}
        self.refresh_cache()

    @property
    def is_linphi(self):
        return self.state._val['ZENO_MOTOR_P'] is True
        
    def get_zeno_scale(self, which):    # specify 'SZ_CW_P', 'SZ_CCW_P', or the _T varieties
        scale = self.state._val[which]
        if scale is None:
            scale = 1.0
        return scale

    def _load_cached_params(self):
        '''Do this *after* refreshing the caches in the axis instances.'''
        self._abs_shaft_speed_cruise_T = abs(self._motor_speed_cruise[pc.T] / self.axis[pc.T].signed_gear_ratio)
        self._abs_shaft_speed_cruise_P = abs(self._motor_speed_cruise[pc.P] / self.axis[pc.P].signed_gear_ratio)
        self._abs_shaft_spinupdown_distance_T = abs(self.axis[pc.T].motor_to_shaft(self._spinupdown_distance(pc.T)))
        if self.is_linphi:
            if self.DEBUG:
                self.printfunc(f'_load_cached_params: LinPhi posid = {self.posid}')  # DEBUG
            speed = pc.P_zeno_speed
            ramp = pc.P_zeno_ramp
            gear_ratio = pc.gear_ratio[self.state._val['GEAR_TYPE_P']]
            self._abs_shaft_spinupdown_distance_P = speed*(speed+1)*ramp/20/gear_ratio # From DESI-1710 Motor Speed Parameters Spreadsheet
            if self.DEBUG:
                self.printfunc(f'Phi Spinupdown = {self._abs_shaft_spinupdown_distance_P}')  # DEBUG
        else:
            self._abs_shaft_spinupdown_distance_P = abs(self.axis[pc.P].motor_to_shaft(self._spinupdown_distance(pc.P)))

    def refresh_cache(self):
        """Reloads state parameters with cached values."""
        for axis in self.axis:
            axis._load_cached_params()
        self._load_cached_params()

    @property
    def _motor_speed_creep(self):
        """Returns motor creep speed (which depends on the creep period) in deg/sec."""
        return self._timer_update_rate * self._stepsize_creep / self.state._val['CREEP_PERIOD']  # deg/sec

    def _spinupdown_distance(self, axisid):
        """Returns distance at the motor shaft in deg over which to spin up to cruise speed or down from cruise speed."""
        if self.is_linphi and axisid == pc.P:
            speed = pc.P_zeno_speed
            ramp = pc.P_zeno_ramp
#           gear_ratio = pc.gear_ratio[self.state._val['GEAR_TYPE_T']]
            sud = speed*(speed+1)*ramp/20 * pc.P_zeno_ramp # From DESI-1710 Motor Speed Parameters Spreadsheet
        elif self.state._val['CURR_SPIN_UP_DOWN'] == 0:
            sud = 0  # special case, where user is trying to prevent FIPOS from doing the physical spin-up down
        else:
            sud = self._spinupdown_dist_per_period[axisid] * self.state._val['SPINUPDOWN_PERIOD']
        return sud

    @property
    def posid(self):
        """Returns the positioner id string."""
        return self.state._val['POS_ID']

    @property
    def canid(self):
        """Returns the positioner hardware id on the CAN bus."""
        return self.state._val['CAN_ID']

    @property
    def busid(self):
        """Returns the name of the can bus that the positioner belongs to."""
        return self.state._val['BUS_ID']

    @property
    def deviceloc(self):
        """Returns the device location id (on the petal) of the positioner."""
        return self.state._val['DEVICE_LOC']

    @property
    def is_enabled(self):
        """Returns whether the positioner has its control enabled or not."""
        return self.state._val['CTRL_ENABLED']

    @property
    def axis_locks(self):
        """Returns 1x2 tuple of booleans, stating whether the (theta, phi) axes are locked."""
        return (self.axis[pc.T].is_locked, self.axis[pc.P].is_locked)

    @property
    def classified_as_retracted(self):
        """Returns whether the positioner has been classified as retracted or not."""
        return self.state._val['CLASSIFIED_AS_RETRACTED']

    @property
    def expected_current_posintTP(self):
        """Returns the internally-tracked expected position of the
        theta and phi shafts at the output of the gearbox."""
        return self.axis[pc.T].pos, self.axis[pc.P].pos

    @property
    def expected_current_poslocTP(self):
        """Returns the expected position of theta and phi bodies, as seen by
        an external observer."""
        posintTP = self.expected_current_posintTP
        return self.trans.posintTP_to_poslocTP(posintTP)

    @property
    def expected_current_position(self):
        """
        Returns a general dictionary of the current expected position in all
        the various coordinate systems.
        The keys are:
            'posintTP'  tuple in deg, independent variable,
                        the internally-tracked expected position of the
                        theta shaft at the output of the gearbox
            'poslocTP'  tuple in deg, pos local coordinates with offsets
            'poslocXY   tuple in, mm, expected local x position
            'QS'        tuple in (deg, mm), expected global Q position
            'flatXY'    tuple in, mm, dependent variable,
                        expected local x in a system where focal surface
                        curvature is flattened out to an approximate plane
            'ptlXY'     tuple in mm, petal local XY projection of ptlXYZ
            'obsXY'     tuple in mm, dependent, expected global x position
        """
        posintTP = self.expected_current_posintTP
        QS = self.trans.posintTP_to_QS(posintTP)
        return {'posintTP': posintTP,
                'poslocTP': self.trans.posintTP_to_poslocTP(posintTP),
                'poslocXY': self.trans.posintTP_to_poslocXY(posintTP),
                'flatXY': self.trans.posintTP_to_flatXY(posintTP),
                'ptlXY': self.trans.posintTP_to_ptlXY(posintTP),
                'obsXY': tuple(
                        self.trans.QS_to_obsXYZ(QS, cast=True).flatten()[:2]),
                'QS': QS}

    @property
    def expected_current_position_str(self):
        """One-line string summarizing current expected position.
        """
        pos = self.expected_current_position
        s = ('Q:{:8.3f}{}, S:{:8.3f}{} | '
             'flatX:{:8.3f}{}, flatY:{:8.3f}{} | '
             'obsX:{:8.3f}{}, obsY:{:8.3f}{} | '
             'ptlX:{:8.3f}{}, ptlY:{:8.3f}{} | '
             'poslocX:{:8.3f}{}, poslocY:{:8.3f}{} | '
             'poslocT:{:8.3f}{}, poslocP:{:8.3f}{} | '
             'posintT:{:8.3f}{}, posintP:{:8.3f}{}'). \
            format(pos['QS'][0],       pc.deg, pos['QS'][1],       pc.mm,
                   pos['flatXY'][0],   pc.mm,  pos['flatXY'][1],   pc.mm,
                   pos['obsXY'][0],    pc.mm,  pos['obsXY'][1],    pc.mm,
                   pos['ptlXY'][0],    pc.mm,  pos['ptlXY'][1],    pc.mm,
                   pos['poslocXY'][0], pc.mm,  pos['poslocXY'][1], pc.mm,
                   pos['poslocTP'][0], pc.deg, pos['poslocTP'][1], pc.deg,
                   pos['posintTP'][0], pc.deg, pos['posintTP'][1], pc.deg)
        return s

    @property
    def targetable_range_posintT(self):
        """Returns a [1x2] array of theta_min, theta_max, after subtracting buffer zones near the hardstops.
        The return is in the posintTP coordinates, not poslocTP. Understand therefore that OFFSET_T is not included in these values.
        """
        return self.axis[pc.T].debounced_range

    @property
    def targetable_range_posintP(self):
        """Returns a [1x2] array of phi_min, phi_max, after subtracting buffer zones near the hardstops.
        The return is in the posintTP coordinates, not poslocTP. Understand therefore that OFFSET_P is not included in these values.
        """
        return self.axis[pc.P].debounced_range

    @property
    def full_range_posintT(self):
        """Returns a [1x2] array of [theta_min, theta_max], from hardstop-to-hardstop.
        The return is in the posintTP coordinates, not poslocTP. Understand therefore that OFFSET_T is not included in these values.
        """
        return self.axis[pc.T].full_range

    @property
    def full_range_posintP(self):
        """Returns a [1x2] array of [phi_min, phi_max], from hardstop-to-hardstop.
        The return is in the posintTP coordinates, not poslocTP. Understand therefore that OFFSET_P is not included in these values.
        """
        return self.axis[pc.P].full_range

    @property
    def abs_shaft_speed_cruise_T(self):
        """Returns the absolute output shaft speed (deg/sec), in cruise mode, of the theta axis.
        """
        return self._abs_shaft_speed_cruise_T

    @property
    def abs_shaft_speed_cruise_P(self):
        """Returns the absolute output shaft speed (deg/sec), in cruise mode, of the phi axis.
        """
        return self._abs_shaft_speed_cruise_P

    @property
    def abs_shaft_spinupdown_distance_T(self):
        '''Acceleration / deceleration distance on theta axis, when spinning up to / down from cruise speed.'''
        return self._abs_shaft_spinupdown_distance_T

    @property
    def abs_shaft_spinupdown_distance_P(self):
        '''Acceleration / deceleration distance on phi axis, when spinning up to / down from cruise speed.'''
        return self._abs_shaft_spinupdown_distance_P

    @property
    def theta_hardstop_ambiguous_zone(self):
        '''Returns a 1x2 tuple of (AMBIG_MIN, AMBIG_MAX), describing the range in which a measured
        a measured theta value could mean it's in fact on *either* side of the hardstop (in which case
        the situation must be resolved by some physical motion). The output values are always returned
        as positive numbers within the range [0, 360]. So like (+170, +190).'''
        full_range = self.full_range_posintT
        ambig_min = (full_range[0] - pc.theta_hardstop_ambig_tol) % 360
        ambig_max = (full_range[1] + pc.theta_hardstop_ambig_tol) % 360
        for key in {'ambig_min', 'ambig_max'}:
            assert 0 <= eval(key) <= 360, f'{self.posid} ambig_min={eval(key)} is not within [0,360]. full_range={full_range}'
        return (ambig_min, ambig_max)

    @property
    def in_theta_hardstop_ambiguous_zone(self):
        '''Returns boolean whether positioner is currently in the ambiguous either-side-of-the-hardstop
        zone.'''
        ambig_range = self.theta_hardstop_ambiguous_zone
        t_test = self.axis[pc.T].pos
        t_test %= 360
        return ambig_range[0] <= t_test <= ambig_range[1]

    def true_move(self, axisid, distance, allow_cruise, limits='debounced', init_posintTP=None):
        """Input move distance on either the theta or phi axis, as seen by the
        observer, in degrees.

        The return values are formatted as a dictionary, with keys:
            'motor_step'   ... integer number of motor steps, signed according to direction
            'move_time'    ... total time the move takes [sec]
            'distance'     ... quantized distance of travel [deg], signed according to direction
            'speed'        ... approximate speed of travel [deg/sec]. unsigned. note 'move_time' is preferred for precision timing calculations
            'speed_mode'   ... 'cruise' or 'creep'

        The argument 'init_posintTP' is used to account for expected
        future shaft position changes. This is necessary for correct checking of
        software travel limits, when a sequence of multiple moves is being planned out.

        The argument 'limits' may take values:
            'debounced'
            'near_full'
            None
        See the debounced_range() and near_full_range() methods of Axis class for their
        specific meanings.
        """
        start = self.expected_current_posintTP if not init_posintTP else init_posintTP
        if self.axis[axisid].is_locked:
            new_distance = 0.0
            if self.is_linphi and axisid == pc.P and distance != new_distance:
                if self.DEBUG > 1:
                    self.printfunc(f'{self.posid} linphi Distance = {distance} changed to {new_distance}')  # DEBUG
            distance = new_distance
        elif limits:
            use_near_full_range = (limits == 'near_full')
            new_distance = self.axis[axisid].truncate_to_limits(distance, start[axisid], use_near_full_range)
            if self.is_linphi and axisid == pc.P and distance != new_distance:
                if self.DEBUG > 1:
                    self.printfunc(f'{self.posid} linphi Distance = {distance} changed to {new_distance}')  # DEBUG
            distance = new_distance
        motor_dist = self.axis[axisid].shaft_to_motor(distance)
        move_data = self.motor_true_move(axisid, motor_dist, allow_cruise)
        move_data['distance'] = self.axis[axisid].motor_to_shaft(move_data['distance'])
        move_data['speed']    = self.axis[axisid].motor_to_shaft(move_data['speed'])
        return move_data

    def motor_true_move(self, axisid, distance, allow_cruise):
        """Calculation of cruise, creep, spinup, and spindown details for a move of
        an argued distance on the axis identified by axisid.
        """
        move_data = {}
        allow_creep = False if self.is_linphi and axisid == pc.P else True
        dist_spinup = 2 * pc.sign(distance) * self._spinupdown_distance(axisid)  # distance over which accel / decel to and from cruise speed
        if allow_creep and ( not(allow_cruise) or abs(distance) <= (abs(dist_spinup) + self.state._val['MIN_DIST_AT_CRUISE_SPEED'])):
            if self.is_linphi and axisid == pc.P and abs(distance) > 0.00001: # in mm, == 10 microns
                ddist = self.axis[axisid].motor_to_shaft(distance)
                self.printfunc(f'{self.posid} linphi Distance = {ddist}, MotDist = {distance}, WARNING: creep on linphi')  # DEBUG
            move_data['motor_step']   = int(round(distance / self._stepsize_creep))
            move_data['distance']     = move_data['motor_step'] * self._stepsize_creep
            move_data['speed_mode']   = 'creep'
            move_data['speed']        = self._motor_speed_creep
            move_data['move_time']    = abs(move_data['distance']) / move_data['speed']
        else:
            dist_cruise = distance - dist_spinup
            move_data['motor_step']   = int(round(dist_cruise / self._stepsize_cruise[axisid]))
            move_data['distance']     = move_data['motor_step'] * self._stepsize_cruise[axisid] + dist_spinup
            move_data['speed_mode']   = 'cruise'
            move_data['speed']        = self._motor_speed_cruise[axisid]
            if move_data['motor_step'] == 0:
                move_data['move_time'] = 0
            else:
                move_data['move_time'] = (abs(move_data['motor_step'])*self._stepsize_cruise[axisid] + 4*self._spinupdown_distance(axisid)) / move_data['speed']
            if self.DEBUG > 1:
                if self.is_linphi and axisid == pc.P and distance != 0.0:
                    ddist = self.axis[axisid].motor_to_shaft(distance)
                    self.printfunc(f'{self.posid} linphi Distance = {ddist}, MotDist = {distance}, Spinupdown = {dist_spinup}, dist_cruise = {dist_cruise}, steps = {move_data["motor_step"]}')  # DEBUG
        return move_data

    def postmove_cleanup(self, cleanup_table):
        """Always perform this after positioner physical moves have been
        completed, to update the internal tracking of shaft positions and
        variables.
        """
        if self.state._val['CTRL_ENABLED'] is False:
            return
        net_distance = {pc.T: cleanup_table['net_dT'][-1],
                        pc.P: cleanup_table['net_dP'][-1]}
        for axis in self.axis:
            cleanup_cmd = cleanup_table['postmove_cleanup_cmds'][axis.axisid]
            if '.pos' not in cleanup_cmd:  # some postmove commands (e.g. limit seeks) force axis position to a particular value
                axis.pos += net_distance[axis.axisid]
            exec(cleanup_cmd)
        move_cmds = {}
        move_cmds['MOVE_CMD'] = [cleanup_table['orig_command']]
        move_cmds['MOVE_CMD'] += cleanup_table['auto_commands']
        has_finite_dist = [i for i in range(cleanup_table['nrows']) if cleanup_table['dT'][i] or cleanup_table['dP'][i]]
        for letter, number in zip(['T', 'P'], ['1', '2']):
            speeds = cleanup_table[f'speed_mode_{letter}']
            dists = cleanup_table[f'd{letter}']
            move_cmds[f'MOVE_VAL{number}'] = [f'{speeds[i]} {dists[i]:.3f}' for i in has_finite_dist]
        for key, sublist in move_cmds.items():
            string = pc.join_notes(*sublist)
            self.state.store(key, string)
        self.state.store('TOTAL_CRUISE_MOVES_T', self.state._val['TOTAL_CRUISE_MOVES_T'] + cleanup_table['TOTAL_CRUISE_MOVES_T'])
        self.state.store('TOTAL_CRUISE_MOVES_P', self.state._val['TOTAL_CRUISE_MOVES_P'] + cleanup_table['TOTAL_CRUISE_MOVES_P'])
        self.state.store('TOTAL_CREEP_MOVES_T', self.state._val['TOTAL_CREEP_MOVES_T'] + cleanup_table['TOTAL_CREEP_MOVES_T'])
        self.state.store('TOTAL_CREEP_MOVES_P', self.state._val['TOTAL_CREEP_MOVES_P'] + cleanup_table['TOTAL_CREEP_MOVES_P'])
        self.state.store('TOTAL_MOVE_SEQUENCES', self.state._val['TOTAL_MOVE_SEQUENCES'] + 1)
        self.state.store('LOG_NOTE', cleanup_table['log_note'])

class Axis(object):
    """Handler for a motion axis. Provides move syntax and keeps tracks of position.
    """

    def __init__(self, posmodel, axisid, printfunc=print):
        self.posmodel = posmodel
        self.axisid = axisid
        self.printfunc = printfunc
        self._load_cached_params()

    def _load_cached_params(self):
        # order of operations here matters
        self.antibacklash_final_move_dir = self._calc_antibacklash_final_move_dir()
        self.principle_hardstop_direction = self._calc_principle_hardstop_direction()
        self.backlash_clearance = self._calc_backlash_clearance()
        self.hardstop_clearance = self._calc_hardstop_clearance()
        self.hardstop_clearance_near_full_range = self._calc_hardstop_clearance_near_full_range()
        self.hardstop_debounce = self._calc_hardstop_debounce()
        motor_props = self.motor_calib_properties
        self.is_locked = motor_props['locked']
        self.signed_gear_ratio = motor_props['ccw_sign'] * motor_props['gear_ratio']
        if not self.is_locked:
            self.signed_gear_ratio /= motor_props['gear_calib'] # avoid divide-by-zero for disabled axis
        self._full_range = self._calc_full_range()
        self._debounced_range = self._calc_debounced_range()
        self._near_full_range = self._calc_near_full_range()

    @property
    def pos(self):
        """Internally-tracked angular position of the axis, at the output of the gear.
        """
        if self.axisid == pc.T:
            return self.posmodel.state._val['POS_T']
        else:
            return self.posmodel.state._val['POS_P']

    @pos.setter
    def pos(self, value):
        full_range = self.full_range
        if value < min(full_range) or value > max(full_range):
            self.printfunc(f'{self.posmodel.posid} axis {self.axisid}: cannot set pos to {value}' +
                  f' (outside allowed range of {full_range}). Keeping old value {self.pos}')
            return
        if self.axisid == pc.T:
            self.posmodel.state.store('POS_T', value)
        else:
            self.posmodel.state.store('POS_P', value)

    @property
    def maxpos(self):
        """Max accessible position, as defined by debounced_range.
        """
        return self.get_maxpos()  # this redirect is for legacy compatibility with external code

    def get_maxpos(self, use_near_full_range=False):
        """Function to return accessible position. By default this is within
        debounced_range. An option exists to get the max as defined by
        near_full_range instead.
        """
        d = self.near_full_range if use_near_full_range else self.debounced_range
        if self.last_primary_hardstop_dir >= 0:
            return max(d)
        else:
            return self.get_minpos(use_near_full_range) + (d[1] - d[0])

    @property
    def minpos(self):
        """Min accessible position, as defined by debounced_range.
        """
        return self.get_minpos()  # this redirect is for legacy compatibility with external code

    def get_minpos(self, use_near_full_range=False):
        """Function to return accessible position. By default this is within
        debounced_range. An option exists to get the min as defined by
        near_full_range instead.
        """
        d = self.near_full_range if use_near_full_range else self.debounced_range
        if self.last_primary_hardstop_dir < 0:
            return min(d)
        else:
            return self.get_maxpos(use_near_full_range) - (d[1] - d[0])

    @property
    def full_range(self):
        '''Returns 1x2 list, [min,max] values for full travel range.'''
        return self._full_range.copy()

    @property
    def debounced_range(self):
        '''Returns 1x2 list, [min,max] values for debounced travel range.'''
        return self._debounced_range.copy()

    @property
    def near_full_range(self):
        '''Returns 1x2 list, [min,max] values for a travel range between "full"
        and "debounced" which allows some small additional margin. In particular
        this is used for robustness when working with travel distances that have
        been discretized according to motor steps.
        '''
        return self._near_full_range.copy()

    @property
    def last_primary_hardstop_dir(self):
        if self.axisid == pc.T:
            return self.posmodel.state._val['LAST_PRIMARY_HARDSTOP_DIR_T']
        else:
            return self.posmodel.state._val['LAST_PRIMARY_HARDSTOP_DIR_P']

    @last_primary_hardstop_dir.setter
    def last_primary_hardstop_dir(self,value):
        if self.axisid == pc.T:
            self.posmodel.state.store('LAST_PRIMARY_HARDSTOP_DIR_T',value)
        else:
            self.posmodel.state.store('LAST_PRIMARY_HARDSTOP_DIR_P',value)

    @property
    def total_limit_seeks(self):
        if self.axisid == pc.T:
            return self.posmodel.state._val['TOTAL_LIMIT_SEEKS_T']
        else:
            return self.posmodel.state._val['TOTAL_LIMIT_SEEKS_P']

    @total_limit_seeks.setter
    def total_limit_seeks(self,value):
        if self.axisid == pc.T:
            self.posmodel.state.store('TOTAL_LIMIT_SEEKS_T',value)
        else:
            self.posmodel.state.store('TOTAL_LIMIT_SEEKS_P',value)

    @property
    def limit_seeking_search_distance(self):
        """A distance magnitude that guarantees hitting a hard limit in either direction.
        """
        full_range = self.full_range
        return abs((full_range[1]-full_range[0])*self.posmodel.state._val['LIMIT_SEEK_EXCEED_RANGE_FACTOR'])

    @property
    def motor_calib_properties(self):
        """Return properties for motor calibration.
        """
        prop = {}
        if self.axisid == pc.T:
            prop['gear_ratio'] = pc.gear_ratio[self.posmodel.state._val['GEAR_TYPE_T']]
            prop['ccw_sign'] = self.posmodel.state._val['MOTOR_CCW_DIR_T']
            prop['gear_calib'] = self.posmodel.state._val['GEAR_CALIB_T']
        else:
            prop['gear_ratio'] = pc.gear_ratio[self.posmodel.state._val['GEAR_TYPE_P']]
            prop['ccw_sign'] = self.posmodel.state._val['MOTOR_CCW_DIR_P']
            prop['gear_calib'] = self.posmodel.state._val['GEAR_CALIB_P']
        prop['locked'] = prop['gear_calib'] == 0.0
        return prop

    def motor_to_shaft(self, distance):
        """
        Convert a distance in motor angle to shaft angle at the gearbox output.
        """
        return distance / self.signed_gear_ratio

    def shaft_to_motor(self, distance):
        """
        Convert a distance in shaft angle to motor angle at the gearbox output.
        """
        return distance * self.signed_gear_ratio

    def truncate_to_limits(self, distance, start_pos=None, use_near_full_range=False):
        """Return distance after truncating it (if necessary) to the software limits.

        An expected starting position can optionally be argued. If None, then the
        internally-tracked position is used as the starting point.

        Can also optionally argue with use_near_full_range=True. This uses a slightly
        wider definiton for the max and min accessible limits at which to truncate.
        See the maxpos and minpos methods, as well as debounced_range vs near_full_range.
        """
        if distance == 0:
            return distance
        maxpos = self.get_maxpos(use_near_full_range)
        minpos = self.get_minpos(use_near_full_range)
        if start_pos == None:
            start_pos = self.pos
        target_pos = start_pos + distance
        if maxpos < minpos:
            distance = 0
        elif target_pos > maxpos:
            new_distance = maxpos - start_pos
            distance = new_distance
        elif target_pos < minpos:
            new_distance = minpos - start_pos
            distance = new_distance
        return distance

    def _calc_full_range(self):
        """Calculated from physical range only, with no subtraction of debounce
        distance.
        Returns [1x2] array of [min,max]
        """
        if self.axisid == pc.T:
            r = abs(self.posmodel.state._val['PHYSICAL_RANGE_T'])
            return [-0.50*r, 0.50*r]  # split theta range such that 0 is essentially in the middle
        else:
            r = abs(self.posmodel.state._val['PHYSICAL_RANGE_P'])
            return [185.0-r, 185.0]  # split phi range such that 0 is essentially at the minimum

    def _calc_debounced_range(self):
        """Calculated from full range, accounting for removal of both the hardstop
        clearance distances and the backlash removal distance.
        Returns [1x2] array of [min,max]
        """
        f = self.full_range
        h = self.hardstop_debounce # includes "backlash" and "clearance"
        return [f[0] + h[0], f[1] + h[1]]

    def _calc_near_full_range(self):
        """In-between full_range and debounced_range. This guarantees sufficient
        backlash clearance, while giving reasonable probability of not contacting
        the hard stop. Intended use case is to allow a small amount of numerical
        margin for making tiny cruise move corrections.
        Returns [1x2] array of [min,max]
        """
        f = self.full_range
        nh = self.hardstop_clearance_near_full_range # does not include "backlash"
        return [f[0] + nh[0], f[1] + nh[1]]

    def _calc_principle_hardstop_direction(self):
        """The "principle" hardstop is the one which is struck during homing.
        (The "secondary" hardstop is only struck when finding the total available travel range.)
        """
        if self.axisid == pc.T:
            return self.posmodel.state._val['PRINCIPLE_HARDSTOP_DIR_T']
        else:
            return self.posmodel.state._val['PRINCIPLE_HARDSTOP_DIR_P']

    def _calc_antibacklash_final_move_dir(self):
        if self.axisid == pc.T:
            return self.posmodel.state._val['ANTIBACKLASH_FINAL_MOVE_DIR_T']
        else:
            return self.posmodel.state._val['ANTIBACKLASH_FINAL_MOVE_DIR_P']

    def _calc_hardstop_debounce(self):
        """This is the amount to debounce off the hardstop after striking it.
        It is the hardstop clearance distance plus the backlash removal distance.
        Returns [1x2] array of [min,max]
        """
        h = self.hardstop_clearance
        b = self.backlash_clearance
        return [h[0] + b[0], h[1] + b[1]]

    def _calc_hardstop_clearance(self):
        """Minimum distance to stay clear from hardstop.
        Returns [1x2] array of [clearance_at_min_limit, clearance_at_max_limit].
        These are DIRECTIONAL quantities (i.e., the sign indicates the direction
        in which one would debounce after hitting the given hardstop, to get
        back into the accessible range).
        """
        if self.axisid == pc.T:
            if self.principle_hardstop_direction < 0:
                return [+self.posmodel.state._val['PRINCIPLE_HARDSTOP_CLEARANCE_T'],
                        -self.posmodel.state._val['SECONDARY_HARDSTOP_CLEARANCE_T']]
            else:
                return [+self.posmodel.state._val['SECONDARY_HARDSTOP_CLEARANCE_T'],
                        -self.posmodel.state._val['PRINCIPLE_HARDSTOP_CLEARANCE_T']]
        else:
            if self.principle_hardstop_direction < 0:
                return [+self.posmodel.state._val['PRINCIPLE_HARDSTOP_CLEARANCE_P'],
                        -self.posmodel.state._val['SECONDARY_HARDSTOP_CLEARANCE_P']]
            else:
                return [+self.posmodel.state._val['SECONDARY_HARDSTOP_CLEARANCE_P'],
                        -self.posmodel.state._val['PRINCIPLE_HARDSTOP_CLEARANCE_P']]

    def _calc_hardstop_clearance_near_full_range(self):
        """Slightly less clearance than provided by the _calc_hardstop_clearance
        function.
        Returns 1x2 array of [clearance_from_low_side_of_range, clearance_from_high_side].
        These are DIRECTIONAL quantities (i.e., the sign indicates the direction
        similarly as the hardstop clearance directions).
        """
        h = self.hardstop_clearance
        return [x * pc.near_full_range_reduced_hardstop_clearance_factor for x in h]

    def _calc_backlash_clearance(self):
        """Minimum clearance distance required for backlash removal moves.
        Returns 1x2 array of [clearance_from_low_side_of_range, clearance_from_high_side].
        These are DIRECTIONAL quantities (i.e., the sign indicates the direction
        similarly as the hardstop clearance directions).
        """
        if self.antibacklash_final_move_dir > 0:
            return [+self.posmodel.state._val['BACKLASH'],0]
        else:
            return [0,-self.posmodel.state._val['BACKLASH']]
