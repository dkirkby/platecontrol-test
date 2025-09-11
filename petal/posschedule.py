import posconstants as pc
import posschedulestage
import posschedstats
import time
import math

# enables debugging code
DEBUG = False

class PosSchedule(object):
    """Generates move table schedules in local (theta,phi) to get positioners
    from starts to finishes. The move tables are instances of the PosMoveTable
    class.

        petal   ... Instance of Petal that this schedule applies to.

        stats   ... Instance of PosSchedStats in which to register scheduling statistics.
                    If stats=None, then no statistics are logged.

        verbose ... Control verbosity at stdout.
    """

    def __init__(self, petal, stats=None, verbose=True):
        self.petal = petal
        if stats:
            schedule_id = pc.timestamp_str()
            self.stats = stats
            self.stats.register_new_schedule(schedule_id, len(self.petal.posids))
        else:
            self.stats = posschedstats.PosSchedStats(enabled=False) # this is really just to get the is_enabled() function available
        self.verbose = verbose
        self.printfunc = self.petal.printfunc
        self._requests = {} # keys: posids, values: target request dictionaries
        self.stage_order = ['direct', 'debounce_polygons', 'retract', 'rotate', 'extend', 'expert', 'final']
        self.RRE_stage_order = ['retract', 'rotate', 'extend']
        self.stages = {name:posschedulestage.PosScheduleStage(
                                collider         = self.collider,
                                stats            = self.stats,
                                power_supply_map = self.petal.power_supply_map,
                                verbose          = self.verbose,
                                printfunc        = self.printfunc,
                                petal            = self.petal
                            ) for name in self.stage_order}
        self.should_check_petal_boundaries = True # allows you to turn off petal-specific boundary checks for non-petal systems (such as positioner test stands)
        self.should_check_sweeps_continuity = False # if True, inspects all quantized sweeps to confirm well-formed. incurs slowdown, and generally is not needed; more for validating if any changes made to quantize function at a lower level
        self.move_tables = {}
        self.extra_log_notes = {} # keys = posid, values = strs --- special collection of extra log notes that should be stored, outside of the usual move_tables data tracking (e.g. for empty or motionless positioner special cases)
        self._expert_added_tables_sequence = []  # copy of original sequence in which expert tables were added (for error-recovery cases)
        self._all_requested_posids = {'regular': set(), 'expert': set()}  # every posid that received a request, whether accepted or not

    @property
    def collider(self):
        return self.petal.collider

    @property
    def regular_requests_accepted(self):
        return {posid for posid in self._requests if self.has_regular_request_already(posid)}

    @property
    def expert_requests_accepted(self):
        return set(self.stages['expert'].move_tables.keys())

    def has_regular_request_already(self, posid):
        if posid in self._requests and not self._requests[posid]['is_dummy']:
            return True
        return False

    def _reinit_stages(self):
        self.stages = {name:posschedulestage.PosScheduleStage(
                                collider         = self.collider,
                                stats            = self.stats,
                                power_supply_map = self.petal.power_supply_map,
                                verbose          = self.verbose,
                                printfunc        = self.printfunc
                            ) for name in self.stage_order}
        return

    def request_target(self, posid, uv_type, u, v, log_note='', allow_initial_interference=True):
        """Adds a request to the schedule for a given positioner to move to the
        target position (u,v) or by the target distance (du,dv) in the
        coordinate system indicated by uv_type.

              posid ... string, unique id of positioner
            uv_type ... string, valid arguments are:
                              ABSOLUTE       RELATIVE
                          ... 'QS'           'dQdS'
                          ... 'obsXY'        'obsdXdY'
                          ... 'poslocXY'     'poslocdXdY'
                          ... 'ptlXY'
                          ... 'posintTP'     'dTdP'
                          ... 'poslocTP'
                  u ... float, value of q, dq, x, dx, t, or dt
                  v ... float, value of s, ds, y, dy, p, or dp
           log_note ... optional string to store alongside the requested move
                        in the log data

        A schedule can only contain 1 target request per positioner at a time.

        The special argument allow_initial_interference takes a boolean. It allows
        us to immediately reject requests to positioners with intially overlapping
        keepout polygons. [JHS] As of 2020-10-29, I think it best in general to *not*
        reject such requests, and let the scheduling code handle the thornier questions
        like "are we barely overlapping, and should move anyway".

        Returns an error value == None if the request was accepted, or an
        explanatory string if the request was denied.
        """
        stats_enabled = self.stats.is_enabled()
        if stats_enabled:
            self.stats.add_request()
        self._all_requested_posids['regular'].add(posid)
        posmodel = self.petal.posmodels[posid]
        trans = posmodel.trans
        current_position = posmodel.expected_current_position
        start_posintTP = current_position['posintTP']
        cmd_target_str = self._make_coord_str(uv_type, [u, v], prefix='user')

        # options used below, for control of t_guess parameter in some coord conversions
        t_guess_OFF = None  # Using this option always puts target poslocP within [0, 180].
        t_guess_START = current_position['poslocTP'][pc.T]  # For "small" moves (where "small" means within
                                                            # t_guess_tol), this picks whichever of the possible
                                                            # poslocTP options that is closer to starting position.
                                                            # Here, poslocP may go outside [0, 180].

        # get into uniform coordinate system
        lims = 'targetable'
        unreachable = False
        if uv_type == 'QS':
            targt_posintTP, unreachable = trans.QS_to_posintTP([u, v], lims, t_guess=t_guess_OFF)
        elif uv_type == 'dQdS':
            start_uv = current_position['QS']
            targt_uv = posmodel.trans.addto_QS(start_uv, [u, v])
            targt_posintTP, unreachable = trans.QS_to_posintTP(targt_uv, lims, t_guess=t_guess_START)
        elif uv_type == 'poslocXY':
            targt_posintTP, unreachable = trans.poslocXY_to_posintTP([u, v], lims, t_guess=t_guess_OFF)
        elif uv_type == 'obsdXdY':
            # global cs5 projected xy, as returned by platemaker
            start_uv = current_position['obsXY']
            targt_uv = posmodel.trans.addto_XY(start_uv, [u, v])
            targt_posintTP, unreachable = posmodel.trans.obsXY_to_posintTP(targt_uv, lims, t_guess=t_guess_START)
        elif uv_type == 'poslocdXdY':
            # in poslocXY coordinates in local tangent plane, not global cs5
            start_uv = current_position['poslocXY']
            targt_uv = posmodel.trans.addto_XY(start_uv, [u, v])
            targt_posintTP, unreachable = posmodel.trans.poslocXY_to_posintTP(targt_uv, lims, t_guess=t_guess_START)
        elif uv_type == 'posintTP':
            targt_posintTP = [u, v]
        elif uv_type == 'dTdP':
            targt_posintTP = trans.addto_posintTP(start_posintTP, [u, v], lims)
        elif uv_type == 'obsXY':
            targt_posintTP, unreachable = trans.obsXY_to_posintTP([u, v], lims, t_guess=t_guess_OFF)
        elif uv_type == 'ptlXY':
            targt_posintTP, unreachable = trans.ptlXY_to_posintTP([u, v], lims, t_guess=t_guess_OFF)
        elif uv_type == 'poslocTP':
            targt_posintTP = trans.poslocTP_to_posintTP([u, v])
        else:
            return self._denied_str(cmd_target_str, 'Bad uv_type')

        # handle locked axes
        # 2021-04-29 [JHS] There may be more sophisiticated things one could do with the coord transformations,
        # to optimize closeness of single-DOF approach to requested target. Here I am currently opting for the
        # simplest approach: zero the target rotation angle for the locked axis.
        lock_msg = ''
        locks = posmodel.axis_locks
        if any(locks) and not all(locks): # all case is handled below in validations section
            locked_axis = 0 if locks[0] else 1
            locked_angle = start_posintTP[locked_axis]
            orig_targt_posintTP_str = self._make_coord_str('posintTP', targt_posintTP, prefix='user')
            lock_msg = pc.join_notes(orig_targt_posintTP_str, f'posint{"T" if locked_axis == 0 else "P"} locked at {locked_angle:.3f}')
            targt_posintTP_mutable = list(targt_posintTP)
            targt_posintTP_mutable[locked_axis] = start_posintTP[locked_axis]
            targt_posintTP = tuple(targt_posintTP_mutable)

        # other standard coordinates for validations and logging
        targt_poslocTP = trans.posintTP_to_poslocTP(targt_posintTP)
        targt_ptlXYZ = trans.poslocTP_to_ptlXYZ(targt_poslocTP)
        target_str_posintTP = self._make_coord_str('posintTP', targt_posintTP, prefix='req')
        target_str_ptlXYZ = self._make_coord_str('ptlXYZ', targt_ptlXYZ, prefix='req')
        target_str = pc.join_notes(lock_msg, target_str_posintTP, target_str_ptlXYZ)

        # validations
        if self._deny_request_because_both_locked(posmodel):
            return self._denied_str(target_str, BOTH_AXES_LOCKED_MSG)
        if unreachable:
            self.petal.pos_flags[posid] |= self.petal.flags.get('UNREACHABLE', self.petal.missing_flag)
            target_str = target_str.replace('req', 'nearest')
            target_str = pc.join_notes(cmd_target_str, target_str)
            return self._denied_str(target_str, 'Target not reachable.')
        if self.has_regular_request_already(posid):
            self.petal.pos_flags[posid] |= self.petal.flags.get('MULTIPLEREQUESTS', self.petal.missing_flag)
            return self._denied_str(target_str, 'Cannot request more than one target per positioner in a given schedule.')
        if self._deny_request_because_disabled(posmodel):
            return self._denied_str(target_str, POS_DISABLED_MSG)
        if self._deny_request_because_starting_out_of_range(posmodel):
            return self._denied_str(target_str, f'Bad initial position (POS_T, POS_P) = {posmodel.expected_current_posintTP} is' +
                                 f' outside allowed range T={posmodel.full_range_posintT} and/or P={posmodel.full_range_posintP}')
        if not allow_initial_interference:
            interfering_neighbors = self._check_init_or_final_neighbor_interference(posmodel)
            if interfering_neighbors:
                return self._denied_str(target_str, f'Interference at initial position with {interfering_neighbors}')
        limit_err = self._deny_request_because_limit(posmodel, targt_poslocTP)
        if limit_err:
            return self._denied_str(target_str, limit_err)
        if self.should_check_petal_boundaries:
            bounds_err = self._deny_request_because_out_of_bounds(posmodel, targt_poslocTP)
            if bounds_err:
                return self._denied_str(target_str, bounds_err)
        interfering_neighbors = self._check_init_or_final_neighbor_interference(posmodel, targt_poslocTP)
        if interfering_neighbors:
            return self._denied_str(target_str, f'Target interferes with existing target(s) of neighbors {interfering_neighbors}')
        if posmodel.is_linphi:
            targXY = trans.posintTP_to_poslocXY(targt_posintTP)
            strtXY = trans.posintTP_to_poslocXY(start_posintTP)
            dist_from_targt = 1000.0 * math.dist(targXY, strtXY)
            LINPHI_DIST_LIMIT = 10.0 # microns
            try:
                if hasattr(self.petal, 'petal_debug'):
                    LINPHI_DIST_LIMIT = float(self.petal.petal_debug.get('linphi_dist_limit'))
            except TypeError:
                pass
            if dist_from_targt < LINPHI_DIST_LIMIT: # 10 microns
                return self._denied_str(target_str, f"Linear phi already close enough, {dist_from_targt} < {LINPHI_DIST_LIMIT} microns to target")
        # form internal request dict
        new_request = {'start_posintTP': start_posintTP,
                       'targt_posintTP': targt_posintTP,
                       'posmodel': posmodel,
                       'posid': posid,
                       'command': uv_type,
                       'cmd_val1': u,
                       'cmd_val2': v,
                       'log_note': log_note,
                       'target_str': target_str,
                       'is_dummy': False,
                       }
        self._requests[posid] = new_request
        if stats_enabled:
            self.stats.add_request_accepted()
        return None

    def _schedule_moves(self, anticollision, should_anneal, scheduling_timer_start):
        if self.expert_mode_is_on():
            self._schedule_expert_tables(anticollision=anticollision, should_anneal=should_anneal)
        else:
            self._fill_enabled_but_nonmoving_with_dummy_requests()
            if anticollision == 'adjust':
                self._schedule_requests_with_path_adjustments(should_anneal=should_anneal, adjust_requested_only=False)
            elif anticollision == 'adjust_requested_only':
                self._schedule_requests_with_path_adjustments(should_anneal=should_anneal, adjust_requested_only=True)
            else:
                self._schedule_requests_with_no_path_adjustments(anticollision=anticollision, should_anneal=should_anneal)
        self._combine_stages_into_final()
        self.printfunc(f'Scheduling calculation done in {time.perf_counter()-scheduling_timer_start:.3f} sec')
        finalcheck_timer_start = time.perf_counter()
        final = self.stages['final']
        if anticollision:
            c, _, p = self._check_final_stage(msg_prefix='Penultimate')
        else:
            c, _, p = self._check_final_stage(msg_prefix='Final')
        colliding_sweeps, collision_pairs = c, p # for readability
        if anticollision:
            if not colliding_sweeps:
                self.printfunc('Final collision check --> skipped (because \'penultimate\' check already succeeded)')
            else:
                adjusted = set()
                frozen = set()
                for posid in colliding_sweeps:
                    these_adjusted, these_frozen = final.adjust_path(posid, freezing='forced_recursive')
                    adjusted.update(these_adjusted)
                    frozen.update(these_frozen)
                prefix = 'Following from \'penultimate\' check:'
                self.printfunc(f'{prefix} adjusted posids --> {adjusted}')
                self.printfunc(f'{prefix} frozen posids --> {frozen}')
                c, _, p = self._check_final_stage(msg_prefix='Final',
                                                  msg_suffix=' (should always be zero)',
                                                  assert_no_unresolved=False)   # Cheanged to False 07/15/2024 cad - now handled by self.schedule_moves()
                colliding_sweeps, collision_pairs = c, p # for readability
        return colliding_sweeps, collision_pairs, finalcheck_timer_start, final

    def _handle_schedule_moves_collision(self, colliding_sweeps, collision_pairs):
        """
        This function handles the case when, after all the collision mitigation strategies have run, there are still collisions
        in the planned moves.

        If some of the colliding posids are from zeno devices, replace the target for those device(s) with a dummy,
        and also temporarily disable that device, then when this function exits, the plan will be re-run (takes 2-4 sec at KPNO).

        If none of the colliding posids are from zeno devices, then resolve_non_zeno determines what happens.

        If resolve_non_zeno is False, then this code does what it originally did, which is to print error messages, then cause an exception.

        If resolve_non_zeno is True, then all colliding devices will have targets replaced with dummies and they will be temporarily disabled.
        Again, after this function exits, the planning for the remaining devices will be re-run.
        """
        zeno_posids = set()
        colliding_posids = [posid for posid in colliding_sweeps]
        for posid in colliding_posids:
            p_state = self.petal.posmodels[posid].state
            if p_state._val['ZENO_MOTOR_P'] is True:
                zeno_posids.add(posid)
        if zeno_posids:
            colliding = set(colliding_sweeps)
            self.printfunc(self.get_details_str(colliding, label=f'Unresolved zeno collision avoided: removed target(s) for {zeno_posids}'))
            for psid in zeno_posids:
                self._make_dummy_request(psid, lognote='target removed due to collision avoidance failure')
#           self.petal.temporary_disable_positioners_reason(zeno_posids,'collision avoidance failure')  # Disable the involved zeno motors so they won't be used until fp_setup is run, likely saves move planning time on subsequent moves
            stats_enabled = self.stats.is_enabled()
            if stats_enabled:
                self.stats.sub_request_accepted()
            self._reinit_stages() # clear out old move tables - starting over
        else:
            colliding = set(colliding_sweeps)
            self.printfunc(self.get_details_str(colliding, label='colliding'))
            resolve_non_zeno = True
            if not resolve_non_zeno:  # blow up PETAL on purpose to call attention to it
                err_str = f'{len(colliding)} collisions were NOT resolved! This indicates a bug that needs to be fixed. See details above.'
                self.printfunc(err_str)
                # self.petal.enter_pdb()  # 2023-07-11 [CAD] This line causes a fault (no function enter_pdb)
                assert False, err_str  # 2020-11-16 [JHS] put a PDB entry point in rather than assert, so I can inspect memory next time this happens online
            else: # temporarily disable the offending robots and press on
                self.printfunc(self.get_details_str(colliding, label=f'Unresolved collision avoided: removed target(s) for {colliding_posids}'))
                for psid in colliding_posids:
                    self._make_dummy_request(psid, lognote='target removed due to collision avoidance failure')
#               self.petal.temporary_disable_positioners_reason(colliding_posids,'collision avoidance failure')  # Disable the involved motors so they won't be used until fp_setup is run, likely saves move planning time on subsequent moves
                stats_enabled = self.stats.is_enabled()
                if stats_enabled:
                    self.stats.sub_request_accepted()
                self._reinit_stages() # clear out old move tables - starting over
        return

    def schedule_moves(self, anticollision='freeze', should_anneal=True):
        """Executes the scheduling algorithm upon the stored list of move requests.

        A single move table is generated for each positioner that has a request
        registered. The resulting tables are stored in the move_tables list.

        There are three options for anticollision behavior during scheduling:

          None      ... Expert use only.

          'freeze'  ... If any collisions are found, the colliding positioner
                        is frozen prior to the requested target position. This
                        setting is suitable for small correction moves.

          'adjust'  ... If any collisions are found, the motion paths of the
                        colliding positioners are adjusted to attempt to avoid
                        each other. If this fails, the colliding positioner
                        is frozen at its original position. This setting is
                        suitable for gross retargeting moves.

                        Occasionally, there is a neighbor positioner with no
                        requested target in the desired path of one with a target.
                        In this case, the 'adjust' algorithm is allowed to move
                        that neighbor out of the way, and then restore it back
                        to its starting position.

          'adjust_requested_only'
                    ... Same as 'adjust', however no automatic
                        "get out of the way" moves will be done on unrequested
                        neighbors.

        If there were ANY pre-existing move tables in the list (for example, hard-
        stop seeking tables directly added by an expert user or expert function),
        then the requests list is ignored. The only changes to move tables are
        for power density annealing. Furthermore, if anticollision='adjust',
        then it reverts to 'freeze' instead. An argument of anticollision=None
        remains as-is.

        The boolean flag should_anneal controls whether or not to spread out
        the move density in time.
        """

        self._schedule_moves_initialize_logging(anticollision)
        if not self._requests and not self.expert_mode_is_on():
            self.printfunc('No requests nor existing move tables found. No move scheduling performed.')
            return
        all_accepted = set()
        for kind in ['regular', 'expert']:
            received = self._all_requested_posids[kind]
            accepted = self.regular_requests_accepted if kind == 'regular' else self.get_posids_with_expert_tables()
            rejected = received - accepted
            if len(received) > 0:
                prefix = f'num {kind} target requests'
                self.printfunc(f'{prefix} received = {len(received)}')
                self.printfunc(f'{prefix} accepted = {len(accepted)}')
                self.printfunc(f'{prefix} rejected = {len(rejected)}')
                if rejected:
                    self.printfunc(f'pos with rejected {kind} request(s): {rejected}')
            all_accepted |= accepted

        # If anticollision mode is None or 'freeze' and number of targets is < limit
        # then there is no need for annealing.  In particular, this will have significant
        # gains in the speed at which fp_setup runs with no additional chance of collisions
        num_targets = len(all_accepted)
        if anticollision in {None, 'freeze'} and num_targets <= 2*pc.max_targets_for_no_anneal:
            count = {'V1': 0, 'V2':0}
            for p in all_accepted:
                for ps in ['V1', 'V2']:
                    if p in self.petal.power_supply_map[ps]:
                        count[ps] += 1
            if count['V1'] <= pc.max_targets_for_no_anneal and count['V2'] <= pc.max_targets_for_no_anneal:
                if hasattr(self.petal, 'petal_debug') and self.petal.petal_debug.get('cancel_anneal_verbose') and should_anneal:
                    self.printfunc(f'Annealing cancelled due to anticollision={anticollision} and number of targets={count} <= max of {pc.max_targets_for_no_anneal}')
                should_anneal = False


        scheduling_timer_start = time.perf_counter()
        do_schedule = True

        while do_schedule:
            colliding_sweeps, collision_pairs, finalcheck_timer_start, final = \
                self._schedule_moves(anticollision, should_anneal, scheduling_timer_start)
#            if DEBUG:
#                colliding_sweeps, collision_pairs = self._possibly_induce_scheduling_error(colliding_sweeps, collision_pairs, anticollision)
            if not collision_pairs or not anticollision:
                do_schedule = False
            else:
                self._handle_schedule_moves_collision(colliding_sweeps, collision_pairs)

        self.printfunc(f'Final collision checks done in {time.perf_counter()-finalcheck_timer_start:.3f} sec')
        self._schedule_moves_check_final_sweeps_continuity()
        self._schedule_moves_store_collisions_and_pairs(colliding_sweeps, collision_pairs)
        self.move_tables = final.rewrite_zeno_move_tables(final.move_tables) # Apply Zeno mods AFTER normal scheduling and anticollision checks -- only possible iff extra moves are well within keepouts
        empties = {posid for posid, table in self.move_tables.items() if not table}
        motionless = {posid for posid, table in self.move_tables.items() if table.is_motionless}
        for posid in empties | motionless:
            # Ignore expert tables because they don't have requests
            if posid in (all_accepted-self.get_posids_with_expert_tables()):
                req = self._requests[posid]
                note = pc.join_notes(req['log_note'], req['target_str']) # because the move tables for these are about to be deleted
                self.extra_log_notes[posid] = note
            del self.move_tables[posid]
        if self.petal.animator_on:
            anim_tables = {posid: table.copy() for posid, table in self.move_tables.items()}
        else:
            anim_tables = {}
        for table in self.move_tables.values():
            table.strip()
        self._schedule_moves_store_requests_info()
        self._schedule_moves_finish_logging(anim_tables)

    def conservative_move_timeout_period(self, safety_factor=4.0):
        """Returns a conservative period of time (in seconds) that one should
        wait for move table execution, before assuming some failure must have
        occurred. An internal time estimate is done. Then this estimate is
        multiplied by the argued safety_factor, and returned.
        """
        times = {table.total_time(suppress_automoves=False) for table in self.move_tables.values()}
        if times:
            return max(times) * safety_factor
        return 0.0

    def expert_add_table(self, move_table):
        """Adds an externally-constructed move table to the schedule. Only simple
        freezing is available as an anticollision method for such externally-made
        tables. If any tables have been added by this method, then any target requests
        will be ignored upon scheduling. Generally, this method should only be used
        by an expert user.

        Returns an error string or None if no error.
        """
        stats_enabled = self.stats.is_enabled()
        if stats_enabled:
            timer_start = time.perf_counter()
        self._all_requested_posids['expert'].add(move_table.posid)
        disabled = self._deny_request_because_disabled(move_table.posmodel)
        both_locked = self._deny_request_because_both_locked(move_table.posmodel)
        if disabled or both_locked:
            sched_table = move_table.for_schedule()
            net_dtdp = [sched_table[key][-1] for key in ['net_dT', 'net_dP']]
            target_str = self._make_coord_str('expert_net_dtdp', net_dtdp, prefix='')
            msg = POS_DISABLED_MSG if disabled else BOTH_AXES_LOCKED_MSG
            return self._denied_str(target_str, msg)
        self.stages['expert'].add_table(move_table)
        self._expert_added_tables_sequence.append(move_table.copy())
        if stats_enabled:
            total_time = time.perf_counter() - timer_start
            self.stats.add_expert_table_time(total_time)
        return None

    def expert_mode_is_on(self):
        """Returns boolean stating whether scheduling is in expert mode. This is
        the case if any calls have been made to expert_add_table(). See that
        function's comments for more detail.
        """
        return self.stages['expert'].is_not_empty()

    def get_requests(self, include_dummies=False):
        """Returns a dict containing copies of all the current requests. Keys
        are posids. Any dummy requests (auto-generated during anticollision
        scheduling) by default are excluded. But can be included with the
        include_dummies boolean.
        """
        requests = {}
        for posid, request in self._requests.items():
            if not request['is_dummy'] or include_dummies:
                requests[posid] = request.copy()
        return requests

    def get_orig_expert_tables_sequence(self):
        """Returns a list containing any tables added by the "expert_add_table()"
        function, in the same order they were added.
        """
        return self._expert_added_tables_sequence

    def get_frozen_posids(self):
        '''Returns set of any posids whose final move tables do not achieve the
        targets in their requests. This is *slower* and *more complete* than simply
        checking the "is_frozen" property of sweeps (which only really work within
        a substage --- may not capture the effects of multiple sequential stages).
        Intended to be called *after* doing schedule_moves().
        '''
        frozen = set()
        SHOW_FROZEN_MT = None
        if hasattr(self.petal, 'petal_debug'):
            SHOW_FROZEN_MT = self.petal.petal_debug.get('show_frozen_mt')
        user_requested = set(self.get_requests(include_dummies=False))
        has_table = set(self.move_tables)
        check = has_table & user_requested # ignores expert tables
        frozen |= user_requested - has_table
        err = {}
        for posid in check:
            request = self._requests[posid]
            sched_table = self.move_tables[posid].for_schedule()
            if SHOW_FROZEN_MT:
                if posid in SHOW_FROZEN_MT:
                    ht = self.move_tables[posid].for_hardware()
                    self.printfunc(f"posid={posid}\nhardware_movetable={str(ht)}\nschedule_movetable={str(sched_table)}")
            net_requested = [request['targt_posintTP'][i] - request['start_posintTP'][i] for i in [0,1]]
            net_scheduled = [sched_table['net_dT'][-1], sched_table['net_dP'][-1]]
            err[posid] = [abs(net_requested[i] - net_scheduled[i]) for i in [0, 1]]
            err[posid] = [min(e, abs(e - 360)) for e in err[posid]]  # wrapped angle cases
            do_set_frozen = False
            if SHOW_FROZEN_MT and posid in SHOW_FROZEN_MT:
                self.printfunc(f'posid={posid}, net_requested={net_requested}, net_scheduled={net_scheduled}, err={err[posid]}')
            if self.petal.posmodels[posid].is_linphi:
                if err[posid][0] > pc.schedule_checking_numeric_angular_tol or\
                   err[posid][1] > pc.schedule_checking_angular_tol_zeno:
                       do_set_frozen = True
            else:
                if max(err[posid]) > pc.schedule_checking_numeric_angular_tol:
                    do_set_frozen = True
            if do_set_frozen:
                frozen.add(posid)
                if SHOW_FROZEN_MT and posid in SHOW_FROZEN_MT:
                    self.printfunc(f'posid {posid} is frozen')
            else:
                if SHOW_FROZEN_MT and posid in SHOW_FROZEN_MT:
                    self.printfunc(f'posid {posid} is not frozen')

        return frozen

    def get_posids_with_expert_tables(self):
        '''Returns set of any posids for which an "expert" table has been added.
        '''
        return {table.posid for table in self._expert_added_tables_sequence}

    def get_overlaps(self, posids):
        '''Returns dict with keys=posids and values=sets of any neighbors with
        overlapping polygons (in the current, initial configuration). Positioners
        with no overlaps are excluded from the returned dict.'''
        overlaps = {}
        for posid in sorted(posids):
            model = self.collider.posmodels[posid]
            these = self._check_init_or_final_neighbor_interference(model, final_poslocTP=None)
            if these:
                overlaps[posid] = these
        return overlaps

    def get_details_str(self, posids, label=''):
        '''Returns string describing details of move tables and sweeps for
        positioners in the final stage. The argued posids should be a set.
        '''
        if pc.is_string(posids):
            posids = {posids}
        else:
            posids = set(posids)
        neighbors = set()
        for posid in posids:
            neighbors |= self.collider.pos_neighbors[posid]
        neighbors -= posids
        label = label if label else 'argued posids'
        s = f'{label}: {posids}\n'
        s += f'neighbors: {neighbors}\n'
        display_posids = (posids | neighbors) & set(self.move_tables)
        for posid in display_posids:
            suffix = label.upper() if posid in posids else 'NEIGHBOR'
            s += f'\n---- {posid} [{suffix}] ----'
            self.move_tables[posid].for_hardware() # dummy call, to ensure inclusion of final creep rows in next display line
            s += '\n' + self.move_tables[posid].display(printfunc=None)
            s += '\n' + self.move_tables[posid].display_for('schedule', printfunc=None)
            s += '\n' + self.move_tables[posid].display_for('hardware', printfunc=None)
            s += f'\nsweep for: {posid}'
            s += '\n' + str(self.stages['final'].sweeps[posid])
            s += '\n'
        return s

    def plot_density(self, path=None):
        '''Bins and plots total power density of motors as-scheduled. Useful for
        checking effects of annealing. Assumes that sweeps have been calculated
        already and stored (c.f. store_collision_finding_results).'''
        import matplotlib.pyplot as plt
        plt.switch_backend('Agg')
        import poscollider
        import os
        plt.ioff()
        quantize_timestep = 0.02 # sec
        plt.figure()
        for supply, posids in self.petal.power_supply_map.items():
            has_table = posids & set(self.move_tables)
            if not any(has_table):
                continue
            num_moving = []
            for posid in has_table:
                table = self.move_tables[posid]
                sweep = poscollider.PosSweep(posid)
                sweep.fill_exact(init_poslocTP=table.init_poslocTP,
                                 table=table.for_collider(suppress_automoves=False),
                                 start_time=0)
                sweep.quantize(timestep=quantize_timestep)
                was_moving = {axis: [sweep.axis_was_moving(i, axis) for i in range(len(sweep.time))] for axis in [pc.T, pc.P]}
                for axis, axis_was_moving in was_moving.items():
                    longer = len(axis_was_moving) - len(num_moving)
                    num_moving.extend([0] * longer)
                    for i in range(len(axis_was_moving)):
                        if axis_was_moving[i]:
                            num_moving[i] += 1
            time = [i * quantize_timestep for i in range(len(num_moving))]
            plt.plot(time, num_moving, label=f'Power supply: {supply}')
        plt.xlabel('time [sec]')
        plt.ylabel('num motors moving')
        plt.title(f'move schedule density - petal id {self.petal.petal_id}\n{pc.timestamp_str()}')
        plt.legend()
        if not path:
            path = pc.dirs['temp_files']
            path = os.path.join(path, f'density_ptlid{self.petal.petal_id:02}_{pc.filename_timestamp_str()}.png')
        plt.tight_layout()
        plt.savefig(path)
        plt.close()
        self.printfunc(f'Saved density plot to {path}')

    def _schedule_expert_tables(self, anticollision, should_anneal):
        """Gathers data from expert-added move tables and populates the 'expert'
        stage. Any move requests are ignored.
        """
        should_freeze = not(not(anticollision))
        if should_freeze:
            if anticollision != 'freeze':
                self.printfunc(f'anticollision method \'{anticollision}\' overriden with \'freeze\', due to presence of expert move tables')
            if self.stats.is_enabled():
                self.stats.set_scheduling_method('freeze')
                self.stats.add_note('expert tables')
        self._direct_stage_conditioning(stage=self.stages['expert'],
                                        should_freeze=should_freeze,
                                        should_anneal=should_anneal)

    def _schedule_requests_with_no_path_adjustments(self, anticollision, should_anneal):
        """Gathers data from requests dictionary and populates the 'direct'
        stage with direct motions from start to finish. The positioners are
        given no path adjustments to avoid each other.
        """
        start_posintTP = {}
        desired_final_posintTP = {}
        dtdp = {}
        for posid, request in self._requests.items():
            start_posintTP[posid] = request['start_posintTP']
            desired_final_posintTP[posid] = request['targt_posintTP']
            trans = self.petal.posmodels[posid].trans
            dtdp[posid] = trans.delta_posintTP(desired_final_posintTP[posid],
                                               start_posintTP[posid],
                                               range_wrap_limits='targetable')
        stage = self.stages['direct']
        stage.initialize_move_tables(start_posintTP, dtdp)
#       stage.move_tables = stage.rewrite_zeno_move_tables(stage.move_tables)
        # double-negative syntax is to be compatible with various
        # False/None/'' negative values
        should_freeze = not(not(anticollision))
        self._direct_stage_conditioning(stage=stage,
                                        should_freeze=should_freeze,
                                        should_anneal=should_anneal)

    def _direct_stage_conditioning(self, stage, should_freeze, should_anneal):
        """Applies annealing and possibly freezing to a 'direct' or 'expert' stage.

            stage         ... instance of PosScheduleStage, needs to already have its move tables initialized
            should_freeze ... boolean, says whether to check for collisions and freeze
            should_anneal ... boolean, enables/disables annealing
        """
        if should_anneal:
            stage.anneal_tables(suppress_automoves=False, mode=self.petal.anneal_mode)
        if should_freeze:
            colliding_sweeps, all_sweeps = stage.find_collisions(stage.move_tables)
            stage.store_collision_finding_results(colliding_sweeps, all_sweeps)
            if self.verbose:
                self.printfunc("initial stage.colliding: " + str(stage.colliding))
            adjustment_performed = False
            for posid in sorted(stage.colliding.copy()): # sort is for repeatability (since stage.colliding is an unordered set, and so path adjustments would otherwise get processed in variable order from run to run). the copy() call is redundant with sorted(), but left there for the sake of clarity, that need to be looping on a copy of *some* kind
                if posid in stage.colliding: # re-check, since earlier path adjustments in loop may have already resolved this posid's collision
                    adjusted, frozen = stage.adjust_path(posid, freezing='forced_recursive')
                    adjustment_performed = True
                    if self.verbose:
                        self.printfunc("remaining stage.colliding " + str(stage.colliding))
            if self.stats.is_enabled() and adjustment_performed:
                self.stats.add_to_num_adjustment_iters(1)

    def _debounce_polygons(self):
        '''Looks for cases where positioners have initial overlap with neighbors.
        In cases where the polygons are just barely touching, this function attempts
        to schedule initial small moves(s) to get the polygons free of one another.
        These moves go into the 'debounce_polygons' stage. Returns a dict with
        keys = posid and values = any adjusted starting POS_T, POS_P (so that later
        stages know where this debounce move puts the robots.)
        '''
        user_requests = self.get_requests(include_dummies=False)
        requested_posids = user_requests.keys()
        overlaps_dict = self.get_overlaps(requested_posids)
        if len(overlaps_dict) == 0:
            return {}
        overlapping = set(overlaps_dict)
        for overlapping_neighbors in overlaps_dict.values():
            overlapping |= overlapping_neighbors
        overlapping -= pc.case.fixed_case_names
        stage = self.stages['debounce_polygons']
        skip = pc.num_timesteps_ignore_overlap
        db = pc.debounce_polys_distance
        models = {posid: self.collider.posmodels[posid] for posid in overlapping}
        spinupdown_distances  = {model.abs_shaft_spinupdown_distance_T for model in models.values()}
        spinupdown_distances |= {model.abs_shaft_spinupdown_distance_P for model in models.values()}
        backlash_distances  = {model.state._val['BACKLASH'] for model in models.values()}
        err_msg_prefix = f'posschedule.py: posconstants.debounce_polys_distance = {db:.3f} is insufficient'
        assert db >= 2*max(spinupdown_distances), f'{err_msg_prefix}, < 2*max(spinupdown) = {2*max(spinupdown_distances):.3f}'
        assert db >= max(backlash_distances), f'{err_msg_prefix}, < max(backlash) = {max(backlash_distances):.3f}'
        delta_options = [(0, db), (db, 0), (-db, 0), (db, db), (-db, db)]
        # delta_options = [(0,-db)]  # uncomment this line for debugging ONLY, to induce more likely failures of debouncing for close polygons
        enabled = self.petal.all_enabled_posids()
        unresolved = overlapping
        for deltas in delta_options:
            if not unresolved:
                break
            posids_to_adjust = unresolved & enabled & overlapping
            start_tp = {posid: self._requests[posid]['start_posintTP'] for posid in posids_to_adjust}
            dtdp = {posid: deltas for posid in start_tp}

            # [JHS] Comments on use of stage here:
            # 1. Each time we initialize_move_tables, only updating the ones for which a new delta is proposed.
            # 2. No annealing allowed! (would mess up the "skip first x timesteps" during collision checking)
            stage.initialize_move_tables(start_tp, dtdp, update_only=True)
#           stage.rewrite_zeno_move_tables(stage.move_tables)  # i.e. a for loop of the single motor similar function
            colliding_sweeps, all_sweeps = stage.find_collisions(stage.move_tables, skip=skip)
            unresolved = set(colliding_sweeps)
        adjustments_failed = unresolved & enabled & overlapping
        if any(adjustments_failed):
            for posid in adjustments_failed:
                stage.del_table(posid)
                if posid in requested_posids & overlapping:
                    req = user_requests[posid]
                    explanation = f'Interference at initial position with {overlaps_dict[posid]}. Could ' + \
                                  f'not resolve within {skip} timestep{"s" if skip != 1 else ""} ({skip*self.collider.timestep:.3f} ' + \
                                  f'sec) by debouncing polygons (jog distances = {db} deg).'
                    deny_msg = self._denied_str(req['target_str'], explanation)
                    self.petal._print_and_store_note(posid, deny_msg)  # since deleting request, must get into log_note now
                del self._requests[posid]
            self._fill_enabled_but_nonmoving_with_dummy_requests()  # to repopulate deletions
        resolved = overlapping - set(colliding_sweeps)
        resolved_overlaps_dict = self.get_overlaps(resolved)
        final_posintTP = {}
        for posid in resolved & enabled:
            sched_table = stage.move_tables[posid].for_schedule()
            dT = sched_table["net_dT"][-1]
            dP = sched_table["net_dP"][-1]
            msg = f'Debounced initial polygon overlap with {resolved_overlaps_dict[posid]} using dtdp=({dT:.3f}, {dP:.3f})'
            self._requests[posid]['log_note'] = pc.join_notes(self._requests[posid]['log_note'], msg)
            self.printfunc(f'{posid}: {msg}')
            final_posintTP[posid] = (stage.start_posintTP[posid][0] + dT,
                                     stage.start_posintTP[posid][1] + dP)
        return final_posintTP

    def _schedule_requests_with_path_adjustments(self, should_anneal=True, adjust_requested_only=False):
        """Gathers data from requests dictionary and populates the 'retract',
        'rotate', and 'extend' stages with motion paths from start to finish.
        The move tables may include adjustments of paths to avoid collisions.
        """
        debounced_start_posintTP = self._debounce_polygons()
        stats_enabled = self.stats.is_enabled()
        start_posintTP = {name: {} for name in self.RRE_stage_order}
        desired_final_posintTP = {name: {} for name in self.RRE_stage_order}
        dtdp = {name: {} for name in self.RRE_stage_order}
        retracted_poslocP = self.collider.Eo_phi  # Ei would also be safe, but unnecessary in most cases. Costs more time and power to get to
        no_auto_adjust = set()
        if adjust_requested_only:
            no_auto_adjust = {posid for posid, req in self._requests.items() if req['is_dummy']}
        for posid, request in self._requests.items():
            # Some care is taken here to use only delta and add functions
            # provided by PosTransforms to ensure that range wrap limits are
            # always safely handled from stage to stage.
            posmodel = self.petal.posmodels[posid]
            trans = posmodel.trans
            if posid in debounced_start_posintTP:
                this_start_posintTP = debounced_start_posintTP[posid]
            else:
                this_start_posintTP = request['start_posintTP']
            this_start_poslocTP = trans.posintTP_to_poslocTP(this_start_posintTP)
            start_posintTP['retract'][posid] = this_start_posintTP
            if this_start_poslocTP[pc.P] > self.collider.Eo_phi or request['start_posintTP'] == request['targt_posintTP']:
                retracted_posintP = start_posintTP['retract'][posid][pc.P]
            else:
                retracted_posintTP = trans.poslocTP_to_posintTP([0, retracted_poslocP])  # poslocT=0 is a dummy value
                retracted_posintP = retracted_posintTP[pc.P]
            desired_final_posintTP['retract'][posid] = [request['start_posintTP'][pc.T], retracted_posintP]
            desired_final_posintTP['rotate'][posid] = [request['targt_posintTP'][pc.T], retracted_posintP]
            desired_final_posintTP['extend'][posid] = request['targt_posintTP']
            def calc_dtdp(name, posid):
                tp_start = start_posintTP[name][posid]
                tp_final = desired_final_posintTP[name][posid]
                return trans.delta_posintTP(tp_final, tp_start, range_wrap_limits='targetable')
            def calc_next_tp(last_name, posid):
                last_tp = start_posintTP[last_name][posid]
                this_dtdp = dtdp[last_name][posid]
                return trans.addto_posintTP(last_tp, this_dtdp, range_wrap_limits='targetable')
            dtdp['retract'][posid] = calc_dtdp('retract', posid)
            start_posintTP['rotate'][posid] = calc_next_tp('retract', posid)
            dtdp['rotate'][posid] = calc_dtdp('rotate', posid)
            start_posintTP['extend'][posid] = calc_next_tp('rotate', posid)
            dtdp['extend'][posid] = calc_dtdp('extend', posid)
        for name in self.RRE_stage_order:
            stage = self.stages[name]
            stage.initialize_move_tables(start_posintTP[name], dtdp[name])
#           stage.move_tables = stage.rewrite_zeno_move_tables(stage.move_tables)
            not_the_last_stage = name != self.RRE_stage_order[-1]
            if should_anneal:
                stage.anneal_tables(suppress_automoves=not_the_last_stage, mode=self.petal.anneal_mode)
            if self.verbose:
                self.printfunc(f'posschedule: finding collisions for {len(stage.move_tables)} positioners, trying {name}')
                self.printfunc('Posschedule first move table: \n' + str(list(stage.move_tables.values())[0].for_collider()))
            colliding_sweeps, all_sweeps = stage.find_collisions(stage.move_tables)
            stage.store_collision_finding_results(colliding_sweeps, all_sweeps)
            attempts_sequence = ['off','on','forced','forced_recursive'] # these are used as freezing arg to adjust_path()
            while stage.colliding and attempts_sequence:
                freezing = attempts_sequence.pop(0)
                for posid in sorted(stage.colliding.copy()): # sort is for repeatability (since stage.colliding is an unordered set, and so path adjustments would otherwise get processed in variable order from run to run). the copy() call is redundant with sorted(), but left there for the sake of clarity, that need to be looping on a copy of *some* kind
                    if posid in stage.colliding: # because it may have been resolved already when a *neighbor* got previously adjusted
                        adjusted, frozen = stage.adjust_path(posid, freezing=freezing, do_not_move=no_auto_adjust)
                        for p in frozen:
                            if not_the_last_stage: # i.e. some next stage exists
                                # must set next stage to begin from the newly-frozen position
                                adjusted_table_data = stage.move_tables[p].for_schedule()
                                adjusted_t = start_posintTP[name][p][pc.T] + adjusted_table_data['net_dT'][-1]
                                adjusted_p = start_posintTP[name][p][pc.P] + adjusted_table_data['net_dP'][-1]
                                next_stage_idx = self.RRE_stage_order.index(name) + 1
                                next_name = self.RRE_stage_order[next_stage_idx]
                                start_posintTP[next_name][p] = [adjusted_t,adjusted_p]
                                dtdp[next_name][p] = calc_dtdp(next_name, p)
                if stats_enabled:
                    self.stats.add_to_num_adjustment_iters(1)
            if stage.colliding:
                self.printfunc('Error: During ' + name.upper() + ' stage of move scheduling (see PosSchedule.py), the positioners ' + str([posid for posid in stage.colliding]) + ' had collision(s) that were NOT resolved. This means there is a bug somewhere in the code that needs to be found and fixed. If this move is executed on hardware, these two positioners will collide!')
                self.printfunc('The move table(s) for these are:')
                for posid in stage.colliding:
                    for n in self.RRE_stage_order:
                        stage_str = str(posid) + ': ' + n.upper()
                        if posid in self.stages[n].move_tables:
                            self.printfunc(stage_str)
                            self.stages[n].move_tables[posid].display(self.printfunc)
                        elif n == name:
                            self.printfunc(stage_str + ' --> no move table found')
                if stats_enabled:
                    sorted_colliding = sorted(stage.colliding) # just for human ease of reading the values
                    colliding_tables = {posid:stage.move_tables[posid] for posid in sorted_colliding}
                    colliding_sweeps = {posid:stage.sweeps[posid] for posid in sorted_colliding}
                    self.stats.add_unresolved_colliding_at_stage(name, sorted_colliding, colliding_tables, colliding_sweeps)

    def _make_dummy_request(self, posid, lognote='generated by path adjustment scheduler for enabled but untargeted positioner'):
        posmodel = self.petal.posmodels[posid]
        current_posintTP = posmodel.expected_current_posintTP
        new_request = {'start_posintTP': current_posintTP,
                       'targt_posintTP': current_posintTP,
                       'posmodel': posmodel,
                       'posid': posid,
                       'command': '(autogenerated)',
                       'cmd_val1': 0,
                       'cmd_val2': 0,
                       'log_note': lognote,
                       'target_str': '',
                       'is_dummy': True,
                       }
        self._requests[posid] = new_request

    def _fill_enabled_but_nonmoving_with_dummy_requests(self):
        enabled = set(self.petal.all_enabled_posids())
        requested = set(self._requests.keys())
        enabled_but_not_requested = enabled - requested
        for posid in enabled_but_not_requested:
            self._make_dummy_request(posid)

    def _deny_request_because_disabled(self, posmodel):
        """Checks enabled status and includes setting flag.
        """
        enabled = posmodel.is_enabled
        if enabled == False:  # this is specifically NOT worded as "if not enabled:", because here we actually do not want a value of None to pass the test, in case the parameter field 'CTRL_ENABLED' has not yet been implemented in the positioner's .conf file
            self.petal.pos_flags[posmodel.posid] |= self.petal.flags.get('NOTCTLENABLED', self.petal.missing_flag)
            return True
        return False

    def _deny_request_because_both_locked(self, posmodel):
        '''Checks for case where both axes are locked.
        '''
        both_locked = all(posmodel.axis_locks)
        return both_locked

    def _check_init_or_final_neighbor_interference(self, posmodel, final_poslocTP=None):
        """Checks for interference of posmodel with any neighbor positioner or fixed
        boundary, at a single  position (i.e. not a whole sweep).

        When final_poslocTP is None, checks the initial positions.

        When final_poslocTP is a coordinate pair, checks final position of posmodel
        against whatever requested final neighbor positions already exist in the
        schedule.

        Returns either empty set (no interference) or the posids (or boundary names)
        of all interfering neighbors.
        """
        use_final = final_poslocTP != None
        posid = posmodel.posid
        interfering_neighbors = set()
        neighbors = self.collider.pos_neighbors[posid]
        if use_final:
            poslocTP = final_poslocTP
            neighbors_with_requests = {n for n in neighbors if n in self._requests}
        else:
            poslocTP = posmodel.expected_current_poslocTP
            neighbors_with_requests = []
        for n in neighbors:
            n_posmodel = self.petal.posmodels[n]
            if use_final:
                if n in neighbors_with_requests:
                    n_posintTP = self._requests[n]['targt_posintTP']
                else:
                    continue
            else:
                n_posintTP = n_posmodel.expected_current_posintTP
            n_poslocTP = n_posmodel.trans.posintTP_to_poslocTP(n_posintTP)
            if self.collider.spatial_collision_between_positioners(posid, n, poslocTP, n_poslocTP):
                interfering_neighbors.add(n)
        fixed_case = self.collider.spatial_collision_with_fixed(posid, poslocTP)
        if fixed_case:
            interfering_neighbors.add(pc.case.names[fixed_case])
        if interfering_neighbors:
            self.petal.pos_flags[posmodel.posid] |= self.petal.flags.get('OVERLAP', self.petal.missing_flag)
        return interfering_neighbors

    def _deny_request_because_out_of_bounds(self, posmodel, target_poslocTP):
        """Checks for case where a target request is definitively unreachable
        due to being beyond a fixed petal or GFA boundary.

        Returns '' if ok otherwise an error string describing the denial reason.
        """
        out_of_bounds = self.collider.spatial_collision_with_fixed(posmodel.posid, target_poslocTP)
        if out_of_bounds:
            self.petal.pos_flags[posmodel.posid] |= self.petal.flags.get('BOUNDARYVIOLATION', self.petal.missing_flag)
            return f'Target exceeds fixed boundary "{pc.case.names[out_of_bounds]}".'
        return ''

    def _deny_request_because_limit(self, posmodel, poslocTP):
        '''Check for cases where angle exceeds phi limit. This limit may either
        be imposed globally across the petal, or due ot a positioner being
        classified as retracted.

        Returns '' if ok otherwise an error string describing the denial reason.
        '''
        err = ''
        limit_angle = None
        if self.petal.limit_angle:
            limit_angle = self.petal.limit_angle
        if posmodel.classified_as_retracted:
            limit_angle = self.collider.Eo_phi if limit_angle == None else max(self.collider.Eo_phi, limit_angle)
        if limit_angle == None or poslocTP[1] >= limit_angle:
            return err
        # [JHS] 2020-10-20 "EXPERTLIMIT" flag below does not capture the distinction between
        # denying due to global limit vs denying due to classified as retracted. But at least
        # if not flagging the precise cause, it does still flag the same effect.
        self.petal.pos_flags[posmodel.posid] |= self.petal.flags.get('EXPERTLIMIT', self.petal.missing_flag)
        err = f'Target poslocP={poslocTP[1]:.3f} is outside phi limit angle={limit_angle:.1f}.'
        return err

    def _deny_request_because_starting_out_of_range(self, posmodel):
        '''Checks for case where an out-of-range initial POS_T or POS_P (the
        internally-tracked angular postion would cause nonsense moves.
        '''
        posintTP = posmodel.expected_current_posintTP
        rangeT = posmodel.full_range_posintT
        if min(rangeT) > posintTP[pc.T] or max(rangeT) < posintTP[pc.T]:
            return True
        rangeP = posmodel.full_range_posintP
        if min(rangeP) > posintTP[pc.P] or max(rangeP) < posintTP[pc.P]:
            return True
        return False

    def _schedule_moves_initialize_logging(self, anticollision):
        """Initial logging tasks for the schedule_moves() function."""
        if self.stats.is_enabled():
            self.__timer_start = time.perf_counter()
            self.stats.set_scheduling_method(str(anticollision))
            self.__original_request_posids = set(self._requests.keys())
            self.__max_net_time = 0

    def _combine_stages_into_final(self):
        """Takes move tables from each individual stage and combines them into
        the "final" stages.
        """
        final = self.stages['final']
        for name in self.stage_order:
            stage = self.stages[name]
            if stage != final:
                stage.equalize_table_times()
                for posid,table in stage.move_tables.items():
                    if posid not in final.move_tables:
                        final.add_table(table)
                    else:
                        final.move_tables[posid].extend(table)

    def _check_final_stage(self, msg_prefix='', msg_suffix='', assert_no_unresolved=False):
        """Checks the special "final" schedule stage for collisions.

        Inputs: Some prefix and suffix text may be argued for printed messages customization.
                assert_all_resolved ... if True, will print details and throw an assert if detect any unresolved collisions

        Outputs: colliding_sweeps ... dictionary of any colliding sweeps (keys are posids)
                 all_sweeps       ... dictionary of all sweeps checked (keys are posids)
                 collision_pairs  ... set of any collision pair id strings
        """
        assert isinstance(msg_prefix, str)
        assert isinstance(msg_suffix, str)
        final = self.stages['final']
        colliding_sweeps, all_sweeps = final.find_collisions(final.move_tables)
        final.store_collision_finding_results(colliding_sweeps, all_sweeps)
        self.printfunc(msg_prefix + ' collision check --> num colliding sweeps = ' + str(len(colliding_sweeps)) + msg_suffix)
        collision_pairs = {final._collision_id(posid,colliding_sweeps[posid].collision_neighbor) for posid in colliding_sweeps}
        self.printfunc(msg_prefix + ' collision pairs: ' + str(collision_pairs))
        colliding = set(colliding_sweeps)
        if assert_no_unresolved and colliding:
            self.printfunc(self.get_details_str(colliding, label='colliding'))
            err_str = f'{len(colliding)} collisions were NOT resolved! This indicates a bug that needs to be fixed. See details above.'
            self.printfunc(err_str)
            # self.petal.enter_pdb()  # 2023-07-11 [CAD] This line causes a fault (no function enter_pdb)
            assert False, err_str  # 2020-11-16 [JHS] put a PDB entry point in rather than assert, so I can inspect memory next time this happens online
        return colliding_sweeps, all_sweeps, collision_pairs

    def _schedule_moves_check_final_sweeps_continuity(self):
        """Helper function for schedule_moves()."""
        final = self.stages['final']
        if final.sweeps: # indicates that results from a collision check do exist
            if self.should_check_sweeps_continuity:
                discontinuous = final.sweeps_continuity_check()
                self.printfunc('Final check of quantized sweeps --> ' + str(len(discontinuous)) + ' discontinuous (should always be zero)')
                if discontinuous:
                    self.printfunc('Discontinous sweeps: ' + str(sorted(discontinuous.keys())))

    def _schedule_moves_store_collisions_and_pairs(self, colliding_sweeps, collision_pairs):
        """Helper function for schedule_moves()."""
        final = self.stages['final']
        if self.stats.is_enabled():
            self.stats.add_final_collision_check(collision_pairs)
            colliding_posids = set(colliding_sweeps.keys())
            colliding_tables = {p:final.move_tables[p] for p in colliding_posids if p in final.move_tables}
            self.stats.add_unresolved_colliding_at_stage('final', colliding_posids, colliding_tables, colliding_sweeps)

    def _schedule_moves_store_requests_info(self):
        """Goes through the move tables, matching up original request information
        and pushing it into the tables. (This is for logging purposes.)
        """
        stats_enabled = self.stats.is_enabled()
        for posid, table in self.move_tables.items():
            if stats_enabled:
                # for_hardware time is the true time to execute the move, including automatic antibacklash and creep moves (unknown to posschedule)
                self.__max_net_time = max(table.for_hardware()['total_time'], self.__max_net_time)
            log_note_addendum = ''
            if posid in self._requests:
                req = self._requests[posid]
                table.store_orig_command(string=req['command'], val1=req['cmd_val1'], val2=req['cmd_val2']) # keep the original commands with move tables
                log_note_addendum = pc.join_notes(req['log_note'], req['target_str']) # keep req target and the original log notes with move tables
                if stats_enabled:
                    match = self._table_matches_request(table.for_schedule(), req)
                    if posid in self.__original_request_posids and match:
                        self.stats.add_table_matching_request()
            elif not self.expert_mode_is_on():
                self.printfunc('Error: ' + str(posid) + ' has a move table despite no request.')
                table.display()
            table.append_log_note(log_note_addendum)

    def _schedule_moves_finish_logging(self, anim_tables=None):
        """Final logging and animation steps for the schedule_moves() function."""
        anim_tables = {} if anim_tables is None else anim_tables
        if self.stats.is_enabled():
            self.stats.set_num_move_tables(len(self.move_tables))
            self.stats.set_max_table_time(self.__max_net_time)
            resolved_by_freeze = self.stats.get_collisions_resolved_by(method='freeze')
            if resolved_by_freeze:
                self.printfunc(f'{len(resolved_by_freeze)} collision(s) prevented by "freeze" method: {resolved_by_freeze}')

            # Patch for occasional corner corner case where two neighbor positioners both must freeze,
            # but during path adjustment, only one of them got the avoidance event registered.
            frozen = self.get_frozen_posids()
            for posid in frozen:
                avoidances = self.stats.get_avoidances(posid)
                if 'freeze' not in str(avoidances):
                    resolved_this_posid_by_freeze = {collision for collision in resolved_by_freeze if posid in collision}
                    for collision_pair_id in resolved_this_posid_by_freeze:
                        self.stats.add_avoidance(posid, 'freeze', collision_pair_id)

        self.printfunc(f'Num move tables in final schedule = {len(self.move_tables)}')
        if self.verbose:
            self.printfunc(f'posids with move tables in final schedule: {sorted(self.move_tables.keys())}')
        total_time = time.perf_counter() - self.__timer_start
        if self.stats.is_enabled():
            self.stats.add_scheduling_time(total_time)
        self.printfunc(f'Total time to calculate and check schedules = {total_time:.3f} sec')
        if self.petal.animator_on and anim_tables:
            final = self.stages['final']
            dummy_stats = posschedstats.PosSchedStats(enabled=False)
            anim_stage = posschedulestage.PosScheduleStage(collider=final.collider,
                                                           stats=dummy_stats,
                                                           power_supply_map=final._power_supply_map,
                                                           verbose=final.verbose,
                                                           printfunc=final.printfunc)
            for table in anim_tables.values():
                anim_stage.add_table(table)
            colliding_sweeps, all_sweeps = anim_stage.find_collisions(anim_stage.move_tables)
            note = f'move {self.petal.animator_move_number}'
            note_time = self.petal.animator_total_time
            self.collider.animator.set_note(note=note, time=note_time)
            self.petal.animator_move_number += 1
            if self.collider.animate_colliding_only:
                sweeps_to_add = {}
                for posid, sweep in colliding_sweeps.items():
                    sweeps_to_add.update({posid:sweep})
                    neighbor_sweeps = {n:all_sweeps[n] for n in self.collider.pos_neighbors[posid]}
                    sweeps_to_add.update(neighbor_sweeps)
            else:
                sweeps_to_add = {posid: all_sweeps[posid] for posid in all_sweeps if posid in self.collider.posids_to_animate}
            if sweeps_to_add:
                self.collider.add_mobile_to_animator(self.petal.animator_total_time, sweeps_to_add)
                for posid in sweeps_to_add:
                    self.collider.add_posid_label(posid)
                self.petal.animator_total_time += max({sweep.time[-1] for sweep in sweeps_to_add.values()})
                if self.collider.animate_colliding_only:
                    self.printfunc('Added ' + str(len(colliding_sweeps)) + ' colliding sweeps (and their neighbors) to the animator.')

    def _table_matches_quantized_sweep(self, move_table, sweep):
        """Takes as input a "for_schedule()" move table and a quantized sweep,
        and then cross-checks whether their total rotations (theta and phi)
        match. Returns a boolean.
        """
        tol = pc.schedule_checking_numeric_angular_tol
        table = move_table.for_schedule()
        end_tp_sweep = [sweep.theta(-1), sweep.phi(-1)]
        end_tp_table = [table['net_dT'][-1] + sweep.theta(0), table['net_dP'][-1] + sweep.phi(0)]
        endpos_diff = [end_tp_sweep[i] - end_tp_table[i] for i in range(2)]
        if abs(endpos_diff[0]) > tol or abs(endpos_diff[1]) > tol:
            self.printfunc(f'table and sweep not matched: {sweep.posid} end_tp: check={end_tp_sweep}, move={end_tp_table}')
            return False
        return True

    def _table_matches_request(self, table_for_schedule, request):
        """Input a move table (for_schedule format) and check whether the total
        motion matches request. Returns a boolean.
        """
        tol = pc.schedule_checking_numeric_angular_tol
        dtdp_request = [request['targt_posintTP'][i] - request['start_posintTP'][i] for i in range(2)]
        dtdp_table = [table_for_schedule['net_dT'][-1], table_for_schedule['net_dP'][-1]]
        diff = [dtdp_request[i] - dtdp_table[i] for i in range(2)]
        diff_abs = [abs(x) for x in diff]
        for i in range(len(diff_abs)):
            if diff_abs[i] > 180:
                diff_abs[i] -= 360
        if diff_abs[0] > tol or diff_abs[1] > tol:
            return False
        return True

    def _denied_str(self, target_str, msg_str):
        return pc.join_notes('Target request denied', target_str, msg_str)

    def _make_coord_str(self, uv_type, uv, prefix=''):
        '''Make a string showing a coordinate pair in a standard format.
        '''
        s = f'{prefix}_' if prefix else ''
        s += f'{uv_type}='
        uv_strs = [f'{x:.3f}' for x in uv]
        uv_str = ", ".join(uv_strs)
        s += f'({uv_str})'
        return s

POS_DISABLED_MSG = 'Positioner is disabled.'
BOTH_AXES_LOCKED_MSG = 'Both theta and phi axes are locked.'
