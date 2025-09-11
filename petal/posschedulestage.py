import posconstants as pc
import posmovetable
import math

class PosScheduleStage(object):
    """This class encapsulates the concept of a 'stage' of the fiber
    positioner motion. The typical usage would be either a direct stage from
    start to finish, or an intermediate stage, used for retraction, rotation,
    or extension.

        collider         ... instance of poscollider for this petal
        stats            ... instance of posschedstats for this petal
        power_supply_map ... dict where key = power supply id, value = set of posids attached to that supply
    """
    def __init__(self, collider, stats, power_supply_map=None, verbose=False, printfunc=None, petal=None):
        self.collider = collider # poscollider instance
        self.move_tables = {} # keys: posids, values: posmovetable instances
        self.start_posintTP = {} # keys: posids, values: initial positions at start of stage
        self.sweeps = {} # keys: posids, values: instances of PosSweep, corresponding to entries in self.move_tables
        self.colliding = set() # positioners currently known to have collisions
        self.stats = stats
        self._power_supply_map = {} if power_supply_map is None else power_supply_map
        self._theta_max_jog_A = 50 # deg, maximum distance to temporarily shift theta when doing path adjustments
        self._theta_max_jog_B = 20
        self._phi_max_jog_A = 45 # deg, maximum distance to temporarily shift phi when doing path adjustments
        self._phi_max_jog_B = 90
        self._max_jog = self._assign_max_jog_values() # collection of all the max jog options above
        self.sweep_continuity_check_stepsize = 4.0 # deg, see PosSweep.check_continuity function
        self.verbose = verbose
        self.printfunc = printfunc
        self.petal_debug = petal.petal_debug if hasattr(petal, 'petal_debug') else {}

    def initialize_move_tables(self, start_posintTP, dtdp, update_only=False):
        """Generates basic move tables for each positioner, starting at position
        start_tp and going to a position a distance dtdp away.

            start_posintTP  ... dict of starting [theta,phi] positions, keys are posids
            dtdp            ... dict of [delta theta, delta phi] from the starting position. keys are posids
            update_only     ... boolean, if True then only update existing start_posintTP and move_tables
                                (The default is False, in which case all move_tables and associated data ---
                                 start_posintTP, sweeps --- are wiped out and replaced with those argued here.

        The user should take care that dtdp vectors have been generated using the
        canonical delta_posintTP function provided by PosTransforms module, with
        range_wrap_limits='targetable'. (The user should *not* generate dtdp by a
        simple vector subtraction or the like, since this would not correctly handle
        physical range limits of the fiber positioner.)
        """
        posids_to_delete = set(start_posintTP) if update_only else set(self.move_tables)
        for posid in posids_to_delete:
            self.del_table(posid)
        for posid in start_posintTP:
            posmodel = self.collider.posmodels[posid]
            final_posintTP = posmodel.trans.addto_posintTP(start_posintTP[posid], dtdp[posid], range_wrap_limits='targetable')
            true_dtdp = posmodel.trans.delta_posintTP(final_posintTP, start_posintTP[posid], range_wrap_limits='targetable')
            table = posmovetable.PosMoveTable(posmodel, start_posintTP[posid])
            table.set_move(0, pc.T, true_dtdp[0])
            table.set_move(0, pc.P, true_dtdp[1])
            table.set_prepause(0, 0.0)
            table.set_postpause(0, 0.0)
            self.move_tables[posid] = table
            self.start_posintTP[posid] = tuple(start_posintTP[posid])

    def _print_table_diff(self, posid, old_table, new_table):
        """
        Print only the differences between two move tables
        """
        o_table = {}
        n_table = {}
        for k, v in old_table.items():
            if k in new_table:
                if v != new_table[k]:
                    o_table[k] = v
                    n_table[k] = new_table[k]
        if o_table:
            self.printfunc(f'zeno movetable change\nold {posid} table: {str(o_table)}\nnew {posid} table: {str(n_table)}')
        return

    def rewrite_zeno_move_tables(self, proposed_tables):
        for posid, table in proposed_tables.items():
            if table.posmodel.is_linphi:
                # self.printfunc(f'Rewriting zeno table for {posid}')
                new_table = self.rewrite_zeno_move_table(table)
                if new_table is not None:
                    vrbose = self.petal_debug.get('linphi_verbose')
                    try:
                        if vrbose and int(vrbose) > 1:
                            self._print_table_diff(posid, table.as_dict(), new_table.as_dict())
                    except TypeError:
                        pass
                    proposed_tables[posid] = new_table
        return proposed_tables

    def rewrite_zeno_move_table(self, table):
        linphi_table = None
        if table.posmodel.is_linphi:
            linphi_table = self._rewrite_linphi_move_table(table)
        return linphi_table

    def _rewrite_linphi_move_table(self, table, verbose=False):
#       last_motor_direction is always > 0
#       last_motor_direction = table.posmodel.linphi_params['LAST_P_DIR']
        last_motor_direction = 1
        vrbose = self.petal_debug.get('linphi_verbose')
        if self.petal_debug.get('compact_linphi'):
            table.compact()
        if table.has_phi_motion:
            new_table = table.copy()
            idx = 0
            l_idx = 0
            if verbose:
                self.printfunc('Proposed table has phi movement') # DEBUG
            for row in table.rows:
                phi_dist = table.get_move(idx, pc.P)
                theta_dist = table.get_move(idx, pc.T)
                postpause = table.get_postpause(idx)
                if phi_dist == 0:
                    if verbose:
                        self.printfunc(f'no phi movement in old row {idx}, new row {l_idx}, skipping') # DEBUG
                    idx += 1
                    l_idx += 1
                else:
                    new_direction = 1 if phi_dist >= 0.0 else -1
#                   scale_ccw = float(table.posmodel.linphi_params['CCW_SCALE_A'])
#                   scale_cw = float(table.posmodel.linphi_params['CW_SCALE_A'])
#                   We don't need the scales here since Petalcontroller will apply them
                    #NOTE: The first and second moves should have abs(move) >= pc.P_zeno_jog
                    if last_motor_direction == 1:  # table.posmodel.linphi_params['LAST_P_DIR']:
                        if new_direction > 0:   # must go negative, then positive
                            first_move = -pc.P_zeno_jog # / scale_cw
                            second_move = (pc.P_zeno_jog + phi_dist) # / scale_ccw
                        else:                       # must go positive, then negative
                            first_move = (-pc.P_zeno_jog + phi_dist) # / scale_cw
                            second_move = pc.P_zeno_jog # / scale_ccw
                    else:
                        if new_direction > 0:   # must go positive, then negative
                            first_move = (pc.P_zeno_jog + phi_dist) # / scale_ccw
                            second_move = -pc.P_zeno_jog # / scale_cw
                        else:                       # must go positive, then negative
                            first_move = (-pc.P_zeno_jog + phi_dist) # / scale_cw
                            second_move = pc.P_zeno_jog # / scale_ccw
#                   Would need next two lines to prevent banging into hard stops, but should never get here
#                   if that were going to happen since additional keepout is > jog size
#                   first_move_limited = self._range_limited_jog(first_move ... and other args)
#                   second_move_limited = self._range_limited_jog(second_move ... and other args)
                    if verbose:
                        self.printfunc(f'original index = {idx}, new indices = {l_idx}, {l_idx+1}') # DEBUG
                    new_table.set_move(l_idx, pc.P, first_move)
                    new_table.set_move(l_idx, pc.T, 0.0)
                    new_table.set_postpause(l_idx, 0)
                    new_table.insert_new_row(l_idx + 1)
                    new_table.set_move(l_idx + 1, pc.P, second_move)
                    new_table.set_move(l_idx + 1, pc.T, theta_dist)
                    if postpause:
                        new_table.set_postpause(l_idx + 1, postpause)
# Second move is always >0, so LAST_P_DIR is always 1, never -1
#                   table.posmodel.linphi_params['LAST_P_DIR'] = 1 if second_move > 0 else -1  # store new direction
                    idx += 1
                    l_idx += 2
            else:
                if verbose:
                    self.printfunc('Proposed table has no phi movement') # DEBUG
            if idx != l_idx:    # table was modified
                return new_table
        return None

    def is_not_empty(self):
        """Returns boolean whether the stage is empty of move_tables.
        """
        return not(not(self.move_tables))

    def add_table(self, move_table):
        """Directly adds a move table to the stage. If a move table representing
        same positioner already exists, then that table is extended.
        """
        this_posid = move_table.posid
        if move_table.posid in self.move_tables:
            self.move_tables[this_posid].extend(move_table)
        else:
            self.move_tables[this_posid] = move_table

    def del_table(self, posid):
        '''Deletes a move table and associated sweep data. This may leave the
        state of self.colliding out of date, until the next find_collisions() call.
        '''
        for d in [self.move_tables, self.start_posintTP, self.sweeps]:
            if posid in d:
                del d[posid]

    def anneal_tables(self, suppress_automoves=False, mode='filled'):
        """Adjusts move table timing, to attempt to reduce peak power consumption
        of the overall array.

            anneal_time ... Time in seconds over which to spread out moves in this stage
                            to reduce overall power density consumed by the array.

            suppress_automoves ... Boolean, if True, don't include auto-generated final
                                   creep and antibacklash moves in the move times. Typically
                                   you want this True only for stages which are not the last
                                   one in the sequence.

            mode ... 'filled' --> try to most efficiently fill time with moves
                     'ramped' --> try to ramp up/down the power (takes more time)

        If anneal_time is less than the time it takes to execute the longest move
        table, then that longer execution time will be used instead of anneal_time.
        """
        assert mode in pc.anneal_density, f'unrecognized anneal mode {mode}'
        times = {posid: table.total_time(suppress_automoves=suppress_automoves)
                 for posid, table in self.move_tables.items()
                 if not(table.is_motionless)}
        if not times:
            return
        sorted_times = sorted(times.values())[::-1]
        sorted_posids = sorted(times, key=lambda k: times[k])[::-1]
        orig_max_time = max(sorted_times)
        anneal_window = sum(sorted_times) / len(sorted_times) / pc.anneal_density[mode]
        anneal_window = max(anneal_window, orig_max_time)  # for case of very large outlier

        if mode == 'filled':
            def first_within(x, vec):
                for i, test in enumerate(vec):
                    if test <= x:
                        return i
                return None
            for map_posids in self._power_supply_map.values():
                posids = [p for p in sorted_posids if p in map_posids]  # maintains sorted-by-time order
                times2 = [sorted_times[i] for i in range(len(sorted_posids)) if sorted_posids[i] in map_posids]
                prepause = 0.0
                while posids:
                    open_time = anneal_window - prepause
                    i = first_within(open_time, times2)
                    if i == None:
                        prepause = 0.0
                    else:
                        posid = posids[i]
                        self.move_tables[posid].insert_new_row(0)
                        self.move_tables[posid].set_prepause(0, prepause)
                        prepause += times2[i]
                        del posids[i]
                        del times2[i]

        elif mode == 'ramped':
            resolution = 0.1 # sec
            for map_posids in self._power_supply_map.values():
                posids = [p for p in sorted_posids if p in map_posids]  # maintains sorted-by-time order
                prepause = 0.0
                which = 0
                while posids:
                    posid = posids.pop(which)
                    prepause = 0.0 if prepause + times[posid] > anneal_window else prepause + resolution
                    self.move_tables[posid].insert_new_row(0)
                    self.move_tables[posid].set_prepause(0, prepause)
                    which = -1 if which == 0 else 0

    def equalize_table_times(self):
        """Makes all move tables in the stage have an equal total time length,
        by adding in post-pauses wherever necessary.
        """
        if not self.move_tables:
            return
        times = {}
        for posid,table in self.move_tables.items():
            postprocessed = table.for_schedule()
            times[posid] = postprocessed['net_time'][-1]
        max_time = max(times.values())
        for posid,table in self.move_tables.items():
            equalizing_pause = max_time - times[posid]
            if equalizing_pause:
                idx = table.n_rows
                table.insert_new_row(idx)
                table.set_postpause(idx,equalizing_pause)
                if self.sweeps: # because no collision checking is performed if anticollsion=None
                    if posid in self.sweeps.keys():
                        self.sweeps[posid].extend(self.collider.timestep, max_time)
        if self.verbose:
            self.printfunc(f'posschedulestage: move time after annealing = {max_time}')
        return max_time

    def adjust_path(self, posid, freezing='on', do_not_move=None):
        """Adjusts move paths for posid to avoid collision. If the positioner
        has no collision, then no adjustment is made.

            posid    ... positioner to adjust. This method assumes that this posid
                         already has a move_table associated with it in the stage,
                         that find_collsions() has already been run, and the results
                         saved via store_collision_finding_results().

            freezing ... string with the following settings for when the "freeze"
                         method of collision resolution shall be applied:

                'on'     ... freeze positioner if the path adjustment options all fail to resolve collisions
                'off'    ... don't freeze, even if the path adjustment options all fail to resolve collisions
                'forced' ... only freeze, and *must* do so
                'forced_recursive' ... like 'forced', then closes any follow-on neighbor collisions

            do_not_move ... set of posids which are *not* allowed to be automatically
                            moved (i.e. like moved out of the way) of posid

        Returns a tuple:

            item 0 ... set containing the posids of any robot(s) whose move
                       tables were adjusted in the course of the function call

            item 1 ... set of only those posids which were specifically "frozen"

        With freezing == 'off' or 'on', the path adjustment algorithm goes through a
        series of options, trying adding various pauses and pre-moves to avoid collision.

        With freezing == 'on', then if all the adjustment options fail, the fallback
        is to attempt to "freeze" the positioner posid. This means its path is simply halted
        prior to the collision, and no attempt is made for it to reach its supposed final
        target. This may still fail, for example in the case where two robots are simultaneously
        sweeping through each other (rather than one simply striking another). In this case, one
        must use freezing == 'forced' if you want to freeze this posid anyway.

        With freezing == 'forced', we skip calculating the various path adjustment options,
        and instead go straight to freezing the positioner before it collides.

        With freezing == 'forced_recursive', it is just like 'forced', plus at the end
        we look for any follow-on collisions among neighbors (due to a neighbor continuing to
        move, or a new side-effect collision). These will be recursively searched out and
        resolved by more forced freezing. This is intended as the final adjustment method,
        to definitively prevent any collisions including side-effects.

        The timing of a neighbor's motion path may be adjusted as well by this
        function, but not the geometric path that the neighbor follows.

        Adjustments to neighbor's timing only go one level deep of neighborliness.
        In other words, a neighbor's neighbor will not be affected by this function.
        """
        stats_enabled = self.stats.is_enabled()
        adjusted = set()
        frozen = set()
        if self.sweeps[posid].collision_case == pc.case.I:
            return adjusted, frozen
        elif self.sweeps[posid].collision_case in pc.case.fixed_cases:
            methods = ['freeze'] if freezing != 'off' else []
        elif freezing in {'forced','forced_recursive'}:
            methods = ['freeze']
        elif freezing == 'off':
            methods = pc.nonfreeze_adjustment_methods
        else:
            methods = pc.all_adjustment_methods
        for method in methods:
            collision_neighbor = self.sweeps[posid].collision_neighbor
            proposed_tables = self._propose_path_adjustment(posid, method, do_not_move)
#           proposed_tables = self.rewrite_zeno_move_tables(proposed_tables)
            colliding_sweeps, all_sweeps = self.find_collisions(proposed_tables)
            should_accept = not(colliding_sweeps) or freezing in {'forced','forced_recursive'}
            should_accept &= any(proposed_tables) # nothing to accept if no proposed tables were generated
            if should_accept:
                self.move_tables.update(proposed_tables)
                adjusted.update(proposed_tables.keys())

                # search for any side effect new collisions
                proposed = set(proposed_tables.keys())
                could_have_changed = {posid, collision_neighbor}
                for p in proposed:
                    could_have_changed.update(self.collider.pos_neighbors[p])
                could_have_changed -= proposed # saves a little time -- no need to recheck these
                tables_to_recheck = {p:self.move_tables[p] for p in could_have_changed if p in self.move_tables}
                col_recheck, all_recheck = self.find_collisions(tables_to_recheck)
                col_to_store = {p:s for p,s in col_recheck.items() if p in tables_to_recheck.keys()}
                all_to_store = {p:s for p,s in all_recheck.items() if p in tables_to_recheck.keys()}
                col_to_store.update({p:s for p,s in colliding_sweeps.items() if p in proposed})
                all_to_store.update({p:s for p,s in all_sweeps.items() if p in proposed})

                # determine if collision was resolved (for statistics tracking)
                if stats_enabled:
                    old_collision_id = self._collision_id(posid,collision_neighbor)
                    if posid in col_to_store:
                        new_collision_id = self._collision_id(posid,col_to_store[posid].collision_neighbor)
                        collision_resolved = new_collision_id != old_collision_id
                    else:
                        collision_resolved = True
                    if collision_resolved:
                        self.stats.add_collisions_resolved(posid, method, {old_collision_id})

                # store results
                old_colliding = self.colliding # note how sequence here emphasizes that this must occur before store_collision_finding_results(), which affects self.colliding. in a perfect world, I would re-factor functionally to remove the state-dependence [JHS]
                self.store_collision_finding_results(col_to_store, all_to_store)
                if method == 'freeze':
                    self.sweeps[posid].register_as_frozen() # needs to occur after storing results above
                    frozen.add(posid)

                # recursively-forced freezing
                if freezing == 'forced_recursive':
                    newly_colliding = set(col_to_store.keys()).difference(old_colliding)
                    remainder = newly_colliding
                    if collision_neighbor in self.colliding:
                        remainder.add(collision_neighbor)
                    if newly_colliding:
                        self.printfunc('Note: adjust_path(' + str(posid) + ', freezing=\'' + str(freezing) + '\') introduced new collisions for ' + str(newly_colliding))
                    for p in sorted(remainder): # sort is for repeatabiity (since 'remainder' is an unordered set, and so path adjustments would otherwise get processed in variable order from run to run)
                        freeze_is_possible = False # starting assumption for this pos
                        if p in self.move_tables: # does p have any move_table to be frozen?
                            if not self.move_tables[p].is_motionless: # does that table have any contents inside to be frozen?
                                freeze_is_possible = True
                        if freeze_is_possible:
                            newly_adjusted, recursed_newly_frozen = self.adjust_path(p, freezing='forced_recursive') # recursively close out any side-effect new collisions
                            adjusted.update(newly_adjusted)
                            frozen.update(recursed_newly_frozen)
                        else:
                            self.printfunc(' --> no further freezing possible on ' + str(p) + ' --- already motionless')
                        verified = p not in self.colliding
                        self.printfunc(' --> recursive forced freeze attempted on ' + str(p) + '. Verified now non-colliding? ' + str(verified))
                break # note indentation level of this return statement is essential. it breaks out of the methods for loop. do not remove again!
        return adjusted, frozen

    def find_collisions(self, move_tables, skip=0):
        """Identifies any collisions that would be induced by executing a collection
        of move tables.

            move_tables ... dict with keys = posids, values = PosMoveTable instances.
            skip ... integer number of initial timesteps for which to skip collision checks

        Two items are returned. They are both dicts with keys = posids, values = PosSweep
        instances (see poscollider.py).

          1st dict: Only contains sweeps for positioners that collide. It will be
                    empty if there are no collisions. These sweeps may indicate
                    collision with another positioner or with a fixed boundary. This
                    information is given internally within each sweep instance.

          2nd dict: Contains all sweeps that were generated during this function call.

        For any pair of positioners that collide, the returned collisions dict will contain
        separate sweeps for each of the pair. These two sweeps are giving you information
        about the same collision event, but from the perspectives of the two different
        positioners. In other words, if there is an entry for posid 'M00001', colliding with
        neighbor 'M00002', then the dict will also contain an entry for posid 'M00002',
        colliding with neighbor 'M0001'.

        If positioner A collides with a fixed boundary, or with a disabled neighbor
        positioner B, then the returned collisions dict only contains the sweep of A.

        If a positioner has collisions with multiple other postioners / fixed boundaries,
        then only the first collision event in time is included in the returned collisions
        dict.

        In the rare event of a three-way exactly simultaneous collision between three
        moving positioners, then all three of those positioners' sweeps would still
        appear in the return dictionary.
        """
        already_checked = {posid:set() for posid in self.collider.posids}
        colliding_sweeps = {posid:set() for posid in self.collider.posids}
        all_sweeps = {}
        for posid in move_tables:
            table_A = move_tables[posid]
            for neighbor in self.collider.pos_neighbors[posid]:
                if neighbor not in already_checked[posid]:
                    table_B = move_tables[neighbor] if neighbor in move_tables else self._get_or_generate_table(neighbor)
                    pospos_sweeps = self.collider.spacetime_collision_between_positioners(
                                            posid, table_A.init_poslocTP, table_A.for_collider(),
                                            neighbor, table_B.init_poslocTP, table_B.for_collider(),
                                            skip=skip)
                    all_sweeps.update({posid:pospos_sweeps[0], neighbor:pospos_sweeps[1]})
                    for sweep in pospos_sweeps:
                        if sweep.collision_case != pc.case.I:
                            colliding_sweeps[sweep.posid].add(sweep)
                    already_checked[posid].add(neighbor)
                    already_checked[neighbor].add(posid)
            for fixed_neighbor in self.collider.fixed_neighbor_cases[posid]:
                posfix_sweep = self.collider.spacetime_collision_with_fixed(
                                       posid, table_A.init_poslocTP, table_A.for_collider(),
                                       skip=skip)[0] # index 0 to immediately retrieve from the one-element list this function returns
                all_sweeps.update({posid:posfix_sweep}) # don't worry --- if pospos colliding sweep takes precedence, this will be appropriately replaced again below
                if posfix_sweep.collision_case != pc.case.I:
                    colliding_sweeps[posid].add(posfix_sweep)
        multiple_collisions = {posid for posid in colliding_sweeps if len(colliding_sweeps[posid]) > 1}
        for posid in multiple_collisions:
            first_collision_time = math.inf
            for sweep in colliding_sweeps[posid]:
                if sweep.collision_time < first_collision_time:
                    first_sweep = sweep
                    first_collision_time = sweep.collision_time
            colliding_sweeps[posid] = {first_sweep}
        colliding_sweeps = {posid:colliding_sweeps[posid].pop() for posid in colliding_sweeps if colliding_sweeps[posid]} # remove set structure from elements, and remove empty elements
        all_sweeps.update(colliding_sweeps)
        return colliding_sweeps, all_sweeps

    def store_collision_finding_results(self, colliding_sweeps, all_sweeps):
        """Stores sweep data as-generated by find_collisions method.
        """
        previously_frozen = {posid for posid, sweep in self.sweeps.items() if sweep.is_frozen}
        should_update = set(all_sweeps.keys()) - previously_frozen
        updates = {posid: all_sweeps[posid] for posid in should_update}
        self.sweeps.update(updates)
        all_checked = {posid for posid in all_sweeps}
        now_colliding = {posid for posid in colliding_sweeps}
        now_not_colliding = all_checked.difference(now_colliding)
        self.colliding = self.colliding.union(now_colliding)
        self.colliding = self.colliding.difference(now_not_colliding)
        # check for special case of lingering incorrect colliding status (no move table therefore wasn't registered as resolved)
        no_table_but_in_colliding = {posid for posid in self.colliding if posid not in self.move_tables}
        for posid in no_table_but_in_colliding:
            n = self.sweeps[posid].collision_neighbor
            if n and self.sweeps[n].collision_neighbor != posid: # i.e., neighbor thinks this collision has been resolved
                self.sweeps[posid].clear_collision()
                self.colliding.remove(posid)
                if posid in colliding_sweeps:
                    colliding_sweeps.pop(posid)
        if self.stats.is_enabled():
            found = {self._collision_id(posid, sweep.collision_neighbor) for posid, sweep in colliding_sweeps.items()}
            self.stats.add_collisions_found(found)

    def sweeps_continuity_check(self):
        """Returns set of posids for any whose sweeps were found to be discontinous.
        """
        discontinuous = {}
        for p,s in self.sweeps.items():
            if not s.check_continuity(self.sweep_continuity_check_stepsize, self.collider.posmodels[p]):
                discontinuous.add(p)
        return discontinuous

    def _propose_path_adjustment(self, posid, method='freeze', do_not_move=None):
        """Generates a proposed alternate move table for the positioner posid
        The alternate table is meant to attempt to avoid collision.

          posid  ... The positioner to propose a path adjustment for.

          method ... The type of adjustment to make. Valid selections are:

             'pause'     ... Pre-delay is added to the positioner's move table in this
                             stage, to wait for the neighbor to possibly move out of the way.

             'extend_X'    ... Positioner phi arm is first extended out, in an attempt to open
                               a clear path for neighbor.

             'retract_X'   ... Positioner phi arm is first retracted in, in an attempt to open
                               a clear path for neighbor.

             'rot_ccw_X'   ... Positioner theta axis is first rotated ccw, in an attempt to
                               open a clear path for neighbor.

             'rot_cw_X'    ... Positioner theta axis is first rotated cw, in an attempt to
                               open a clear path for neighbor.

             'repel_ccw_X' ... Positioner theta axis is first rotated ccw, while simultaneously
                               the neighbor axis is rotated the opposite direction (cw).

             'repel_cw_X'  ... Positioner theta axis is first rotated cw, while simultaneously
                               the neighbor axis is rotated the opposite direction (ccw).

             'freeze'      ... Positioner is halted prior to the collision, and no attempt
                               is made for its final target.

          do_not_move ... see comments in adjust_path() docstr

        The subscript 'X' in many of the adjustment methods above refers to the
        size of the maximum distance jog step to attempt, labeled 'A', 'B', etc.

        The return value is a dict of proposed new move tables. Key = posid, value =
        new move table. If no change is proposed, the dict is empty. If the proposal
        includes changing both the argued positioner and its neighbor, than they will
        both have tables in the dict.

        No positioner or neighbor will be changed if it is disabled, or if it has
        already been "frozen".

        No new collision checking is performed by this method. So it is important
        after receiving a new proposal, to re-check for collisions against all the
        neighbors of all the proposed new move tables.
        """
        # [JHS] 2020-10-29 Additional methods to consider implementing are *combined*
        # rotation followed immediately by retraction or extension. It's a fancier
        # move, which I believe may solve a few corner cases. Might significantly affect
        # calculation time? Because once you get into combinations, you have a whole lot
        # more possible path adjustment cases. I.e. 4 rot x 4 extension/retraction options
        # --> 16 additional methods! If go this route, would be perhaps cleanest to
        # generate them in code as combo of these options. (At initialization time ---
        # one time.) Anyway, worth considering and doing a few timing tests. Really, most
        # time is spent doing collision checks --- that's where the real cost of additional
        # path adjustments has to be considered.
        if self.sweeps[posid].collision_case == pc.case.I: # no collision
            return {}
        if self.sweeps[posid].is_frozen:
            return {}  # can't do anything with an already-frozen positioner
        if not self.collider.posmodels[posid].is_enabled:
            return {}  # can't do anything with a disabled positioner
        fixed_collision = self.sweeps[posid].collision_case in pc.case.fixed_cases
        if fixed_collision and method in pc.useless_with_fixed_boundary:
            return {}
        neighbor = self.sweeps[posid].collision_neighbor
        do_not_move = set() if do_not_move == None else do_not_move
        unmoving_neighbor = neighbor not in self.move_tables or neighbor in do_not_move
        if unmoving_neighbor and method in pc.useless_with_unmoving_neighbor:
            return {}

        # table that will be adjusted
        tables = {posid: self._get_or_generate_table(posid,should_copy=True)}
        tables_data = {}
        if method == 'freeze' or 'repel' in method:
            tables_data[posid] = tables[posid].for_schedule()

        # freeze method
        if method == 'freeze':
            for row_idx in reversed(range(tables[posid].n_rows)):
                net_time = tables_data[posid]['net_time'][row_idx]
                collision_time = self.sweeps[posid].collision_time
                if math.fmod(net_time, self.collider.timestep): # i.e., if not exactly matching a quantized timestep (almost always the case)
                    collision_time -= self.collider.timestep # handles coarseness of discrete time by treating the collision time as one timestep earlier in the sweep
                if net_time >= collision_time:
                    if self.verbose:
                        self.printfunc("net time, collision_time " + str(net_time) + ' ' + str(collision_time))
                    tables[posid].delete_row(row_idx)
                else:
                    break
            if tables[posid].n_rows == 0:
                tables[posid].set_move(0, 0, 0)
            return tables

        # get neighbor table
        posmodels = {posid:self.collider.posmodels[posid], neighbor:self.collider.posmodels[neighbor]}
        neighbor_can_move = posmodels[neighbor].is_enabled and not self.sweeps[neighbor].is_frozen and neighbor not in do_not_move
        if neighbor_can_move:
            tables[neighbor] = self._get_or_generate_table(neighbor,should_copy=True)
            tables_data[neighbor] = tables[neighbor].for_schedule()
            for neighbor_clearance_time in tables_data[neighbor]['net_time']:
                if neighbor_clearance_time > self.sweeps[posid].collision_time:
                    break
            neighbor_clearance_time += pc.num_timesteps_clearance_margin * self.collider.timestep
        elif method == 'pause':
            return {} # no point in pausing if neighbor never moves
        else:
            neighbor_clearance_time = 0

        # pause method
        if method == 'pause':
            tables[posid].insert_new_row(0)
            tables[posid].set_prepause(0,neighbor_clearance_time)
            return {posid:tables[posid]} # exclude neighbor here, since nothing being done to it

        # retract, rot, extend, repel methods: calculate jog distances
        max_abs_jog = abs(self._max_jog[method])
        jogs = {} # will hold distance(s) to move away from target and then back toward target
        jog_times = {} # will hold corresponding move time for each jog
        if 'extend' in method or 'retract' in method:
            axis = pc.P
            limits = posmodels[posid].targetable_range_posintP
            if 'retract' in method:
                limits[1] = min(limits[1], self.collider.Ei_phi) # note deeper retraction than Eo, to give better chance of avoidance
                direction = +1
            else:
                direction = -1
            jogs[posid] = self._range_limited_jog(nominal=max_abs_jog, direction=direction, start=tables[posid].init_posintTP[1], limits=limits)
        elif 'rot' in method:
            axis = pc.T
            direction = +1 if 'ccw' in method else -1
            jogs[posid] = self._range_limited_jog(nominal=max_abs_jog, direction=direction, start=tables[posid].init_posintTP[0], limits=posmodels[posid].targetable_range_posintT)
        elif 'repel' in method:
            axis = pc.T
            for p,t in tables.items():
                start = t.init_posintTP[0]
                limits = posmodels[p].targetable_range_posintT
                direction = 1 if p == posid else -1 # neighbor repels from primary
                if 'ccw' not in method:
                    direction *= -1
                if direction > 0:
                    furthest_existing_excursion = max(tables_data[p]['net_dT'] + [0]) # zero element prevents accidentally adding margin (if all net_dT < 0) that temporally might not exist when you hoped it would
                    limits[1] -= furthest_existing_excursion
                else:
                    furthest_existing_excursion = min(tables_data[p]['net_dT'] + [0]) # zero element prevents accidentally adding margin (if all net_dT > 0) that temporally might not exist when you hoped it would
                    limits[0] -= furthest_existing_excursion # minus sign here is correct, since furthest_existing_excursion <= 0
                jogs[p] = self._range_limited_jog(nominal=max_abs_jog, direction=direction, start=start, limits=limits)
        else:
            self.printfunc('Error: unknown path adjustment method \'' + str(method) + '\'')
            return {}
        for p in jogs:
            jog_times[p] = posmodels[p].true_move(axisid=axis, distance=jogs[p], allow_cruise=True, limits=None, init_posintTP=None)['move_time']

        # apply jogs and associated pauses to tables
        if 'repel' in method:
            wait_for_posid_to_jog = 0
            wait_for_neighbor_to_jog = 0
            if neighbor in tables:
                jog_times_diff = jog_times[neighbor] - jog_times[posid]
                wait_for_neighbor_to_jog = max(jog_times_diff, 0)
                wait_for_posid_to_jog = max(-jog_times_diff, 0)
            for p,t in tables.items():
                t.insert_new_row(0)
                t.set_move(0, axis, jogs[p])
                if p == posid:
                    t.set_postpause(0, wait_for_neighbor_to_jog + neighbor_clearance_time)
                else:
                    t.set_postpause(0, wait_for_posid_to_jog)
                t.set_move(t.n_rows, axis, -jogs[p]) # note how this is happening in a new final row
        elif neighbor_can_move:
            tables[posid].insert_new_row(0)
            tables[posid].insert_new_row(0)
            tables[posid].set_move(0, axis, jogs[posid])
            tables[posid].set_postpause(0, neighbor_clearance_time)
            tables[posid].set_move(1, axis, -jogs[posid])
            if neighbor in tables:
                tables[neighbor].insert_new_row(0)
                tables[neighbor].set_prepause(0, jog_times[posid])
        else:
            tables[posid].insert_new_row(0)
            tables[posid].set_move(0, axis, jogs[posid])
            new_final_idx = tables[posid].n_rows
            tables[posid].set_move(new_final_idx, axis, -jogs[posid])
        return tables

    @staticmethod
    def _range_limited_jog(nominal, direction, start, limits):
        """Returns a range-limited jog distance.
            nominal   ... nominal jog distance, prior to applying range limits
            direction ... +1 means counter-clockwise, -1 means clockwise
            start     ... starting angle position
            limits    ... 2-element list or tuple stating max and min accessible angles
        """
        if direction >= 0:
            limited = min(start + nominal, max(limits))
        else:
            limited = max(start - nominal, min(limits))
        return limited - start

    def _get_or_generate_table(self, posid, should_copy=False):
        """Fetches move table for posid from self.move_tables. If no such table
        exists, generates a new one. The should_copy flag allows you to request
        that the returned table be a duplicate of the original table (not a pointer
        to it.)
        """
        if posid in self.move_tables:
            table = self.move_tables[posid]
            if should_copy:
                table = table.copy()
        else:
            start = self.start_posintTP[posid] if posid in self.start_posintTP else None
            table = posmovetable.PosMoveTable(self.collider.posmodels[posid], start)
            table.insert_new_row(0)
        return table

    def _assign_max_jog_values(self):
        """Makes convenience dict containing theta and phi max jog values for
        all adjustment methods.
        """
        max_jogs = {}
        for method in pc.nonfreeze_adjustment_methods:
            if '_A' in method:
                if 'extend' in method or 'retract' in method:
                    max_jogs[method] = self._phi_max_jog_A
                else:
                    max_jogs[method] = self._theta_max_jog_A
            else:
                if 'extend' in method or 'retract' in method:
                    max_jogs[method] = self._phi_max_jog_B
                else:
                    max_jogs[method] = self._theta_max_jog_B
        return max_jogs

    @staticmethod
    def _collision_id(A,B):
        """Returns an id string combining string A and string B. The returned
        string will be the same regardless of whether A or B is argued first.
        """
        s = sorted({str(A),str(B)})
        return s[0] + '-' + s[1]
