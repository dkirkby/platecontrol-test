import posmodel
import posconstants as pc
import copy as copymodule

class PosMoveTable(object):
    """A move table contains the information for a single positioner's move
    sequence, in both axes. This object defines the move table structure, and
    presents functions to view or convert the data in the several formats
    required through the move scheduling pipeline.

    The internal representation of the table data should not be directly accessed
    or modified, since there is some extra logic associated with correctly
    appending, inserting, validating, etc. Use the provided setters / getters
    instead.

    The initial starting position (in positioner local [theta,phi] coordinates)
    of the move table may be specified upon initialization (1x2 list). If it is
    not provided, then the move table will automatically look up the current
    expected position from the posmodel's state object.
    """

    def __init__(self, this_posmodel=None, init_posintTP=None):
        if not(this_posmodel):
            this_posmodel = posmodel.PosModel()
        self.posmodel = this_posmodel    # the particular positioner this table applies to
        self.posid = self.posmodel.posid # the string ID number of the positioner this table applies to
        self._log_note = ''              # optional note string which user can associate with this table, to be stored in any logging
        self.rows = []                   # internal representation of the move data
        self._rows_extra = []            # auto-generated backlash and final creep rows get internally stored here
        if init_posintTP:
            self.init_posintTP = init_posintTP  # initial [theta,phi] position (positioner local coordinates)
        else:
            self.init_posintTP = self.posmodel.expected_current_posintTP
        self.init_poslocTP = self.posmodel.trans.posintTP_to_poslocTP(self.init_posintTP)  # initial theta, phi position (in petal CS)
        self.should_antibacklash = self.posmodel.state._val['ANTIBACKLASH_ON']
        self.should_final_creep  = self.posmodel.state._val['FINAL_CREEP_ON']
        self.allow_cruise = not(self.posmodel.state._val['ONLY_CREEP'])
        self.allow_exceed_limits = self.posmodel.state._val['ALLOW_EXCEED_LIMITS']
        self._is_required = True
        self._postmove_cleanup_cmds = {pc.T: '', pc.P: ''}
        self._orig_command = ''
        self._warning_flag = 'WARNING'
        self._error_flag = 'ERROR'
        self._not_yet_calculated = '(not yet calculated)'

    def _set_zeno_dict(self, d):
        if self.posmodel.is_linphi:
            if 'zeno' in d:
                d['zeno'] += 'P'
            else:
                d['zeno'] = 'P'    # Denotes a movetable for a linear phi positioner
            d['PCCWA'] = float(self.posmodel.get_zeno_scale('SZ_CCW_P'))
            d['PCWA'] = float(self.posmodel.get_zeno_scale('SZ_CW_P'))
        return d

    def as_dict(self):
        """Returns a dictionary containing copies of all the table data."""
        c = self.copy()
        d = {'posid':                 c.posid,
             'log_note':              c.log_note,
             'rows':                  c.rows,
             '_rows_extra':           c._rows_extra,
             'init_posintTP':         c.init_posintTP,
             'init_poslocTP':         c.init_poslocTP,
             'should_antibacklash':   c.should_antibacklash,
             'should_final_creep':    c.should_final_creep,
             'allow_exceed_limits':   c.allow_exceed_limits,
             'allow_cruise':          c.allow_cruise,
             'postmove_cleanup_cmds': c._postmove_cleanup_cmds,
             'orig_command':          c._orig_command,
             'total_time':            self.total_time(suppress_automoves=False),
             'is_required':           c._is_required,
             }
        d = self._set_zeno_dict(d)
        return d

    def __repr__(self):
        return str(self.as_dict())

    def __str__(self):
        return str(self.as_dict())

    def display(self, printfunc=print, show_posid=True):
        '''Pretty-prints the table.  To return a string, instead of printing
        immediately, argue printfunc=None.
        '''
        def fmt(x):
            if x == None:
                x = str(x)
            if pc.is_string(x):
                return format(x,'>11s')
            elif pc.is_integer(x) or pc.is_float(x):
                return format(x,'>11g')
        tab = '  '
        output = f'move table for: {self.posid}\n' if show_posid else ''
        output += f'{tab}Original command: {self._orig_command}\n'
        output += f'{tab}Initial posintTP: {self.init_posintTP}\n'
        output += f'{tab}Initial poslocTP: {self.init_poslocTP}\n'
        for axisid, cmd_str in self._postmove_cleanup_cmds.items():
            output += f'{tab}Axis {axisid} postmove cmds: {repr(cmd_str)}\n'
        d = self.as_dict()
        keys_elsewhere = {'posid', 'rows', '_rows_extra', 'init_posintTP',
                          'init_poslocTP', 'postmove_cleanup_cmds', 'orig_command'}
        other_keys = set(d.keys()) - keys_elsewhere
        for key in other_keys:
            output += f'{tab}{key}: {d[key]}\n'
        if self.rows or self._rows_extra:
            output += fmt('row_type')
            headers = PosMoveRow().data.keys()
            for header in headers:
                output += fmt(header)
            for row in self.rows:
                output += '\n' + fmt('normal')
                for header in headers:
                    output += fmt(row.data[header])
            for extra_row in self._rows_extra:
                output += '\n' + fmt('extra')
                for header in headers:
                    output += fmt(extra_row.data[header])
        else:
            output += ' (empty: contains no row data)'
        if printfunc:
            printfunc(output)
        else:
            return output

    def display_for(self, output_type='hardware', printfunc=print):
        '''Pretty-prints the version that gets sent to hardware. To return a
        string, instead of printing immediately, argue printfunc=None.
        '''
        tab = '  '
        if output_type == 'hardware':
            t = self.for_hardware()
            t['row_time'] = [t['move_time'][i] + t['postpause'][i]/1000 for i in range(t['nrows'])]
        elif output_type == 'collider':
            t = self.for_collider()
        elif output_type == 'schedule':
            t = self.for_schedule()
        else:
            assert False, f'pretty printing of output_type {output_type} not yet defined'
        if 'row_time' not in t:
            t['row_time'] = [t['move_time'][i] + t['prepause'][i] + t['postpause'][i] for i in range(t['nrows'])]
        if 'net_time' not in t:
            t['net_time'] = [sum(t['row_time'][:i]) for i in range(1, t['nrows'] + 1)]
        output = f'move table for: {self.posid} ({output_type} version)'
        singletons = [k for k,v in t.items() if not isinstance(v, (list, tuple))]
        newline = f'\n{tab}'
        for key in singletons:
            output += f'\n{tab}{key}: {t[key]}'
        multiples = [k for k in t if k not in singletons]
        headers = [str(m) for m in multiples]
        widths = [max(8, len(h)) for h in headers]
        output += newline + tab.join([format(headers[i], f'^{widths[i]}s') for i in range(len(headers))])
        output += newline + tab.join(['-' * w for w in widths])
        lengths = {len(t[key]) for key in multiples}
        if len(lengths) != 1:
            output += '\n{tab}ERROR! NOT ALL COLUMNS HAVE SAME NUMBER OF ROWS!'
        for i in range(lengths.pop()):
            formats = []
            for key in headers:
                formats += [f'>{widths[headers.index(key)]}']
                value = t[key][i]
                if pc.is_integer(value):
                    formats[-1] += 'd'
                elif pc.is_float(value):
                    formats[-1] += '.3f'
            output += newline + tab.join([format(t[key][i], formats[multiples.index(key)]) for key in multiples])
        if printfunc:
            printfunc(output)
        else:
            return output

    def copy(self):
        new = copymodule.copy(self) # intentionally shallow, then will deep-copy just the row instances as needed below
        new.rows = [row.copy() for row in self.rows]
        new._rows_extra = [row.copy() for row in self._rows_extra]
        return new

    # getters
    def for_schedule(self, suppress_automoves=True):
        """Version of the table suitable for move scheduling. Distances are given at
        the output shafts, in degrees. Times are given in seconds.

        An option is provided to select whether auto-generated final creep and antibacklash
        moves should be suppressed from the schedule-formatted output table. For typical
        anticollision calcs, this option is left as True, since margin space for these
        moves is included in the geometric envelope of the positioner.
        """
        if suppress_automoves:
            return self._format_while_suppressing_automoves('schedule')
        return self._for_output_type('schedule')

    def for_collider(self, suppress_automoves=True):
        """Version of the table that is same as for_schedule, except with reduced
        amount of data returned (only what the poscollider requires and no more.)
        """
        if suppress_automoves:
            return self._format_while_suppressing_automoves('collider')
        return self._for_output_type('collider')

    def for_hardware(self):
        """Version of the table suitable for the hardware side.
        Distances are given at the motor shafts, in discrete steps.
        Times are given in milliseconds.
        """
        return self._for_output_type('hardware')

    def for_cleanup(self):
        """Version of the table suitable for updating the software internal
        position tracking after the physical move has been performed.
        """
        return self._for_output_type('cleanup')

    def angles(self):
        """Reduced version of the table giving just the theta and phi angles and
        deltas.
        """
        return self._for_output_type('angles')

    def full_table(self):
        """Version of the table with all data.
        """
        return self._for_output_type('full')

    def timing(self, suppress_automoves=False):
        '''Version of the table with just the time data.
        '''
        if suppress_automoves:
            return self._format_while_suppressing_automoves('timing')
        return self._for_output_type('timing')

    def _format_while_suppressing_automoves(self, output_type):
        '''Calculate version of table according to output_type, but don't include
        final creep or antibacklash moves.'''
        old_finalcreep = self.should_final_creep
        old_antibacklash = self.should_antibacklash
        self.should_final_creep = False
        self.should_antibacklash = False
        output = self._for_output_type(output_type)
        self.should_final_creep = old_finalcreep
        self.should_antibacklash = old_antibacklash
        return output

    @property
    def n_rows(self):
        """Number of rows in table.
        """
        return len(self.rows)

    @property
    def is_motionless(self):
        """Boolean saying whether the move table contains no motion at all, on
        neither the theta nor phi axis, in any row.
        """
        for row in self.rows:
            if row.has_motion:
                return False
        return True

    @property
    def has_phi_motion(self):
        """Boolean saying whether the move table contains any phi motion at all in any row.
        """
        for row in self.rows:
            if row.has_phi_motion:
                return True
        return False

    @property
    def log_note(self):
        '''Returns a copy of property log_note.'''
        return self._log_note

    @log_note.setter
    def log_note(self, note):
        '''Sets property log_note. The argument will be converted to str.'''
        self._log_note = str(note)

    def append_log_note(self, note):
        '''Appends a note to the current log_note.'''
        self._log_note = pc.join_notes(self._log_note, note)

    @property
    def error_str(self):
        '''Returns a string that is either empty or contains human-readable
        error messages, suitable for printout at a console or in a log.'''
        msg = ''
        i = 0
        for row in self.rows + self._rows_extra:
            auto_cmd = row.data['auto_cmd']
            is_err = self._warning_flag in auto_cmd or self._error_flag in auto_cmd
            if is_err:
                msg += f'\nRow {i}: {auto_cmd}'
            i += 1
        if msg:
            msg = f'{self.posid} move table contains errors/warnings:' + msg
        return msg

    def total_time(self, suppress_automoves=False):
        '''Returns total time to execute the table.'''
        times_table = self.timing(suppress_automoves=suppress_automoves)
        return times_table['net_time'][-1]

    def get_move(self, rowidx, axisid):
        ''' Returns distance for specified axis in row of move table '''
        dist_label = {pc.T:'dT_ideal', pc.P:'dP_ideal'}
        if rowidx >= len(self.rows):
            return None
        if axisid not in dist_label:
            return None
        return self.rows[rowidx].data[dist_label[axisid]]

    def get_prepause(self, rowidx):
        ''' Returns prepause (integer msec) for specified axis in row of move table '''
        if rowidx >= len(self.rows):
            return None
        return self.rows[rowidx].data['prepause']

    def get_postpause(self, rowidx):
        ''' Returns postpause (integer msec) for specified axis in row of move table '''
        if rowidx >= len(self.rows):
            return None
        return self.rows[rowidx].data['postpause']

    # setters
    def set_move(self, rowidx, axisid, distance):
        """Put or update a move distance into the table.
        If row index does not exist yet, then it will be added, and any blank
        filler rows will be generated in-between.
        """
        dist_label = {pc.T:'dT_ideal', pc.P:'dP_ideal'}
        if rowidx >= len(self.rows):
            self.insert_new_row(rowidx)
        self.rows[rowidx].data[dist_label[axisid]] = distance

    def store_orig_command(self, string, val1=None, val2=None):
        '''To keep a note of original move command associated with this move
        table. If any such notes already exist, the action is like an append.
        args val1 and val2 are optional, and are intended to represent a pair
        of coordinates.
        '''
        if val1 != None or val2 != None:
            string = f'{string}=[{val1}, {val2}]'
        self._orig_command = pc.join_notes(self._orig_command, string)

    def append_postmove_cleanup_cmd(self, axisid, cmd_str):
        """Add a posmodel cleanup command for execution after the move has
        been completed.
        """
        if cmd_str:
            separator = '\n'
            existing = self._postmove_cleanup_cmds[axisid]
            if existing and existing[-1] != separator:
                self._postmove_cleanup_cmds[axisid] += separator
            self._postmove_cleanup_cmds[axisid] += str(cmd_str)

    def set_prepause(self, rowidx, prepause):
        """Put or update a prepause into the table.
        If row index does not exist yet, then it will be added, and any blank
        filler rows will be generated in-between.
        """
        if rowidx >= len(self.rows):
            self.insert_new_row(rowidx)
        self.rows[rowidx].data['prepause'] = prepause

    def set_postpause(self, rowidx, postpause):
        """Put or update a postpause into the table.
        If row index does not exist yet, then it will be added, and any blank
        filler rows will be generated in-between.
        """
        if rowidx >= len(self.rows):
            self.insert_new_row(rowidx)
        self.rows[rowidx].data['postpause'] = postpause

    def set_required(self, boolean):
        '''Set flag which will be ultimately passed along to hardware. It tells
        the petalcontroller whether hardware *must* execute this table (to avoid
        collisions etc). If petalcontroller fails to send this table to its
        positioner, while required == True, then petalcontroller must return an
        error or reset power supplies, or some other such fallback plan, rather
        than just firing off the "execute tables" signal.
        '''
        self._is_required = pc.boolean(boolean)

    def strip(self):
        '''Removes two things from table:

            1. Any "zero" rows, i.e. with no motion and no pauses.
            2. Any pauses that come after the last finite move.

        Stripping is performed only on the user-defined rows, *not* on any
        internally auto-generated _rows_extra. In particular, case (2) means
        that any auto-creep moves will be pushed earlier in time, so that they
        occur immediately upon completion of the final user-defined motion.
        '''
        remove_pause = True
        for i in reversed(range(len(self.rows))):
            has_motion = self.rows[i].has_motion
            if remove_pause:
                if has_motion:
                    self.rows[i].data['postpause'] = 0
                    remove_pause = False
                else:
                    self.rows[i].data['prepause'] = 0
                    self.rows[i].data['postpause'] = 0
            has_prepause = self.rows[i].has_prepause
            has_postpause = self.rows[i].has_postpause
            if not has_motion and not has_prepause and not has_postpause:
                del self.rows[i]

    def compact(self):
        '''Removes two things from table:

            1. Any "zero" rows, i.e. with no motion and no pauses.
            2. Any pause-only rows (except the first)

        Compacting is performed only on the user-defined rows, *not* on any
        internally auto-generated _rows_extra. 
        '''
        for i in reversed(range(len(self.rows))):
            if i != 0:
                if not self.rows[i].has_motion:
                    prepause = self.get_prepause(i)
                    postpause = self.get_postpause(i)
                    if prepause or postpause:
                        prev_postpause = self.get_postpause(i-1)
                        prev_postpause += prepause + postpause
                        self.set_postpause(i-1, prev_postpause)
                    del self.rows[i]

    # row manipulations
    def insert_new_row(self, index):
        newrow = PosMoveRow()
        self.rows.insert(index, newrow)
        if index > len(self.rows):
            self.insert_new_row(index) # to fill in any blanks up to index

    def delete_row(self, index):
        del self.rows[index]

    def extend(self, other_move_table):
        """Extend one move table with another.
        Flags for antibacklash, creep, and limits remain the same as self, so if
        you want other_move_table to override these flags, this must be explicitly
        done separately.
        """
        if self == other_move_table:
            return
        for otherrow in other_move_table.rows:
            self.rows.append(otherrow.copy())
        for axisid, cmd_str in other_move_table._postmove_cleanup_cmds.items():
            self.append_postmove_cleanup_cmd(axisid=axisid, cmd_str=cmd_str)
        self.append_log_note(other_move_table.log_note)
        self.store_orig_command(string=other_move_table._orig_command)

        # Second table (since it comes last) takes precedence  when determining
        # whether to do automatic final moves.
        self.should_final_creep = other_move_table.should_final_creep
        self.should_antibacklash = other_move_table.should_antibacklash

        # Any disabling of range limits takes precedence.
        self.allow_exceed_limits |= other_move_table.allow_exceed_limits

        # Any forcing of creep-only takes precedence.
        self.allow_cruise &= other_move_table.allow_cruise

        # Any required table takes precedence.
        self._is_required |= other_move_table._is_required

    # internal methods
    def _calculate_true_moves(self):
        """Uses PosModel instance to get the real, quantized, calibrated values.
        Anti-backlash and final creep moves are added as necessary.
        """
        latest_TP = [x for x in self.init_posintTP]
        backlash = [0, 0]
        true_moves = [[], []]
        new_moves = [[], []]
        true_and_new = [[], []]
        has_moved = [False, False]
        normal_row_limits = None if self.allow_exceed_limits else 'debounced'
        extra_row_limits = None if self.allow_exceed_limits else 'near_full'
        axis_idxs = [pc.T, pc.P]
        for row in self.rows:
            ideal_dist = [row.data['dT_ideal'], row.data['dP_ideal']]
            for i in [pc.T,pc.P]:
                if self.posmodel.is_linphi and i == pc.P:
                    my_allow_cruise = True
                else:
                    my_allow_cruise = self.allow_cruise
                move = self.posmodel.true_move(axisid=i,
                                               distance=ideal_dist[i],
                                               allow_cruise=my_allow_cruise,
                                               limits=normal_row_limits,
                                               init_posintTP=latest_TP)
                true_moves[i].append(move)
                latest_TP[i] += true_moves[i][-1]['distance']
                if true_moves[i][-1]['distance']:
                    has_moved[i] = True
        if self.should_antibacklash and any(has_moved):
            backlash_dir = [self.posmodel.state._val['ANTIBACKLASH_FINAL_MOVE_DIR_T'],
                            self.posmodel.state._val['ANTIBACKLASH_FINAL_MOVE_DIR_P']]
            backlash_mag = self.posmodel.state._val['BACKLASH']
            for i in axis_idxs:
                if self.posmodel.is_linphi and i == pc.P:
                    backlash[i] = 0.0
                    auto_cmd_msg = ''
                    my_allow_cruise = False
                else:
                    backlash[i] = -backlash_dir[i] * backlash_mag * has_moved[i]
                    auto_cmd_msg = '(auto backlash backup)'
                    my_allow_cruise = self.allow_cruise
                move = self.posmodel.true_move(axisid=i,
                                               distance=backlash[i],
                                               allow_cruise=my_allow_cruise,
                                               limits=extra_row_limits,
                                               init_posintTP=latest_TP)
                new_moves[i].append(move)
                new_moves[i][-1]['auto_cmd'] = auto_cmd_msg
                latest_TP[i] += new_moves[i][-1]['distance']
        if self.should_final_creep or any(backlash):
            ideal_total = [0, 0]
            ideal_total[pc.T] = sum([row.data['dT_ideal'] for row in self.rows])
            ideal_total[pc.P] = sum([row.data['dP_ideal'] for row in self.rows])
            actual_total = [0, 0]
            err_dist = [0, 0]
            for i in axis_idxs:
                if self.posmodel.is_linphi and i == pc.P:
                    err_dist[i] = 0.0
                    auto_cmd_warning = ''
                else:
                    actual_total[i] = latest_TP[i] - self.init_posintTP[i]
                    if not self.posmodel.axis[i].is_locked:
                        err_dist[i] = ideal_total[i] - actual_total[i]
                    if abs(err_dist[i]) > pc.max_auto_creep_distance:
                        auto_cmd_warning = f' - {self._warning_flag}: auto creep distance={err_dist[i]:.3f} deg was ' + \
                                           f'truncated to pc.max_auto_creep_distance={pc.max_auto_creep_distance}, ' + \
                                            'which indicates likely upstream problem in the move schedule!'
                        err_dist[i] = pc.sign(err_dist[i]) * pc.max_auto_creep_distance
                    else:
                        auto_cmd_warning = ''
                move = self.posmodel.true_move(axisid=i,
                                               distance=err_dist[i],
                                               allow_cruise=False,
                                               limits=extra_row_limits,
                                               init_posintTP=latest_TP)
                new_moves[i].append(move)
                new_moves[i][-1]['auto_cmd'] = f'(auto final creep{auto_cmd_warning})'
                latest_TP[i] += new_moves[i][-1]['distance']
        self._rows_extra = []
        for i in range(len(new_moves[0])):
            self._rows_extra.append(PosMoveRow())
            self._rows_extra[i].data = {'dT_ideal': new_moves[pc.T][i]['distance'],
                                        'dP_ideal': new_moves[pc.P][i]['distance'],
                                        'prepause': 0,
                                        'postpause': 0,
                                        'auto_cmd': new_moves[pc.T][i]['auto_cmd'],
                                        }
        for i in [pc.T,pc.P]:
            true_and_new[i].extend(true_moves[i])
            true_and_new[i].extend(new_moves[i])
        return true_and_new

    def _for_output_type(self, output_type):
        """Internal function that calculates the various output table formats and
        passes them up to the wrapper functions above.
        """
        true_moves = self._calculate_true_moves()
        rows = self.rows.copy()
        rows.extend(self._rows_extra)
        row_range = range(len(rows))
        table = {}
        lock_note = ''
        if output_type in {'collider', 'schedule', 'full', 'cleanup', 'angles'}:
            table['dT'] = [true_moves[pc.T][i]['distance'] for i in row_range]
            table['dP'] = [true_moves[pc.P][i]['distance'] for i in row_range]
        if output_type in {'collider', 'schedule', 'full'}:
            table['Tdot'] = [true_moves[pc.T][i]['speed'] for i in row_range]
            table['Pdot'] = [true_moves[pc.P][i]['speed'] for i in row_range]
        if output_type in {'collider', 'schedule', 'full', 'timing'}:
            table['prepause'] = [rows[i].data['prepause'] for i in row_range]
            table['postpause'] = [rows[i].data['postpause'] for i in row_range]
        if output_type in {'hardware', 'full'}:
            table['motor_steps_T'] = [true_moves[pc.T][i]['motor_step'] for i in row_range]
            table['motor_steps_P'] = [true_moves[pc.P][i]['motor_step'] for i in row_range]
        if output_type in {'hardware', 'full', 'cleanup'}:
            table['speed_mode_T'] = [true_moves[pc.T][i]['speed_mode'] for i in row_range]
            table['speed_mode_P'] = [true_moves[pc.P][i]['speed_mode'] for i in row_range]
        if output_type in {'full', 'cleanup'}:
            table['orig_command'] = self._orig_command
            table['auto_commands'] = [rows[i].data['auto_cmd'] for i in row_range]
        if output_type in {'collider', 'schedule', 'full', 'hardware', 'timing'}:
            table['move_time'] = [max(true_moves[pc.T][i]['move_time'],
                                      true_moves[pc.P][i]['move_time']) for i in row_range]
        if output_type in {'full', 'cleanup', 'hardware'}:
            locked_axes = [name for num, name in {pc.T: 'T', pc.P: 'P'}.items() if self.posmodel.axis[num].is_locked]
            for axis in locked_axes:
                dkey = f'motor_steps_{axis}' if output_type == 'hardware' else f'd{axis}'
                if any(table[dkey]):
                    pc.printfunc(f'ERROR: {self.posid}: lock calculation failed on axis {axis}, non-zero {dkey} detected in' + \
                                   ' posmovetable.py. Move command to positioner will be altered to prevent any motion, though' + \
                                   ' could still result in timing errors or collisions.')
                    table[dkey] = [0 for i in row_range]
                    lock_note = pc.join_notes(lock_note, f'lock error on axis {axis}')
                if output_type in {'full', 'cleanup'}:
                    ideal_deltas = [rows[i].data[f'{dkey}_ideal'] for i in row_range]
                    zeroed_deltas = [table[dkey][i] == 0.0 and ideal_deltas[i] != 0.0 for i in row_range]
                    if any(zeroed_deltas):
                        lock_note = pc.join_notes(lock_note, f'locked {dkey}=0.0')  # ensures some indication of locking event gets into log for "expert" tables

        if output_type == 'hardware':
            table['posid'] = self.posmodel.posid
            table['canid'] = self.posmodel.canid
            table['busid'] = self.posmodel.busid
            table = self._set_zeno_dict(table)

            # interior rows
            table['postpause'] = [rows[i].data['postpause'] + rows[i+1].data['prepause'] for i in range(len(rows) - 1)]

            # last row
            table['postpause'].append(rows[-1].data['postpause'])

            # new first row, if necessary (because hardware only supports postpauses)
            leading_prepause = rows[0].data['prepause']
            if leading_prepause:
                table['postpause'].insert(0, leading_prepause)
                for key in ['motor_steps_T','motor_steps_P','move_time']:
                    table[key].insert(0, 0)
                for key in ['speed_mode_T','speed_mode_P']:
                    table[key].insert(0, 'creep') # speed mode doesn't matter here

            table['nrows'] = len(table['move_time'])
            table['total_time'] = sum(table['move_time'] + table['postpause']) # in seconds
            table['postpause'] = [int(round(x*1000)) for x in table['postpause']] # hardware postpause in integer milliseconds
            table['required'] = self._is_required
            self._latest_total_time_est = table['total_time']
            return table
        table['nrows'] = len(table['move_time']) if 'move_time' in table else len(table['dT'])
        if output_type == 'collider':
            return table
        table['posid'] = self.posmodel.posid
        if output_type in {'schedule', 'full', 'timing'}:
            table['net_time'] = [table['move_time'][i] + table['prepause'][i] + table['postpause'][i] for i in row_range]
            for i in range(1,table['nrows']):
                table['net_time'][i] += table['net_time'][i-1]
        if output_type == 'timing':
            return table
        if output_type in {'schedule', 'cleanup', 'full', 'angles'}:
            table['net_dT'] = table['dT'].copy()
            table['net_dP'] = table['dP'].copy()
            for i in range(1, table['nrows']):
                table['net_dT'][i] += table['net_dT'][i-1]
                table['net_dP'][i] += table['net_dP'][i-1]
        if output_type in {'cleanup', 'full'}:
            table.update({'TOTAL_CRUISE_MOVES_T':0,'TOTAL_CRUISE_MOVES_P':0,'TOTAL_CREEP_MOVES_T':0,'TOTAL_CREEP_MOVES_P':0})
            for i in row_range:
                table['TOTAL_CRUISE_MOVES_T'] += int(table['speed_mode_T'][i] == 'cruise' and table['dT'][i] != 0)
                table['TOTAL_CRUISE_MOVES_P'] += int(table['speed_mode_P'][i] == 'cruise' and table['dP'][i] != 0)
                table['TOTAL_CREEP_MOVES_T'] += int(table['speed_mode_T'][i] == 'creep' and table['dT'][i] != 0)
                table['TOTAL_CREEP_MOVES_P'] += int(table['speed_mode_P'][i] == 'creep' and table['dP'][i] != 0)
            table['postmove_cleanup_cmds'] = self._postmove_cleanup_cmds
            linphi_note = ''
            if self.posmodel.is_linphi:
                for s in [('CCW_SCALE_A','SZ_CCW_P'),('CW_SCALE_A','SZ_CW_P')]:
                    linphi_note = pc.join_notes(linphi_note, f'p{s[0]}={self.posmodel.get_zeno_scale(s[1])}')
            table['log_note'] = pc.join_notes(self.log_note, lock_note, linphi_note)
        if output_type in {'full', 'angles'}:
            trans = self.posmodel.trans
            posintT = [self.init_posintTP[pc.T] + table['net_dT'][i] for i in row_range]
            posintP = [self.init_posintTP[pc.P] + table['net_dP'][i] for i in row_range]
            table['posintTP'] = [[posintT[i], posintP[i]] for i in row_range]
            table['poslocTP'] = [trans.posintTP_to_poslocTP(tp) for tp in table['posintTP']]
        if output_type in {'full'}:
            table['poslocXY'] = [trans.posintTP_to_poslocXY(tp) for tp in table['posintTP']]
            table['QS'] = [trans.poslocXY_to_QS(xy) for xy in table['poslocXY']]
        return table

class PosMoveRow(object):
    """The general user does not directly use the internal values of a
    PosMoveRow instance, but rather should rely on the higher level table
    formats that are exported by PosMoveTable.
    """
    def __init__(self):
        self.data = {'prepause':  0,  # [sec] delay for this number of seconds before executing the move
                     'dT_ideal':  0,  # [deg] ideal theta distance to move (as seen by external observer)
                     'dP_ideal':  0,  # [deg] ideal phi distance to move (as seen by external observer)
                     'postpause': 0,  # [sec] delay for this number of seconds after the move has completed
                     'auto_cmd': '',  # [string] auto-generated command info corresponding to this row
                     }

    def __repr__(self):
        return repr(self.data)

    def copy(self):
        return copymodule.deepcopy(self)

    @property
    def has_motion(self):
        return self.data['dP_ideal'] != 0 or self.data['dT_ideal'] != 0

    @property
    def has_phi_motion(self):
        return self.data['dP_ideal'] != 0

    @property
    def has_prepause(self):
        return self.data['prepause'] != 0

    @property
    def has_postpause(self):
        return self.data['postpause'] != 0

