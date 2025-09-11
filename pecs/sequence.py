# -*- coding: utf-8 -*-
"""
Represents a sequence of positioner moves, with detailed control over motor
and move scheduling parameters.
"""

import os
from astropy.table import Table, vstack
import numpy as np
import datetime
import pandas
import sys
sys.path.append(os.path.abspath('../petal'))
import posconstants as pc

pos_defaults = {'CURR_SPIN_UP_DOWN': 100,
                'CURR_CRUISE': 100,
                'CURR_CREEP': 100,
                'CREEP_PERIOD': 2,
                'SPINUPDOWN_PERIOD': 12,
                'FINAL_CREEP_ON': True,
                'ANTIBACKLASH_ON': True,
                'ONLY_CREEP': False,
                'MIN_DIST_AT_CRUISE_SPEED': 180.0,
                'ALLOW_EXCEED_LIMITS': False,
                'BACKLASH': 3.0,
                'PRINCIPLE_HARDSTOP_CLEARANCE_T': 3.0,
                'PRINCIPLE_HARDSTOP_CLEARANCE_P': 3.0,
                'MOTOR_CCW_DIR_P': -1,
                'MOTOR_CCW_DIR_T': -1,
                }

pos_comments = {'CURR_SPIN_UP_DOWN': 'int, 0-100, spin up / spin down current',
                'CURR_CRUISE': 'int, 0-100, cruise current',
                'CURR_CREEP': 'int, 0-100, creep current',
                'CREEP_PERIOD': 'int, number of timer intervals corresponding to a creep step. a higher value causes slower creep',
                'SPINUPDOWN_PERIOD': 'int, number of 55 us periods to repeat each displacement during spin up to cruise speed or spin down from cruise speed. a higher value causes slower acceleration, over a longer travel distance',
                'FINAL_CREEP_ON': 'bool, if true do a finishing creep move after cruising',
                'ANTIBACKLASH_ON': 'boolean, if true do an antibacklash sequence at end of a move',
                'ONLY_CREEP': 'bool, if true disable cruising speed',
                'MIN_DIST_AT_CRUISE_SPEED': 'float, minimum rotor distance in deg to travel when at cruise speed before slowing back down',
                'ALLOW_EXCEED_LIMITS': 'bool, flag to allow positioner to go past software limits or not. exercise some caution if setting',
                'BACKLASH': 'deg, backlash removal distance',
                'PRINCIPLE_HARDSTOP_CLEARANCE_T': 'float, minimum distance in deg to stay clear of theta principle hardstop',
                'PRINCIPLE_HARDSTOP_CLEARANCE_P': 'float, minimum distance in deg to stay clear of phi principle hardstop',
                'MOTOR_CCW_DIR_P': '+1 or -1, defining as-wired motor counter-clockwise direction',
                'MOTOR_CCW_DIR_T': '+1 or -1, defining as-wired motor counter-clockwise direction',
                }

nominals = {}
nominals['gear_ratio'] = ((46+14)/14)**4
nominals['timer_update_rate'] = 18e3 # Hz
nominals['stepsize_creep'] = 0.1 # deg
nominals['stepsize_cruise'] = 3.3 # deg
nominals['motor_speed_cruise'] = 9900.0 * 360.0 / 60.0 # deg/sec (= RPM *360/60)
nominals['spinupdown_dist_per_period'] = sum(range(round(nominals['stepsize_cruise']/nominals['stepsize_creep']) + 1))*nominals['stepsize_creep']
nominals['spinupdown_distance'] = nominals['spinupdown_dist_per_period'] * pos_defaults['SPINUPDOWN_PERIOD']
nominals['spinupdown_distance_output'] = nominals['spinupdown_distance'] / nominals['gear_ratio']

global_commands = {'QS', 'obsXY', 'ptlXY'}
abs_local_commands = {'poslocXY', 'poslocTP', 'posintTP'}
delta_commands = {'dQdS', 'obsdXdY', 'poslocdXdY', 'dTdP'}
abs_commands = global_commands | abs_local_commands
general_commands = abs_commands | delta_commands
homing_commands = {'home_and_debounce', 'home_no_debounce'}
local_commands = abs_local_commands | delta_commands | homing_commands
valid_commands = general_commands.copy()
valid_commands.update(homing_commands)

is_float = lambda x: isinstance(x, (float, np.floating))
is_int = lambda x: isinstance(x, (int, np.integer))
is_number = lambda x: is_float(x) or is_int(x)
is_bool = lambda x: isinstance(x, (bool))
is_string = lambda x: isinstance(x, (str))

move_idx_key = 'move_idx'
get_datestr = lambda: datetime.datetime.now().isoformat(sep=' ', timespec='seconds')
sequence_note_prefix = 'sequence: '

class Sequence(object):
    '''Iterable structure that defines a positioner test, as a sequence of Move instances.
    Typical operations you can do on a list work here (indexing, slices, del, etc).
    However, all elements must have type Move (see class lower in this module).
    
        short_name   ... string, brief name for the test
        long_name    ... string, optional longer descriptive name for the test
        details      ... string, optional additional string to explain the test's purpose, etc
        pos_settings ... optional dict of positioner settings to apply during the sequence
        
    After initialization, populate the sequence using the "add_move" function.
    '''
    def __init__(self, short_name, long_name='', details='', pos_settings=None):
        self.moves = []
        self.short_name = short_name
        self.long_name = str(long_name)
        self.details = str(details)
        self.creation_date = get_datestr()
        pos_settings = {} if pos_settings is None else pos_settings
        self.pos_settings = pos_settings
        self._max_print_lines = np.inf
    
    save_as_metadata = ['short_name', 'long_name', 'details', 'creation_date', 'pos_settings']

    @property
    def normalized_short_name(self):
        '''A particular format, used for example when saving files to disk.'''
        return self.short_name.replace(' ','_').upper()
    
    @staticmethod
    def read(path):
        '''Reads in and validates format for a saved Sequence from a file. E.g.
        sequence = Sequence.read(path)
        '''
        path = os.path.abspath(path)
        kwargs = {'fill_values': None} if os.path.splitext(path)[1] == '.ecsv' else {}
        table = Table.read(path, **kwargs)
        table = Sequence._validate_table(table)
        basename = os.path.basename(path)
        fallbacks = {'short_name': os.path.splitext(basename)[0],
                     'long_name': '',
                     'details': f'Read in from {path} on {get_datestr()}',
                     'pos_settings': {},
                     }
        seq = Sequence(**fallbacks)
        seq.creation_date = 'unknown_date'
        for key, value in table.meta.items():
            setattr(seq, key, value)  # overrides the fallback values where possible
        move_idxs = sorted(set(table[move_idx_key]))
        for m in move_idxs:
            subtable = table[m == table[move_idx_key]]
            kwargs = {key: subtable[key][0] for key in Move.single_keys}
            kwargs.update({key: subtable[key].tolist() for key in Move.multi_keys})
            if len(subtable) == 1:
                for key, val in kwargs.items():
                    if isinstance(val, list) and len(val) == 1:
                        kwargs[key] = val[0]
            move = Move(**kwargs)
            seq.append(move)
        return seq
    
    def save(self, directory='.', basename=''):
        '''Saves an ecsv file representing the sequence to directory/basename.ecsv
        Returns the path of the saved file.
        '''
        if not basename:
            basename = f'seq_{self.normalized_short_name}'
        table = self.to_table()
        path = os.path.join(directory, basename + '.ecsv')
        table.write(path, overwrite=True, delimiter=',')
        return path
    
    def to_table(self):
        '''Generates an astropy table which completely summarizes the sequence.
        '''
        tables = []
        for i in range(len(self)):
            move = self.moves[i]
            this_dict = {move_idx_key: [i] * len(move)}
            this_dict.update(move.to_dict(sparse=False))
            this_table = Table(this_dict)
            tables.append(this_table)
        table = vstack(tables)
        table.meta = {key: getattr(self, key) for key in Sequence.save_as_metadata}
        return table
    
    def get_posids(self):
        '''Returns set of all posid strings specified in the sequence.
        '''
        posids = set()
        for move in self:
            posids |= set(move.posids)
        if 'any' in posids:
            posids.remove('any')
        return posids
    
    @property
    def pos_settings(self):
        '''Copy of internal dict, containing all pos settings. Will be some
        subset of key/value pairs like those in pos_defaults.'''
        return self._pos_settings.copy()
    
    @pos_settings.setter
    def pos_settings(self, settings):
        '''Settings dict, may contain any subset of pos_defaults keys.'''
        self._pos_settings = self._validate_pos_settings(settings)
    
    def merge(self, other):
        '''Return a new sequnce, which merges this one with another. A number of
        restrictions apply. The two sequences must be:
            - same number of moves
            - same pos_settings
            - within each move, must have:
                - unique posids
                - same command
                - same allow_corr
        '''
        assert isinstance(other, Sequence)
        new = Sequence(short_name=pc.join_notes(self.short_name, other.short_name),
                       long_name=pc.join_notes(self.long_name, other.long_name),
                       details=pc.join_notes(self.details, other.details),
                       )
        assert len(self) == len(other), 'merge of non-equal length sequences is not defined'
        assert self.pos_settings == other.pos_settings, 'merge of differing pos_settings is not defined'
        for i in range(len(self)):
            A = self[i]
            B = other[i]
            assert not A.is_uniform and not B.is_uniform, 'merge of move with undifferentiated posids is not defined'
            for key in Move.single_keys:
                assert getattr(A, key) == getattr(B, key), f'merge of moves with differing values for "{key}" is not defined'
            assert len(set(A.posids) & set(B.posids)) == 0, 'merge of moves with overlapping posids is not defined'
            move = Move(command=A.command,
                        target0=A.target0 + B.target0,
                        target1=A.target1 + B.target1,
                        posids=A.posids + B.posids,
                        log_note=A._log_note + B._log_note,
                        allow_corr=A.allow_corr
                        )
            new.append(move)
        return new
    
    def str(self, max_lines=np.inf):
        '''Returns same as __str__ except limits printed table length to max_lines.'''
        old_max_print_lines = self._max_print_lines
        self._max_print_lines = max_lines
        s = self.__str__()
        self._max_print_lines = old_max_print_lines
        return s
    
    def __str__(self):
        s = self._meta_str()
        s += '\n'
        table = self.to_table()
        align = []
        for col in table.columns:
            if col.lower() in {'log_note'}:
                align += ['<']
            elif col.lower() in {'command', 'allow_corr'}:
                align += ['^']
            else:
                align += ['>']
        lines = table.pformat_all(align=align)
        overage = len(lines) - self._max_print_lines
        if overage > 0:
            tail = int(np.floor(self._max_print_lines/2))
            head = int(self._max_print_lines - tail)
            lines = lines[:head] + ['...'] + lines[-tail:]
        s += '\n'.join(lines)
        if len(lines) > 50:
            # repeat headers and metadata for convenience with long tables
            s += f'\n{lines[1]}\n{lines[0]}\n'
            s += f'\n{self._meta_str()}'
        return s
    
    def __repr__(self):
        s = object.__repr__(self)
        s += '\n'
        s += self._meta_str()
        return s
        
    def _meta_str(self):
        '''Returns a few lines of descriptive meta data.'''
        s = f'{self.normalized_short_name}'
        if self.long_name:
            s += f': {self.long_name}'
        if self.creation_date:
            s += f'\nCreated: {self.creation_date}'
        if self.details:
            s += f'\n{self.details}'
        s += f'\nSettings: {self.pos_settings}'
        s += f'\nNumber of rows = {len(self.to_table())}'
        s += f'\nNumber of moves = {len(self)}'
        return s
    
    def _validate_move(self, move):
        assert isinstance(move, Move), f'{move} must be an instance of Move class'
        move.register_sequence_name_getter(lambda: self.normalized_short_name)
    
    def _validate_table(table):
        '''Validates an astropy table representing a sequence. Returns a new
        table, which may be slightly modified. For example, adding any new columns
        which have been added to sequence since that table was saved to disk, or
        putting in placeholder metadata.
        '''
        new = table.copy()
        example_move = Move(command='dTdP', target0=0.0, target1=0.0)
        example_dict = example_move.to_dict()
        example_table = Table(example_dict)
        exclude_columns = {move_idx_key}
        test_columns = set(new.columns) - exclude_columns
        for col in test_columns:
            assert col in example_table.columns
            for i in range(len(new)):
                val = new[col][i]
                test = example_table[col][0]
                if is_number(val) and is_number(test) or is_bool(val) and is_bool(test):
                    continue
                assert isinstance(val, type(test))
            try:
                np.isfinite(example_table[col][0])
                isnumber = True
            except:
                isnumber = False
            if isnumber:
                assert all(np.isfinite(new[col]))
        for row in new:
            assert row['command'] in valid_commands
        missing_col = set(example_table.columns) - set(new.columns)
        defaults = example_move.to_dict(sparse=True)
        for k, v in defaults.items():
            if isinstance(v, list):
                defaults[k] = v[0]
        defaults.update(pos_defaults)
        missing_col_default_available = missing_col & set(defaults)
        for col in missing_col_default_available:
            values = [defaults[col]] * len(new)
            new[col] = values
        if move_idx_key not in new.columns:
            new[move_idx_key] = range(len(new))
        return new
    
    @staticmethod
    def _validate_pos_settings(settings):
        '''Returns a new, validated pos_settings dict, with possible type casts
        of values.'''
        new = settings.copy()
        for key, value in settings.items():
            assert key in pos_defaults
            example = pos_defaults[key]
            expected_type =  type(example)
            if is_bool(value):
                new[key] = bool(value)
            elif is_number(value) and is_float(example):
                    new[key] = float(value)
            elif is_int(value):
                new[key] = int(value)
            assert isinstance(new[key], expected_type)
        return new
    
    # Standard sequence functions
    def __iter__(self):
        self.__idx = 0
        return self
        
    def __next__(self):
        if self.__idx < len(self):
            move = self.moves[self.__idx]
            self.__idx += 1
            return move
        else:
            raise StopIteration
    
    def __len__(self):
        return len(self.moves)
    
    def __getitem__(self, key):
        return self.moves[key]
        
    def __setitem__(self, key, value):
        self._validate_move(value)
        self.moves[key] = value

    def __delitem__(self, key):
        del self.moves[key]
        
    def __contains__(self, value):
        return value in self.moves
    
    def index(self, value):
        for i in range(len(self)):
            if self[i] == value:
                return i
        raise ValueError(f'{object.__repr__(self)} is not in sequence')
    
    def count(self, value):
        matches = [True for x in self if value == x]
        return len(matches)
    
    def append(self, value):
        self._validate_move(value)
        self.moves.append(value)

    
class Move(object):
    '''Encapsulates the command data for a simultaneous movement of one or more
    fiber positioners.
    
        command      ... move command string as described in petal.request_targets()
        target0      ... 1st target coordinate(s) or delta(s), as described in petal.request_targets()
        target1      ... 2nd target coordinate(s) or delta(s), as described in petal.request_targets()
        posids       ... optional sequence of which positioner ids should be commanded
        log_note     ... optional string (uniformly applied) or sequence of strings (per device) to store alongside log data in this move
        allow_corr   ... optional boolean, whether correction moves are allowed to be performed after the primary ("blind") move
        
    Multiple positioners, same targets:
        By default, posids='any', which means that the move command should be applied
        to any positioners selected at runtime (presumably by the PECS initialization
        process). In this case, all positioners get the same target0 and target1. The
        command string must be one of local_commands, as defined in this module.
        
    Multiple positioners, multiple targets:
        By arguing ordered sequences for target0, target1, and posids, the user
        can specify different targets for different positioners. The command will be
        the same in all cases, and must be one of general_commands, as defined in this
        module.
        
    Homing commands:
        Set target0 = 1 if you want to home theta axis, target1 = 1 to home phi axis,
        or set both to home both axes.
    '''
    def __init__(self, command, target0, target1, posids='any', log_note='', allow_corr=True):
        if is_string(posids):
            assert command in local_commands, f'cannot apply a non-local command {command} to multiple positioners'
            for x in [target0, target1]:
                assert is_number(x), f'target {x} is not a number'
            self.target0 = [target0]
            self.target1 = [target1]
            self.posids = [posids]           
        else:
            assert command in general_commands, f'command {command} not recognized, see general_commands collection'
            assert (len(target0) == len(target1) == len(posids)), 'args target0, target1, and posid are not of equal length'
            self.target0 = list(target0)
            self.target1 = list(target1)
            self.posids = list(posids)
        self.log_note = log_note
        for X in [self.target0, self.target1]:
            assert all([is_number(x) for x in X]), 'all elements of target0 or target1 must be numbers'
        if command in homing_commands:
            assert target0 or target1, 'for a homing command, need to set target0 to 1 for theta ' + \
                                       'homing, target1 to 1 for phi homing, or both simultaneously'
        self.command = command
        self.allow_corr = bool(allow_corr)
        self._sequence_name_getter = None
    
    # properties of Move with single values vs multiple values
    single_keys = {'command', 'allow_corr'}
    multi_keys = {'target0', 'target1', 'posids', 'log_note'}
    
    def register_sequence_name_getter(self, function):
        '''Register function that can get name of the sequence that contains this move.'''
        self._sequence_name_getter = function
    
    @property
    def is_uniform(self):
        '''Boolean whether the same uniform move command and target is to be done on
        any positioner.'''
        return self.posids == ['any'] or len(self.posids) <= 1
    
    @property
    def log_note(self):
        '''Log note strings, including sequence name if available.'''
        notes = []
        for note in self._log_note:
            existing_parts = [part.strip() for part in note.split(';')]
            new_parts = []
            has_name_yet = False
            for part in existing_parts:
                is_name = sequence_note_prefix in part
                not_a_name = sequence_note_prefix not in part
                use_old_name = is_name and not has_name_yet and not self._sequence_name_getter
                if not_a_name or use_old_name:
                    new_parts.append(part)
                if use_old_name:
                    has_name_yet = True
            if not has_name_yet and self._sequence_name_getter:
                seq_name = str(sequence_note_prefix) + str(self._sequence_name_getter())
                new_parts = [seq_name] + new_parts
            note = pc.join_notes(*new_parts)
            notes.append(note)
        return notes
    
    @log_note.setter
    def log_note(self, note):
        '''Set with a single string (will be uniformly applied in case of multiple targets)
        or a list of strings, of same length as number of targets.'''
        if is_string(note):
            self._log_note = [note] * len(self.posids)
        else:
            assert isinstance(note, (list, tuple)), f'invalid collection type {type(note)} for log note'
            assert len(note) == len(self.posids)
            self._log_note = [str(x) for x in note]
    
    def get_log_notes(self, posids='any'):
        '''Returns list of log note values for the collection posids, in
        the same order as posids.'''
        notes = self.log_note
        if posids != 'any':
            assert isinstance(posids, (list, tuple)), f'posids={posids} of type {type(posids)} is not supported for get_log_notes()'
        if posids == 'any':
            posids = self.posids
        elif self.is_uniform:
            return [notes[0]] * len(posids)
        return [notes[self.posids.index(posid)] for posid in posids]        
    
    def __len__(self):
        return len(self.posids)
    
    @property
    def has_multiple_targets(self):
        '''Boolean whether this move has mulitple targets.'''
        return len(self) > 1
    
    def is_defined_for_all_positioners(self, posids):
        '''Returns boolean whether the move has valid command definitions for
        the argued collection of posids.
        '''
        if 'any' in self.posids:
            return True
        missing = set(posids) - set(self.posids)
        if any(missing):
            return False
        return True
    
    def to_dict(self, sparse=False, posids='any'):
        '''Returns a dict, which is ready for direct conversion to an astropy
        table or pandas dataframe. Keys = column names and values = equal-length
        lists of column data.
        
        INPUTS:  sparse ... Optional boolean, if True, will make a more sparse
                            representation (non-arrays not expanded). This form
                            is not suitable for direct conversion to table.
                            
                 posids ... Optional set, tuple, or list of posids. Any
                            arrays in the returned dictionary will only contain
                            entries corresponding to these positioners.
        '''
        data = {}
        key_order = ['command', 'target0', 'target1', 'posids', 'allow_corr', 'log_note']
        for key in self.single_keys:
            data[key] = [getattr(self, key)] * len(self)
        for key in self.multi_keys:
            data[key] = getattr(self, key)
        key_order = key_order + sorted(set(data) - set(key_order))
        data = {key: data[key] for key in key_order}
        if posids != 'any':
            ids = sorted(posids)
            if self.is_uniform:
                selection = [0] if sparse else [0 for i in range(len(ids))]
            else:
                selection = [i for i in range(len(ids)) if ids[i] in self.posids]
            for key, value in data.items():
                data[key] = [value[i] for i in selection]
        if sparse:
            for key in self.single_keys:
                data[key] = data[key][0] if len(data[key]) > 0 else []
            for key in self.multi_keys:
                if len(set(data[key])) == 1:
                    data[key] = data[key][0]
        return data
    
    def make_request(self, posids, log_note=''):
        '''Make a move request data structure, ready for sending to the online control system.
        Any positioners known to this move instance, but not included in posids, will be
        skipped.
        
        INPUT:  posids ... collection of positioner ids
                log_note ... optional string, will be appended to any existing log note
        
        OUTPUT: pandas dataframe with columns 'DEVICE_ID', 'COMMAND', 'X1', 'X2', 'LOG_NOTE'
        '''
        assert self.command not in homing_commands, 'Cannot make a request for homing. Try make_homing_kwargs() instead.'
        posids = list(set(posids))
        sorted_posids, target0, target1, final_log_note = [], [], [], []
        possible_notes = [pc.join_notes(note, log_note) for note in self.log_note]
        if self.has_multiple_targets:
            for i in range(len(self)):
                posid = self.posids[i]
                if posid in posids:
                    sorted_posids += [posid]
                    target0 += [self.target0[i]]
                    target1 += [self.target1[i]]
                    final_log_note += [possible_notes[i]]
        else:
            sorted_posids = posids
            target0 = self.target0[0]
            target1 = self.target1[0]
            final_log_note = possible_notes[0]
        request_data = {'DEVICE_ID': sorted_posids,
                        'COMMAND': self.command,
                        'X1': target0,
                        'X2': target1,
                        'LOG_NOTE': final_log_note,
                        }
        request = pandas.DataFrame(request_data)
        return request
        
    def make_homing_kwargs(self, posids, log_note=''):
        '''Make a kwargs dictionary for homing moves, ready for sending to the online
        control system.
        
        INPUT:  posids ... positioners to home
                log_note ... optional string, will be appended to any existing log note
        
        OUTPUT: kwargs dictionary
        '''
        assert self.command in homing_commands, 'Cannot make homing kwargs for non-homing. Try make_request() instead.'
        should_debounce = self.command == 'home_and_debounce'
        axis = 'theta_only' if not self.target1[0] else 'phi_only' if not self.target0[0] else 'both'
        possible_notes = [pc.join_notes(note, log_note) for note in self.log_note]
        if self.has_multiple_targets:
            log_note = possible_notes
        else:
            log_note = possible_notes[0]
        kwargs = {'posids': posids,
                  'axis': axis,
                  'debounce': should_debounce,
                  'log_note': log_note,
                  }
        return kwargs
        
    def __str__(self):
        s = object.__repr__(self)
        s += '\n'
        s += str(self.to_dict(sparse=True))
        return s
    
    def __repr__(self):
        return self.__str__()
