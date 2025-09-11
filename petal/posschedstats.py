import os
import io
import numpy as np
import pandas as pd
import csv
import posconstants as pc
import copy


# separated out here because so long
final_checks_str = 'num collisions found by final check (should be zero, or None if no check done)'
found_not_registered_str = 'found but not directly resolved (useful for debugging, not necessarily empty -- some collisions may be indirectly resolved)'
_unregistered_schedule_str = 'unregistered'
_blank_str = '-'

class PosSchedStats(object):
    """Collects statistics from runs of the PosSchedule.
    
    By default, the module starts out disabled. You can override this with a
    boolean value True to initialization argument enabled. Or use the enable()
    and disable() methods separately.
    """
    def __init__(self, enabled=False):
        self._init_data_structures()
        self.disable() # called here in all cases, to ensure consistent initialization of state
        if enabled:
            self.enable()
        self.clear_cache_after_save = True
        self.filename_suffix = ''
    
    def is_enabled(self):
        '''Boolean whether statistics tracking is currently turned on.'''
        if self._is_enabled and self.schedule_ids: # 2nd expression is to avoid storing to improperly empty structure (i.e. no schedule yet registered) 
            return True
        return False
    
    def enable(self):
        '''Turn module on, ready for statistics tracking. Any existing
        data in the cache will be cleared.'''
        carryover_npos = None
        if self.n_rows > 0:
            carryover_npos = self.numbers['n pos'][-1] # for cases where enabling stats while some posschedule instance (unregistered) already exists
        self._init_data_structures(num_pos=carryover_npos)
        self._is_enabled = True
    
    def disable(self):
        '''Turn module off, disabled for any statistics tracking.'''
        self._is_enabled = False
    
    @property
    def n_rows(self):
        '''Number of rows present in the data structure, considering each
        registered schedule to be a 'row'.'''
        return len(self.schedule_ids)
    
    @property
    def latest(self):
        '''Get latest 'row', i.e. latest registered schedule id, which acts
        as an index to that row.'''
        return self.schedule_ids[-1]
    
    def _init_data_structures(self, num_pos=None):
        '''Initialize data structures.
        
        Note the argument num_pos may be entered in special cases where the
        number of positioners is already known. This value will be put into the
        first (unregistered) row of the data structure. In most cases no need
        to worry about this --- it occurs whenever registering a new schedule.
        '''
        self.schedule_ids = []
        self.collisions = {}
        self.unresolved = {}
        self.unresolved_tables = {}
        self.unresolved_sweeps = {}
        self.final_checked_collision_pairs = {}
        self.strings = {'method':[], 'note':[]}
        self._strings_to_print_first = ['method']
        self._strings_to_print_last = ['note']
        self.numbers = {'n pos':[],
                        'n requests':[],
                        'n requests accepted':[],
                        'n move tables':[],
                        'n tables achieving requested-and-accepted targets':[],
                        final_checks_str:[],
                        'max table move time':[],
                        'num path adjustment iters':[],
                        'request_target calc time':[],
                        'schedule_moves calc time':[],
                        'request + schedule calc time':[],
                        'expert_add_table calc time':[],
                        'max num table rows':[],
                        'avg num table rows':[],
                        'std num table rows':[],
                        }
        self.avoidances = {}
        self._latest_saved_row = None
        dummy_id = pc.timestamp_str() + ' (' + _unregistered_schedule_str + ')'
        self.register_new_schedule(schedule_id=dummy_id, num_pos=num_pos)
    
    def register_new_schedule(self, schedule_id, num_pos=None):
        """Register a new schedule object to track statistics of.
        
           schedule_id ... unique string
           num_pos     ... number of positioners
        """
        self.schedule_ids.append(schedule_id)
        self.collisions[self.latest] = {'found':set(), 'resolved':{}}
        self.unresolved[self.latest] = {}
        self.unresolved_tables[self.latest] = {}
        self.unresolved_sweeps[self.latest] = {}
        self.final_checked_collision_pairs[self.latest] = {}
        for key in self.strings:
            self.strings[key].append(_blank_str)
        for key in self.numbers:
            val = None if key == final_checks_str else 0
            self.numbers[key].append(val)
        self.numbers['n pos'][-1] = num_pos
        self.avoidances[self.latest] = {}

    def set_num_move_tables(self, n):
        """Record number of move tables in the current schedule."""
        self.numbers['n move tables'][-1] = n
        
    def add_hardware_move_tables(self, tables):
        """Records stats (not necessarily each whole table) on tables."""
        if not tables:
            return
        lengths = [t['nrows'] for t in tables]
        self.numbers['max num table rows'][-1] = max(lengths)
        self.numbers['avg num table rows'][-1] = np.mean(lengths)
        self.numbers['std num table rows'][-1] = np.std(lengths)
        
    def add_collisions_found(self, collision_pair_ids):
        """Add collisions that have been found in the current schedule.
            collision_pair_ids ... set of colliding posids
        """
        these_collisions = self.collisions[self.latest]
        these_collisions['found'] = these_collisions['found'].union(collision_pair_ids)
        
    def add_collisions_resolved(self, posid, method, collision_pair_ids):
        """Add collisions that have been resolved in the current schedule.
            posid ... string, the positioner that was adjusted
            method ... string, collision resolution method
            collision_pair_ids ... set of collision_pair_ids resolved by method
        """
        latest = self.latest
        this_dict = self.collisions[latest]['resolved']
        if method not in this_dict:
            this_dict[method] = set()
        this_dict[method] = this_dict[method].union(collision_pair_ids)
        for pair in collision_pair_ids:
            self.add_avoidance(posid, method, pair)
            
    def add_avoidance(self, posid, method, collision_pair_id):
        '''Add an avoidance method for a collision that has been resolved in the
        current schedule.
            posid ... string, the positioner that was adjusted
            method ... string, collision resolution method
            collision_pair_id ... collision_pair_id resolved by method, posid must be one of the pair!
        '''
        split = set(collision_pair_id.split('-'))
        assert posid in split, f'posschedstats: {posid} not found in collision pair {collision_pair_id}!'
        other = (split - {posid}).pop()
        latest = self.latest
        if posid not in self.avoidances[latest]:
            self.avoidances[latest][posid] = []
        self.avoidances[latest][posid] += [f'{method}/{other}']
        
    def get_collisions_resolved_by(self, method='freeze'):
        '''Returns set of all collisions resolved by method in the latest
        schedule.'''
        this_dict = self.collisions[self.latest]['resolved']
        if method in this_dict:
            return this_dict[method].copy()
        return {}
    
    def get_avoidances(self, posid):
        '''Returns collection of all collision avoidances for that positioner
        in the latest schedule. Will be empty if none found.'''
        latest = self.latest
        out = []
        if posid in self.avoidances[latest]:
            out = self.avoidances[latest][posid]
        return out
    
    @property
    def total_unresolved(self):
        '''Returns count of how many collisions were recorded as unresolved after
        final collision check. For the latest schedule_id.
        '''
        value = self.numbers[final_checks_str][-1]
        if value == None:
            return 0
        return value
    
    @property
    def total_resolved(self):
        '''Returns count of how many collisions were recorded as resolved.
        '''
        this_dict = self.collisions[self.latest]['resolved']
        count = 0
        for collision_pairs in this_dict.values():
            count += len(collision_pairs)
        return count
    
    @property
    def unresolved_posids(self):
        '''Returns set of posids that were recorded as unresolved after final
        collision check. For the latest schedule_id.
        '''
        unresolved = self.unresolved[self.latest]
        if 'final' in unresolved:
            return unresolved['final']
        return {}
        
    def add_request(self):
        """Increment requests count."""
        self.numbers['n requests'][-1] += 1
        
    def add_request_accepted(self):
        """Increment requests accepted count."""
        self.numbers['n requests accepted'][-1] += 1
    
    def sub_request_accepted(self):
        """Deccrement requests accepted count."""
        if self.numbers['n requests accepted'][-1] > 0:
            self.numbers['n requests accepted'][-1] -= 1
    
    def add_table_matching_request(self):
        """Increment number of tables matching their original requests."""
        self.numbers['n tables achieving requested-and-accepted targets'][-1] += 1
            
    def add_requesting_time(self, time):
        """Add in time spent processing target requests in the current schedule."""
        self.numbers['request_target calc time'][-1] += time
        self.numbers['request + schedule calc time'][-1] += time

    def add_scheduling_time(self, time):
        """Add in time spent processing move schedules in the current schedule."""
        self.numbers['schedule_moves calc time'][-1] += time
        self.numbers['request + schedule calc time'][-1] += time

    def add_expert_table_time(self, time):
        """Add in time spent putting expert tables into the current schedule."""
        self.numbers['expert_add_table calc time'][-1] += time

    def set_scheduling_method(self, method):
        """Set scheduling method (a string) for the current schedule."""
        self.strings['method'][-1] = str(method)
    
    def add_note(self, note):
        """Add a note string for the current schedule. If one already exists,
        then the argued note will be appended to it, with a standard separator.
        """
        if not note:
            return
        old = self.strings['note'][-1]
        if old == _blank_str:
            new = note
        else:
            new = pc.join_notes(old, note)
        self.strings['note'][-1] = new
        
    def set_max_table_time(self, time):
        """Set the maximum move table time in the current schedule."""
        self.numbers['max table move time'][-1] = time

    def add_to_num_adjustment_iters(self, iterations):
        """Add data recording number of iterations of path adjustment were made."""
        self.numbers['num path adjustment iters'][-1] += iterations
        
    def add_final_collision_check(self, collision_pairs):
        """Add data recording if there were still any bots colliding after a
        final check."""
        if self.numbers[final_checks_str][-1] == None:
            self.numbers[final_checks_str][-1] = 0
        self.numbers[final_checks_str][-1] += len(collision_pairs)
        self.final_checked_collision_pairs[self.latest] = collision_pairs
    
    def add_unresolved_colliding_at_stage(self, stage_name, colliding_set, colliding_tables, colliding_sweeps):
        """Add data recording ids of any bots colliding in a given stage. This
        is distinct from "final" check data."""
        if stage_name not in self.unresolved[self.latest]:
            self.unresolved[self.latest][stage_name] = set()
            self.unresolved_tables[self.latest][stage_name] = {}
            self.unresolved_sweeps[self.latest][stage_name] = {}
        colliding_tables_copies = {} # stores copy-able / pickle-able versions, rather thatn the complete posmovetable instances
        for key,table in colliding_tables.items():
            colliding_tables_copies[key] = table.as_dict()
        self.unresolved[self.latest][stage_name].update(colliding_set)
        self.unresolved_tables[self.latest][stage_name].update(colliding_tables_copies)
        self.unresolved_sweeps[self.latest][stage_name].update(colliding_sweeps)
    
    def summarize_collision_resolutions(self):
        """Returns a summary dictionary of the collisions and resolutions data."""
        summary = {}
        summary['found total collisions (before anticollision)'] = [len(self.collisions[sched]['found']) for sched in self.schedule_ids]
        summary['collisions with PTL'] = [len({c for c in self.collisions[sched]['found'] if 'PTL' in c}) for sched in self.schedule_ids]
        summary['collisions with GFA'] = [len({c for c in self.collisions[sched]['found'] if 'GFA' in c}) for sched in self.schedule_ids]
        summary['resolved total collisions'] = []
        summary['found set'] = []
        summary['resolved set'] = []
        summary[found_not_registered_str] = []
        for sched in self.schedule_ids:
            coll = self.collisions[sched]
            summary['resolved total collisions'].append(0)
            for method in pc.all_adjustment_methods:
                if method not in summary:
                    summary[method] = []
                if method in coll['resolved']:
                    summary[method].append(len(coll['resolved'][method]))
                else:
                    summary[method].append(0)
                summary['resolved total collisions'][-1] += summary[method][-1]
            summary['found set'].append(coll['found'])
            summary['resolved set'].append(coll['resolved'])
            summary[found_not_registered_str].append(self.found_but_not_resolved(coll['found'],coll['resolved']))
        return summary
    
    def summarize_unresolved_colliding(self):
        """Returns a summary dictionary of the unresolved colliding items."""
        summary = {'unresolved colliding':[], 'unresolved colliding tables':[], 'unresolved colliding sweeps':[], 'final check collision pairs found':[]}
        for sched in self.schedule_ids:
            summary['unresolved colliding'].append(self.unresolved[sched])
            summary['unresolved colliding tables'].append(self.unresolved_tables[sched])
            summary['unresolved colliding sweeps'].append(self.unresolved_sweeps[sched])
            summary['final check collision pairs found'].append(self.final_checked_collision_pairs[sched])
        return summary
        
    def summarize_all(self):
        """Returns a summary dictionary of all data."""
        data = {'schedule id':self.schedule_ids}
        for key in self._strings_to_print_first:
            data.update({key:self.strings[key]})
        data.update(self.numbers)
        data.update(self.summarize_collision_resolutions())
        data.update(self.summarize_unresolved_colliding())
        nrows = len(next(iter(data.values())))
        safe_divide = lambda a,b: a / b if b else np.inf # avoid divide-by-zero errors
        data['calc: fraction of target requests accepted'] = [safe_divide(data['n requests accepted'][i], data['n requests'][i]) for i in range(nrows)]
        data['calc: fraction of targets achieved (of those accepted)'] = [safe_divide(data['n tables achieving requested-and-accepted targets'][i], data['n requests accepted'][i]) for i in range(nrows)]
        for key in self._strings_to_print_last:
            data.update({key:self.strings[key]})
        stripped_data, stripped_nrows = self._copy_and_strip_null_rows(data)
        return stripped_data, stripped_nrows
    
    def _copy_and_strip_null_rows(self, data):
        '''Copies summary data structure and strips out any rows that are
        considered "null". Very low-level implementation-specifc stuff, broken
        out into a separate function here just for readability.'''
        null_rows = [row for row in range(self.n_rows) if not self._row_contains_real_data(row)]
        stripped = copy.deepcopy(data) # to ensure nothing screwy happens to underlying data before deleting any null rows below
        stripped_nrows = self.n_rows
        for row in null_rows:
            stripped_nrows -= 1
            for somelist in stripped.values():
                if 0 <= row < len(somelist):
                    del somelist[row]
        return stripped, stripped_nrows

    def _row_contains_real_data(self, row=-1):
        '''Checks whether a row contains any actual data, or is instead
        completely filled with "null" values. Arguments to input row:
            
            -1 ... last row in data structure
                   (function returns False if data structure is completely empty)
            
            integer >= 0 ... data in that row index
                   (function returns False if that's not a valid index)
        '''
        if row < -1 or row >= self.n_rows or self.n_rows == 0:
            return False
        sched_id = self.schedule_ids[row]
        for dictionary in [self.collisions,
                           self.unresolved,
                           self.unresolved_tables,
                           self.unresolved_sweeps,
                           self.final_checked_collision_pairs]:
            is_empty = self._recursive_is_collection_empty(dictionary[sched_id])
            if not is_empty:
                return True
        for strlist in self.strings.values():
            if strlist[row] != _blank_str:
                return True
        for key, numlist in self.numbers.items():
            value = numlist[row]
            if key == final_checks_str:
                if value != None:
                    return True
            elif key == 'n pos':
                pass # n pos is specifically NOT used as a criterion for real data, due to how it may be pre-set upon initializion
            elif value != 0:
                return True
        return False
    
    def _recursive_is_collection_empty(self, item):
        '''Recursively search a collection and any subcollections it contains,
        returning 1 if any of them contain data, and 0 if not. Here, "data"
        means anything other than an empty collection. And "collection" is only
        guaranteed to work for lists, dicts, and sets.'''
        if isinstance(item, dict) or isinstance(item, list) or isinstance(item, set):
            if len(item) == 0:
                return 1
            else:
                boolean = 1
                for sub in item:
                    if isinstance(item, dict):
                        sub = item[sub]
                    boolean *= self._recursive_is_collection_empty(sub)
                return boolean
        else:
            return 0

    def generate_table(self, footers=False):
        """Returns a pandas dataframe representing a complete report table of
        the current stats.
        
        If the argument footers=True, then a number of rows will be returned
        instead, which match up to the normal table. These special rows will
        give max, min, mean, rms, and median for each data column, collated by
        anticollision method (e.g. 'adjust' vs 'freeze'). That's no new
        information --- just a convenience when for example evaluating sims of
        hundreds of targets in a row.
        """
        data, nrows = self.summarize_all()
        if footers:
            blank_row = {'method':''}
            rms = lambda X: (sum([x**2 for x in X])/len(X))**0.5
            unique_methods = sorted(set(data['method']) - {_blank_str})
            categories = ['overall'] + unique_methods
            calcs = {'max':max, 'min':min, 'rms':rms, 'avg':np.mean, 'med':np.median}
            stats = {}
            for category in categories:
                stats[category] = {}
                for calc in calcs:
                    stats[category][calc] = {'method':calc} # puts the name of this calc in the method column of csv file
                stats[category]['max']['schedule id'] = category # puts the category of this group of calcs in the schedule id column of csv file, in same row as maxes
            for key in data:
                if len(data[key]) > 0:
                    type_test_val = data[key][0]
                    if isinstance(type_test_val,int) or isinstance(type_test_val,float):
                        this_data = {'overall': [data[key][i] for i in range(nrows)]}
                        for i in range(nrows):
                            if this_data['overall'][i] == None:
                                this_data['overall'][i] = 0
                        for category in unique_methods:
                            this_data[category] = [this_data['overall'][i] for i in range(nrows) if data['method'][i] == category]
                        for category in categories:
                            for calc,stat_function in calcs.items():
                                if this_data[category]:
                                    stats[category][calc][key] = stat_function(this_data[category])
        file = io.StringIO(newline='\n')
        writer = csv.DictWriter(file, fieldnames=data.keys())
        writer.writeheader()
        if not footers:
            for i in range(nrows):
                row = {key: val[i] for key, val in data.items()}
                writer.writerow(row)
        else:
            for category in categories:
                writer.writerow(blank_row)
                for calc in calcs:
                    writer.writerow(stats[category][calc])
        file.seek(0)  # go back to the beginning after finishing write
        return pd.read_csv(file)  # returns a pandas dataframe

    def save(self, path=None, footers=False):
        """Saves stats results to disk. If no path was specified, the return
        value is the path that was generated. If path is specified, it should
        have '.csv' file extension. Repeated calls to save() with the same path
        will generally append rows to that csv file.
        
        Boolean argument "footers" instead appends some extra summary statistics
        to the bottom of the output table. These are helpful for  debugging ---
        saves time processing the table --- but may look weird if you are going
        to keep appending new results after them. In other words, this option is
        best applied only when you know you're about to stop appending new data
        to the file at the current path. N.B. Footers will only really work
        right if the clear_cache_after_save property is set to False at the
        beginning of the statistics collection. (Otherwise, the data in memory
        will be flushed after every save to disk.)
        """
        dir_name = os.path.dirname(str(path))
        dir_exists = os.path.isdir(dir_name)
        if not dir_exists:
            parent_dir = os.path.dirname(dir_name)
            parent_dir_exists = os.path.isdir(parent_dir)
            if parent_dir_exists:
                try:
                    os.mkdir(dir_name)
                except FileExistsError: #multiple petals could be trying to do this at the same time
                    pass
                dir_exists = os.path.isdir(dir_name)
        if path is None or not dir_exists:
            suffix = str(self.filename_suffix)
            suffix = '_' + suffix if suffix else ''
            filename = f'{pc.filename_timestamp_str()}_schedstats{suffix}.csv'
            path = os.path.join(pc.dirs['temp_files'], filename)
        include_headers = os.path.exists(path) == False
        frame = self.generate_table(footers=footers)
        if footers:
            save_frame = frame
        else:
            start_row = 0 if self._latest_saved_row == None else self._latest_saved_row + 1
            n_rows_to_save = len(frame) - start_row
            save_frame = frame.tail(n_rows_to_save)
            self._latest_saved_row = len(frame) - 1
        save_frame.to_csv(path, mode='a', header=include_headers, index=False)
        if self.clear_cache_after_save:
            self._init_data_structures()
        return path
    
    @staticmethod
    def found_but_not_resolved(found, resolved):
        """Searches through the dictionary of resolved collision pairs, and
        eliminates these from the set of found collision pairs.
        
        Inputs:
            found ... set of strings
            resolved ... dict with keys = strings, values = sets of strings
            
        Output:
            set of strings
        """
        output = found.copy()
        for f_pair in found:
            for method in resolved:
                for r_pair in resolved[method]:
                    if f_pair == r_pair and f_pair in output: # second check since pair may repeat in multiple resolved[method] sets, for example if collision was fixed first in a retract stage, and then again in an extend stage
                        output.remove(f_pair)
        return output