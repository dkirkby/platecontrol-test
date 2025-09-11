# -*- coding: utf-8 -*-
"""
Reads in a schedule stats file representing simulation results from a sequence
of lots of targets. Reduces this to smaller sequences, with preference for
those that are more "stressful" on the algorithm.

Created on Mon Mar  9 15:30:16 2020
@author: jhsilber
"""

import os
import sys
import random
import pandas as pd
sys.path.append(os.path.abspath('../../petal/'))
import posconstants as pc


# load stats data
statsfile_dir = 'C:\\Users\\joe\\fp_temp_files'
statsfile_name = 'sched_stats_075RQO_01000-01999.csv'
statsfile_path = os.path.join(statsfile_dir, statsfile_name)
stats = pd.read_csv(statsfile_path)
sequence_prefix = '01002'

# clear out footers
not_string = [i for i in range(len(stats)) if not isinstance(stats['method'].iloc[i], str)]
footer_start = not_string[0]
footer_rows = [i for i in range(footer_start, len(stats))]
stats.drop(index=footer_rows, inplace=True)

# extract move request ids
requests_key = 'MOVE_REQUESTS_ID'
requests_ids = [s.split(': ')[-1] for s in stats['note']]
stats[requests_key] = requests_ids

# desired shape of output target request sequence
output_shapes = [{'n_targets': 1,  'n_sets':1},
                 {'n_targets': 5,  'n_sets':2},
                 {'n_targets': 10, 'n_sets':3},
                 {'n_targets': 50, 'n_sets':1}]

# weighting towards "stressful" cases
found_total_key = 'found total collisions (before anticollision)'
stress_indicators = {}
stress_indicators[0] = (2 * stats[found_total_key]).to_list()
stress_indicators[1] = stats['num path adjustment iters'].to_list()
stress_indicators[2] = (3 * (stats['n move tables'] - stats['n requests accepted'])).to_list()
stress_indicators[3] = (stats['n requests'] - stats['n tables achieving requested-and-accepted targets']).to_list()
weights = [0]*len(stress_indicators[0])
for s in stress_indicators.values():
    weights = [weights[i] + s[i] for i in range(len(weights))]
weight_scaler = 10.0 # another enhancement of stressful cases probability
weights = [w * weight_scaler for w in weights]
min_weight = 1.0
weights = [w + min_weight for w in weights]

# generate sequences
sequences = {}
num_collisions_to_avoid = {}
for shape in output_shapes:
    n_targets = shape['n_targets']
    n_sets = shape['n_sets']
    for s in range(n_sets):
        sequence_name = sequence_prefix + '_ntarg' + format(n_targets,'03d') + '_set' + format(s,'03d')
        sequence = []
        n_missing = n_targets
        while n_missing: # while loop here ensures no duplicate entries   
            new = random.choices(population = stats[requests_key].to_list(),
                                 weights = weights,
                                 k = n_missing)
            sequence.extend(new)
            n_missing = n_targets - len(set(sequence)) # set conversion ensures uniqueness
        num_coll = [int(stats[found_total_key][stats[requests_key] == request]) for request in sequence]
        sequences[sequence_name] = [int(request) for request in sequence]
        num_collisions_to_avoid[sequence_name] = num_coll

# generate sequences and print to an output text file
output_name = 'stress_test_sequences.txt'
output_path = os.path.join(statsfile_dir, output_name)
with open(output_path, 'w', newline='') as file:
    writeline = lambda string: file.write(string + '\n')
    writeline('Anticollision target sequences for "stress-tests"')
    writeline('Date generated: ' + pc.timestamp_str())
    writeline('Input stats file: ' + statsfile_name)
    writeline('')
    writeline('Generated sequences:')
    file.write(str(sequences))
    writeline('')
    writeline('')
    writeline('Number of collisions to avoid in each sequence:')
    file.write(str(num_collisions_to_avoid))
