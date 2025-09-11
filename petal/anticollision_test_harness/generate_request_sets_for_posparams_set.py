"""Creates a set of simulated move requests for testing the anticollision code.
"""

import os
import sys
sys.path.append(os.path.abspath('../../petal'))
import csv
import random
import poscollider
import posmodel
import posstate
import harness_constants as hc
import sequences

# input parameters
petal_id = int(input('Enter petal id (max 2 digits): '))
assert 0 <= petal_id < 100
params_id = int(input('Enter params id (max 5 digits): '))
assert 0 <= params_id < 100000
set_number = int(input('Enter set number (max 2 digits): '))
assert 0 <= set_number < 100
num_req_to_make = int(input('Enter number of requests to make (max 5 digits): '))
assert 0 <= num_req_to_make < 100000
params_prefix = hc.make_params_prefix(petal_id)
posparams = sequences._read_data(data_id=params_id, prefix=params_prefix)
posids =  'all' # 'all' for all posids in posparams, otherwise a set of selected ones

# initialize a poscollider instance (for fair target checking)
collider = poscollider.PosCollider()
posmodels = {}
keys_to_copy = ['POS_ID','DEVICE_LOC','CTRL_ENABLED',
                'LENGTH_R1','LENGTH_R2',
                'OFFSET_T','OFFSET_P','OFFSET_X','OFFSET_Y',
                'PHYSICAL_RANGE_T','PHYSICAL_RANGE_P',
                'KEEPOUT_EXPANSION_PHI_RADIAL', 'KEEPOUT_EXPANSION_PHI_ANGULAR',
                'KEEPOUT_EXPANSION_THETA_RADIAL', 'KEEPOUT_EXPANSION_THETA_ANGULAR',
                'CLASSIFIED_AS_RETRACTED']
for posid,data in posparams.items():
    if posid in posids or posids == 'all':
        state = posstate.PosState(posid)
        for key in keys_to_copy:
            state.store(key,data[key], register_if_altered=False)  
        model = posmodel.PosModel(state)
        posmodels[posid] = model
collider.add_positioners(posmodels.values())

# generate random targets
all_targets = []
for i in range(num_req_to_make):
    targets_obsTP = {}
    targets_posXY = {}
    for posid,model in posmodels.items():
        attempts_remaining = 10000 # just to prevent infinite loop if there's a bug somewhere
        while posid not in targets_obsTP and attempts_remaining:
            rangeT = model.targetable_range_posintT
            rangeP = model.targetable_range_posintP
            max_patrol = posparams[posid]['LENGTH_R1'] + posparams[posid]['LENGTH_R2']
            min_patrol = abs(posparams[posid]['LENGTH_R1'] - posparams[posid]['LENGTH_R2'])
            x = random.uniform(-max_patrol,max_patrol)
            y = random.uniform(-max_patrol,max_patrol)
            this_radius = (x**2 + y**2)**0.5
            if min_patrol > this_radius or this_radius > max_patrol:
                attempts_remaining = max(attempts_remaining - 1, 0)
            else:
                this_posXY = [x,y]
                this_posTP = model.trans.poslocXY_to_posintTP(this_posXY)[0]
                this_obsTP = model.trans.posintTP_to_poslocTP(this_posTP)
                target_interference = False
                for neighbor in collider.pos_neighbors[posid]:
                    if neighbor in targets_obsTP:
                        if collider.spatial_collision_between_positioners(posid, neighbor, this_obsTP, targets_obsTP[neighbor]):
                            target_interference = True
                            break
                out_of_bounds = collider.spatial_collision_with_fixed(posid, this_obsTP)
                if not target_interference and not out_of_bounds:
                    targets_obsTP[posid] = this_obsTP
                    targets_posXY[posid] = this_posXY
                else:
                    attempts_remaining = max(attempts_remaining - 1, 0)
        if not attempts_remaining:
            v = model.state._val
            print(f'Warning: no valid target found for posid: {posid} at location {v["DEVICE_LOC"]}' +
                  f' (x0, y0) = ({v["OFFSET_X"]:.3f}, {v["OFFSET_Y"]:.3f}). Target set {i} will not' +
                  ' include an entry for this positioner.')
    all_targets.append(targets_posXY)

# gather device locations and sort rows
device_loc_list = sorted(collider.devicelocs.keys())
all_targets_by_loc = []
for targets in all_targets:
    targets_by_loc = {}
    for loc in device_loc_list:
        posid = collider.devicelocs[loc]
        if posid in targets:
            targets_by_loc[loc] = targets[posid]
    all_targets_by_loc.append(targets_by_loc)

# save set files
i = 0
prefix = hc.make_request_prefix(petal_id, set_number)
overwrite = False
existing = os.listdir(hc.req_dir)
for target in all_targets_by_loc:
    save_path = hc.filepath(hc.req_dir, prefix, i)
    name = os.path.basename(save_path)
    if name in existing:
        if not overwrite:
            yesno = input(f'Some files (e.g. {name}) already exist. Overwrite them? (y/n): ')
            overwrite = yesno.lower() in {'true', 'y', 't', 'yes', '1'}
        ok_to_write = overwrite
    else:
        ok_to_write = True
    if ok_to_write:
        with open(save_path,'w',newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['DEVICE_LOC','command','u','v'])
            for loc,uv in target.items():
                writer.writerow([loc,'poslocXY',uv[0],uv[1]])
    i += 1