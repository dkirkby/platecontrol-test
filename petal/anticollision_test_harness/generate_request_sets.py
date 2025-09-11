"""Creates a set of simulated move requests for testing the anticollision code.
"""

import os
import sys
sys.path.append(os.path.abspath('../../../petal'))
import csv
import random
from read_nominal_pos_locations import locations
import poscollider
import posmodel
import posstate
import harness_constants as hc

# input parameters
num_sets_to_make = 75
nom_R1 = 3.0 # mm, nominal LENGTH_R1
nom_R2 = 3.0 # mm, nominal LENGTH_R2
min_patrol = abs(nom_R1 - nom_R2)
max_patrol = nom_R1 + nom_R2

# initialize a nominal poscollider instance (for reasonable target checking)
collider = poscollider.PosCollider()
posmodels = {}
for posid,data in locations.items():
    state = posstate.PosState(posid)
    state.store('DEVICE_LOC',data['DEVICE_LOC'])
    state.store('LENGTH_R1',nom_R1)
    state.store('LENGTH_R2',nom_R2)
    state.store('OFFSET_T',0.0)
    state.store('OFFSET_P',0.0)
    state.store('OFFSET_X',data['nomX'])
    state.store('OFFSET_Y',data['nomY'])
    state.store('PHYSICAL_RANGE_T',380.0)
    state.store('PHYSICAL_RANGE_P',200.0)    
    model = posmodel.PosModel(state)
    posmodels[posid] = model
collider.add_positioners(posmodels.values())

# generate random targets
all_targets = []
for i in range(num_sets_to_make):
    targets_obsTP = {}
    targets_posXY = {}
    for posid,model in posmodels.items():
        attempts_remaining = 100 # just to prevent infinite loop if there's a bug somewhere
        while posid not in targets_obsTP:
            rangeT = model.targetable_range_posintT
            rangeP = model.targetable_range_posintP
            while attempts_remaining:
                x = random.uniform(-max_patrol,max_patrol)
                y = random.uniform(-max_patrol,max_patrol)
                this_radius = (x**2 + y**2)**0.5
                if min_patrol <= this_radius and this_radius <= max_patrol:
                    break
                attempts_remaining -= 1
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
                attempts_remaining -= 1
        if not attempts_remaining:
            v = state._val
            print('Warning: no valid target found for posid: ' + posid + ' at location ' + str(v['DEVICE_LOC']) + ' (x,y) = (' + format(v['OFFSET_X'],'.3f') + ',' + format(v['OFFSET_Y'],'.3f') + ')')
    all_targets.append(targets_posXY)
    
# save set files
existing = os.listdir(hc.req_dir)
existing_filenumbers = [hc.filenumber(name,hc.req_prefix) for name in existing]
next_filenumber = max(existing_filenumbers) + 1 if existing_filenumbers else 0
for target in all_targets:
    save_path = hc.filepath(hc.req_dir, hc.req_prefix, next_filenumber)
    with open(save_path,'w',newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['DEVICE_LOC','command','u','v'])
        for posid,uv in target.items():
            writer.writerow([posmodels[posid].state._val['DEVICE_LOC'],'posXY',uv[0],uv[1]])
    next_filenumber += 1