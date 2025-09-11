"""Creates a set of simulated move requests for testing the anticollision code.
This particular script makes commands that exercise all of the different move
command options:
    
    ABSOLUTE       RELATIVE
    'QS'           'dQdS'
    'obsXY'        'obsdXdY'
    'poslocXY'     'poslocdXdY'
    'ptlXY'
    'posintTP'     'dTdP'
    'poslocTP'
    
For ease of verification, and to minimize opportunities for odd calibration
or anticollision conditions to interfere, in all these coordinate systems the
positioners are simply made to move back and forth between two easy-to-reach
keypoints.

This means that this sequence is not a complete test inspecting all the corners
of the move scheduling system. Rather it is focused on demonstrating proper
function just of these various coordinate systems.

To get some variability into the system, the simulated calibration parameters
are randomly selected.
"""



import os
import sys
sys.path.append(os.path.abspath('../../../petal'))
import csv
from read_nominal_pos_locations import locations
import poscollider
import posmodel
import posstate
import harness_constants as hc
import numpy

# save options
reqests_first_file_num = 90000
params_file_num = 90000
should_save = True

# key points
keypoints = [[0,150],
             [30,120]]

# uv_types
uv_types = ['QS', 'dQdS', 'obsXY', 'obsdXdY', 'poslocXY', 'poslocdXdY',
            'ptlXY', 'posintTP', 'dTdP', 'poslocTP']

# input parameters
params = {}
params['LENGTH_R1'] = {'nom':3.0, 'min':2.7, 'max':3.3}
params['LENGTH_R2'] = {'nom':3.0, 'min':2.7, 'max':3.3}
params['OFFSET_T'] = {'nom':0.0, 'min':-179.0, 'max':179.0}
params['OFFSET_P'] = {'nom':0.0, 'min':-5.0, 'max':15.0}

# initialize a nominal poscollider instance (for reasonable target checking)
collider = poscollider.PosCollider()
posmodels = {}
for posid,data in locations.items():
    state = posstate.PosState(posid)
    state.store('DEVICE_LOC',data['DEVICE_LOC'])
    state.store('OFFSET_X',data['nomX'])
    state.store('OFFSET_Y',data['nomY'])
    state.store('PHYSICAL_RANGE_T',380.0)
    state.store('PHYSICAL_RANGE_P',200.0)
    for key,defs in params.items():
        random_val = numpy.random.uniform(low=defs['min'], high=defs['max'])
        state.store(key, random_val)
    model = posmodel.PosModel(state)
    posmodels[posid] = model
collider.add_positioners(posmodels.values())

# generate random targets
all_requests = []
start_posintTP = {posid:None for posid in posmodels}
final_posintTP = {posid:None for posid in posmodels}
for cmd in uv_types:
    for kp in keypoints:
        these_requests = []
        for posid,model in posmodels.items():
            trans = model.trans
            start_posintTP[posid] = final_posintTP[posid] # therefore first target in sequence should be an absolute target, not relative, so that something will exist here to assign
            final_posintTP[posid] = trans.poslocTP_to_posintTP(kp)
            if cmd == 'QS':
                new_uv = trans.posintTP_to_QS(final_posintTP[posid])
            elif cmd == 'dQdS':
                start_qs = trans.posintTP_to_QS(start_posintTP[posid])
                final_qs = trans.posintTP_to_QS(final_posintTP[posid])
                new_uv = trans.delta_QS(final_qs, start_qs)
            elif cmd == 'obsXY':
                new_uv = trans.posintTP_to_obsXY(final_posintTP[posid])
            elif cmd == 'obsdXdY':
                start_xy = trans.posintTP_to_obsXY(start_posintTP[posid])
                final_xy = trans.posintTP_to_obsXY(final_posintTP[posid])
                new_uv = trans.delta_XY(final_xy, start_xy)                
            elif cmd == 'poslocXY':
                new_uv = trans.posintTP_to_poslocXY(final_posintTP[posid])
            elif cmd == 'poslocdXdY':
                start_xy = trans.posintTP_to_poslocXY(start_posintTP[posid])
                final_xy = trans.posintTP_to_poslocXY(final_posintTP[posid])
                new_uv = trans.delta_XY(final_xy, start_xy)      
            elif cmd == 'ptlXY':
                new_uv = trans.posintTP_to_ptlXY(final_posintTP[posid])
            elif cmd ==  'posintTP':
                new_uv = final_posintTP[posid]
            elif cmd == 'dTdP':
                new_uv = trans.delta_posintTP(final_posintTP[posid], start_posintTP[posid])
            elif cmd == 'poslocTP':
                new_uv = trans.posintTP_to_poslocTP(final_posintTP[posid])
            else:
                print(f'Error {cmd} not matched to a defined transformation.')
            new_request = {'DEVICE_LOC': model.state._val['DEVICE_LOC'],
                           'command': cmd,
                           'u': new_uv[0],
                           'v': new_uv[1]}
            these_requests.append(new_request)
        all_requests.append(these_requests)
    
# save request files
if should_save:
    next_filenumber = reqests_first_file_num
    for requests in all_requests:
        save_path = hc.filepath(hc.req_dir, hc.req_prefix, next_filenumber)
        with open(save_path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=requests[0].keys())
            writer.writeheader()
            for request in requests:
                writer.writerow(request)
        next_filenumber += 1
        
# save pos param files
if should_save:
    save_path = hc.filepath(hc.pos_dir, hc.pos_prefix, params_file_num)
    with open(save_path, 'w', newline='') as csvfile:
        headers = ['POS_ID', 'DEVICE_LOC', 'CTRL_ENABLED', 'LENGTH_R1', 'LENGTH_R2',
                   'OFFSET_T', 'OFFSET_P', 'OFFSET_X', 'OFFSET_Y', 'PHYSICAL_RANGE_T',
                   'PHYSICAL_RANGE_P']
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        for model in posmodels.values():
            state = model.state
            data = {key:state._val[key] for key in headers}
            writer.writerow(data)
