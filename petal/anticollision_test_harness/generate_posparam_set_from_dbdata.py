# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 16:01:48 2020
@author: jhsilber

Customized script for merging data into posparam format. Data is presumed
pulled from online databases:
    
    "calib" --> http://web.replicator.dev-cattle.stable.spin.nersc.org:60040/DESIPositioners/Positioner_calibration.html
    "moves" --> http://web.replicator.dev-cattle.stable.spin.nersc.org:60040/DESIPositioners/Positioner_moves.html
    
Alternate DB server, for petals 0 and 1 at LBNL: http://beyonce.lbl.gov:8080

Data is assumed saved as csv using the browser tool.

Sample SQL queries:
    
    "calib" --> select petal_id,pos_id,time_recorded,device_loc,length_r1,length_r2,offset_x,offset_y,offset_t,offset_p,physical_range_t,physical_range_p,gear_calib_t,gear_calib_p,keepout_expansion_phi_radial,keepout_expansion_phi_angular,keepout_expansion_theta_radial,keepout_expansion_theta_angular,classified_as_retracted from posmovedb.positioner_calibration_p3 PM where EXISTS (select * from posmovedb.positioner_calibration_index_p3 where id = PM.pos_id and last_index=PM.pos_calib_index)
    "moves" --> select petal_id,pos_id,time_recorded,ctrl_enabled from posmovedb.positioner_moves_p3 PM where exists (select * from posmovedb.positioner_moves_index_p3 where id = PM.pos_id and last_index=PM.pos_move_index)

Take care when copy/pasting these queries into the web tool that you replace the "_p3"
suffixes with ones that identify the petal you are interested in. Like "_p1" for
petal ID 1, etcetera. And note there are two instances of "_p3" in each query line.

"""

import os
import pandas

posids = 'all' # 'all', for all posids on that petal, or else a limited set
output_file = 'posparams_00001.csv'

data_folder = 'C:/Users/joe/Desktop/DESI-5834-v1/online_calib_data_retrieved_2020-09-18'
calibdb_name = 'PTL00 online calib db, retrieved 2020-09-18.csv'
movesdb_name = 'PTL00 online moves db, retrieved 2020-09-18.csv'
calibpath = os.path.join(data_folder, calibdb_name)
movespath = os.path.join(data_folder, movesdb_name)
calib = pandas.read_csv(calibpath)
moves = pandas.read_csv(movespath)

if posids == 'all':
    posids = set(calib['pos_id']).intersection(moves['pos_id'])
posids = sorted(posids)

calib.sort_values('time_recorded', ascending=False, inplace=True)
moves.sort_values('time_recorded', ascending=False, inplace=True)

columns_calib = ['DEVICE_LOC', 'LENGTH_R1', 'LENGTH_R2', 'OFFSET_T',
                 'OFFSET_P', 'OFFSET_X', 'OFFSET_Y', 'PHYSICAL_RANGE_T',
                 'PHYSICAL_RANGE_P', 'KEEPOUT_EXPANSION_PHI_RADIAL',
                 'KEEPOUT_EXPANSION_PHI_ANGULAR', 'KEEPOUT_EXPANSION_THETA_RADIAL',
                 'KEEPOUT_EXPANSION_THETA_ANGULAR', 'CLASSIFIED_AS_RETRACTED']
columns_moves = ['CTRL_ENABLED']
columns_calib_upperlower = {C:C.lower() for C in columns_calib}
columns_moves_upperlower = {C:C.lower() for C in columns_moves}

result = pandas.DataFrame(index=posids, columns=columns_calib + columns_moves)

def grab(dataframe, columns_upperlower):
    remaining = set(posids)
    for i in dataframe.index:
        posid = dataframe['pos_id'][i]
        if posid in remaining:
            for C,c in columns_upperlower.items():
                result[C][posid] = dataframe[c][i]
            remaining.remove(posid)
            
grab(calib, columns_calib_upperlower)
grab(moves, columns_moves_upperlower)

result.to_csv(output_file, index_label='POS_ID')
      