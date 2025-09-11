"""Creates a set of simulated positioners for testing the anticollision code.
"""

import csv
import random
import os
import harness_constants as hc

# input paramters
percent_disabled = 0.05 # fraction
armlength_max_err = 0.2 # mm
offset_phi_max_err = 10.0 # deg
offset_xy_max_err = 0.1 # mm
range_max_err = 10.0 # deg

# gather nominal positioner locations from file
from read_nominal_pos_locations import locations

# gather / randomize positioner config file parameters
conf_data = {posid:{} for posid in locations}
for posid,loc_data in locations.items():
    data = conf_data[posid]
    data['POS_ID']           = loc_data['POS_ID']
    data['DEVICE_LOC']        = loc_data['DEVICE_LOC']
    data['CTRL_ENABLED']     = random.choices([True,False],[1-percent_disabled,percent_disabled])[0]
    data['LENGTH_R1']        = 3.0 + random.uniform(-armlength_max_err,armlength_max_err)
    data['LENGTH_R2']        = 3.0 + random.uniform(-armlength_max_err,armlength_max_err)
    data['OFFSET_T']         = random.uniform(-180,180)
    data['OFFSET_P']         = random.uniform(-offset_phi_max_err,offset_phi_max_err)
    data['OFFSET_X']         = locations[posid]['nomX'] + random.uniform(-offset_xy_max_err,offset_xy_max_err)
    data['OFFSET_Y']         = locations[posid]['nomY'] + random.uniform(-offset_xy_max_err,offset_xy_max_err)
    data['PHYSICAL_RANGE_T'] = 380.0 + random.uniform(-range_max_err,range_max_err)
    data['PHYSICAL_RANGE_P'] = 200.0 + random.uniform(-range_max_err,range_max_err)
    
# save set file
existing = os.listdir(hc.pos_dir)
existing_filenumbers = [hc.filenumber(name,hc.pos_prefix) for name in existing]
new_filenumber = max(existing_filenumbers) + 1 if existing_filenumbers else 0
save_path = hc.filepath(hc.pos_dir, hc.pos_prefix, new_filenumber)
with open(save_path,'w',newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=conf_data[next(iter(conf_data))])
    writer.writeheader()
    for data in conf_data.values():
        writer.writerow(data)