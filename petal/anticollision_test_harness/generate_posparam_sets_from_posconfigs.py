"""Reads in config files from pos_settings folder, and generates sets of
corresponding simulated positioners, for testing the anticollision code.
"""

import os
import configobj
import harness_constants as hc
import csv

petals_to_generate = [2,3,4,5,6,7,8,9,10,11] # petal IDs
pos_settings_dir_rel = '../../../../../fp_settings/pos_settings/'
pos_settings_dir_abs = os.path.abspath(pos_settings_dir_rel)
all_files = os.listdir(pos_settings_dir_abs)
m_files = [file for file in all_files if 'unit_M' in file]
data = {petal_id:{} for petal_id in petals_to_generate}
keys_to_copy = ['POS_ID','DEVICE_LOC','CTRL_ENABLED',
                'LENGTH_R1','LENGTH_R2',
                'OFFSET_T','OFFSET_P','OFFSET_X','OFFSET_Y',
                'PHYSICAL_RANGE_T','PHYSICAL_RANGE_P']
for file in m_files:
    filepath = os.path.join(pos_settings_dir_abs,file)
    conf = configobj.ConfigObj(filepath,unrepr=True,encoding='utf-8')
    petal_id = conf['PETAL_ID']
    pos_id = conf['POS_ID']
    if petal_id in data:
        data[petal_id][pos_id] = {}
        for key in keys_to_copy:
            data[petal_id][pos_id][key] = conf[key] 
    
# save set files
for petal_id,this_data in data.items():
    filenumber = petal_id + 10000 # make filenames be like 'posparams_100XX.csv', where XX is the petal_id, and compatible with existing anticollision harness code
    save_path = hc.filepath(hc.pos_dir, hc.pos_prefix, filenumber)
    with open(save_path,'w',newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=this_data[next(iter(this_data))])
        writer.writeheader()
        for values in this_data.values():
            writer.writerow(values)