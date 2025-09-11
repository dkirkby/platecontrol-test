'''
Set armlengths R1 and R2 to nominal values. Needs running DOS instance
'''
import os
import pandas as pd
from pecs import PECS
import posconstants as pc

seed = PECS(interactive=True, no_expid=True)
print(f'Seeding armlengths for PCIDs: {seed.pcids}')
updates = []
for posid, row in seed.posinfo.iterrows():
    update = {'DEVICE_ID': posid, 'MODE': 'seed_offsets_tp'}
    petal_loc = row['PETAL_LOC']
    if max(seed.pcids) > 899:
        role = seed._pcid2role(900+petal_loc)
    else:
        role = seed._pcid2role(petal_loc)
    update = seed.ptlm.collect_calib(update, tag='OLD_',
                                     participating_petals=role)[role]
    seed.ptlm.set_posfid_val(posid, 'LENGTH_R1',
                             pc.nominals['LENGTH_R1']['value'],
                             participating_petals=role)
    seed.ptlm.set_posfid_val(posid, 'LENGTH_R2',
                             pc.nominals['LENGTH_R2']['value'],
                             participating_petals=role)
    update = seed.ptlm.collect_calib(update, tag='',
                                     participating_petals=role)[role]
    updates.append(update)
seed.ptlm.commit(mode='calib', calib_note='seed_armlengths')
updates = pd.DataFrame(updates).set_index('DEVICE_ID').sort_index()
path = os.path.join(pc.dirs['calib_logs'],
                    f'{pc.filename_timestamp_str()}-seed_armlengths.csv')
updates.to_csv(path)
# preview calibration updates
print(updates[['POS_T', 'POS_P', 'LENGTH_R1', 'LENGTH_R2',
               'OFFSET_X', 'OFFSET_Y', 'OFFSET_T', 'OFFSET_P']])
print(f'Seed offsets TP data saved to: {path}')
