'''
Sets offsetsXY to nominal values. Needs running DOS instance

Also called seed xy offsets, this script establishes initial values for
positioner offsets using nominal specs or
metrology data, right after petal assemblies are mounted to the ring
on-mountain, so that calibration can be done.

Calibration, or particularly spotmatch, relies on knowing a priori
approximately where things are, and searching within a search radius.
After calibration, the estimatd parameters will be updated with newly
measured values.

The x, y offsets of positioenrs for the central axis can be estimated
from the designed focal plate layout. So we'll use the theoretical
positioner x, y offsets as a starting point.

In the new local coordinate transformation scheme, xy offsets are
in the petal's local CS, from petal's nominal origin to positioner
centre (of patrol area).
We are only doing this for positioners, no fiducials.

Returns two pandas Dataframes with columns
DEVICE_ID, PETAL_LOC, DEVICE_LOC, old and new calibration vals, with
MODE for updates and old values.
'''
import os
import pandas as pd
import posconstants as pc
from pecs import PECS

seed = PECS(interactive=True, no_expid=True)
print('Seeding offsets XY...')
nominals = pd.read_csv(pc.positioner_locations_file)
nominals.set_index('device_location_id', inplace=True)
updates = []
for posid, row in seed.posinfo.iterrows():
    update = {'DEVICE_ID': posid, 'MODE': 'seed_offsets_xy'}
    device_loc, petal_loc = row['DEVICE_LOC'], row['PETAL_LOC']
    if max(seed.pcids) > 899:
        role = seed._pcid2role(900+petal_loc) #sim mode
    else:
        role = seed._pcid2role(petal_loc)
    x = nominals['FLAT_X'][device_loc]
    y = nominals['FLAT_Y'][device_loc]
    update = seed.ptlm.collect_calib(update, tag='OLD_',
                                     participating_petals=role)[role]
    seed.ptlm.set_posfid_val(posid, 'OFFSET_X', x, participating_petals=role)
    seed.ptlm.set_posfid_val(posid, 'OFFSET_Y', y, participating_petals=role)
    update = seed.ptlm.collect_calib(update, tag='',
                                     participating_petals=role)[role]
    updates.append(update)
seed.ptlm.commit(mode='calib', calib_note='seed_offsets_xy')
updates = pd.DataFrame(updates).set_index('DEVICE_ID').sort_index()
path = os.path.join(pc.dirs['calib_logs'],
                    f'{pc.filename_timestamp_str()}-seed_offsets_xy.csv')
updates.to_csv(path)
# preview calibration updates
print(updates[['POS_T', 'POS_P', 'LENGTH_R1', 'LENGTH_R2',
               'OFFSET_X', 'OFFSET_Y', 'OFFSET_T', 'OFFSET_P']])
print(f'Seed offsets XY data saved to: {path}')

