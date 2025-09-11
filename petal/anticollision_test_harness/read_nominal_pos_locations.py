"""Reads in the data for nominal device_location_id.
"""

import csv
import harness_constants as hc

# gather nominal positioner locations from file
locations = {}
id_prefix = 'S'
with open(hc.device_locations_path,'r',newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if row['device_type'] in {'POS','OPT'}:
            posid = id_prefix + format(int(row['device_location_id']),'05d')
            locations[posid] = {'POS_ID':posid,
                                'DEVICE_LOC':int(row['device_location_id']),
                                'nomX':float(row['X']),
                                'nomY':float(row['Y'])}