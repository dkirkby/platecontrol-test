# -*- coding: utf-8 -*-
"""
Created on Fri May 31 18:09:54 2019

@author: Duan Yutong
"""


from DOSlib.positioner_index import PositionerIndex
PETAL_ID = 900
tablepath = '/software/products/PositionerIndexTable-trunk/index_files/mock_observing_20180619.csv'
p = PositionerIndex(tablepath)
ret = p.find_by_arbitrary_keys(PETAL_ID=PETAL_ID, DEVICE_TYPE='POS')
posids = [item['DEVICE_ID'] for item in ret]
ret = p.find_by_arbitrary_keys(PETAL_ID=PETAL_ID, DEVICE_TYPE='FIF')
fidids = [item['DEVICE_ID'] for item in ret]
#print('posids[0]', posids[0])
posidpath = f'/home/msdos/focalplane/fp_settings/ptl_settings/unit_{PETAL_ID}_posids.txt'
fididpath = f'/home/msdos/focalplane/fp_settings/ptl_settings/unit_{PETAL_ID}_fidids.txt'
with open(posidpath, 'w') as h:
    h.write(', '.join(posids))
with open(fididpath, 'w') as h:
    h.write(', '.join(fidids))
print(f'PTL {PETAL_ID}, localtion {ret[0]["PETAL_LOC"]}',
      f'{len(posids)} posids and {len(fidids)} fidids written')
