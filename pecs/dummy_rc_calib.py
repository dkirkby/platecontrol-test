from astropy.table import Table
from pecs import PECS
import pandas as pd
import traceback

rc_calib = './seq_RC_CALIB.ecsv'
tab = Table.read(rc_calib, fill_values=None)

cs = PECS(interactive=True)
posinfo = cs.ptlm.get_positioners(posids=cs.posids, enabled_only=True)
posinfo = pd.concat([p for p in posinfo.values()])
base_req = posinfo.drop(columns='BUS_ID')

def make_req(base, x1, x2, command, log=''):
    req = base.copy()
    req['X1'] = x1
    req['X2'] = x2
    req['COMMAND'] = command
    req['LOG_NOTE'] = log
    return req

limit_angles = cs.ptlm.get_phi_limit_angle()
cs.ptlm.set_phi_limit_angle(None)

tablength = len(tab)
i = 0
try:
    for row in tab:
        command = row['command']
        x1 = row['target0']
        x2 = row['target1']
        log = row['log_note']
        print(f'Move {i+1} of {tablength}')
        print(f'Moving to ({x1}, {x2}) with command {command}.')
        print(log)
        if command not in ['posintTP', 'poslocXY']:
            #pass
            cs.rehome_and_measure(cs.posids, log_note=log)
        else:
            req = make_req(base_req, x1, x2, command, log=log)
            cs.move_measure(request=req)
        i += 1
except:
    traceback.print_exc()

cs.fvc_collect()
for ptl, value in limit_angles.items():
    cs.ptlm.set_phi_limit_angle(value, participating_petals=ptl)
