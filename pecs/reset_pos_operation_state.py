import argparse
from astropy.table import Table
from DOSlib.proxies import PetalMan
from DOSlib.positioner_index import PositionerIndex
import glob
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i','--infile', type=str, help='Input calibration file to determine what to enable/disable. Defaults to latest file.', default='')
parser.add_argument('-dir', '--directory', type=str, help='Directory where calibration files are located. Default /data/focalplane/calibration/', default='/data/focalplane/calibration/')
parser.add_argument('-c', '--comment', required=True, type=str, help='Comment to go with commit')

parser.add_argument('-d', '--ctrl_enable', action='store_true', help='Use file to set ctrl_enable')
parser.add_argument('-n', '--nonfunctional', action='store_true', help='Use file to set device_classified_nonfunctional.')
parser.add_argument('-b', '--brokenfiber', action='store_true', help='Use file to set fiber_intact.')
parser.add_argument('-r', '--retracted', action='store_true', help='Use file to set classified_as_retracted.')

uargs = parser.parse_args()

infile = uargs.infile
if infile == '':
    files = glob.glob(f'{uargs.directory}/*.ecsv')
    infile = max(files, key=os.path.getctime).split('/')[-1]

df = Table.read(uargs.directory + infile)

assert uargs.comment.strip(' ') != '', 'Please provide a comment!'
comment = f'reset_pos_operation_state file={infile}: {uargs.comment}'

if uargs.ctrl_enable | uargs.nonfunctional | uargs.brokenfiber | uargs.retracted:
    ptlm = PetalMan()
    change_flags = ''
    if uargs.ctrl_enable:
        change_flags += 'ctrl_enable '
    if uargs.nonfunctional:
        change_flags += 'device_classified_nonfunctional '
    if uargs.brokenfiber:
        change_flags += 'fiber_intact '
    if uargs.retracted:
        change_flags += 'classified_as_retracted '
    input(f'WARNING: this script will set the following flags ({change_flags}) to match those found in the file ({infile}). Are you sure this is what you want to do? (enter to continue)')



if uargs.ctrl_enable:
    print('Setting ctrl_enable flags...')
    nonfunc = set(df[~df['CTRL_ENABLE']]['POS_ID'])
    func = set(df[df['CTRL_ENABLE']]['POS_ID'])
    ptlm.set_operability(nonfunc, 'CTRL_ENABLE', False, comment=comment)
    ptlm.set_operability(func, 'CTRL_ENABLE', True, comment=comment)

if uargs.nonfunctional:
    print('Setting device_classified_nonfunctional flags...')
    nonfunc = set(df[df['DEVICE_CLASSIFIED_NONFUNCTIONAL']]['POS_ID'])
    func = set(df[~df['DEVICE_CLASSIFIED_NONFUNCTIONAL']]['POS_ID'])
    ptlm.set_operability(nonfunc, 'DEVICE_CLASSIFIED_NONFUNCTIONAL', True, comment=comment)
    ptlm.set_operability(func, 'DEVICE_CLASSIFIED_NONFUNCTIONAL', False, comment=comment)

if uargs.brokenfiber:
    print('Setting fiber_intact flags...')
    intact = set(df[df['FIBER_INTACT']]['POS_ID'])
    broke = set(df[~df['FIBER_INTACT']]['POS_ID'])
    ptlm.set_operability(intact, 'FIBER_INTACT', True, comment=comment)
    ptlm.set_operability(broke, 'FIBER_INTACT', False, comment=comment)

if uargs.retracted:
    print('Setting classified_as_retracted flags...')
    retract = set(df[df['CLASSIFIED_AS_RETRACTED']]['POS_ID'])
    ok = set(df[~df['CLASSIFIED_AS_RETRACTED']]['POS_ID'])
    ptlm.set_operability(retract, 'CLASSIFIED_AS_RETRACTED', True, comment=comment)
    ptlm.set_operability(ok, 'CLASSIFIED_AS_RETRACTED', False, comment=comment)

print('Done! Please double check posDB.')
