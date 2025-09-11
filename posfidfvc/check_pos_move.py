import numpy as np
from astropy.io import fits
from lmfit import minimize, Parameters
import os
import sys
sys.path.append(os.path.abspath('../petal/'))
sys.path.append(os.path.abspath('../posfidfvc/'))

import petal
import posmovemeasure
import fvchandler
import posconstants as pc
import instrmaker
import tkinter
import tkinter.filedialog
import tkinter.messagebox
import configobj
import csv
import getpass

    
    
    """This code load a petal and fvc (SBIG), take a image, move one positioner, take another image, 
       and subtract the two images. This is to check if the positioner really moves. 
    """


# start set of new and changed files
new_and_changed_files = set()

# unique timestamp and fire up the gui
start_filename_timestamp = pc.filename_timestamp_str()
gui_root = tkinter.Tk()

# get the station config info
message = 'Pick hardware setup file.'
hwsetup_conf = tkinter.filedialog.askopenfilename(initialdir=pc.dirs['hwsetups'], filetypes=(("Config file","*.conf"),("All Files","*")), title=message)
hwsetup = configobj.ConfigObj(hwsetup_conf,unrepr=True)
new_and_changed_files.add(hwsetup.filename)

# are we in simulation mode?
sim = hwsetup['fvc_type'] == 'simulator'
# update log and settings files from the SVN
if not(sim):
    print("")
    svn_user=input("Please enter your SVN username: ")
    svn_pass=getpass.getpass("Please enter your SVN password: ")
    svn_auth_err=False
    svn_update_dirs = [pc.dirs[key] for key in ['pos_logs','pos_settings','xytest_logs','xytest_summaries']]
    should_update_from_svn = tkinter.messagebox.askyesno(title='Update from SVN?',message='Overwrite any existing local positioner log and settings files to match what is currently in the SVN?')
    if should_update_from_svn:
        if svn_auth_err:
            print('Could not validate svn user/password.')
        else:
            for d in svn_update_dirs:
                os.system('svn update --username ' + svn_user + ' --password ' + svn_pass + ' --non-interactive ' + d)
                os.system('svn revert --username ' + svn_user + ' --password ' + svn_pass + ' --non-interactive ' + d + '*')

# software initialization and startup
fvc = fvchandler.FVCHandler(fvc_type=hwsetup['fvc_type'],save_sbig_fits=hwsetup['save_sbig_fits'])
fvc.rotation = hwsetup['rotation'] # this value is used in setups without fvcproxy / platemaker
fvc.scale = hwsetup['scale'] # this value is used in setups without fvcproxy / platemaker
posids = hwsetup['pos_ids']
fidids = hwsetup['fid_ids']
ptl = petal.Petal(hwsetup['ptl_id'], posids, fidids, simulator_on=sim, user_interactions_enabled=True, anticollision=None)
m = posmovemeasure.PosMoveMeasure([ptl],fvc)
m.make_plots_during_calib = True
print('Automatic generation of calibration plots is turned ' + ('ON' if m.make_plots_during_calib else 'OFF') + '.')
for ptl in m.petals:
    for posmodel in ptl.posmodels:
        new_and_changed_files.add(posmodel.state.conf.filename)
        new_and_changed_files.add(posmodel.state.log_path)
    for fidstate in ptl.fidstates.values():
        new_and_changed_files.add(fidstate.conf.filename)
        new_and_changed_files.add(fidstate.log_path)
m.n_extradots_expected = hwsetup['num_extra_dots']


m.request_home()

# Read dots identification result from ptl and store a dictionary
posids=ptl.posids
fidids=ptl.fidids
n_pos=len(posids)
n_fid=len(fidids)
pos_fid_dots={}
obsX_arr,obsY_arr=[],[]
metro_X_arr,metro_Y_arr=[],[]

# read the Metrology data firest, then match positioners to DEVICE_LOC 
pos_locs = '../positioner_locations_0530v14.csv'#os.path.join(allsetdir,'positioner_locations_0530v12.csv')
positioners = Table.read(pos_locs,format='ascii.csv',header_start=0,data_start=1)
for row in positioners:
    idnam, typ, xloc, yloc, zloc, qloc, rloc, sloc = row
    device_loc_file.append(idnam)
    metro_X_file.append(xloc)
    metro_Y_file.append(yloc)

PFS_file='PFS_ID_Map.csv'
PFS=Table.read(PFS_file,format='ascii.csv',header_start=0,data_start=1)
for row in PFS:
    posids_PFS.append(row['DEVICE_ID'])
    device_loc_PFS.append(row['DEVICE_LOC'])
        
for i in range(n_pos):
    posid=posids[i]
    obsX_arr.append(float(ptl.get(posid=posid,key=['LAST_MEAS_OBS_X']))
    obsY_arr.append(float(ptl.get(posid=posid,key=['LAST_MEAS_OBS_Y']))        
    index=posids_PFS.index(posid.strip('M'))
    device_loc_this=device_loc_PFS[index]
    index2=divice_loc_arr.index(device_loc_this)
    metro_X_arr.append(metro_X_file[index2])
    metro_Y_arr.append(metro_Y_file[index2])

for i in range(n_pos):
    #   Take an image
    obsXY,peaks,fwhms,imgfile1 = m.fvc.measure(num_objects=1)
    #  Move this positioner
    ptl.quick_direct_dtdp(posids, [5.,5.], log_note='')
    # Take another image
    obsXY,peaks,fwhms,imgfile1 = m.fvc.measure(num_objects=1)
    img1=fits.open(imgfile1)
    img2=fits.open(imgfile2)
    pdb.set_trace()



            
