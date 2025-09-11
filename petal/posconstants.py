# -*- coding: utf-8 -*-
import os
import inspect
import numpy as np
import math
from datetime import datetime, timedelta
import time
import pytz
from collections import OrderedDict
import csv
try:
    from DOSlib.flags import POSITIONER_FLAGS_BITS, POSITIONER_FLAGS_MASKS, POSITIONER_FLAGS_VERBOSE
    flags_imported = True
except:
    POSITIONER_FLAGS_BITS = {}
    POSITIONER_FLAGS_MASK = {}
    POSITIONER_FLAGS_VERSOSE = {}
    flags_imported = False

"""Constants, environment variables, and convenience methods used in the
plate_control code for fiber positioners.
"""

# Required environment variables
POSITIONER_LOGS_PATH = os.environ.get('POSITIONER_LOGS_PATH') # corresponds to https://desi.lbl.gov/svn/code/focalplane/positioner_logs
FP_SETTINGS_PATH = os.environ.get('FP_SETTINGS_PATH') # corresponds to https://desi.lbl.gov/svn/code/focalplane/fp_settings

# Interpreter settings
np.set_printoptions(suppress=True) # suppress auto-scientific notation when printing np arrays

# Verson of plate_control code we are running
# Determined by reading what path in the directory structure we are running from.
# This works on assumption of the particular directory structure that we have set up in the SVN, true as of 2017-02-07.
# Should result in either 'trunk' or 'v0.31' or the like.
petal_directory = os.path.dirname(os.path.abspath(inspect.getframeinfo(inspect.currentframe()).filename))
code_version = petal_directory.split(os.path.sep)[-2]

# Directory locations
dirs = {}
if 'FOCALPLANE_HOME' in os.environ:
    home = os.environ.get('FOCALPLANE_HOME')
elif 'DESI_HOME' in os.environ:
    home = os.environ.get('DESI_HOME')
elif 'HOME' in os.environ:
    home = os.environ.get('HOME')
else:
    home = petal_directory
dirs['temp_files'] = os.path.join(home, 'fp_temp_files')
dir_keys_logs = ['pos_logs', 'fid_logs', 'ptl_logs', 'xytest_data',
                 'xytest_logs', 'xytest_plots', 'xytest_summaries',
                 'calib_logs', 'kpno', 'sequence_logs']
dir_keys_settings = ['pos_settings', 'fid_settings', 'test_settings',
                     'collision_settings', 'hwsetups', 'ptl_settings',
                     'other_settings']
for key in dir_keys_logs:
    dirs[key] = os.path.join(POSITIONER_LOGS_PATH, key)
for key in dir_keys_settings:
    dirs[key] = os.path.join(FP_SETTINGS_PATH, key)
for directory in dirs.values():
    if not os.path.isfile(directory):
        os.makedirs(directory, exist_ok=True)

# File locations
positioner_locations_file = os.path.join(petal_directory, 'positioner_locations_0530v18.csv')
small_array_locations_file = os.path.join(dirs['hwsetups'], 'SWIntegration_XY.csv')
default_collider_filename = '_collision_settings_DEFAULT.conf'
def get_keepouts_cache_path(petal_id):
    filename = f'keepouts_cache_petal_id_{petal_id}.ecsv'
    return os.path.join(dirs['temp_files'], filename)
fp_calibs_path_cache = os.path.join(dirs['temp_files'], 'latest_fp_calibs.txt')

# Lookup tables for focal plane coordinate conversions
R_lookup_path = petal_directory + os.path.sep + 'focal_surface_lookup.csv'
R_lookup_data = np.genfromtxt(R_lookup_path, comments="#", delimiter=",")

# Mapping of radial coordinate R to pseudo-radial coordinate S
# (distance along focal surface from optical axis)
def R2S_lookup(R):
    return np.interp(R,R_lookup_data[:,0],R_lookup_data[:,2],left=float('nan'))
def S2R_lookup(S):
    return np.interp(S,R_lookup_data[:,2],R_lookup_data[:,0],left=float('nan'))
def R2Z_lookup(R):
    return np.interp(R,R_lookup_data[:,0],R_lookup_data[:,1],left=float('nan'))
def Z2R_lookup(S):
    return np.interp(S,R_lookup_data[:,1],R_lookup_data[:,0],left=float('nan'))
def R2N_lookup(R):
    return np.interp(R,R_lookup_data[:,0],R_lookup_data[:,3],left=float('nan'))
def N2R_lookup(S):
    return np.interp(S,R_lookup_data[:,3],R_lookup_data[:,0],left=float('nan'))

# composite focal surface lookup methods for 3D out-of-shell QST transforms
def S2N_lookup(S):
    return R2N_lookup(S2R_lookup(S))  # return nutation angles in degrees

def S2Z_lookup(S):
    return R2Z_lookup(S2R_lookup(S))

def Z2S_lookup(Z):
    return R2S_lookup(Z2R_lookup(Z))

def N2S_lookup(N):
    return R2S_lookup(N2R_lookup(N))  # takes nutation angles in degrees

# Generic map of positioner device locs to their neighboring locs
generic_pos_neighbor_locs_path = os.path.join(petal_directory, 'generic_pos_neighbor_locs.csv')
generic_pos_neighbor_locs = {}
if os.path.exists(generic_pos_neighbor_locs_path):
    with open(generic_pos_neighbor_locs_path, 'r', newline='') as file:
        reader = csv.DictReader(file)
        neighbor_fields = [f for f in reader.fieldnames if f != 'DEVICE_LOC']
        for row in reader:
            neighbors = {int(row[f]) for f in neighbor_fields if row[f] not in ['']}
            generic_pos_neighbor_locs[int(row['DEVICE_LOC'])] = neighbors
else:
    generic_pos_neighbor_locs = 'not found: ' + generic_pos_neighbor_locs_path

# Mapping of positioner power supplies to can channels
#  2020-06-25 [JHS] For usage of new petalboxes with 20 can channels, an alternate
#  map will need to be provided here. Selection of which map to use will need to
#  be given by some configuration argument during petal initialization.
#  2024-11-11 [CAD] Since the 10 and 20 channel canids are disjoint sets, can use one map.
power_supply_canbus_map = {'V1': {'can10', 'can11', 'can13', 'can22', 'can23', 'can30', 'can31', 'can32', 'can33', 'can34', 'can35', 'can39', 'can40', 'can41', 'can42'},
                           'V2': {'can12', 'can14', 'can15', 'can16', 'can17', 'can36', 'can37', 'can38', 'can43', 'can44', 'can45', 'can46', 'can47', 'can48', 'can49'}}

# Constants
deg = '\u00b0'
mm = 'mm'
um_per_mm = 1000
deg_per_rad = 180./math.pi
rad_per_deg = math.pi/180.
timestamp_format = '%Y-%m-%dT%H:%M:%S%z' # see strftime documentation
filename_timestamp_format = '%Y%m%dT%H%M%S%z'
gear_ratio = {}
gear_ratio['namiki'] = (46.0/14.0+1)**4  # namiki    "337:1", output rotation/motor input
gear_ratio['maxon'] = 4100625.0/14641.0  # maxon     "280:1", output rotation/motor input
gear_ratio['faulhaber'] = 256.0  		 # faulhaber "256:1", output rotation/motor input
T = 0  # theta axis idx -- NOT the motor axis ID!!
P = 1  # phi axis idx -- NOT the motor axis ID!!
axis_labels = ('theta', 'phi')

# Zeno motor parameters for "linear phi"
P_zeno_speed = 100  # 30,000 RPM - see DESI 1710, Motor Speed Parameters spreadsheet
P_zeno_ramp  = 1    # 1.497 deg
# P_zeno_speed = 66   # 19,800 RPM - see DESI 1710, Motor Speed Parameters spreadsheet
# P_zeno_ramp  = 2    # 1.311 deg
# P_zeno_speed = 33   # 9,900 RPM - see DESI 1710, Motor Speed Parameters spreadsheet
# P_zeno_ramp  = 12    # 1.995 deg
P_zeno_jog = 5.0    # degrees, must be greater than ramp and backlash size

# common print function
# note the implementation may be replaced at runtime by petal.py, for logging
printfunc = print

# handle_fvc_feedback defaults
err_thresh = 0.5 # tracking error over which to disable, False means none
up_tol = 0.065 # mm over which to apply tp updates
up_tol_disabled = 0.2 #mm over which to apply tp updates to disabled positioners
up_frac = 1.0 # amount of update to apply to posTP
err_disable = False # automatically disable positioners above error threshhold
unreachable_margin = 0.2 # Additional margin beyond max patrol radius which an unreachable measurement may still be considered valid for a TP update
didnotmove_tol = 0.080 # Margin for how close to old move a positioner needs to be to be considered "nonmoving"
didnotmove_check_tol = 0.080 # Margin for how far off tracking needs to be to check for non motion

# some numeric tolerances for scheduling moves
schedule_checking_numeric_angular_tol = 0.01 # deg, equiv to about 1 um at full extension of both arms
schedule_checking_angular_tol_zeno = 0.1174 # 0.009782 deg/motor_step * 6 rows * 2 zeno jogs/normal row
near_full_range_reduced_hardstop_clearance_factor = 0.75 # applies to hardstop clearance values in special case of "near_full_range" (c.f. Axis class in posmodel.py)
max_auto_creep_distance = 10.0 # deg, fallback value to prevent huge / long creep moves in case of error in distance calculation -- only affects auto-generated creep moves
theta_hardstop_ambig_tol = 8.0 # deg, for determining when within ambiguous zone of theta hardstops
theta_hardstop_ambig_exit_margin = 5.0 # deg, additional margin to ensure getting out of ambiguous zone

# Annealing spreads out motor power consumption in time, as well as naturally reducing
# potential collision frequency. See posschedulestage.py for more info. Do *not* increase
# the anneal_density values below on a whim. These values were chosen to broadly ensure
# not too many motors spinning simultaneously.
anneal_density = {'filled': 0.5,
                  'ramped': 0.35,
                  }
max_targets_for_no_anneal = 123 # If the number of targets for each power supply is below or equal to this number, and anticollision is None or 'freeze', there is no need for annealing

# Nominal and tolerance calibration values
nominals = OrderedDict()
nominals['LENGTH_R1']        = {'value':   3.0, 'tol':    1.0}
nominals['LENGTH_R2']        = {'value':   3.0, 'tol':    1.0}
nominals['OFFSET_T']         = {'value':   0.0, 'tol':  200.0}
nominals['OFFSET_P']         = {'value':   0.0, 'tol':   50.0}
nominals['OFFSET_X']         = {'value':   0.0, 'tol': 1000.0}
nominals['OFFSET_Y']         = {'value':   0.0, 'tol': 1000.0}
nominals['PHYSICAL_RANGE_T'] = {'value': 370.0, 'tol':   50.0}
nominals['PHYSICAL_RANGE_P'] = {'value': 190.0, 'tol':   50.0}
nominals['GEAR_CALIB_T']     = {'value':   1.0, 'tol':    1.0}
nominals['GEAR_CALIB_P']     = {'value':   1.0, 'tol':    1.0}

# Tolerance for theta guesses when performing xy2tp transform
default_t_guess_tol = 30.0  # deg

# Nominal value for when positioner is sufficiently off-center, to ensure that
# a theta measurement by the FVC will be valid. Also a tol value meant for identifying
# when the off-center phi would cause an intolerably large theta jump.
_off_center_threshold_mm = 0.5  # radial distance off-center
_ctrd_phi_theta_change_tol_mm = 0.1 # allowable max positioning error induced by sudden theta change while centered phi
_nom_max_r = nominals['LENGTH_R1']['value'] + nominals['LENGTH_R2']['value']
_min_length_r2 = nominals['LENGTH_R2']['value'] - nominals['LENGTH_R2']['tol']
phi_off_center_threshold = 180 - math.floor(_off_center_threshold_mm / _min_length_r2 * deg_per_rad)
ctrd_phi_theta_change_tol = math.ceil(_ctrd_phi_theta_change_tol_mm / _nom_max_r * deg_per_rad)

# Hardware (operations) States
PETAL_OPS_STATES = {'INITIALIZED' : OrderedDict({#'CAN_EN':(['on','on'], 1.0), #CAN Power ON
                                                 'GFA_FAN':({'inlet':['off',0],'outlet':['off',0]}, 1.0), #GFA Fan Power OFF
                                                 'GFAPWR_EN':('off', 60.0),  #GFA Power Enable OFF
                                                 'TEC_CTRL':('off', 15.0), #TEC Power EN OFF
                                                 'BUFFERS':(['on','on'], 1.0), #SYNC Buffer EN ON
                                                 #GFA CCD OFF
                                                 #GFA CCD Voltages EN OFF
                                                 #TEC Control EN OFF - handeled by camera.py
                                                 #PetalBox Power ON - controlled by physical raritan switch
                                                 'PS1_EN':('off', 1.0), #Positioner Power EN OFF
                                                 'PS2_EN':('off', 1.0)}),
                    'STANDBY' : OrderedDict({#'CAN_EN':(['on','on'], 1.0), #CAN Power ON
                                             'GFAPWR_EN':('off', 60.0), #GFA Power Enable OFF
                                             'GFA_FAN':({'inlet':['off',0],'outlet':['off',0]}, 1.0), #GFA Fan Power OFF
                                             'TEC_CTRL': ('off', 15.0), #TEC Power EN OFF
                                             'BUFFERS':(['on','on'], 1.0), #SYNC Buffer EN ON
                                             #GFA CCD OFF
                                             #GFA CCD Voltages EN OFF
                                             #TEC Control EN OFF - handeled by camera.py
                                             #PetalBox Power ON - controlled by physical raritan switch
                                             'PS1_EN':('off', 1.0), #Positioner Power EN OFF
                                             'PS2_EN':('off', 1.0)}),
                    'READY' : OrderedDict({#'CAN_EN':(['on','on'], 1.0), #CAN Power ON
                                           'GFA_FAN':({'inlet':['on',15],'outlet':['on',15]}, 1.0), #GFA Fan Power ON
                                           'GFAPWR_EN':('on', 60.0), #GFA Power Enable ON
                                           'TEC_CTRL': ('off', 15.0), #TEC Power EN OFF for now
                                           'BUFFERS':(['on','on'], 1.0), #SYNC Buffer EN ON
                                           #GFA CCD OFF
                                           #GFA CCD Voltages EN OFF
                                           #TEC Control EN ON - controlled by camera.py
                                           #PetalBox Power ON - controlled by physical raritan switch
                                           'PS1_EN': ('off', 1.0), #Positioner Power EN OFF
                                           'PS2_EN': ('off', 1.0)}),
                    'OBSERVING' : OrderedDict({#'CAN_EN':(['on','on'], 1.0), #CAN Power ON
                                               'GFA_FAN':({'inlet':['on',15],'outlet':['on',15]}, 1.0), #GFA Fan Power ON
                                               'GFAPWR_EN':('on', 60.0), #GFA Power Enable ON
                                               'TEC_CTRL':('off', 15.0), #TEC Power EN OFF for now
                                               'BUFFERS':(['on','on'], 1.0), #SYNC Buffer EN ON
                                               #GFA CCD ON
                                               #GFA CCD Voltages EN ON
                                               #TEC Control EN ON - controlled by camera.py
                                               #PetalBox Power ON - controlled by physical raritan switch
                                               'PS1_EN':('on', 1.0), #Positioner Power EN ON
                                               'PS2_EN':('on', 1.0)})}


# Conservatively accessible angle ranges (intended to be valid for any basically
# functional postioner, and for which a seed calibration is at least roughly known).
# All angles in deg.
conservative_min_posintP_extended = 10.0
conservative_min_posintP_retracted = 120.0
conservative_max_posintP = 170.0
conservative_min_posintT = -170.0
conservative_max_posintT = 170.0

# keepout envelope expansion parameter keys
keepout_expansion_keys = ['KEEPOUT_EXPANSION_PHI_RADIAL',
                          'KEEPOUT_EXPANSION_PHI_ANGULAR',
                          'KEEPOUT_EXPANSION_THETA_RADIAL',
                          'KEEPOUT_EXPANSION_THETA_ANGULAR']
keepout_keys = keepout_expansion_keys + ['CLASSIFIED_AS_RETRACTED']

# other "calib" keys, meaning keys that for whatever historical or other reason,
# are kept in the "calib" tables in the online db rather than "moves"
other_pos_calib_keys = {'TOTAL_LIMIT_SEEKS_T', 'TOTAL_LIMIT_SEEKS_P',
                        'LAST_PRIMARY_HARDSTOP_DIR_T', 'LAST_PRIMARY_HARDSTOP_DIR_P',
                        'CALIB_NOTE', 'DEVICE_CLASSIFIED_NONFUNCTIONAL', 'FIBER_INTACT'}
fiducial_calib_keys = {'DUTY_STATE', 'DUTY_DEFAULT_ON', 'DUTY_DEFAULT_OFF'}
zenomotor_keys = {'ZENO_MOTOR_P', 'ZENO_MOTOR_T', 'SZ_CW_P', 'SZ_CCW_P',  'SZ_CW_T', 'SZ_CCW_T'}
posmodel_range_names = {'targetable_range_posintT', 'targetable_range_posintP',
                        'full_range_posintT', 'full_range_posintP',
                        'theta_hardstop_ambiguous_zone'}
posmodel_keys = {'in_theta_hardstop_ambiguous_zone',
                 'abs_shaft_speed_cruise_T', 'abs_shaft_speed_cruise_P',
                 'abs_shaft_spinupdown_distance_T', 'abs_shaft_spinupdown_distance_P'}

for name in posmodel_range_names:
    posmodel_keys |= {f'max_{name}', f'min_{name}'}

# test for whether certain posstate keys are classified as "calibration" vals
calib_keys = set(nominals.keys()) | set(keepout_keys) | other_pos_calib_keys | fiducial_calib_keys | zenomotor_keys
def is_calib_key(key):
    return key.upper() in calib_keys

# keys / default types in constants DB
constants_keys = {"ALLOW_EXCEED_LIMITS": False,
                  "ANTIBACKLASH_FINAL_MOVE_DIR_P":  1,
                  "ANTIBACKLASH_FINAL_MOVE_DIR_T": -1,
                  "ANTIBACKLASH_ON": True,
                  "BACKLASH": 3.0,
                  "BUMP_CCW_FLG": False,
                  "BUMP_CW_FLG": False,
                  "CREEP_PERIOD": 2,
                  "CREEP_TO_LIMITS": False,
                  "CURR_CREEP": 70,
                  "CURR_CRUISE": 70,
                  "CURR_HOLD": 0,
                  "CURR_SPIN_UP_DOWN": 70,
                  "FINAL_CREEP_ON": True,
                  "GEAR_TYPE_P": "namiki",
                  "GEAR_TYPE_T": "namiki",
                  "LIMIT_SEEK_EXCEED_RANGE_FACTOR": 1.3,
                  "MIN_DIST_AT_CRUISE_SPEED": 180.0,
                  "MOTOR_CCW_DIR_P": -1,
                  "MOTOR_CCW_DIR_T": -1,
                  "MOTOR_ID_P": 0,
                  "MOTOR_ID_T": 1,
                  "ONLY_CREEP": False,
                  "PRINCIPLE_HARDSTOP_CLEARANCE_P": 3.0,
                  "PRINCIPLE_HARDSTOP_CLEARANCE_T": 3.0,
                  "PRINCIPLE_HARDSTOP_DIR_P": 1,
                  "PRINCIPLE_HARDSTOP_DIR_T": -1,
                  "SECONDARY_HARDSTOP_CLEARANCE_P": 3.0,
                  "SECONDARY_HARDSTOP_CLEARANCE_T": 3.0,
                  "SPINUPDOWN_PERIOD": 12,
                  }
def is_constants_key(key):
    return key.upper() in constants_keys

# data fields which must *always* be accompanied by a comment
require_comment_to_store = {'CTRL_ENABLED',
                            'DEVICE_CLASSIFIED_NONFUNCTIONAL',
                            'FIBER_INTACT'}

# state data fields associated with "late" committing to database
late_commit_defaults = {'OBS_X': None,
                        'OBS_Y': None,
                        'PTL_X': None,
                        'PTL_Y': None,
                        'PTL_Z': None,
                        'FLAGS': None,
                        'POSTSCRIPT': None}

# ordered lists of motor parameters (in order that petalcontroller expects them)
# interface is funky, for example the doubling in some cases for the two motors
ordered_motor_current_keys = ['CURR_SPIN_UP_DOWN', 'CURR_CRUISE', 'CURR_CREEP', 'CURR_HOLD'] * 2
ordered_motor_period_keys = ['CREEP_PERIOD'] * 2 + ['SPINUPDOWN_PERIOD']
motor_param_keys = set(ordered_motor_current_keys + ordered_motor_period_keys)

# state keys which posmodel keeps a cache of (for speed)
keys_cached_in_posmodel = {'BACKLASH',
                           'PRINCIPLE_HARDSTOP_CLEARANCE_T',
                           'PRINCIPLE_HARDSTOP_CLEARANCE_P',
                           'SECONDARY_HARDSTOP_CLEARANCE_T',
                           'SECONDARY_HARDSTOP_CLEARANCE_P',
                           'ANTIBACKLASH_FINAL_MOVE_DIR_T',
                           'ANTIBACKLASH_FINAL_MOVE_DIR_P',
                           'PRINCIPLE_HARDSTOP_DIR_T',
                           'PRINCIPLE_HARDSTOP_DIR_P',
                           'GEAR_TYPE_T',
                           'GEAR_TYPE_P',
                           'MOTOR_CCW_DIR_T',
                           'MOTOR_CCW_DIR_P',
                           'GEAR_CALIB_T',
                           'GEAR_CALIB_P',
                           'PHYSICAL_RANGE_T',
                           'PHYSICAL_RANGE_P',
                           }
def is_cached_in_posmodel(key):
    return key.upper() in keys_cached_in_posmodel

# performance grade letters
grades = ['A', 'B', 'C', 'D', 'F', 'N/A']

# move command strings
delta_move_commands = {'dQdS', 'obsdXdY', 'poslocdXdY', 'dTdP'}
abs_move_commands = {'QS', 'obsXY', 'ptlXY', 'poslocXY', 'poslocTP', 'posintTP'}
valid_move_commands = abs_move_commands | delta_move_commands

# common formatting / split-up names
coord_formats = {key: '6.1f' for key in ['posintTP', 'poslocTP']}
coord_formats.update({key: '7.3f' for key in ['poslocXY']})
coord_formats.update({key: '8.3f' for key in ['QS', 'flatXY', 'obsXY', 'ptlXY']})
coord_pair2single = {c: (c[:-1], c[:-2] + c[-1]) for c in coord_formats}
coord_formats.update({s[i]: coord_formats[c] for i in [0,1] for c, s in coord_pair2single.items()})
single_coords = set(coord_formats) - set(coord_pair2single)
coord_single2pair = {}
for pair_key, coord_tuple in coord_pair2single.items():
    for single_key in coord_tuple:
        coord_single2pair[single_key] = pair_key

def decipher_posflags(flags, sep=';', verbose=True):
    '''translates posflag to readable reasons, bits taken from petal.py
    simple problem of locating the leftmost set bit, always 0b100 on the
    right input flags. presence of non-positioner bits from FVC/PM is OK
    input flags must be an array-like or list-like object'''

    def decipher_flag(flag, sep, verbose):
            if flags_imported:
                status_list = []
                for key, val in POSITIONER_FLAGS_MASKS.items():
                    if (flag & val) != 0:
                        if verbose:
                            status_list.append(POSITIONER_FLAGS_VERBOSE[key])
                        else:
                            status_list.append(key)
                return sep.join(status_list)
            else:
                return 'Flags not imported'
    flags = np.array(flags).reshape(-1,).astype(int)  # 1d to enable indexing

    return [decipher_flag(flag, sep, verbose) for flag in flags]


class collision_case(object):
    """Enumeration of collision cases. The I, II, and III cases are described in
    detail in DESI-0899.
    """
    def __init__(self):
        self.I    = 0  # no collision
        self.II   = 1  # phi arm against neighboring phi arm
        self.III  = 2  # phi arm against neighboring central body
        self.IV   = 3  # phi arm against neighboring circular keepout
        self.GFA  = 4  # phi arm against the GFA fixed keepout envelope
        self.PTL  = 5  # phi arm against the Petal edge keepout envelope
        self.fixed_cases = {self.PTL, self.GFA}
        self.names = {self.I:   'no collision',
                      self.II:  'phi',
                      self.III: 'central body',
                      self.IV:  'circular keepout',
                      self.GFA: 'GFA',
                      self.PTL: 'PTL'}
        self.fixed_case_names = {self.names[c] for c in self.fixed_cases}
case = collision_case()

# Collision resolution methods
nonfreeze_adjustment_methods = ['pause',
                                'extend_A', 'retract_A',
                                'extend_B', 'retract_B',
                                'rot_ccw_A', 'rot_cw_A',
                                'rot_ccw_B', 'rot_cw_B',
                                'repel_ccw_A', 'repel_cw_A',
                                'repel_ccw_B', 'repel_cw_B']
all_adjustment_methods = nonfreeze_adjustment_methods + ['freeze']
useless_with_unmoving_neighbor = {'pause'} | {m for m in nonfreeze_adjustment_methods if 'repel' in m}
useless_with_fixed_boundary = useless_with_unmoving_neighbor | {m for m in nonfreeze_adjustment_methods if 'extend' in m}
num_timesteps_clearance_margin = 2  # this value * PosCollider.timestep --> small extra wait for a neighbor to move out of way

# Initial polygon debouncing settings
num_timesteps_ignore_overlap = 1  # ignore collisions during these first few timesteps (just during debounce stage)
debounce_polys_distance = 5  # deg, for attempts to slightly step one polygon off another when barely touching


# Convenience methods
rotmat2D = lambda angle: [math.cos(angle*rad_per_deg), - math.sin(angle*rad_per_deg), math.sin(angle*rad_per_deg), math.cos(angle*rad_per_deg)]
transpose = lambda matrix: [list(x) for x in zip(*matrix)] # slower than numpy.transpose() if you *already* have a numpy array, but much faster on small python lists

def sign(x):
    """Return the sign of the value x as +1, -1, or 0."""
    if x > 0.:
        return 1
    elif x < 0.:
        return -1
    else:
        return 0

def join_notes(*args):
    '''Concatenate items into a "note" string with standard format. A list or
    tuple arg is treated as a single "item". So for example if you want the
    subelements of a list "joined", then argue it expanded, like *mylist
    '''
    separator = '; '
    if len(args) == 0:
        return ''
    elif len(args) == 1:
        return str(args[0])
    strings = (str(x) for x in args if x != '')
    return separator.join(strings)

def linspace(start,stop,num):
    """Return a list of floats linearly spaced from start to stop (inclusive).
    List has num elements."""
    return [i*(stop-start)/(num-1)+start for i in range(num)]

# Functions for handling mixes of [M][N] vs [M] dimension lists
def listify(uv, keep_flat=False):
    """Turn [u,v] into [[u],[v]], if it isn't already.
    Turn uv into a list [uv] if it isn't already.
    In the special case where uv is a single item list, it remains so in the return.
    A boolean is returned saying whether the item was modified.
    If optional argument 'flat' is true, then [u,v] remains [u,v].
    """
    if not(isinstance(uv,list)):
        return [uv], True
    if len(uv) == 1 or keep_flat:
        return uv, False
    new_uv = []
    was_listified = False
    for i in range(len(uv)):
        if not(isinstance(uv[i],list)):
            new_uv.append([uv[i]])
            was_listified = True
        else:
            new_uv.append(uv[i].copy())
    return new_uv, was_listified

def delistify(uv):
    """Turn [[u],[v]] into [u,v], if it isn't already.
    For a non-list, there is no modification.
    For a list with one item only, the element inside is returned.
    For a list with multiple items, the return is a flat list of multiple items.
    """
    if not(isinstance(uv,list)):
        return uv
    if len(uv) == 1:
        return uv[0]
    new_uv = []
    for i in range(len(uv)):
        if isinstance(uv[i],list):
            new_uv.extend(uv[i][:])
        else:
            new_uv.append(uv[i])
    return new_uv

def listify2d(uv):
    """If uv is [u,v], return [[u,v]]. Otherwise return uv.
    I.e., so [[u1,v1],[u2,v2],...] would stay the same
    """
    if isinstance(uv,list):
        if not(isinstance(uv[0],list)) and len(uv) == 2:
            return [uv]
        else:
            return uv

def concat_lists_of_lists(L1, L2):
    """Always output an [N][M] list, for M dimensional coordinate inputs initial and new.
    E.g.
    [] + [u,v] --> [[u,v]]
    [] + [[u1,v1],[u2,v2],...] --> same
    [[u1,v1],[u2,v2]] + [u3,v3] --> [[u1,v1],[u2,v2],[u3,v3]]
    [[u1,v1],[u2,v2]] + [[u3,v3],[u4,v4],...] --> [[u1,v1],[u2,v2],[u3,v3],[u4,v4],...]
    """
    if not(L1):
        L1 = []
    elif not(isinstance(L1[0],list)):
        L1 = [L1]
    if not(L2):
        L2 = []
    elif not(isinstance(L2[0],list)):
        L2 = [L2]
    return L1 + L2

# Enumeration of verbosity level to stdout
not_verbose = 0
verbose = 1
very_verbose = 2


def is_verbose(verbosity_enum):
    boole = True
    if verbosity_enum == not_verbose:
        boole = False
    return boole


def is_very_verbose(verbosity_enum):
    boole = False
    if verbosity_enum == very_verbose:
        boole = True
    return boole


# timestamp functions
def now():
    # current TZ unaware local time (), then TZ aware in local timezone
    return datetime.now().astimezone()


def utcnow():
    return datetime.now().astimezone(pytz.timezone('UTC'))


def timestamp_str(t=None):
    if t is None:
        t = utcnow()
    return t.strftime(timestamp_format)


def filename_timestamp_str(t=None):
    if t is None:
        t = utcnow()
    return t.strftime(filename_timestamp_format)


def dir_date_str(t=None):
    '''returns date string for the directory name, changes at noon Arizona'''
    if t is None:
        t = now()
    t = t.astimezone(pytz.timezone('America/Phoenix')) - timedelta(hours=12)
    return f'{t.year:04}{t.month:02}{t.day:02}'

def compact_timestamp(nowtime=None, basetime=1582915648):
    '''Compact, readable time code. Default return string is six characters
    in length; will exceed this length at basetime + 69 years. Precision is
    rounded to seconds. Default argument baselines it at a recent time on
    Feb 28, 2020, 10:47 AM PST. The argument nowtime is just there for testing.
    '''
    maxchar = 6
    nowtime = time.time() if not nowtime else nowtime
    relative_now = math.floor(nowtime - basetime)
    converted = np.base_repr(relative_now, base=36)
    padded = converted.rjust(maxchar,'0') if len(converted) < maxchar else converted
    return padded

# other misc functions
def ordinal_str(number):
    '''
    Returns a string of the number plus 'st', 'nd', 'rd', 'th' as appropriate.
    '''
    numstr = str(number)
    last_digit = numstr[-1]
    if last_digit == '1':
        return numstr + 'st'
    if last_digit == '2':
        return numstr + 'nd'
    if last_digit == '3':
        return numstr + 'rd'
    return numstr + 'th'

def is_integer(x):
    return isinstance(x, (int, np.integer))

def is_float(x):
    return isinstance(x, (float, np.floating))

def is_string(x):
    return isinstance(x, (str))

def is_boolean(x):
    return x in [True, False, 0, 1] or str(x).lower() in ['true', 'false', '0', '1']

def is_none(x):
    return x in [None, 'None', 'none', 'NONE']

def is_collection(x):
    if isinstance(x, (dict, list, tuple, set)):
        return True
    if is_integer(x) or is_float(x) or is_string(x) or is_boolean(x):
        return False
    return '__len__' in dir(x)

def boolean(x):
    '''Cast input to boolean.'''
    if x in [True, False]:
        return x
    if x == None or is_integer(x) or is_float(x):
        return bool(x)
    if is_string(x):
        return x.lower() not in {'false', '0', 'none', 'null', 'no', 'n'}
    if is_collection(x):
        return len(x) > 0
    assert False, f'posconstants.boolean(): undefined interpretation for {x}'

# style info for plotting positioners
plot_styles = {
    'ferrule':
        {'linestyle' : '-',
         'linewidth' : 0.5,
         'edgecolor' : 'green',
         'facecolor' : 'none'},

    'phi arm':
        {'linestyle' : '-',
         'linewidth' : 1,
         'edgecolor' : 'green',
         'facecolor' : 'none'},

    'phi arc':
        {'linestyle' : '--',
         'linewidth' : 1,
         'edgecolor' : 'green',
         'facecolor' : 'none'},

    'central body':
        {'linestyle' : '-',
         'linewidth' : 1,
         'edgecolor' : 'green',
         'facecolor' : 'none'},

    'positioner element unbold':
        {'linestyle' : '--',
         'linewidth' : 0.5,
         'edgecolor' : '0.6',
         'facecolor' : 'none'},

    'collision':
        {'linestyle' : '-',
         'linewidth' : 2,
         'edgecolor' : 'red',
         'facecolor' : 'none'},

    'frozen':
        {'linestyle' : '-',
         'linewidth' : 2,
         'edgecolor' : 'blue',
         'facecolor' : 'none'},

    'line t0':
        {'linestyle' : '-.',
         'linewidth' : 0.5,
         'edgecolor' : 'gray',
         'facecolor' : 'none'},

    'arm lines':
        {'linestyle' : '--',
         'linewidth' : 0.7,
         'edgecolor' : 'black',
         'facecolor' : 'none'},

    'Eo':
        {'linestyle' : '-',
         'linewidth' : 0.5,
         'edgecolor' : '0.9',
         'facecolor' : 'none'},

    'Eo bold':
        {'linestyle' : '-',
         'linewidth' : 1,
         'edgecolor' : 'green',
         'facecolor' : 'none'},

    'Ei':
        {'linestyle' : '-',
         'linewidth' : 0.5,
         'edgecolor' : '0.9',
         'facecolor' : 'none'},

    'Ee':
        {'linestyle' : '-',
         'linewidth' : 0.5,
         'edgecolor' : '0.9',
         'facecolor' : 'none'},

    'PTL':
        {'linestyle' : '--',
         'linewidth' : 1,
         'edgecolor' : '0.5',
         'facecolor' : 'none'},

    'GFA':
        {'linestyle' : '--',
         'linewidth' : 1,
         'edgecolor' : '0.5',
         'facecolor' : 'none'}
    }
