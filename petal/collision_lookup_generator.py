"""This module generates a set of hashable codes for collision lookup.

[These notes are current as of 9/12/2018]

POS-POS CASE
To check for collision between any two positioners:
    1. Form the code based on:
        ... center-to-center distance between positioners
        ... theta of each positioner
        ... phi w.r.t. theta for each positioner
    2. See if the code is in the collisions lookup set.

    CODE FORMAT: 'rrr-ttt-ppp-qqq-nnn', all digits are integers.
        rrr ... [um]  rrr = ctr_to_ctr_dist - 10.1 mm
                      (10.1 mm is the nominal distance between the centers of two positioners)
        ttt ... [deg] ttt = obsT1
        ppp ... [deg] ppp = posP1
        qqq ... [deg] qqq = obsT2 (neighbor)
        nnn ... [deg] nnn = posP2
        [WHAT ABOUT VARIATION OF R1 ON EACH POSITIONER?]

POS-FIXED CASE
To check for collision between a positioner and fixed:
    1. Form the code based on:
        ... device_location_id (i.e., hole number in petal)
        ... dx,dy center variation from nominal
        ... obsT, obsP
    2. See if the code is in the collisions lookup set.

    CODE FORMAT: 'hhh-xxx-yyy-ttt-ppp', all digits are integers.
        hhh ... integer device location (hole) id number
        xxx ... [um]  xxx = actual_obsX - nominal_X
        yyy ... [um]  yyy = actual_obsY - nominal_Y
        ttt ... [deg] ttt = obsT
        ppp ... [deg] ppp = posP
        [WHAT ABOUT VARIATION OF R1?]
"""

# PSEUDO-CODE
# Divest much of polygon-handling code from poscollider.py into this module.
# Initialize, bringing in the usual poscollider stuff like polygon definitions etc.
#
# POS-POS calculations:
# for r in the range [0, 50, 100 ... 950] um
#   for t in the range [0, 1, 2 ... 360] deg
#       for p in the range [60, 61, 62 ... 210] deg
#           for n in the range [60, 61, 62 ... 210] deg
#               if spatial collision
#                   form code
#                   add to set
#
# Start much **much** coarser than this and see how long it takes to generate / how big is the set.
# Issue: I really want 0.5 deg ~ 53 um resolution @ fully-extended phi. But that may be unreasonable.
#
# POS-FIXED calculations:
# Import petal layout data (device_location_id vs. obsX,obsY).
# Generate a set of simulation positioners filling all locations.
# for h in the range of device location ids
#   for x in the range [0, 50, 100 ... 950] um
#       for y in the range [0, 50, 100 ... 950] um
#           etc...
#

import time
from astropy.io import ascii
from matplotlib import pyplot as plt
from matplotlib.pyplot import show
import poscollider
import numpy as np

# Global variables for making pos-pos table
nominal_r = 10.1
R1 = 3.0
pos = poscollider.PosCollider()
pos._load_keepouts()

# For making pos-fixed table
posloc = ascii.read('positioner_locations_0530v14.csv')
pos_ind = np.where(posloc['device_type'] == 'POS')[0] # get only the positioners
pos_loc_id = posloc['device_location_id'][pos_ind]    # location id of positioners
xc, yc = posloc['X'][pos_ind], posloc['Y'][pos_ind] # center locations of positioners

ptl = pos.fixed_neighbor_keepouts[5] # petal keepout
gfa = pos.fixed_neighbor_keepouts[4] # gfa keepout

def place_phi_arm(poslocT, poslocP, x0=0, y0=0):

    """
    copied from poscollider.py
    """

    poly = pos.keepout_P.rotated(poslocP)
    poly = poly.translated(R1, 0)
    poly = poly.rotated(poslocT)
    poly = poly.translated(x0, y0)
    return poly

def place_central_body(poslocT, x0=0, y0=0):

    """
    copied from poscollider.py
    """

    poly = pos.keepout_T.rotated(poslocT)
    poly = poly.translated(x0, y0)
    return poly

def plot_placements(x1, y1, x2, y2, poslocT1, posintP1, poslocT2, posintP2):

    """
    plot configurations of 2 positioners (arm+central body)

    x1, y1 = (x, y) positions for pos1
    x2, y2 = (x, y) positions for pos2
    obsT1, posP1 = theta and phi angles for pos1
    obsT2, posP2 = theta and phi angles for pos2
    """
    # create pos1 with angle obsT1 and obsP1 angle
    central1 = place_central_body(poslocT1, x0=x1, y0=y1)
    arm1 = place_phi_arm(poslocT1, posintP1, x0=x1, y0=y1)

    # create pos2 with angle obsT2 and posP2 angle
    central2 = place_central_body(poslocT2, x0=x2, y0=y2)
    arm2 = place_phi_arm(poslocT2, posintP2, x0=x2, y0=y2)

    # plot the 2 positioners
    plt.plot(central1.points[0], central1.points[1], 'k-')
    plt.plot(arm1.points[0], arm1.points[1], 'k-')

    plt.plot(central2.points[0], central2.points[1], 'r-')
    plt.plot(arm2.points[0], arm2.points[1], 'r-')

def around_pos1(dr, theta1, phi1, theta2, phi2):

    pos2_circle = np.arange(0, 360, 60)
    unitX = np.cos(np.pi*pos2_circle/180)
    unitY = np.sin(np.pi*pos2_circle/180)

    rnew = 10.1 + dr/1000.
    xnew = rnew*unitX
    ynew = rnew*unitY

    central0 = place_central_body(theta1)
    central1 = place_central_body(theta2, x0=xnew[0], y0=ynew[0])
    central2 = place_central_body(theta2, x0=xnew[1], y0=ynew[1])
    central3 = place_central_body(theta2, x0=xnew[2], y0=ynew[2])
    central4 = place_central_body(theta2, x0=xnew[3], y0=ynew[3])
    central5 = place_central_body(theta2, x0=xnew[4], y0=ynew[4])
    central6 = place_central_body(theta2, x0=xnew[5], y0=ynew[5])

    arm0 = place_phi_arm(theta1, phi1)
    arm1 = place_phi_arm(theta2, phi2, x0=xnew[0], y0=ynew[0])
    arm2 = place_phi_arm(theta2, phi2, x0=xnew[1], y0=ynew[1])
    arm3 = place_phi_arm(theta2, phi2, x0=xnew[2], y0=ynew[2])
    arm4 = place_phi_arm(theta2, phi2, x0=xnew[3], y0=ynew[3])
    arm5 = place_phi_arm(theta2, phi2, x0=xnew[4], y0=ynew[4])
    arm6 = place_phi_arm(theta2, phi2, x0=xnew[5], y0=ynew[5])

    plt.text(0, 0, '1', fontsize=15, color='r')
    plt.text(xnew[0], ynew[0], '2', fontsize=15, color='r')
    plt.text(xnew[1], ynew[1], '3', fontsize=15, color='r')
    plt.text(xnew[2], ynew[2], '4', fontsize=15, color='r')
    plt.text(xnew[3], ynew[3], '5', fontsize=15, color='r')
    plt.text(xnew[4], ynew[4], '6', fontsize=15, color='r')
    plt.text(xnew[5], ynew[5], '7', fontsize=15, color='r')

    plt.plot(central0.points[0], central0.points[1], 'k-')
    plt.plot(central1.points[0], central1.points[1], 'k-')
    plt.plot(central2.points[0], central2.points[1], 'k-')
    plt.plot(central3.points[0], central3.points[1], 'k-')
    plt.plot(central4.points[0], central4.points[1], 'k-')
    plt.plot(central5.points[0], central5.points[1], 'k-')
    plt.plot(central6.points[0], central6.points[1], 'k-')

    plt.plot(arm0.points[0], arm0.points[1], 'k-')
    plt.plot(arm1.points[0], arm1.points[1], 'k-')
    plt.plot(arm2.points[0], arm2.points[1], 'k-')
    plt.plot(arm3.points[0], arm3.points[1], 'k-')
    plt.plot(arm4.points[0], arm4.points[1], 'k-')
    plt.plot(arm5.points[0], arm5.points[1], 'k-')
    plt.plot(arm6.points[0], arm6.points[1], 'k-')

def make_code_pp(dr, poslocT1, posintP1, poslocT2, posintP2):

    # ~ 7 us to run on my laptop

    """
    dr    = [um] center-to-center distance between 2 positioners
    obsT1 = [deg] theta angle of pos1
    posP1 = [deg] phi angle of pos1
    obsT2 = [deg] theta angle of pos2
    posP2 = [deg] phi angle of pos2
    """

    if poslocT1 < 0: poslocT1 = 360 + poslocT1
    if poslocT2 < 0: poslocT2 = 360 + poslocT2

    rrr = str(int(dr)).zfill(3)     # ~1 us
    ttt1 = str(int(poslocT1)).zfill(3) # ...
    ppp = str(int(posintP1)).zfill(3)
    ttt2 = str(int(poslocT2)).zfill(3) # ...
    nnn = str(int(posintP2)).zfill(3)

    # join() is faster than using '+'
    code = '-'.join([rrr, ttt1, ppp, ttt2, nnn])
    return code

def make_table_pp(obsT1_start=0, obsT1_end=360, obsT1_step=5, \
                  obsT2_start=0, obsT2_end=360, obsT2_step=5,
                  posP1_start=-20, posP1_end=210, posP1_step=5, \
                  posP2_start=-20, posP2_end=210, posP2_step=5,
                  dr_step=50):

    """
    dr_step in um (for now, dr_start=0 and dr_end=1000 hardcoded).
    5 deg res: ~18 hrs (on lear)

    For now pickle save the output table separately.
    File naming convention: table_pp_[dr_step]_[obsT1_step]_[obsT2_step]_[posP1_step]_[posP2_step].out
    """

    table_pp = {}
    n1, n2, n3 = 0, 0, 0

    start = time.time()

    # 60 deg increment is chosen because a closer look at the petal geometry
    # suggests that each positioner is surrounded (evenly) by at most 6 other positioners,
    # therefore giving 360/6 = 60 deg spacing.
    pos2_circle = np.arange(0, 360, 60)
    unitX = np.cos(np.pi*pos2_circle/180)
    unitY = np.sin(np.pi*pos2_circle/180)

    # rotating pos1 by (theta, phi)
    for obsT1 in range(obsT1_start, obsT1_end, obsT1_step):
        for posP1 in range(posP1_start, posP1_end, posP1_step):

            central1 = place_central_body(obsT1)
            arm1 = place_phi_arm(obsT1, posP1)

            # creating pos2 at a range of (theta, phi) angles
            # displaced (dr, dx, dy) away from pos1
            for dr in range(0, 1000, dr_step):
                rnew = nominal_r + dr/1000.
                xnew = rnew*unitX # dx distance of pos2 from pos1 at [0, 60, 120, .., 300] deg
                ynew = rnew*unitY # dy distance ...

                # rotating pos2 by (theta, phi)
                for obsT2 in range(obsT2_start, obsT2_end, obsT2_step):
                    for posP2 in range(posP2_start, posP2_end, posP2_step):

                        # translating pos2 by (dx, dy)
                        for i in range(len(xnew)):
                            central2 = place_central_body(obsT2, x0=xnew[i], y0=ynew[i])
                            arm2 = place_phi_arm(obsT2, posP2, x0=xnew[i], y0=ynew[i])

                            # case II
                            if arm1.collides_with(central2):

                                # Remove the relative angle between pos1 and pos2 by subtracting out this angle
                                # from obsT1 and obsT2:
                                # obsT1 = obsT1-pos2_circle[i]
                                # obsT2 = obsT2-pos2_circle[i]

                                code_pp = make_code_pp(dr, obsT1-pos2_circle[i], posP1, obsT2-pos2_circle[i], posP2)
                                table_pp[code_pp] = 2
                                n2 += 1

                            # case III
                            elif arm2.collides_with(central1):
                                code_pp = make_code_pp(dr, obsT1-pos2_circle[i], posP1, obsT2-pos2_circle[i], posP2)
                                table_pp[code_pp] = 3
                                n3 += 1

                            # case I
                            elif arm1.collides_with(arm2):
                                code_pp = make_code_pp(dr, obsT1-pos2_circle[i], posP1, obsT2-pos2_circle[i], posP2)
                                table_pp[code_pp] = 1
                                n1 += 1

    end = time.time()
    print((end-start)/60.)
    print(len(table_pp))
    print(n1, n2, n3)

    return table_pp

################################### POS-FIXED ###################################
def identify_id_near_fixed(start, end, write=True):

    """
    One-time run: manually identity the location id of positioners (really the hole ids)
    near the petal or a gfa. Outputs a file containing these location ids.

    start = index corresponding to xc/yc, from 0-499
    end = [same]
    """

    f = open('ids_near_fixed.dat', 'a') # where the loc id will be written

    # loop through positioner in a petal and plot where it is on the petal
    # requests for user input after each plot to determine whether to record down
    # the center positions of that positioner or not.
    for i in range(start, end):
        plt.plot(xc, yc, 'ko')

        plt.plot(ptl.points[0], ptl.points[1])
        plt.plot(gfa.points[0], gfa.points[1])
        plt.plot(xc[i], yc[i], 'ro')
        show(block=False)
        user = input('record?')
        if write:
            if user == 'y':
                f.write(str(pos_loc_id[i]) + '\n')
    f.close()

def plot_id_near_fixed(file, newout=None):

    """
    Also one-time run; plots the recorded loc ids and saves the center positions of
    the positioners at these loc ids to "newout".

    file = "ids_near_fixed.dat"
    newout = "pos_near_fixed.csv"
    """

    # extract the loc ids from "ids_near_fixed.dat"
    f = open(file, 'r')
    id_near_fixed = []

    for line in f:
        id_near_fixed.append(int(line.strip('\n')))

    f.close()

    # get the center positions (xc, yc) of the loc ids
    x_near_fixed, y_near_fixed = [], []
    for i in range(len(pos_loc_id)):
        if pos_loc_id[i] in id_near_fixed:
            x_near_fixed.append(xc[i])
            y_near_fixed.append(yc[i])

    plt.plot(ptl.points[0], ptl.points[1])
    plt.plot(gfa.points[0], gfa.points[1])
    plt.plot(xc, yc, 'ko')
    plt.plot(x_near_fixed, y_near_fixed, 'ro')
    plt.savefig('pos_near_fixed.png')

    # write out the (loc_id, xc, yc) of the loc ids if an output file is provided
    if newout != None:
        f = open(newout, 'w')
        for i in range(len(id_near_fixed)):
            f.write("%d,%0.6f,%0.6f\n" % (id_near_fixed[i], x_near_fixed[i], y_near_fixed[i]))
        f.close()

def make_code_pf(loc_id, dx, dy, poslocT, poslocP):

    """
    loc_id = 3-digit location id
    dx =
    dy =
    obsT = theta angle of positioner
    posP = phi angle of positioner
    """

    hhh = str(int(loc_id)).zfill(3)
    xxx = str(dx).zfill(3)
    yyy = str(dy).zfill(3)
    ttt = str(poslocT).zfill(3)
    ppp = str(poslocP).zfill(3)

    #code = 'PF-' + hhh + '-' + xxx + '-' + yyy + '-' + ttt + '-' + ppp
    code = hhh + '-' + xxx + '-' + yyy + '-' + ttt + '-' + ppp
    return code

def make_table_pf(dx_step=50, dy_step=50, \
                  obsT_start=0, obsT_end=360, obsT_step=5, \
                  posP_start=-20, posP_end=150, posP_step=5):

    """
    dx_step, dy_step in um.
    phi angle ranges from -20 to 150 deg (phi arm fully tucked in >150 deg and so should never
    collides with a fixed neighbor; checked)

    For now pickle save the output table separately.
    File naming convention: table_pf_[dx_step]_[dy_step]_[obsT_step]_[posP_step].out
    """

    # retrieve the location id and center positions of holes near fixed neighbors
    pos_near_fixed = ascii.read('pos_near_fixed.csv')
    near_fixed_id = pos_near_fixed['device_location_id']
    xc_near_fixed = pos_near_fixed['X']
    yc_near_fixed = pos_near_fixed['Y']

    table = dict()
    #case4, case5 = [], []

    start = time.time()

    # loop through all hole id near fixed boundaries
    for i in range(len(near_fixed_id)):
        # "dx" and "dy" are xy variation of phi arm from their nominal locations (x0, y0)
        for dx in range(0, 1000, dx_step):
            for dy in range(0, 1000, dy_step):

                new_x = xc_near_fixed[i]+dx/1000.
                new_y = yc_near_fixed[i]+dy/1000.

                # rotate phi arm by (obsT, posP)
                # translate phi arm to (new_x, new_y)
                for obsT in range(obsT_start, obsT_end, obsT_step):
                    for posP in range(posP_start, posP_end, posP_step):
                        arm = place_phi_arm(obsT, posP, x0=new_x, y0=new_y)

                        # loop through ptl and gfa keepouts and check for collisions
                        for ind in pos.fixed_neighbor_keepouts:
                            fixed_neighbor = pos.fixed_neighbor_keepouts[ind]

                            if arm.collides_with(fixed_neighbor):
                                code = make_code_pf(near_fixed_id[i], dx, dy, obsT, posP)
                                table[code] = ind  # ind either 4 or 5 as ptl = fixed_neighbor_keepouts[5]
                                                   # and gfa = fixed_neighbor_keepout[4]

    end = time.time()
    print((end-start)/60.)
    print(len(table))
    return table

def nearest(x, spacing, min_allowed, max_allowed):
    lower = (x // spacing) * spacing
    if x - lower > spacing / 2:
        val = lower + spacing
    else:
        val = lower
    val = max(val,min_allowed)
    val = min(val,max_allowed)
    return val

############################### VALIDATION SETS (POS-POS) ###################################
def random_collide(n):

    # find possible collision cases by randomly drawing angles and distances
    # n = number of random draws

    collide = []

    pos2_angle = np.arange(0, 360, 60) # possible pos2 angles
    unitX = np.cos(np.pi*pos2_angle/180)
    unitY = np.sin(np.pi*pos2_angle/180)

    for k in range(n):

        # randomly draw theta1, theta2, phi1, and phi2 values within the allowed range
        # for now drawing integers rather than floats.
        obsT1 = np.random.randint(0, 360)
        posP1 = np.random.randint(-20, 210)
        obsT2 = np.random.randint(0, 360)
        posP2 = np.random.randint(-20, 210)

        # randomly draw values of dr, then compute the corresponding (pos2_x, pos2_y)
        # locations at the 6 angles
        dr = np.random.randint(0, 950)
        rnew = dr/1000. + nominal_r
        xnew = rnew*unitX
        ynew = rnew*unitY

        # place pos1 at the randomly-drawn angles
        central1 = place_central_body(obsT1)
        arm1 = place_phi_arm(obsT1, posP1)

        # loop through each pos2 angle
        for i in range(len(xnew)):
            central2 = place_central_body(obsT2, x0=xnew[i], y0=ynew[i])
            arm2 = place_phi_arm(obsT2, posP2, x0=xnew[i], y0=ynew[i])

            # if you use random.uniform() for drawing theta and phi angles,
            # you need to bin these values before checking for collision

            # obsT1 = nearest(obsT1, 1, 0, 359)
            # obsT2 = nearest(obsT2, 1, 0, 359)
            # posP1 = nearest(posP1, 1, -20, 209)
            # posP2 = nearest(posP2, 1, -20, 209)

            # record down the configuration of the collision
            if arm1.collides_with(central2):
                collide.append([pos2_angle[i], dr, obsT1, posP1, obsT2, posP2])

            elif arm2.collides_with(central1):
                collide.append([pos2_angle[i], dr, obsT1, posP1, obsT2, posP2])

            elif arm1.collides_with(arm2):
                collide.append([pos2_angle[i], dr, obsT1, posP1, obsT2, posP2])

    print(len(collide))
    return collide

def create_validation_never_collides(dr_step, theta_step, phi_step, n):

    """
    create validation cases where collision should NEVER happen
    (both phi and phi2 within inner envelopes), by randomly-drawing
    phi1 and phi2 within 145-210 degrees. Theta angles are left free.

    dr_step = hash table ctr-to-ctr resolution
    theta_step = hash table theta angle resolution
    phi_step = hash table phi angle resolution
    n = number of validation case to create
    """

    # set dr_step, theta_step, phi_step to the resolution of the hash table
    # this ensures that the random drawing sample is consistent with the table

    dr = np.arange(0, 1000, dr_step)
    theta = np.arange(0, 360, theta_step)
    phi = np.arange(-20, 210, phi_step)
    angle = np.arange(0, 360, 60) # only 6 positions available for pos2 around pos1

    phi_min, phi_max = 145, 210 # boundaries of phi angle within inner envelope
    ind = np.where((phi >= phi_min) & (phi <= phi_max))[0]
    phi = phi[ind]

    # random draws
    angle_random = angle[np.random.randint(0, len(angle), n)]
    dr_random = dr[np.random.randint(0, len(dr), n)]
    theta1_random = theta[np.random.randint(0, len(theta), n)]
    theta2_random = theta[np.random.randint(0, len(theta), n)]

    # randomly draw phi angles between 145-210 deg
    phi1_random = phi[np.random.randint(0, len(phi), n)]
    phi2_random = phi[np.random.randint(0, len(phi), n)]

    # put the cases in the right format
    nevercollide_sets = np.stack((angle_random, dr_random, theta1_random, phi1_random, theta2_random, phi2_random), axis=1)
    return nevercollide_sets

def create_validation_others():

    # examples of validation sets for collide and no-collide
    # these sets have been vosually inspected.
    # created using random_collide()
    collide_sets = [[0, 845, 329, 36, 181, 195], \
                    [0, 950, 25, 29, 119, 44], \
                    [0, 672, 315, -2, 196, 4], \
                    [60, 118, 61, 36, 266, 2], \
                    [60, 128, 130, 15, 221, 19], \
                    [60, 127, 84, 116, 245, 4], \
                    [120, 204, 122, 172, 258, 74], \
                    [120, 55, 170, 116, 286, 32], \
                    [120, 556, 134, 93, 318, -16], \
                    [180, 565, 157, 17, 48, 169], \
                    [180, 409, 144, 97, 312, 82], \
                    [180, 433, 149, 117, 328, 38], \
                    [180, 180, 155, 9, 65, 201], \
                    [240, 783, 244, -15, 27, 90], \
                    [240, 357, 206, 104, 44, 1], \
                    [240, 268, 258, 155, 30, 33], \
                    [300, 950, 320, 120, 80, 30], \
                    [300, 724, 301, 51, 53, 83], \
                    [300, 546, 323, -10, 112, 15], \
                    # 2-way collision
                    [0, 0, 100, 10, 200, -15], \
                    [60, 0, 100, 10, 200, -15], \
                    # another 2-way collision
                    [60, 19, 61, -15, 301, 10], \
                    [120, 19, 61, -15, 301, 10]]

    # constructed based on subset of collide_sets by changing the angle in each case
    # to other non-colliding angles, as well as some manually inserted cases
    nocollide_sets = [[120, 845, 329, 36, 181, 195], \
                        [60, 950, 25, 29, 119, 44], \
                        [0, 118, 61, 36, 266, 2], \
                        [300, 128, 130, 15, 221, 19], \
                        [180, 127, 84, 116, 245, 4], \
                        [0, 204, 122, 172, 258, 74], \
                        [180, 55, 170, 116, 286, 32], \
                        [0, 565, 157, 17, 48, 169], \
                        [240, 409, 144, 97, 312, 82], \
                        [60, 783, 244, -15, 27, 90], \
                        [0, 357, 206, 104, 44, 1], \
                        [0, 950, 320, 120, 80, 30], \
                        [120, 724, 301, 51, 53, 83], \
                        [180, 0, 100, 10, 200, -15], \
                        [240, 0, 100, 10, 200, -15], \
                        [240, 950, 215, 10, 120, 95], \
                        [60, 350, 130, 15, 221, 19], \
                        [180, 180, 155, 50, 65, 201]]

    return collide_sets, nocollide_sets

def plot_validation(valid_set, ind):

    """
    Visually check that the validation set is "valid".
    In the figure:
    left panel: pos1 always at (0,0) with pos2 configuration duplicated at all 6 positions
    right panel: after removing the angle between pos1 and pos2

    For checking the "nocollide" set, left panel could be confusing, so just look at the
    the right panel.
    For checking the "collide" set, both panels should show a collision, with the left
    panel showing the pos2 angle where the collision happens.

    valid_set = list of list, with [angle, dr, obsT1, posP1, obsT2, posP2].
                Same output format from create_validation...() and same input format
                for validate_table...()

    ind = index in valid_set to plot
    """

    angle = valid_set[ind][0]
    dr = valid_set[ind][1]
    obsT1 = valid_set[ind][2]
    posP1 = valid_set[ind][3]
    obsT2 = valid_set[ind][4]
    posP2 = valid_set[ind][5]
    print(angle, dr, obsT1, posP1, obsT2, posP2)

    plt.figure(figsize=(10, 5))
    plt.subplot(121)
    around_pos1(dr, obsT1, posP1, obsT2, posP2)

    plt.subplot(122)
    obsT1 -= angle
    obsT2 -= angle
    if obsT1 < 0: obsT1 += 360
    if obsT2 < 0: obsT2 += 360

    plot_placements(0, 0, nominal_r + dr/1000., 0, obsT1, posP1, obsT2, posP2)
    #plot_placements(nominal_r + dr/1000., 0, obsT1, posP1, obsT2, posP2)
    plt.tight_layout()

def validate_table(table, cases_to_validate, check_collide, bin_dr, bin_theta, bin_phi):

    """
    1) Go through each case in the validation set
    2) Bin the case to the resolution of the hash table
    3) Assert that the resultant code is found (for check_collide=True) or not found
    (for check_collide=False) in the table.

    This validates the table by checking for the existence/non-existence of a certain
    collision in the table, rather than validating the type of collision.

    table = hash table
    cases_to_validate =  list of list, with [angle, dr, obsT1, posP1, obsT2, posP2]
    check_collide = True or False
    bin_dr = [um] hash table ctr-to-ctr resolution
    bin_theta = [deg] hash table theta angle resolution
    bin_phi = [deg] hash table phi angle resolution
    """

    for case in cases_to_validate:

        angle = case[0]

        # bin the cases to the same resolution as the table and then
        # create the code for the binned case
        dr = nearest(case[1], bin_dr, 0, 950)
        obsT1 = nearest(case[2], bin_theta, 0, 360)
        posP1 = nearest(case[3], bin_phi, -20, 210)
        obsT2 = nearest(case[4], bin_theta, 0, 360)
        posP2 = nearest(case[5], bin_phi, -20, 210)

        new_case = [angle, dr, obsT1, posP1, obsT2, posP2]
        code = make_code_pp(dr, obsT1-angle, posP1, obsT2-angle, posP2)

        # if want to check for collision, assert that code is found
        if check_collide:
            try:
                assert code in table
            except AssertionError:
                print("AssertionError!")
                print("Original case: ", case)
                print("Binned case:   ", new_case)
                print("Resulting code:", code)
                print()

        # if want to check for non-collision, assert that code is NOT found
        else:
            try:
                assert code not in table
            except AssertionError:
                print("AssertionError!")
                print("Original case: ", case)
                print("Binned case:   ", new_case)
                print("Resulting code:", code)
                print()

def validate_table_all(table, nevercollide, collide, nocollide, bin_dr, bin_theta, bin_phi):
    # run all validation sets at once
    # refer to validate_table() for explanation of input arguments

    print("Now validating nevercollide")
    validate_table(table, nevercollide, False, bin_dr, bin_theta, bin_phi)
    print("Now validating collide")
    validate_table(table, collide, True, bin_dr, bin_theta, bin_phi)
    print("Now validating nocollide")
    validate_table(table, nocollide, False, bin_dr, bin_theta, bin_phi)

############################### VALIDATION SETS (POS-FIX) ###################################
def create_validation_nevercollides_pf(n):

    """
    create validation cases for pos-fixed where collision should NEVER happen,
    by random drawing hole ids that are not located directly to a fixed neighbor.

    n = number of random draws
    """

    # hole ids near fixed neighbors
    f = open('ids_near_fixed.dat', 'r')
    near_fixed = []

    for line in f:
        near_fixed.append(int(line.strip('\n')))
    f.close()

    # get the hole ids that are not near fixed neigbors and their xy-centers
    not_near_fixed = list(set(pos_loc_id) - set(near_fixed))
    xc_notnearfixed, yc_notnearfixed = [], []
    for id in not_near_fixed:
        ind = np.where(pos_loc_id == id)[0]
        xc_notnearfixed.append(xc[ind])
        yc_notnearfixed.append(yc[ind])

    # randomly draw n holes from these
    ind_random = np.random.randint(0, len(not_near_fixed), n)
    random_id = np.array(not_near_fixed)[ind_random]
    random_xc = np.array(xc_notnearfixed)[ind_random]
    random_yc = np.array(yc_notnearfixed)[ind_random]

    random_xc = random_xc[:,0]
    random_yc = random_yc[:,0]

    # randomly draw dx and dy values
    # dx and dy are variations from the nominal positions of the positioner
    dx = np.random.randint(0, 950, n)
    dy = np.random.randint(0, 950, n)

    # positions to translate the phi arms
    new_x = random_xc + dx/1000.
    new_y = random_yc + dy/1000.

    # randomly draw theta and phi angles for the phi arm
    obsT = np.random.randint(0, 360, n)
    posP = np.random.randint(-20, 210, n)

    # create the validation set
    nevercollide = []
    for i in range(n):
        nevercollide.append([random_id[i], new_x[i], new_y[i], dx[i], dy[i], obsT[i], posP[i]])

    return nevercollide

def create_validation_collide_pf():

    """
    validation set that collides, set manually; visually inspection done.
    format is [hole-id, new_x (mm), new_y (mm), dx (um), dy (um), obsT, posP]
    """

    collide = [[100, 142.616722, 95.751246, 0, 0, 90, 0], \
               [162, 179.106358, 122.765412, 0, 0, 90, 0], \
               [459, 304.656415, 194.714323, 0, 0, 90, 0], \
               [459, 304.656415, 194.714323, 0, 0, 60, 0], \
               [8, 69.800289, 5.201506, 0, 0, 270, 0], \
               [54, 153.117665, 5.201978, 0, 0, 270, 0], \
               [124, 226.004616, 5.20282, 0, 0, 270, 0], \
               [416, 404.081201, 5.208477, 0, 0, 270, 0]]

    return collide

def plot_validation_pf(valid_set, ind):

    """
    visually inspect a validation set
    """

    new_x = valid_set[ind][1]
    new_y = valid_set[ind][2]
    obsT = valid_set[ind][5]
    posP = valid_set[ind][6]

    arm = place_phi_arm(obsT, posP, x0=new_x, y0=new_y)
    plt.plot(gfa.points[0], gfa.points[1])
    plt.plot(ptl.points[0], ptl.points[1])
    plt.plot(arm.points[0], arm.points[1])


def validate_table_pf(table, cases_to_validate, check_collide, bin_dx, bin_dy, bin_theta, bin_phi):

    """
    same as validate_table() for pos-pos.
    """

    for case in cases_to_validate:
        id = case[0]

        # bin to the same resolution as hash table
        dx = nearest(case[3], bin_dx, 0, 950)
        dy = nearest(case[4], bin_dy, 0, 950)
        obsT = nearest(case[5], bin_theta, 0, 360)
        posP = nearest(case[6], bin_phi, -20, 150)

        new_case = [id, dx, dy, obsT, posP]
        code = make_code_pf(id, dx, dy, obsT, posP) # create code using binned values

        # if want to check for collision, assert that code is found
        if check_collide:
            try:
                assert code in table
            except AssertionError:
                print("AssertionError!")
                print("Original case: ", case)
                print("Binned case:   ", new_case)
                print("Resulting code:", code)
                print()

        # if want to check for non-collision, assert that code is NOT found
        else:
            try:
                assert code not in table
            except AssertionError:
                print("AssertionError!")
                print("Original case: ", case)
                print("Binned case:   ", new_case)
                print("Resulting code:", code)
                print()