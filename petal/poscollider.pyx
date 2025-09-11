# cython: profile=False
# cython: language_level=3
cimport cython
import numpy as np
import posconstants as pc
import posanimator
import configobj
import os
import copy as copymodule
import math

class PosCollider(object):
    """PosCollider contains geometry definitions for mechanical components of the
    fiber positioner, GFA camera, and petal. It provides the methods to check for
    collisions between neighboring positioners, and to check for crossing the
    boundaries of the GFA or petal envelope. It maintains a list of all the particular
    positioners in its petal.

    See DESI-0899 for geometry specifications, illustrations, and kinematics.
    """
    def __init__(self, configfile='',
                 use_neighbor_loc_dict=False,
                 config=None,
                 animator_label_type='both',
                 printfunc=print):
        self.printfunc = printfunc
        if not config:
            if not configfile:
                filename = pc.default_collider_filename
            else:
                filename = configfile
            filepath = os.path.join(pc.dirs['collision_settings'], filename)
            self.config = configobj.ConfigObj(filepath, unrepr=True)
        else:
            self.config = config
        self.posids = set() # posid strings for all the positioners
        self.posindexes = {} # key: posid string, value: index number for positioners in animations
        self.posmodels = {} # key: posid string, value: posmodel instance
        self.devicelocs = {} # key: device loc, value: posid string
        self.pos_neighbors = {} # all the positioners that surround a given positioner. key is a posid, value is a set of neighbor posids
        self.fixed_neighbor_cases = {} # all the fixed neighbors that apply to a given positioner. key is a posid, value is a set of the fixed neighbor cases
        self.R1, self.R2, self.x0, self.y0, self.t0, self.p0 = {}, {}, {}, {}, {}, {}  # keys are posids
        self.keepout_expansions = {} # keys are posids, values are subdicts, containing the 4x keepout expansion params for each positioner
        self.plotting_on = True
        self.timestep = self.config['TIMESTEP']
        self.animator = posanimator.PosAnimator(fignum=0, timestep=self.timestep)
        self.animate_colliding_only = False # overrides posids_to_animate and fixed_items_to_animate
        self.posids_to_animate = set() # set of posids to be animated (i.e., allows restricting which ones get plotted)
        self.fixed_items_to_animate = {'PTL','GFA'} # may contain nothing, 'PTL', and/or 'GFA'
        self.labeled_posids = set() # tracks those already labeled in animator
        self.animator_label_type = animator_label_type # 'loc' --> device location ids, 'posid' --> posids, None --> no labels
        self.use_neighbor_loc_dict = use_neighbor_loc_dict
        self.keepouts_T = {} # key: posid, value: central body keepout of type PosPoly
        self.keepouts_P = {} # key: posid, value: phi arm keepout of type PosPoly
        self.keepouts_arcP = {} # key: posid, value: phi arm keepout swept through its full range, of type PosPoly
        self.keepouts_arcP_resolution = 25 # number of points to add when generating polygonal full-range phi arc
        self.classified_as_retracted = set() # posids of robots classified as retracted. overrides polygonal keepout calcs

        # load fixed dictionary containing locations of neighbors for each positioner DEVICE_LOC (if this option has been selected)
        if self.use_neighbor_loc_dict:
            self.neighbor_locs = pc.generic_pos_neighbor_locs

    def refresh_calibrations(self, verbose=True):
        """Reloads positioner parameters.
        """
        self._load_config_data(verbose=verbose)

    def add_positioners(self, posmodels, verbose=True):
        """Add a collection of positioners to the collider object.
        """
        for posmodel in posmodels:
            posid = posmodel.posid
            if posid not in self.posids:
                self.posids.add(posid)
                self.posindexes[posid] = len(self.posindexes) + 1
                self.posmodels[posid] = posmodel
                self.devicelocs[posmodel.deviceloc] = posid
                self.pos_neighbors[posid] = set()
                self.fixed_neighbor_cases[posid] = set()
        self._load_config_data(verbose=verbose)
        for p in self.posids:
            self._identify_neighbors(p)
        self.posids_to_animate.update(self.posids) # this doesn't turn the animator on --- just gathering the set in case

    def add_fixed_to_animator(self, start_time=0):
        """Add unmoving polygon shapes to the animator.

            start_time ... seconds, global time when the move begins
        """
        if not self.animate_colliding_only:
            if 'GFA' in self.fixed_items_to_animate:
                self.animator.add_or_change_item('GFA', '', start_time, self.keepout_GFA.points)
            if 'PTL' in self.fixed_items_to_animate:
                self.animator.add_or_change_item('PTL', '', start_time, self.keepout_PTL.points)
            for posid in self.posids_to_animate:
                Eo_style_override = 'Eo bold' if posid in self.classified_as_retracted else ''
                self.animator.add_or_change_item('Eo', self.posindexes[posid], start_time, self.Eo_polys[posid].points, Eo_style_override)
                # self.animator.add_or_change_item('Ei', self.posindexes[posid], start_time, self.Ei_polys[posid].points)
                # self.animator.add_or_change_item('Ee', self.posindexes[posid], start_time, self.Ee_polys[posid].points)
                self.animator.add_or_change_item('line t0', self.posindexes[posid], start_time, self.line_t0_polys[posid].points)
                self.add_posid_label(posid)
            
    def add_posid_label(self, posid):
        """Adds label to animator for posid (if not already there)."""
        if self.animator_label_type == 'posid':
            label = str(posid)
        elif self.animator_label_type == 'loc':
            label = format(self.posmodels[posid].deviceloc,'03d')
        elif self.animator_label_type == 'both':
            label = f'{posid}\n{self.posmodels[posid].deviceloc:03d}'
        else:
            label = ''
        if posid not in self.labeled_posids:
            self.animator.add_label(label, self.x0[posid], self.y0[posid])
            self.labeled_posids.add(posid)

    def add_mobile_to_animator(self, start_time, sweeps):
        """Add a collection of PosSweeps to the animator, describing positioners'
        real-time motions.

            start_time ... seconds, global time when the move begins
            sweeps     ... dict with keys = posids, values = PosSweep instances
        """
        for posid,s in sweeps.items():
            if posid in self.posids_to_animate or self.animate_colliding_only:
                posidx = self.posindexes[posid]
                for i in range(len(s.time)):
                    style_override = ''
                    collision_has_occurred = s.time[i] >= s.collision_time
                    freezing_has_occurred = s.time[i] >= s.frozen_time
                    if posid in self.classified_as_retracted or not self.posmodels[posid].is_enabled:
                        style_override = 'positioner element unbold'
                    if freezing_has_occurred:
                        style_override = 'frozen'
                    if collision_has_occurred:
                        style_override = 'collision'
                    time = start_time + s.time[i]
                    self.animator.add_or_change_item('central body', posidx, time, self.place_central_body(posid, s.tp[i][0]).points, style_override)
                    self.animator.add_or_change_item('phi arm',      posidx, time, self.place_phi_arm(     posid, s.tp[i]).points, style_override)
                    self.animator.add_or_change_item('ferrule',      posidx, time, self.place_ferrule(     posid, s.tp[i]).points, style_override)
                    if collision_has_occurred and s.collision_case == pc.case.GFA:
                        self.animator.add_or_change_item('GFA', '', time, self.keepout_GFA.points, style_override)
                    elif collision_has_occurred and s.collision_case == pc.case.PTL:
                        self.animator.add_or_change_item('PTL', '', time, self.keepout_PTL.points, style_override)

    def spacetime_collision_between_positioners(self, posid_A, init_poslocTP_A, tableA,
                                                      posid_B, init_poslocTP_B, tableB,
                                                      skip=0):
        """Wrapper for spacetime_collision method, specifically for checking
        two positioners against each other."""
        return self.spacetime_collision(posid_A, init_poslocTP_A, tableA,
                                        posid_B, init_poslocTP_B, tableB,
                                        skip=skip)

    def spacetime_collision_with_fixed(self, posid, init_poslocTP, table, skip=0):
        """Wrapper for spacetime_collision method, specifically for checking
        one positioner against the fixed keepouts.
        """
        return self.spacetime_collision(posid, init_poslocTP, table, skip=skip)

    def spacetime_collision(self, posid_A, init_poslocTP_A, tableA,
                                  posid_B=None, init_poslocTP_B=None, tableB=None,
                                  skip=0):
        """Searches for collisions in time and space between two positioners
        which are rotating according to the argued tables.

            posid_A, posid_B ... posid strings of the two positioners to check against each other
            init_poslocTP_A, init_poslocTP_B  ... starting (theta,phi) positions, in the poslocTP coordinate systems
            tableA, table  ... dictionaries defining rotation schedules as described below
            skip ... integer number of initial timesteps for which to skip collision check
        
        ** Subtlety regarding "skip" **
            As of 2020-10-29, I *always* skip the first timestep, regardless of the
            value of skip. In effect, arguing skip=0 is like saying skip=1. This is
            due to some details of how I most efficiently process the quantized sweeps
            (in particular, the "was_moving" flags). So far in practice, I believe this
            is fine --- but it may not be what you expect at first glance. You can
            still use "skip" to skip more than one timestep, although at present with
            our timestep values of 0.02 sec there's no use case I can think of to do so.

        If no arguments are provided for the "B" positioner (i.e. no args for idxB, init_poslocTP_B, tableB)
        then the method checks the "A" positioner against the fixed keepout envelopes.

        The input table dictionaries must contain the following fields:

            'nrows'     : number of rows in the lists below (all must be the same length)
            'dT'        : list of theta rotation distances in degrees
            'dP'        : list of phi rotation distances in degrees
            'prepause'  : list of prepause (before the rotations begin) values in seconds
            'move_time' : list of durations of rotations in seconds, approximately equals max(dT/Tdot,dP/Pdot), but calculated more exactly for the physical hardware
            'postpause' : list of postpause (after the rotations end) values in seconds

        The return is a list of instances of PosSweep. (These contain the theta and phi rotations
        in real time, when if any collision, and the collision type and neighbor.)
        """
        pospos = posid_B is not None
        if pospos:
            init_poslocTPs = [init_poslocTP_A, init_poslocTP_B]
            tables = [tableA,tableB]
            sweeps = [PosSweep(posid_A),PosSweep(posid_B)]
            steps_remaining = [0,0]
            step = [0,0]
            pos_range = [0,1]
        else:
            init_poslocTPs = [init_poslocTP_A]
            tables = [tableA]
            sweeps = [PosSweep(posid_A)]
            steps_remaining = [0]
            step = [0]
            pos_range = [0]
        rev_pos_range = pos_range[::-1]
        for i in pos_range:
            sweeps[i].fill_exact(init_poslocTPs[i], tables[i])
            sweeps[i].quantize(self.timestep)
            steps_remaining[i] = len(sweeps[i].time)
        while any(steps_remaining):
            check_collision_this_loop = False
            for i in pos_range:
                if sweeps[i].was_moving_cached[step[i]] and step[i] >= skip:
                    check_collision_this_loop = True
            if check_collision_this_loop:
                if pospos:
                    collision_case = self.spatial_collision_between_positioners(posid_A, posid_B, sweeps[0].tp[step[0]], sweeps[1].tp[step[1]])
                else:
                    collision_case = self.spatial_collision_with_fixed(posid_A, sweeps[0].tp[step[0]])
                if collision_case != pc.case.I:
                    for i, j in zip(pos_range, rev_pos_range):
                        sweeps[i].collision_case = collision_case
                        if pospos:
                            sweeps[i].collision_neighbor = sweeps[j].posid
                            possible_collision_times = {sweeps[i].time[step[i]]:step[i], sweeps[j].time[step[j]]:step[j]} # key: time value, value: step index
                            sweeps[i].collision_time = max(possible_collision_times)
                            sweeps[i].collision_idx = possible_collision_times[sweeps[i].collision_time]
                        else:
                            sweeps[i].collision_neighbor = 'PTL' if collision_case == pc.case.PTL else 'GFA'
                            sweeps[i].collision_time = sweeps[i].time[step[i]]
                            sweeps[i].collision_idx = step[i]
                        steps_remaining[i] = 0 # halt the sweep here
            steps_remaining = [max(step-1,0) for step in steps_remaining]
            for i in pos_range:
                if steps_remaining[i]:
                    step[i] += 1
                else:
                    pass
        return sweeps

    def spatial_collision_between_positioners(self, posid_A, posid_B, poslocTP_A, poslocTP_B):
        """Searches for collisions in space between two fiber positioners.

            posid_A, posid_B  ...  posid strings of the two positioners to check against each other
            poslocTP_A, poslocTP_B  ...  (theta,phi) positions of the axes for the two positioners

        poslocTP_A and poslocTP_B are in the (poslocT,poslocP) coordinate system, as defined in
        PosTransforms.

        The return is an enumeration of type "case", indicating what kind of collision
        was first detected, if any.
        """
        A_is_within_Eo = poslocTP_A[1] >= self.Eo_phi or posid_A in self.classified_as_retracted 
        B_is_within_Eo = poslocTP_B[1] >= self.Eo_phi or posid_B in self.classified_as_retracted
        if A_is_within_Eo and B_is_within_Eo:
            return pc.case.I
        elif not A_is_within_Eo and posid_B in self.classified_as_retracted: # check case IV, A upon B
            if self._case_IV_collision(posid_A, posid_B, poslocTP_A):
                return pc.case.IV
            else:
                return pc.case.I
        elif not B_is_within_Eo and posid_A in self.classified_as_retracted: # check case IV, B upon B
            if self._case_IV_collision(posid_B, posid_A, poslocTP_B):
                return pc.case.IV
            else:
                return pc.case.I
        elif poslocTP_A[1] < self.Eo_phi and poslocTP_B[1] >= self.Ei_phi: # check case III, A upon B
            if self._case_III_collision(posid_A, posid_B, poslocTP_A, poslocTP_B[0]):
                return pc.case.III
            else:
                return pc.case.I
        elif poslocTP_B[1] < self.Eo_phi and poslocTP_A[1] >= self.Ei_phi: # check case III, B upon A
            if self._case_III_collision(posid_B, posid_A, poslocTP_B, poslocTP_A[0]):
                return pc.case.III
            else:
                return pc.case.I
        else: # check cases II and III
            if self._case_III_collision(posid_A, posid_B, poslocTP_A, poslocTP_B[0]):
                return pc.case.III
            elif self._case_III_collision(posid_B, posid_A, poslocTP_B, poslocTP_A[0]):
                return pc.case.III
            elif self._case_II_collision(posid_A, posid_B, poslocTP_A, poslocTP_B):
                return pc.case.II
            else:
                return pc.case.I

    def spatial_collision_with_fixed(self, posid, poslocTP, use_phi_arc=False):
        """Searches for collisions in space between a fiber positioner and all
        fixed keepout envelopes.

            posid ... positioner to check
            poslocTP ... (theta,phi) position of the axes of the positioner
            use_phi_arc ... boolean, uses phi full-range arc instead of typical arm keepout

        poslocTP is in the (poslocT,poslocP) coordinate system, as defined in
        PosTransforms.

        The return is an enumeration of type "case", indicating what kind of collision
        was first detected, if any.
        """
        cdef PosPoly poly1
        cdef PosPoly poly2
        if self.fixed_neighbor_cases[posid]:
            if use_phi_arc:
                poly1 = self.place_phi_arc(posid, poslocTP[0])
            else:
                poly1 = self.place_phi_arm(posid, poslocTP)
            for fixed_case in self.fixed_neighbor_cases[posid]:
                poly2 = self.fixed_neighbor_keepouts[fixed_case]
                if poly1.collides_with(poly2):
                    return fixed_case
        return pc.case.I
    
    def phi_range_collision(self, posid_A, poslocT_A, posid_B, poslocTP_B=None):
        """Special function to search for collisions in space between a positioner
        and its neighbors for all possible positions of phi. Depending on the arguments,
        this function can also check the full-range phi against neighbors' full-range
        phi, or against fixed boundaries.
        
        The intended use case is for identifying when it is safe to move an
        unpredictable phi arm outside the Eo retracted envelope.
        
            posid_A ... positioner with fixed full-range phi arc
            poslocT_A ... theta position of A
            posid_B ... neighbor of A, which may be either a posid or the string
                        'fixed' (in which case check against any fixed boundaries)
            poslocTP_B ... position of B, where:
                           ... (theta, phi) pair --> use neighbor's normal phi arm keepout
                           ... scalar value --> treat this scalar as the neighbor's theta, and
                                                use full-range phi arc for neighbor's phi keepout
                           ... arg is ignored if posid_B == 'fixed'
                           
        The return is an enumeration of type "case", indicating what kind of collision
        was first detected, if any.
        """
        cdef PosPoly poly1
        cdef PosPoly poly2
        if posid_B == 'fixed':
            dummyP = 0.0
            return self.spatial_collision_with_fixed(posid_A, [poslocT_A, dummyP], use_phi_arc=True)
        poly1 = self.place_phi_arc(posid_A, poslocT_A)
        if posid_B in self.classified_as_retracted:
            if poly1.collides_with_circle(self.x0[posid_B], self.y0[posid_B], self.Eo_radius_with_margin):
                return pc.case.IV
        use_neighbor_arc = not isinstance(poslocTP_B, (list, tuple))
        neighbor_place_phi = self.place_phi_arc if use_neighbor_arc else self.place_phi_arm
        poly2 = neighbor_place_phi(posid_B, poslocTP_B)
        if poly1.collides_with(poly2):
            return pc.case.II
        poslocT_B = poslocTP_B if use_neighbor_arc else poslocTP_B[0]
        poly2 = self.place_central_body(posid_B, poslocT_B)
        if poly1.collides_with(poly2):
            return pc.case.III
        return pc.case.I

    def place_phi_arm(self, posid, poslocTP):
        """Rotates and translates the phi arm to position defined by the positioner's
        (x0,y0) and the argued poslocTP (theta,phi) angles.
        """
        return self.keepouts_P[posid].place_as_phi_arm(theta=poslocTP[0],
                                                       phi=poslocTP[1],
                                                       x0=self.x0[posid],
                                                       y0=self.y0[posid],
                                                       r1=self.R1[posid])

    def place_phi_arc(self, posid, poslocT):
        """Rotates and translates the full-range phi arc to position defined by the
        positioner's (x0,y0) and the argued poslocT (theta) angle.
        
        This function is like place_phi_arm(), but rather than the typical phi arm
        keepout, it uses a larger swept-arc polygon. That polygon represents the
        total area the phi arm can possibly inhabit, in its full range of motion,
        when theta is held constant at the argued value. 
        """
        return self.keepouts_arcP[posid].rotated(poslocT).translated(self.x0[posid], self.y0[posid])

    def place_central_body(self, posid, poslocT):
        """Rotates and translates the central body of positioner
        to its (x0,y0) and the argued poslocT theta angle.
        """
        return self.keepouts_T[posid].place_as_central_body(theta=poslocT,
                                                            x0=self.x0[posid],
                                                            y0=self.y0[posid])
    
    def place_arm_lines(self, posid, poslocTP):
        '''Generates, rotates, and translates a PosPoly with two lines through the
        central body and phi arms. These graphically indicate the theta dnd phi
        angles more precisely, when making plots.
        '''
        cdef PosPoly poly_t
        cdef PosPoly poly_p
        cdef PosPoly combo
        t_line_x = [0, self.R1[posid]]
        p_line_x = [0, self.R2[posid]]
        line_y = [0,0]
        poly_t = PosPoly([t_line_x, line_y], close_polygon=False)
        poly_p = PosPoly([p_line_x, line_y], close_polygon=False)
        poly_t = poly_t.place_as_central_body(theta=poslocTP[0], x0=self.x0[posid], y0=self.y0[posid])
        poly_p = poly_p.place_as_phi_arm(theta=poslocTP[0], phi=poslocTP[1], x0=self.x0[posid], y0=self.y0[posid], r1=self.R1[posid])
        t_points = poly_t.points
        p_points = poly_p.points
        combo_x_pts = t_points[0] + p_points[0]
        combo_y_pts = t_points[1] + p_points[1]
        combo = PosPoly([combo_x_pts, combo_y_pts], close_polygon=False)
        return combo

    def place_ferrule(self, posid, poslocTP):
        """Rotates and translates the ferrule to position defined by the positioner's
        (x0,y0) and the argued poslocTP (theta,phi) angles.
        """
        cdef PosPoly poly
        poly = self.ferrule_poly.translated(self.R2[posid], 0)
        poly = poly.rotated(poslocTP[1])
        poly = poly.translated(self.R1[posid],0)
        poly = poly.rotated(poslocTP[0])
        poly = poly.translated(self.x0[posid], self.y0[posid])
        return poly

    def _case_II_collision(self, posid1, posid2, tp1, tp2):
        """Search for case II collision, positioner 1 arm against positioner 2 arm."""
        cdef PosPoly poly1
        cdef PosPoly poly2
        poly1 = self.place_phi_arm(posid1, tp1)
        poly2 = self.place_phi_arm(posid2, tp2)
        return poly1.collides_with(poly2)

    def _case_III_collision(self, posid1, posid2, tp1, t2):
        """Search for case III collision, positioner 1 arm against positioner 2 central body."""
        cdef PosPoly poly1
        cdef PosPoly poly2
        poly1 = self.place_phi_arm(posid1, tp1)
        poly2 = self.place_central_body(posid2, t2)
        return poly1.collides_with(poly2)
    
    def _case_IV_collision(self, posid1, posid2, tp1):
        """Search or case IV collision, positioner 1 arm against positioner 2 circular envelope."""
        cdef PosPoly poly1
        poly1 = self.place_phi_arm(posid1, tp1)
        return poly1.collides_with_circle(self.x0[posid2], self.y0[posid2], self.Eo_radius_with_margin)

    def _load_config_data(self, verbose=True):
        """Reads latest versions of all configuration file offset, range, and
        polygon definitions, and updates stored values accordingly.
        """
        self.timestep = self.config['TIMESTEP']
        self._load_positioner_params(verbose=verbose)
        self._load_keepouts()
        self._adjust_keepouts()
        self._load_circle_envelopes()
        self._load_keepouts_arcP()

    def _load_positioner_params(self, verbose=True):
        """Read latest versions of all positioner parameters."""
        for posid, posmodel in self.posmodels.items():
            self.R1[posid] = posmodel.state.read('LENGTH_R1')
            self.R2[posid] = posmodel.state.read('LENGTH_R2')
            self.x0[posid] = posmodel.state.read('OFFSET_X')
            self.y0[posid] = posmodel.state.read('OFFSET_Y')
            self.t0[posid] = posmodel.state.read('OFFSET_T')
            self.p0[posid] = posmodel.state.read('OFFSET_P')
            self.keepout_expansions[posid] = {key:posmodel.state.read(key) for key in pc.keepout_expansion_keys}
            if posmodel.is_linphi:
                angP = self.keepout_expansions[posid].get('KEEPOUT_EXPANSION_PHI_ANGULAR', 0.0)
                self.keepout_expansions[posid]['KEEPOUT_EXPANSION_PHI_ANGULAR'] = max(pc.P_zeno_jog, angP)
            classified_retracted = posmodel.state.read('CLASSIFIED_AS_RETRACTED')
            disabled = not posmodel.state.read('CTRL_ENABLED')
            if classified_retracted:
                self.classified_as_retracted.add(posid)

    def _load_keepouts(self):
        """Read latest versions of all keepout geometries."""
        self.general_keepout_P = PosPoly(self.config['KEEPOUT_PHI'])
        self.general_keepout_T = PosPoly(self.config['KEEPOUT_THETA'])
        self.keepout_PTL = PosPoly(self.config['KEEPOUT_PTL'])
        self.keepout_GFA = PosPoly(self.config['KEEPOUT_GFA'])
        self.fixed_neighbor_keepouts = {pc.case.PTL : self.keepout_PTL, pc.case.GFA : self.keepout_GFA}

    def _adjust_keepouts(self):
        """Expand/contract, and pre-shift the theta and phi keepouts for each positioner."""
        R1_nom = pc.nominals['LENGTH_R1']['value']
        R2_nom = pc.nominals['LENGTH_R2']['value']
        for posid in self.posids:
            R1_error = self.R1[posid] - R1_nom  # true R1 err desired, since it is kinematically real
            R2_error = max(self.R2[posid] - R2_nom, 0.0)  # only expand phi (not contract) it, since this is just distance to fiber, and contraction might not represent true mechanical shape
            expansions = self.keepout_expansions[posid]
            keepout_P = self.general_keepout_P.translated(0,0) # effectively just a copy operation
            keepout_T = self.general_keepout_T.translated(0,0) # effectively just a copy operation
            keepout_P = keepout_P.expanded_radially(expansions['KEEPOUT_EXPANSION_PHI_RADIAL'])
            keepout_P = keepout_P.expanded_angularly(expansions['KEEPOUT_EXPANSION_PHI_ANGULAR'])
            keepout_P = keepout_P.translated(R1_error,0)
            keepout_P = keepout_P.expanded_x(left_shift=R1_error, right_shift=R2_error)
            keepout_T = keepout_T.expanded_radially(expansions['KEEPOUT_EXPANSION_THETA_RADIAL'])
            keepout_T = keepout_T.expanded_angularly(expansions['KEEPOUT_EXPANSION_THETA_ANGULAR'])
            self.keepouts_P[posid] = keepout_P
            self.keepouts_T[posid] = keepout_T

    def _load_circle_envelopes(self):
        """Read latest versions of all circular envelopes, including outer clear rotation
        envelope (Eo), inner clear rotation envelope (Ei) and extended-phi clear rotation
        envelope (Ee).
        """
        self.Eo_phi = self.config['PHI_EO']   # poslocP angle above which phi is guaranteed to be within envelope Eo
        self.Ei_phi = self.config['PHI_EI']   # poslocP angle above which phi is guaranteed to be within envelope Ei
        self.Eo = self.config['ENVELOPE_EO']  # outer clear rotation envelope
        self.Eo_with_margin = self.Eo + 2 * self.config['EO_RADIAL_TOL'] # outer clear rotation envelope for collision checks (diameter)
        self.Eo_radius_with_margin = self.Eo_with_margin / 2 # outer clear rotation envelope for collision checks (radius)
        self.Ei = self.config['ENVELOPE_EI']  # inner clear rotation envelope
        self.Ee = self._max_extent() * 2      # extended-phi clear rotation envelope
        self.Eo_poly = PosPoly(self._circle_poly_points(self.Eo, self.config['RESOLUTION_EO']))
        self.Eo_poly_with_margin = PosPoly(self._circle_poly_points(self.Eo_with_margin, self.config['RESOLUTION_EO']))
        self.Ei_poly = PosPoly(self._circle_poly_points(self.Ei, self.config['RESOLUTION_EI']))
        self.Ee_poly = PosPoly(self._circle_poly_points(self.Ee, self.config['RESOLUTION_EE']))
        line_x = [0, self.Eo/2]
        line_y = [0, 0]
        self.line_t0_poly = PosPoly([line_x, line_y], close_polygon=False)
        self.Eo_polys = {}
        self.Ei_polys = {}
        self.Ee_polys = {}
        self.line_t0_polys = {}
        for posid in self.posids:
            x = self.x0[posid]
            y = self.y0[posid]
            if posid in self.classified_as_retracted:
                Eo_poly = self.Eo_poly_with_margin
            else:
                Eo_poly = self.Eo_poly
            self.Eo_polys[posid] = Eo_poly.translated(x,y)
            self.Ei_polys[posid] = self.Ei_poly.translated(x,y)
            self.Ee_polys[posid] = self.Ee_poly.translated(x,y)
            self.line_t0_polys[posid] = self.line_t0_poly.rotated(self.t0[posid]).translated(x,y)
        self.ferrule_diam = self.config['FERRULE_DIAM']
        self.ferrule_poly = PosPoly(self._circle_poly_points(self.ferrule_diam, self.config['FERRULE_RESLN']))

    def _load_keepouts_arcP(self):
        '''Generate latest full-range phi arc keepouts.'''
        dummy_T = 0
        gpts = self.general_keepout_P.points
        gpts = [gpts[0][:-1], gpts[1][:-1]] # take off the last polygon closure point
        idx0 = gpts[1].index(0.0) # specific to the phi polygon definition as of 2021-06-21
        remaining_idxs = list(range(idx0 + 1, len(gpts[0]))) + list(range(idx0))
        for posid, posmodel in self.posmodels.items():
            range_posintP = posmodel.full_range_posintP
            angular_range = max(range_posintP) - min(range_posintP)
            expanded = self.keepouts_P[posid].expanded_angularly(angular_range/2)
            pts = expanded.points
            arc_radius = math.hypot(pts[0][idx0], pts[1][idx0])
            n = self.keepouts_arcP_resolution
            dA = angular_range / n
            arc_angles = [i*dA - angular_range/2 for i in range(n+1)]
            arc_x = [arc_radius*math.cos(math.radians(a)) for a in arc_angles]
            arc_y = [arc_radius*math.sin(math.radians(a)) for a in arc_angles]
            old_x = [pts[0][i] for i in remaining_idxs]
            old_y = [pts[1][i] for i in remaining_idxs]
            new = PosPoly([arc_x + old_x, arc_y + old_y])
            center_posintP = sum(range_posintP)/2
            center_poslocP = posmodel.trans.posintTP_to_poslocTP([dummy_T, center_posintP])[1]
            rotated = new.rotated(center_poslocP)
            self.keepouts_arcP[posid] = rotated.translated(self.R1[posid], 0)

    def _identify_neighbors(self, posid):
        """Find all neighbors which can possibly collide with a given positioner."""
        Ee = self.Ee_poly.translated(self.x0[posid], self.y0[posid])
        if self.use_neighbor_loc_dict:
            deviceloc = self.posmodels[posid].deviceloc
            neighbors = self.neighbor_locs[deviceloc].intersection(self.devicelocs.keys())
            self.pos_neighbors[posid] = {self.devicelocs[loc] for loc in neighbors}
        else:
            for possible_neighbor in self.posids:
                Ee_neighbor = self.Ee_poly.translated(self.x0[possible_neighbor], self.y0[possible_neighbor])
                if not(posid == possible_neighbor) and Ee.collides_with(Ee_neighbor):
                    self.pos_neighbors[posid].add(possible_neighbor)
        for possible_neighbor in self.fixed_neighbor_keepouts:
            EE_neighbor = self.fixed_neighbor_keepouts[possible_neighbor]
            if Ee.collides_with(EE_neighbor):
                self.fixed_neighbor_cases[posid].add(possible_neighbor)
        assert len(self.pos_neighbors[posid]) <= 6, f'{posid}: num neighbors > 6 is geometrically invalid. This indicates a problem ' + \
                                                     ' with calibration values or polygon geometry. Must be fixed before proceeding.'

    def _max_extent(self):
        """Calculation of max radius of keepout for a positioner with fully-extended phi arm."""
        extended_phi = self.general_keepout_P.translated(max(self.R1.values()),0) # assumption here that phi arm polygon defined at 0 deg angle
        return max(np.sqrt(np.sum(np.array(extended_phi.points)**2, axis=0)))

    @staticmethod
    def _circle_poly_points(diameter, npts, outside=True):
        """Constuct a polygon approximating a circle of a given diameter, with a
        given number of points. The polygon's segments are tangent to the circle if
        the optional argument 'outside' is true. If 'outside' is false, then the
        segment points lie on the circle.
        """
        assert diameter > 0
        assert npts > 2
        alpha = np.linspace(0, 2 * math.pi, npts + 1)[0:-1]
        if outside:
            half_angle = (alpha[1] - alpha[0])/2
            points_radius = diameter/2 / math.cos(half_angle)
        else:
            points_radius = diameter/2
        x = points_radius * np.cos(alpha)
        y = points_radius * np.sin(alpha)
        return [x,y]


class PosSweep(object):
    """Contains a real-time description of the sweep of positioner mechanical
    geometries through space.
    """
    def __init__(self, posid=None):
        self.posid = posid                  # unique posid string of the positioner
        self.time = []                      # time at which each TP position value occurs
        self.tp = []                        # theta,phi angles (poslocTP coordinates) as function of time (sign indicates direction)
        self.was_moving_cached = []         # cached boolean values, corresponding to timesteps, as computed by "was_moving() method
        self.collision_case = pc.case.I     # enumeration of type "case", indicating what kind of collision first detected, if any
        self.collision_time = math.inf      # time at which collision occurs. if no collision, the time is inf
        self.collision_idx = None           # index in time and theta,phi lists at which collision occurs
        self.collision_neighbor = ''        # id string (posid, 'PTL', or 'GFA') of neighbor it collides with, if any
        self.frozen_time = math.inf         # time at which positioner is frozen in place. if no freezing, the time is inf

    def copy(self):
        return copymodule.deepcopy(self)
    
    def as_dict(self):
        """Returns a dictionary containing copies of all the sweep data."""
        c = self.copy()
        d = {'posid':               c.posid,
             'time':                c.time,
             'tp':                  c.tp,
             'was_moving':          c.was_moving_cached,
             'collision_case':      c.collision_case,
             'collision_time':      c.collision_time,
             'collision_idx':       c.collision_idx,
             'collision_neighbor':  c.collision_neighbor,
             'frozen_time':         c.frozen_time}
        return d
    
    def __repr__(self):
        d = self.as_dict()
        d['time_start'] = d['time'][0] if d['time'] else None
        d['time_final'] = d['time'][-1] if d['time'] else None
        d['tp_start'] = d['tp'][0] if d['tp'] else None
        d['tp_final'] = d['tp'][-1] if d['tp'] else None
        del d['time']
        del d['tp']
        del d['was_moving']
        return str(d)

    def __str__(self):
        return self.__repr__()
    
    def fill_exact(self, init_poslocTP, table, start_time=0):
        """Fills in a sweep object based on the input table. Time and position
        are handled continuously and exactly (i.e. not yet quantized).
        """
        time = [start_time]
        tp = [list(init_poslocTP)]
        for i in range(table['nrows']):
            if table['prepause'][i]:
                time.append(table['prepause'][i] + time[-1])
                tp.append(tp[-1].copy())
            if table['move_time'][i]:
                time.append(table['move_time'][i] + time[-1])
                tp.append([table['dT'][i] + tp[-1][0],
                           table['dP'][i] + tp[-1][1]])
            if table['postpause'][i]:
                time.append(table['postpause'][i] + time[-1])
                tp.append(tp[-1].copy())
        self.time = time
        self.tp = tp
        self.was_moving_cached = [self.was_moving(step) for step in range(len(self.time))]
        
    def quantize(self, timestep):
        """Converts itself from exact, continuous time to quantized, discrete time.
        The result has approximate intermediate (theta,phi) positions and speeds,
        all as a function of discrete time. The quantization is according to the
        parameter 'timestep'.
        """
        axes = [pc.T,pc.P]
        discrete_time = [self.time[0]]
        discrete_position = [self.tp[0].copy()]
        for i in range(1,len(self.time)):
            time_diff = self.time[i] - discrete_time[-1]
            n_steps = int(time_diff / timestep)
            distance_diffs = [self.tp[i][ax] - self.tp[i-1][ax] for ax in axes]
            if n_steps == 0 and any(distance_diffs):
                n_steps = 1
            if n_steps:                
                distance_steps = [distance_diffs[ax] / n_steps for ax in axes]
                for j in range(n_steps):
                    discrete_time.append(discrete_time[-1] + timestep)
                    discrete_position.append([discrete_position[-1][ax] + distance_steps[ax] for ax in axes])
        self.time = discrete_time
        self.tp = discrete_position
        self.was_moving_cached = [self.was_moving(step) for step in range(len(self.time))]

    def extend(self, timestep, max_time):
        """Extends a sweep object to max_time to reflect the postpauses inserted into the move table
        in equalize_table_times() in posschedulestage.py, ensuring that the sweep object
        is in sync with the move table so that the animator is reflecting true moves. """

        starttime_extension = self.time[-1] + timestep
        n_steps = int((max_time + timestep)/timestep)
        time_extension = [starttime_extension + k * timestep for k in range(n_steps)]
        self.time += time_extension
        tp_extension = [self.tp[-1].copy() for k in range(n_steps)]
        self.tp += tp_extension
        self.was_moving_cached += [self.was_moving(step) for step in range(len(self.was_moving_cached),len(self.time))]

    def register_as_frozen(self):
        """Sets an indicator that the sweep has been frozen at the end."""
        self.frozen_time = self.time[-1]
        
    def clear_collision(self):
        """Resets collision state to default (non-colliding) values."""
        blank_sample = PosSweep()
        self.collision_case = blank_sample.collision_case
        self.collision_time = blank_sample.collision_time
        self.collision_idx = blank_sample.collision_idx
        self.collision_neighbor = blank_sample.collision_neighbor

    @property
    def is_frozen(self):
        """Returns boolean value whether the sweep has a "freezing" event."""
        return self.frozen_time < math.inf

    def was_moving(self, step):
        """Returns boolean value whether the sweep is moving in the most recent
        timestep. 'Most recent' here means the period from (step-1) to step.
        By definition, when step == 0 this function returns False."""
        if step <= 0 or step >= len(self.tp):
            return False
        if self.tp[step][0] == self.tp[step-1][0] and self.tp[step][1] == self.tp[step-1][1]:
            return False
        return True
    
    def axis_was_moving(self, step, axis):
        '''Like was_moving, but for axis = 0 (means theta) or 1 (means phi).
        Separately defined from was_moving, to avoid some overhead when that
        function is repeatedly called during anticollision calcs.'''
        if step <= 0 or step >= len(self.tp):
            return False
        if self.tp[step][axis] == self.tp[step-1][axis]:
            return False
        return True

    def theta(self, step):
        """Returns theta position of the sweep at the specified timestep index."""
        return self.tp[step][0]

    def phi(self, step):
        """Returns phi position of the sweep at the specified timestep index."""
        return self.tp[step][1]
    
    def check_continuity(self, stepsize, posmodel):
        """Checks that no abs delta theta or abs delta phi is greater than
        stepsize. Returns True if continous by this definition, False if not.
        Note that this check only makes sense to do after the sweep has been
        quantized. Returns True if only 0 or 1 elements exist in the sweep.
        
        Note that a jump e.g. from 180 deg to -180 deg (as defined in posintTP
        coordinates) would intentionally fail this test, because mechanically
        these do in fact represent an enormous discontinuity (due to theta
        hardstop).
        """
        L = len(self.tp)
        if L < 2:
            return True
        posintTP = [posmodel.trans.poslocTP_to_posintTP(self.tp[i]) for i in range(0,L)]
        abs_delta_t = [abs(posintTP[i][0] - posintTP[i-1][0]) for i in range(1,L)]
        abs_delta_p = [abs(posintTP[i][1] - posintTP[i-1][1]) for i in range(1,L)]
        if max(abs_delta_t) > stepsize or max(abs_delta_p) > stepsize:
            return False
        return True

from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport sin as c_sin
from libc.math cimport cos as c_cos
from libc.math cimport pi as c_pi
from libc.math cimport fmax as c_fmax
from libc.math cimport fmin as c_fmin
from libc.math cimport atan2 as c_atan2
cdef double rad_per_deg = c_pi / 180.0
cdef double deg_per_rad = 180.0 / c_pi

cdef class PosPoly:
    """Represents a collidable polygonal envelope definition for a mechanical component
    of the fiber positioner.

        points ... normal usage: 2 x N list of x and y coordinates of polygon vertices
               ... expert usage: ignored

        close_polygon ... normal usage: whether to make an identical last point matching the first
                      ... expert usage: ignored

        allocation ... normal usage: leave False, and it is ignored
                   ... expert usage: positive integer, stating just a size of arrays to allocate.
                                     (overrides both points and close_polygon)
    """
    cdef double* x
    cdef double* y
    cdef unsigned int n_pts

    def __cinit__(self, points, close_polygon=True, allocation=False):
        cdef size_t xlen
        cdef size_t ylen
        cdef size_t extra = 0
        cdef size_t xlen_plus
        cdef unsigned int fill_zeros = False
        cdef unsigned int i
        cdef unsigned int allocation_qty = allocation

        if allocation_qty > 0:
                xlen = allocation_qty
                fill_zeros = True
        else:
            xlen = len(points[0])
            ylen = len(points[1])
            if close_polygon:
                extra = 1
            if xlen != ylen:
                print('Error: different lengths x (' + str(xlen) + ') and y (' + str(ylen) + ')')
                return
        xlen_plus = xlen + extra
        self.x = <double*> PyMem_Malloc((xlen_plus) * sizeof(double))
        self.y = <double*> PyMem_Malloc((xlen_plus) * sizeof(double))
        if fill_zeros:
            for i in range(xlen_plus):
                self.x[i] = 0.0
                self.y[i] = 0.0
        else:
            for i in range(xlen):
                self.x[i] = points[0][i]
                self.y[i] = points[1][i]
            if extra:
                self.x[xlen] = self.x[0]
                self.y[xlen] = self.y[0]
        self.n_pts = xlen_plus
        if not self.x or not self.y:
            raise MemoryError()

    def __dealloc__(self):
        PyMem_Free(self.x)
        PyMem_Free(self.y)

    def __repr__(self):
        points = self.points
        output = 'PosPoly ['
        for i in range(len(points)):
            output += '['
            output += ', '.join([format(x,'.3f') for x in points[i]])
            output += '],'
        output = output[:-1]
        output += ']'
        return output

    @property
    def points(self):
        xy = [[],[]]
        for i in range(self.n_pts):
            xy[0].append(self.x[i])
            xy[1].append(self.y[i])
        return xy

    cpdef PosPoly copy(self):
        cdef PosPoly new = PosPoly(points=0, close_polygon=False, allocation=self.n_pts)
        cdef unsigned int i
        for i in range(new.n_pts):
            new.x[i] = self.x[i]
            new.y[i] = self.y[i]
        return new

    cpdef PosPoly rotated(self, angle):
        """Returns a copy of the polygon object, with points rotated by angle (unit degrees)."""
        cdef PosPoly new = self.copy()
        cdef double a = angle
        a *= rad_per_deg
        cdef double c = c_cos(a)
        cdef double s = c_sin(a)
        cdef double this_x
        cdef unsigned int i
        for i in range(new.n_pts):
            this_x = new.x[i]
            new.x[i] = c*this_x + -s*new.y[i]
            new.y[i] = s*this_x +  c*new.y[i]
        return new

    cpdef PosPoly translated(self, dx, dy):
        """Returns a copy of the polygon object, with points translated by distance (dx,dy)."""
        cdef PosPoly new = self.copy()
        cdef unsigned int i
        cdef double delta_x = dx
        cdef double delta_y = dy
        for i in range(new.n_pts):
            new.x[i] += delta_x
            new.y[i] += delta_y
        return new

    cpdef PosPoly expanded_radially(self, dR):
        """Returns a copy of the polygon object, with points expanded radially by distance dR.
        The linear expansion is made along lines from the polygon's [0,0] center.
        A value dR < 0 is allowed, causing contraction of the polygon.
        """
        cdef PosPoly new = self.copy()
        cdef unsigned int i
        cdef double delta_r = dR
        cdef double angle
        for i in range(new.n_pts):
            angle = c_atan2(new.y[i], new.x[i])
            new.x[i] += delta_r * c_cos(angle)
            new.y[i] += delta_r * c_sin(angle)
        return new

    cpdef PosPoly expanded_x(self, left_shift, right_shift):
        """Returns a copy of the polygon object, with points expanded along the x direction only.
        Points leftward of the line x=0 are shifted further left by an amount left_shift > 0.
        Points rightward of the line x=0 are shifted further right by an amount right_shift > 0.
        Points upon the line x=0 are not shifted.
        Negative values for left_dx and right_dx are allowed. They contract the polygon rightward and leftward, respectively.
        """
        cdef PosPoly new = self.copy()
        cdef unsigned int i
        cdef double L = left_shift
        cdef double R = right_shift
        for i in range(new.n_pts):
            if new.x[i] > 0:
                new.x[i] += R
            elif new.x[i] < 0:
                new.x[i] -= L
        return new

    cpdef PosPoly expanded_angularly(self, dA):
        """Returns a copy of the polygon object, with points expanded rotationally by angle dA (in degrees).
        The rotational expansion is made about the polygon's [0,0] center.
        Expansion is made in both the clockwise and counter-clockwise directions from the line y=0.
        A value dA < 0 is allowed, causing contraction of the polygon.
        """
        cdef PosPoly new = self.copy()
        cdef unsigned int i
        cdef double delta_angle = dA
        cdef double this_angle
        cdef double radius
        delta_angle *= rad_per_deg
        for i in range(new.n_pts):
            this_angle = c_atan2(new.y[i], new.x[i])
            if this_angle > 0:
                this_angle += delta_angle
            elif this_angle < 0:
                this_angle -= delta_angle
            radius = (new.x[i]**2 + new.y[i]**2)**0.5
            new.x[i] = radius * c_cos(this_angle)
            new.y[i] = radius * c_sin(this_angle)
        return new

    cpdef PosPoly place_as_phi_arm(self, theta, phi, x0, y0, r1):
        """Treating polygon as a phi arm, rotates and translates it to position defined
        by angles theta and phi (deg) and calibration values x0, y0, r1.
        """
        cdef PosPoly new = self.rotated(theta + phi) # units deg
        cdef double X0 = x0
        cdef double Y0 = y0
        cdef double R1 = r1
        cdef double T = theta
        cdef unsigned int i
        T *= rad_per_deg
        cdef double delta_x = X0 + R1 * c_cos(T) # now radians
        cdef double delta_y = Y0 + R1 * c_sin(T)
        for i in range(new.n_pts):
            new.x[i] += delta_x
            new.y[i] += delta_y
        return new

    cpdef PosPoly place_as_central_body(self, theta, x0, y0):
        """Treating polygon as a central body, rotates and translates it to postion
        defined by angle theta (deg) and calibration values x0, y0.
        """
        cdef PosPoly new = self.rotated(theta)
        cdef double delta_x = x0
        cdef double delta_y = y0
        cdef unsigned int i
        for i in range(new.n_pts):
            new.x[i] += delta_x
            new.y[i] += delta_y
        return new

    cpdef unsigned int collides_with(self, PosPoly other):
        """Searches for collisions in space between this polygon and
        another PosPoly object. Returns a bool, where true indicates a
        collision.
        """
        if _bounding_boxes_collide(self.x, self.y, self.n_pts, other.x, other.y, other.n_pts):
            return _polygons_collide(self.x, self.y, self.n_pts, other.x, other.y, other.n_pts)
        else:
            return False
        
    cpdef unsigned int collides_with_circle(self, x, y, radius):
        """Searches for collisions in space between this polygon and
        a circle, defined by center (x,y) and radius. Returns a bool,
        where true indicates a collision.
        """
        cdef double X = x
        cdef double Y = y
        cdef double R = radius
        cdef double distance
        cdef unsigned int i
        for i in range(self.n_pts):
            distance = ((self.x[i] - X)**2 + (self.y[i] - Y)**2)**0.5
            if distance < R:
                return True
        return False

cdef unsigned int _bounding_boxes_collide(double x1[], double y1[], unsigned int len1, double x2[], double y2[], unsigned int len2):
    """Check whether the rectangular bounding boxes of two polygons collide.

       x1,y1 ... 1xN c-arrays of the 1st polygon's vertices
       len1  ... length of x1 and y1
       x2,y2,len ... similarly for polygon 2

    Returns True if the bounding boxes collide, False if they do not. Intended as
    a fast check, to enhance speed of other collision detection functions.
    """
    if   _c_max(x1,len1) < _c_min(x2,len2):
        return False
    elif _c_max(y1,len1) < _c_min(y2,len2):
        return False
    elif _c_max(x2,len2) < _c_min(x1,len1):
        return False
    elif _c_max(y2,len2) < _c_min(y1,len1):
        return False
    else:
        return True

cdef unsigned int _polygons_collide(double x1[], double y1[], unsigned int len1, double x2[], double y2[], unsigned int len2):
    """Check whether two closed polygons collide.

       x1,y1 ... 1xN c-arrays of the 1st polygon's vertices
       len1  ... length of x1 and y1
       x2,y2,len ... similarly for polygon 2

    Returns True if the polygons intersect, False if they do not.

    The algorithm is by detecting intersection of line segments, therefore the case of
    a small polygon completely enclosed by a larger polygon will return False. Not checking
    for this condition admittedly breaks some conceptual logic, but this case is not
    anticipated to occur given the DESI petal geometry, and speed is at a premium.
    """
    cdef double[2] A1,A2,B1,B2
    cdef unsigned int i,j
    for i in range(len1 - 1):
        A1 = [x1[i],   y1[i]]
        A2 = [x1[i+1], y1[i+1]]
        for j in range(len2 - 1):
            B1 = [x2[j],   y2[j]]
            B2 = [x2[j+1], y2[j+1]]
            if _segments_intersect(A1,A2,B1,B2):
                return True
    return False

@cython.cdivision(True)
cdef unsigned int _segments_intersect(double A1[2], double A2[2], double B1[2], double B2[2]):
    """Checks whether two 2d line segments intersect. The endpoints for segments
    A and B are each a pair of (x,y) coordinates.
    """
    dx_A = A2[0] - A1[0]
    dy_A = A2[1] - A1[1]
    dx_B = B2[0] - B1[0]
    dy_B = B2[1] - B1[1]
    delta = dx_B * dy_A - dy_B * dx_A
    if delta == 0.0:
        return False  # parallel segments
    s = (dx_A * (B1[1] - A1[1]) + dy_A * (A1[0] - B1[0])) / delta
    t = (dx_B * (A1[1] - B1[1]) + dy_B * (B1[0] - A1[0])) / (-delta)
    return (0 <= s <= 1) and (0 <= t <= 1)

cdef double _c_max(double x[], unsigned int length):
    """Get the max of a c-array of doubles."""
    cdef unsigned int i
    cdef double max_val
    if length >= 1:
        max_val = x[0]
        for i in range(1,length):
            max_val = c_fmax(x[i],max_val)
        return max_val
    else:
        return 0.0

cdef double _c_min(double x[], unsigned int length):
    """Get the min of a c-array of doubles."""
    cdef unsigned int i
    cdef double min_val
    if length >= 1:
        min_val = x[0]
        for i in range(1,length):
            min_val = c_fmin(x[i],min_val)
        return min_val
    else:
        return 0.0

cpdef test():
    polys = []
    polys += [PosPoly([[0,1,1],[0,0,1]])]
    polys += [polys[0].translated(0.5,0)]
    polys += [polys[0].translated(10,0)]
    polys += [polys[1].rotated(30)]
    polys += [PosPoly([[0,1,2,3,4,5,6,7,8,9],[10,11,12,13,14,15,16,17,18,19]])]
    polys += [polys[4].rotated(45)]
    polys += [PosPoly([[20,10,43.2,-74],[1,2,3,4]])]
    bools = []
    bools += [polys[0].collides_with(polys[0])]
    bools += [polys[0].collides_with(polys[1])]
    bools += [polys[0].collides_with(polys[2])]
    bools += [polys[0].collides_with(polys[3])]
    bools += [polys[0].collides_with(polys[4])]
    bools += [polys[4].collides_with(polys[5])]
    bools += [polys[6].collides_with(polys[0])]
    t = 20
    p = -100
    x = 10
    y = -4
    r1 = 3
    a = polys[0].rotated(p)
    a = a.translated(r1,0)
    a = a.rotated(t)
    a = a.translated(x,y)
    b = polys[0].place_as_phi_arm(t,p,x,y,r1)
    return polys,bools, a, b

cpdef test2():
    t = {'nrows':5,
              'dT':[10,-20,0,0,0],
              'dP':[ 0,  0,-10,20,-10],
              'Tdot':[10,10,1,10,20],
              'Pdot':[5,5,5,5,5],
              'prepause':[0,1,0,0,0],
              'postpause':[0,0,0,0,1]}
    t['move_time'] = [max(abs(t['dT'][i]/t['Tdot'][i]), abs(t['dP'][i]/t['Pdot'][i])) for i in range(t['nrows'])]
    s = PosSweep('posid')
    s.fill_exact([100,-100], t, start_time=10)
    q = s.copy()
    q.quantize(0.1)
    return t, s, q

cpdef test3():
    keepout_phi = [[3.967, 3.918, 3.269, -1.172, -1.172,  3.269,  3.918],
                  [0.000, 1.014, 1.583,  1.037, -1.037, -1.583, -1.014]]
    p = PosPoly(keepout_phi)
    return p
