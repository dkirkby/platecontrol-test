# -*- coding: utf-8 -*-
"""
Generates xy target grids.

Created on Fri May  8 13:54:31 2020
@author: joe
"""

import numpy as np

def filled_circle(n_points=24, radius=3.3, random=False, verbose=False):
    '''Make a cartesian grid of n_points. Wrapper for filled_annulus.
    '''
    return filled_annulus(n_points, 0, radius, random, verbose)

def filled_annulus(n_points=24, r_min=0.0, r_max=3.3, random=False, verbose=False):
    '''Make a cartesian grid of n_points, filling an annular region.
    
    mirror symmetric about both axes,
    and filling a circle of size radius. Returns a list of size Nx2, where
    N is <= n_points (and as close to n_points as possible).
    '''
    if random:
        return _random_filled_annulus(n_points, r_min, r_max, verbose)
    else:
        return _regular_filled_annulus(n_points, r_min, r_max, verbose)

def _regular_filled_annulus(n_points, r_min, r_max, verbose, _is_first_pass=True):
    '''Unimplemented, to-do. As of 2020-05-14.'''
    n_init = int(n_points**0.5) # initial grid radius in integer num points
    odd = n_points % 2
    line = np.linspace(-1, 1, 2*n_init + odd) / np.sqrt(2)
    grid = [[x,y] for x in line for y in line]
    grid_r = [(xy[0]**2 + xy[1]**2)**0.5 for xy in grid]
    def n_points_within(inner_radius, outer_radius):
        within = [r for r in grid_r if inner_radius <= r <= outer_radius]
        return len(within)
    rate = 1 / n_points / 10 # a bit arbitrary, but seems to work fine
    cutoff = 0.0
    while cutoff <= 1.0:
        cutoff += rate
        inner_cutoff = r_min/r_max * cutoff
        n = n_points_within(inner_cutoff, cutoff)
        if verbose:
            print(f'count={n:5d},  normalized radii: inner={inner_cutoff:8.6f}, outer={cutoff:8.6f}')
        if n > n_points:
            cutoff -= rate
            break
        if n == n_points:
            break
    selected_r = [r for r in grid_r if r <= cutoff]
    selected_xy = [grid[i] for i in range(len(grid)) if grid_r[i] <= cutoff]
    scale = r_max / max(selected_r)
    scaled_xy = [[xy[0]*scale, xy[1]*scale] for xy in selected_xy]
    if _is_first_pass and n_points > 4:
        retry = _regular_filled_annulus(n_points=n_points-1, r_min=r_min, r_max=r_max,
                              verbose=verbose, _is_first_pass=False)
        current_err = len(scaled_xy) - n_points
        retry_err = len(retry) - n_points
        if abs(retry_err) < abs(current_err):
            scaled_xy = retry
        if verbose:
            scaled_r = [(xy[0]**2 + xy[1]**2)**0.5 for xy in scaled_xy]
            print(f'min radius = {min(scaled_r)}')
            print(f'max radius = {max(scaled_r)}')
            print(f'num points = {len(scaled_r)}')
            print(f'(x,y) locations: {scaled_xy}')
    return scaled_xy

def _random_filled_annulus(n_points, r_min, r_max, verbose):
    '''to-do'''
    return

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    reg_circle = filled_circle(n_points=24, radius=3.0, random=False, verbose=True)
    reg_annulus = filled_annulus(n_points=48, r_min=2.0, r_max=3.3, random=False, verbose=True)
    plt.clf()
    results = [reg_circle, reg_annulus]
    for result, plotnum in zip(results, range(1,len(results)+1)):
        xy = np.transpose(result)
        plt.subplot(1, len(results), plotnum)
        plt.plot(xy[0], xy[1], '-o')
        plt.axis('equal')
        plt.title(f'n_points = {len(result)}')