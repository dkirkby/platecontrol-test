'''
fvc_measure.py

Loops over fvc_measure from PECS
'''

import time
import argparse
import traceback
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-r', '--match_radius', type=int, default=None, help='int, specify a particular match radius, other than default')
parser.add_argument('-t', '--exposure_time', type=float, default=2.0, help='float, specify exposure time for the fvc.')
parser.add_argument('-l', '--loops', type=int, default=1, help='int, number of times to loop over fvc_measure.')
parser.add_argument('--cadence', type=float, default=0.0, help='float, loop at this cadence in seconds or with no delay if 0 (default)')
max_fvc_iter = 100
parser.add_argument('-nm', '--num_meas', type=int, default=1, help=f'int, number of measurements by the FVC per move (default is 1, max is {max_fvc_iter})')
uargs = parser.parse_args()
assert 1 <= uargs.num_meas <= max_fvc_iter, f'out of range argument {uargs.num_meas} for num_meas parameter'
assert uargs.loops > 0, f'out of range argument, {uargs.loops} cannot be negative!'

from pecs import PECS
cs = PECS(interactive=True, test_name='pecs_fvc_measure')
cs.exptime = uargs.exposure_time
try:
	last = 0
	for i in range(uargs.loops):
		while True:
			now = time.time()
			remaining = last + uargs.cadence - now
			if remaining <= 0:
				break
			print(f'{remaining:.0f}s until next image...')
			time.sleep(1)
		last = now
		print(f'Taking image {i} of {uargs.loops}....')
		cs.fvc_measure(match_radius=uargs.match_radius, num_meas=uargs.num_meas)
except:
	traceback.print_exc()
cs.fvc_collect()
