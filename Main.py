from amuse.lab import *
from amuse.units.units import *
from amuse.units.optparse import OptionParser

from InitialConditions import *
from GalacticCenter import *

def main(Nstars=50, Ntest=5, a=2e-3 | parsec, ecc=1-1e-2, pert='None', runname=None, integ='Hermite', dtparam=0.01, tend=1 | Myr, numthreads=1, dt=100 | yr):
	print(Nstars, Ntest, a, ecc, pert, runname, integ, dtparam, tend, numthreads, dt)
	initial = get_initial_conditions(Nstars-Ntest, Ntest, a, ecc)
	
	if runname == None:
		print('Doing a timed run.')
		time, Einit, Epost, dummy = timed_run(initial, numthreads, dtparam, integ, pert, tend)
		print('Simulation took', time, 'seconds.')
		print('Energy error was', (Epost - Einit) / Einit)
	else:
		saved_run(initial, tend, numthreads, dtparam, integ, pert, runname, dt=dt)

def new_option_parser():
	result = OptionParser()
	
	result.add_option('--Nstars', dest='Nstars', type='int', default=50, help='number of stars [%default]')
	result.add_option('--Ntest', dest='Ntest', type='int', default=5, help='number of test particles [%default]')
	result.add_option('--a', dest='a', type='float', unit=parsec, default=2e-3 | parsec, help='initial semi-major axis for the test particles in parsec [%default]')
	result.add_option('--ecc', dest='ecc', type='float', default = 1-1e-2, help='initial eccentricity for the test particles [%default]')
	result.add_option('--pert', dest='pert', default='None', help='the perturbation to use [%default]')
	result.add_option('--runname', dest='runname', help='the runname to save the simulation to [not saved]')
	result.add_option('--integ', dest='integ', default='Hermite', help='the integrator to use [%default]')
	result.add_option('--dtparam', dest='dtparam', type='float', default=0.01, help='the time step parameter to use [%default]')
	result.add_option('--tend', dest='tend', type='float', unit=Myr, default=1 | Myr, help='the end time in Myr [%default]')
	result.add_option('--numthreads', dest='numthreads', type='int', default=1, help='the number of threads to run simultaniously [%default]')
	result.add_option('--dt', dest='dt', type='float', unit=yr, default=100 | yr, help='the time between snapshot saves in yr [%default]')
	
	return result

if __name__ == '__main__':
	o, arguments = new_option_parser().parse_args()
	main(**o.__dict__)
