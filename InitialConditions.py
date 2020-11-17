from amuse.lab import *
from amuse.units.units import *
from amuse.units import constants
import orbital_elements as orb_elem
import numpy as np
from math import *
import random

def kepler_eq(E, M, e):
	return E - e * np.sin(E) - M

def kepler_eq_prime(E, M, e):
	return 1 - e * np.cos(E)

def true_anomaly_from_mean_anomaly(M, ecc):
	E = orb_elem.newton(kepler_eq, M, kepler_eq_prime, (M, ecc), maxiter=100)[0]
	true_anom = orb_elem.true_anomaly_from_eccentric_anomaly(E, ecc)
	return degrees(true_anom)

def get_initial_conditions(Nfield, Ntest, a_fixed, ecc_fixed):
	N = Nfield + Ntest
	
	M_mbh = 1e6 | MSun
	M_bh = (2500.0 / N) | MSun
	
	a_min = 0.1e-3 | parsec
	a_max = 10e-3 | parsec
	
	converter = nbody_system.nbody_to_si(N*M_bh+M_mbh, parsec)

	mbh = Particle()
	mbh.position = [0,0,0] | parsec
	mbh.velocity = [0,0,0] | parsec / yr
	mbh.mass = M_mbh
	mbh.radius = 8 * constants.G * M_mbh / constants.c**2

	bodies = Particles()
	bodies.add_particle(mbh)
	
	a_arr = np.random.uniform(a_min.value_in(parsec), a_max.value_in(parsec), N) | parsec
	ecc_arr = np.sqrt(np.random.uniform(0, 1, N))
	incl_arr = np.array([-1,1])[np.random.randint(0, 1, N)] * np.degrees(np.arccos(np.random.uniform(-1,1, N)))
	long_asc_arr = np.random.uniform(-180,180, N)
	arg_per_arr = np.random.uniform(-180,180, N)
	mean_anom_arr = np.random.uniform(0, 2*np.pi, N)
	
	if Ntest > 0:
		a_arr[0:Ntest] = a_fixed
		ecc_arr[0:Ntest] = ecc_fixed

	for a, ecc, incl, long_asc, arg_per, mean_anom in zip(a_arr, ecc_arr, incl_arr, long_asc_arr, arg_per_arr, mean_anom_arr):
		true_anom = true_anomaly_from_mean_anomaly(mean_anom, ecc)
		
		binary = orb_elem.new_binary_from_orbital_elements(M_mbh, M_bh, a, ecc, true_anom, incl, long_asc, arg_per, constants.G)
		
		bh = binary[1]
		bh.position -= binary[0].position
		bh.velocity -= binary[0].velocity
		bh.radius = 0 | parsec
		bodies.add_particle(bh)
	
	bodies.move_to_center()
	return bodies[0], bodies[1:], converter
