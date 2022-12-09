from amuse.lab import *
from amuse.units.units import *
from amuse.units import constants

import orbital_elements as orb_elem
from amuse.community.hermite_grx.interface import *

import numpy as np
import time
from math import *
import subprocess
import os

def evolve_sim(grav, out_particles, perturbation, t_end):
	start = time.time()
	# Try evolving
	grav.evolve_model(t_end)
	
	E_delta = 0 | kg * m**2 / s**2
	
	while grav.stopping_conditions.collision_detection.is_set():		
		# Collision between black hole and star
		p = grav.stopping_conditions.collision_detection.particles
		
		colliding_particles = Particles()
		if p(1).mass < p(0).mass:
			colliding_particles.add_particle(p(0))
			colliding_particles.add_particle(p(1))
		else:
			colliding_particles.add_particle(p(1))
			colliding_particles.add_particle(p(0))
		
		# Print orbital param
		try:
			print(orb_elem.orbital_elements_from_binary(colliding_particles, constants.G))
		except Exception:
			pass
		
		if perturbation != None:
			E_begin = grav.get_total_energy_with(perturbation)
			
			# Merge the star into the black hole, taking into account momentum and center of mass
			grav.large_particles[0].position = colliding_particles.center_of_mass()
			grav.large_particles[0].velocity = colliding_particles.center_of_mass_velocity()
			grav.large_particles[0].mass = colliding_particles.total_mass()
		
			# Remove the star from the simulation
			grav.small_particles.remove_particle(colliding_particles[1])
			
			E_end = grav.get_total_energy_with(perturbation)
		else:
			E_begin = grav.get_kinetic_energy() + grav.get_potential_energy()
			
			# Remove star particle and save mbh
			grav.particles.remove_particle(colliding_particles[1])
			mbh = grav.particles[np.nonzero(grav.particles.key == colliding_particles[0].key)[0]]
			
			# Add new particle which is the merged of the colliding particles
			mbh.position = colliding_particles.center_of_mass()
			mbh.velocity = colliding_particles.center_of_mass_velocity()
			mbh.mass = colliding_particles.total_mass()
			
			E_end = grav.get_kinetic_energy() + grav.get_potential_energy()
		
		# Record energy difference
		E_delta += E_begin - E_end
		
		# Setting its mass zero in out_particles
		index = np.nonzero(out_particles.key == colliding_particles[1].key)[0]
		out_particles[index].mass = 0 | MSun
		
		grav.particles.new_channel_to(out_particles).copy()
		
		print(out_particles.mass)
		
		print('RESOLVED A COLLISION at t =', grav.get_time().in_(yr))
		
		# Reduce timeout by number of seconds already evolved
		if grav.stopping_conditions.timeout_detection.is_enabled():
			end = time.time()
			timeevolved = (end - start) | s
			start = time.time()
			
			timeout = grav.parameters.stopping_conditions_timeout
			grav.parameters.stopping_conditions_timeout = timeout - timeevolved
		
		# Try evolving again
		grav.evolve_model(t_end)
	
	return E_delta

def setup_sim(initial, num_threads, dt_param, integrator, perturbation):
	mbh, bhs, converter = initial
	
	if perturbation != None:
		grav = HermiteGRX(converter)
		grav.parameters.light_speed = constants.c
		grav.parameters.perturbation = perturbation
		grav.parameters.integrator = integrator
		grav.parameters.num_threads = num_threads
		
		grav.large_particles.add_particle(mbh)
		grav.small_particles.add_particles(bhs)
	else:
		grav = Hermite(converter, number_of_workers = num_threads)
		grav.particles.add_particle(mbh)
		grav.particles.add_particles(bhs)
	
	grav.parameters.dt_param = dt_param
	grav.stopping_conditions.collision_detection.enable()
	
	out_particles = Particles()
	out_particles.add_particle(mbh)
	out_particles.add_particles(bhs)
	
	channel = grav.particles.new_channel_to(out_particles)
	channel.copy()
	
	return grav, out_particles, channel

# Load the latest snapshot of a saved run.
# Create a new saved run if it doesn't exist.
def load_snapshot(initial, runname):
	mbh, bhs, converter = initial
	
	snapshot_number = 0
	
	# Check if run already exists
	directory = '../../Data/%s/' % runname
	if os.path.isdir(directory):
		if len(os.listdir(directory)) != 0:
			# Read the last snapshot
			files = np.sort(os.listdir(directory))
			lastfile = files[-1]
			while not lastfile.endswith('.hdf'):
				files = files[:-1]
				lastfile = files[-1]
			
			particles = read_set_from_file(directory + lastfile, 'amuse', close_file = True)
			mbh = particles[0]
			bhs = particles[1:]
			
			# Remove particles with zero mass as they were captured by the mbh
			bhs.remove_particles(bhs[bhs.mass == 0 | MSun])
			
			snapshot_number = int(files[-1][:-4])
			
			print('Loaded snapshot #%d of %s.' % (snapshot_number, runname))
	else:
		os.mkdir(directory)
		
		print('Started run %s.' % runname)
	
	return snapshot_number, directory, (mbh, bhs, converter)

def save_snapshot(directory, out_particles, snapshot_number, E, E_delta, w):
	filename = '%05d.hdf' % snapshot_number
	write_set_to_file(out_particles, directory + filename, 'amuse')
	
	if snapshot_number != 0:
		dat = np.loadtxt(directory + 'info.txt')
		if dat.ndim == 1:
			# Only one entry
			dat = [dat]
		else:
			dat = list(dat)
		
		E_init = dat[0][1] | kg * m**2 / s**2
		E_delta += dat[-1][3] | kg * m**2 / s**2
		w += dat[-1][4]
	else:
		E_init = E
		dat = []
	
	dat.append(np.array([snapshot_number, E.value_in(kg * m**2 / s**2), (E - E_init + E_delta) / E_init, E_delta.value_in(kg * m**2 / s**2), w]))
	np.savetxt(directory + 'info.txt', dat)
	
	print('Saved snapshot #%d.' % snapshot_number)

def saved_run(initial, t_end, num_threads, dt_param, integrator, perturbation, runname, dt=100|yr):
	i, directory, initial = load_snapshot(initial, runname)
	grav, out_particles, channel = setup_sim(initial, num_threads, dt_param, integrator, perturbation)
	
	if i == 0:
		if perturbation != None:
			E = grav.get_total_energy_with(perturbation)[0]
		else:
			E = grav.get_kinetic_energy() + grav.get_potential_energy()
		save_snapshot(directory, out_particles, i, E, 0 | kg * m**2 / s**2, 0)
	
	t = i * dt
	t_start = t
	
	while t < t_end:
		t += dt
		i += 1
		
		start = time.time()
		E_delta = evolve_sim(grav, out_particles, perturbation, t - t_start)
		end = time.time()
		w = end - start
		
		if perturbation != None:
			E = grav.get_total_energy_with(perturbation)[0]
		else:
			E = grav.get_kinetic_energy() + grav.get_potential_energy()
		
		print('t        =', (grav.get_time() + t_start).in_(yr))
		print('remaining=', w * (t_end - grav.get_time()) / dt / 60, 'min')
		print('E_tot    =', E)
		
		channel.copy()
		save_snapshot(directory, out_particles, i, E, E_delta, w)
	
	grav.stop()

def timed_run(initial, num_threads, dt_param, integrator, perturbation, t_end = 1e6 | yr, timeout= -1 | s):
	grav, out_particles, channel = setup_sim(initial, num_threads, dt_param, integrator, perturbation)
	
	evolve_sim(grav, out_particles, perturbation, 0.0001 | yr)
	
	if perturbation == None:
		E_init = grav.get_kinetic_energy() + grav.get_potential_energy()
	else:
		E_init = grav.get_total_energy_with(perturbation)
	
	if timeout > 0 | s:
		grav.stopping_conditions.timeout_detection.enable()
		grav.parameters.stopping_conditions_timeout = timeout
	
	start = time.time()
	E_delta = evolve_sim(grav, out_particles, perturbation, t_end)
	end = time.time()
	
	if perturbation == None:
		E_post = grav.get_kinetic_energy() + grav.get_potential_energy()
	else:
		E_post = grav.get_total_energy_with(perturbation)
	
	gtime = grav.get_time()
	
	grav.stop()
	return end - start, E_init, E_post+E_delta, gtime
