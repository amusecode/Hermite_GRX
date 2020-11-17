from amuse.lab import *
from amuse.units.units import *
from amuse.units import constants
import numpy as np
import matplotlib.pyplot as plt
from hermitepn.interface import *

import orbital_elements as orb_elem

def get_initial_conditions_pythagorean(unit_mass, unit_length):
	converter = nbody_system.nbody_to_si(12*unit_mass, unit_length)
	
	particles = Particles(3)
	particles.position = np.array([[1,3,0], [-2,-1,0], [1,-1,0]]) * unit_length
	particles.velocity = [[0,0,0],[0,0,0],[0,0,0]] | units.m / units.s
	particles.mass = np.array([3,4,5]) * unit_mass
	particles.radius = 0 * unit_length
	
	particles.move_to_center()
	
	return particles, converter

def get_initial_conditions_figure_eight(unit_mass, unit_length):
	converter = nbody_system.nbody_to_si(unit_mass, unit_length)
	
	particles = Particles(3)
	particles.position = np.array([[0.9700436, -0.24308753, 0], [-0.9700436, 0.24308753, 0], [0,0,0]]) * unit_length
	particles.velocity = np.array([[0.466203685, 0.43236573, 0], [0.466203685, 0.43236573, 0], [-2*0.466203685, -2*0.43236573, 0]]) * unit_length / converter.to_si(nbody_system.time)
	particles.mass = np.array([1,1,1]) * unit_mass
	particles.radius = 0 * unit_length
	
	particles.move_to_center()
	
	return particles, converter
	
def get_energy(grav, pert=None):
	if pert != None:
		return grav.get_total_energy_with(pert)[0]
	else:
		return grav.particles.kinetic_energy() + grav.particles.potential_energy()

def get_trajectories(initial, grav, t_end, pert=None):
	grav.particles.add_particles(initial[0])
	
	pre_energy = get_energy(grav, pert)
	x1 = [] | units.AU
	x2 = [] | units.AU
	x3 = [] | units.AU
	y1 = [] | units.AU
	y2 = [] | units.AU
	y3 = [] | units.AU
	z1 = [] | units.AU
	z2 = [] | units.AU
	z3 = [] | units.AU
	E = [] | units.J
	
	t = 0 * t_end
	dt = t_end / 10000.0
	i = 0
	
	plt.ion()
	while t < t_end:
		x1.append(grav.particles[0].position.x)
		y1.append(grav.particles[0].position.y)
		z1.append(grav.particles[0].position.z)
		x2.append(grav.particles[1].position.x)
		y2.append(grav.particles[1].position.y)
		z2.append(grav.particles[1].position.z)
		x3.append(grav.particles[2].position.x)
		y3.append(grav.particles[2].position.y)
		z3.append(grav.particles[2].position.z)
		
		E.append(get_energy(grav, pert))
		
		if i == 50:
			plt.clf()
			plt.plot(x1[-1].value_in(AU), y1[-1].value_in(AU), 'ro')
			plt.plot(x1.value_in(AU), y1.value_in(AU), 'r-')
			plt.plot(x2[-1].value_in(AU), y2[-1].value_in(AU), 'go')
			plt.plot(x2.value_in(AU), y2.value_in(AU), 'g-')
			plt.plot(x3[-1].value_in(AU), y3[-1].value_in(AU), 'bo')
			plt.plot(x3.value_in(AU), y3.value_in(AU), 'b-')
			plt.title('t=%f yr, (E-E[0])/E[0]=%g' % (grav.get_time().value_in(units.yr), (E[-1]-E[0])/E[0]))
			plt.xlim([-5,5])
			plt.ylim([-5,5])
			plt.draw()
			i = 0
		
		t += dt
		i += 1
		grav.evolve_model(t)
		
	plt.ioff()
	
	return E

initial = get_initial_conditions_pythagorean(1 | units.MSun, 1 | AU)
converter = initial[1]
grav = HermitePN(converter)
grav.parameters.perturbation = 'None'
grav.parameters.integrator = 'RegularizedHermite'
grav.parameters.dt_param = 0.001
grav.parameters.light_speed = constants.c
get_trajectories(initial, grav, 11 | units.yr, 'None')
