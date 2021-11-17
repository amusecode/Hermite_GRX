from amuse.lab import *
from amuse.units.units import *
from amuse.units import constants
import numpy as np
from hermitepn.interface import *

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
    if pert == 'None':
        return grav.particles.kinetic_energy() + \
            grav.particles.potential_energy()
    else:
        return grav.get_total_energy_with(pert)[0]
    
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
    time = [] |units.yr
    E = [] | units.J 
	
    t = 0 * t_end
    dt = t_end / 10000.0
    i = 0.
	
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
	
        time.append(t)
        E.append(get_energy(grav, pert))
		
        t += dt
        i += 1
        grav.evolve_model(t)
		
    return x1, y1, z1, x2, y2, z2, x3, y3, z3, time, E
    
initial = get_initial_conditions_figure_eight(1 | units.MSun,
                                              1 | units.RSun)
#initial[0].velocity += 10 |units.kms
bodies = initial[0]
converter = initial[1]
Nbody_code = HermitePN
#Nbody_code = Hermite
#grav = HermitePN(converter)
grav = Nbody_code(converter)
pert = '1PN_Pairwise'
if Nbody_code=="HermitePN":
    print("perturbation=", grav.parameters.perturbation)
    grav.parameters.perturbation = pert
    grav.parameters.integrator = 'RegularizedHermite'
    grav.parameters.dt_param = 0.1
    grav.parameters.light_speed = 0.001*constants.c
else:
    grav.parameters.dt_param = 0.1
print("v=", bodies[0].velocity.length()/grav.parameters.light_speed)
    
x1, y1, z1, x2, y2, z2, x3, y3, z3, time, E = get_trajectories(initial,
                                                      grav,
                                                      0.001 | units.yr,
                                                      pert)

print("n=", len(x1))
from matplotlib import pyplot
pyplot.plot(x1.value_in(units.au), y1.value_in(units.au), c='r', lw=1)
pyplot.plot(x2.value_in(units.au), y2.value_in(units.au), c='g', lw=1)
pyplot.plot(x3.value_in(units.au), y3.value_in(units.au), c='b', lw=1)
pyplot.show()

pyplot.scatter(time.value_in(units.yr), E/E[0], c='g')
aylo
