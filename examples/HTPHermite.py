from amuse.lab import *
from amuse.units.units import *
from amuse.units import constants
from amuse.community.hermite_grx.interface import *

def HTpulsar(m1, m2, a, e):
    converter = nbody_system.nbody_to_si(m1+m2, a)
    from amuse.ext.orbital_elements import new_binary_from_orbital_elements
    HTbinary = new_binary_from_orbital_elements(m1, m2, a, e,
                                                G=constants.G)
    HTbinary.radius = 10 | units.km
    return HTbinary, converter

    
def get_trajectories(initial, grav, t_end, pert=None):
    grav.particles.add_particles(initial[0])
	
    pre_energy = get_energy(grav, pert)
    x1 = [] | units.AU
    x2 = [] | units.AU
    y1 = [] | units.AU
    y2 = [] | units.AU
    z1 = [] | units.AU
    z2 = [] | units.AU
    time = [] |units.yr
    E = [] | units.J 
	
    t = 0 * t_end
    dt = t_end / 100.0
    i = 0.
	
    while t < t_end:
        x1.append(grav.particles[0].position.x)
        y1.append(grav.particles[0].position.y)
        z1.append(grav.particles[0].position.z)
        x2.append(grav.particles[1].position.x)
        y2.append(grav.particles[1].position.y)
        z2.append(grav.particles[1].position.z)
	
        time.append(t)
        E.append(get_energy(grav, pert))
		
        t += dt
        i += 1
        grav.evolve_model(t)
		
    return x1, y1, z1, x2, y2, z2, time, E

m1 = 1.441 | units.MSun
m2 = 1.397 | units.MSun
a = 1950100 | units.km
e = 0.6171334
initial = HTpulsar(m1, m2, a, e)
bodies = initial[0]
Porb = 7.751938773864 | units.hour

converter = initial[1]
Nbody_code = Hermite
grav = Nbody_code(converter)
grav.parameters.dt_param = 0.1
    
x1, y1, z1, x2, y2, z2, time, E = get_trajectories(initial,
                                                   grav,
                                                   10*Porb)
                                                   
print("n=", len(x1))
from matplotlib import pyplot
pyplot.plot(x1.value_in(units.au), y1.value_in(units.au), c='r', lw=1)
pyplot.plot(x2.value_in(units.au), y2.value_in(units.au), c='g', lw=1)
pyplot.show()

pyplot.scatter(time.value_in(units.yr), E/E[0], c='g')
pyplot.sh
