from amuse.lab import *
from amuse.units.units import *
from matplotlib import pyplot
from amuse.units import constants
import numpy as np
from hermitepn.interface import *
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.orbital_elements import orbital_elements_from_binary
           
def HTpulsar(m1, m2, a, e):
    converter = nbody_system.nbody_to_si(m1+m2, a)
    HTbinary = new_binary_from_orbital_elements(m1, m2, a, e,
                                                G=constants.G)
    HTbinary.radius = 10 | units.km
    return HTbinary, converter

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
    
def get_trajectories(initial, grav, t_end, dt, pert=None):
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

    sma = [] | units.RSun
    ecc = []
	
    t = 0 * t_end
	
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
        grav.evolve_model(t)

        M, m, a, e, ta_out, inc, lan_out, aop_out = orbital_elements_from_binary(grav.particles, G=constants.G)
        sma.append(a)
        ecc.append(e)
	
    return x1, y1, z1, x2, y2, z2, time, E, sma, ecc

def run_nbody_code(Nbody_code, label, dt, tend):
    grav = Nbody_code(converter)
    if Nbody_code==HermitePN:
        pert = '1PN_Pairwise'
        grav.parameters.integrator = 'RegularizedHermite'
        if "EIH" in label:
            pert = '1PN_EIH'
            grav.parameters.integrator = 'SymmetrizedHermite'
        print("v=", bodies[0].velocity.length()/grav.parameters.light_speed)
        grav.parameters.perturbation = pert
        grav.parameters.dt_param = 0.02
        grav.parameters.light_speed = 0.05*constants.c
        print("v=", bodies[0].velocity.length()/grav.parameters.light_speed)
    else:
        pert = 'None'
        grav.parameters.dt_param = 0.02

    x1, y1, z1, x2, y2, z2, time, E, sma, ecc = get_trajectories(initial,
                                                                 grav,
                                                                 t_end,
                                                                 dt,
                                                                 pert)
    return x1, y1, z1, x2, y2, z2, time, E, sma, ecc
m1 = 1.441 | units.MSun
m2 = 1.397 | units.MSun
a = 1950100 | units.km
e = 0.6171334
initial = HTpulsar(m1, m2, a, e)
bodies = initial[0]
Porb = 7.751938773864 | units.hour
converter = initial[1]
Nbody_codes = [Hermite, HermitePN, HermitePN]
colors = ['r', 'b', 'g']
labels = ["Hermite", "hermitePN", "HermitePN\_EIH"]
dt = 0.1*Porb
t_end = 10*Porb

figure = pyplot.figure(figsize = (10, 10))

for Nbody_code, ci, li in zip(Nbody_codes, colors, labels):
    x1, y1, z1, x2, y2, z2, time, E, sma, ecc = run_nbody_code(Nbody_code, li, dt, t_end)
    print("n=", len(x1))
    print("t=", time/Porb)

    subplot = figure.add_subplot(2, 2, 1)
    pyplot.plot(x1.value_in(units.au), y1.value_in(units.au), c=ci, lw=1)
    pyplot.plot(x2.value_in(units.au), y2.value_in(units.au), c=ci, lw=1)
    pyplot.xlabel('x')
    pyplot.ylabel('y')

    subplot = figure.add_subplot(2, 2, 2)
    pyplot.plot(sma.value_in(units.RSun), ecc, c=ci, label=li)
    pyplot.xlabel('a')
    pyplot.ylabel('e')

    subplot = figure.add_subplot(2, 2, 3)
    pyplot.plot(time.value_in(units.yr), sma.value_in(units.RSun), c=ci)
    pyplot.xlabel('t')
    pyplot.ylabel('a')

    subplot = figure.add_subplot(2, 2, 4)
    pyplot.plot(time.value_in(units.yr), ecc, c=ci)
    pyplot.xlabel('t')
    pyplot.ylabel('e')
    
subplot = figure.add_subplot(2, 2, 2)
pyplot.legend()

#pyplot.show()
pyplot.savefig("HTpulsar.pdf", fontsize=10)
