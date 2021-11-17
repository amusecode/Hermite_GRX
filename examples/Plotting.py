from amuse.lab import *
from amuse.units.units import *
from amuse.units import constants

import matplotlib.pyplot as plt
import orbital_elements as orb_elem

import numpy as np
from numpy.linalg import inv
from math import *
import os
import cairo
import scipy.misc

import GalacticCenter as GC

def get_orbit(binary, lengthunit, G=nbody_system.G):
	try:
		m1, m2, a, ecc, true_anom, incl, long_asc, arg_per = orb_elem.orbital_elements_from_binary(binary, G)
	except Exception:
		return None
	
	arg_per = radians(arg_per)
	long_asc = radians(long_asc)
	incl = radians(incl)
	x1 = cos(long_asc) * cos(arg_per) - sin(long_asc) * cos(incl) * sin(arg_per)
	x2 = sin(long_asc) * cos(arg_per) + cos(long_asc) * cos(incl) * sin(arg_per)
	x3 = sin(incl) * sin(arg_per)
	y1 = -cos(long_asc) * sin(arg_per) - sin(long_asc) * cos(incl) * cos(arg_per)
	y2 = -sin(long_asc) * sin(arg_per) + cos(long_asc) * cos(incl) * cos(arg_per)
	y3 = sin(incl) * cos(arg_per)
	z1 = sin(incl) * sin(long_asc)
	z2 = -sin(incl) * cos(long_asc)
	z3 = cos(incl)
	
	M = inv(np.array([[x1, x2, x3], [y1, y2, y3], [z1, z2, z3]]))
	
	N = 500
	f = np.linspace(0, 2*np.pi, N) + true_anom
	r = a * (1-ecc**2) / (1 + ecc * np.cos(f))
	x = np.cos(f) * r
	y = np.sin(f) * r
	z = np.zeros_like(f) | lengthunit
	
	pos = np.array([x.value_in(lengthunit),y.value_in(lengthunit),z.value_in(lengthunit)])
	x, y, z = np.dot(M, pos) | lengthunit
	
	cm = binary.center_of_mass()
	x += cm.x
	y += cm.y
	z += cm.z
	
	P = 2*np.pi*(a**3/(G*(m1+m2))).sqrt().value_in(yr)
	alpha = 0.2 + 0.8 * np.exp(-np.linspace(0,1,N)*P/10)
	
	return x, y, z, alpha

def plot_bh_stars(grav):
	plt.clf()
	plt.plot(grav.large_particles.x.value_in(parsec), grav.large_particles.y.value_in(parsec), 'r+')
	plt.plot(grav.small_particles.x.value_in(parsec), grav.small_particles.y.value_in(parsec), 'b+')
	plt.axis('equal')
	plt.xlim([-0.005, 0.005])
	plt.ylim([-0.005, 0.005])
	bh = Particle()
	bh.position = grav.large_particles[0].position
	bh.velocity = grav.large_particles[0].velocity
	bh.mass = grav.large_particles[0].mass
	
	for p in grav.small_particles:
		p1 = Particle()
		p1.position = p.position
		p1.velocity = p.velocity
		p1.mass = p.mass
		
		binary = Particles()
		binary.add_particle(bh)
		binary.add_particle(p1)
		
		res = get_orbit(binary, parsec, constants.G)
		if res != None:
			x, y, z, a = res
			plt.plot(x.value_in(parsec), y.value_in(parsec), 'b')
		
	plt.title('t=%2f yr' % grav.get_time().value_in(yr))
	
def cairo_bh_stars(runname, snapshot_number):
	filename = '%05d.hdf' % snapshot_number
	directory = '../../Data/%s/' % runname
	
	particles = read_set_from_file(directory + filename, 'amuse', close_file = True)
	
	mbh = Particle()
	mbh.position = particles[0].position
	mbh.velocity = particles[0].velocity
	mbh.mass = particles[0].mass
	
	bhs = particles[1:]
	bhs.remove_particles(bhs[bhs.mass == 0 | MSun])
	
	bhs.position -= mbh.position
	mbh.position -= mbh.position
	
	image = np.zeros((800,800,4), dtype=np.uint8)
	surface = cairo.ImageSurface.create_for_data(image, cairo.FORMAT_ARGB32, 800, 800)
	cr = cairo.Context(surface)
	
	cr.set_source_rgb(0,0,0)
	cr.paint()
	
	for p in bhs:
		p1 = Particle()
		p1.position = p.position
		p1.velocity = p.velocity
		p1.mass = p.mass
		
		binary = Particles()
		binary.add_particle(mbh)
		binary.add_particle(p1)
		
		res = get_orbit(binary, parsec, constants.G)
		if res == None:
			continue
		
		x, y, z, a = res
		x = x.value_in(parsec)*25000
		y = y.value_in(parsec)*25000
		z = z.value_in(parsec)*25000
		
		last = np.array([x[-1], y[-1]])
		
		for i,j,k,l in np.array([x,y,z,a]).T:
			cr.move_to(last[0]+400, last[1]+400)
			cr.line_to(i+400, j+400)
			cr.set_source_rgb(l,l,l)
			cr.set_line_width(1)
			cr.stroke()
			last = np.array([i, j])
	
	directory = '../../Movies/%s/' % runname
	filename = '%d.png' % snapshot_number
	
	if not os.path.isdir(directory):
		os.mkdir(directory)
	
	scipy.misc.imsave(directory + filename, image)

def cairo_whole_run(runname):
	for i in range(1000):
		cairo_bh_stars(runname, i+9000)
	
	direct = 'cd ../../Movies/%s/; ' % runname
	os.system(direct + 'ffmpeg -y -f image2 -i %d.png -vcodec mpeg4 -framerate 30 -q:v 5 movie.avi')
	os.system(direct + 'rm *.png')

def gnuplot_bh_stars(runname, snapshot_number):
	filename = '%05d.hdf' % snapshot_number
	directory = '../../Data/%s/' % runname
	
	particles = read_set_from_file(directory + filename, 'amuse', close_file = True)
	
	mbh = Particle()
	mbh.position = particles[0].position
	mbh.velocity = particles[0].velocity
	mbh.mass = particles[0].mass
	
	bhs = particles[1:]
	bhs.remove_particles(bhs[bhs.mass == 0 | MSun])
	
	file_output = ''
	
	for p in bhs:
		p1 = Particle()
		p1.position = p.position
		p1.velocity = p.velocity
		p1.mass = p.mass
		
		binary = Particles()
		binary.add_particle(mbh)
		binary.add_particle(p1)
		
		res = get_orbit(binary, parsec, constants.G)
		if res == None:
			continue
		
		x, y, z, a = res
		x = x.value_in(parsec)
		y = y.value_in(parsec)
		z = z.value_in(parsec)
		
		for i,j,k,l in np.array([x,y,z,a]).T:
			file_output += '%f\t%f\t%f\t%f\n' % (i,j,k,l)
		file_output += '\n\n'
	
	filename = '%d.dat' % snapshot_number
	directory = '../../Movies/%s/' % runname
	
	if not os.path.isdir(directory):
		os.mkdir(directory)
	
	f = open(directory + filename, 'wb')
	f.write(file_output)
	f.close()

def gnuplot_whole_run(runname):
	for i in range(10000):
		gnuplot_bh_stars(runname, i)
	
	os.system('cd ../../Movies/%s/; gnuplot plot.gp')

if __name__ == '__main__':
	cairo_whole_run('1PN_Pairwise02')
