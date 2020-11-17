""" 

Routines to convert Newtonian Keplerian elements to Cartesian coordinates 
(new_binary_from_orbital_elements) and vice versa (orbital_elements_from_binary).
Note that the latter routine cannot give the true anomaly and argument of
periapse if the orbit is circular.

Last modification December 13 2013

"""

import numpy

from amuse.units import units,nbody_system,constants
from amuse.datamodel import Particles,rotation

def newton(f,x0,fprime=None,args=(),tol=1.48e-8,maxiter=50):
    if fprime is None:
        print("provide fprime")
        return x0
    i=0
    x=x0
    while (i<maxiter):
        fv=f(x,*args)
        dfv=fprime(x,*args)
        if(dfv==0):
            return x0,-2
        delta=-fv/dfv
        if(abs(delta)<tol):
            return x+delta,0
        x=x+delta
        i=i+1
    return x,-1    


def true_anomaly_from_eccentric_anomaly(E,e):
  return 2*numpy.arctan2((1+e)**0.5*numpy.sin(E/2),(1-e)**0.5*numpy.cos(E/2))

# E from M,e
# newton solver for M=E-e sin E

def new_binary_from_orbital_elements(
        mass1,
        mass2,
        semimajor_axis, 
        eccentricity = 0,
        true_anomaly = 0, 
        inclination = 0,
        longitude_of_the_ascending_nodes = 0,
        argument_of_periapsis = 0,
        G=nbody_system.G
  ):
  """ 

  Function that returns two-particle Particle set, with the second 
  particle position and velocities computed from the input orbital 
  elements. angles in degrees, inclination between 0 and 180

  """
    
  if semimajor_axis.value_in(units.m) <= 0.0:
    raise Exception("expects positive semimajor axes a > 0 (i.e. bound orbits)")    

  if eccentricity < 0.0 or eccentricity >= 1.0:
    raise Exception("expects eccentricities 0 <= e < 1 (i.e. bound orbits)")    
    
  inclination = numpy.radians(inclination)
  argument_of_periapsis = numpy.radians(argument_of_periapsis)
  longitude_of_the_ascending_nodes = numpy.radians(longitude_of_the_ascending_nodes)
  true_anomaly = numpy.radians(true_anomaly)

  cos_true_anomaly = numpy.cos(true_anomaly)
  sin_true_anomaly = numpy.sin(true_anomaly)    

  cos_inclination = numpy.cos(inclination)
  sin_inclination = numpy.sin(inclination)    

  cos_arg_per = numpy.cos(argument_of_periapsis)
  sin_arg_per = numpy.sin(argument_of_periapsis)

  cos_long_asc_nodes = numpy.cos(longitude_of_the_ascending_nodes)
  sin_long_asc_nodes = numpy.sin(longitude_of_the_ascending_nodes)

  ### alpha is a unit vector directed along the line of node ###
  alphax = cos_long_asc_nodes*cos_arg_per - sin_long_asc_nodes*sin_arg_per*cos_inclination
  alphay = sin_long_asc_nodes*cos_arg_per + cos_long_asc_nodes*sin_arg_per*cos_inclination
  alphaz = sin_arg_per*sin_inclination
  alpha = [alphax,alphay,alphaz]

  ### beta is a unit vector perpendicular to alpha and the orbital angular momentum vector ###
  betax = -cos_long_asc_nodes*sin_arg_per - sin_long_asc_nodes*cos_arg_per*cos_inclination
  betay = -sin_long_asc_nodes*sin_arg_per + cos_long_asc_nodes*cos_arg_per*cos_inclination
  betaz = cos_arg_per*sin_inclination
  beta = [betax,betay,betaz]

#  print 'alpha',alphax**2+alphay**2+alphaz**2 # For debugging; should be 1
#  print 'beta',betax**2+betay**2+betaz**2 # For debugging; should be 1

  ### Relative position and velocity ###
  separation = semimajor_axis*(1.0 - eccentricity**2)/(1.0 + eccentricity*cos_true_anomaly) # Compute the relative separation
  position_vector = separation*cos_true_anomaly*alpha + separation*sin_true_anomaly*beta
  velocity_tilde = (G*(mass1 + mass2)/(semimajor_axis*(1.0 - eccentricity**2))).sqrt() # Common factor
  velocity_vector = -1.0*velocity_tilde*sin_true_anomaly*alpha + velocity_tilde*(eccentricity + cos_true_anomaly)*beta

  result = Particles(2)
  result[0].mass = mass1
  result[1].mass = mass2
    
  result[1].position = position_vector
  result[1].velocity = velocity_vector
    
  result.move_to_center()
  return result
    
def orbital_elements_from_binary( binary, G=nbody_system.G):
  """ 

  Function that computes orbital elements from given two-particle set. 
  Elements are computed for the second particle in this set and the 
  return values are: mass1, mass2, semimajor axis, eccentricity, 
  true anomaly, inclination, longitude of the ascending nodes and the
  argument of pericenter. All angles are returned in degrees.
  In case of a perfectly circular orbit the true anomaly 
  and argument of pericenter cannot be determined; in this case the 
  return values are 0.0 for both angles. 

  """
  if len(binary)>2:
    raise Exception("expects binary or single part")

  if len(binary)==2:
    mass1=binary[0].mass
    mass2=binary[1].mass
    position = binary[1].position-binary[0].position
    velocity = binary[1].velocity-binary[0].velocity
    total_mass = mass1 + mass2
  else:
    mass1=binary[0].mass
    mass2=0.*mass1
    position = binary[0].position
    velocity = binary[0].velocity
    total_mass = mass1
      
  ### specific energy ###
  specific_energy = (1.0/2.0)*velocity.lengths_squared() - G*total_mass/position.lengths()
  if specific_energy.value_in(units.m**2/(units.s**2)) >= 0.0:
    raise Exception("not a bound orbit")
  semimajor_axis = -G*total_mass/(2.0*specific_energy)
    
  ### specific angular momentum ###    
  specific_angular_momentum = position.cross(velocity)
  specific_angular_momentum_norm = specific_angular_momentum.lengths()    
  specific_angular_momentum_unit=specific_angular_momentum/specific_angular_momentum_norm

  maximum_specific_angular_momentum_norm = G*total_mass/((-2.0*specific_energy).sqrt())
  ell = specific_angular_momentum_norm/maximum_specific_angular_momentum_norm ### specific AM in units of maximum AM

  ### for e = 0 or e nearly 0, ell can be slightly larger than unity due to numerical reasons ###
  ell_epsilon = 1e-15

  completely_or_nearly_circular = False
  
  if ell > 1.0:
    if 1.0 < ell <= ell + ell_epsilon: ### still unity within numerical precision
      print('orbit is completely or nearly circular; in this case the LRL vector cannot be used to reliably obtain the argument of pericenter and true anomaly; the output values of the latter will be set to zero; output e will be e = 0')
      ell = 1.0
      completely_or_nearly_circular = True
    else: ### larger than unity within numerical precision
      raise Exception("angular momentum larger than maximum angular momentum for bound orbit")

  eccentricity = numpy.sqrt(1.0 - ell**2)

  ### Orbital inclination ###
  inclination = numpy.degrees(numpy.arccos(specific_angular_momentum.z/specific_angular_momentum_norm))

  ### Longitude of ascending nodes, with reference direction along x-axis ###
  z_vector = [0.,0.,1.] | units.none
  ascending_node_vector = z_vector.cross(specific_angular_momentum)
  if ascending_node_vector.lengths().number==0:
    ascending_node_vector_unit= numpy.array([1.,0.,0.]) 
  else:
    ascending_node_vector_unit = ascending_node_vector/ascending_node_vector.lengths()
  long_asc_nodes = numpy.degrees(numpy.arctan2(ascending_node_vector_unit[1],ascending_node_vector_unit[0]))

  ### Argument of periapsis and true anomaly, using eccentricity a.k.a. Laplace-Runge-Lenz (LRL) vector ###
  mu = G*total_mass
  position_unit = position/position.lengths()
  e_vector = ( (1.0/mu)*velocity.cross(specific_angular_momentum) - position_unit ) | units.none ### Laplace-Runge-Lenz vector

  if completely_or_nearly_circular == True: ### orbit is completely or nearly circular; in this case the LRL vector cannot be used to reliably obtain the argument of pericenter and true anomaly; the output values of the latter will be set to zero; output e will be e = 0
    arg_per=0.
    true_anomaly=0.
  else:
    e_vector_unit = e_vector/e_vector.lengths()
    
    e_vector_unit_cross_AM_unit = numpy.cross(e_vector_unit,specific_angular_momentum_unit)
    sin_arg_per = ascending_node_vector_unit.dot(e_vector_unit_cross_AM_unit)
    cos_arg_per = e_vector_unit.dot(ascending_node_vector_unit)
    arg_per=numpy.degrees(numpy.arctan2(sin_arg_per,cos_arg_per))

    sin_true_anomaly = position_unit.dot(-1.0*e_vector_unit_cross_AM_unit)
    cos_true_anomaly = position_unit.dot(e_vector_unit)
    true_anomaly=numpy.degrees(numpy.arctan2(sin_true_anomaly,cos_true_anomaly))

  return mass1, mass2, semimajor_axis, eccentricity, true_anomaly, inclination, long_asc_nodes, arg_per
