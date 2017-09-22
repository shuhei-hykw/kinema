// -*- C++ -*-

#include "Particle.hh"

#include <iomanip>
#include <iostream>

//______________________________________________________________________________
Particle::Particle( std::string name, double s,
		    double range, double range_e,
		    double theta, double theta_e,
		    double phi,   double phi_e )
  : m_is_database(false),
    m_name(name), m_a(-1.), m_z(-1.), m_s(s),
    m_mass(-1.), m_mass_e(-1.),
    m_energy(-1.), m_energy_e(-1.),
    m_range(range), m_range_e(range_e),
    m_theta(theta), m_theta_e(theta_e),
    m_phi(phi),     m_phi_e(phi_e)
{
}

//______________________________________________________________________________
Particle::Particle( std::string name, double a,
		    double z, double s, double mass )
  : m_is_database(true),
    m_name(name), m_a(a), m_z(z), m_s(s),
    m_mass(mass), m_mass_e(0.),
    m_energy(0.), m_energy_e(0.),
    m_range(0.), m_range_e(0.),
    m_theta(0.), m_theta_e(0.),
    m_phi(0.),     m_phi_e(0.)
{
}

//______________________________________________________________________________
Particle::~Particle( void )
{
}

//______________________________________________________________________________
void
Particle::Print( const std::string& arg ) const
{
  if( m_is_database ){
    std::cout << " name = " << std::setw(8) << std::left << m_name
	      << " A = " << std::setw(2) << m_a
	      << " Z = " << std::setw(2) << m_z
	      << " S = " << std::setw(2) << m_s
	      << " Mass = " << std::setw(6) << m_mass
	      << std::endl;
  } else {
    std::cout << " name = " << std::setw(8) << std::left << m_name
	      << " A = " << std::setw(2) << m_a
	      << " Z = " << std::setw(2) << m_z
	      << " S = " << std::setw(2) << m_s
	      << " Mass = " << std::setw(6) << m_mass
	      << " +/- " << std::setw(6) << m_mass_e
	      << " Energy = " << std::setw(6) << m_energy
	      << " +/- " << std::setw(6) << m_energy_e
	      << " Range = " << std::setw(6) << m_range
	      << " +/- " << std::setw(6) << m_range_e
	      << " Theta = " << std::setw(6) << m_theta
	      << " +/- " << std::setw(6) << m_theta_e
	      << " Phi = " << std::setw(6) << m_phi
	      << " +/- " << std::setw(6) << m_phi_e
	      << std::endl;
  }
}
