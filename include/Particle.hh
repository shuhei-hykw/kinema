// -*- C++ -*-

#ifndef PARTICLE_HH
#define PARTICLE_HH

#include <string>

//______________________________________________________________________________
class Particle
{
public:
  // for input
  Particle( std::string name, double s,
	    double range, double range_e,
	    double theta, double theta_e,
	    double phi,   double phi_e );
  // for database
  Particle( std::string name, double a,
	    double z, double s, double mass );
  ~Particle( void );

private:
  bool        m_is_database;
  std::string m_name;
  int         m_a; // A : mass number
  int         m_z; // Z : proton number
  int         m_s; // S : strangeness
  double      m_mass;
  double      m_mass_e;
  double      m_energy;
  double      m_energy_e;
  // input parameters
  double      m_range;
  double      m_range_e;
  double      m_theta;
  double      m_theta_e;
  double      m_phi;
  double      m_phi_e;

public:
  std::string        Name( void ) const { return m_name; }
  int                A( void ) const { return m_a; }
  int                Z( void ) const { return m_z; }
  int                S( void ) const { return m_s; }
  double             Mass( void ) const { return m_mass; }
  double             MassE( void ) const { return m_mass_e; }
  double             Energy( void ) const { return m_energy; }
  double             EnergyE( void ) const { return m_energy_e; }
  double             Range( void ) const { return m_range; }
  double             RangeE( void ) const { return m_range_e; }
  double             Theta( void ) const { return m_theta; }
  double             ThetaE( void ) const { return m_theta_e; }
  double             Phi( void ) const { return m_phi; }
  double             PhiE( void ) const { return m_phi_e; }
  void               Print( const std::string& arg="" ) const;
  static std::string ClassName( void ){ return "Particle"; }
};

#endif
