// -*- C++ -*-

#ifndef KINEMATICS_ANALYZER_HH
#define KINEMATICS_ANALYZER_HH

#include <string>
#include <vector>

#include "Particle.hh"

//______________________________________________________________________________
class KinematicsAnalyzer
{
public:
  static KinematicsAnalyzer& GetInstance( void );
  ~KinematicsAnalyzer( void );

private:
  KinematicsAnalyzer( void );
  KinematicsAnalyzer( const KinematicsAnalyzer& );
  KinematicsAnalyzer& operator =( const KinematicsAnalyzer& );

private:
  std::vector<Particle*> m_particle_array;

public:
  void               Add( Particle* p ) { m_particle_array.push_back(p); }
  bool               DoSearch( void );
  Particle*          Get( int i ) const { return m_particle_array.at(i); }
  Particle* operator []( int i ) const { return Get(i); }
  std::string        Name( int i ) const { return Get(i)->Name(); }
  int                A( int i ) const { return Get(i)->A(); }
  int                Z( int i ) const { return Get(i)->Z(); }
  int                S( int i ) const { return Get(i)->S(); }
  double             Mass( int i ) const { return Get(i)->Mass(); }
  double             MassE( int i ) const { return Get(i)->MassE(); }
  double             Energy( int i ) const { return Get(i)->Energy(); }
  double             EnergyE( int i ) const { return Get(i)->EnergyE(); }
  double             Range( int i ) const { return Get(i)->Range(); }
  double             RangeE( int i ) const { return Get(i)->RangeE(); }
  double             Theta( int i ) const { return Get(i)->Theta(); }
  double             ThetaE( int i ) const { return Get(i)->ThetaE(); }
  double             Phi( int i ) const { return Get(i)->Phi(); }
  double             PhiE( int i ) const { return Get(i)->PhiE(); }
  void               Print( const std::string& arg="" ) const;
  static std::string ClassName( void ){ return "KinematicsAnalyzer"; }
};

//______________________________________________________________________________
inline KinematicsAnalyzer&
KinematicsAnalyzer::GetInstance( void )
{
  static KinematicsAnalyzer g_instance;
  return g_instance;
}

#endif
