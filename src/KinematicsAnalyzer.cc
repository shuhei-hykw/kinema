// -*- C++ -*-

#include "KinematicsAnalyzer.hh"

#include <iostream>

#include "Calculator.hh"
#include "charged.hh"
#include "n_charged.hh"
#include "P.hh"
#include "KE.hh"

//______________________________________________________________________________
KinematicsAnalyzer::KinematicsAnalyzer( void )
{
}

//______________________________________________________________________________
KinematicsAnalyzer::~KinematicsAnalyzer( void )
{
  for( auto& p : m_particle_array ){
    delete p; p = 0;
  }
}

//______________________________________________________________________________
bool
KinematicsAnalyzer::DoSearch( void )
{
  int n_charged = m_particle_array.size();
  int n_neutron = 0;
  double Range1, Range2, Range3, errRange1, errRange2, errRange3;
  double theta1, theta2, theta3, errtheta1, errtheta2, errtheta3;
  double phi1, phi2, phi3, errphi1, errphi2, errphi3;

  if( n_charged==3 && n_neutron==0 )//charged particles が３個の場合
    {
      bool answer1 = calc::ThreeCharged( (int)S(0), (int)S(1), (int)S(2),
					 Range(0),  Range(1),  Range(2),
					 RangeE(0), RangeE(1), RangeE(2),
					 Theta(0),  Theta(1),  Theta(2),
					 ThetaE(0), ThetaE(1), ThetaE(2),
					 Phi(0),    Phi(1),    Phi(2),
					 PhiE(0),   PhiE(1),   PhiE(2) );

    }

  else if( n_charged==2 && n_neutron==0 )//charged partibles が２個の場合
    {
      printf("\n\ncharged particles が２本の計算を行います。\n\n");

      double answer2 = two_charged(Range1, Range2, errRange1, errRange2, theta1, theta2, errtheta1, errtheta2, phi1, phi2, errphi1, errphi2);

      printf("\n\n計算が終了しました!");
      getchar();
    }


  else if( n_charged==1 && n_neutron==1 )//charged particles が１個と中世粒子が１個の場合
    {
      printf("\n\ncharged particles が1本と、中性粒子が１個（中性子の場合、２個も含む）の計算を行います。\n\n");

      double answer3 = n_one(Range1, errRange1, theta1, errtheta1, phi1, errphi1);

      printf("\n\n計算が終了しました！");
      getchar();

    }

  else if( n_charged==2 && n_neutron==1 )//charged particles が２個と中世粒子が１個の場合
    {
      printf("\n\ncharged particles が２本と、中性粒子が１個（中性子の場合、２個も含む）の計算を行います。\n\n");

      double answer5 = n_two(Range1, Range2, errRange1, errRange2, theta1, theta2, errtheta1, errtheta2, phi1, phi2, errphi1, errphi2);

      printf("\n\n計算が終了しました！");
      getchar();

    }


  else if( n_charged==3 && n_neutron==1 )//charged particles が３個と中世粒子が１個の場合
    {
      printf("\n\ncharged particles が３本と、中性粒子が１個（中性子の場合、２個も含む）の計算を行います。\n\n");

      double answer5 = n_three(Range1, Range2, Range3, errRange1, errRange2, errRange3, theta1, theta2, theta3, errtheta1, errtheta2, errtheta3, phi1, phi2, phi3, errphi1, errphi2, errphi3);

      printf("\n\n計算が終了しました！");
      getchar();

    }

  return true;
}

//______________________________________________________________________________
void
KinematicsAnalyzer::Print( const std::string& arg ) const
{
  static const std::string func_name("["+ClassName()+"::"+__func__+"()]");

  std::cout << "#D " << func_name << " " << arg << std::endl
	    << " NParticle = " << m_particle_array.size() << std::endl;
  for( const auto& p : m_particle_array ){
    p->Print();
  }
}
