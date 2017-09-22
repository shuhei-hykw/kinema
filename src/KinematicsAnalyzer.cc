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

  //キーボードから実験値を入力する。
  printf("\n\n");
  printf("（１）このプログラムは、キーボードより打ち込まれたtrackのRange、θ、φから\n");
  printf("　　　当てはまる崩壊過程を見つけだすプログラムです。\n");
  printf("（２）このプログラムでは、以下の崩壊パターンより、数値入力にて選択し、\n");
  printf("      使用することができます。\n");
  printf("（３）あらかじめ決めておいた、粒子１から粒子３の順に打ち込んでください。\n");
  printf("（４）計算が完了しましたら、このウインドウに「計算されました！」と表記されます。");
  printf("　　　表記されるまで、しばらくお待ちください。\n");
  printf("（５）打ち込む数値は、半角数字で入力してください。\n");
  printf("（６）数値は、打ち込むごとにenterキーを押してください。\n\n\n");

  printf("はじめに,崩壊パターンを決定します。\n");
  printf("「charged particle の数」は、「１」〜「３」を入力してください。\n");
  printf("「中性粒子の数」は、「０」か「１」を入力してください。\n\n\n");

  printf("\n\nここからは、Range、θ、φの入力を行います。\n\n\n");

  //動いているのを確認。
  printf("\n\n与えられた値をもとに、以下に示した内容で計算します。\n\n");

  if( n_charged==3 && n_neutron==0 )//charged particles が３個の場合
    {
      double answer1 = calc::ThreeCharged( (int)S(0), (int)S(1), (int)S(2),
					   Range(0),  Range(1),  Range(2),
					   RangeE(0), RangeE(1), RangeE(2),
					   Theta(0),  Theta(1),  Theta(2),
					   ThetaE(0), ThetaE(1), ThetaE(2),
					   Phi(0),    Phi(1),    Phi(2),
					   PhiE(0),   PhiE(1),   PhiE(2) );

      printf("\n\n計算が終了しました!");
      getchar();
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
