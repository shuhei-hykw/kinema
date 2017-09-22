// -*- C++ -*-

//このプログラムは、実験結果より得られたRange、theta、phiより、存在が確認されているすべての粒子においてtrackを当てはめて、親粒子の質量を計算する。
//それが誤差の３倍以内で成り立つものを見つけ出すものである。表示する際は、どの粒子を当てはめて、計算したかを、Z,A,Sを用いて、表示しています。
//出力したresult.txtには、反応式を表記して、どのような反応が仮定されたかを見やすくしてあります。

#include <stdio.h>
#include <math.h>

#include "charged.hh"
#include "n_charged.hh"
#include "P.hh"
#include "KE.hh"

//______________________________________________________________________________
int
main( int argc, char *argv[] )
{
  int    charged_kind, n_kind;
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
  printf("「charged particle の数」は、「１」～「３」を入力してください。\n");
  printf("「中性粒子の数」は、「０」か「１」を入力してください。\n\n\n");

  printf("charged particle の数 = ");
  scanf("%d", &charged_kind);

  printf("\n中性粒子の数 = ");
  scanf("%d", &n_kind);



  printf("\n\nここからは、Range、θ、φの入力を行います。\n\n\n");

  printf("Range１= ");
  scanf("%lf", &Range1);
  printf("σRange１= ");
  scanf("%lf", &errRange1);
  printf("θ１ = ");
  scanf("%lf", &theta1);
  printf("σθ１ = ");
  scanf("%lf", &errtheta1);
  printf("φ１= ");
  scanf("%lf", &phi1);
  printf("σφ１= ");
  scanf("%lf", &errphi1);
  printf("\nRange２ = ");
  scanf("%lf", &Range2);
  printf("σRange２ = ");
  scanf("%lf", &errRange2);
  printf("θ２ = ");
  scanf("%lf", &theta2);
  printf("σθ２ = ");
  scanf("%lf", &errtheta2);
  printf("φ２= ");
  scanf("%lf", &phi2);
  printf("σφ２= ");
  scanf("%lf", &errphi2);
  printf("\nRange３= ");
  scanf("%lf", &Range3);
  printf("σRange３= ");
  scanf("%lf", &errRange3);
  printf("θ３ = ");
  scanf("%lf", &theta3);
  printf("σθ３ = ");
  scanf("%lf", &errtheta3);
  printf("φ３= ");
  scanf("%lf", &phi3);
  printf("σφ３= ");
  scanf("%lf", &errphi3);




  //Range1 = 68.5;
  //Range2 = 37.4;
  //Range3 = 0;
  //errRange1 = 1.6;
  //errRange2 = 1.3;
  //errRange3 = 0;
  //theta1 = 117.9;
  //theta2 = 61.2;
  //theta3 = 0;
  //errtheta1 = 3;
  //errtheta2 = 7;
  //errtheta3 = 0;
  //phi1 = 179.9;
  //phi2 = 0.1;
  //phi3 = 0;
  //errphi1 = 0.2;
  //errphi2 = 7;
  //errphi3 = 0;


  //動いているのを確認。
  printf("\n\n与えられた値をもとに、以下に示した内容で計算します。\n\n");


  if((charged_kind==3)&&(n_kind==0))//charged particles が３個の場合
    {
      printf("\n\ncharged particles が３本の計算を行います。\n\n");

      double answer1 = three_charged(Range1, Range2, Range3, errRange1, errRange2, errRange3, theta1, theta2, theta3, errtheta1, errtheta2, errtheta3, phi1, phi2, phi3, errphi1, errphi2, errphi3);

      printf("\n\n計算が終了しました!");
      getchar();
    }

  else if((charged_kind==2)&&(n_kind==0))//charged partibles が２個の場合
    {
      printf("\n\ncharged particles が２本の計算を行います。\n\n");

      double answer2 = two_charged(Range1, Range2, errRange1, errRange2, theta1, theta2, errtheta1, errtheta2, phi1, phi2, errphi1, errphi2);

      printf("\n\n計算が終了しました!");
      getchar();
    }


  else if((charged_kind==1)&&(n_kind==1))//charged particles が１個と中世粒子が１個の場合
    {
      printf("\n\ncharged particles が1本と、中性粒子が１個（中性子の場合、２個も含む）の計算を行います。\n\n");

      double answer3 = n_one(Range1, errRange1, theta1, errtheta1, phi1, errphi1);

      printf("\n\n計算が終了しました！");
      getchar();

    }


  else if((charged_kind==2)&&(n_kind==1))//charged particles が２個と中世粒子が１個の場合
    {
      printf("\n\ncharged particles が２本と、中性粒子が１個（中性子の場合、２個も含む）の計算を行います。\n\n");

      double answer5 = n_two(Range1, Range2, errRange1, errRange2, theta1, theta2, errtheta1, errtheta2, phi1, phi2, errphi1, errphi2);

      printf("\n\n計算が終了しました！");
      getchar();

    }


  else if((charged_kind==3)&&(n_kind==1))//charged particles が３個と中世粒子が１個の場合
    {
      printf("\n\ncharged particles が３本と、中性粒子が１個（中性子の場合、２個も含む）の計算を行います。\n\n");

      double answer5 = n_three(Range1, Range2, Range3, errRange1, errRange2, errRange3, theta1, theta2, theta3, errtheta1, errtheta2, errtheta3, phi1, phi2, phi3, errphi1, errphi2, errphi3);

      printf("\n\n計算が終了しました！");
      getchar();

    }

}
