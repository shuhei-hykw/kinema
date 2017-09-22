// -*- C++ -*-

#include <cstdio>
#include <TMath.h>

#include "charged.hh"
#include "P.hh"
#include "KE.hh"

double three_charged(double Range1,double Range2,double Range3,double errRange1,double errRange2,double errRange3,
		     double theta1,double theta2,double theta3,double errtheta1,double errtheta2,double errtheta3,
		     double phi1,double phi2,double phi3,double errphi1,double errphi2,double errphi3)

{


  //Rangeのコピー
  double Range1_1,Range1_1_1,Range2_1,Range2_1_1,Range3_1,Range3_1_1;
  Range1_1 = Range1;
  Range1_1_1 = Range1;
  Range2_1 = Range2;
  Range2_1_1 = Range2;
  Range3_1 = Range3;
  Range3_1_1 = Range3;


  double Range1_2,Range1_3,Range2_2,Range2_3,Range3_2,Range3_3;


  int Z1,Z2,Z3,A1,A2,A3,S1,S2,S3;//粒子の原子番号、質量数、ストレンジ数（S=3は、中間子など）
  int Z1_1,Z2_1,Z3_1,A1_1,A2_1,A3_1,S1_1,S2_1,S3_1;
  double Mass1,Mass2,Mass3;//粒子の質量
  double pi = 6*asin(0.5);//円周率の定義

  int SS,ZZ,AA;//親粒子の種類
  double Estimated_Mass,Mass_gap;//推定されるsingleおよびdoubleの質量、推定される親粒子と計算された質量との差



  double M[4][19][8] ={{{938.272,0.,0.,0.,0.,0.,0.,0.},//考えられる粒子の質量表
			{1875.613,0.,0.,0.,0.,0.,0.,0.},
			{2808.922,2808.392,0.,0.,0.,0.,0.,0.},
			{0.,3727.380,0.,0.,0.,0.,0.,0.},
			{0.,4667.845,4667.624,0.,0.,0.,0.,0.},
			{0.,5605.541,5601.520,5605.298,0.,0.,0.,0.},
			{0.,6545.550,6533.836,6534.186,0.,0.,0.,0.},
			{0.,7482.542,7471.369,7454.852,7472.321,0.,0.,0.},
			{0.,0.,8406.871,8392.753,8393.310,8409.295,0.,0.},
			{0.,0.,0.,9325.507,9324.440,9327.580,0.,0.},
			{0.,0.,10285.846,10264.570,10252.550,10254.022,0.,0.},
			{0.,0.,0.,11200.922,11188.746,11174.866,11191.693,0.},
			{0.,0.,0.,12142.542,12123.434,12109.485,12111.195,12128.443},
			{0.,0.,0.,0.,13062.023,13040.874,13040.207,13044.841},
			{0.,0.,0.,0.,0.,13979.222,13968.939,13971.182},
			{0.,0.,0.,0.,0.,14914.536,14906.014,14895.084},
			{0.,0.,0.,0.,0.,0.,0.,15830.506},
			{139.570,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.}},
		       {{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{2991.166,0.,0.,0.,0.,0.,0.,0.},
			{3922.565,3921.685,0.,0.,0.,0.,0.,0.},
			{0.,4839.943,0.,0.,0.,0.,0.,0.},
			{0.,5779.348,5778.807,0.,0.,0.,0.,0.},
			{0.,0.,6711.623,6715.821,0.,0.,0.,0.},
			{0.,7654.073,7642.719,7643.029,0.,0.,0.,0.},
			{0.,0.,8578.552,8563.825,8579.714,0.,0.,0.},
			{0.,0.,0.,9499.326,9500.103,0.,0.,0.},
			{0.,0.,0.,0.,10429.883,0.,0.,0.},
			{0.,0.,0.,0.,11356.863,11358.905,0.,0.},
			{0.,0.,0.,0.,12293.059,12278.859,0.,0.},
			{0.,0.,0.,0.,0.,13212.998,13214.708,0.},
			{0.,0.,0.,0.,0.,0.,14142.300,0.},
			{0.,0.,0.,0.,0.,0.,0.,15074.365},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,16931.689},
			{0.,0.,0.,0.,0.,0.,0.,0.}},
		       {{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{4106.719,0.,0.,0.,0.,0.,0.,0.},
			{5036.208,5034.978,0.,0.,0.,0.,0.,0.},
			{0.,5952.506,0.,0.,0.,0.,0.,0.},
			{0.,6890.851,6889.990,0.,0.,0.,0.,0.},
			{0.,0.,7821.726,7826.344,0.,0.,0.,0.},
			{0.,8762.596,8751.602,8751.872,0.,0.,0.,0.},
			{0.,0.,9685.735,9672.798,9687.107,0.,0.,0.},
			{0.,0.,0.,10605.899,10606.896,0.,0.,0.},
			{0.,0.,0.,11538.653,11535.326,0.,0.,0.},
			{0.,0.,0.,0.,12461.176,12463.788,0.,0.},
			{0.,0.,0.,0.,13397.372,13382.852,0.,0.},
			{0.,0.,0.,0.,0.,14316.511,14318.221,0.},
			{0.,0.,0.,0.,0.,15247.900,15244.393,0.},
			{0.,0.,0.,0.,0.,0.,0.,16177.548},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,18032.872}},
		       {{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,12496.576,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,14361.917,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,16216.794,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.}}};


  char sort[4][19][8][10] = {{{"p","0.","0.","0.","0.","0.","0.","0."},//上の配列に対応した、考えられる粒子の表記（３次元文字配列）
			      {"d","0.","0.","0.","0.","0.","0.","0."},
			      {"t","He3","0.","0.","0.","0.","0.","0."},
			      {"0","He4","0.","0.","0.","0.","0.","0."},
			      {"0.","He5","Li5","0.","0.","0.","0.","0."},
			      {"0.","He6","Li6","Be6","0.","0.","0.","0."},
			      {"0.","He7","Li7","Be7","0.","0.","0.","0."},
			      {"0.","He8","Li8","Be8","B8","0.","0.","0."},
			      {"0.","0.","Li9","Be9","B9","C9","0.","0."},
			      {"0.","0.","0.","Be10","B10","C10","0.","0."},
			      {"0.","0.","Li11","Be11","B11","C11","0.","0."},
			      {"0.","0.","0.","Be12","B12","C12","N12","0."},
			      {"0.","0.","0.","Be13","B13","C13","N13","O13"},
			      {"0.","0.","0.","0.","B14","C14","N14","O14"},
			      {"0.","0.","0.","0.","0.","C15","N15","O15"},
			      {"0.","0.","0.","0.","0.","C16","N16","O16"},
			      {"0.","0.","0.","0.","0.","0.","0.","O17"},
			      {"π-","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."}},
			     {{"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"H3L","0.","0.","0.","0.","0.","0.","0."},
			      {"H4L","He4L","0.","0.","0.","0.","0.","0."},
			      {"0.","He5L","0.","0.","0.","0.","0.","0."},
			      {"0.","He6L","Li6L","0.","0.","0.","0.","0."},
			      {"0.","0.","Li7L","Be7L","0.","0.","0.","0."},
			      {"0.","He8L","Li8L","Be8L","0.","0.","0.","0."},
			      {"0.","0.","Li9L","Be9L","B9L","0.","0.","0."},
			      {"0.","0.","0.","Be10L","B10L","0.","0.","0."},
			      {"0.","0.","0.","0.","B11L","0.","0.","0."},
			      {"0.","0.","0.","0.","B12L","C12L","0.","0."},
			      {"0.","0.","0.","0.","B13L","C13L","0.","0."},
			      {"0.","0.","0.","0.","0.","C14L","N14L","0."},
			      {"0.","0.","0.","0.","0.","0.","N15L","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","O16L"},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","O18L"},
			      {"0.","0.","0.","0.","0.","0.","0.","0."}},
			     {{"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"H4LL","0.","0.","0.","0.","0.","0.","0."},
			      {"H5LL","He5LL","0.","0.","0.","0.","0.","0."},
			      {"0.","He6LL","0.","0.","0.","0.","0.","0."},
			      {"0.","He7LL","Li7LL","0.","0.","0.","0.","0."},
			      {"0.","0.","Li8LL","Be8LL","0.","0.","0.","0."},
			      {"0.","He9LL","Li9LL","Be9LL","0.","0.","0.","0."},
			      {"0.","0.","Li10LL","Be10LL","B10LL","0.","0.","0."},
			      {"0.","0.","0.","Be11LL","B11LL","0.","0.","0."},
			      {"0.","0.","0.","Be12LL","B12LL","0.","0.","0."},
			      {"0.","0.","0.","0.","B13LL","C13LL","0.","0."},
			      {"0.","0.","0.","0.","B14LL","C14LL","0.","0."},
			      {"0.","0.","0.","0.","0.","C15LL","N15LL","0."},
			      {"0.","0.","0.","0.","0.","C16LL","N16LL","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","O17LL"},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","O19LL"}},
			     {{"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","Ξ- + C12","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","Ξ- + N14","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","Ξ- + O16","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."},
			      {"0.","0.","0.","0.","0.","0.","0.","0."}}};

  printf("\n\n最後に、親粒子と娘粒子のストレンジ数を決定します。\n");
  printf("ここでは、「０」～「３」を入力します。\n");
  printf("通常原子核（π-も含む ）＝  「０」\n");
  printf("    single hyper     ＝  「1」\n");
  printf("    double hyper     ＝  「2」\n");
  printf("    Ξ- + C(N,O)      ＝  「3」\n\n");

  printf("親粒子　= ");
  scanf("%d",&SS);
  printf("track #1 = ");
  scanf("%d",&S1);
  printf("track #2 = ");
  scanf("%d",&S2);
  printf("track #3 = ");
  scanf("%d",&S3);




  for(ZZ=1;ZZ<=8;ZZ++){
    for(AA=1;AA<=19;AA++){
      for(Z1=1;Z1<=8;Z1++){
	for(A1=1;A1<=19;A1++){
	  for(Z2=1;Z2<=8;Z2++){
	    for(A2=1;A2<=19;A2++){
	      for(Z3=1;Z3<=8;Z3++){
		for(A3=1;A3<=19;A3++){     //すべての粒子に対してRange,theta,phiを当てはめていく。




		  Mass1 = M[S1][A1-1][Z1-1];
		  Mass2 = M[S2][A2-1][Z2-1];
		  Mass3 = M[S3][A3-1][Z3-1];//粒子の質量
		  Estimated_Mass = M[SS][AA-1][ZZ-1];//推定されるsingleおよびdouble、Ξが吸収された原子核の質量




		  //反応式の構成のための原子番号、質量数のコピー
		  Z1_1=Z1;
		  Z2_1=Z2;
		  Z3_1=Z3;
		  A1_1=A1;
		  A2_1=A2;
		  A3_1=A3;
		  S1_1=S1;
		  S2_1=S2;
		  S3_1=S3;


		  double Q = (Estimated_Mass - Mass1 - Mass2 - Mass3);//Q値

		  if(Estimated_Mass <= 0.)//親粒子の質量が０になってしまうところは飛ばす。
		    {
		      continue;
		    }


		  if((Mass1 <= 0.)||(Mass2 <= 0.)||(Mass3 <= 0.))//粒子の質量が０になってしまうところは飛ばす。
		    {
		      continue;
		    }

		  if(Q <= 0.0000)//Q値が正になる反応式以外は、飛ばす。
		    {
		      continue;
		    }

		  double kinetic_energy1 = KE(Mass1,Range1,Z1,A1,S1) ;//粒子１の運動エネルギー
		  double momentum1 = P(Mass1,Range1,Z1,A1,S1);//粒子１の運動量
		  Range1_2 = Range1_1 + errRange1;
		  double kinetic_energy1_up = KE(Mass1,Range1_2,Z1,A1,S1);
		  double momentum1_up = P(Mass1,Range1_2,Z1,A1,S1);
		  Range1_3 = Range1_1_1 - errRange1;
		  double kinetic_energy1_down = KE(Mass1,Range1_3,Z1,A1,S1);
		  double momentum1_down = P(Mass1,Range1_3,Z1,A1,S1);
		  double err_kinetic_energy1 = fabs(kinetic_energy1_up - kinetic_energy1_down)/2;//粒子１のRangeから考えられる運動エネルギーの誤差
		  double err_momentum1 = fabs(momentum1_up - momentum1_down)/2;//粒子１のRangeから考えられる運動量の誤差



		  double kinetic_energy2 = KE(Mass2,Range2,Z2,A2,S2);//粒子２の運動エネルギー
		  double momentum2 = P(Mass2,Range2,Z2,A2,S2);//粒子２の運動量
		  Range2_2 = Range2_1 + errRange2;
		  double kinetic_energy2_up = KE(Mass2,Range2_2,Z2,A2,S2);
		  double momentum2_up = P(Mass2,Range2_2,Z2,A2,S2);
		  Range2_3 = Range2_1_1 - errRange2;
		  double kinetic_energy2_down = KE(Mass2,Range2_3,Z2,A2,S2);
		  double momentum2_down = P(Mass2,Range2_3,Z2,A2,S2);
		  double err_kinetic_energy2 = fabs(kinetic_energy2_up - kinetic_energy2_down)/2;//粒子２のRangeから考えられる運動エネルギーの誤差
		  double err_momentum2 = fabs(momentum2_up - momentum2_down)/2;//粒子２のRangeから考えられる運動量の誤差


		  double kinetic_energy3 = KE(Mass3,Range3,Z3,A3,S3);//粒子３の運動エネルギー
		  double momentum3 = P(Mass3,Range3,Z3,A3,S3);//粒子３の運動量
		  Range3_2 = Range3_1 + errRange3;
		  double kinetic_energy3_up = KE(Mass3,Range3_2,Z3,A3,S3);
		  double momentum3_up = P(Mass3,Range3_2,Z3,A3,S3);
		  Range3_3 = Range3_1_1 - errRange3;
		  double kinetic_energy3_down = KE(Mass3,Range3_3,Z3,A3,S3);
		  double momentum3_down = P(Mass3,Range3_3,Z3,A3,S3);
		  double err_kinetic_energy3 = fabs(kinetic_energy3_up - kinetic_energy3_down)/2;//粒子３のRangeから考えられる運動エネルギーの誤差
		  double err_momentum3 = fabs(momentum3_up - momentum3_down)/2;//粒子３のRangeから考えられる運動量の誤差


		  if((Z1==1)&&(A1==18)&&(S1==0))
		    {
		      Z1_1=-1;
		      A1_1=0;
		    }

		  if((Z2==1)&&(A2==18)&&(S2==0))
		    {
		      Z2_1=-1;
		      A2_1=0;
		    }

		  if((Z3==1)&&(A3==18)&&(S3==0))
		    {
		      Z3_1=-1;
		      A3_1=0;
		    }



		  if((ZZ!=Z1_1+Z2_1+Z3_1)||(AA!=A1_1+A2_1+A3_1))
		    {
		      continue;
		    }



		  double p1x = momentum1*sin(theta1*pi/180)*cos(phi1*pi/180);//粒子１の運動量の各成分
		  double p1y = momentum1*sin(theta1*pi/180)*sin(phi1*pi/180);
		  double p1z = momentum1*cos(theta1*pi/180);

		  double p2x = momentum2*sin(theta2*pi/180)*cos(phi2*pi/180);//粒子２の運動量の各成分
		  double p2y = momentum2*sin(theta2*pi/180)*sin(phi2*pi/180);
		  double p2z = momentum2*cos(theta2*pi/180);

		  double p3x = momentum3*sin(theta3*pi/180)*cos(phi3*pi/180);//粒子３の運動量の各成分
		  double p3y = momentum3*sin(theta3*pi/180)*sin(phi3*pi/180);
		  double p3z = momentum3*cos(theta3*pi/180);


		  //各粒子の運動量ベクトルの角度をA,B,Cとした時のcos成分
		  //粒子１と粒子２の角度A
		  //粒子２と粒子３の角度B
		  //粒子１と粒子３の角度C
		  double cosA = (p1x*p2x+p1y*p2y+p1z*p2z)/((sqrt(p1x*p1x+p1y*p1y+p1z*p1z))*(sqrt(p2x*p2x+p2y*p2y+p2z*p2z)));
		  double cosB = (p2x*p3x+p2y*p3y+p2z*p3z)/((sqrt(p2x*p2x+p2y*p2y+p2z*p2z))*(sqrt(p3x*p3x+p3y*p3y+p3z*p3z)));
		  double cosC = (p1x*p3x+p1y*p3y+p1z*p3z)/((sqrt(p1x*p1x+p1y*p1y+p1z*p1z))*(sqrt(p3x*p3x+p3y*p3y+p3z*p3z)));


		  double errp1x,errp1y,errp1z,errp2x,errp2y,errp2z,errp3x,errp3y,errp3z;//各粒子の運動量成分の誤差
		  double errp1,errp2,errp3;//各粒子の運動量の誤差
		  double errcosA,errcosB,errcosC;//各粒子の運動量ベクトルの角度をA,B,Cとした時のcos成分の誤差



		  //各粒子の運動量成分の誤差
		  errp1x = sqrt(pow((sin(theta1*pi/180)*cos(phi1*pi/180)),2)*pow(err_momentum1,2)
				+pow(momentum1*cos(theta1*pi/180)*cos(phi1*pi/180),2)*pow(errtheta1*pi/180,2)+pow(momentum1*sin(theta1*pi/180)*sin(phi1*pi/180),2)*pow(errphi1*pi/180,2));
		  errp1y = sqrt(pow((sin(theta1*pi/180)*sin(phi1*pi/180)),2)*pow(err_momentum1,2)
				+pow(momentum1*cos(theta1*pi/180)*sin(phi1*pi/180),2)*pow(errtheta1*pi/180,2)+pow(momentum1*sin(theta1*pi/180)*cos(phi1*pi/180),2)*pow(errphi1*pi/180,2));
		  errp1z = sqrt(pow((cos(theta1*pi/180)),2)*pow(err_momentum1,2)+pow(momentum1*sin(theta1*pi/180),2)*pow(errtheta1*pi/180,2));

		  errp2x = sqrt(pow((sin(theta2*pi/180)*cos(phi2*pi/180)),2)*pow(err_momentum2,2)
				+pow(momentum2*cos(theta2*pi/180)*cos(phi2*pi/180),2)*pow(errtheta2*pi/180,2)+pow(momentum2*sin(theta2*pi/180)*sin(phi2*pi/180),2)*pow(errphi2*pi/180,2));
		  errp2y = sqrt(pow((sin(theta2*pi/180)*sin(phi2*pi/180)),2)*pow(err_momentum2,2)
				+pow(momentum2*cos(theta2*pi/180)*sin(phi2*pi/180),2)*pow(errtheta2*pi/180,2)+pow(momentum2*sin(theta2*pi/180)*cos(phi2*pi/180),2)*pow(errphi2*pi/180,2));
		  errp2z = sqrt(pow((cos(theta2*pi/180)),2)*pow(err_momentum2,2)+pow(momentum2*sin(theta2*pi/180),2)*pow(errtheta2*pi/180,2));

		  errp3x = sqrt(pow((sin(theta3*pi/180)*cos(phi3*pi/180)),2)*pow(err_momentum3,2)
				+pow(momentum3*cos(theta3*pi/180)*cos(phi3*pi/180),2)*pow(errtheta3*pi/180,2)+pow(momentum3*sin(theta3*pi/180)*sin(phi3*pi/180),2)*pow(errphi3*pi/180,2));
		  errp3y = sqrt(pow((sin(theta3*pi/180)*sin(phi3*pi/180)),2)*pow(err_momentum3,2)
				+pow(momentum3*cos(theta3*pi/180)*sin(phi3*pi/180),2)*pow(errtheta3*pi/180,2)+pow(momentum3*sin(theta3*pi/180)*cos(phi3*pi/180),2)*pow(errphi3*pi/180,2));
		  errp3z = sqrt(pow((cos(theta3*pi/180)),2)*pow(err_momentum3,2)+pow(momentum3*sin(theta3*pi/180),2)*pow(errtheta3*pi/180,2));


		  //各粒子の運動量の誤差
		  errp1 = sqrt(pow(p1x/sqrt(pow(p1x,2)+pow(p1y,2)+pow(p1z,2)),2)*pow(errp1x,2)
			       +pow(p1y/sqrt(pow(p1x,2)+pow(p1y,2)+pow(p1z,2)),2)*pow(errp1y,2)+pow(p1z/sqrt(pow(p1x,2)+pow(p1y,2)+pow(p1z,2)),2)*pow(errp1z,2));
		  errp2 = sqrt(pow(p2x/sqrt(pow(p2x,2)+pow(p2y,2)+pow(p2z,2)),2)*pow(errp2x,2)
			       +pow(p2y/sqrt(pow(p2x,2)+pow(p2y,2)+pow(p2z,2)),2)*pow(errp2y,2)+pow(p2z/sqrt(pow(p2x,2)+pow(p2y,2)+pow(p2z,2)),2)*pow(errp2z,2));
		  errp3 = sqrt(pow(p3x/sqrt(pow(p3x,2)+pow(p3y,2)+pow(p3z,2)),2)*pow(errp3x,2)
			       +pow(p3y/sqrt(pow(p3x,2)+pow(p3y,2)+pow(p3z,2)),2)*pow(errp3y,2)+pow(p3z/sqrt(pow(p3x,2)+pow(p3y,2)+pow(p3z,2)),2)*pow(errp3z,2));


		  //各粒子の運動量ベクトルの角度をA,B,Cとした時のcos成分の誤差
		  errcosA = sqrt(pow((p2x*sqrt((pow(p1x,2)+pow(p1y,2)+pow(p1z,2)))-p1x*(p1x*p2x+p1y*p2y+p1z*p2z)*pow(pow(p1x,2)+pow(p1y,2)+pow(p1z,2),-0.5))
				     /((pow(p1x,2)+pow(p1y,2)+pow(p1z,2))*sqrt(pow(p2x,2)+pow(p2y,2)+pow(p2z,2))),2)*pow(errp1x,2)+
				 pow((p1x*sqrt((pow(p2x,2)+pow(p2y,2)+pow(p2z,2)))-p2x*(p1x*p2x+p1y*p2y+p1z*p2z)*pow(pow(p2x,2)+pow(p2y,2)+pow(p2z,2),-0.5))
				     /((pow(p2x,2)+pow(p2y,2)+pow(p2z,2))*sqrt(pow(p1x,2)+pow(p1y,2)+pow(p1z,2))),2)*pow(errp2x,2)+
				 pow((p2y*sqrt((pow(p1x,2)+pow(p1y,2)+pow(p1z,2)))-p1y*(p1x*p2x+p1y*p2y+p1z*p2z)*pow(pow(p1x,2)+pow(p1y,2)+pow(p1z,2),-0.5))
				     /((pow(p1x,2)+pow(p1y,2)+pow(p1z,2))*sqrt(pow(p2x,2)+pow(p2y,2)+pow(p2z,2))),2)*pow(errp1y,2)+
				 pow((p1y*sqrt((pow(p2x,2)+pow(p2y,2)+pow(p2z,2)))-p2y*(p1x*p2x+p1y*p2y+p1z*p2z)*pow(pow(p2x,2)+pow(p2y,2)+pow(p2z,2),-0.5))
				     /((pow(p2x,2)+pow(p2y,2)+pow(p2z,2))*sqrt(pow(p1x,2)+pow(p1y,2)+pow(p1z,2))),2)*pow(errp2y,2)+
				 pow((p2z*sqrt((pow(p1x,2)+pow(p1y,2)+pow(p1z,2)))-p1z*(p1x*p2x+p1y*p2y+p1z*p2z)*pow(pow(p1x,2)+pow(p1y,2)+pow(p1z,2),-0.5))
				     /((pow(p1x,2)+pow(p1y,2)+pow(p1z,2))*sqrt(pow(p2x,2)+pow(p2y,2)+pow(p2z,2))),2)*pow(errp1z,2)+
				 pow((p1z*sqrt((pow(p2x,2)+pow(p2y,2)+pow(p2z,2)))-p2z*(p1x*p2x+p1y*p2y+p1z*p2z)*pow(pow(p2x,2)+pow(p2y,2)+pow(p2z,2),-0.5))
				     /((pow(p2x,2)+pow(p2y,2)+pow(p2z,2))*sqrt(pow(p1x,2)+pow(p1y,2)+pow(p1z,2))),2)*pow(errp2z,2));

		  errcosB = sqrt(pow((p3x*sqrt((pow(p2x,2)+pow(p2y,2)+pow(p2z,2)))-p2x*(p2x*p3x+p2y*p3y+p2z*p3z)*pow(pow(p2x,2)+pow(p2y,2)+pow(p2z,2),-0.5))
				     /((pow(p2x,2)+pow(p2y,2)+pow(p2z,2))*sqrt(pow(p3x,2)+pow(p3y,2)+pow(p3z,2))),2)*pow(errp2x,2)+
				 pow((p2x*sqrt((pow(p3x,2)+pow(p3y,2)+pow(p3z,2)))-p3x*(p2x*p3x+p2y*p3y+p2z*p3z)*pow(pow(p3x,2)+pow(p3y,2)+pow(p3z,2),-0.5))
				     /((pow(p3x,2)+pow(p3y,2)+pow(p3z,2))*sqrt(pow(p2x,2)+pow(p2y,2)+pow(p2z,2))),2)*pow(errp3x,2)+
				 pow((p3y*sqrt((pow(p2x,2)+pow(p2y,2)+pow(p2z,2)))-p2y*(p2x*p3x+p2y*p3y+p2z*p3z)*pow(pow(p2x,2)+pow(p2y,2)+pow(p2z,2),-0.5))
				     /((pow(p2x,2)+pow(p2y,2)+pow(p2z,2))*sqrt(pow(p3x,2)+pow(p3y,2)+pow(p3z,2))),2)*pow(errp2y,2)+
				 pow((p2y*sqrt((pow(p3x,2)+pow(p3y,2)+pow(p3z,2)))-p3y*(p2x*p3x+p2y*p3y+p2z*p3z)*pow(pow(p3x,2)+pow(p3y,2)+pow(p3z,2),-0.5))
				     /((pow(p3x,2)+pow(p3y,2)+pow(p3z,2))*sqrt(pow(p2x,2)+pow(p2y,2)+pow(p2z,2))),2)*pow(errp3y,2)+
				 pow((p3z*sqrt((pow(p2x,2)+pow(p2y,2)+pow(p2z,2)))-p2z*(p2x*p3x+p2y*p3y+p2z*p3z)*pow(pow(p2x,2)+pow(p2y,2)+pow(p2z,2),-0.5))
				     /((pow(p2x,2)+pow(p2y,2)+pow(p2z,2))*sqrt(pow(p3x,2)+pow(p3y,2)+pow(p3z,2))),2)*pow(errp2z,2)+
				 pow((p2z*sqrt((pow(p3x,2)+pow(p3y,2)+pow(p3z,2)))-p3z*(p2x*p3x+p2y*p3y+p2z*p3z)*pow(pow(p3x,2)+pow(p3y,2)+pow(p3z,2),-0.5))
				     /((pow(p3x,2)+pow(p3y,2)+pow(p3z,2))*sqrt(pow(p2x,2)+pow(p2y,2)+pow(p2z,2))),2)*pow(errp3z,2));

		  errcosC = sqrt(pow((p3x*sqrt((pow(p1x,2)+pow(p1y,2)+pow(p1z,2)))-p1x*(p1x*p3x+p1y*p3y+p1z*p3z)*pow(pow(p1x,2)+pow(p1y,2)+pow(p1z,2),-0.5))
				     /((pow(p1x,2)+pow(p1y,2)+pow(p1z,2))*sqrt(pow(p3x,2)+pow(p3y,2)+pow(p3z,2))),2)*pow(errp1x,2)+
				 pow((p1x*sqrt((pow(p3x,2)+pow(p3y,2)+pow(p3z,2)))-p3x*(p1x*p3x+p1y*p3y+p1z*p3z)*pow(pow(p3x,2)+pow(p3y,2)+pow(p3z,2),-0.5))
				     /((pow(p3x,2)+pow(p3y,2)+pow(p3z,2))*sqrt(pow(p1x,2)+pow(p1y,2)+pow(p1z,2))),2)*pow(errp3x,2)+
				 pow((p3y*sqrt((pow(p1x,2)+pow(p1y,2)+pow(p1z,2)))-p1y*(p1x*p3x+p1y*p3y+p1z*p3z)*pow(pow(p1x,2)+pow(p1y,2)+pow(p1z,2),-0.5))
				     /((pow(p1x,2)+pow(p1y,2)+pow(p1z,2))*sqrt(pow(p3x,2)+pow(p3y,2)+pow(p3z,2))),2)*pow(errp1y,2)+
				 pow((p1y*sqrt((pow(p3x,2)+pow(p3y,2)+pow(p3z,2)))-p3y*(p1x*p3x+p1y*p3y+p1z*p3z)*pow(pow(p3x,2)+pow(p3y,2)+pow(p3z,2),-0.5))
				     /((pow(p3x,2)+pow(p3y,2)+pow(p3z,2))*sqrt(pow(p1x,2)+pow(p1y,2)+pow(p1z,2))),2)*pow(errp3y,2)+
				 pow((p3z*sqrt((pow(p1x,2)+pow(p1y,2)+pow(p1z,2)))-p1z*(p1x*p3x+p1y*p3y+p1z*p3z)*pow(pow(p1x,2)+pow(p1y,2)+pow(p1z,2),-0.5))
				     /((pow(p1x,2)+pow(p1y,2)+pow(p1z,2))*sqrt(pow(p3x,2)+pow(p3y,2)+pow(p3z,2))),2)*pow(errp1z,2)+
				 pow((p1z*sqrt((pow(p3x,2)+pow(p3y,2)+pow(p3z,2)))-p3z*(p1x*p3x+p1y*p3y+p1z*p3z)*pow(pow(p3x,2)+pow(p3y,2)+pow(p3z,2),-0.5))
				     /((pow(p3x,2)+pow(p3y,2)+pow(p3z,2))*sqrt(pow(p1x,2)+pow(p1y,2)+pow(p1z,2))),2)*pow(errp3z,2));



		  //計算された親粒子の質量の誤差
		  double err_Caluculated_Mass = sqrt(pow(0.5*pow(pow(sqrt(pow(Mass1,2)+pow(momentum1,2))+sqrt(pow(Mass2,2)+pow(momentum2,2))+sqrt(pow(Mass3,2)+pow(momentum3,2)),2)
								 -(pow(momentum1,2)+pow(momentum2,2)+pow(momentum3,2)+2*momentum1*momentum2*cosA+2*momentum2*momentum3*cosB+2*momentum1*momentum3*cosC),-0.5)
							 *(2*momentum1*pow(pow(Mass1,2)+pow(momentum1,2),-0.5)*(sqrt(pow(Mass1,2)+pow(momentum1,2))+sqrt(pow(Mass2,2)+pow(momentum2,2))
														+sqrt(pow(Mass3,2)+pow(momentum3,2)))-(2*momentum1+2*momentum2*cosA+2*momentum3*cosC)),2)*pow(errp1,2)+
						     pow(0.5*pow(pow(sqrt(pow(Mass1,2)+pow(momentum1,2))+sqrt(pow(Mass2,2)+pow(momentum2,2))+sqrt(pow(Mass3,2)+pow(momentum3,2)),2)
								 -(pow(momentum1,2)+pow(momentum2,2)+pow(momentum3,2)+2*momentum1*momentum2*cosA+2*momentum2*momentum3*cosB+2*momentum1*momentum3*cosC),-0.5)
							 *(2*momentum2*pow(pow(Mass2,2)+pow(momentum2,2),-0.5)*(sqrt(pow(Mass1,2)+pow(momentum1,2))+sqrt(pow(Mass2,2)+pow(momentum2,2))
														+sqrt(pow(Mass3,2)+pow(momentum3,2)))-(2*momentum2+2*momentum1*cosA+2*momentum3*cosB)),2)*pow(errp2,2)+
						     pow(0.5*pow(pow(sqrt(pow(Mass1,2)+pow(momentum1,2))+sqrt(pow(Mass2,2)+pow(momentum2,2))+sqrt(pow(Mass3,2)+pow(momentum3,2)),2)
								 -(pow(momentum1,2)+pow(momentum2,2)+pow(momentum3,2)+2*momentum1*momentum2*cosA+2*momentum2*momentum3*cosB+2*momentum1*momentum3*cosC),-0.5)
							 *(2*momentum3*pow(pow(Mass3,2)+pow(momentum3,2),-0.5)*(sqrt(pow(Mass1,2)+pow(momentum1,2))+sqrt(pow(Mass2,2)+pow(momentum2,2))
														+sqrt(pow(Mass3,2)+pow(momentum3,2)))-(2*momentum3+2*momentum2*cosB+2*momentum1*cosC)),2)*pow(errp3,2)+
						     pow(0.5*pow(pow(sqrt(pow(Mass1,2)+pow(momentum1,2))+sqrt(pow(Mass2,2)+pow(momentum2,2))+sqrt(pow(Mass3,2)+pow(momentum3,2)),2)
								 -(pow(momentum1,2)+pow(momentum2,2)+pow(momentum3,2)+2*momentum1*momentum2*cosA+2*momentum2*momentum3*cosB+2*momentum1*momentum3*cosC),-0.5)
							 *(2*momentum1*momentum2),2)*pow(errcosA,2)+
						     pow(0.5*pow(pow(sqrt(pow(Mass1,2)+pow(momentum1,2))+sqrt(pow(Mass2,2)+pow(momentum2,2))+sqrt(pow(Mass3,2)+pow(momentum3,2)),2)
								 -(pow(momentum1,2)+pow(momentum2,2)+pow(momentum3,2)+2*momentum1*momentum2*cosA+2*momentum2*momentum3*cosB+2*momentum1*momentum3*cosC),-0.5)
							 *(2*momentum2*momentum3),2)*pow(errcosB,2)+
						     pow(0.5*pow(pow(sqrt(pow(Mass1,2)+pow(momentum1,2))+sqrt(pow(Mass2,2)+pow(momentum2,2))+sqrt(pow(Mass3,2)+pow(momentum3,2)),2)
								 -(pow(momentum1,2)+pow(momentum2,2)+pow(momentum3,2)+2*momentum1*momentum2*cosA+2*momentum2*momentum3*cosB+2*momentum1*momentum3*cosC),-0.5)
							 *(2*momentum1*momentum3),2)*pow(errcosC,2));



		  //娘粒子の運動量の合計
		  double total_momentum = sqrt(pow(p1x+p2x+p3x,2)+pow(p1y+p2y+p3y,2)+pow(p1z+p2z+p3z,2));


		  //娘粒子の運動量の合計の誤差
		  double err_total_momentum = sqrt(pow((p1x+p2x+p3x)/sqrt(pow(p1x+p2x+p3x,2)+pow(p1y+p2y+p3y,2)+pow(p1z+p2z+p3z,2)),2)*pow(errp1x,2)+
						   pow((p1x+p2x+p3x)/sqrt(pow(p1x+p2x+p3x,2)+pow(p1y+p2y+p3y,2)+pow(p1z+p2z+p3z,2)),2)*pow(errp2x,2)+
						   pow((p1x+p2x+p3x)/sqrt(pow(p1x+p2x+p3x,2)+pow(p1y+p2y+p3y,2)+pow(p1z+p2z+p3z,2)),2)*pow(errp3x,2)+
						   pow((p1y+p2y+p3y)/sqrt(pow(p1x+p2x+p3x,2)+pow(p1y+p2y+p3y,2)+pow(p1z+p2z+p3z,2)),2)*pow(errp1y,2)+
						   pow((p1y+p2y+p3y)/sqrt(pow(p1x+p2x+p3x,2)+pow(p1y+p2y+p3y,2)+pow(p1z+p2z+p3z,2)),2)*pow(errp2y,2)+
						   pow((p1y+p2y+p3y)/sqrt(pow(p1x+p2x+p3x,2)+pow(p1y+p2y+p3y,2)+pow(p1z+p2z+p3z,2)),2)*pow(errp3y,2)+
						   pow((p1z+p2z+p3z)/sqrt(pow(p1x+p2x+p3x,2)+pow(p1y+p2y+p3y,2)+pow(p1z+p2z+p3z,2)),2)*pow(errp1z,2)+
						   pow((p1z+p2z+p3z)/sqrt(pow(p1x+p2x+p3x,2)+pow(p1y+p2y+p3y,2)+pow(p1z+p2z+p3z,2)),2)*pow(errp2z,2)+
						   pow((p1z+p2z+p3z)/sqrt(pow(p1x+p2x+p3x,2)+pow(p1y+p2y+p3y,2)+pow(p1z+p2z+p3z,2)),2)*pow(errp3z,2));



		  //親粒子の質量
		  double Caluculated_Mass;

		  Caluculated_Mass = sqrt(pow(sqrt(pow(Mass1,2)+pow(momentum1,2))+sqrt(pow(Mass2,2)+pow(momentum2,2))+sqrt(pow(Mass3,2)+pow(momentum3,2)),2)
					  -(pow(momentum1,2)+pow(momentum2,2)+pow(momentum3,2)+2*momentum1*momentum2*cosA+2*momentum2*momentum3*cosB+2*momentum1*momentum3*cosC));


		  Mass_gap = fabs(Caluculated_Mass - Estimated_Mass);//推定される親粒子と計算された質量との差の絶対値



		  double total_kinetic_energy = kinetic_energy1 + kinetic_energy2 + kinetic_energy3;//娘粒子の運動エネルギーの合計
		  double err_total_kinetic_energy = sqrt(pow(err_kinetic_energy1,2) + pow(err_kinetic_energy2,2) + pow(err_kinetic_energy3,2));//娘粒子の運動エネルギーの合計の誤差

		  double QQ = fabs(Q-total_kinetic_energy);//Q値と娘粒子の運動エネルギーの合計の差




		  if((Mass_gap <= err_Caluculated_Mass*3)&&(total_momentum <= err_total_momentum*3)&&(QQ <= err_total_kinetic_energy*3))//条件の合う崩壊過程のものだけresult.txtに書き出す
		    {
		      FILE *file1;
		      file1 = fopen("C:/Users/A-cham/Documents/OK.txt","a");
		      fprintf(file1,"ーーーーーーーーーーーーーーーーこの崩壊過程が考えられます！！！－－－－－－－－－－－－－\n\n");
		      fprintf(file1,"          反応式　　　\n");
		      fprintf(file1,"          %s  →  %s  +  %s  +  %s  \n\n",sort[SS][AA-1][ZZ-1],sort[S1][A1-1][Z1-1],sort[S2][A2-1][Z2-1],sort[S3][A3-1][Z3-1]);
		      fprintf(file1,"p1 = %9.3f ± %9.3f      P1x = %9.3f ± %9.3f\n",momentum1,errp1,p1x,errp1x);
		      fprintf(file1,"                                 p1y = %9.3f ± %9.3f\n",p1y,errp1y);
		      fprintf(file1,"                                 p1z = %9.3f ± %9.3f\n\n",p1z,errp1z);
		      fprintf(file1,"p2 = %9.3f ± %9.3f      p2x = %9.3f ± %9.3f\n",momentum2,errp2,p2x,errp2x);
		      fprintf(file1,"                                 p2y = %9.3f ± %9.3f\n",p2y,errp2y);
		      fprintf(file1,"                                 p2z = %9.3f ± %9.3f\n\n",p2z,errp2z);
		      fprintf(file1,"p3 = %9.3f ± %9.3f      p3x = %9.3f ± %9.3f\n",momentum3,errp3,p3x,errp3x);
		      fprintf(file1,"                                 p3y = %9.3f ± %9.3f\n",p3y,errp3y);
		      fprintf(file1,"                                 p3z = %9.3f ± %9.3f\n\n",p3z,errp3z);
		      fprintf(file1,"KE = %9.3f ± %9.3f\n",kinetic_energy1,err_kinetic_energy1);
		      fprintf(file1,"KE2 = %9.3f ± %9.3f\n",kinetic_energy2,err_kinetic_energy2);
		      fprintf(file1,"KE3 = %9.3f ± %9.3f\n\n",kinetic_energy3,err_kinetic_energy3);
		      fprintf(file1,"     Calculated_Mass = %9.3f ± %9.3f\n",Caluculated_Mass,err_Caluculated_Mass);
		      fprintf(file1,"       Estimated_Mass = %9.3f\n",Estimated_Mass);
		      fprintf(file1,"             Mass_gap = %9.3f\n",Mass_gap);
		      fprintf(file1,"                 Rate = %9.3f\n",Mass_gap/err_Caluculated_Mass);
		      fprintf(file1,"       Total_momentum = %9.3f ± %9.3f\n",total_momentum,err_total_momentum);
		      fprintf(file1,"                    Q = %9.3f\n",Q);
		      fprintf(file1," Total_Kinetic_Energy = %9.3f ± %9.3f\n\n",total_kinetic_energy,err_total_kinetic_energy);
		      fprintf(file1,"ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー\n\n\n\n\n");
		      fclose(file1);
		    }
		  else
		    {
		      FILE *file1;
		      file1 = fopen("C:/Users/A-cham/Documents/result.txt","a");
		      fprintf(file1,"ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー\n\n");
		      fprintf(file1,"         反応式　　　\n");
		      fprintf(file1,"          %s  →  %s  +  %s  +  %s  \n\n",sort[SS][AA-1][ZZ-1],sort[S1][A1-1][Z1-1],sort[S2][A2-1][Z2-1],sort[S3][A3-1][Z3-1]);
		      fprintf(file1,"p1 = %9.3f ± %9.3f      P1x = %9.3f ± %9.3f\n",momentum1,errp1,p1x,errp1x);
		      fprintf(file1,"                                 p1y = %9.3f ± %9.3f\n",p1y,errp1y);
		      fprintf(file1,"                                 p1z = %9.3f ± %9.3f\n\n",p1z,errp1z);
		      fprintf(file1,"p2 = %9.3f ± %9.3f      p2x = %9.3f ± %9.3f\n",momentum2,errp2,p2x,errp2x);
		      fprintf(file1,"                                 p2y = %9.3f ± %9.3f\n",p2y,errp2y);
		      fprintf(file1,"                                 p2z = %9.3f ± %9.3f\n\n",p2z,errp2z);
		      fprintf(file1,"p3 = %9.3f ± %9.3f      p3x = %9.3f ± %9.3f\n",momentum3,errp3,p3x,errp3x);
		      fprintf(file1,"                                 p3y = %9.3f ± %9.3f\n",p3y,errp3y);
		      fprintf(file1,"                                 p3z = %9.3f ± %9.3f\n\n",p3z,errp3z);
		      fprintf(file1,"KE1 = %9.3f ± %9.3f\n",kinetic_energy1,err_kinetic_energy1);
		      fprintf(file1,"KE2 = %9.3f ± %9.3f\n",kinetic_energy2,err_kinetic_energy2);
		      fprintf(file1,"KE3 = %9.3f ± %9.3f\n\n",kinetic_energy3,err_kinetic_energy3);
		      fprintf(file1,"     Calculated_Mass = %9.3f ± %9.3f\n",Caluculated_Mass,err_Caluculated_Mass);
		      fprintf(file1,"       Estimated_Mass = %9.3f\n",Estimated_Mass);
		      fprintf(file1,"             Mass_gap = %9.3f\n",Mass_gap);
		      fprintf(file1,"                 Rate = %9.3f\n",Mass_gap/err_Caluculated_Mass);
		      fprintf(file1,"       Total_momentum = %9.3f ± %9.3f\n",total_momentum,err_total_momentum);
		      fprintf(file1,"                    Q = %9.3f\n",Q);
		      fprintf(file1," Total_Kinetic_Energy = %9.3f ± %9.3f\n\n",total_kinetic_energy,err_total_kinetic_energy);
		      fprintf(file1,"ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー\n\n\n\n\n");
		      fclose(file1);
		    }



		}
	      }
	    }
	  }


	}

      }



    }
  }

  return 0;

}
