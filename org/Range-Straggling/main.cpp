#include<stdio.h>
#include"math.h"
double function2(double KEM);
double function3(double KEM);

int main()
{
	double KE,KEM,LKEM,E,P,B,Rs,Mass,M,Mp;//,Range_straggling1,Range_straggling2,Range_straggling3;

	int Z,A,S;

	loop:

	double MM[4][19][8] ={{{938.272,0.,0.,0.,0.,0.,0.,0.},
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
						{0.,0.,0.,0.,0.,0.,0.,0.},
						{0.,0.,0.,0.,0.,0.,0.,0.}},
					   {{0.,0.,0.,0.,0.,0.,0.,0.},
						{0.,0.,0.,0.,0.,0.,0.,0.},
						{2991.166,0.,0.,0.,0.,0.,0.,0.},
						{3922.565,3921.685,0.,0.,0.,0.,0.,0.},
						{0.,4839.942,0.,0.,0.,0.,0.,0.},
						{0.,5779.348,5778.807,0.,0.,0.,0.,0.},
						{0.,0.,6711.623,6715.821,0.,0.,0.,0.},
						{0.,7654.073,7642.719,7643.029,0.,0.,0.,0.},
						{0.,0.,8578.552,8563.825,8579.714,0.,0.,0.},
						{0.,0.,0.,9499.323,9500.103,0.,0.,0.},
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
					   {{139.570,0.,0.,0.,0.,0.,0.,0.},
					    {134.977,0.,0.,0.,0.,0.,0.,0.},
						{493.677,0.,0.,0.,0.,0.,0.,0.},
						{497.648,0.,0.,0.,0.,0.,0.,0.},
						{1869.4,0.,0.,0.,0.,0.,0.,0.},
						{1115.683,0.,0.,0.,0.,0.,0.,0.},
						{1189.37,0.,0.,0.,0.,0.,0.,0.},
						{1192.449,0.,0.,0.,0.,0.,0.,0.},
						{1314.83,0.,0.,0.,0.,0.,0.,0.},
						{1321.31,0.,0.,0.,0.,0.,0.,0.},
						{1672.45,0.,0.,0.,0.,0.,0.,0.},
						{0.,0.,0.,0.,0.,0.,0.,0.},
						{0.,0.,0.,0.,0.,0.,0.,0.},
						{0.,0.,0.,0.,0.,0.,0.,0.},
						{0.,0.,0.,0.,0.,0.,0.,0.},
						{0.,0.,0.,0.,0.,0.,0.,0.},
						{0.,0.,0.,0.,0.,0.,0.,0.},
						{0.,0.,0.,0.,0.,0.,0.,0.},
						{0.,0.,0.,0.,0.,0.,0.,0.}}};

	




 
 printf("\n");
 printf(" Z = ");
 scanf_s("%d",&Z);
 printf(" A = ");
 scanf_s("%d",&A);
 printf(" S = ");
 scanf_s("%d",&S);
 printf(" Mass = ");
 scanf_s("%lf",&Mass);
 
 
 printf(" KE = ");
 scanf_s("%lf",&KE);

//	Mass = MM[S][A-1][Z-1];






	Mp = 938.272;
	//Mass = 12278.859;//10429.883;//8563.825;//3922.565;//3727.380;//6533.836;//139.570;//3727.380;//139.570;//938.272;//139.570;//3727.380;//938.272;//139.570;//139.570;//938.272;//3727.380;//938.272;//139.570;//938.272;//139.570;//139.570;//3727.380;//139.570;//938.272;
	//Z = 6;

	//KE = 5.573;










//for(KE=0.01;KE<=50;KE+=0.01)
//{
    
	KEM = KE/Mass;
    E = Mass + KE;
    P = sqrt(E*E-Mass*Mass);
    B = P/E;
    LKEM = log10(KEM);
	M = Mass/Mp;



 if(KEM>0 && KEM<0.0001)
	{
	    Rs = 479.210*pow(KEM,0.675899);
	}
 else
	{
		Rs = function2(KEM);
	}

 double Rw,F,RR,D,D0,Rp,FX,Cz,R1,R2,R;

 D = 3.815;
 D0 = 3.815;
 RR = 0.885;


  if(KEM<0.0001)
    {
        Rw = 6.2813+1.5342*LKEM-0.15997*LKEM*LKEM;
        Rw = Rw-0.025245*LKEM*LKEM*LKEM;
        Rw = Rs+pow(10.0,Rw);
    }
    else
    {
		Rw = function3(KEM);
    }


    F = (D)/(D0)+(RR*(D0-D)*Rs)/((RR*D0-1.0)*Rw);
    Rp = Rs/F;
    
   
   if(Z>1.0)
   {
       FX = 137.0*B/Z;
       
       if(FX<=0.5)
		{
		Cz = 0.168550736771407*pow(FX,1.90707106569386);
		}
	   else if(FX<=2.51)
	   {
	   Cz =  0.002624371
		    -0.081622520*FX
		    +0.643381535*FX*FX
		    -0.903648583*FX*FX*FX
			+0.697505012*FX*FX*FX*FX
			-0.302935572*FX*FX*FX*FX*FX
			+0.067662990*FX*FX*FX*FX*FX*FX
			-0.006004180*FX*FX*FX*FX*FX*FX*FX;
	   }
		else
		{
		Cz = 0.217598079611354;//0.218764318175282;//0.220443959563487;
		}
   }
   else
   {
       Cz = 0.0;
   }
    
    R1 = Rp/(Z*Z)/Mp*Mass;
    R2 = Mass/Mp*Cz*pow(Z,2.0/3.0);
    R = R1+R2;
	


double rate_Range_straggling;
 double X = log10(KE/M);

if(X<log10(2000/M))
{
rate_Range_straggling = 1.0164402329033*pow(10.0,0.30754491034543-0.110592840462518*X)/100;
}
else
{
rate_Range_straggling = pow(10.0,-0.461423281494685+0.124489776442783*X)/100;
}


double Rpp;
if(KEM<0.0001)
{
Rpp = 479.210*pow(((KE/M)/938.272),0.675899);
}
else
{
//Rpp = pow(10.0, 1.147272863
//			        +1.481654835*X
//			        +0.156395018*X*X
//			        -0.078039243*X*X*X
//			        +0.065765281*X*X*X*X
//			        -0.030414179*X*X*X*X*X
//			        +0.005482754*X*X*X*X*X*X
//			        -0.000336044*X*X*X*X*X*X*X);

Rpp = function2((KE/M)/938.272);
}

double Range_straggling = (sqrt(M) / (Z*Z) ) * Rpp * rate_Range_straggling;


 //double Range_straggling = (sqrt(M) / (Z*Z) ) * pow(10,-0.501981130936234110213068180662+1.49146841766292474741422639802*X);

 double rate = Range_straggling/R*100;



// -0.501981130936228
//1.49146841766292


 printf("Range_straggling = %.5f [um]\n\n",Range_straggling);
 getchar();





// Range_straggling1 = (sqrt(M)/(Z*Z))*pow(10.0,-0.531508790
//                                             +1.362619405*X
//                                             +0.104614482*X*X
//                                             -0.038676012*X*X*X
//                                             +0.068558525*X*X*X*X
//                                             -0.039271090*X*X*X*X*X
//                                             +0.007828660*X*X*X*X*X*X
//                                             -0.000503500*X*X*X*X*X*X*X);

 //Range_straggling2 = (sqrt(M)/(Z*Z))*pow(10, -0.63227200185628797256055926851
//										   +1.705307505164731952882456538*X
//										   -0.05345977187545166980105768927*X*X);
 
 
 
 //Range_straggling3 = (sqrt(M)/(Z*Z))*pow(10,-0.501981130936234110213068180662+1.49146841766292474741422639802*X);






//double rate1 = Range_straggling1/R*100;
//double rate2 = Range_straggling2/R*100;
//double rate3 = Range_straggling3/R*100;


// FILE *file1;
// fopen_s(&file1,"result.txt","a");
// //fprintf(file1,"%.1f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f\n",KE,R,Range_straggling1,Range_straggling2,Range_straggling3,rate1,rate2,rate3);
// //fprintf(file1,"%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  \n",KE,R,rate_Range_straggling,Rpp,Range_straggling,rate);
// fprintf(file1,"%.2f  %.10f\n",KE,rate);
// fclose(file1);



 goto loop;
end:


 //}

	return 0;
}



double function2(double KEM)
{
	double LKEM = log10(KEM);

	double Rs;
	double Mp = 938.272;

	if(LKEM<log10(0.20000/Mp))
{
Rs =         0.8463768109 * LKEM +         3.3577221391;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(0.20000/Mp)<=LKEM && LKEM<log10(0.40000/Mp))
{
Rs =         1.1352913663 * LKEM +         4.4184137885;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(0.40000/Mp)<=LKEM && LKEM<log10(0.60000/Mp))
{
Rs =         1.3245936319 * LKEM +         5.0564133011;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(0.60000/Mp)<=LKEM && LKEM<log10(0.80000/Mp))
{
Rs =         1.4180699098 * LKEM +         5.3549931254;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(0.80000/Mp)<=LKEM && LKEM<log10(1.00000/Mp))
{
Rs =         1.4553837130 * LKEM +         5.4695180967;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(1.00000/Mp)<=LKEM && LKEM<log10(1.20000/Mp))
{
Rs =         1.4885031975 * LKEM +         5.5679600930;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(1.20000/Mp)<=LKEM && LKEM<log10(1.40000/Mp))
{
Rs =         1.5140162855 * LKEM +         5.6417732199;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(1.40000/Mp)<=LKEM && LKEM<log10(1.60000/Mp))
{
Rs =         1.5334393343 * LKEM +         5.6966666544;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(1.60000/Mp)<=LKEM && LKEM<log10(1.80000/Mp))
{
Rs =         1.5429425311 * LKEM +         5.7229734873;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(1.80000/Mp)<=LKEM && LKEM<log10(2.00000/Mp))
{
Rs =         1.5592671817 * LKEM +         5.7673284811;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(2.00000/Mp)<=LKEM && LKEM<log10(2.50000/Mp))
{
Rs =         1.5871926911 * LKEM +         5.8419258600;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(2.50000/Mp)<=LKEM && LKEM<log10(3.00000/Mp))
{
Rs =         1.6066746768 * LKEM +         5.8920800647;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(3.00000/Mp)<=LKEM && LKEM<log10(3.50000/Mp))
{
Rs =         1.6085809913 * LKEM +         5.8968367149;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(3.50000/Mp)<=LKEM && LKEM<log10(4.00000/Mp))
{
Rs =         1.6419695971 * LKEM +         5.9779129546;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(4.00000/Mp)<=LKEM && LKEM<log10(4.50000/Mp))
{
Rs =         1.6471916876 * LKEM +         5.9902907127;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(4.50000/Mp)<=LKEM && LKEM<log10(5.00000/Mp))
{
Rs =         1.6391660282 * LKEM +         5.9716782755;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(5.00000/Mp)<=LKEM && LKEM<log10(5.50000/Mp))
{
Rs =         1.6573310151 * LKEM +         6.0129738075;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(5.50000/Mp)<=LKEM && LKEM<log10(6.00000/Mp))
{
Rs =         1.6546600399 * LKEM +         6.0070122814;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(6.00000/Mp)<=LKEM && LKEM<log10(6.50000/Mp))
{
Rs =         1.6505154049 * LKEM +         5.9979182166;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(6.50000/Mp)<=LKEM && LKEM<log10(7.00000/Mp))
{
Rs =         1.6933142061 * LKEM +         6.0903386069;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(7.00000/Mp)<=LKEM && LKEM<log10(7.50000/Mp))
{
Rs =         1.6788954138 * LKEM +         6.0596665091;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(7.50000/Mp)<=LKEM && LKEM<log10(8.00000/Mp))
{
Rs =         1.6848885408 * LKEM +         6.0722356995;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(8.00000/Mp)<=LKEM && LKEM<log10(8.50000/Mp))
{
Rs =         1.6873221594 * LKEM +         6.0772714375;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(8.50000/Mp)<=LKEM && LKEM<log10(9.00000/Mp))
{
Rs =         1.6980049481 * LKEM +         6.0990954115;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(9.00000/Mp)<=LKEM && LKEM<log10(9.50000/Mp))
{
Rs =         1.6936788523 * LKEM +         6.0903649771;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(9.50000/Mp)<=LKEM && LKEM<log10(10.00000/Mp))
{
Rs =         1.7048590250 * LKEM +         6.1126650072;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(10.00000/Mp)<=LKEM && LKEM<log10(11.00000/Mp))
{
Rs =         1.7088879928 * LKEM +         6.1206114562;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(11.00000/Mp)<=LKEM && LKEM<log10(12.00000/Mp))
{
Rs =         1.7234004255 * LKEM +         6.1486340360;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(12.00000/Mp)<=LKEM && LKEM<log10(13.00000/Mp))
{
Rs =         1.7112246690 * LKEM +         6.1255835329;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(13.00000/Mp)<=LKEM && LKEM<log10(14.00000/Mp))
{
Rs =         1.7212888429 * LKEM +         6.1442866468;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(14.00000/Mp)<=LKEM && LKEM<log10(15.00000/Mp))
{
Rs =         1.7296610495 * LKEM +         6.1595759765;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(15.00000/Mp)<=LKEM && LKEM<log10(16.00000/Mp))
{
Rs =         1.7255660184 * LKEM +         6.1522203280;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(16.00000/Mp)<=LKEM && LKEM<log10(17.00000/Mp))
{
Rs =         1.7353019506 * LKEM +         6.1694354888;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(17.00000/Mp)<=LKEM && LKEM<log10(18.00000/Mp))
{
Rs =         1.7331479287 * LKEM +         6.1656834415;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(18.00000/Mp)<=LKEM && LKEM<log10(19.00000/Mp))
{
Rs =         1.7323708098 * LKEM +         6.1643490846;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(19.00000/Mp)<=LKEM && LKEM<log10(20.00000/Mp))
{
Rs =         1.7326317288 * LKEM +         6.1647909705;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(20.00000/Mp)<=LKEM && LKEM<log10(22.50000/Mp))
{
Rs =         1.7488913069 * LKEM +         6.1919655832;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(22.50000/Mp)<=LKEM && LKEM<log10(25.00000/Mp))
{
Rs =         1.7456856695 * LKEM +         6.1867719819;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(25.00000/Mp)<=LKEM && LKEM<log10(27.50000/Mp))
{
Rs =         1.7530353820 * LKEM +         6.1983432866;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(27.50000/Mp)<=LKEM && LKEM<log10(30.00000/Mp))
{
Rs =         1.7518538573 * LKEM +         6.1965320138;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(30.00000/Mp)<=LKEM && LKEM<log10(32.50000/Mp))
{
Rs =         1.7560114139 * LKEM +         6.2027484237;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(32.50000/Mp)<=LKEM && LKEM<log10(35.00000/Mp))
{
Rs =         1.7583207257 * LKEM +         6.2061210474;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(35.00000/Mp)<=LKEM && LKEM<log10(37.50000/Mp))
{
Rs =         1.7591183559 * LKEM +         6.2072602714;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(37.50000/Mp)<=LKEM && LKEM<log10(40.00000/Mp))
{
Rs =         1.7611310988 * LKEM +         6.2100746846;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(40.00000/Mp)<=LKEM && LKEM<log10(42.50000/Mp))
{
Rs =         1.7615988832 * LKEM +         6.2107156750;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(42.50000/Mp)<=LKEM && LKEM<log10(45.00000/Mp))
{
Rs =         1.7630519452 * LKEM +         6.2126685029;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(45.00000/Mp)<=LKEM && LKEM<log10(50.00000/Mp))
{
Rs =         1.7552508154 * LKEM +         6.2023779059;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(50.00000/Mp)<=LKEM && LKEM<log10(55.00000/Mp))
{
Rs =         1.7610046077 * LKEM +         6.2097045477;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(55.00000/Mp)<=LKEM && LKEM<log10(60.00000/Mp))
{
Rs =         1.7551403040 * LKEM +         6.2024799245;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(60.00000/Mp)<=LKEM && LKEM<log10(65.00000/Mp))
{
Rs =         1.7571380298 * LKEM +         6.2048655638;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(65.00000/Mp)<=LKEM && LKEM<log10(70.00000/Mp))
{
Rs =         1.7524528048 * LKEM +         6.1994334417;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(70.00000/Mp)<=LKEM && LKEM<log10(75.00000/Mp))
{
Rs =         1.7580457486 * LKEM +         6.2057379797;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(75.00000/Mp)<=LKEM && LKEM<log10(80.00000/Mp))
{
Rs =         1.7492530118 * LKEM +         6.1960899955;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(80.00000/Mp)<=LKEM && LKEM<log10(85.00000/Mp))
{
Rs =         1.7433695160 * LKEM +         6.1897991337;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(85.00000/Mp)<=LKEM && LKEM<log10(90.00000/Mp))
{
Rs =         1.7463804290 * LKEM +         6.1929392444;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(90.00000/Mp)<=LKEM && LKEM<log10(100.00000/Mp))
{
Rs =         1.7426147338 * LKEM +         6.1891054420;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(100.00000/Mp)<=LKEM && LKEM<log10(110.00000/Mp))
{
Rs =         1.7356763738 * LKEM +         6.1823590750;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(110.00000/Mp)<=LKEM && LKEM<log10(120.00000/Mp))
{
Rs =         1.7332450424 * LKEM +         6.1800956610;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(120.00000/Mp)<=LKEM && LKEM<log10(130.00000/Mp))
{
Rs =         1.7311283907 * LKEM +         6.1782051788;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(130.00000/Mp)<=LKEM && LKEM<log10(140.00000/Mp))
{
Rs =         1.7241936365 * LKEM +         6.1722524870;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(140.00000/Mp)<=LKEM && LKEM<log10(150.00000/Mp))
{
Rs =         1.7181126882 * LKEM +         6.1672284031;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(150.00000/Mp)<=LKEM && LKEM<log10(160.00000/Mp))
{
Rs =         1.7125184994 * LKEM +         6.1627741002;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(160.00000/Mp)<=LKEM && LKEM<log10(170.00000/Mp))
{
Rs =         1.7050531046 * LKEM +         6.1570391184;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(170.00000/Mp)<=LKEM && LKEM<log10(180.00000/Mp))
{
Rs =         1.7000689585 * LKEM +         6.1533414809;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(180.00000/Mp)<=LKEM && LKEM<log10(190.00000/Mp))
{
Rs =         1.6950059414 * LKEM +         6.1497110128;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(190.00000/Mp)<=LKEM && LKEM<log10(200.00000/Mp))
{
Rs =         1.6822039782 * LKEM +         6.1408318892;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(200.00000/Mp)<=LKEM && LKEM<log10(220.00000/Mp))
{
Rs =         1.6857293490 * LKEM +         6.1431984663;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(220.00000/Mp)<=LKEM && LKEM<log10(240.00000/Mp))
{
Rs =         1.6731756467 * LKEM +         6.1352908130;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(240.00000/Mp)<=LKEM && LKEM<log10(260.00000/Mp))
{
Rs =         1.6503273813 * LKEM +         6.1217619548;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(260.00000/Mp)<=LKEM && LKEM<log10(280.00000/Mp))
{
Rs =         1.6494833686 * LKEM +         6.1212915398;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(280.00000/Mp)<=LKEM && LKEM<log10(300.00000/Mp))
{
Rs =         1.6362826607 * LKEM +         6.1143589145;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(300.00000/Mp)<=LKEM && LKEM<log10(320.00000/Mp))
{
Rs =         1.6271350251 * LKEM +         6.1098289367;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(320.00000/Mp)<=LKEM && LKEM<log10(340.00000/Mp))
{
Rs =         1.6143136585 * LKEM +         6.1038390663;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(340.00000/Mp)<=LKEM && LKEM<log10(360.00000/Mp))
{
Rs =         1.6048262097 * LKEM +         6.0996565260;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(360.00000/Mp)<=LKEM && LKEM<log10(380.00000/Mp))
{
Rs =         1.5665186610 * LKEM +         6.0837195800;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(380.00000/Mp)<=LKEM && LKEM<log10(400.00000/Mp))
{
Rs =         1.6200196566 * LKEM +         6.1047211368;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(400.00000/Mp)<=LKEM && LKEM<log10(420.00000/Mp))
{
Rs =         1.5724030714 * LKEM +         6.0870902027;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(420.00000/Mp)<=LKEM && LKEM<log10(440.00000/Mp))
{
Rs =         1.5663374900 * LKEM +         6.0849728327;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(440.00000/Mp)<=LKEM && LKEM<log10(460.00000/Mp))
{
Rs =         1.5617734162 * LKEM +         6.0834718180;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(460.00000/Mp)<=LKEM && LKEM<log10(480.00000/Mp))
{
Rs =         1.5473922253 * LKEM +         6.0790198194;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(480.00000/Mp)<=LKEM && LKEM<log10(500.00000/Mp))
{
Rs =         1.5406126320 * LKEM +         6.0770463645;
Rs = pow(10.0,Rs);
return Rs;
}




else if(log10(500.00000/Mp)<=LKEM && LKEM<log10(520.00000/Mp))
{
Rs =         1.5298257351 * LKEM +         6.0740976718;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(520.00000/Mp)<=LKEM && LKEM<log10(540.00000/Mp))
{
Rs =         1.5258439380 * LKEM +         6.0730770360;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(540.00000/Mp)<=LKEM && LKEM<log10(560.00000/Mp))
{
Rs =         1.5125308185 * LKEM +         6.0698827527;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(560.00000/Mp)<=LKEM && LKEM<log10(580.00000/Mp))
{
Rs =         1.5058951365 * LKEM +         6.0683954261;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(580.00000/Mp)<=LKEM && LKEM<log10(600.00000/Mp))
{
Rs =         1.5002655728 * LKEM +         6.0672194060;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(600.00000/Mp)<=LKEM && LKEM<log10(620.00000/Mp))
{
Rs =         1.4906428611 * LKEM +         6.0653508918;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(620.00000/Mp)<=LKEM && LKEM<log10(640.00000/Mp))
{
Rs =         1.4869492586 * LKEM +         6.0646862758;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(640.00000/Mp)<=LKEM && LKEM<log10(660.00000/Mp))
{
Rs =         1.4744430938 * LKEM +         6.0626083918;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(660.00000/Mp)<=LKEM && LKEM<log10(680.00000/Mp))
{
Rs =         1.4678356470 * LKEM +         6.0615988742;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(680.00000/Mp)<=LKEM && LKEM<log10(700.00000/Mp))
{
Rs =         1.4666031861 * LKEM +         6.0614265517;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(700.00000/Mp)<=LKEM && LKEM<log10(720.00000/Mp))
{
Rs =         1.4521155942 * LKEM +         6.0595832850;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(720.00000/Mp)<=LKEM && LKEM<log10(740.00000/Mp))
{
Rs =         1.4523128717 * LKEM +         6.0596059712;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(740.00000/Mp)<=LKEM && LKEM<log10(760.00000/Mp))
{
Rs =         1.4394396565 * LKEM +         6.0582787808;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(760.00000/Mp)<=LKEM && LKEM<log10(780.00000/Mp))
{
Rs =         1.4363691909 * LKEM +         6.0579977867;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(780.00000/Mp)<=LKEM && LKEM<log10(800.00000/Mp))
{
Rs =         1.4293632224 * LKEM +         6.0574356687;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(800.00000/Mp)<=LKEM && LKEM<log10(820.00000/Mp))
{
Rs =         1.4229701478 * LKEM +         6.0569930201;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(820.00000/Mp)<=LKEM && LKEM<log10(840.00000/Mp))
{
Rs =         1.4171317355 * LKEM +         6.0566513860;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(840.00000/Mp)<=LKEM && LKEM<log10(860.00000/Mp))
{
Rs =         1.4160172008 * LKEM +         6.0565978332;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(860.00000/Mp)<=LKEM && LKEM<log10(880.00000/Mp))
{
Rs =         1.4026000502 * LKEM +         6.0560902583;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(880.00000/Mp)<=LKEM && LKEM<log10(900.00000/Mp))
{
Rs =         1.3900213974 * LKEM +         6.0557399920;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(900.00000/Mp)<=LKEM && LKEM<log10(920.00000/Mp))
{
Rs =         1.4193240127 * LKEM +         6.0562699664;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(920.00000/Mp)<=LKEM && LKEM<log10(940.00000/Mp))
{
Rs =         1.3658666593 * LKEM +         6.0558133909;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(940.00000/Mp)<=LKEM && LKEM<log10(960.00000/Mp))
{
Rs =         1.3959002423 * LKEM +         6.0557893912;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(960.00000/Mp)<=LKEM && LKEM<log10(980.00000/Mp))
{
Rs =         1.3845957829 * LKEM +         6.0559017855;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(980.00000/Mp)<=LKEM && LKEM<log10(1000.00000/Mp))
{
Rs =         1.3739190279 * LKEM +         6.0561035475;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(1000.00000/Mp)<=LKEM && LKEM<log10(1200.00000/Mp))
{
Rs =         1.3444513218 * LKEM +         6.0569189556;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(1200.00000/Mp)<=LKEM && LKEM<log10(1400.00000/Mp))
{
Rs =         1.3062405101 * LKEM +         6.0610018760;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(1400.00000/Mp)<=LKEM && LKEM<log10(1600.00000/Mp))
{
Rs =         1.2676490394 * LKEM +         6.0677090458;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(1600.00000/Mp)<=LKEM && LKEM<log10(1800.00000/Mp))
{
Rs =         1.2415009122 * LKEM +         6.0737699523;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(1800.00000/Mp)<=LKEM && LKEM<log10(2000.00000/Mp))
{
Rs =         1.2137228593 * LKEM +         6.0816295787;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(2000.00000/Mp)<=LKEM && LKEM<log10(2200.00000/Mp))
{
Rs =         1.1925202398 * LKEM +         6.0885989060;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(2200.00000/Mp)<=LKEM && LKEM<log10(2400.00000/Mp))
{
Rs =         1.1758892561 * LKEM +         6.0947539320;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(2400.00000/Mp)<=LKEM && LKEM<log10(2600.00000/Mp))
{
Rs =         1.1564978913 * LKEM +         6.1026633301;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(2600.00000/Mp)<=LKEM && LKEM<log10(2800.00000/Mp))
{
Rs =         1.1402053137 * LKEM +         6.1098751514;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(2800.00000/Mp)<=LKEM && LKEM<log10(3000.00000/Mp))
{
Rs =         1.1262730542 * LKEM +         6.1164905961;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(3000.00000/Mp)<=LKEM && LKEM<log10(3200.00000/Mp))
{
Rs =         1.1171576427 * LKEM +         6.1210919875;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(3200.00000/Mp)<=LKEM && LKEM<log10(3400.00000/Mp))
{
Rs =         1.1033423617 * LKEM +         6.1284530624;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(3400.00000/Mp)<=LKEM && LKEM<log10(3600.00000/Mp))
{
Rs =         1.0939037498 * LKEM +         6.1337306637;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(3600.00000/Mp)<=LKEM && LKEM<log10(3800.00000/Mp))
{
Rs =         1.0854289308 * LKEM +         6.1386797355;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(3800.00000/Mp)<=LKEM && LKEM<log10(4000.00000/Mp))
{
Rs =         1.0748159452 * LKEM +         6.1451266450;
Rs = pow(10.0,Rs);
return Rs;
}



else if(log10(4000.00000/Mp)<=LKEM && LKEM<=log10(4200.00000/Mp))
{
Rs =         1.0679731269 * LKEM +         6.1494357814;
Rs = pow(10.0,Rs);
return Rs;
}







else
{
return Rs = 0;
}




}



double function3(double KEM)
{

double LKEM = log10(KEM);

double Rw;
	double Mp = 938.272;

	if(LKEM<log10(0.20000/Mp))
{
Rw =         1.0641303374 * LKEM +         4.2684682249;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(0.20000/Mp)<=LKEM && LKEM<log10(0.40000/Mp))
{
Rw =         1.3344190391 * LKEM +         5.2607788003;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(0.40000/Mp)<=LKEM && LKEM<log10(0.60000/Mp))
{
Rw =         1.4637938697 * LKEM +         5.6968067508;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(0.60000/Mp)<=LKEM && LKEM<log10(0.80000/Mp))
{
Rw =         1.5287356869 * LKEM +         5.9042424427;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(0.80000/Mp)<=LKEM && LKEM<log10(1.00000/Mp))
{
Rw =         1.6203511705 * LKEM +         6.1854322366;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(1.00000/Mp)<=LKEM && LKEM<log10(1.20000/Mp))
{
Rw =         1.6129297913 * LKEM +         6.1633734579;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(1.20000/Mp)<=LKEM && LKEM<log10(1.40000/Mp))
{
Rw =         1.6669415534 * LKEM +         6.3196374531;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(1.40000/Mp)<=LKEM && LKEM<log10(1.60000/Mp))
{
Rw =         1.6637124771 * LKEM +         6.3105114351;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(1.60000/Mp)<=LKEM && LKEM<log10(1.80000/Mp))
{
Rw =         1.6945926719 * LKEM +         6.3959942613;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(1.80000/Mp)<=LKEM && LKEM<log10(2.00000/Mp))
{
Rw =         1.6946093363 * LKEM +         6.3960395396;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(2.00000/Mp)<=LKEM && LKEM<log10(2.40000/Mp))
{
Rw =         1.7060814371 * LKEM +         6.4266849483;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(2.40000/Mp)<=LKEM && LKEM<log10(2.80000/Mp))
{
Rw =         1.7168039891 * LKEM +         6.4544790630;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(2.80000/Mp)<=LKEM && LKEM<log10(3.20000/Mp))
{
Rw =         1.7255261180 * LKEM +         6.4765039274;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(3.20000/Mp)<=LKEM && LKEM<log10(3.60000/Mp))
{
Rw =         1.7362754752 * LKEM +         6.5030245134;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(3.60000/Mp)<=LKEM && LKEM<log10(4.00000/Mp))
{
Rw =         1.7397868632 * LKEM +         6.5115081191;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(4.00000/Mp)<=LKEM && LKEM<log10(4.40000/Mp))
{
Rw =         1.7491051283 * LKEM +         6.5335949117;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(4.40000/Mp)<=LKEM && LKEM<log10(4.80000/Mp))
{
Rw =         1.7482588592 * LKEM +         6.5316240559;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(4.80000/Mp)<=LKEM && LKEM<log10(5.20000/Mp))
{
Rw =         1.7562757963 * LKEM +         6.5499915604;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(5.20000/Mp)<=LKEM && LKEM<log10(5.60000/Mp))
{
Rw =         1.7581487428 * LKEM +         6.5542175373;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(5.60000/Mp)<=LKEM && LKEM<log10(6.00000/Mp))
{
Rw =         1.7611973274 * LKEM +         6.5609980184;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(6.00000/Mp)<=LKEM && LKEM<log10(6.80000/Mp))
{
Rw =         1.7644854207 * LKEM +         6.5682126787;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(6.80000/Mp)<=LKEM && LKEM<log10(7.60000/Mp))
{
Rw =         1.7718308775 * LKEM +         6.5839306329;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(7.60000/Mp)<=LKEM && LKEM<log10(8.40000/Mp))
{
Rw =         1.7746224991 * LKEM +         6.5897693519;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(8.40000/Mp)<=LKEM && LKEM<log10(9.20000/Mp))
{
Rw =         1.7834295773 * LKEM +         6.6078066837;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(9.20000/Mp)<=LKEM && LKEM<log10(10.00000/Mp))
{
Rw =         1.7726008365 * LKEM +         6.5860567146;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(10.00000/Mp)<=LKEM && LKEM<log10(12.00000/Mp))
{
Rw =         1.7878536403 * LKEM +         6.6161402582;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(12.00000/Mp)<=LKEM && LKEM<log10(14.00000/Mp))
{
Rw =         1.7915745050 * LKEM +         6.6231844040;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(14.00000/Mp)<=LKEM && LKEM<log10(16.00000/Mp))
{
Rw =         1.7967384454 * LKEM +         6.6326147956;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(16.00000/Mp)<=LKEM && LKEM<log10(18.00000/Mp))
{
Rw =         1.7994331731 * LKEM +         6.6373796367;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(18.00000/Mp)<=LKEM && LKEM<log10(20.00000/Mp))
{
Rw =         1.8001247827 * LKEM +         6.6385671694;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(20.00000/Mp)<=LKEM && LKEM<log10(24.00000/Mp))
{
Rw =         1.8036351188 * LKEM +         6.6444339898;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(24.00000/Mp)<=LKEM && LKEM<log10(28.00000/Mp))
{
Rw =         1.8054841239 * LKEM +         6.6473778232;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(28.00000/Mp)<=LKEM && LKEM<log10(32.00000/Mp))
{
Rw =         1.7991286159 * LKEM +         6.6376845884;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(32.00000/Mp)<=LKEM && LKEM<log10(36.00000/Mp))
{
Rw =         1.8165636373 * LKEM +         6.6632648817;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(36.00000/Mp)<=LKEM && LKEM<log10(40.00000/Mp))
{
Rw =         1.8033018386 * LKEM +         6.6444858266;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(40.00000/Mp)<=LKEM && LKEM<log10(44.00000/Mp))
{
Rw =         1.8056552130 * LKEM +         6.6477105821;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(44.00000/Mp)<=LKEM && LKEM<log10(48.00000/Mp))
{
Rw =         1.8095619228 * LKEM +         6.6529021153;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(48.00000/Mp)<=LKEM && LKEM<log10(52.00000/Mp))
{
Rw =         1.7986494693 * LKEM +         6.6388131828;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(52.00000/Mp)<=LKEM && LKEM<log10(56.00000/Mp))
{
Rw =         1.8018492308 * LKEM +         6.6428331245;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(56.00000/Mp)<=LKEM && LKEM<log10(60.00000/Mp))
{
Rw =         1.8010518994 * LKEM +         6.6418570786;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(60.00000/Mp)<=LKEM && LKEM<log10(64.00000/Mp))
{
Rw =         1.8016457418 * LKEM +         6.6425662319;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(64.00000/Mp)<=LKEM && LKEM<log10(68.00000/Mp))
{
Rw =         1.7946975435 * LKEM +         6.6344635988;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(68.00000/Mp)<=LKEM && LKEM<log10(72.00000/Mp))
{
Rw =         1.7940890411 * LKEM +         6.6337700158;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(72.00000/Mp)<=LKEM && LKEM<log10(76.00000/Mp))
{
Rw =         1.7944622296 * LKEM +         6.6341861195;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(76.00000/Mp)<=LKEM && LKEM<log10(80.00000/Mp))
{
Rw =         1.7880496185 * LKEM +         6.6271866572;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(80.00000/Mp)<=LKEM && LKEM<log10(84.00000/Mp))
{
Rw =         1.7906154549 * LKEM +         6.6299301490;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(84.00000/Mp)<=LKEM && LKEM<log10(88.00000/Mp))
{
Rw =         1.7829704792 * LKEM +         6.6219178363;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(88.00000/Mp)<=LKEM && LKEM<log10(92.00000/Mp))
{
Rw =         1.7837917847 * LKEM +         6.6227620119;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(92.00000/Mp)<=LKEM && LKEM<log10(96.00000/Mp))
{
Rw =         1.7817871748 * LKEM +         6.6207402808;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(96.00000/Mp)<=LKEM && LKEM<log10(100.00000/Mp))
{
Rw =         1.7995675231 * LKEM +         6.6383438484;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(100.00000/Mp)<=LKEM && LKEM<log10(120.00000/Mp))
{
Rw =         1.7680784633 * LKEM +         6.6077261301;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(120.00000/Mp)<=LKEM && LKEM<log10(140.00000/Mp))
{
Rw =         1.7559734176 * LKEM +         6.5969145386;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(140.00000/Mp)<=LKEM && LKEM<log10(160.00000/Mp))
{
Rw =         1.7477293018 * LKEM +         6.5901032442;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(160.00000/Mp)<=LKEM && LKEM<log10(180.00000/Mp))
{
Rw =         1.7355408314 * LKEM +         6.5807399543;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(180.00000/Mp)<=LKEM && LKEM<log10(200.00000/Mp))
{
Rw =         1.7202384788 * LKEM +         6.5697673067;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(200.00000/Mp)<=LKEM && LKEM<log10(240.00000/Mp))
{
Rw =         1.7035326827 * LKEM +         6.5585527265;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(240.00000/Mp)<=LKEM && LKEM<log10(280.00000/Mp))
{
Rw =         1.6793087147 * LKEM +         6.5442092908;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(280.00000/Mp)<=LKEM && LKEM<log10(320.00000/Mp))
{
Rw =         1.6570650075 * LKEM +         6.5325275469;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(320.00000/Mp)<=LKEM && LKEM<log10(360.00000/Mp))
{
Rw =         1.6343839720 * LKEM +         6.5219314485;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(360.00000/Mp)<=LKEM && LKEM<log10(400.00000/Mp))
{
Rw =         1.6136589250 * LKEM +         6.5133092848;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(400.00000/Mp)<=LKEM && LKEM<log10(500.00000/Mp))
{
Rw =         1.5798762100 * LKEM +         6.5008006006;
Rw = pow(10.0,Rw);
return Rw;
}




else if(log10(500.00000/Mp)<=LKEM && LKEM<log10(600.00000/Mp))
{
Rw =         1.5367803894 * LKEM +         6.4890199809;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(600.00000/Mp)<=LKEM && LKEM<log10(700.00000/Mp))
{
Rw =         1.4983953991 * LKEM +         6.4815664792;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(700.00000/Mp)<=LKEM && LKEM<log10(800.00000/Mp))
{
Rw =         1.4612210659 * LKEM +         6.4768367622;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(800.00000/Mp)<=LKEM && LKEM<log10(900.00000/Mp))
{
Rw =         1.4309603543 * LKEM +         6.4747415478;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(900.00000/Mp)<=LKEM && LKEM<log10(1000.00000/Mp))
{
Rw =         1.4035628718 * LKEM +         6.4742460301;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(1000.00000/Mp)<=LKEM && LKEM<log10(1200.00000/Mp))
{
Rw =         1.3690552944 * LKEM +         6.4752008977;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(1200.00000/Mp)<=LKEM && LKEM<log10(1400.00000/Mp))
{
Rw =         1.3289360139 * LKEM +         6.4794877427;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(1400.00000/Mp)<=LKEM && LKEM<log10(1600.00000/Mp))
{
Rw =         1.2947919415 * LKEM +         6.4854219579;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(1600.00000/Mp)<=LKEM && LKEM<log10(1800.00000/Mp))
{
Rw =         1.2667309749 * LKEM +         6.4919262438;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(1800.00000/Mp)<=LKEM && LKEM<log10(2000.00000/Mp))
{
Rw =         1.2402098497 * LKEM +         6.4994302304;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(2000.00000/Mp)<=LKEM && LKEM<log10(4000.00000/Mp))
{
Rw =         1.1606433452 * LKEM +         6.5255838390;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(4000.00000/Mp)<=LKEM && LKEM<log10(6000.00000/Mp))
{
Rw =         1.0696044442 * LKEM +         6.5829138786;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(6000.00000/Mp)<=LKEM && LKEM<log10(8000.00000/Mp))
{
Rw =         1.0295175541 * LKEM +         6.6152167963;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(8000.00000/Mp)<=LKEM && LKEM<log10(10000.00000/Mp))
{
Rw =         1.0045449299 * LKEM +         6.6384603467;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(10000.00000/Mp)<=LKEM && LKEM<log10(14000.00000/Mp))
{
Rw =         0.9848303768 * LKEM +         6.6587204261;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(14000.00000/Mp)<=LKEM && LKEM<log10(18000.00000/Mp))
{
Rw =         0.9681720937 * LKEM +         6.6782739067;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(18000.00000/Mp)<=LKEM && LKEM<log10(22000.00000/Mp))
{
Rw =         0.9591833749 * LKEM +         6.6898059274;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(22000.00000/Mp)<=LKEM && LKEM<log10(26000.00000/Mp))
{
Rw =         0.9531865910 * LKEM +         6.6980220845;
Rw = pow(10.0,Rw);
return Rw;
}



else if(log10(26000.00000/Mp)<=LKEM && LKEM<=log10(30000.00000/Mp))
{
Rw =         0.9493480684 * LKEM +         6.7035597085;
Rw = pow(10.0,Rw);
return Rw;
}





else
{
return Rw = 0;
}



}