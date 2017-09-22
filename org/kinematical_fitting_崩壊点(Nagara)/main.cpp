#pragma warning(disable : 4996)

#include<stdio.h>
#include"opencv2\\opencv.hpp"
#include<string>
#include<math.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>
#include "Math/GSLMinimizer.h"
#include "Minuit2/Minuit2Minimizer.h"
using namespace std; 


#ifdef _DEBUG
    //DebugÉÇÅ[ÉhÇÃèÍçá
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_core246d.lib")
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_imgproc246d.lib")
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_highgui246d.lib")
 
#else
    //ReleaseÉÇÅ[ÉhÇÃèÍçá
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_core246.lib")
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_imgproc246.lib")
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_highgui246.lib") 
#endif


double chi(const double *xx )
{


	double m1,m2,m3;//mass of single Hyper or nuclei
	double initial_m_D,m_D; //mass of double hyper nuclei 
	
	double momentum1,initial_momentum2,momentum2,momentum3;//momentum of single Hyper
	double err_momentum1,err_momentum2,err_momentum3;

	

	double theta1,theta2,theta3;
	double err_theta1,err_theta2,err_theta3;
	double theta1_rad,theta2_rad,theta3_rad;
	double err_theta1_rad,err_theta2_rad,err_theta3_rad;
	
	double phi1,phi2,phi3;
	double err_phi1,err_phi2,err_phi3;
	double phi1_rad,phi2_rad,phi3_rad;
	double err_phi1_rad,err_phi2_rad,err_phi3_rad;
	
	
	double err_theta1_phi1_rad,err_theta2_phi2_rad,err_theta3_phi3_rad;


	initial_m_D = m_D = xx[0];//5951.989;
	initial_momentum2 = momentum2 = xx[1];//86.179;



	m1 = 4839.943;
	m2 = 938.272;
	m3 = 139.570;

	

	momentum1 = 162.577;//162.532;//162.532;//170.961;
  //momentum2 = 82.661;
	momentum3 = 91.984;//91.994;//91.994;//166.049;

	theta1 = 78.651;//77.7;
	theta2 = 123.014;//122.8;
	theta3 = 81.010;//81.0;

	phi1 = 116.118;//115.9;
	phi2 = 284.138;//284.2;
	phi3 = 305.499;//305.5;

	err_momentum1 = 2.493;//3.506;
	//err_momentum2 = 5.855;
	err_momentum3 = 0.872;//1.333;

	err_theta1 = 4.48;//1.6;
	err_theta2 = 2.78;//1.0;
	err_theta3 = 0.69;//0.5;

	err_phi1 = 4.26;//0.8;
	err_phi2 = 2.69;//0.7;
	err_phi3 = 0.51;//0.2;

	err_theta1_phi1_rad = -5.13066*pow(10,-6);//0;//1.74231*pow(10,-5);//0;
	err_theta2_phi2_rad = -9.41478*pow(10,-6);//0;//7.87901*pow(10,-4);//0;
	err_theta3_phi3_rad = -7.44537*pow(10,-7);//0;//5.28957*pow(10,-6);//0;



	double pi = 6*asin(0.5);//â~é¸ó¶ÇÃíËã`


	//ÅãÇradÇ…ïœä∑
	theta1_rad = theta1*pi/180;
	theta2_rad = theta2*pi/180;
	theta3_rad = theta3*pi/180;

	err_theta1_rad = err_theta1*pi/180;
	err_theta2_rad = err_theta2*pi/180;
	err_theta3_rad = err_theta3*pi/180;

	phi1_rad = phi1*pi/180;
	phi2_rad = phi2*pi/180;
	phi3_rad = phi3*pi/180;

	err_phi1_rad = err_phi1*pi/180;
	err_phi2_rad = err_phi2*pi/180;
	err_phi3_rad = err_phi3*pi/180;

	double a0_11 = momentum1;      //ì¸óÕÉfÅ[É^ÇÃÇ†ÇÈÉpÉâÉÅÅ[É^
	double a0_21 = theta1_rad;
	double a0_31 = phi1_rad;
	double a0_41 = theta2_rad;
	double a0_51 = phi2_rad;
	double a0_61 = momentum3;
	double a0_71 = theta3_rad;
	double a0_81 = phi3_rad;


	cv::Mat a0 = (cv::Mat_<double>(8,1) <<  a0_11,a0_21,a0_31,a0_41,a0_51,a0_61,a0_71,a0_81);  //data from emulsion matrix
								

	
	
	double D11 = momentum1/sqrt(pow(m1,2)+pow(momentum1,2));
	double D21 = sin(theta1_rad)*cos(phi1_rad);
	double D31 = sin(theta1_rad)*sin(phi1_rad);
	double D41 = cos(theta1_rad);
	
	double D12 = 0;
	double D22 = momentum1*cos(theta1_rad)*cos(phi1_rad);
	double D32 = momentum1*cos(theta1_rad)*sin(phi1_rad);
	double D42 = -momentum1*sin(theta1_rad);
	
	double D13 = 0;
	double D23 = -momentum1*sin(theta1_rad)*sin(phi1_rad);
	double D33 = momentum1*sin(theta1_rad)*cos(phi1_rad);
	double D43 = 0;
	
	double D14 = 0;
	double D24 = momentum2*cos(theta2_rad)*cos(phi2_rad);
	double D34 = momentum2*cos(theta2_rad)*sin(phi2_rad);
	double D44 = -momentum2*sin(theta2_rad);
	
	double D15 = 0;
	double D25 = -momentum2*sin(theta2_rad)*sin(phi2_rad);
	double D35 = momentum2*sin(theta2_rad)*cos(phi2_rad);
	double D45 = 0;

	double D16 = momentum3/sqrt(pow(m3,2)+pow(momentum3,2));
	double D26 = sin(theta3_rad)*cos(phi3_rad);
	double D36 = sin(theta3_rad)*sin(phi3_rad);
	double D46 = cos(theta3_rad);
	
	double D17 = 0;
	double D27 = momentum3*cos(theta3_rad)*cos(phi3_rad);
	double D37 = momentum3*cos(theta3_rad)*sin(phi3_rad);
	double D47 = -momentum3*sin(theta3_rad);
	
	double D18 = 0;
	double D28 = -momentum3*sin(theta3_rad)*sin(phi3_rad);
	double D38 = momentum3*sin(theta3_rad)*cos(phi3_rad);
	double D48 = 0;
	

	cv::Mat D = (cv::Mat_<double>(4,8) <<   D11,D12,D13,D14,D15,D16,D17,D18,						//constraints derivation matrix
										    D21,D22,D23,D24,D25,D26,D27,D28,
										    D31,D32,D33,D34,D35,D36,D37,D38,
										    D41,D42,D43,D44,D45,D46,D47,D48);


	//covariance  momentum-theta,momentum-phi,,particle-particle = 0

	double Va11,Va12,Va13,Va14,Va15,Va16,Va17,Va18;
	double Va21,Va22,Va23,Va24,Va25,Va26,Va27,Va28;
	double Va31,Va32,Va33,Va34,Va35,Va36,Va37,Va38;
	
	double Va41,Va42,Va43,Va44,Va45,Va46,Va47,Va48;
	double Va51,Va52,Va53,Va54,Va55,Va56,Va57,Va58;
	
	double Va61,Va62,Va63,Va64,Va65,Va66,Va67,Va68;
	double Va71,Va72,Va73,Va74,Va75,Va76,Va77,Va78;
	double Va81,Va82,Va83,Va84,Va85,Va86,Va87,Va88;
	
	
	Va14 = Va15 = Va16 = Va17 = Va18 = 0;    //äeó±éqä‘ÇÃã§ï™éUÅÅÇO
	Va24 = Va25 = Va26 = Va27 = Va28 = 0;
	Va34 = Va35 = Va36 = Va37 = Va38 = 0;
	
	Va41 = Va42 = Va43 = Va46 = Va47 = Va48 = 0;
	Va51 = Va52 = Va53 = Va56 = Va57 = Va58 = 0;

	Va61 = Va62 = Va63 = Va64 = Va65 = 0;
	Va71 = Va72 = Va73 = Va74 = Va75 = 0;
	Va81 = Va82 = Va83 = Va84 = Va85 = 0;


	Va21 = Va12 = Va31 = Va13 = 0;     //ÇPÇ¬ÇÃó±éqÇÃâ^ìÆó Ç∆É∆ÅAÉ”ä‘ÇÃã§ï™éUÅÅÇO
	Va76 = Va67 = Va86 = Va68 = 0;


	Va11 = pow(err_momentum1,2);     //äeÉpÉâÉÅÅ[É^ÇÃï™éUÅAÉ∆ÅAÉ”ä‘ÇÃã§ï™éU
	Va22 = pow(err_theta1_rad,2);
	Va33 = pow(err_phi1_rad,2);
	Va32 = Va23 = err_theta1_phi1_rad;

	Va44 = pow(err_theta2_rad,2);
	Va55 = pow(err_phi2_rad,2); 
	Va54 = Va45 = err_theta2_phi2_rad;

	Va66 = pow(err_momentum3,2);
	Va77 = pow(err_theta3_rad,2);
	Va88 = pow(err_phi3_rad,2);
	Va87 = Va78 = err_theta3_phi3_rad;

	cv::Mat Va = (cv::Mat_<double>(8,8) <<     Va11,Va12,Va13,Va14,Va15,Va16,Va17,Va18,		//covariance matrix
										       Va21,Va22,Va23,Va24,Va25,Va26,Va27,Va28,
										       Va31,Va32,Va33,Va34,Va35,Va36,Va37,Va38,
										       Va41,Va42,Va43,Va44,Va45,Va46,Va47,Va48,
											   Va51,Va52,Va53,Va54,Va55,Va56,Va57,Va58,
											   Va61,Va62,Va63,Va64,Va65,Va66,Va67,Va68,
											   Va71,Va72,Va73,Va74,Va75,Va76,Va77,Va78,
											   Va81,Va82,Va83,Va84,Va85,Va86,Va87,Va88);





//	double z0_11 = m_D;      //unknwonÉpÉâÉÅÅ[É^
//	double z0_21 = momentum2;
//
//	cv::Mat z0 = (cv::Mat_<double>(2,1) <<  z0_11,z0_21);  //matrix of unknown data
								



	double E11 = momentum2/sqrt(pow(m2,2)+pow(momentum2,2));
	double E21 = sin(theta2_rad)*cos(phi2_rad);
	double E31 = sin(theta2_rad)*sin(phi2_rad);
	double E41 = cos(theta2_rad);

	double E12 = -1;
	double E22 = 0;
	double E32 = 0;
	double E42 = 0;

	cv::Mat E = (cv::Mat_<double>(4,2) << E11,E12,
										  E21,E22,
										  E31,E32,
										  E41,E42);  //unknown parameter 





	
	double H1,H2,H3,H4;//constraints 
	H1 = sqrt(pow(m1,2)+pow(momentum1,2))+sqrt(pow(m2,2)+pow(momentum2,2))+sqrt(pow(m3,2)+pow(momentum3,2))-m_D;
	H2 = momentum1*sin(theta1_rad)*cos(phi1_rad)+momentum2*sin(theta2_rad)*cos(phi2_rad)+momentum3*sin(theta3_rad)*cos(phi3_rad);
	H3 = momentum1*sin(theta1_rad)*sin(phi1_rad)+momentum2*sin(theta2_rad)*sin(phi2_rad)+momentum3*sin(theta3_rad)*sin(phi3_rad);
	H4 = momentum1*cos(theta1_rad)+momentum2*cos(theta2_rad)+momentum3*cos(theta3_rad);


	cv::Mat d = (cv::Mat_<double>(4,1) <<  H1,H2,H3,H4);





	cv::Mat Dt = D.t();

	
	cv::Mat Vd0 = D*Va*Dt; 
	cv::Mat Vd = Vd0.inv();
	
//	cv::Mat test = Vd0*Vd;



	cv::Mat E_t =  E.t();
	cv::Mat VE0 = E_t*Vd*E;
	cv::Mat VE = VE0.inv();

	

//	std::cout << "" << E << std::endl << std::endl;
//	getchar();
//
//	std::cout << "" << E_t << std::endl << std::endl;
//	getchar();



	cv::Mat lambda0 = Vd*d;

	cv::Mat z = -VE*E_t*lambda0;

	cv::Mat lambda = lambda0 + Vd*E*z;


	cv::Mat lambda_t = lambda.t();

	cv::Mat Vd_i = Vd.inv();

	cv::Mat chi_square = lambda_t*Vd_i*lambda;

		

	cv::Mat a = a0 - Va*Dt*lambda;//caluclated value

	cv::Mat Vz = VE;
	
	cv::Mat V_lambda = Vd - Vd*E*VE*E_t*Vd;
	cv::Mat V =  Va - Va*Dt*V_lambda*D*Va;


//	std::cout << "" << z << std::endl << std::endl;
//	getchar();


	
	double chi_square_value = chi_square.at<double>(0,0);


	double new_momentum1 = a.at<double>(0,0);
	double new_theta1_rad = a.at<double>(1,0);
	double new_phi1_rad = a.at<double>(2,0);
	double new_theta2_rad = a.at<double>(3,0);
	double new_phi2_rad = a.at<double>(4,0);
	double new_momentum3 = a.at<double>(5,0);
	double new_theta3_rad = a.at<double>(6,0);
	double new_phi3_rad = a.at<double>(7,0);


	double new_err_momentum1 = sqrt(V.at<double>(0,0));
	double new_err_theta1_do = sqrt(V.at<double>(1,1)) * 180 / pi;
	double new_err_phi1_do = sqrt(V.at<double>(2,2)) * 180 / pi;
	double new_err_theta2_do = sqrt(V.at<double>(3,3)) * 180 / pi;
	double new_err_phi2_do = sqrt(V.at<double>(4,4)) * 180 / pi;
	double new_err_momentum3 = sqrt(V.at<double>(5,5));
	double new_err_theta3_do = sqrt(V.at<double>(6,6)) * 180 / pi;
	double new_err_phi3_do = sqrt(V.at<double>(7,7)) * 180 / pi;


	double new_momentum2 = initial_momentum2 + z.at<double>(0,0);
	double new_m_D = initial_m_D + z.at<double>(1,0);

	double new_err_momentum2 = sqrt(Vz.at<double>(0,0));
	double new_err_m_D = sqrt(Vz.at<double>(1,1));



//	printf("%.3f\n",new_momentum2);
//	printf("%.3f\n\n\n",new_m_D);
//	getchar();

	double new_H1,new_H2,new_H3,new_H4;//constraints 
	new_H1 = sqrt(pow(m1,2)+pow(new_momentum1,2))+sqrt(pow(m2,2)+pow(new_momentum2,2))+sqrt(pow(m3,2)+pow(new_momentum3,2))-new_m_D;
	new_H2 = new_momentum1*sin(new_theta1_rad)*cos(new_phi1_rad)+new_momentum2*sin(new_theta2_rad)*cos(new_phi2_rad)+new_momentum3*sin(new_theta3_rad)*cos(new_phi3_rad);
	new_H3 = new_momentum1*sin(new_theta1_rad)*sin(new_phi1_rad)+new_momentum2*sin(new_theta2_rad)*sin(new_phi2_rad)+new_momentum3*sin(new_theta3_rad)*sin(new_phi3_rad);
	new_H4 = new_momentum1*cos(new_theta1_rad)+new_momentum2*cos(new_theta2_rad)+new_momentum3*cos(new_theta3_rad);



	
	double delta_H1 = d.at<double>(0,0)/sqrt(Vd_i.at<double>(0,0));
	double delta_H2 = d.at<double>(1,0)/sqrt(Vd_i.at<double>(1,1));
	double delta_H3 = d.at<double>(2,0)/sqrt(Vd_i.at<double>(2,2));
	double delta_H4 = d.at<double>(3,0)/sqrt(Vd_i.at<double>(3,3));


	
		
//	printf("new_H1      %.4f   %.4f \n",new_H1,delta_H1);
//	printf("new_H2      %.4f   %.4f \n",new_H2,delta_H2);
//	printf("new_H3      %.4f   %.4f \n",new_H3,delta_H3);
//	printf("new_H4      %.4f   %.4f\n\n\n",new_H4,delta_H4);
//
//	printf("chi_square         %.4f\n",chi_square_value);
//	printf("momentum1 = %.3f Å} %.3f\n",new_momentum1,new_err_momentum1);
//	printf("momentum2 = %.3f Å} %.3f\n",new_momentum2,new_err_momentum2);
//	printf("momentum3 = %.3f Å} %.3f\n",new_momentum3,new_err_momentum3);
//	printf("theta1    = %.3f Å} %.3f\n",new_theta1_rad*180/pi,new_err_theta1_do);
//	printf("phi1      = %.3f Å} %.3f\n",new_phi1_rad*180/pi,new_err_phi1_do);
//	printf("theta2    = %.3f Å} %.3f\n",new_theta2_rad*180/pi,new_err_theta2_do);
//	printf("phi2      = %.3f Å} %.3f\n",new_phi2_rad*180/pi,new_err_phi2_do);
//	printf("theta3    = %.3f Å} %.3f\n",new_theta3_rad*180/pi,new_err_theta3_do);
//	printf("phi3      = %.3f Å} %.3f\n",new_phi3_rad*180/pi,new_err_phi3_do);
//	printf("m_D       = %.3f\n",m_D);
//	printf("BLL       = %.3f\n",1116.483*2+3727.380-m_D);
//	printf("m_D(new)  = %.3f Å} %.3f\n",new_m_D,new_err_m_D);
//	printf("BLL(new)  = %.3f",1115.683*2+3727.380-new_m_D);
//
//	getchar();





	return chi_square_value;//0;
}



double chi2(double mass,double momentum)
{


	double m1,m2,m3;//mass of single Hyper or nuclei
	double initial_m_D,m_D; //mass of double hyper nuclei 
	
	double momentum1,initial_momentum2,momentum2,momentum3;//momentum of single Hyper
	double err_momentum1,err_momentum2,err_momentum3;

	

	double theta1,theta2,theta3;
	double err_theta1,err_theta2,err_theta3;
	double theta1_rad,theta2_rad,theta3_rad;
	double err_theta1_rad,err_theta2_rad,err_theta3_rad;
	
	double phi1,phi2,phi3;
	double err_phi1,err_phi2,err_phi3;
	double phi1_rad,phi2_rad,phi3_rad;
	double err_phi1_rad,err_phi2_rad,err_phi3_rad;
	
	
	double err_theta1_phi1_rad,err_theta2_phi2_rad,err_theta3_phi3_rad;


	initial_m_D = m_D = mass;//5951.989;
	initial_momentum2 = momentum2 = momentum;//86.179;



	m1 = 4839.943;
	m2 = 938.272;
	m3 = 139.570;

	

	momentum1 = 162.577;//162.532;//162.532;//170.961;
  //momentum2 = 82.661;
	momentum3 = 91.984;//91.994;//91.994;//166.049;

	theta1 = 78.651;//77.7;
	theta2 = 123.014;//122.8;
	theta3 = 81.010;//81.0;

	phi1 = 116.118;//115.9;
	phi2 = 284.138;//284.2;
	phi3 = 305.499;//305.5;

	err_momentum1 = 2.493;//3.506;
	//err_momentum2 = 5.855;
	err_momentum3 = 0.872;//1.333;

	err_theta1 = 4.48;//1.6;
	err_theta2 = 2.78;//1.0;
	err_theta3 = 0.69;//0.5;

	err_phi1 = 4.26;//0.8;
	err_phi2 = 2.69;//0.7;
	err_phi3 = 0.51;//0.2;

	err_theta1_phi1_rad = -5.13066*pow(10,-6);//0;//1.74231*pow(10,-5);//0;
	err_theta2_phi2_rad = -9.41478*pow(10,-6);//0;//7.87901*pow(10,-4);//0;
	err_theta3_phi3_rad = -7.44537*pow(10,-7);//0;//5.28957*pow(10,-6);//0;



	double pi = 6*asin(0.5);//â~é¸ó¶ÇÃíËã`


	//ÅãÇradÇ…ïœä∑
	theta1_rad = theta1*pi/180;
	theta2_rad = theta2*pi/180;
	theta3_rad = theta3*pi/180;

	err_theta1_rad = err_theta1*pi/180;
	err_theta2_rad = err_theta2*pi/180;
	err_theta3_rad = err_theta3*pi/180;

	phi1_rad = phi1*pi/180;
	phi2_rad = phi2*pi/180;
	phi3_rad = phi3*pi/180;

	err_phi1_rad = err_phi1*pi/180;
	err_phi2_rad = err_phi2*pi/180;
	err_phi3_rad = err_phi3*pi/180;

	double a0_11 = momentum1;      //ì¸óÕÉfÅ[É^ÇÃÇ†ÇÈÉpÉâÉÅÅ[É^
	double a0_21 = theta1_rad;
	double a0_31 = phi1_rad;
	double a0_41 = theta2_rad;
	double a0_51 = phi2_rad;
	double a0_61 = momentum3;
	double a0_71 = theta3_rad;
	double a0_81 = phi3_rad;


	cv::Mat a0 = (cv::Mat_<double>(8,1) <<  a0_11,a0_21,a0_31,a0_41,a0_51,a0_61,a0_71,a0_81);  //data from emulsion matrix
								

	
	
	double D11 = momentum1/sqrt(pow(m1,2)+pow(momentum1,2));
	double D21 = sin(theta1_rad)*cos(phi1_rad);
	double D31 = sin(theta1_rad)*sin(phi1_rad);
	double D41 = cos(theta1_rad);
	
	double D12 = 0;
	double D22 = momentum1*cos(theta1_rad)*cos(phi1_rad);
	double D32 = momentum1*cos(theta1_rad)*sin(phi1_rad);
	double D42 = -momentum1*sin(theta1_rad);
	
	double D13 = 0;
	double D23 = -momentum1*sin(theta1_rad)*sin(phi1_rad);
	double D33 = momentum1*sin(theta1_rad)*cos(phi1_rad);
	double D43 = 0;
	
	double D14 = 0;
	double D24 = momentum2*cos(theta2_rad)*cos(phi2_rad);
	double D34 = momentum2*cos(theta2_rad)*sin(phi2_rad);
	double D44 = -momentum2*sin(theta2_rad);
	
	double D15 = 0;
	double D25 = -momentum2*sin(theta2_rad)*sin(phi2_rad);
	double D35 = momentum2*sin(theta2_rad)*cos(phi2_rad);
	double D45 = 0;

	double D16 = momentum3/sqrt(pow(m3,2)+pow(momentum3,2));
	double D26 = sin(theta3_rad)*cos(phi3_rad);
	double D36 = sin(theta3_rad)*sin(phi3_rad);
	double D46 = cos(theta3_rad);
	
	double D17 = 0;
	double D27 = momentum3*cos(theta3_rad)*cos(phi3_rad);
	double D37 = momentum3*cos(theta3_rad)*sin(phi3_rad);
	double D47 = -momentum3*sin(theta3_rad);
	
	double D18 = 0;
	double D28 = -momentum3*sin(theta3_rad)*sin(phi3_rad);
	double D38 = momentum3*sin(theta3_rad)*cos(phi3_rad);
	double D48 = 0;
	

	cv::Mat D = (cv::Mat_<double>(4,8) <<   D11,D12,D13,D14,D15,D16,D17,D18,						//constraints derivation matrix
										    D21,D22,D23,D24,D25,D26,D27,D28,
										    D31,D32,D33,D34,D35,D36,D37,D38,
										    D41,D42,D43,D44,D45,D46,D47,D48);


	//covariance  momentum-theta,momentum-phi,,particle-particle = 0

	double Va11,Va12,Va13,Va14,Va15,Va16,Va17,Va18;
	double Va21,Va22,Va23,Va24,Va25,Va26,Va27,Va28;
	double Va31,Va32,Va33,Va34,Va35,Va36,Va37,Va38;
	
	double Va41,Va42,Va43,Va44,Va45,Va46,Va47,Va48;
	double Va51,Va52,Va53,Va54,Va55,Va56,Va57,Va58;
	
	double Va61,Va62,Va63,Va64,Va65,Va66,Va67,Va68;
	double Va71,Va72,Va73,Va74,Va75,Va76,Va77,Va78;
	double Va81,Va82,Va83,Va84,Va85,Va86,Va87,Va88;
	
	
	Va14 = Va15 = Va16 = Va17 = Va18 = 0;    //äeó±éqä‘ÇÃã§ï™éUÅÅÇO
	Va24 = Va25 = Va26 = Va27 = Va28 = 0;
	Va34 = Va35 = Va36 = Va37 = Va38 = 0;
	
	Va41 = Va42 = Va43 = Va46 = Va47 = Va48 = 0;
	Va51 = Va52 = Va53 = Va56 = Va57 = Va58 = 0;

	Va61 = Va62 = Va63 = Va64 = Va65 = 0;
	Va71 = Va72 = Va73 = Va74 = Va75 = 0;
	Va81 = Va82 = Va83 = Va84 = Va85 = 0;


	Va21 = Va12 = Va31 = Va13 = 0;     //ÇPÇ¬ÇÃó±éqÇÃâ^ìÆó Ç∆É∆ÅAÉ”ä‘ÇÃã§ï™éUÅÅÇO
	Va76 = Va67 = Va86 = Va68 = 0;


	Va11 = pow(err_momentum1,2);     //äeÉpÉâÉÅÅ[É^ÇÃï™éUÅAÉ∆ÅAÉ”ä‘ÇÃã§ï™éU
	Va22 = pow(err_theta1_rad,2);
	Va33 = pow(err_phi1_rad,2);
	Va32 = Va23 = err_theta1_phi1_rad;

	Va44 = pow(err_theta2_rad,2);
	Va55 = pow(err_phi2_rad,2); 
	Va54 = Va45 = err_theta2_phi2_rad;

	Va66 = pow(err_momentum3,2);
	Va77 = pow(err_theta3_rad,2);
	Va88 = pow(err_phi3_rad,2);
	Va87 = Va78 = err_theta3_phi3_rad;

	cv::Mat Va = (cv::Mat_<double>(8,8) <<     Va11,Va12,Va13,Va14,Va15,Va16,Va17,Va18,		//covariance matrix
										       Va21,Va22,Va23,Va24,Va25,Va26,Va27,Va28,
										       Va31,Va32,Va33,Va34,Va35,Va36,Va37,Va38,
										       Va41,Va42,Va43,Va44,Va45,Va46,Va47,Va48,
											   Va51,Va52,Va53,Va54,Va55,Va56,Va57,Va58,
											   Va61,Va62,Va63,Va64,Va65,Va66,Va67,Va68,
											   Va71,Va72,Va73,Va74,Va75,Va76,Va77,Va78,
											   Va81,Va82,Va83,Va84,Va85,Va86,Va87,Va88);





//	double z0_11 = m_D;      //unknwonÉpÉâÉÅÅ[É^
//	double z0_21 = momentum2;
//
//	cv::Mat z0 = (cv::Mat_<double>(2,1) <<  z0_11,z0_21);  //matrix of unknown data
								



	double E11 = momentum2/sqrt(pow(m2,2)+pow(momentum2,2));
	double E21 = sin(theta2_rad)*cos(phi2_rad);
	double E31 = sin(theta2_rad)*sin(phi2_rad);
	double E41 = cos(theta2_rad);

	double E12 = -1;
	double E22 = 0;
	double E32 = 0;
	double E42 = 0;

	cv::Mat E = (cv::Mat_<double>(4,2) << E11,E12,
										  E21,E22,
										  E31,E32,
										  E41,E42);  //unknown parameter 





	
	double H1,H2,H3,H4;//constraints 
	H1 = sqrt(pow(m1,2)+pow(momentum1,2))+sqrt(pow(m2,2)+pow(momentum2,2))+sqrt(pow(m3,2)+pow(momentum3,2))-m_D;
	H2 = momentum1*sin(theta1_rad)*cos(phi1_rad)+momentum2*sin(theta2_rad)*cos(phi2_rad)+momentum3*sin(theta3_rad)*cos(phi3_rad);
	H3 = momentum1*sin(theta1_rad)*sin(phi1_rad)+momentum2*sin(theta2_rad)*sin(phi2_rad)+momentum3*sin(theta3_rad)*sin(phi3_rad);
	H4 = momentum1*cos(theta1_rad)+momentum2*cos(theta2_rad)+momentum3*cos(theta3_rad);


	cv::Mat d = (cv::Mat_<double>(4,1) <<  H1,H2,H3,H4);





	cv::Mat Dt = D.t();

	
	cv::Mat Vd0 = D*Va*Dt; 
	cv::Mat Vd = Vd0.inv();
	
//	cv::Mat test = Vd0*Vd;



	cv::Mat E_t =  E.t();
	cv::Mat VE0 = E_t*Vd*E;
	cv::Mat VE = VE0.inv();

	

//	std::cout << "" << E << std::endl << std::endl;
//	getchar();
//
//	std::cout << "" << E_t << std::endl << std::endl;
//	getchar();



	cv::Mat lambda0 = Vd*d;

	cv::Mat z = -VE*E_t*lambda0;

	cv::Mat lambda = lambda0 + Vd*E*z;


	cv::Mat lambda_t = lambda.t();

	cv::Mat Vd_i = Vd.inv();

	cv::Mat chi_square = lambda_t*Vd_i*lambda;

		

	cv::Mat a = a0 - Va*Dt*lambda;//caluclated value

	cv::Mat Vz = VE;
	
	cv::Mat V_lambda = Vd - Vd*E*VE*E_t*Vd;
	cv::Mat V =  Va - Va*Dt*V_lambda*D*Va;


//	std::cout << "" << z << std::endl << std::endl;
//	getchar();


	
	double chi_square_value = chi_square.at<double>(0,0);


	double new_momentum1 = a.at<double>(0,0);
	double new_theta1_rad = a.at<double>(1,0);
	double new_phi1_rad = a.at<double>(2,0);
	double new_theta2_rad = a.at<double>(3,0);
	double new_phi2_rad = a.at<double>(4,0);
	double new_momentum3 = a.at<double>(5,0);
	double new_theta3_rad = a.at<double>(6,0);
	double new_phi3_rad = a.at<double>(7,0);


	double new_err_momentum1 = sqrt(V.at<double>(0,0));
	double new_err_theta1_do = sqrt(V.at<double>(1,1)) * 180 / pi;
	double new_err_phi1_do = sqrt(V.at<double>(2,2)) * 180 / pi;
	double new_err_theta2_do = sqrt(V.at<double>(3,3)) * 180 / pi;
	double new_err_phi2_do = sqrt(V.at<double>(4,4)) * 180 / pi;
	double new_err_momentum3 = sqrt(V.at<double>(5,5));
	double new_err_theta3_do = sqrt(V.at<double>(6,6)) * 180 / pi;
	double new_err_phi3_do = sqrt(V.at<double>(7,7)) * 180 / pi;


	double new_momentum2 = initial_momentum2 + z.at<double>(0,0);
	double new_m_D = initial_m_D + z.at<double>(1,0);

	double new_err_momentum2 = sqrt(Vz.at<double>(0,0));
	double new_err_m_D = sqrt(Vz.at<double>(1,1));



//	printf("%.3f\n",new_momentum2);
//	printf("%.3f\n\n\n",new_m_D);
//	getchar();

	double new_H1,new_H2,new_H3,new_H4;//constraints 
	new_H1 = sqrt(pow(m1,2)+pow(new_momentum1,2))+sqrt(pow(m2,2)+pow(new_momentum2,2))+sqrt(pow(m3,2)+pow(new_momentum3,2))-new_m_D;
	new_H2 = new_momentum1*sin(new_theta1_rad)*cos(new_phi1_rad)+new_momentum2*sin(new_theta2_rad)*cos(new_phi2_rad)+new_momentum3*sin(new_theta3_rad)*cos(new_phi3_rad);
	new_H3 = new_momentum1*sin(new_theta1_rad)*sin(new_phi1_rad)+new_momentum2*sin(new_theta2_rad)*sin(new_phi2_rad)+new_momentum3*sin(new_theta3_rad)*sin(new_phi3_rad);
	new_H4 = new_momentum1*cos(new_theta1_rad)+new_momentum2*cos(new_theta2_rad)+new_momentum3*cos(new_theta3_rad);



	
	double delta_H1 = d.at<double>(0,0)/sqrt(Vd_i.at<double>(0,0));
	double delta_H2 = d.at<double>(1,0)/sqrt(Vd_i.at<double>(1,1));
	double delta_H3 = d.at<double>(2,0)/sqrt(Vd_i.at<double>(2,2));
	double delta_H4 = d.at<double>(3,0)/sqrt(Vd_i.at<double>(3,3));

	
	double new_err_BLL = sqrt(4*pow(0.006,2)+pow(new_err_m_D,2));
	
	
		
	printf("new_H1      %.4f   %.4f \n",new_H1,delta_H1);
	printf("new_H2      %.4f   %.4f \n",new_H2,delta_H2);
	printf("new_H3      %.4f   %.4f \n",new_H3,delta_H3);
	printf("new_H4      %.4f   %.4f\n\n\n",new_H4,delta_H4);

	printf("chi_square         %.4f\n",chi_square_value);
	printf("momentum1 = %.3f Å} %.3f\n",new_momentum1,new_err_momentum1);
	printf("momentum2 = %.3f Å} %.3f\n",new_momentum2,new_err_momentum2);
	printf("momentum3 = %.3f Å} %.3f\n",new_momentum3,new_err_momentum3);
	printf("theta1    = %.3f Å} %.3f\n",new_theta1_rad*180/pi,new_err_theta1_do);
	printf("theta2    = %.3f Å} %.3f\n",new_theta2_rad*180/pi,new_err_theta2_do);
	printf("theta3    = %.3f Å} %.3f\n",new_theta3_rad*180/pi,new_err_theta3_do);
	printf("phi1      = %.3f Å} %.3f\n",new_phi1_rad*180/pi,new_err_phi1_do);
	printf("phi2      = %.3f Å} %.3f\n",new_phi2_rad*180/pi,new_err_phi2_do);
	printf("phi3      = %.3f Å} %.3f\n",new_phi3_rad*180/pi,new_err_phi3_do);
	printf("m_D       = %.3f\n",m_D);
	printf("BLL       = %.3f\n",1116.483*2+3727.380-m_D);
	printf("m_D(new)  = %.4f Å} %.4f\n",new_m_D,new_err_m_D);
	printf("BLL(new)  = %.4f Å} %.4f\n",1115.683*2+3727.380-new_m_D,new_err_m_D);

	getchar();





	return chi_square_value;//0;
}





int NumericalMinimization()
{
	ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );

	min.SetMaxFunctionCalls(1000000);
	min.SetMaxIterations(100000);
	min.SetTolerance(0.001);

	ROOT::Math::Functor f(&chi,2);
	double step[2] = {0.001,0.001};
	double variable[2] = {5952.600,80.000};
	
	min.SetFunction(f);

	std::string x;
	std::string y;

	min.SetLimitedVariable(0,x,variable[0],step[0],5950.000,5955.000);
	min.SetLimitedVariable(1,y,variable[1],step[1],70.000,90.000);

//	min.SetVariable(0,"X",-1.0,0.001);
//	min.SetVariable(1,"y",1.00,0.001);

	min.Minimize();

	const double *xs = min.X();

	double new_m_D = xs[0];
	double new_momentum2 = xs[1];
	double result = chi(xs);

	double result2 = chi2(new_m_D,new_momentum2);

//	printf("      new_m_D = %.3f\n",new_m_D);
//	printf("new_momentum2 = %.3f\n",new_momentum2);
//	printf("       result = %.3f\n",result);

   return 0;
}




int main()
{
	NumericalMinimization();
	getchar();
	return 0;

}