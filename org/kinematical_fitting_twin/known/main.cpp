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
#include "Math/BrentMinimizer1D.h"
using namespace std; 


#ifdef _DEBUG
    //Debugƒ‚[ƒh‚Ìê‡
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_core246d.lib")
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_imgproc246d.lib")
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_highgui246d.lib")
 
#else
    //Releaseƒ‚[ƒh‚Ìê‡
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_core246.lib")
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_imgproc246.lib")
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_highgui246.lib") 
#endif



int main()
{


	double m1,m2;
	double B_xi,initial_B_xi;

	double momentum1,momentum2;
	double err_momentum1,err_momentum2;

	

	double theta1,theta2;
	double err_theta1,err_theta2;
	double theta1_rad,theta2_rad;
	double err_theta1_rad,err_theta2_rad;
	
	double phi1,phi2;
	double err_phi1,err_phi2;
	double phi1_rad,phi2_rad;
	double err_phi1_rad,err_phi2_rad;
	
	
	double err_theta1_phi1_rad,err_theta2_phi2_rad;


	B_xi = initial_B_xi = 3.000;
	
	double m_C = 13040.207;//11174.866;//14895.084;//
	double m_xi = 1321.71;
	double err_m_xi = 0.07;

	m1 = 9499.326;
	m2 = 4839.943;
	

	

	momentum1 = 363.004;
    momentum2 = 344.294;
	

	theta1 = 136.2;
	theta2 = 41.7;
	

	phi1 = 13.3;
	phi2 = 193.5;
	

	err_momentum1 = 7.859;
	err_momentum2 = 0.953;
	

	err_theta1 = 11.87;
	err_theta2 = 5.29;
	

	err_phi1 = 12.00;
	err_phi2 = 5.34;
	

	err_theta1_phi1_rad = 0;//-5.13066*pow(10,-6);
	err_theta2_phi2_rad = 0;//-9.41478*pow(10,-6);
	



	double pi = 6*asin(0.5);//‰~ü—¦‚Ì’è‹`


	//‹‚ğrad‚É•ÏŠ·
	theta1_rad = theta1*pi/180;
	theta2_rad = theta2*pi/180;
	

	err_theta1_rad = err_theta1*pi/180;
	err_theta2_rad = err_theta2*pi/180;
	

	phi1_rad = phi1*pi/180;
	phi2_rad = phi2*pi/180;


	err_phi1_rad = err_phi1*pi/180;
	err_phi2_rad = err_phi2*pi/180;
	

	double a0_11 = momentum1;      //“ü—Íƒf[ƒ^‚Ì‚ ‚éƒpƒ‰ƒ[ƒ^
	double a0_21 = theta1_rad;
	double a0_31 = phi1_rad;
	double a0_41 = momentum2;
	double a0_51 = theta2_rad;
	double a0_61 = phi2_rad;
	double a0_71 = m_xi;


	cv::Mat a0 = (cv::Mat_<double>(7,1) <<  a0_11,a0_21,a0_31,a0_41,a0_51,a0_61,a0_71);  //data from emulsion matrix
								

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
	
	double D14 = momentum2/sqrt(pow(m2,2)+pow(momentum2,2));
	double D24 = sin(theta2_rad)*cos(phi2_rad);
	double D34 = sin(theta2_rad)*sin(phi2_rad);
	double D44 = cos(theta2_rad);
	
	double D15 = 0;
	double D25 = momentum2*cos(theta2_rad)*cos(phi2_rad);
	double D35 = momentum2*cos(theta2_rad)*sin(phi2_rad);
	double D45 = -momentum2*sin(theta2_rad);
	
	double D16 = 0;
	double D26 = -momentum2*sin(theta2_rad)*sin(phi2_rad);
	double D36 = momentum2*sin(theta2_rad)*cos(phi2_rad);
	double D46 = 0;

	double D17 = -1;
	double D27 = 0;
	double D37 = 0;
	double D47 = 0;
	
	

	cv::Mat D = (cv::Mat_<double>(4,7) <<   D11,D12,D13,D14,D15,D16,D17,						//constraints derivation matrix
										    D21,D22,D23,D24,D25,D26,D27,
										    D31,D32,D33,D34,D35,D36,D37,
										    D41,D42,D43,D44,D45,D46,D47);


	//covariance  momentum-theta,momentum-phi,,particle-particle = 0

	double Va11,Va12,Va13,Va14,Va15,Va16,Va17;
	double Va21,Va22,Va23,Va24,Va25,Va26,Va27;
	double Va31,Va32,Va33,Va34,Va35,Va36,Va37;
	double Va41,Va42,Va43,Va44,Va45,Va46,Va47;
	double Va51,Va52,Va53,Va54,Va55,Va56,Va57;
	double Va61,Va62,Va63,Va64,Va65,Va66,Va67;
	double Va71,Va72,Va73,Va74,Va75,Va76,Va77;
	
	
	
	Va14 = Va15 = Va16 = Va17 = 0;    //Še—±qŠÔ‚Ì‹¤•ªU‚O
	Va24 = Va25 = Va26 = Va27 = 0;
	Va34 = Va35 = Va36 = Va37 = 0;
	
	Va41  = Va42 = Va43 = 0;
	Va51  = Va52 = Va53 = 0;
	Va61  = Va62 = Va63 = 0;
	Va71  = Va72 = Va73 = 0;
	
	Va47 = 0;
	Va57 = 0;
	Va67 = 0;

	Va74 = Va75 = Va76 = 0;   //Še—±qŠÔ‚Ì‹¤•ªU‚O


	Va21 = Va12 = Va31 = Va13 = 0;     //‚P‚Â‚Ì—±q‚Ì‰^“®—Ê‚ÆƒÆAƒÓŠÔ‚Ì‹¤•ªU‚O
	Va45 = Va54 = Va46 = Va64 = 0;


	Va11 = pow(err_momentum1,2);     //Šeƒpƒ‰ƒ[ƒ^‚Ì•ªUAƒÆAƒÓŠÔ‚Ì‹¤•ªU
	Va22 = pow(err_theta1_rad,2);
	Va33 = pow(err_phi1_rad,2);
	Va32 = Va23 = err_theta1_phi1_rad;

	Va44 = pow(err_momentum2,2);
	Va55 = pow(err_theta2_rad,2);
	Va66 = pow(err_phi2_rad,2); 
	Va56 = Va65 = err_theta2_phi2_rad;

	Va77 = pow(err_m_xi,2);

	cv::Mat Va = (cv::Mat_<double>(7,7) <<     Va11,Va12,Va13,Va14,Va15,Va16,Va17,		//covariance matrix
										       Va21,Va22,Va23,Va24,Va25,Va26,Va27,
										       Va31,Va32,Va33,Va34,Va35,Va36,Va37,
										       Va41,Va42,Va43,Va44,Va45,Va46,Va47,
											   Va51,Va52,Va53,Va54,Va55,Va56,Va57,
											   Va61,Va62,Va63,Va64,Va65,Va66,Va67,
											   Va71,Va72,Va73,Va74,Va75,Va76,Va77);





//	double z0_11 = m_D;      //unknwonƒpƒ‰ƒ[ƒ^
//	double z0_21 = momentum2;
//
//	cv::Mat z0 = (cv::Mat_<double>(2,1) <<  z0_11,z0_21);  //matrix of unknown data
								



	double E11 = 1;
	double E21 = 0;
	double E31 = 0;
	double E41 = 0;


	cv::Mat E = (cv::Mat_<double>(4,1) << E11,
										  E21,
										  E31,
										  E41);  //unknown parameter 





	
	double H1,H2,H3,H4;//constraints 
	H1 = sqrt(pow(m1,2)+pow(momentum1,2))+sqrt(pow(m2,2)+pow(momentum2,2))-(m_C+m_xi-B_xi);
	H2 = momentum1*sin(theta1_rad)*cos(phi1_rad)+momentum2*sin(theta2_rad)*cos(phi2_rad);
	H3 = momentum1*sin(theta1_rad)*sin(phi1_rad)+momentum2*sin(theta2_rad)*sin(phi2_rad);
	H4 = momentum1*cos(theta1_rad)+momentum2*cos(theta2_rad);


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
	double new_momentum2 = a.at<double>(3,0);
	double new_theta2_rad = a.at<double>(4,0);
	double new_phi2_rad = a.at<double>(5,0);
	double new_m_xi = a.at<double>(6,0);


	double new_err_momentum1 = sqrt(V.at<double>(0,0));
	double new_err_theta1_do = sqrt(V.at<double>(1,1)) * 180 / pi;
	double new_err_phi1_do = sqrt(V.at<double>(2,2)) * 180 / pi;
	double new_err_momentum2 = sqrt(V.at<double>(3,3));
	double new_err_theta2_do = sqrt(V.at<double>(4,4)) * 180 / pi;
	double new_err_phi2_do = sqrt(V.at<double>(5,5)) * 180 / pi;
	double new_err_m_xi = sqrt(V.at<double>(6,6));


	double new_B_xi = initial_B_xi + z.at<double>(0,0);
	

	double new_err_B_xi = sqrt(Vz.at<double>(0,0));
	



//	printf("%.3f\n",new_momentum2);
//	printf("%.3f\n\n\n",new_m_D);
//	getchar();

	double new_H1,new_H2,new_H3,new_H4;//constraints 
	new_H1 = sqrt(pow(m1,2)+pow(new_momentum1,2))+sqrt(pow(m2,2)+pow(new_momentum2,2))-(m_C+new_m_xi-new_B_xi);
	new_H2 = new_momentum1*sin(new_theta1_rad)*cos(new_phi1_rad)+new_momentum2*sin(new_theta2_rad)*cos(new_phi2_rad);
	new_H3 = new_momentum1*sin(new_theta1_rad)*sin(new_phi1_rad)+new_momentum2*sin(new_theta2_rad)*sin(new_phi2_rad);
	new_H4 = new_momentum1*cos(new_theta1_rad)+new_momentum2*cos(new_theta2_rad);



	
	double delta_H1 = d.at<double>(0,0)/sqrt(Vd_i.at<double>(0,0));
	double delta_H2 = d.at<double>(1,0)/sqrt(Vd_i.at<double>(1,1));
	double delta_H3 = d.at<double>(2,0)/sqrt(Vd_i.at<double>(2,2));
	double delta_H4 = d.at<double>(3,0)/sqrt(Vd_i.at<double>(3,3));


	
		
	printf("new_H1      %.4f   %.4f \n",new_H1,delta_H1);
	printf("new_H2      %.4f   %.4f \n",new_H2,delta_H2);
	printf("new_H3      %.4f   %.4f \n",new_H3,delta_H3);
	printf("new_H4      %.4f   %.4f\n\n\n",new_H4,delta_H4);

	printf("chi_square         %.4f\n",chi_square_value);
	printf("momentum1 = %.3f } %.3f\n",new_momentum1,new_err_momentum1);
	printf("momentum2 = %.3f } %.3f\n",new_momentum2,new_err_momentum2);
	printf("theta1    = %.3f } %.3f\n",new_theta1_rad*180/pi,new_err_theta1_do);
	printf("phi1      = %.3f } %.3f\n",new_phi1_rad*180/pi,new_err_phi1_do);
	printf("theta2    = %.3f } %.3f\n",new_theta2_rad*180/pi,new_err_theta2_do);
	printf("phi2      = %.3f } %.3f\n",new_phi2_rad*180/pi,new_err_phi2_do);
	printf("B_xi(new)  = %.3f } %.3f\n",new_B_xi,new_err_B_xi);

	getchar();





	return 0;
}





