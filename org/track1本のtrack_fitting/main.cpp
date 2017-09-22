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
    //Debugモードの場合
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_core246d.lib")
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_imgproc246d.lib")
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_highgui246d.lib")
 
#else
    //Releaseモードの場合
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_core246.lib")
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_imgproc246.lib")
    #pragma comment(lib,"C:\\OpenCV\\opencv_2.4.6\\build\\x86\\vc11\\lib\\opencv_highgui246.lib") 
#endif
struct coordinate{


	double x,y,z;
};



bool ReadMCData(char* filename, std::vector<coordinate>& ve){
	coordinate e;
	FILE *data_in;
	fopen_s(&data_in,filename,"r");
	
	double	x,y,z;//track_cordinates
	double S;//shrinkage factor

	S = 2.1750;//2.089;//1.993;//2.262;//2.146;//2.262;//2.126;//2.246;//2.262;//1;//2.311;//1.85;//1.84712/1.8375;//2.125;//1.84712/1.8375;//2;

	while(fscanf_s(data_in,"%lf %lf %lf",
				&x,&y,&z)!=EOF){
					e.x = x;			
					e.y = y;
					e.z = S*z;//*1.84712/1.8375;
					ve.push_back(e);
	}
	return true;
}

int main()
{
	std::vector<coordinate> ve;
	ReadMCData("10.txt", ve);
	
	coordinate sum;
	
		sum.x=0.0;
		sum.y=0.0;
		sum.z=0.0;
		
	double sum_x_y=0.0;
	double sum_x_z=0.0;
	double sum_x_x=0.0;
	double sum_y_y=0.0;
	double sum_z_z=0.0;
	
	double err_x = 0.00025;
	double err_y = 0.00025;
	double err_z = 0.00036;


	for(int i= 0;i<ve.size();i++)
	{
			sum.x += ve[i].x;
			sum.y += ve[i].y;
			sum.z += ve[i].z;
			sum_x_y += ve[i].x * ve[i].y;
			sum_x_z += ve[i].x * ve[i].z;
			sum_x_x += ve[i].x * ve[i].x;
			sum_y_y += ve[i].y * ve[i].y;
			sum_z_z += ve[i].z * ve[i].z;
	}

	int n = ve.size();//データの数

	double alfa_0 = (n*sum_x_y-sum.y*sum.x)/(n*sum_x_x-pow(sum.x,2));
	double beta_0 = (sum.y*sum_x_x-sum_x_y*sum.x)/(n*sum_x_x-pow(sum.x,2));
	double gamma_0 = (n*sum_x_z-sum.z*sum.x)/(n*sum_x_x-pow(sum.x,2));
	double delta_0 = (sum.z*sum_x_x-sum_x_z*sum.x)/(n*sum_x_x-pow(sum.x,2));

	cv::Mat z0 = (cv::Mat_<double>(4,1) << alfa_0,beta_0,gamma_0,delta_0);

	
	// n x n の単位行列
   cv::Mat n_n_tani = cv::Mat::eye(n, n, CV_64F);

   // 要素がすべて 0 の n x n 行列
   cv::Mat n_n_zero = cv::Mat::zeros(n, n, CV_64F);

   //Ｄ行列の左上部分行列
   cv::Mat up_left = -alfa_0*n_n_tani;
   
   //Ｄ行列の左下部分行列
   cv::Mat down_left = -gamma_0*n_n_tani;

   //known parameter のヤコビ行列
   cv::Mat D = cv::Mat(2*n,3*n,CV_64F);
   cv::Mat roi1,roi2,roi3,roi4,roi5,roi6;

  roi1 = D(cv::Rect(0,0,n,n));//Ｄ行列の左上
  up_left.copyTo(roi1);
  roi2 = D(cv::Rect(0,n,n,n));//Ｄ行列の左下
  down_left.copyTo(roi2);
  roi3 = D(cv::Rect(n,0,n,n));//Ｄ行列の真ん中上
  n_n_tani.copyTo(roi3);
  roi4 = D(cv::Rect(n,n,n,n));//Ｄ行列の真ん中下
  n_n_zero.copyTo(roi4);
  roi5 = D(cv::Rect(2*n,0,n,n));//Ｄ行列の右上
  n_n_zero.copyTo(roi5);
  roi6 = D(cv::Rect(2*n,n,n,n));//Ｄ行列の右下
  n_n_tani.copyTo(roi6);



  //unknown parameter のヤコビ行列
  cv::Mat E = cv::Mat(2*n,4,CV_64F);


  for(int j= 0;j<ve.size();j++)
  {
  E.at<double>(j,0) = -1*ve[j].x;
  E.at<double>(j,1) = -1;
  E.at<double>(j,2) = 0;
  E.at<double>(j,3) = 0;
  E.at<double>(j+n,0) = 0;
  E.at<double>(j+n,1) = 0;
  E.at<double>(j+n,2) = -1*ve[j].x;
  E.at<double>(j+n,3) = -1;
  }



//  //初期constraint
//  cv::Mat d0 = cv::Mat(2*n,1,CV_64F);
//  for(int k = 0;k<ve.size();k++)
//  {
//	  d0.at<double>(k,0) = ve[k].y-alfa_0-beta_0*ve[k].x;
//	  d0.at<double>(k+n,0) = ve[k].z-gamma_0-delta_0*ve[k].x; 
//  }

  //constraint
  cv::Mat d = cv::Mat(2*n,1,CV_64F);
  for(int k = 0;k<ve.size();k++)
  {
	  d.at<double>(k,0) = alfa_0*ve[k].x;
	  d.at<double>(k+n,0) = gamma_0*ve[k].x;
  }



  //parameter
  cv::Mat a0 = cv::Mat(3*n,1,CV_64F);
  for(int h = 0;h<ve.size();h++)
  {
	  a0.at<double>(h,0) = ve[h].x;
	  a0.at<double>(h+n,0) = ve[h].y;
	  a0.at<double>(h+2*n,0) = ve[h].z;
  }

  //covariance matrix
  cv::Mat Va = cv::Mat::zeros(3*n,3*n,CV_64F);
  for(int g = 0;g<ve.size();g++)
  {
	  Va.at<double>(g,g) = pow(err_x,2);
	  Va.at<double>(g+n,g+n) = pow(err_y,2);
	  Va.at<double>(g+2*n,g+2*n) = pow(err_z,2);
  }


//    cv::Mat d = d0 - D*a0 - E*z0;
    cv::Mat Dt = D.t();
    cv::Mat Vd0 = D*Va*Dt; 
    cv::Mat Vd = Vd0.inv();
    
    
    cv::Mat E_t =  E.t();
	cv::Mat VE0 = E_t*Vd*E;
	cv::Mat VE = VE0.inv();

	cv::Mat lambda0 = Vd*(D*a0+d);


	cv::Mat z = -VE*E_t*lambda0;

	cv::Mat lambda = lambda0 + Vd*E*z;


	cv::Mat lambda_t = lambda.t();

	cv::Mat Vd_i = Vd.inv();

	cv::Mat chi_square = lambda_t*Vd_i*lambda;

		

	cv::Mat a = a0 - Va*Dt*lambda;

	cv::Mat Vz = VE;
	
	cv::Mat V_lambda = Vd - Vd*E*VE*E_t*Vd;
	cv::Mat V =  Va - Va*Dt*V_lambda*D*Va;

	double pi = 6*asin(0.5);//円周率の定義

	double alfa = z.at<double>(0,0);
	double beta = z.at<double>(1,0);
	double gamma = z.at<double>(2,0);
	double delta = z.at<double>(3,0);

	double err_alfa = Vz.at<double>(0,0);
	double err_beta = Vz.at<double>(1,1);
	double err_alfa_beta = Vz.at<double>(0,1);
	double err_gamma = Vz.at<double>(2,2);
	double err_delta = Vz.at<double>(3,3);
	double err_gamma_delta = Vz.at<double>(2,3);
	double err_alfa_gamma = Vz.at<double>(2,0);


	printf("alfa = %.7f +/- %.7f\n",alfa,sqrt(err_alfa));
	printf("beta = %.7f +/- %.7f\n",beta,sqrt(err_beta));
	printf("gamma = %.7f +/- %.7f\n",gamma,sqrt(err_gamma));
	printf("delta = %.7f +/- %.7f\n",delta,sqrt(err_delta));
	printf("err_alfa_beta = %.7f\n",err_alfa_beta);
	printf("err_gamma_delta = %.7f\n",err_gamma_delta);




	
	double x1 = ve[0].x;
	double x2 = ve[ve.size()-1].x;
	double y1 = beta + alfa*x1;
	double y2 = beta + alfa*x2;
	double z1 = delta + gamma*x1;
	double z2 = delta + gamma*x2;

	double range_theta = sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));
	double range_phi = sqrt(pow(x2-x1,2)+pow(y2-y1,2));

	double sin_phi = (y2-y1)/range_phi;

	double theta = acos((z2-z1)/range_theta)*180/pi;

	double phi = acos((x2-x1)/range_phi)*180/pi;
	if(sin_phi<0)
	{
	phi = 360 - phi;
	}

	


	double err_y1 = sqrt(pow(x1,2)*err_alfa+err_beta+pow(alfa,2)*pow(err_x,2)+2*x1*err_alfa_beta);
	double err_y2 = sqrt(pow(x2,2)*err_alfa+err_beta+pow(alfa,2)*pow(err_x,2)+2*x2*err_alfa_beta);

	double err_z1 = sqrt(pow(x1,2)*err_gamma+err_delta+pow(gamma,2)*pow(err_x,2)+2*x1*err_gamma_delta);
	double err_z2 = sqrt(pow(x2,2)*err_gamma+err_delta+pow(gamma,2)*pow(err_x,2)+2*x2*err_gamma_delta);


	double err_R = sqrt(2*pow(err_x,2)*pow((x2-x1)/range_theta,2)
					+(pow(err_y1,2)+pow(err_y2,2))*pow((y2-y1)/range_theta,2)
					+(pow(err_z1,2)+pow(err_z2,2))*pow((z2-z1)/range_theta,2));
	
	double err_delta_z = sqrt(pow(err_z2,2)+pow(err_z1,2));

	



	//fitting parameter を使った角度の算出と場合分け

	double phi2 = atan(alfa)*180/pi;


	if(phi2<0)
	{
		if(sin((y2-y1)/range_phi)<0)
		{
		phi2 = phi2 + 360;
		}
		else
		{
		phi2 = phi2 + 180;
		}

	}

	
	double theta2 = atan(1/(gamma*cos(phi2*pi/180)))*180/pi;

//	if(theta2<0)
//	{
//	theta2 = 180 + theta2;
//	}





	double err_phi2 = (sqrt(err_alfa)/(1+pow(alfa,2)))*180/pi;

	//tanθ = 1/(gamma*cosφ) = A

	double err_A = sqrt(pow(alfa/(gamma*sqrt(1+pow(alfa,2))),2)*err_alfa
					+pow(sqrt(1+pow(alfa,2))/pow(gamma,2),2)*err_gamma
					-2*(alfa/(gamma*sqrt(1+pow(alfa,2))))*(sqrt(1+pow(alfa,2))/pow(gamma,2))*err_alfa_gamma);

	double err_theta2 = (1/(1+pow(tan(theta2*pi/180),2)))*err_A*180/pi;


	double A = sqrt(1+pow(alfa,2))/gamma;

	double err_theta_phi = (1/(1+pow(alfa,2)))*(1/(1+pow(A,2)))*(alfa/(gamma*sqrt(1+pow(alfa,2))))*err_alfa
							+(1/(1+pow(alfa,2)))*(1/(1+pow(A,2)))*(-sqrt(1+pow(alfa,2))/gamma)*err_alfa_gamma;






//角度算出のためのRange
	double sum_Range = 0;
	for(int i=0;i<ve.size()-1;i++)
	{
		sum_Range += sqrt(pow(ve[i+1].x-ve[i].x,2)
						 +pow(ve[i+1].y-ve[i].y,2)
						 +pow(ve[i+1].z-ve[i].z,2));
	}
	sum_Range = sum_Range * 1000;


	
	printf("θ = %.1f +/- %.1f\n\n",theta,err_theta2);
	printf("φ = %.1f +/- %.1f\n\n",phi,err_phi2); 

	printf("θ = %.1f +/- %.1f\n\n",theta2,err_theta2);
	printf("φ = %.1f +/- %.1f\n\n",phi2,err_phi2);

	printf("%.3f um\n\n\n",sum_Range);
	printf("%.9f\n\n\n",err_theta_phi);


//	std::cout << z << std::endl << std::endl;
//	std::cout << z0 << std::endl << std::endl;
//	std::cout << a << std::endl << std::endl;
//
//	std::cout << Vz << std::endl << std::endl;
//	std::cout << D << std::endl << std::endl;
//	std::cout << E << std::endl << std::endl;
//	std::cout << Va << std::endl << std::endl;
//	std::cout << d << std::endl << std::endl;

//	std::cout << d0 << std::endl << std::endl;

//	std::cout << chi_square << std::endl << std::endl;

getchar();
return 0;
}











