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
#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"

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

//測定点からフィットした最終結果
double aa1;
double bb1;
double aa2;
double bb2;
double aa3;
double bb3;
double xx0;
double yy0;
double zz0;
double result;

//各座標値につく誤差
double err_x = 0.25;
double err_y = 0.25;
double err_z = 0.36;

double wx = pow(1/err_x,2);
double wy = pow(1/err_y,2);
double wz = pow(1/err_z,2);


int monte_number = 10000;

int number = 0;

struct coordinate{
	double x,y,z;
};


struct results{
	double a1,b1,a2,b2,a3,b3,x0,y0,z0,chi_2;
};



std::vector<coordinate> track1;
std::vector<coordinate> track2;
std::vector<coordinate> track3;


std::vector<std::vector<coordinate>> track1_monte;
std::vector<std::vector<coordinate>> track2_monte;
std::vector<std::vector<coordinate>> track3_monte;

bool ReadMCData(char* filename, std::vector<coordinate>& track){
	track.clear();
	coordinate tr;
	FILE *data_in;
	fopen_s(&data_in,filename,"r");
	
	double	x,y,z;//track_cordinates
	double S;//shrinkage factor

	S = 2.21;//1.84712/1.8375;

	while(fscanf_s(data_in,"%lf %lf %lf",
				&x,&y,&z)!=EOF){
					tr.x = x;			
					tr.y = y;
					tr.z = z*S;
					track.push_back(tr);
	}
	fclose(data_in);
	
	return true;
}









double chi(const double *xx)
{
//	std::vector<coordinate1> track1;
//	ReadMCData1("track2.txt", track1);
//	std::vector<coordinate2> track2;
//	ReadMCData2("track6.txt", track2);
//	std::vector<coordinate3> track3;
//	ReadMCData3("track8.txt", track3);

	
	double a1 = xx[0];
	double b1 = xx[1];
	double a2 = xx[2];
	double b2 = xx[3];
	double a3 = xx[4];
	double b3 = xx[5];
	double x0 = xx[6];
	double y0 = xx[7];
	double z0 = xx[8];


	
		
	double chi_square1=0.0;
	double chi_square2=0.0;
	double chi_square3=0.0;
	
	


   //chi_square1 = pow(a1-6,2)+pow(b1-2,2)+pow(a2-3,2)+pow(b2-4,2)+pow(a3-5,2)+pow(b3-6,2)+pow(x0-7,2)+pow(y0-8,2)+pow(z0-9,2);

	

	for( int m1= 0;m1<track1.size();m1++)
	{
		double term1 = (wx*(track1[m1].x-x0)+wy*a1*(track1[m1].y-y0)+wz*b1*(track1[m1].z-z0))/(wx+wy*a1*a1+wz*b1*b1);

		chi_square1 += wx*pow((track1[m1].x-x0)-term1,2)+wy*pow(track1[m1].y-a1*term1-y0,2)+wz*pow(track1[m1].z-b1*term1-z0,2);



	//chi_square1 +=  wx*pow((track1[m1].x-x0)-(((wx*(track1[m1].x-x0)+wy*a1*(track1[m1].y-y0)+wz*b1*(track1[m1].z-z0))/(wx+wy*a1*a1+wz*b1*b1))),2)
	//		+ wy*pow(track1[m1].y-a1*(((wx*(track1[m1].x-x0)+wy*a1*(track1[m1].y-y0)+wz*b1*(track1[m1].z-z0))/(wx+wy*a1*a1+wz*b1*b1)))-y0,2)
	//		+ wz*pow(track1[m1].z-b1*(((wx*(track1[m1].x-x0)+wy*a1*(track1[m1].y-y0)+wz*b1*(track1[m1].z-z0))/(wx+wy*a1*a1+wz*b1*b1)))-z0,2);
	}

	for( int m2= 0;m2<track2.size();m2++)
	{

		double term2 = (wx*(track2[m2].x-x0)+wy*a2*(track2[m2].y-y0)+wz*b2*(track2[m2].z-z0))/(wx+wy*a2*a2+wz*b2*b2);

		chi_square2 += wx*pow((track2[m2].x-x0)-term2,2)+wy*pow(track2[m2].y-a2*term2-y0,2)+wz*pow(track2[m2].z-b2*term2-z0,2);




//	chi_square2 +=  wx*pow((track2[m2].x-x0)-(((wx*(track2[m2].x-x0)+wy*a2*(track1[m2].y-y0)+wz*b2*(track2[m2].z-z0))/(wx+wy*a2*a2+wz*b2*b2))),2)
//			+ wy*pow(track2[m2].y-a2*(((wx*(track2[m2].x-x0)+wy*a2*(track2[m2].y-y0)+wz*b2*(track2[m2].z-z0))/(wx+wy*a2*a2+wz*b2*b2)))-y0,2)
//			+ wz*pow(track2[m2].z-b2*(((wx*(track2[m2].x-x0)+wy*a2*(track2[m2].y-y0)+wz*b2*(track2[m2].z-z0))/(wx+wy*a2*a2+wz*b2*b2)))-z0,2);
	}

	for( int m3= 0;m3<track3.size();m3++)
	{

		double term3 = (wx*(track3[m3].x-x0)+wy*a3*(track3[m3].y-y0)+wz*b3*(track3[m3].z-z0))/(wx+wy*a3*a3+wz*b3*b3);

		chi_square3 += wx*pow((track3[m3].x-x0)-term3,2)+wy*pow(track3[m3].y-a3*term3-y0,2)+wz*pow(track3[m3].z-b3*term3-z0,2);





//	chi_square3 +=  wx*pow((track3[m3].x-x0)-(((wx*(track3[m3].x-x0)+wy*a3*(track3[m3].y-y0)+wz*b3*(track3[m3].z-z0))/(wx+wy*a3*a3+wz*b3*b3))),2)
//			+ wy*pow(track3[m3].y-a3*(((wx*(track3[m3].x-x0)+wy*a3*(track3[m3].y-y0)+wz*b3*(track3[m3].z-z0))/(wx+wy*a3*a3+wz*b3*b3)))-y0,2)
//			+ wz*pow(track3[m3].z-b3*(((wx*(track3[m3].x-x0)+wy*a3*(track3[m3].y-y0)+wz*b3*(track3[m3].z-z0))/(wx+wy*a3*a3+wz*b3*b3)))-z0,2);
	}


	//   ((wx*(track1[m].x-x0)+wy*a1*(track1[m].y-y0)+wz*b1*(track1[m].z-z0))/(wx+wy*a1*a1+wz*b1*b1))


	double chi_square = chi_square1 + chi_square2 + chi_square3;







	return chi_square;
}



double chi_monte(const double *xx)
{
//	std::vector<coordinate1> track1;
//	ReadMCData1("track2.txt", track1);
//	std::vector<coordinate2> track2;
//	ReadMCData2("track6.txt", track2);
//	std::vector<coordinate3> track3;
//	ReadMCData3("track8.txt", track3);

	
	double a1 = xx[0];
	double b1 = xx[1];
	double a2 = xx[2];
	double b2 = xx[3];
	double a3 = xx[4];
	double b3 = xx[5];
	double x0 = xx[6];
	double y0 = xx[7];
	double z0 = xx[8];


	
		
	double chi_square1=0.0;
	double chi_square2=0.0;
	double chi_square3=0.0;
	
	


   //chi_square1 = pow(a1-6,2)+pow(b1-2,2)+pow(a2-3,2)+pow(b2-4,2)+pow(a3-5,2)+pow(b3-6,2)+pow(x0-7,2)+pow(y0-8,2)+pow(z0-9,2);

	

	for( int m1= 0;m1<track1.size();m1++)
	{
		double term1 = (wx*(track1_monte[number][m1].x-x0)+wy*a1*(track1_monte[number][m1].y-y0)+wz*b1*(track1_monte[number][m1].z-z0))/(wx+wy*a1*a1+wz*b1*b1);

		chi_square1 += wx*pow((track1_monte[number][m1].x-x0)-term1,2)+wy*pow(track1_monte[number][m1].y-a1*term1-y0,2)+wz*pow(track1_monte[number][m1].z-b1*term1-z0,2);



	//chi_square1 +=  wx*pow((track1[m1].x-x0)-(((wx*(track1[m1].x-x0)+wy*a1*(track1[m1].y-y0)+wz*b1*(track1[m1].z-z0))/(wx+wy*a1*a1+wz*b1*b1))),2)
	//		+ wy*pow(track1[m1].y-a1*(((wx*(track1[m1].x-x0)+wy*a1*(track1[m1].y-y0)+wz*b1*(track1[m1].z-z0))/(wx+wy*a1*a1+wz*b1*b1)))-y0,2)
	//		+ wz*pow(track1[m1].z-b1*(((wx*(track1[m1].x-x0)+wy*a1*(track1[m1].y-y0)+wz*b1*(track1[m1].z-z0))/(wx+wy*a1*a1+wz*b1*b1)))-z0,2);
	}

	for( int m2= 0;m2<track2.size();m2++)
	{

		double term2 = (wx*(track2_monte[number][m2].x-x0)+wy*a2*(track2_monte[number][m2].y-y0)+wz*b2*(track2_monte[number][m2].z-z0))/(wx+wy*a2*a2+wz*b2*b2);

		chi_square2 += wx*pow((track2_monte[number][m2].x-x0)-term2,2)+wy*pow(track2_monte[number][m2].y-a2*term2-y0,2)+wz*pow(track2_monte[number][m2].z-b2*term2-z0,2);




//	chi_square2 +=  wx*pow((track2[m2].x-x0)-(((wx*(track2[m2].x-x0)+wy*a2*(track1[m2].y-y0)+wz*b2*(track2[m2].z-z0))/(wx+wy*a2*a2+wz*b2*b2))),2)
//			+ wy*pow(track2[m2].y-a2*(((wx*(track2[m2].x-x0)+wy*a2*(track2[m2].y-y0)+wz*b2*(track2[m2].z-z0))/(wx+wy*a2*a2+wz*b2*b2)))-y0,2)
//			+ wz*pow(track2[m2].z-b2*(((wx*(track2[m2].x-x0)+wy*a2*(track2[m2].y-y0)+wz*b2*(track2[m2].z-z0))/(wx+wy*a2*a2+wz*b2*b2)))-z0,2);
	}

	for( int m3= 0;m3<track3.size();m3++)
	{

		double term3 = (wx*(track3_monte[number][m3].x-x0)+wy*a3*(track3_monte[number][m3].y-y0)+wz*b3*(track3_monte[number][m3].z-z0))/(wx+wy*a3*a3+wz*b3*b3);

		chi_square3 += wx*pow((track3_monte[number][m3].x-x0)-term3,2)+wy*pow(track3_monte[number][m3].y-a3*term3-y0,2)+wz*pow(track3_monte[number][m3].z-b3*term3-z0,2);





//	chi_square3 +=  wx*pow((track3[m3].x-x0)-(((wx*(track3[m3].x-x0)+wy*a3*(track3[m3].y-y0)+wz*b3*(track3[m3].z-z0))/(wx+wy*a3*a3+wz*b3*b3))),2)
//			+ wy*pow(track3[m3].y-a3*(((wx*(track3[m3].x-x0)+wy*a3*(track3[m3].y-y0)+wz*b3*(track3[m3].z-z0))/(wx+wy*a3*a3+wz*b3*b3)))-y0,2)
//			+ wz*pow(track3[m3].z-b3*(((wx*(track3[m3].x-x0)+wy*a3*(track3[m3].y-y0)+wz*b3*(track3[m3].z-z0))/(wx+wy*a3*a3+wz*b3*b3)))-z0,2);
	}


	//   ((wx*(track1[m].x-x0)+wy*a1*(track1[m].y-y0)+wz*b1*(track1[m].z-z0))/(wx+wy*a1*a1+wz*b1*b1))


	double chi_square = chi_square1 + chi_square2 + chi_square3;







	return chi_square;
}


int NumericalMinimization()
{
	ReadMCData("track1.txt", track1);
	ReadMCData("track4.txt", track2);
	ReadMCData("track5.txt", track3);
	
	
	
	//track1//
	coordinate sum1;
	
		sum1.x=0.0;
		sum1.y=0.0;
		sum1.z=0.0;
		
	double sum1_x_y=0.0;
	double sum1_x_z=0.0;
	double sum1_x_x=0.0;
	double sum1_y_y=0.0;
	double sum1_z_z=0.0;

	for(unsigned int i= 0; i<track1.size(); i++)
	{
			sum1.x += track1[i].x;
			sum1.y += track1[i].y;
			sum1.z += track1[i].z;
			sum1_x_y += track1[i].x * track1[i].y;
			sum1_x_z += track1[i].x * track1[i].z;
			sum1_x_x += track1[i].x * track1[i].x;
			sum1_y_y += track1[i].y * track1[i].y;
			sum1_z_z += track1[i].z * track1[i].z;
	}

	int n1 = track1.size();//データの数
	//trackの傾き
	double a1_0 = (n1*sum1_x_y-sum1.y*sum1.x)/(n1*sum1_x_x-pow(sum1.x,2));
	double b1_0 = (n1*sum1_x_z-sum1.z*sum1.x)/(n1*sum1_x_x-pow(sum1.x,2));
	//trackの切片
	double c1_0 = (sum1.y*sum1_x_x-sum1_x_y*sum1.x)/(n1*sum1_x_x-pow(sum1.x,2));
	double d1_0 = (sum1.z*sum1_x_x-sum1_x_z*sum1.x)/(n1*sum1_x_x-pow(sum1.x,2));






	//track2//
	coordinate sum2;
	
		sum2.x=0.0;
		sum2.y=0.0;
		sum2.z=0.0;
		
	double sum2_x_y=0.0;
	double sum2_x_z=0.0;
	double sum2_x_x=0.0;
	double sum2_y_y=0.0;
	double sum2_z_z=0.0;

	for(unsigned int j= 0;j<track2.size();j++)
	{
			sum2.x += track2[j].x;
			sum2.y += track2[j].y;
			sum2.z += track2[j].z;
			sum2_x_y += track2[j].x * track2[j].y;
			sum2_x_z += track2[j].x * track2[j].z;
			sum2_x_x += track2[j].x * track2[j].x;
			sum2_y_y += track2[j].y * track2[j].y;
			sum2_z_z += track2[j].z * track2[j].z;
	}

	int n2 = track2.size();//データの数
	//trackの傾き
	double a2_0 = (n2*sum2_x_y-sum2.y*sum2.x)/(n2*sum2_x_x-pow(sum2.x,2));
	double b2_0 = (n2*sum2_x_z-sum2.z*sum2.x)/(n2*sum2_x_x-pow(sum2.x,2));
	//trackの切片
	double c2_0 = (sum2.y*sum2_x_x-sum2_x_y*sum2.x)/(n2*sum2_x_x-pow(sum2.x,2));
	double d2_0 = (sum2.z*sum2_x_x-sum2_x_z*sum2.x)/(n2*sum2_x_x-pow(sum2.x,2));


	//track3//
	coordinate sum3;
	
		sum3.x=0.0;
		sum3.y=0.0;
		sum3.z=0.0;
		
	double sum3_x_y=0.0;
	double sum3_x_z=0.0;
	double sum3_x_x=0.0;
	double sum3_y_y=0.0;
	double sum3_z_z=0.0;

	for(unsigned int k= 0;k<track3.size();k++)
	{
			sum3.x += track3[k].x;
			sum3.y += track3[k].y;
			sum3.z += track3[k].z;
			sum3_x_y += track3[k].x * track3[k].y;
			sum3_x_z += track3[k].x * track3[k].z;
			sum3_x_x += track3[k].x * track3[k].x;
			sum3_y_y += track3[k].y * track3[k].y;
			sum3_z_z += track3[k].z * track3[k].z;
	}

	int n3 = track3.size();//データの数
	//trackの傾き
	double a3_0 = (n3*sum3_x_y-sum3.y*sum3.x)/(n3*sum3_x_x-pow(sum3.x,2));
	double b3_0 = (n3*sum3_x_z-sum3.z*sum3.x)/(n3*sum3_x_x-pow(sum3.x,2));
	//trackの切片
	double c3_0 = (sum3.y*sum3_x_x-sum3_x_y*sum3.x)/(n3*sum3_x_x-pow(sum3.x,2));
	double d3_0 = (sum3.z*sum3_x_x-sum3_x_z*sum3.x)/(n3*sum3_x_x-pow(sum3.x,2));



	//vertex point の初期値
	
	cv::Mat A1 = (cv::Mat_<double>(2,2) << pow(a1_0,2)+pow(b1_0,2)+1,-(a1_0*a2_0+b1_0*b2_0+1),
										(a1_0*a2_0+b1_0*b2_0+1),-(pow(a2_0,2)+pow(b2_0,2)+1));
	cv::Mat A2 = (cv::Mat_<double>(2,1) << a1_0*(c2_0-c1_0)+b1_0*(d2_0-d1_0),
										  a2_0*(c2_0-c1_0)+b2_0*(d2_0-d1_0));
	cv::Mat A1_i = A1.inv();
	cv::Mat cordinate1 = A1_i*A2;
	double tr12_x1 = cordinate1.at<double>(0,0);
	double tr12_x2 = cordinate1.at<double>(1,0);
	double tr12_y1 = a1_0*tr12_x1+c1_0;
	double tr12_y2 = a2_0*tr12_x2+c2_0;
	double tr12_z1 = b1_0*tr12_x1+d1_0;
	double tr12_z2 = b2_0*tr12_x2+d2_0;
	double tr12_x1_true = (tr12_x1+tr12_x2)/2;//track1とtrack2の最短距離を構成する座標の中点
	double tr12_y1_true = (tr12_y1+tr12_y2)/2;
	double tr12_z1_true = (tr12_z1+tr12_z2)/2;
	
	

	cv::Mat B1 = (cv::Mat_<double>(2,2) << pow(a1_0,2)+pow(b1_0,2)+1,-(a1_0*a3_0+b1_0*b3_0+1),
										(a1_0*a3_0+b1_0*b3_0+1),-(pow(a3_0,2)+pow(b3_0,2)+1));
	cv::Mat B2 = (cv::Mat_<double>(2,1) << a1_0*(c3_0-c1_0)+b1_0*(d3_0-d1_0),
										  a3_0*(c3_0-c1_0)+b3_0*(d3_0-d1_0));
	cv::Mat B1_i = B1.inv();
	cv::Mat cordinate2 = B1_i*B2;
	double tr13_x1 = cordinate2.at<double>(0,0);
	double tr13_x3 = cordinate2.at<double>(1,0);
	double tr13_y1 = a1_0*tr13_x1+c1_0;
	double tr13_y3 = a3_0*tr13_x3+c3_0;
	double tr13_z1 = b1_0*tr13_x1+d1_0;
	double tr13_z3 = b3_0*tr13_x3+d3_0;
	double tr13_x2_true = (tr13_x1+tr13_x3)/2;//track1とtrack3の最短距離を構成する座標の中点
	double tr13_y2_true = (tr13_y1+tr13_y3)/2;
	double tr13_z2_true = (tr13_z1+tr13_z3)/2;
	
	
	
	cv::Mat C1 = (cv::Mat_<double>(2,2) << pow(a2_0,2)+pow(b2_0,2)+1,-(a2_0*a3_0+b2_0*b3_0+1),
										(a2_0*a3_0+b2_0*b3_0+1),-(pow(a3_0,2)+pow(b3_0,2)+1));
	cv::Mat C2 = (cv::Mat_<double>(2,1) << a2_0*(c3_0-c2_0)+b2_0*(d3_0-d2_0),
										  a3_0*(c3_0-c2_0)+b3_0*(d3_0-d2_0));
	cv::Mat C1_i = C1.inv();
	cv::Mat cordinate3 = C1_i*C2;
	double tr23_x2 = cordinate3.at<double>(0,0);
	double tr23_x3 = cordinate3.at<double>(1,0);
	double tr23_y2 = a2_0*tr23_x2+c2_0;
	double tr23_y3 = a3_0*tr23_x3+c3_0;
	double tr23_z2 = b2_0*tr23_x2+d2_0;
	double tr23_z3 = b3_0*tr23_x3+d3_0;
	double tr23_x3_true = (tr23_x2+tr23_x3)/2;//track2とtrack3の最短距離を構成する座標の中点
	double tr23_y3_true = (tr23_y2+tr23_y3)/2;
	double tr23_z3_true = (tr23_z2+tr23_z3)/2;

	double center_x = (tr12_x1_true + tr13_x2_true + tr23_x3_true)/3;
	double center_y = (tr12_y1_true + tr13_y2_true + tr23_y3_true)/3;
	double center_z = (tr12_z1_true + tr13_z2_true + tr23_z3_true)/3;//vertex_pointの初期値

	//-35.3689;//
	// 25.7359;//
	// 15.3208;//

//	a1_0 = -0.75;
//	b1_0 = 2.83;
//	a2_0 = -0.08;
//	b2_0 = 0.63;
//	a3_0 = -0.41;
//	b3_0 = 1.08;


//	printf("%.6f  ,  %.6f  ,  %.6f  ,  %.6f\n",a1_0,b1_0,c1_0,d1_0);
//	printf("%.6f  ,  %.6f  ,  %.6f  ,  %.6f\n",a2_0,b2_0,c2_0,d2_0);
//	printf("%.6f  ,  %.6f  ,  %.6f  ,  %.6f\n",a3_0,b3_0,c3_0,d3_0);
	
	
//	printf("%.3f  ,  %.3f  ,  %.3f\n",center_x,center_y,center_z);
	
//	std::cout << "" << cordinate1 << std::endl << std::endl;
//	getchar();
	
	
	ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );//最小値を探す方法

	min.SetMaxFunctionCalls(10000000);
	min.SetMaxIterations(10000000);
	min.SetTolerance(0.00000001);

	ROOT::Math::Functor f(&chi,9);
	double step[9] = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.00001,0.00001,0.00001};//unknown parameter の動かす間隔
	//double step[9] = {1,1,1,1,1,1,1,1,1};//unknown parameter の動かす間隔
	//double variable[9] = {a1_0,b1_0,a2_0,b2_0,a3_0,b3_0,x0_0,y0_0,z0_0};//unknown parameter の初期値
	double variable[9] = {a1_0,b1_0,a2_0,b2_0,a3_0,b3_0,center_x,center_y,center_z};//unknown parameter の初期値
	//double variable[9] = {1,1,-1,1,-1,-1,0.0,0.0,0.0};//unknown parameter の初期値


	min.SetFunction(f);

	std::string A;
	std::string B;
	std::string C;
	std::string D;
	std::string E;
	std::string F;
	std::string G;
	std::string H;
	std::string I;

	min.SetLimitedVariable(0,A,variable[0],step[0],a1_0-2.0,a1_0+2.0);//番号,名前,初期値,間隔,最小値,最大値
	min.SetLimitedVariable(1,B,variable[1],step[1],b1_0-2.0,b1_0+2.0);
	min.SetLimitedVariable(2,C,variable[2],step[2],a2_0-2.0,a2_0+2.0);
	min.SetLimitedVariable(3,D,variable[3],step[3],b2_0-2.0,b2_0+2.0);
	min.SetLimitedVariable(4,E,variable[4],step[4],a3_0-2.0,a3_0+2.0);
	min.SetLimitedVariable(5,F,variable[5],step[5],b3_0-2.0,b3_0+2.0);
	min.SetLimitedVariable(6,G,variable[6],step[6],center_x-2.0,center_x+2.0);
	min.SetLimitedVariable(7,H,variable[7],step[7],center_y-2.0,center_y+2.0);
	min.SetLimitedVariable(8,I,variable[8],step[8],center_z-2.0,center_z+2.0);



	min.Minimize();

	const double *xs = min.X();

	aa1 = xs[0];
	bb1 = xs[1];
	aa2 = xs[2];
	bb2 = xs[3];
	aa3 = xs[4];
	bb3 = xs[5];
	xx0 = xs[6];
	yy0 = xs[7];
	zz0 = xs[8];
	result = chi(xs);

	



	/*
	cv::Mat Tx1 = cv::Mat(9,n1,CV_64F);
	cv::Mat Ty1 = cv::Mat(9,n1,CV_64F);
	cv::Mat Tz1 = cv::Mat(9,n1,CV_64F);

	for(int i=0;i<n1;i++)
	{
		Tx1.at<double>(0,i) = 1.0;
		Tx1.at<double>(1,i) = -aa1;
		Tx1.at<double>(2,i) = -bb1;
		Tx1.at<double>(3,i) = -(track1[i].y-yy0)/pow(track1[i].x-xx0,2);
		Tx1.at<double>(4,i) = -(track1[i].z-zz0)/pow(track1[i].x-xx0,2);
		Tx1.at<double>(5,i) = 0;
		Tx1.at<double>(6,i) = 0;
		Tx1.at<double>(7,i) = 0;
		Tx1.at<double>(8,i) = 0;
	}

	for(int i=0;i<n1;i++)
	{
		Ty1.at<double>(0,i) = -1.0/aa1;
		Ty1.at<double>(1,i) = 1;
		Ty1.at<double>(2,i) = 0;
		Ty1.at<double>(3,i) = 1/(track1[i].x-xx0);
		Ty1.at<double>(4,i) = 0;
		Ty1.at<double>(5,i) = 0;
		Ty1.at<double>(6,i) = 0;
		Ty1.at<double>(7,i) = 0;
		Ty1.at<double>(8,i) = 0;
	}

	for(int i=0;i<n1;i++)
	{
		Tz1.at<double>(0,i) = 0;
		Tz1.at<double>(1,i) = 0;
		Tz1.at<double>(2,i) = 1;
		Tz1.at<double>(3,i) = 0;
		Tz1.at<double>(4,i) = 1/(track1[i].x-xx0);
		Tz1.at<double>(5,i) = 0;
		Tz1.at<double>(6,i) = 0;
		Tz1.at<double>(7,i) = 0;
		Tz1.at<double>(8,i) = 0;
	}

	cv::Mat Tx2 = cv::Mat(9,n2,CV_64F);
	cv::Mat Ty2 = cv::Mat(9,n2,CV_64F);
	cv::Mat Tz2 = cv::Mat(9,n2,CV_64F);
	
	for(int i=0;i<n2;i++)
	{
		Tx2.at<double>(0,i) = 1.0;
		Tx2.at<double>(1,i) = -aa2;
		Tx2.at<double>(2,i) = -bb2;
		Tx2.at<double>(3,i) = 0;
		Tx2.at<double>(4,i) = 0;
		Tx2.at<double>(5,i) = -(track2[i].y-yy0)/pow(track2[i].x-xx0,2);
		Tx2.at<double>(6,i) = -(track2[i].z-zz0)/pow(track2[i].x-xx0,2);
		Tx2.at<double>(7,i) = 0;
		Tx2.at<double>(8,i) = 0;
	}

	for(int i=0;i<n2;i++)
	{
		Ty2.at<double>(0,i) = -1.0/aa2;
		Ty2.at<double>(1,i) = 1;
		Ty2.at<double>(2,i) = 0;
		Ty2.at<double>(3,i) = 0;
		Ty2.at<double>(4,i) = 0;
		Ty2.at<double>(5,i) = 1/(track2[i].x-xx0);
		Ty2.at<double>(6,i) = 0;
		Ty2.at<double>(7,i) = 0;
		Ty2.at<double>(8,i) = 0;
	}

	for(int i=0;i<n2;i++)
	{
		Tz2.at<double>(0,i) = 0;
		Tz2.at<double>(1,i) = 0;
		Tz2.at<double>(2,i) = 1;
		Tz2.at<double>(3,i) = 0;
		Tz2.at<double>(4,i) = 0;
		Tz2.at<double>(5,i) = 0;
		Tz2.at<double>(6,i) = 1/(track2[i].x-xx0);
		Tz2.at<double>(7,i) = 0;
		Tz2.at<double>(8,i) = 0;
	}




	cv::Mat Tx3 = cv::Mat(9,n3,CV_64F);
	cv::Mat Ty3 = cv::Mat(9,n3,CV_64F);
	cv::Mat Tz3 = cv::Mat(9,n3,CV_64F);
	
	for(int i=0;i<n3;i++)
	{
		Tx3.at<double>(0,i) = 1.0;
		Tx3.at<double>(1,i) = -aa3;
		Tx3.at<double>(2,i) = -bb3;
		Tx3.at<double>(3,i) = 0;
		Tx3.at<double>(4,i) = 0;
		Tx3.at<double>(5,i) = 0;
		Tx3.at<double>(6,i) = 0;
		Tx3.at<double>(7,i) = -(track3[i].y-yy0)/pow(track3[i].x-xx0,2);
		Tx3.at<double>(8,i) = -(track3[i].z-zz0)/pow(track3[i].x-xx0,2);
	}

	for(int i=0;i<n3;i++)
	{
		Ty3.at<double>(0,i) = -1.0/aa3;
		Ty3.at<double>(1,i) = 1;
		Ty3.at<double>(2,i) = 0;
		Ty3.at<double>(3,i) = 0;
		Ty3.at<double>(4,i) = 0;
		Ty3.at<double>(5,i) = 0;
		Ty3.at<double>(6,i) = 0;
		Ty3.at<double>(7,i) = 1/(track3[i].x-xx0);
		Ty3.at<double>(8,i) = 0;
	}

	for(int i=0;i<n3;i++)
	{
		Tz3.at<double>(0,i) = 0;
		Tz3.at<double>(1,i) = 0;
		Tz3.at<double>(2,i) = 1;
		Tz3.at<double>(3,i) = 0;
		Tz3.at<double>(4,i) = 0;
		Tz3.at<double>(5,i) = 0;
		Tz3.at<double>(6,i) = 0;
		Tz3.at<double>(7,i) = 0;
		Tz3.at<double>(8,i) = 1/(track3[i].x-xx0);
	}

	
	cv::Mat covMat = pow(err_x,2)*Tx1*Tx1.t()+pow(err_y,2)*Ty1*Ty1.t()+pow(err_z,2)*Tz1*Tz1.t()
					+pow(err_x,2)*Tx2*Tx2.t()+pow(err_y,2)*Ty2*Ty2.t()+pow(err_z,2)*Tz2*Tz2.t()
					+pow(err_x,2)*Tx3*Tx3.t()+pow(err_y,2)*Ty3*Ty3.t()+pow(err_z,2)*Tz3*Tz3.t();
	







	std::cout << covMat << std::endl << std::endl;

	*/


	
//	printf("       result = %.3f\n",result);

   return 0;
}




/*
bool Detamake1(std::vector<std::vector<coordinate>>& track_monte)
{
	track_monte.clear();
	coordinate tr;
	double	x,y,z;
	for(int j=0;j<monte_number;j++)//乱数の数
	{
		std::vector<coordinate> vcoord;
		for(int i=0;i<track1.size();i++)//点の数
		{
			tr.x = gRandom->Gaus(track1[i].x,err_x);
			tr.y = gRandom->Gaus(aa1*(track1[i].x-xx0)+yy0,err_y);
			tr.z = gRandom->Gaus(bb1*(track1[i].x-xx0)+zz0,err_z);
			vcoord.push_back(tr);
		}
		track_monte.push_back(vcoord);
	}
return true;
}
*/

bool Detamake1(std::vector<std::vector<coordinate>>& track_monte)
{
	track_monte.clear();
	coordinate tr;
	double	x,y,z;
	for(int j=0;j<monte_number;j++)//乱数の数
	{
		std::vector<coordinate> vcoord;
		for(int i=0;i<track1.size();i++)//点の数
		{
			double real_track1_x = (wx*(track1[i].x-xx0)+wy*aa1*(track1[i].y-yy0)+wz*bb1*(track1[i].z-zz0))/(wx+wy*pow(aa1,2)+wz*pow(bb1,2))+xx0;
			tr.x = gRandom->Gaus(real_track1_x,err_x);
			tr.y = gRandom->Gaus(aa1*(real_track1_x-xx0)+yy0,err_y);
			tr.z = gRandom->Gaus(bb1*(real_track1_x-xx0)+zz0,err_z);
			vcoord.push_back(tr);
		}
		track_monte.push_back(vcoord);
	}
return true;
}


bool Detamake2(std::vector<std::vector<coordinate>>& track_monte)
{
	track_monte.clear();
	coordinate tr;
	double	x,y,z;
	for(int j=0;j<monte_number;j++)//乱数の数
		{std::vector<coordinate> vcoord;
	for(int i=0;i<track2.size();i++)//点の数
	{
		double real_track2_x = (wx*(track2[i].x-xx0)+wy*aa2*(track2[i].y-yy0)+wz*bb2*(track2[i].z-zz0))/(wx+wy*pow(aa2,2)+wz*pow(bb2,2))+xx0;
		tr.x = gRandom->Gaus(real_track2_x,err_x);
		tr.y = gRandom->Gaus(aa2*(real_track2_x-xx0)+yy0,err_y);
		tr.z = gRandom->Gaus(bb2*(real_track2_x-xx0)+zz0,err_z);
		vcoord.push_back(tr);
		}
		track_monte.push_back(vcoord);
	}
return true;
}


bool Detamake3(std::vector<std::vector<coordinate>>& track_monte)
{
	track_monte.clear();
	coordinate tr;
	double	x,y,z;
	for(int j=0;j<monte_number;j++)//乱数の数
		{std::vector<coordinate> vcoord;
	for(int i=0;i<track3.size();i++)//点の数
	{
		double real_track3_x = (wx*(track3[i].x-xx0)+wy*aa3*(track3[i].y-yy0)+wz*bb3*(track3[i].z-zz0))/(wx+wy*pow(aa3,2)+wz*pow(bb3,2))+xx0;
		tr.x = gRandom->Gaus(real_track3_x,err_x);
		tr.y = gRandom->Gaus(aa3*(real_track3_x-xx0)+yy0,err_y);
		tr.z = gRandom->Gaus(bb3*(real_track3_x-xx0)+zz0,err_z);
		vcoord.push_back(tr);
		}
		track_monte.push_back(vcoord);
	}
return true;
}












int NumericalMinimization_Monte(std::vector<coordinate>& track1_monte_data,std::vector<coordinate>& track2_monte_data,std::vector<coordinate>& track3_monte_data,std::vector<results>& element,int h)
{	
	//element.clear();
	results re;






	//track1//
	coordinate sum1;
	
		sum1.x=0.0;
		sum1.y=0.0;
		sum1.z=0.0;
		
	double sum1_x_y=0.0;
	double sum1_x_z=0.0;
	double sum1_x_x=0.0;
	double sum1_y_y=0.0;
	double sum1_z_z=0.0;

	for(unsigned int i= 0; i<track1_monte_data.size(); i++)
	{
			sum1.x += track1_monte_data[i].x;
			sum1.y += track1_monte_data[i].y;
			sum1.z += track1_monte_data[i].z;
			sum1_x_y += track1_monte_data[i].x * track1_monte_data[i].y;
			sum1_x_z += track1_monte_data[i].x * track1_monte_data[i].z;
			sum1_x_x += track1_monte_data[i].x * track1_monte_data[i].x;
			sum1_y_y += track1_monte_data[i].y * track1_monte_data[i].y;
			sum1_z_z += track1_monte_data[i].z * track1_monte_data[i].z;
	}

	int n1 = track1_monte_data.size();//データの数
	//trackの傾き
	double a1_0 = (n1*sum1_x_y-sum1.y*sum1.x)/(n1*sum1_x_x-pow(sum1.x,2));
	double b1_0 = (n1*sum1_x_z-sum1.z*sum1.x)/(n1*sum1_x_x-pow(sum1.x,2));
	//trackの切片
	double c1_0 = (sum1.y*sum1_x_x-sum1_x_y*sum1.x)/(n1*sum1_x_x-pow(sum1.x,2));
	double d1_0 = (sum1.z*sum1_x_x-sum1_x_z*sum1.x)/(n1*sum1_x_x-pow(sum1.x,2));






	//track2//
	coordinate sum2;
	
		sum2.x=0.0;
		sum2.y=0.0;
		sum2.z=0.0;
		
	double sum2_x_y=0.0;
	double sum2_x_z=0.0;
	double sum2_x_x=0.0;
	double sum2_y_y=0.0;
	double sum2_z_z=0.0;

	for(unsigned int j= 0;j<track2_monte_data.size();j++)
	{
			sum2.x += track2_monte_data[j].x;
			sum2.y += track2_monte_data[j].y;
			sum2.z += track2_monte_data[j].z;
			sum2_x_y += track2_monte_data[j].x * track2_monte_data[j].y;
			sum2_x_z += track2_monte_data[j].x * track2_monte_data[j].z;
			sum2_x_x += track2_monte_data[j].x * track2_monte_data[j].x;
			sum2_y_y += track2_monte_data[j].y * track2_monte_data[j].y;
			sum2_z_z += track2_monte_data[j].z * track2_monte_data[j].z;
	}

	int n2 = track2_monte_data.size();//データの数
	//trackの傾き
	double a2_0 = (n2*sum2_x_y-sum2.y*sum2.x)/(n2*sum2_x_x-pow(sum2.x,2));
	double b2_0 = (n2*sum2_x_z-sum2.z*sum2.x)/(n2*sum2_x_x-pow(sum2.x,2));
	//trackの切片
	double c2_0 = (sum2.y*sum2_x_x-sum2_x_y*sum2.x)/(n2*sum2_x_x-pow(sum2.x,2));
	double d2_0 = (sum2.z*sum2_x_x-sum2_x_z*sum2.x)/(n2*sum2_x_x-pow(sum2.x,2));







	//track3//
	coordinate sum3;
	
		sum3.x=0.0;
		sum3.y=0.0;
		sum3.z=0.0;
		
	double sum3_x_y=0.0;
	double sum3_x_z=0.0;
	double sum3_x_x=0.0;
	double sum3_y_y=0.0;
	double sum3_z_z=0.0;

	for(unsigned int k= 0;k<track3_monte_data.size();k++)
	{
			sum3.x += track3_monte_data[k].x;
			sum3.y += track3_monte_data[k].y;
			sum3.z += track3_monte_data[k].z;
			sum3_x_y += track3_monte_data[k].x * track3_monte_data[k].y;
			sum3_x_z += track3_monte_data[k].x * track3_monte_data[k].z;
			sum3_x_x += track3_monte_data[k].x * track3_monte_data[k].x;
			sum3_y_y += track3_monte_data[k].y * track3_monte_data[k].y;
			sum3_z_z += track3_monte_data[k].z * track3_monte_data[k].z;
	}

	int n3 = track3_monte_data.size();//データの数
	//trackの傾き
	double a3_0 = (n3*sum3_x_y-sum3.y*sum3.x)/(n3*sum3_x_x-pow(sum3.x,2));
	double b3_0 = (n3*sum3_x_z-sum3.z*sum3.x)/(n3*sum3_x_x-pow(sum3.x,2));
	//trackの切片
	double c3_0 = (sum3.y*sum3_x_x-sum3_x_y*sum3.x)/(n3*sum3_x_x-pow(sum3.x,2));
	double d3_0 = (sum3.z*sum3_x_x-sum3_x_z*sum3.x)/(n3*sum3_x_x-pow(sum3.x,2));



	//vertex point の初期値
	
	cv::Mat A1 = (cv::Mat_<double>(2,2) << pow(a1_0,2)+pow(b1_0,2)+1,-(a1_0*a2_0+b1_0*b2_0+1),
										(a1_0*a2_0+b1_0*b2_0+1),-(pow(a2_0,2)+pow(b2_0,2)+1));
	cv::Mat A2 = (cv::Mat_<double>(2,1) << a1_0*(c2_0-c1_0)+b1_0*(d2_0-d1_0),
										  a2_0*(c2_0-c1_0)+b2_0*(d2_0-d1_0));
	cv::Mat A1_i = A1.inv();
	cv::Mat cordinate1 = A1_i*A2;
	double tr12_x1 = cordinate1.at<double>(0,0);
	double tr12_x2 = cordinate1.at<double>(1,0);
	double tr12_y1 = a1_0*tr12_x1+c1_0;
	double tr12_y2 = a2_0*tr12_x2+c2_0;
	double tr12_z1 = b1_0*tr12_x1+d1_0;
	double tr12_z2 = b2_0*tr12_x2+d2_0;
	double tr12_x1_true = (tr12_x1+tr12_x2)/2;//track1とtrack2の最短距離を構成する座標の中点
	double tr12_y1_true = (tr12_y1+tr12_y2)/2;
	double tr12_z1_true = (tr12_z1+tr12_z2)/2;
	
	

	cv::Mat B1 = (cv::Mat_<double>(2,2) << pow(a1_0,2)+pow(b1_0,2)+1,-(a1_0*a3_0+b1_0*b3_0+1),
										(a1_0*a3_0+b1_0*b3_0+1),-(pow(a3_0,2)+pow(b3_0,2)+1));
	cv::Mat B2 = (cv::Mat_<double>(2,1) << a1_0*(c3_0-c1_0)+b1_0*(d3_0-d1_0),
										  a3_0*(c3_0-c1_0)+b3_0*(d3_0-d1_0));
	cv::Mat B1_i = B1.inv();
	cv::Mat cordinate2 = B1_i*B2;
	double tr13_x1 = cordinate2.at<double>(0,0);
	double tr13_x3 = cordinate2.at<double>(1,0);
	double tr13_y1 = a1_0*tr13_x1+c1_0;
	double tr13_y3 = a3_0*tr13_x3+c3_0;
	double tr13_z1 = b1_0*tr13_x1+d1_0;
	double tr13_z3 = b3_0*tr13_x3+d3_0;
	double tr13_x2_true = (tr13_x1+tr13_x3)/2;//track1とtrack3の最短距離を構成する座標の中点
	double tr13_y2_true = (tr13_y1+tr13_y3)/2;
	double tr13_z2_true = (tr13_z1+tr13_z3)/2;
	
	
	
	cv::Mat C1 = (cv::Mat_<double>(2,2) << pow(a2_0,2)+pow(b2_0,2)+1,-(a2_0*a3_0+b2_0*b3_0+1),
										(a2_0*a3_0+b2_0*b3_0+1),-(pow(a3_0,2)+pow(b3_0,2)+1));
	cv::Mat C2 = (cv::Mat_<double>(2,1) << a2_0*(c3_0-c2_0)+b2_0*(d3_0-d2_0),
										  a3_0*(c3_0-c2_0)+b3_0*(d3_0-d2_0));
	cv::Mat C1_i = C1.inv();
	cv::Mat cordinate3 = C1_i*C2;
	double tr23_x2 = cordinate3.at<double>(0,0);
	double tr23_x3 = cordinate3.at<double>(1,0);
	double tr23_y2 = a2_0*tr23_x2+c2_0;
	double tr23_y3 = a3_0*tr23_x3+c3_0;
	double tr23_z2 = b2_0*tr23_x2+d2_0;
	double tr23_z3 = b3_0*tr23_x3+d3_0;
	double tr23_x3_true = (tr23_x2+tr23_x3)/2;//track1とtrack3の最短距離を構成する座標の中点
	double tr23_y3_true = (tr23_y2+tr23_y3)/2;
	double tr23_z3_true = (tr23_z2+tr23_z3)/2;

	double center_x = (tr12_x1_true + tr13_x2_true + tr23_x3_true)/3;
	double center_y = (tr12_y1_true + tr13_y2_true + tr23_y3_true)/3;
	double center_z = (tr12_z1_true + tr13_z2_true + tr23_z3_true)/3;//vertex_pointの初期値

	
	
	ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );//最小値を探す方法

	min.SetMaxFunctionCalls(10000000);
	min.SetMaxIterations(10000000);
	min.SetTolerance(0.00000001);

	ROOT::Math::Functor f(&chi_monte,9);
	double step[9] = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.00001,0.00001,0.00001};//unknown parameter の動かす間隔
	double variable[9] = {a1_0,b1_0,a2_0,b2_0,a3_0,b3_0,center_x,center_y,center_z};//unknown parameter の初期値

	min.SetFunction(f);

	std::string A;
	std::string B;
	std::string C;
	std::string D;
	std::string E;
	std::string F;
	std::string G;
	std::string H;
	std::string I;

	min.SetLimitedVariable(0,A,variable[0],step[0],a1_0-2.0,a1_0+2.0);//番号,名前,初期値,間隔,最小値,最大値
	min.SetLimitedVariable(1,B,variable[1],step[1],b1_0-2.0,b1_0+2.0);
	min.SetLimitedVariable(2,C,variable[2],step[2],a2_0-2.0,a2_0+2.0);
	min.SetLimitedVariable(3,D,variable[3],step[3],b2_0-2.0,b2_0+2.0);
	min.SetLimitedVariable(4,E,variable[4],step[4],a3_0-2.0,a3_0+2.0);
	min.SetLimitedVariable(5,F,variable[5],step[5],b3_0-2.0,b3_0+2.0);
	min.SetLimitedVariable(6,G,variable[6],step[6],center_x-2.0,center_x+2.0);
	min.SetLimitedVariable(7,H,variable[7],step[7],center_y-2.0,center_y+2.0);
	min.SetLimitedVariable(8,I,variable[8],step[8],center_z-2.0,center_z+2.0);



	min.Minimize();

	const double *xs = min.X();

	re.a1 = xs[0];
	re.b1 = xs[1];
	re.a2 = xs[2];
	re.b2 = xs[3];
	re.a3 = xs[4];
	re.b3 = xs[5];
	re.x0 = xs[6];
	re.y0 = xs[7];
	re.z0 = xs[8];
	re.chi_2 = chi(xs);
	element.push_back(re);
	
	number = number + 1;
	


   return 0;



}



int main()
{
	

	std::vector<results> element;
	//element.reserve(1000);

	//最小値アルゴリズムへ挿入
	NumericalMinimization();

	Detamake1(track1_monte);
	Detamake2(track2_monte);
	Detamake3(track3_monte);

	printf("Data OK\n");



	TCanvas* c1 = new TCanvas("c1");
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat("neiuoMR"); 
	gStyle->SetPalette(1);
	double x_max=50.0;
	double x_min=20.0;
	TH1D *h1;
	h1= new TH1D("h1","h1 sample", 50, x_min, x_max);
	h1->SetXTitle("X-axis");
	h1->SetYTitle("Y-axis");
	
	

	for(int h=0;h<track1_monte.size();h++)
	{
		NumericalMinimization_Monte(track1_monte[h],track2_monte[h],track3_monte[h],element,h);
		
		if(h==5000)
		{
		printf("%d OK\n",h);
		}
	}

	double sum_a1 = 0;
	double sum_b1 = 0;
	double sum_a2 = 0;
	double sum_b2 = 0;
	double sum_a3 = 0;
	double sum_b3 = 0;
	double sum_x0 = 0;
	double sum_y0 = 0;
	double sum_z0 = 0;
	double sum_chi = 0;

	double dvi_a1 = 0;
	double dvi_b1 = 0;
	double dvi_a2 = 0;
	double dvi_b2 = 0;
	double dvi_a3 = 0;
	double dvi_b3 = 0;
	double dvi_x0 = 0;
	double dvi_y0 = 0;
	double dvi_z0 = 0;
	double dvi_chi = 0;

	for(int i=0;i<element.size();i++)
	{
		sum_a1 += element[i].a1;
		sum_b1 += element[i].b1;
		sum_a2 += element[i].a2;
		sum_b2 += element[i].b2;
		sum_a3 += element[i].a3;
		sum_b3 += element[i].b3;
		sum_x0 += element[i].x0;
		sum_y0 += element[i].y0;
		sum_z0 += element[i].z0;
		sum_chi += element[i].chi_2;
	}

	double average_a1 = sum_a1/element.size();
	double average_b1 = sum_b1/element.size();
	double average_a2 = sum_a2/element.size();
	double average_b2 = sum_b2/element.size();
	double average_a3 = sum_a3/element.size();
	double average_b3 = sum_b3/element.size();
	double average_x0 = sum_x0/element.size();
	double average_y0 = sum_y0/element.size();
	double average_z0 = sum_z0/element.size();
	double average_chi = sum_chi/element.size();


	for(int i=0;i<element.size();i++)
	{
		dvi_a1 += pow(element[i].a1-average_a1,2);
		dvi_b1 += pow(element[i].b1-average_b1,2);
		dvi_a2 += pow(element[i].a2-average_a2,2);
		dvi_b2 += pow(element[i].b2-average_b2,2);
		dvi_a3 += pow(element[i].a3-average_a3,2);
		dvi_b3 += pow(element[i].b3-average_b3,2);
		dvi_x0 += pow(element[i].x0-average_x0,2);
		dvi_y0 += pow(element[i].y0-average_y0,2);
		dvi_z0 += pow(element[i].z0-average_z0,2);
		dvi_chi += pow(element[i].chi_2-average_chi,2);
	}

	dvi_a1 = sqrt(dvi_a1/element.size());
	dvi_b1 = sqrt(dvi_b1/element.size());
	dvi_a2 = sqrt(dvi_a2/element.size());
	dvi_b2 = sqrt(dvi_b2/element.size());
	dvi_a3 = sqrt(dvi_a3/element.size());
	dvi_b3 = sqrt(dvi_b3/element.size());
	dvi_x0 = sqrt(dvi_x0/element.size());
	dvi_y0 = sqrt(dvi_y0/element.size());
	dvi_z0 = sqrt(dvi_z0/element.size());
	dvi_chi = sqrt(dvi_chi/element.size());

	printf("a1 = %.9f\n",average_a1);
	printf("b1 = %.9f\n",average_b1);
	printf("a2 = %.9f\n",average_a2);
	printf("b2 = %.9f\n",average_b2);
	printf("a3 = %.9f\n",average_a3);
	printf("b3 = %.9f\n",average_b3);
	printf("x0 = %.9f\n",average_x0);
	printf("y0 = %.9f\n",average_y0);
	printf("z0 = %.9f\n",average_z0);
	printf("chi = %.9f\n\n\n\n",average_chi);



	printf("      a1 = %.9f  ± %.9f\n",aa1,dvi_a1);
	printf("      b1 = %.9f  ± %.9f\n",bb1,dvi_b1);
	printf("      a2 = %.9f  ± %.9f\n",aa2,dvi_a2);
	printf("      b2 = %.9f  ± %.9f\n",bb2,dvi_b2);
	printf("      a3 = %.9f  ± %.9f\n",aa3,dvi_a3);
	printf("      b3 = %.9f  ± %.9f\n",bb3,dvi_b3);
	printf("      x0 = %.9f  ± %.9f\n",xx0,dvi_x0);
	printf("      y0 = %.9f  ± %.9f\n",yy0,dvi_y0);
	printf("      z0 = %.9f  ± %.9f\n",zz0,dvi_z0);
	printf("	 χ^2 = %.20f\n",result);




	for(int f=0;f<track1_monte.size();f++)
	{
		h1->Fill(element[f].chi_2);
		//h1->Fill(track1_monte[f][0].x);
	}
	
//	double chisquared_pdf(double x, double r, double x0) {
//      
//      if ((x-x0) <  0) {
//         return 0.0;
//      }
//      double a = r/2 -1.; 
//      // let return inf for case x  = x0 and treat special case of r = 2 otherwise will return nan
//      if (x == x0 && a == 0) return 0.5;
//
//      return std::exp ((r/2 - 1) * std::log((x-x0)/2) - (x-x0)/2 - ROOT::Math::lgamma(r/2))/2;
//      
//   }

		
	
	
	h1->Draw();
	c1->Print("chi.png");

	



	printf("計算できた");
	getchar();
	return 0;
}