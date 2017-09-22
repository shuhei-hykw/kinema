#include<stdio.h>
#include"opencv2\\opencv.hpp"
#include<string>
#include<math.h>
#include <iostream>




struct coordinate {
	double x, y, z;
};



bool ReadXYZData(char* filename, std::vector<coordinate>& vp) {
	FILE *data_in;
	fopen_s(&data_in, filename, "r");
	double	x, y, z;
	while (fscanf_s(data_in, "%lf %lf %lf", &x, &y, &z) != EOF) {
		coordinate p;
		p.x = x;
		p.y = y;
		p.z = z;
		vp.emplace_back(p);
	}
	return true;
}


int main()
{
	std::vector<coordinate> vp;
	ReadXYZData("10.txt", vp);


	double sumX = 0.0;
	double sumY = 0.0;
	double sumZ = 0.0;
	double sumXY = 0.0;
	double sumXZ = 0.0;
	double sumXX = 0.0;
	double sumYY = 0.0;
	double sumZZ = 0.0;
	double errX = 0.00025;
	double errY = 0.00025;
	double errZ = 0.00036;



	for (int i = 0; i < vp.size(); i++)
	{
		sumX += vp[i].x;
		sumY += vp[i].y;
		sumZ += vp[i].z;
		sumXY += vp[i].x * vp[i].y;
		sumXZ += vp[i].x * vp[i].z;
		sumXX += vp[i].x * vp[i].x;
		sumYY += vp[i].y * vp[i].y;
		sumZZ += vp[i].z * vp[i].z;
	}
	int n = vp.size();//データの数


	//Y = alpha*X + beta
	//Z = gamma*X + delta
	//least square method
	//http://hooktail.sub.jp/computPhys/least-square/
	double alpha0 = (n*sumXY - sumY*sumX) / (n*sumXX - pow(sumX, 2));
	double beta0 = (sumY*sumXX - sumXY*sumX) / (n*sumXX - pow(sumX, 2));
	double gamma0 = (n*sumXZ - sumZ*sumX) / (n*sumXX - pow(sumX, 2));
	double delta0 = (sumZ*sumXX - sumXZ*sumX) / (n*sumXX - pow(sumX, 2));
	cv::Mat z0 = (cv::Mat_<double>(4, 1) << alpha0, beta0, gamma0, delta0);



	cv::Mat matEye = cv::Mat::eye(n, n, CV_64F);
	cv::Mat matZero = cv::Mat::zeros(n, n, CV_64F);
	cv::Mat up_left = -alpha0*matEye;
	cv::Mat down_left = -gamma0*matEye;

	//known parameter のヤコビ行列
	cv::Mat D = cv::Mat(2 * n, 3 * n, CV_64F);
	cv::Mat roi1 = D(cv::Rect(0, 0, n, n));//Ｄ行列の左上
	up_left.copyTo(roi1);
	cv::Mat roi2 = D(cv::Rect(0, n, n, n));//Ｄ行列の左下
	down_left.copyTo(roi2);
	cv::Mat roi3 = D(cv::Rect(n, 0, n, n));//Ｄ行列の真ん中上
	matEye.copyTo(roi3);
	cv::Mat roi4 = D(cv::Rect(n, n, n, n));//Ｄ行列の真ん中下
	matZero.copyTo(roi4);
	cv::Mat roi5 = D(cv::Rect(2 * n, 0, n, n));//Ｄ行列の右上
	matZero.copyTo(roi5);
	cv::Mat roi6 = D(cv::Rect(2 * n, n, n, n));//Ｄ行列の右下
	matEye.copyTo(roi6);


	//unknown parameter のヤコビ行列
	cv::Mat E = cv::Mat(2 * n, 4, CV_64F);


	for (int j = 0; j < vp.size(); j++)
	{
		E.at<double>(j, 0) = -1 * vp[j].x;
		E.at<double>(j, 1) = -1;
		E.at<double>(j, 2) = 0;
		E.at<double>(j, 3) = 0;
		E.at<double>(j + n, 0) = 0;
		E.at<double>(j + n, 1) = 0;
		E.at<double>(j + n, 2) = -1 * vp[j].x;
		E.at<double>(j + n, 3) = -1;
	}




	//constraint
	cv::Mat d = cv::Mat(2 * n, 1, CV_64F);
	for (int k = 0; k < vp.size(); k++)
	{
		d.at<double>(k, 0) = alpha0*vp[k].x;
		d.at<double>(k + n, 0) = gamma0*vp[k].x;
	}



	//parameter
	cv::Mat a0 = cv::Mat(3 * n, 1, CV_64F);
	for (int h = 0; h < vp.size(); h++)
	{
		a0.at<double>(h, 0) = vp[h].x;
		a0.at<double>(h + n, 0) = vp[h].y;
		a0.at<double>(h + 2 * n, 0) = vp[h].z;
	}


	//covariance matrix
	cv::Mat Va = cv::Mat::zeros(3 * n, 3 * n, CV_64F);
	for (int g = 0; g < vp.size(); g++)
	{
		Va.at<double>(g, g) = errX*errX;
		Va.at<double>(g + n, g + n) = errY*errY;
		Va.at<double>(g + 2 * n, g + 2 * n) = errZ*errZ;
	}


	//    cv::Mat d = d0 - D*a0 - E*z0;
	cv::Mat Dt = D.t();
	cv::Mat Vd0 = D*Va*Dt;
	cv::Mat Vd = Vd0.inv();

	cv::Mat E_t = E.t();
	cv::Mat VE0 = E_t*Vd*E;
	cv::Mat VE = VE0.inv();
	cv::Mat lambda0 = Vd*(D*a0 + d);

	cv::Mat z = -VE*E_t*lambda0;
	cv::Mat lambda = lambda0 + Vd*E*z;


	cv::Mat lambda_t = lambda.t();
	cv::Mat Vd_i = Vd.inv();
	cv::Mat chi_square = lambda_t*Vd_i*lambda;


	cv::Mat a = a0 - Va*Dt*lambda;
	cv::Mat Vz = VE;
	cv::Mat V_lambda = Vd - Vd*E*VE*E_t*Vd;
	cv::Mat V = Va - Va*Dt*V_lambda*D*Va;

	double pi = 6 * asin(0.5);

	double alpha = z.at<double>(0, 0);
	double beta = z.at<double>(1, 0);
	double gamma = z.at<double>(2, 0);
	double delta = z.at<double>(3, 0);

	double Variance_alpha = Vz.at<double>(0, 0);
	double Variance_beta = Vz.at<double>(1, 1);
	double Variance_gamma = Vz.at<double>(2, 2);
	double Variance_delta = Vz.at<double>(3, 3);
	double Variance_alpha_gamma = Vz.at<double>(2, 0);
	double Variance_alpha_beta = Vz.at<double>(0, 1);
	double Variance_gamma_delta = Vz.at<double>(2, 3);


	printf("alpha = %.7f +/- %.7f\n", alpha, sqrt(Variance_alpha));
	printf("beta = %.7f +/- %.7f\n", beta, sqrt(Variance_beta));
	printf("gamma = %.7f +/- %.7f\n", gamma, sqrt(Variance_gamma));
	printf("delta = %.7f +/- %.7f\n", delta, sqrt(Variance_delta));
	printf("err_alpha_beta = %.7f\n", Variance_alpha_beta);
	printf("err_gamma_delta = %.7f\n", Variance_gamma_delta);





	double x1 = vp[0].x;
	double x2 = vp[vp.size() - 1].x;
	double y1 = beta + alpha*x1;
	double y2 = beta + alpha*x2;
	double z1 = delta + gamma*x1;
	double z2 = delta + gamma*x2;
	double dx = x2 - x1;
	double dy = y2 - y1;
	double dz = z2 - z1;


	double range = sqrt(dx*dx + dy*dy + dz*dz);
	double rangeXYproj = sqrt(dx*dx + dy*dy);

	double sin_phi = dy / rangeXYproj;

	double theta = acos(dz / range) * 180 / pi;
	double phi = acos(dx / rangeXYproj) * 180 / pi;
	if (sin_phi<0) phi = 360 - phi;





	double errY1 = sqrt(x1*x1*Variance_alpha + Variance_beta + alpha*alpha*errX*errX + 2 * x1*Variance_alpha_beta);
	double errY2 = sqrt(x2*x2*Variance_alpha + Variance_beta + alpha*alpha*errX*errX + 2 * x2*Variance_alpha_beta);

	double errZ1 = sqrt(x1*x1*Variance_gamma + Variance_delta + gamma*gamma*errX*errX + 2 * x1*Variance_gamma_delta);
	double errZ2 = sqrt(x2*x2*Variance_gamma + Variance_delta + gamma*gamma*errX*errX + 2 * x2*Variance_gamma_delta);


	double err_R = sqrt(
		2*errX*errX*(dx/range)*(dx/range)
		+ (errY1*errY1 + errY2*errY2)*(dy/range)*(dy/range)
		+ (errZ1*errZ1 + errZ2*errZ2)*(dz/range)*(dz/range)
	);

	double err_delta_z = sqrt(errZ2*errZ2 + errZ1*errZ1);





	//fitting parameter を使った角度の算出と場合分け
	double phi2 = atan(alpha) * 180 / pi;


	if (phi2<0)
	{
		(sin(dy / rangeXYproj)<0 ) ? phi2 += 360 : phi2 += 180;
	}


	double theta2 = atan(1 / (gamma*cos(phi2*pi / 180))) * 180 / pi;



	double err_phi2 = (sqrt(Variance_alpha) / (1 + alpha*alpha)) * 180 / pi;

	//tanθ = 1/(gamma*cosφ) = A
	double err_A = sqrt(
		pow(alpha / (gamma*sqrt(1 + alpha*alpha)), 2)*Variance_alpha
		+ pow(sqrt(1 + alpha*alpha) / pow(gamma, 2), 2)*Variance_gamma
		- 2 * (alpha / (gamma*sqrt(1 + alpha*alpha)))*(sqrt(1 + alpha*alpha) / pow(gamma, 2))*Variance_alpha_gamma
	);

	double err_theta2 = (1 / (1 + pow(tan(theta2*pi / 180), 2)))*err_A * 180 / pi;


	double A = sqrt(1 + alpha*alpha) / gamma;

	double err_theta_phi = (1 / (1 + alpha*alpha))*(1 / (1 + pow(A, 2)))*(alpha / (gamma*sqrt(1 + alpha*alpha)))*Variance_alpha
		+ (1 / (1 + alpha*alpha))*(1 / (1 + pow(A, 2)))*(-sqrt(1 + alpha*alpha) / gamma)*Variance_alpha_gamma;






	//角度算出のためのRange
	double sum_Range = 0;
	for (int i = 0; i<vp.size() - 1; i++)
	{
		sum_Range += sqrt(
			pow(vp[i + 1].x - vp[i].x, 2)
			+ pow(vp[i + 1].y - vp[i].y, 2)
			+ pow(vp[i + 1].z - vp[i].z, 2)
		);
	}



	printf("θ = %.1f +/- %.1f\n", theta, err_theta2);
	printf("φ = %.1f +/- %.1f\n", phi, err_phi2);

	printf("θ = %.1f +/- %.1f\n", theta2, err_theta2);
	printf("φ = %.1f +/- %.1f\n", phi2, err_phi2);

	printf("%.3f mm\n", sum_Range);
	printf("%.9f\n", err_theta_phi);


	return 0;
}
