double chi(double m1)
{

	double m_Nuclear, m_xi, err_m_xi;//mass of Carbon , xi
	double B_xi;//binding energy of xi
	double initial_m1, m2, m3;//mass of single Hyper or nuclei


	double momentum1, initial_momentum1, momentum2, momentum3;//momentum of single Hyper
	double err_momentum1, err_momentum2, err_momentum3;



	double theta1, theta2, theta3;
	double err_theta1, err_theta2, err_theta3;
	double theta1_rad, theta2_rad, theta3_rad;
	double err_theta1_rad, err_theta2_rad, err_theta3_rad;

	double phi1, phi2, phi3;
	double err_phi1, err_phi2, err_phi3;
	double phi1_rad, phi2_rad, phi3_rad;
	double err_phi1_rad, err_phi2_rad, err_phi3_rad;


	double err_theta1_phi1_rad, err_theta2_phi2_rad, err_theta3_phi3_rad;

	B_xi = 0.0;
	m_Nuclear = 11174.866;
	m_xi = 1321.71;
	err_m_xi = 0.07;

	initial_m1 = m1;
	m2 = 3727.380;
	m3 = 2808.922;

	momentum1 = 170.961;//173.65;//173.050;//173.118;//;//173.65;//170.453;//170.961;
	momentum2 = 82.661;//  86.19;//75.955;//75.920;//;//86.19;//82.142;Å@Å@Å@Å@82.661;
	momentum3 = 166.049;//166.53;//166.713;//166.629;//;//166.53;//164.454;//166.049;

	theta1 = 44.9;//44.91;//47.928;//47.926;//44.9;
	theta2 = 57.7;//57.74;//61.089;//60.937;//57.7;
	theta3 = 156.2;//156.18;//156.318;//156.324;//156.2;

	phi1 = 337.5;//337.47;//338.492;//338.524;//337.5;
	phi2 = 174.9;//174.91;//174.275;//174.236;//174.9;
	phi3 = 143.0;//143.01;//142.824;//142.808;//143.0;

	err_momentum1 = 3.506;//3.14;//1.083;//1.525;//3.14;//3.216;//3.506;//
	err_momentum2 = 5.855;//6.27;//1.244;//1.790;//6.27;//6.169;//5.855;//
	err_momentum3 = 0.666;//0.80;//1.244;//1.790;//0.80;//0.659;//0.666;//

	err_theta1 = 4.72;//2.94;//4.72;//2.01;//0.665;//1.218;//2.0;
	err_theta2 = 8.41;//6.19;//8.41;//5.24;//1.306;//1.990;//5.2;
	err_theta3 = 1.86;//1.01;//1.86;//0.46;//0.458;//0.490;//0.5;

	err_phi1 = 4.64;//2.80;//4.64;//1.80;//0.887;//1.282;//1.8;
	err_phi2 = 7.21;//4.38;//7.21;//2.87;//1.677;//2.139;//2.9;
	err_phi3 = 2.05;//1.37;//2.05;//1.03;//0.936;//0.979;//1.0;


	err_theta1_phi1_rad = 1.74231*pow(10, -5);//0;
	err_theta2_phi2_rad = 7.87901*pow(10, -4);//0;
	err_theta3_phi3_rad = 5.28957*pow(10, -6);//0;

	B_xi = 0;//0.13;





	double pi = 6 * asin(0.5);//â~é¸ó¶ÇÃíËã`


	//ÅãÇradÇ…ïœä∑
	theta1_rad = theta1*pi / 180;
	theta2_rad = theta2*pi / 180;
	theta3_rad = theta3*pi / 180;

	err_theta1_rad = err_theta1*pi / 180;
	err_theta2_rad = err_theta2*pi / 180;
	err_theta3_rad = err_theta3*pi / 180;

	phi1_rad = phi1*pi / 180;
	phi2_rad = phi2*pi / 180;
	phi3_rad = phi3*pi / 180;

	err_phi1_rad = err_phi1*pi / 180;
	err_phi2_rad = err_phi2*pi / 180;
	err_phi3_rad = err_phi3*pi / 180;




	double a0_11 = momentum1;
	double a0_21 = theta1_rad;
	double a0_31 = phi1_rad;
	double a0_41 = momentum2;
	double a0_51 = theta2_rad;
	double a0_61 = phi2_rad;
	double a0_71 = momentum3;
	double a0_81 = theta3_rad;
	double a0_91 = phi3_rad;
	double a0_101 = m_xi;


	cv::Mat a0 = (cv::Mat_<double>(10, 1) << a0_11, a0_21, a0_31, a0_41, a0_51, a0_61, a0_71, a0_81, a0_91, a0_101);  //data from emulsion matrix




	double D11 = momentum1 / sqrt(pow(m1, 2) + pow(momentum1, 2));
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

	double D14 = momentum2 / sqrt(pow(m2, 2) + pow(momentum2, 2));
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

	double D17 = momentum3 / sqrt(pow(m3, 2) + pow(momentum3, 2));
	double D27 = sin(theta3_rad)*cos(phi3_rad);
	double D37 = sin(theta3_rad)*sin(phi3_rad);
	double D47 = cos(theta3_rad);

	double D18 = 0;
	double D28 = momentum3*cos(theta3_rad)*cos(phi3_rad);
	double D38 = momentum3*cos(theta3_rad)*sin(phi3_rad);
	double D48 = -momentum3*sin(theta3_rad);

	double D19 = 0;
	double D29 = -momentum3*sin(theta3_rad)*sin(phi3_rad);
	double D39 = momentum3*sin(theta3_rad)*cos(phi3_rad);
	double D49 = 0;

	double D110 = -1;
	double D210 = 0;
	double D310 = 0;
	double D410 = 0;



	cv::Mat D = (cv::Mat_<double>(4, 10) << D11, D12, D13, D14, D15, D16, D17, D18, D19, D110,						//constraints derivation matrix
		D21, D22, D23, D24, D25, D26, D27, D28, D29, D210,
		D31, D32, D33, D34, D35, D36, D37, D38, D39, D310,
		D41, D42, D43, D44, D45, D46, D47, D48, D49, D410);




	//covariance  momentum-theta,momentum-phi,theta-phi ,particle1-particle2 = 0

	double Va11, Va12, Va13, Va14, Va15, Va16, Va17, Va18, Va19, Va110;
	double Va21, Va22, Va23, Va24, Va25, Va26, Va27, Va28, Va29, Va210;
	double Va31, Va32, Va33, Va34, Va35, Va36, Va37, Va38, Va39, Va310;
	double Va41, Va42, Va43, Va44, Va45, Va46, Va47, Va48, Va49, Va410;
	double Va51, Va52, Va53, Va54, Va55, Va56, Va57, Va58, Va59, Va510;
	double Va61, Va62, Va63, Va64, Va65, Va66, Va67, Va68, Va69, Va610;
	double Va71, Va72, Va73, Va74, Va75, Va76, Va77, Va78, Va79, Va710;
	double Va81, Va82, Va83, Va84, Va85, Va86, Va87, Va88, Va89, Va810;
	double Va91, Va92, Va93, Va94, Va95, Va96, Va97, Va98, Va99, Va910;
	double Va101, Va102, Va103, Va104, Va105, Va106, Va107, Va108, Va109, Va1010;

	Va14 = Va15 = Va16 = Va17 = Va18 = Va19 = Va110 = 0;    //äeó±éqä‘ÇÃã§ï™éUÅÅÇO
	Va24 = Va25 = Va26 = Va27 = Va28 = Va29 = Va210 = 0;
	Va34 = Va35 = Va36 = Va37 = Va38 = Va39 = Va310 = 0;
	Va47 = Va48 = Va49 = Va410 = 0;
	Va57 = Va58 = Va59 = Va510 = 0;
	Va67 = Va68 = Va69 = Va610 = 0;
	Va710 = 0;
	Va810 = 0;
	Va910 = 0;
	Va41 = Va42 = Va43 = 0;
	Va51 = Va52 = Va53 = 0;
	Va61 = Va62 = Va63 = 0;
	Va71 = Va72 = Va73 = 0;
	Va81 = Va82 = Va83 = 0;
	Va91 = Va92 = Va93 = 0;
	Va101 = Va102 = Va103 = 0;
	Va74 = Va75 = Va76 = 0;
	Va84 = Va85 = Va86 = 0;
	Va94 = Va95 = Va96 = 0;
	Va104 = Va105 = Va106 = 0;
	Va107 = Va108 = Va109 = 0;


	Va21 = Va12 = Va31 = Va13 = 0;     //ÇPÇ¬ÇÃó±éqÇÃâ^ìÆó Ç∆É∆ÅAÉ”ä‘ÇÃã§ï™éUÅÅÇO
	Va54 = Va45 = Va64 = Va46 = 0;
	Va87 = Va78 = Va97 = Va79 = 0;


	Va11 = pow(err_momentum1, 2);     //äeÉpÉâÉÅÅ[É^ÇÃï™éUÅAÉ∆ÅAÉ”ä‘ÇÃã§ï™éU
	Va22 = pow(err_theta1_rad, 2);
	Va33 = pow(err_phi1_rad, 2);
	Va32 = Va23 = err_theta1_phi1_rad;

	Va44 = pow(err_momentum2, 2);
	Va55 = pow(err_theta2_rad, 2);
	Va66 = pow(err_phi2_rad, 2);
	Va65 = Va56 = err_theta2_phi2_rad;

	Va77 = pow(err_momentum3, 2);
	Va88 = pow(err_theta3_rad, 2);
	Va99 = pow(err_phi3_rad, 2);
	Va98 = Va89 = err_theta3_phi3_rad;

	Va1010 = pow(err_m_xi, 2);


	cv::Mat Va = (cv::Mat_<double>(10, 10) << Va11, Va12, Va13, Va14, Va15, Va16, Va17, Va18, Va19, Va110,		//covariance matrix
		Va21, Va22, Va23, Va24, Va25, Va26, Va27, Va28, Va29, Va210,
		Va31, Va32, Va33, Va34, Va35, Va36, Va37, Va38, Va39, Va310,
		Va41, Va42, Va43, Va44, Va45, Va46, Va47, Va48, Va49, Va410,
		Va51, Va52, Va53, Va54, Va55, Va56, Va57, Va58, Va59, Va510,
		Va61, Va62, Va63, Va64, Va65, Va66, Va67, Va68, Va69, Va610,
		Va71, Va72, Va73, Va74, Va75, Va76, Va77, Va78, Va79, Va710,
		Va81, Va82, Va83, Va84, Va85, Va86, Va87, Va88, Va89, Va810,
		Va91, Va92, Va93, Va94, Va95, Va96, Va97, Va98, Va99, Va910,
		Va101, Va102, Va103, Va104, Va105, Va106, Va107, Va108, Va109, Va1010);



	double H1, H2, H3, H4;//constraints 
	H1 = sqrt(pow(m1, 2) + pow(momentum1, 2)) + sqrt(pow(m2, 2) + pow(momentum2, 2)) + sqrt(pow(m3, 2) + pow(momentum3, 2)) - (m_Nuclear + m_xi - B_xi);
	H2 = momentum1*sin(theta1_rad)*cos(phi1_rad) + momentum2*sin(theta2_rad)*cos(phi2_rad) + momentum3*sin(theta3_rad)*cos(phi3_rad);
	H3 = momentum1*sin(theta1_rad)*sin(phi1_rad) + momentum2*sin(theta2_rad)*sin(phi2_rad) + momentum3*sin(theta3_rad)*sin(phi3_rad);
	H4 = momentum1*cos(theta1_rad) + momentum2*cos(theta2_rad) + momentum3*cos(theta3_rad);


	cv::Mat d = (cv::Mat_<double>(4, 1) << H1, H2, H3, H4);


	double E11 = m1 / sqrt(pow(m1, 2) + pow(momentum1, 2));
	double E21 = 0;
	double E31 = 0;
	double E41 = 0;

	cv::Mat E = (cv::Mat_<double>(4, 1) << E11,
		E21,
		E31,
		E41);  //unknown parameter 






	cv::Mat Dt = D.t();


	cv::Mat Vd0 = D*Va*Dt;
	cv::Mat Vd = Vd0.inv();



	cv::Mat E_t = E.t();
	cv::Mat VE0 = E_t*Vd*E;
	cv::Mat VE = VE0.inv();





	cv::Mat lambda0 = Vd*d;

	cv::Mat z = -VE*E_t*lambda0;

	cv::Mat lambda = lambda0 + Vd*E*z;


	cv::Mat lambda_t = lambda.t();

	cv::Mat Vd_i = Vd.inv();

	cv::Mat chi_square = lambda_t*Vd_i*lambda;



	cv::Mat a = a0 - Va*Dt*lambda;//caluclated value

	cv::Mat Vz = VE;

	cv::Mat V_lambda = Vd - Vd*E*VE*E_t*Vd;
	cv::Mat V = Va - Va*Dt*V_lambda*D*Va;






	double chi_square_value = chi_square.at<double>(0, 0);


	double new_momentum1 = a.at<double>(0, 0);
	double new_theta1_rad = a.at<double>(1, 0);
	double new_phi1_rad = a.at<double>(2, 0);
	double new_momentum2 = a.at<double>(3, 0);
	double new_theta2_rad = a.at<double>(4, 0);
	double new_phi2_rad = a.at<double>(5, 0);
	double new_momentum3 = a.at<double>(6, 0);
	double new_theta3_rad = a.at<double>(7, 0);
	double new_phi3_rad = a.at<double>(8, 0);
	double new_m_xi = a.at<double>(9, 0);


	double new_err_momentum1 = sqrt(V.at<double>(0, 0));
	double new_err_theta1_do = sqrt(V.at<double>(1, 1)) * 180 / pi;
	double new_err_phi1_do = sqrt(V.at<double>(2, 2)) * 180 / pi;
	double new_err_momentum2 = sqrt(V.at<double>(3, 3));
	double new_err_theta2_do = sqrt(V.at<double>(4, 4)) * 180 / pi;
	double new_err_phi2_do = sqrt(V.at<double>(5, 5)) * 180 / pi;
	double new_err_momentum3 = sqrt(V.at<double>(6, 6));
	double new_err_theta3_do = sqrt(V.at<double>(7, 7)) * 180 / pi;
	double new_err_phi3_do = sqrt(V.at<double>(8, 8)) * 180 / pi;
	double new_err_m_xi = sqrt(V.at<double>(9, 9));


	double new_m1 = initial_m1 + z.at<double>(0, 0);

	double new_err_m1 = sqrt(Vz.at<double>(0, 0));



	double new_H1, new_H2, new_H3, new_H4;//constraints 
	new_H1 = sqrt(pow(new_m1, 2) + pow(new_momentum1, 2)) + sqrt(pow(m2, 2) + pow(new_momentum2, 2)) + sqrt(pow(m3, 2) + pow(new_momentum3, 2)) - (m_Nuclear + new_m_xi - B_xi);
	new_H2 = new_momentum1*sin(new_theta1_rad)*cos(new_phi1_rad) + new_momentum2*sin(new_theta2_rad)*cos(new_phi2_rad) + new_momentum3*sin(new_theta3_rad)*cos(new_phi3_rad);
	new_H3 = new_momentum1*sin(new_theta1_rad)*sin(new_phi1_rad) + new_momentum2*sin(new_theta2_rad)*sin(new_phi2_rad) + new_momentum3*sin(new_theta3_rad)*sin(new_phi3_rad);
	new_H4 = new_momentum1*cos(new_theta1_rad) + new_momentum2*cos(new_theta2_rad) + new_momentum3*cos(new_theta3_rad);



	double delta_H1 = d.at<double>(0, 0) / sqrt(Vd_i.at<double>(0, 0));
	double delta_H2 = d.at<double>(1, 0) / sqrt(Vd_i.at<double>(1, 1));
	double delta_H3 = d.at<double>(2, 0) / sqrt(Vd_i.at<double>(2, 2));
	double delta_H4 = d.at<double>(3, 0) / sqrt(Vd_i.at<double>(3, 3));





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
	//	printf("m_D       = %.3f\n",m1);
	//	printf("BLL       = %.3f\n",1116.483*2+3727.380-m1);
	//	printf("m_D(new)  = %.3f Å} %.3f\n",new_m1,new_err_m1);
	//	printf("BLL(new)  = %.3f",1115.683*2+3727.380-new_m1);
	//
	//	getchar();



	return chi_square_value;
}