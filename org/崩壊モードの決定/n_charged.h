#ifndef N_CHARGED_H
#define N_CHARGED_H

double n_three(double Range1,double Range2,double Range3,double errRange1,double errRange2,double errRange3,
					 double theta1,double theta2,double theta3,double errtheta1,double errtheta2,double errtheta3,
					 double phi1,double phi2,double phi3,double errphi1,double errphi2,double errphi3);


double n_two(double Range1,double Range2,double errRange1,double errRange2,double theta1,double theta2,double errtheta1,double errtheta2,
			double phi1,double phi2,double errphi1,double errphi2);

double n_one(double Range1,double errRange1,double theta1,double errtheta1,double phi1,double errphi1);


#endif
