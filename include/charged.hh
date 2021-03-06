#ifndef CHARGED_H
#define CHARGED_H

//______________________________________________________________________________
double
ThreeCharged( double Range1,  double Range2,  double Range3,
	      double RangeE1, double RangeE2, double RangeE3,
	      double theta1,  double theta2,  double theta3,
	      double errtheta1,double errtheta2,double errtheta3,
	      double phi1,double phi2,double phi3,
	      double errphi1,double errphi2,double errphi3);

//______________________________________________________________________________
double
two_charged( double Range1, double Range2,
	     double errRange1, double errRange2,
	     double theta1, double theta2,
	     double errtheta1,double errtheta2,
	     double phi1,double phi2,
	     double errphi1,double errphi2 );

#endif
