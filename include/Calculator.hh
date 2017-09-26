// -*- C++ -*-

//______________________________________________________________________________
namespace calc
{

//______________________________________________________________________________
double
KineticEnergy( double mass, double range, int a, int z, int s);

//______________________________________________________________________________
bool
ThreeCharged( int S1, int S2, int S3,
	      double Range1, double Range2, double Range3,
	      double errRange1,double errRange2,double errRange3,
	      double theta1,double theta2,double theta3,
	      double errtheta1,double errtheta2,double errtheta3,
	      double phi1,double phi2,double phi3,
	      double errphi1,double errphi2,double errphi3 );
}
