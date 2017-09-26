// -*- C++ -*-

#include <libgen.h>
#include <iostream>

#include <TTimeStamp.h>

#include "Database.hh"
#include "DebugTimer.hh"
#include "UserParamMan.hh"
#include "KinematicsAnalyzer.hh"

namespace
{
  enum EArg { kProcess, kParamFile, kArgc };
  UserParamMan&       gParam  = UserParamMan::GetInstance();
  KinematicsAnalyzer& gKinema = KinematicsAnalyzer::GetInstance();
  Database&           gData   = Database::GetInstance();
}

//______________________________________________________________________________
int
main( int argc, char *argv[] )
{
  debug::Timer s;

  if( argc != kArgc ){
    std::cout << std::endl
	      << "Usage : " << basename(argv[kProcess])
	      << " [ParamFile]" << std::endl
	      << std::endl;
    return EXIT_SUCCESS;
  }

  if( !gParam.Initialize( argv[kParamFile] ) )
    return EXIT_FAILURE;

  // gData.Print();

  gKinema.Print("Before");
  gKinema.DoSearch();
  gKinema.Print("After");
}
