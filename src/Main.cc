// -*- C++ -*-

//このプログラムは、実験結果より得られたRange、theta、phiより、存在が確認されているすべての粒子においてtrackを当てはめて、親粒子の質量を計算する。
//それが誤差の３倍以内で成り立つものを見つけ出すものである。表示する際は、どの粒子を当てはめて、計算したかを、Z,A,Sを用いて、表示しています。
//出力したresult.txtには、反応式を表記して、どのような反応が仮定されたかを見やすくしてあります。

#include <libgen.h>
#include <iostream>
#include <vector>

#include "Database.hh"
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
  if( argc != kArgc ){
    std::cout << std::endl
	      << "Usage : " << basename(argv[kProcess])
	      << " [ParamFile]" << std::endl
	      << std::endl;
    return EXIT_SUCCESS;
  }

  if( !gParam.Initialize( argv[kParamFile] ) )
    return EXIT_FAILURE;

  gData.Print();

  gKinema.Print("Before");
  gKinema.DoSearch();
  gKinema.Print("After");
}
