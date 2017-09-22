// -*- C++ -*-

//���̃v���O�����́A�������ʂ�蓾��ꂽRange�Atheta�Aphi���A���݂��m�F����Ă��邷�ׂĂ̗��q�ɂ�����track�𓖂Ă͂߂āA�e���q�̎��ʂ��v�Z����B
//���ꂪ�덷�̂R�{�ȓ��Ő��藧���̂������o�����̂ł���B�\������ۂ́A�ǂ̗��q�𓖂Ă͂߂āA�v�Z���������AZ,A,S��p���āA�\�����Ă��܂��B
//�o�͂���result.txt�ɂ́A��������\�L���āA�ǂ̂悤�Ȕ��������肳�ꂽ�������₷�����Ă���܂��B

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
