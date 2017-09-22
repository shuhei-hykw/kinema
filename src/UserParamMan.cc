// -*- C++ -*-

#include "UserParamMan.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include "KinematicsAnalyzer.hh"
#include "RootHelper.hh"

namespace
{
  const double default_value = -9999.;
  KinematicsAnalyzer& gKinema = KinematicsAnalyzer::GetInstance();
}

// if no parameter,
//   0: throw exception
//   1: return default value
#define ReturnDefaultValue 0

//______________________________________________________________________________
UserParamMan::UserParamMan( void )
  : m_is_ready(false), m_file_name("")
{
}

//______________________________________________________________________________
UserParamMan::~UserParamMan( void )
{
}

//______________________________________________________________________________
Bool_t
UserParamMan::Initialize( void )
{
  static const std::string func_name("["+ClassName()+"::"+__func__+"()]");

  if( m_is_ready ){
    std::cerr << "#W " << func_name << " "
	      << "Already initialized" << std::endl;
    return false;
  }

  std::ifstream ifs( m_file_name );
  if( !ifs.is_open() ){
    std::cerr << "#E " << func_name << " "
	      << "No such parameter file : " << m_file_name << std::endl;
    return false;
  }

  std::string line;
  while( ifs.good() && std::getline(ifs,line) ){
    if( line.empty() || line[0]=='#' ) continue;

    std::istringstream input_line( line );

    std::string first_param;
    input_line >> first_param;

    std::string key = first_param;
    ParamArray param_array;
    double   param;
    while( input_line >> param ){
      param_array.push_back( TMath::Abs(param) );
    }

    if( param_array.size() == NTrackParam ){
      gKinema.Add( new Particle( key, param_array[S],
				 param_array[Range], param_array[RangeE],
				 param_array[Theta], param_array[ThetaE],
				 param_array[Phi], param_array[PhiE] ) );
    }

    m_param_map[key] = param_array;
    m_key_list.push_back( key );
  }

  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
Bool_t
UserParamMan::Initialize( const std::string& filename )
{
  m_file_name = filename;
  return Initialize();
};

//______________________________________________________________________________
int
UserParamMan::GetSize( const std::string& key ) const
{
  static const std::string func_name("["+ClassName()+"::"+__func__+"()]");

  PIterator itr = m_param_map.find(key);
  if( itr==m_param_map.end() ){
    Print(m_file_name);
    std::cerr << "#E " << func_name << " "
	      << "No such key : " << key << std::endl;
    return 0;
  }

  return itr->second.size();
}

//______________________________________________________________________________
double
UserParamMan::GetParameter( const std::string& key, int i ) const
{
  static const std::string func_name("["+ClassName()+"::"+__func__+"()]");

  std::stringstream param;
  param << key << "(" << i << ")";

  PIterator itr = m_param_map.find(key);

  if( itr==m_param_map.end() ){
    Print( m_file_name );
#if ReturnDefaultValue
    std::cerr << "#E " << func_name
	      << "set default value : " << param.str() << " -> "
	      << default_value << std::endl;
    return default_value;
#else
    throw std::out_of_range( func_name+" No such key : "+key );
#endif
  }

  if( i+1 > (int)itr->second.size() ){
    Print( m_file_name );
#if ReturnDefaultValue
    std::cerr << "#E " << func_name
	      << "set default value : " << param.str() << " -> "
	      << default_value << std::endl;
    return default_value;
#else
    throw std::out_of_range( func_name+" No such key : "+key );
#endif
  }

  return itr->second.at(i);
}

//______________________________________________________________________________
void
UserParamMan::Print( const std::string& arg ) const
{
  static const std::string func_name("["+ClassName()+"::"+__func__+"()]");

  std::cout << "#D " << func_name << " " << arg << std::endl;

  const int w = 20;
  PIterator itr, end=m_param_map.end();
  for( itr=m_param_map.begin(); itr!=end; ++itr ){
    std::cout << " key = " << std::setw(w) << std::left
	      << itr->first << itr->second.size() << " : ";
    for( int i=0, n=itr->second.size(); i<n; ++i ){
      std::cout << std::setw(5) << std::right
		<< itr->second.at(i) << " ";
    }
    std::cout << std::endl;
  }
}
