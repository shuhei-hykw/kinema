// -*- C++ -*-

#ifndef USER_PARAM_MAN_HH
#define USER_PARAM_MAN_HH

#include <string>
#include <map>
#include <vector>

//______________________________________________________________________________
class UserParamMan
{
public:
  static UserParamMan&  GetInstance( void );
  ~UserParamMan( void );

private:
  UserParamMan( void );
  UserParamMan( const UserParamMan&  );
  UserParamMan& operator =( const UserParamMan& );

private:
  typedef std::vector<double>               ParamArray;
  typedef std::map<std::string, ParamArray> ParamMap;
  typedef ParamMap::const_iterator          PIterator;
  bool                     m_is_ready;
  std::string              m_file_name;
  ParamMap                 m_param_map;
  std::vector<std::string> m_key_list;

public:
  enum ETrackParam { S, Range, RangeE, Theta, ThetaE, Phi, PhiE, NTrackParam };

public:
  bool               Initialize( void );
  bool               Initialize( const std::string& filename );
  bool               IsReady( void ) const { return m_is_ready; }
  int                GetSize( const std::string& key ) const;
  double             GetParameter( const std::string& key, int i=0 ) const;
  void               Print( const std::string& arg="" ) const;
  static std::string ClassName( void ){ return "UserParamMan"; }
};

//______________________________________________________________________________
inline UserParamMan&
UserParamMan::GetInstance( void )
{
  static UserParamMan g_instance;
  return g_instance;
}

#endif
