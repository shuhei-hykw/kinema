// -*- C++ -*-

#ifndef DATABASE_HH
#define DATABASE_HH

#include <string>
#include <vector>

class Particle;

//______________________________________________________________________________
class Database
{
public:
  static Database& GetInstance( void );
  ~Database( void );

private:
  Database( void );
  Database( const Database& );
  Database& operator =( const Database& );

private:
  std::vector<Particle*> m_database;

public:
  void                   Print( const std::string& arg="" ) const;
  std::vector<Particle*> Get( void ) const { return m_database; }
  std::size_t            Size( void ) const { return m_database.size(); }
  static std::string     ClassName( void ){ return "Database"; }

};

//______________________________________________________________________________
inline Database&
Database::GetInstance( void )
{
  static Database g_instance;
  return g_instance;
}

#endif
