// -*- C++ -*-

#ifndef DEBUG_TIMER_HH
#define DEBUG_TIMER_HH

#include <iostream>
#include <string>

//_____________________________________________________________________
namespace debug
{
  class Timer
  {
  public:
    explicit Timer( const std::string& msg="",
		    bool verbose=true );
    ~Timer( void );

  private:
    Timer( const Timer& );
    Timer& operator =( const Timer& );

  private:
    ::timespec* m_start;
    ::timespec* m_stop;
    std::string m_cstart;
    std::string m_cstop;
    std::string m_msg;
    bool        m_verbose;

  public:
    double sec( void ) const;
    double nsec( void ) const;
    void stop( void );
    void print( const std::string& arg="",
		std::ostream& ost=std::cout ) const;
    void print_verbose( const std::string& arg="",
			std::ostream& ost=std::cout ) const;
  };
}
#endif
