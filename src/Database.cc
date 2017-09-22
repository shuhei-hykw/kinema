// -*- C++ -*-

#include "Database.hh"

#include <iostream>
#include <string>
#include <vector>

#include "Particle.hh"

//______________________________________________________________________________
Database::Database( void )
{
  //                        Particle( Name,      A,  Z,  S,  Mass[MeV/c2] )
  m_database.push_back( new Particle( "n",       1,  0,  0,   939.565 ) );
  m_database.push_back( new Particle( "p",       1,  1,  0,   938.272 ) );
  m_database.push_back( new Particle( "d",       2,  1,  0,  1875.613 ) );
  m_database.push_back( new Particle( "t",       3,  1,  0,  2808.922 ) );
  m_database.push_back( new Particle( "2n",      2,  0,  0,  1879.128 ) );
  m_database.push_back( new Particle( "He3",     3,  2,  0,  2808.392 ) );
  m_database.push_back( new Particle( "He4",     4,  2,  0,  3727.380 ) );
  m_database.push_back( new Particle( "He5",     5,  2,  0,  4667.845 ) );
  m_database.push_back( new Particle( "Li5",     5,  3,  0,  4667.624 ) );
  m_database.push_back( new Particle( "He6",     6,  2,  0,  5605.541 ) );
  m_database.push_back( new Particle( "Li6",     6,  3,  0,  5601.520 ) );
  m_database.push_back( new Particle( "Be6",     6,  4,  0,  5605.298 ) );
  m_database.push_back( new Particle( "He7",     7,  2,  0,  6545.550 ) );
  m_database.push_back( new Particle( "Li7",     7,  3,  0,  6533.836 ) );
  m_database.push_back( new Particle( "Be7",     7,  4,  0,  6534.186 ) );
  m_database.push_back( new Particle( "He8",     8,  2,  0,  7482.542 ) );
  m_database.push_back( new Particle( "Li8",     8,  3,  0,  7471.369 ) );
  m_database.push_back( new Particle( "Be8",     8,  4,  0,  7454.852 ) );
  m_database.push_back( new Particle( "B8",      8,  5,  0,  7472.321 ) );
  m_database.push_back( new Particle( "Li9",     9,  3,  0,  8406.871 ) );
  m_database.push_back( new Particle( "Be9",     9,  4,  0,  8392.753 ) );
  m_database.push_back( new Particle( "B9",      9,  5,  0,  8393.310 ) );
  m_database.push_back( new Particle( "C9",      9,  6,  0,  8409.295 ) );
  m_database.push_back( new Particle( "Be10",   10,  4,  0,  9325.507 ) );
  m_database.push_back( new Particle( "B10",    10,  5,  0,  9324.440 ) );
  m_database.push_back( new Particle( "C10",    10,  6,  0,  9327.580 ) );
  m_database.push_back( new Particle( "Li11",   11,  3,  0, 10285.846 ) );
  m_database.push_back( new Particle( "Be11",   11,  4,  0, 10264.570 ) );
  m_database.push_back( new Particle( "B11",    11,  5,  0, 10252.550 ) );
  m_database.push_back( new Particle( "C11",    11,  6,  0, 10254.022 ) );
  m_database.push_back( new Particle( "Be12",   12,  4,  0, 11200.922 ) );
  m_database.push_back( new Particle( "B12",    12,  5,  0, 11188.746 ) );
  m_database.push_back( new Particle( "C12",    12,  6,  0, 11174.866 ) );
  m_database.push_back( new Particle( "N12",    12,  7,  0, 11191.693 ) );
  m_database.push_back( new Particle( "Be13",   13,  4,  0, 12142.542 ) );
  m_database.push_back( new Particle( "B13",    13,  5,  0, 12123.434 ) );
  m_database.push_back( new Particle( "C13",    13,  6,  0, 12109.485 ) );
  m_database.push_back( new Particle( "N13",    13,  7,  0, 12111.195 ) );
  m_database.push_back( new Particle( "O13",    13,  8,  0, 12128.443 ) );
  m_database.push_back( new Particle( "B14",    14,  5,  0, 13062.023 ) );
  m_database.push_back( new Particle( "C14",    14,  6,  0, 13040.874 ) );
  m_database.push_back( new Particle( "N14",    14,  7,  0, 13040.207 ) );
  m_database.push_back( new Particle( "O14",    14,  8,  0, 13044.841 ) );
  m_database.push_back( new Particle( "C15",    15,  6,  0, 13979.222 ) );
  m_database.push_back( new Particle( "N15",    15,  7,  0, 13968.939 ) );
  m_database.push_back( new Particle( "O15",    15,  8,  0, 13971.182 ) );
  m_database.push_back( new Particle( "C16",    16,  6,  0, 14914.536 ) );
  m_database.push_back( new Particle( "N16",    16,  7,  0, 14906.014 ) );
  m_database.push_back( new Particle( "O16",    16,  8,  0, 14895.084 ) );
  m_database.push_back( new Particle( "O17",    17,  8,  0, 15830.506 ) );
  m_database.push_back( new Particle( "pi-",     0, -1,  0,   139.570 ) );
  m_database.push_back( new Particle( "pi0",     0,  0,  0,   134.977 ) );
  m_database.push_back( new Particle( "L",       1,  0, -1,  1115.683 ) );
  m_database.push_back( new Particle( "H3L",     3,  1, -1,  2991.166 ) );
  m_database.push_back( new Particle( "H4L",     4,  1, -1,  3922.565 ) );
  m_database.push_back( new Particle( "He4L",    4,  2, -1,  3921.685 ) );
  m_database.push_back( new Particle( "He5L",    5,  2, -1,  4839.943 ) );
  m_database.push_back( new Particle( "He6L",    6,  2, -1,  5779.348 ) );
  m_database.push_back( new Particle( "Li6L",    6,  3, -1,  5778.807 ) );
  m_database.push_back( new Particle( "Li7L",    7,  3, -1,  6711.623 ) );
  m_database.push_back( new Particle( "Be7L",    7,  4, -1,  6715.821 ) );
  m_database.push_back( new Particle( "He8L",    8,  2, -1,  7654.073 ) );
  m_database.push_back( new Particle( "Li8L",    8,  3, -1,  7642.719 ) );
  m_database.push_back( new Particle( "Be8L",    8,  4, -1,  7643.029 ) );
  m_database.push_back( new Particle( "Li9L",    9,  3, -1,  8578.552 ) );
  m_database.push_back( new Particle( "Be9L",    9,  4, -1,  8563.825 ) );
  m_database.push_back( new Particle( "B9L",     9,  5, -1,  8579.714 ) );
  m_database.push_back( new Particle( "Be10L",  10,  4, -1,  9499.326 ) );
  m_database.push_back( new Particle( "B10L",   10,  5, -1,  9500.103 ) );
  m_database.push_back( new Particle( "B11L",   11,  5, -1, 10429.883 ) );
  m_database.push_back( new Particle( "B12L",   12,  5, -1, 11356.863 ) );
  m_database.push_back( new Particle( "C12L",   12,  6, -1, 11358.905 ) );
  m_database.push_back( new Particle( "B13L",   13,  5, -1, 12293.059 ) );
  m_database.push_back( new Particle( "C13L",   13,  6, -1, 12278.859 ) );
  m_database.push_back( new Particle( "C14L",   14,  6, -1, 13212.998 ) );
  m_database.push_back( new Particle( "N14L",   14,  7, -1, 13214.708 ) );
  m_database.push_back( new Particle( "N15L",   15,  7, -1, 14142.300 ) );
  m_database.push_back( new Particle( "O16L",   16,  8, -1, 15074.365 ) );
  m_database.push_back( new Particle( "O18L",   18,  8, -1, 16931.689 ) );
  m_database.push_back( new Particle( "H4LL",    4,  1, -2,  4106.719 ) );
  m_database.push_back( new Particle( "H5LL",    5,  1, -2,  5036.208 ) );
  m_database.push_back( new Particle( "He5LL",   5,  2, -2,  5034.978 ) );
  m_database.push_back( new Particle( "He6LL",   6,  2, -2,  5952.506 ) );
  m_database.push_back( new Particle( "He7LL",   7,  2, -2,  6890.851 ) );
  m_database.push_back( new Particle( "Li7LL",   7,  3, -2,  6889.990 ) );
  m_database.push_back( new Particle( "Li8LL",   8,  3, -2,  7821.726 ) );
  m_database.push_back( new Particle( "Be8LL",   8,  4, -2,  7826.344 ) );
  m_database.push_back( new Particle( "He9LL",   9,  2, -2,  8762.596 ) );
  m_database.push_back( new Particle( "Li9LL",   9,  3, -2,  8751.602 ) );
  m_database.push_back( new Particle( "Be9LL",   9,  4, -2,  8751.872 ) );
  m_database.push_back( new Particle( "Li10LL", 10,  3, -2,  9685.735 ) );
  m_database.push_back( new Particle( "Be10LL", 10,  4, -2,  9672.798 ) );
  m_database.push_back( new Particle( "B10LL",  10,  5, -2,  9687.107 ) );
  m_database.push_back( new Particle( "Be11LL", 11,  4, -2, 10605.899 ) );
  m_database.push_back( new Particle( "B11LL",  11,  5, -2, 10606.896 ) );
  m_database.push_back( new Particle( "Be12LL", 12,  4, -2, 11538.653 ) );
  m_database.push_back( new Particle( "B12LL",  12,  5, -2, 11535.326 ) );
  m_database.push_back( new Particle( "B13LL",  13,  5, -2, 12461.176 ) );
  m_database.push_back( new Particle( "C13LL",  13,  6, -2, 12463.788 ) );
  m_database.push_back( new Particle( "B14LL",  14,  5, -2, 13397.372 ) );
  m_database.push_back( new Particle( "C14LL",  14,  6, -2, 13382.852 ) );
  m_database.push_back( new Particle( "C15LL",  15,  6, -2, 14316.511 ) );
  m_database.push_back( new Particle( "N15LL",  15,  7, -2, 14318.221 ) );
  m_database.push_back( new Particle( "C16LL",  16,  6, -2, 15247.900 ) );
  m_database.push_back( new Particle( "N16LL",  16,  7, -2, 15244.393 ) );
  m_database.push_back( new Particle( "O17LL",  17,  8, -2, 16177.548 ) );
  m_database.push_back( new Particle( "O19LL",  19,  8, -2, 18032.872 ) );
  m_database.push_back( new Particle( "Xi-C12", 13,  5, -2, 12496.576 ) );
  m_database.push_back( new Particle( "Xi-N14", 15,  6, -2, 14361.917 ) );
  m_database.push_back( new Particle( "Xi-O16", 17,  7, -2, 16216.794 ) );

}

//______________________________________________________________________________
Database::~Database( void )
{
  for( auto& p : m_database ){
    delete p; p = 0;
  }
}

//______________________________________________________________________________
void
Database::Print( const std::string& arg ) const
{
  static const std::string func_name("["+ClassName()+"::"+__func__+"()]");

  std::cout << "#D " << func_name << " " << arg << std::endl
	    << " NParticle = " << m_database.size() << std::endl;
  for( const auto& p : m_database ){
    p->Print();
  }
}
