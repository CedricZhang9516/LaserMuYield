#ifndef __CWLASER_HH__
#define __CWLASER_HH__


#include <iostream>
#include <fstream>
#include <cmath>



class CWLaser
{

 public:
  CWLaser( double x, double y, double z, double w, double p, double onTime, double offTime );
  CWLaser( double x, double y, double z, double w, double p );
  void   SetWaist0( double w ){ fWaist0 = w; };
  double GetWaist0( void ) { return fWaist0; };
  void   SetPower( double p ){ fPower = p; };
  double GetPower( void ) { return fPower; };
  void   SetDetuning( double d /*[rad/s]*/ ){ fDetuning = d; };
  double GetDetuning( void ) { return fDetuning; };
  double GetWavelength( void ) { return fWavelength; };
  double GetOnTime( void ) { return fONTime; };
  double GetOffTime( void ) { return fOFFTime; };
  std::array<double,3> GetCenterPosition( void ) { return fCenterPos; };
  std::array<double,3> GetCenterPosition( double &x, double &y, double &z );

  double GetIntensity( double x, double y, double z ); // [W/mm^2] = [MW/m^2]
  double GetIntensity( std::array<double,3> pos ) { return GetIntensity( pos[0], pos[1], pos[2] ); };
  double GetIntensity( double x, double y, double z, double t ); // [W/mm^2] = [MW/m^2] (On/OFF is considered)
  double GetIntensity( std::array<double,3> pos, double t ) { return GetIntensity( pos[0], pos[1], pos[2], t ); };
  double GetPeakIntensity( double x ) { return GetIntensity( x, fCenterPos[1], fCenterPos[2] ); };
  double GetPeakIntensity( void ) { return GetIntensity( fCenterPos ); };



  double GetRayleighLength( void ) { return M_PI * pow(fWaist0,2) / fWavelength; };
  double GetWaist( double x );
  void DumpSetting( void );

 private:
  static constexpr double fWavelength = 244.17749E-6; // [mm] (PRA187, 1994, p.247)
  std::array<double,3> fCenterPos; // (x,y,z) [mm]
  double fWaist0; // 1/e^2 radius [mm] (D2sigma)
  double fPower; // W (one direction)
  double fONTime; // [s]
  double fOFFTime; // [s]

  double fDetuning; // [rad/s]
};



CWLaser::CWLaser( double x, double y, double z, double w, double p, double onTime, double offTime )
{
  fCenterPos[0] = x;
  fCenterPos[1] = y;
  fCenterPos[2] = z;
  fWaist0 = w;
  fPower = p;
  fONTime = onTime;
  fOFFTime = offTime;
  fDetuning = 0;
}

CWLaser::CWLaser( double x, double y, double z, double w, double p )
{
  fCenterPos[0] = x;
  fCenterPos[1] = y;
  fCenterPos[2] = z;
  fWaist0 = w;
  fPower = p;
  fONTime =  std::numeric_limits<double>::lowest();
  fOFFTime = std::numeric_limits<double>::max();
  fDetuning = 0;
}


std::array<double,3> CWLaser::GetCenterPosition( double &x, double &y, double &z )
{
  x = fCenterPos[0];
  y = fCenterPos[1];
  z = fCenterPos[2];

  return fCenterPos;
}


double CWLaser::GetWaist( double x ){
  // Yariv (2.5-14)
  return fWaist0 * sqrt( 1 + pow( (x-fCenterPos[0])/GetRayleighLength(), 2) );
}


double CWLaser::GetIntensity( double x, double y, double z )
{
  // ON/OFF is not considered.
  // Demtroder, Laser spectroscopy Vol. 1 (5.147)
  double r = sqrt( pow(y-fCenterPos[1],2) + pow(z-fCenterPos[2],2) );
  double w = GetWaist( x );
  return 2*fPower / M_PI / pow(w,2) * exp( -2*pow(r/w,2) ); 
}


double CWLaser::GetIntensity( double x, double y, double z, double t )
{
  // ON/OFF is considered.
  if( t < fONTime ) return 0;
  if( fOFFTime < t ) return 0;
  return GetIntensity( x, y, z);
}


void CWLaser::DumpSetting( void )
{
  std::cout << "Beam center : " << fCenterPos[0] << " " << fCenterPos[1] << " " << fCenterPos[2] << " mm" << std::endl;
  std::cout << "Beam power  : " << fPower << " W" << std::endl;
  std::cout << "Beam waist  : " << fWaist0 << " mm (radius at the intensity of 1/e^2)" << std::endl;
  std::cout << "RayleighLength : " << GetRayleighLength() << " mm (M^2==1)" << std::endl;
  std::cout << "Peak intensity : " << GetPeakIntensity() << " W/mm^2" << std::endl;
  std::cout << "Time duration : " << fONTime << " - " << fOFFTime << " s" << std::endl;
  std::cout << "Detuning : " << fDetuning << " [rad/s] (" << fDetuning/(2*M_PI) << " [Hz])" << std::endl;
}


#endif //__CWLASER_HH__
