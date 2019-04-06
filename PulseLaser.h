#ifndef __PULSELASER_HH__
#define __PULSELASER_HH__



#include <iostream>
#include <fstream>
#include <cmath>



class PulseLaser
{

 public:
  PulseLaser( double x, double y, double z, double t, double d, double w, double e );
  void   SetWaist0( double w ){ fWaist0 = w; };
  double GetWaist0( void ) { return fWaist0; };
  void   SetEnergy( double e ){ fEnergy = e; };
  double GetEnergy( void ) { return fEnergy; };

  std::array<double,3> GetCenterPosition( void ) { return fCenterPos; };
  std::array<double,3> GetCenterPosition( double &x, double &y, double &z );

  double GetPeakPower( void ){ return fEnergy / sqrt(2*M_PI) / fTSigma; };
  double GetPeakTime( double x ){ return (x-fCenterPos[0])/C + fT0; };
  double GetIntensity( double x, double y, double z, double t ); // [W/mm^2] = [MW/m^2]
  double GetIntensity( std::array<double,3> pos, double t ) { return GetIntensity( pos[0], pos[1], pos[2], t); };


  double GetRayleighLength( void ) { return M_PI * pow(fWaist0,2) / fWavelength; };
  double GetWaist( double x );

  static double FWHM2Sigma( double fwhm ){ return fwhm / (2.*sqrt(2*log(2))); };
  static double Sigma2FWHM( double sigma ){ return sigma * (2.*sqrt(2*log(2))); };

  void DumpSetting( void );

 private:
  static constexpr double fWavelength = 355E-6; // [mm]
  static constexpr double C = 299792458E3; // light velocity [mm/s]
  std::array<double,3> fCenterPos; // (x,y,z) [mm]
  double fT0; // [s] Pulse exist at (fCenterX, fCenterY, fCenterZ)
  double fWaist0; // 1/e^2 radius [mm]
  double fDuration; // FWHM time duration [s]
  double fTSigma; // time duration in sigma [s] 
  double fEnergy; // J


};



PulseLaser::PulseLaser( double x, double y, double z, double t, double d, double w, double e )
{
  fCenterPos[0] = x;
  fCenterPos[1] = y;
  fCenterPos[2] = z;
  fT0 = t;
  fDuration = d;
  fTSigma = FWHM2Sigma( d );
  fWaist0 = w;
  fEnergy = e;
}


std::array<double,3> PulseLaser::GetCenterPosition( double &x, double &y, double &z )
{
  x = fCenterPos[0];
  y = fCenterPos[1];
  z = fCenterPos[2];

  return fCenterPos;
}

double PulseLaser::GetWaist( double x ){
  // Yariv (2.5-14)
  return fWaist0 * sqrt( 1 + pow( (x-fCenterPos[0])/GetRayleighLength(), 2) );
}


double PulseLaser::GetIntensity( double x, double y, double z, double t )
{
  double powerX = GetPeakPower() * exp( - pow( (t-GetPeakTime(x))/fTSigma, 2)/2 );

  // Demtroder, Laser spectroscopy Vol. 1 (5.147)
  double r = sqrt( pow(y-fCenterPos[1],2) + pow(z-fCenterPos[2],2) );
  double w = GetWaist( x );

  return 2*powerX / M_PI / pow(w,2) * exp( -2*pow(r/w,2) );
}


void PulseLaser::DumpSetting( void )
{
  std::cout << "Beam center : " << fCenterPos[0] << " " << fCenterPos[1] << " " << fCenterPos[2] << " mm" << std::endl;
  std::cout << "Beam center timing : " << fT0 << " s" << std::endl;
  std::cout << "Beam energy  : " << fEnergy << " J" << std::endl;
  std::cout << "Beam waist  : " << fWaist0 << " mm (radius at the intensity of 1/e^2)" << std::endl;
  std::cout << "Beam duration : " << fDuration << " s (FWHM), " << fTSigma << " s (sigma)" << std::endl;
  std::cout << "RayleighLength : " << GetRayleighLength() << " mm (M^2==1)" << std::endl;
  std::cout << "Peak power : " << GetPeakPower() << " W" << std::endl;
  std::cout << "Peak intensity at the waist : " << GetIntensity( fCenterPos, fT0 ) << " W/mm^2" << std::endl;
}

#endif //__PULSELASER_HH__
