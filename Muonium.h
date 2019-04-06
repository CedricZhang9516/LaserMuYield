#ifndef __MUONIUM_HH__
#define __MUONIUM_HH__

#include <iostream>
#include <fstream>
#include <cmath>


class Muonium
{
 public:
  Muonium( std::array<double,3> pos, double t, std::array<double,3> vel );
  ~Muonium();

  std::array<double,3> GetPosition( double t, double& x, double& y, double& z );
  std::array<double,3> GetPosition( double t );
  std::array<double,3> GetVelocity( double t, double& vx, double& vy, double& vz );
  std::array<double,3> GetVelocity( double t );
  double GetAbsVelocity( double t );
  double GetStartTime( void ) { return fStartTime; };
  void SetTemperature( double temp /*[K]*/ );
  double GetTemperature( void ) { return fTemperature; };  // [K]
  
  static double GetOmega( double intensity /*[W/mm^2]*/) { return 4*M_PI * beta_ge * intensity; };             // [rad/s] (PRA73,052501,2006 Eq.6)
  static double GetIonizationRate( double intensity /*[W/mm^2]*/ ) { return 2*M_PI * beta_ioni * intensity; }; // [rad/s] (PRA73,052501,2006 Eq.9)
  static double Get2ndDoppler( double v/*[mm/s]*/ ) { return Omega_eg/2 * pow(v/2.99792458e11,2); };           // [rad/s] (PRA73,052501,2006 Eq.12)
  static double GetACStarkShift1S( double intensity /*[W/mm^2]*/) { return 2*M_PI * beta_ac_1S * intensity; }; // [rad/s] for 244 nm (tentatively used for also 355 nm)
  static double GetACStarkShift2S( double intensity /*[W/mm^2]*/) { return 2*M_PI * beta_ac_2S * intensity; }; // [rad/s] for 244 nm (tentatively used for also 355 nm)
  static double GetDCStarkShift( double efield /*[V/mm]*/ ) { return 2*M_PI * 4.7E5 * pow(efield,2); };        // [rad/s] (PLA187,247,1994)

  // constants
  static constexpr double ElectronMass = 0.510998928; // MeV/c^2 (PDG2014)
  static constexpr double MuonMass     = 105.6583715; // MeV/c^2 (PDG2014)
  static constexpr double Omega_eg = 2*M_PI * 2455529.002e9; // rad/s (PLA187,247,1994)
  static constexpr double Lambda_eg = 122.0887467E-6; // mm (PLA187,247,1994)

  static constexpr double GammaS     = 2*M_PI * 1.31; // rad/s (spontaneous decay, PRA73,052501,2006, Table1)
  static constexpr double GammaDecay = 2*M_PI * 72.442e3; // rad/s (FWHM)
  static constexpr double beta_ge    =  3.68111e1 * (1+ElectronMass/MuonMass) * (1+ElectronMass/MuonMass) * (1+ElectronMass/MuonMass); // [Hz/(W/mm^2)] (PRA73,052501,2006, Table2 & Eq.43)
  static constexpr double beta_ac_1S = -2.67827e1 * (1+ElectronMass/MuonMass) * (1+ElectronMass/MuonMass) * (1+ElectronMass/MuonMass); // [Hz/(W/mm^2)] for 244 nm (PRA73,052501,2006, Table4 & Eq.43)
  static constexpr double beta_ac_2S =  1.39927e2 * (1+ElectronMass/MuonMass) * (1+ElectronMass/MuonMass) * (1+ElectronMass/MuonMass); // [Hz/(W/mm^2)] for 244 nm (PRA73,052501,2006, Table4 & Eq.43)
  static constexpr double beta_ioni  =  1.20208e2 * (1+ElectronMass/MuonMass) * (1+ElectronMass/MuonMass) * (1+ElectronMass/MuonMass); // [Hz/(W/mm^2)] (PRA73,052501,2006, Table4 & Eq.43)

 private:
  static constexpr int nRhos = 5;
  std::array<double,3> fStartPos; // (X_sf, Y_sf, Z_sf)
  double fStartTime;   // ( t )
  std::array<double,3> fStartVel; // (VX_sf, XY_sf, VZ_sf)
  std::array<double,nRhos> fRho;  // (rho_gg, Re(rho_ge), Im(rho_ge), rho_ee, rho_ion)

  double fTemperature;
};



Muonium::Muonium( std::array<double,3> pos, double t, std::array<double,3> vel )
{
  fStartPos = pos;
  fStartVel = vel;
  fStartTime = t;
  for( auto& r: fRho ){ r = 0.0; }
  fRho[0] = 1.0;
  fTemperature = 322.; // (default value)
}

Muonium::~Muonium()
{}


void Muonium::SetTemperature( double temp )
{
  for( auto &v: fStartVel ){
    v *= sqrt( temp / fTemperature );
  }
  fTemperature = temp;
}


std::array<double,3> Muonium::GetPosition( double t, double& x, double& y, double& z )
{
  if( t < fStartTime ){
    // The Mu does not exist yet
    x = std::numeric_limits<double>::lowest();
    y = std::numeric_limits<double>::lowest();
    z = std::numeric_limits<double>::lowest();
  }else{
    x = fStartPos[0] + fStartVel[0]*(t-fStartTime);
    y = fStartPos[1] + fStartVel[1]*(t-fStartTime);
    z = fStartPos[2] + fStartVel[2]*(t-fStartTime);
  }
  return {x, y, z};

}

std::array<double,3> Muonium::GetPosition( double t )
{
  double x=0, y=0, z=0;
  return GetPosition( t, x, y, z );
}


std::array<double,3> Muonium::GetVelocity( double t, double& vx, double& vy, double& vz )
{
  vx = fStartVel[0];
  vy = fStartVel[1];
  vz = fStartVel[2];

  return {vx, vy, vz};
}

std::array<double,3> Muonium::GetVelocity( double t )
{
  return fStartVel;
}


double Muonium::GetAbsVelocity( double t )
{
  double v = 0;
  for( auto vi: fStartVel ){
    v += pow( vi, 2 );
  }
  return sqrt( v );
}



#endif //__MUONIUM_HH__
