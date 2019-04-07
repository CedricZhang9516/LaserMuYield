#ifndef __MUONIUM0_HH__
#define __MUONIUM0_HH__

#include <iostream>
#include <fstream>
#include <cmath>
#include "TMath.h"

double D = 87000;// diffussion coefficient mm^2/s
double T = 322;
double light = 299792458; // m/s
double massMu = 106.16/light/light; // MeV/c2
double lifeMu = 2.2e-6; //s
double k = 8.62e-11 ; //Boltzmann constant, MeV/K
double PI = 3.1415926;
double Thick = 8.8;
double vel0_avrg = 1000*sqrt(8*k*T/(PI*massMu));

class Muonium0
{
 public:
  Muonium0( double pos_x, double pos_y, double pos_z, double t );
  ~Muonium0();

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


  void Diffusion();

  double Get_decayT(){return decayT;};
  double Get_vel0(){return vel0;};
  double Get_theta0(){return theta0;};
  double Get_phi0(){return phi0;};
  double Get_Vx0(){return Vx0;};// = vel0 * sin(theta0) * cos(phi0);
  double Get_Vy0(){return Vy0;};// = vel0 * sin(theta0) * sin(phi0);
  double Get_Vz0(){return Vz0;};// = vel0 * cos(theta0);
  double Get_Lmfp(){return Lmfp;};// = 12*D/(PI * vel0_avrg);//mm

  double Get_X_sf(){return   X_sf;};
  double Get_Y_sf(){return   Y_sf;};
  double Get_Z_sf(){return   Z_sf;};
  double Get_VX_sf(){return  VX_sf;};
  double Get_VY_sf(){return  VY_sf;};
  double Get_VZ_sf(){return  VZ_sf;};
  double Get_T_sf(){return   T_sf;};
  double Get_theta_sf(){return   theta_sf;};
  double Get_phi_sf(){return     phi_sf;};

  bool Get_flag_sf(){return flag_sf;};

 private:
  static constexpr int nRhos = 5;
  std::array<double,3> fStartPos; // (X_sf, Y_sf, Z_sf)
  double fStartTime;   // ( t )
  std::array<double,3> fStartVel; // (VX_sf, XY_sf, VZ_sf)
  std::array<double,nRhos> fRho;  // (rho_gg, Re(rho_ge), Im(rho_ge), rho_ee, rho_ion)

  double fTemperature;

  double X_0 ;//= fMu.Get_X_sf();
  double Y_0 ;//= fMu.Get_Y_sf();
  double Z_0 ;//= fMu.Get_Z_sf();
  double T_0 ;//= fMu.Get_Z_sf();

  double vel0;
  double theta0, phi0;
  double Vx0;// = vel0 * sin(theta0) * cos(phi0);
  double Vy0;// = vel0 * sin(theta0) * sin(phi0);
  double Vz0;// = vel0 * cos(theta0);
  double Lmfp;// = 12*D/(PI * vel0_avrg);//mm
  double decayT;

  bool flag_sf;
  double X_sf ;//= fMu.Get_X_sf();
  double Y_sf ;//= fMu.Get_Y_sf();
  double Z_sf ;//= fMu.Get_Z_sf();
  double VX_sf ;//= fMu.Get_VX_sf();
  double VY_sf ;//= fMu.Get_VY_sf();
  double VZ_sf ;//= fMu.Get_VZ_sf();
  double T_sf ;//= fMu.Get_T_sf();
  double theta_sf ;//= fMu.Get_theta_sf();
  double phi_sf ;//= fMu.Get_phi_sf();
  

};



Muonium0::Muonium0( double pos_x, double pos_y, double pos_z, double t )
{
  X_0 = pos_x;
  Y_0 = pos_y;
  Z_0 = pos_z;
  T_0 = t;
  flag_sf = 0;

}

Muonium0::~Muonium0()
{}


void Muonium0::Diffusion(){

  double tempX, tempY;

  double t = 0;
  double z = Z_0;
  double x = X_0;
  double y = Y_0;
  double vx = 0;// = Vx0;
  double vy = 0;// = Vy0;
  double vz = 0;// = Vz0;

  double L = 0;
  double theta, phi;

  do
  {
    tempX = ((double) rand() / (RAND_MAX)) * 3 * sqrt( 2 * k * T/massMu );
    tempY = ((double) rand() / (RAND_MAX)) * 2 * sqrt( 2 * massMu/PI/k/T) * exp(-1.0);
  }while(tempY>sqrt(2/PI*pow(massMu/k/T,3))*pow(tempX,2)*exp(-massMu*tempX*tempX/(2*k*T)));
  
  vel0 = tempX*1000 ; // mm/s

  // generate time of Mu decay event
  do {
    tempX = ((double) rand() / (RAND_MAX)) * 10e-5;  // t=0~10^-5 s
    tempY = ((double) rand() / (RAND_MAX)) * 1/lifeMu;
  } while (tempY>1/lifeMu*exp(-tempX/lifeMu) );

  decayT = tempX;
  
  // Random walk at begining
  // Generate theta and phi
  tempX = TMath::ACos(-1 + 2 * ((double)rand()/(RAND_MAX)));
  tempY = ((double) rand()/(RAND_MAX)) * 2* PI;

  theta0 = tempX;  // in rad unit: [0,PI]
  phi0 = tempY; //[0,+2PI] radians
  
  // initial velocity input in MCmfp random walk
  Vx0 = vel0 * sin(theta0) * cos(phi0);
  Vy0 = vel0 * sin(theta0) * sin(phi0);
  Vz0 = vel0 * cos(theta0);

  vx = Vx0;
  vy = Vy0;
  vz = Vz0;

  Lmfp = 12*D/(PI * vel0_avrg);//mm
  //hLmfp->Fill(Lmfp * 1000);//um

  do{

    //if(flag_DiffusionTrack==1){DiffusionTrack->SetPoint(N,z,y);}//cout<<"filling track"<<N<<endl;}

    // generate the L step lenth according to Lmfp exp. distri.
    do {
      tempX = ((double) rand() / (RAND_MAX)) * 5;  // l=0~0.01 mm
      tempY = ((double) rand() / (RAND_MAX)) * 1/Lmfp;
    } while (tempY>1/Lmfp*exp(-tempX/Lmfp) );

    L = tempX;

    t = t + L/vel0;
    x = x + vx * (L/vel0);
    y = y + vy * (L/vel0);
    z = z + vz * (L/vel0);

    if( (z <= 0 && z >= -Thick) )//(z <= 14.12 && z >= 7) )//(z <= -9 && z >= -11) || (z <= 9 && z >= 7) )
    {

      tempX = TMath::ACos(-1 + 2 * ((double)rand()/(RAND_MAX)) );
      tempY = ((double) rand()/(RAND_MAX)) * 2* PI;

      theta = tempX;  // [0,+PI] radians
      phi = tempY; //[0,+2PI] radians

      vx = vel0 * sin(theta) * cos(phi);
      vy = vel0 * sin(theta) * sin(phi);
      vz = vel0 * cos(theta);
    }

    if (t>decayT) t = -10e-6;


  }while( (z <= 0 && z >= -Thick && t<decayT));

  if(z>=0) flag_sf = 1;
  
  X_sf = x;//= fMu.Get_X_sf();
  Y_sf = y;//= fMu.Get_Y_sf();
  Z_sf = z;//= fMu.Get_Z_sf();
  VX_sf = vx ;//= fMu.Get_VX_sf();
  VY_sf = vy;//= fMu.Get_VY_sf();
  VZ_sf = vz;//= fMu.Get_VZ_sf();
  T_sf = T_0 + t;//= fMu.Get_T_sf();
  theta_sf = theta;//= fMu.Get_theta_sf();
  phi_sf = phi;//= fMu.Get_phi_sf();
  

}

void Muonium0::SetTemperature( double temp )
{
  for( auto &v: fStartVel ){
    v *= sqrt( temp / fTemperature );
  }
  fTemperature = temp;
}


std::array<double,3> Muonium0::GetPosition( double t, double& x, double& y, double& z )
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

std::array<double,3> Muonium0::GetPosition( double t )
{
  double x=0, y=0, z=0;
  return GetPosition( t, x, y, z );
}


std::array<double,3> Muonium0::GetVelocity( double t, double& vx, double& vy, double& vz )
{
  vx = fStartVel[0];
  vy = fStartVel[1];
  vz = fStartVel[2];

  return {vx, vy, vz};
}

std::array<double,3> Muonium0::GetVelocity( double t )
{
  return fStartVel;
}


double Muonium0::GetAbsVelocity( double t )
{
  double v = 0;
  for( auto vi: fStartVel ){
    v += pow( vi, 2 );
  }
  return sqrt( v );
}



#endif //__MUONIUM0_HH__
