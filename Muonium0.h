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



  void Diffusion();
  std::array<double,3> Emission();
  std::array<double,3> EmissionPosition(double);
  bool InsideLaserRegion_cw(double, double, double);
  bool InsideLaserRegion_pulse(double, double, double);

  double Get_decayT(){return decayT;};
  double Get_vel0(){return vel0;};
  double Get_Vx0(){return Vx0;};// = vel0 * sin(theta0) * cos(phi0);
  double Get_Vy0(){return Vy0;};// = vel0 * sin(theta0) * sin(phi0);
  double Get_Vz0(){return Vz0;};// = vel0 * cos(theta0);
  double Get_Lmfp(){return Lmfp;};// = 12*D/(PI * vel0_avrg);//mm
  double Get_X0(){return X_0;};
  double Get_Y0(){return Y_0;};
  double Get_Z0(){return Z_0;};
  double Get_T0(){return T_0;};
  double Get_theta0(){return theta0;};
  double Get_phi0(){return phi0;};


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
  bool Get_flag_laser_cw(){return flag_laser_cw;};
  bool Get_flag_laser_pulse(){return flag_laser_pulse;};

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
  bool flag_laser_cw;
  bool flag_laser_pulse;
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
  T_0 = t; //Tbeam
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

  tempX = 0;tempY = 0;
  do
  {
    tempX = ((double) rand() / (RAND_MAX)) * 3 * sqrt( 2 * k * T/massMu );
    tempY = ((double) rand() / (RAND_MAX)) * 2 * sqrt( 2 * massMu/PI/k/T) * exp(-1.0);
  }while(tempY>sqrt(2/PI*pow(massMu/k/T,3))*pow(tempX,2)*exp(-massMu*tempX*tempX/(2*k*T)));
  
  vel0 = tempX*1000 ; // mm/s

  tempX = 0;tempY = 0;
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
    tempX = 0;tempY = 0;
    do {
      tempX = ((double) rand() / (RAND_MAX)) * 5;  // l=0~0.01 mm
      tempY = ((double) rand() / (RAND_MAX)) * 1/Lmfp;
    } while (tempY>1/Lmfp*exp(-tempX/Lmfp) );

    L = tempX;

    t = t + L/vel0;
    x = x + vx * (L/vel0);
    y = y + vy * (L/vel0);
    z = z + vz * (L/vel0);

    if( (z >= 0 || z <= -Thick) )continue;
    //if( (z <= 0 && z >= -Thick) )//(z <= 14.12 && z >= 7) )//(z <= -9 && z >= -11) || (z <= 9 && z >= 7) )
    //{

    tempX = TMath::ACos(-1 + 2 * ((double)rand()/(RAND_MAX)) );
    tempY = ((double) rand()/(RAND_MAX)) * 2* PI;

    theta = tempX;  // [0,+PI] radians
    phi = tempY; //[0,+2PI] radians

    vx = vel0 * sin(theta) * cos(phi);
    vy = vel0 * sin(theta) * sin(phi);
    vz = vel0 * cos(theta);
    //}

  }while( (z <= 0 && z >= -Thick && t<decayT));

  if(z>=0 && t<decayT) flag_sf = 1;
  
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

std::array<double,3> Muonium0::EmissionPosition( double t){

    double delT = t - (T_sf);
    double tempX, tempY;
    
    double x = X_sf;
    double y = Y_sf;
    double z = Z_sf;
    double vx = VX_sf;
    double vy = VY_sf;
    double vz = VZ_sf;

    //if( t < T_sf) cout<<"wrong"<<endl;
    //if( t > decayT + T_0) cout<<"wrong"<<endl;
  
  if(t < T_sf || t > decayT + T_0){
    x = std::numeric_limits<double>::lowest();
    y = std::numeric_limits<double>::lowest();
    z = std::numeric_limits<double>::lowest();
  }
  else{
    x = x + vx * (delT);
    y = y + vy * (delT);
    z = z + vz * (delT);
  }
  return {x, y, z};

}

std::array<double,3> Muonium0::Emission(){

    double delT = decayT - (T_sf - T_0);
    double tempX, tempY;

    double x = X_sf;
    double y = Y_sf;
    double z = Z_sf;
    double vx = VX_sf;
    double vy = VY_sf;
    double vz = VZ_sf;

    //std::array<double,3> position = EmissionPosition (T_sf + Tstep);
    std::array<double,3> position = EmissionPosition (T_sf);
/*
    for(int i = 0; i < nbinT; i++){
      
      if(Tstep*i >= delT)break;
      x = x + vx * (Tstep);
      y = y + vy * (Tstep);
      z = z + vz * (Tstep);
      t = t + Tstep;

      hZT2D->Fill(TBeam + t, z);
      if(fabs(y)<=20 && Tstep*i <= delT){
        if( (MCtype == 1 || MCtype == 3) && flag_xfree == 0 && fabs(x)>20)continue;
        if( z >= 1 && z <= 6) {hTlaserR->Fill(TBeam + t);}
        if( z >= (-6-Thick) && z <= (-1-Thick)) {hTlaserL->Fill(TBeam + t);}
      }
    }
*/
/*
    DecayX = X_sf + VX_sf * (decayT - t0);
    DecayY = Y_sf + VY_sf * (decayT - t0);
    DecayZ = Z_sf + VZ_sf * (decayT - t0);
    
    LaserX = X_sf + VX_sf * (tLaser - TBeam - t0 );
    LaserY = Y_sf + VY_sf * (tLaser - TBeam - t0 );
    LaserZ = Z_sf + VZ_sf * (tLaser - TBeam - t0 );
    LaserXp = VX_sf/VZ_sf;
    LaserYp = VY_sf/VZ_sf;
    LaserE = 0.5 * massMu * 1e-6 * (VX_sf*VX_sf + VY_sf*VY_sf + VZ_sf*VZ_sf);//v:mm/s, Ek: MeV
*/    
    return {x,y,z};

}

bool Muonium0::InsideLaserRegion_cw(double x, double y, double z){
  
  if( sqrt((z-2)*(z-2)+y*y)<=0.3 )return true;
  return false;


}

bool Muonium0::InsideLaserRegion_pulse(double x, double y, double z){

  if( sqrt((z-3)*(z-3)+y*y)<=2 )return true;
  return false;
  
}
/*

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

*/

#endif //__MUONIUM0_HH__
