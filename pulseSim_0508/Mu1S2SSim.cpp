#include <boost/serialization/array.hpp> //necessary depending on environment
#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include "CWLaser.h"
#include "PulseLaser.h"
#include "Muonium.h"
#include "TreeManager.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"


using namespace boost::numeric::odeint;

static constexpr int Nstates = Muonium::GetnRhoparams();
typedef std::array< double, Nstates > state_type; 

Muonium* mu = 0;
CWLaser cw( 0, 0, 2, 0.3, 0 ); // position (x,y,z) [mm], beam waist (1/e^2 radius) [mm], beam power [W] // essentially do not use but set detuning value function here
//PulseLaser pulse( 0, 0, 3.0, 11e-7, 2e-9, 1.0, 0.005 ); // ( position (x,y,z) [mm], time t [s] Pulse exist at fCenter*, FWHM pulse duration [s], beam waist [mm], pulse energy [J] )
//PulseLaser pulse( 0, 0, 3.0, 13e-7, 2e-9, 1.0, 0.005 ); // ( position (x,y,z) [mm], time t [s] Pulse exist at fCenter*, FWHM pulse duration [s], beam waist [mm], pulse energy [J] )
PulseLaser pulse( 0, 0, 3.0, 14e-7, 2e-9, 1.0, 0.005 ); // ( position (x,y,z) [mm], time t [s] Pulse exist at fCenter*, FWHM pulse duration [s], beam waist [mm], pulse energy [J] )

void OpticalBloch( const state_type &x , state_type &dxdt , double t )
{
  // current muonium position
  std::array<double,3> pos = mu->GetPosition( t );
  if( pos[2] < 0 ){
    // During the muonium is in the aerogel, the density matrix does not change.
    for(int i=0; i<Nstates; ++i) dxdt[i] = 0;
    return;
  }

  // laser intensity at the position
  //double cIntensity = cw.GetIntensity( pos, t );
  constexpr double cIntensity = 0.; // we do not use CW laser
  double pIntensity = pulse.GetIntensity( pos, t );
  double Linewidth = pulse.GetLinewidth(); //pulse laser linewidth
  
  // Rabi frequency at the position
  //double Omega = Muonium::GetOmega( cIntensity );
  double Omega = Muonium::GetOmega( pIntensity ); // pulse excitation scheme

  // Ionization rate at the position due to 244 (& 355)
  double GIon = Muonium::GetIonizationRate( cIntensity + pIntensity ); // 2S, ionization cross section is considered?

  // actual detuning
  
  // Optical-Bloch for the Muonium

  // original version: no stark mixing due to static electric field and assume Fourier-transform-limited laser linewidth
  //x[0]:rho_gg, x[1]:Re(rho_ge), x[2]:Im(rho_ge), x[3]:rho_ee, x[4]:rho_ion
  /*
  double dOmega = cw.GetDetuning();                                 // (2omega_L - omega_eg)
  dOmega -= Muonium::GetACStarkShift2S( cIntensity + pIntensity );  // 2S AC Stark shift
  dOmega += Muonium::GetACStarkShift1S( cIntensity + pIntensity );  // 1S AC Stark shift
  dOmega += Muonium::Get2ndDoppler( mu->GetAbsVelocity(t) );        // second-order Doppler shift
  //dOmega -= Muonium::Get1stDoppler( mu->GetVelocity(t) );           // first-order Doppler shift
  //dOmega -= 117437433*2*M_PI // Muonium recoil shift [rad/s]

  dxdt[0] = -  Omega * x[2] + (Muonium::GammaS * x[3]) - (Muonium::GammaDecay * x[0]);
  dxdt[1] =   dOmega * x[2] - (Muonium::GammaDecay + GIon/2 + Muonium::GammaS/2) * x[1]; 
  dxdt[2] = - dOmega * x[1] + Omega/2 * (x[0]-x[3]) - (Muonium::GammaDecay + GIon/2 + Muonium::GammaS/2) * x[2];
  dxdt[3] =    Omega * x[2] - (GIon + Muonium::GammaS + Muonium::GammaDecay) * x[3];
  dxdt[4] =     GIon * x[3] - (Muonium::GammaDecay * x[4]);
  */
  
  //current version, no static electric field
  //Reference:(Laser Physics, Vol. 12, No. 7, 2002, pp. 1–8.)
  //x[0]:rho_gg, x[1]:Re(rho_ge), x[2]:Im(rho_ge), x[3]:rho_ee, x[4]:rho_ion
  /*
  double dOmega = cw.GetDetuning();                                 // (2omega_L - omega_eg)
  dOmega -= 3.0 * Muonium::GetACStarkShift2S( cIntensity + pIntensity );  // 2S AC Stark shift
  dOmega += 3.0 * Muonium::GetACStarkShift1S( cIntensity + pIntensity );  // 1S AC Stark shift
  dOmega += Muonium::Get2ndDoppler( mu->GetAbsVelocity(t) );        // second-order Doppler shift
  dOmega -= Muonium::Get1stDoppler( mu->GetVelocity(t) );           // first-order Doppler shift
  //dOmega -= 117437433*2*M_PI // Muonium recoil shift [rad/s]
  
  dxdt[0] = - sqrt(2.0) * Omega * x[2] + (Muonium::GammaS * x[3]) - (Muonium::GammaDecay * x[0]);
  dxdt[1] =   dOmega * x[2] - ( (Muonium::GammaDecay + 3.0 * GIon/2 + Muonium::GammaS/2) + Linewidth ) * x[1];
  dxdt[2] = - dOmega * x[1] + Omega/sqrt(2.0) * (x[0]-x[3]) - ( (Muonium::GammaDecay + 3.0 * GIon/2 + Muonium::GammaS/2) + Linewidth ) * x[2];
  dxdt[3] =   sqrt(2.0) * Omega * x[2] - (GIon + Muonium::GammaS + Muonium::GammaDecay) * x[3];
  dxdt[4] =   GIon * x[3] - (Muonium::GammaDecay * x[4]);
  */

  //current version, with static electric field
  //Reference:(Laser Physics, Vol. 12, No. 7, 2002, pp. 1–8.)
  //1:1S, 2:2p1/2, 3:2S, 4:2p3/2 (4-level system + ionization)
  //4 population levels, 6*2 coherence parameters (Real and Imaginary), 1 ionization probability parameter
  //11, 22, 33, 44, Re12, Im12, Re13, Im13, Re14, Im14, Re23, Im23, Re24, Im24, Re34, Im34, ion
  
  double Detuning = cw.GetDetuning();// (2omega_L - omega_eg)
  Detuning += Muonium::Get2ndDoppler( mu->GetAbsVelocity(t) ); // second-order Doppler shift
  Detuning -= Muonium::Get1stDoppler( mu->GetVelocity(t) );    // first-order Doppler shift
  //dOmega -= 117437433*2*M_PI // Muonium recoil shift [rad/s]
  double AC2S = Muonium::GetACStarkShift2S( cIntensity + pIntensity );
  double AC1S = Muonium::GetACStarkShift1S( cIntensity + pIntensity );
  double AC2P = Muonium::GetACStarkShift2P( cIntensity + pIntensity );
  double GIon2P = Muonium::GetIonizationRate2P( cIntensity + pIntensity );; // 2P1/2, 2P3/2, (Laser Physics, Vol. 12, No. 7, 2002, pp. 1–8.)

  static constexpr double DCfield  = 0.060e3;//1.0e3; //[V/cm]
  static constexpr double DCmixing = 2.0 * M_PI * 2.2162e6 * DCfield; // sqrt(3)*e*a0*E (Laser Physics, Vol. 12, No. 7, 2002, pp. 1–8.)

  dxdt[0] = - sqrt(2.0) * Omega * x[7] + (Muonium::GammaS * x[2]) + (Muonium::GammaP * (x[1] + x[3])) - (Muonium::GammaDecay * x[0]);
  dxdt[1] = - (Muonium::GammaP + GIon2P + Muonium::GammaDecay) * x[1] - 2.0 * DCmixing * x[11];
  dxdt[2] =   sqrt(2.0) * Omega * x[7] - (GIon + Muonium::GammaS + Muonium::GammaDecay) * x[2] + 2.0 * DCmixing * x[11] + 2.0 * sqrt(2.0) * DCmixing * x[15] ;
  dxdt[3] = - (Muonium::GammaP + GIon2P + Muonium::GammaDecay) * x[3] - 2.0 * sqrt(2.0) * DCmixing * x[15];
  
  dxdt[4] =   (Muonium::Omega_2S2P1_2 + Detuning - 3.0 * AC2P + 3.0 * AC1S) * x[5] - (Muonium::GammaP + 3.0 * GIon2P + 2.0 * Muonium::GammaDecay + 2.0 * Linewidth)/2.0 * x[4] - Omega/sqrt(2.0) * x[11] - DCmixing * x[7]; 
  dxdt[5] = - (Muonium::Omega_2S2P1_2 + Detuning - 3.0 * AC2P + 3.0 * AC1S) * x[4] - (Muonium::GammaP + 3.0 * GIon2P + 2.0 * Muonium::GammaDecay + 2.0 * Linewidth)/2.0 * x[5] - Omega/sqrt(2.0) * x[10] + DCmixing * x[6]; 
  
  dxdt[6] =   (Detuning - 3.0 * AC2S + 3.0 * AC1S) * x[7] - (Muonium::GammaDecay + 3.0 * GIon/2.0 + Muonium::GammaS/2.0 + Linewidth) * x[6] - DCmixing * x[5] + sqrt(2.0) * DCmixing * x[9];
  dxdt[7] = - (Detuning - 3.0 * AC2S + 3.0 * AC1S) * x[6] - (Muonium::GammaDecay + 3.0 * GIon/2.0 + Muonium::GammaS/2.0 + Linewidth) * x[7] + DCmixing * x[4] - sqrt(2.0) * DCmixing * x[8] + Omega/sqrt(2.0) * (x[0]-x[2]);
  
  dxdt[8] =   (-Muonium::Omega_2S2P3_2 + Detuning - 3.0 * AC2P + 3.0 * AC1S) * x[9] - (Muonium::GammaP + 3.0 * GIon2P + 2.0 * Muonium::GammaDecay + 2.0 * Linewidth)/2.0 * x[8] - Omega/sqrt(2.0) * x[15] + DCmixing * x[7];
  dxdt[9] = - (-Muonium::Omega_2S2P3_2 + Detuning - 3.0 * AC2P + 3.0 * AC1S) * x[8] - (Muonium::GammaP + 3.0 * GIon2P + 2.0 * Muonium::GammaDecay + 2.0 * Linewidth)/2.0 * x[9] - Omega/sqrt(2.0) * x[14] - DCmixing * x[6];
  
  dxdt[10] =   (-Muonium::Omega_2S2P1_2 - AC2S + AC2P) * x[11] - (Muonium::GammaP + GIon2P + GIon + 2.0 * Muonium::GammaDecay)/2.0 * x[10] + Omega/sqrt(2.0) * x[5] + sqrt(2.0) * DCmixing * x[13];
  dxdt[11] = - (-Muonium::Omega_2S2P1_2 - AC2S + AC2P) * x[10] - (Muonium::GammaP + GIon2P + GIon + 2.0 * Muonium::GammaDecay)/2.0 * x[11] + Omega/sqrt(2.0) * x[4] - sqrt(2.0) * DCmixing * x[12] + DCmixing * (x[1]-x[2]);
  
  dxdt[12] =   (-Muonium::Omega_2S2P3_2-Muonium::Omega_2S2P1_2) * x[13] - (Muonium::GammaP + GIon2P + Muonium::GammaDecay) * x[12] + DCmixing * x[15] + sqrt(2.0) * DCmixing * x[11];
  dxdt[13] = - (-Muonium::Omega_2S2P3_2-Muonium::Omega_2S2P1_2) * x[12] - (Muonium::GammaP + GIon2P + Muonium::GammaDecay) * x[13] - DCmixing * x[14] - sqrt(2.0) * DCmixing * x[10];
  
  dxdt[14] =   (-Muonium::Omega_2S2P3_2 - AC2P + AC2S) * x[15] - (Muonium::GammaP + GIon2P + GIon + 2.0 * Muonium::GammaDecay)/2.0 * x[14] + Omega/sqrt(2.0) * x[9] + DCmixing * x[13];
  dxdt[15] = - (-Muonium::Omega_2S2P3_2 - AC2P + AC2S) * x[14] - (Muonium::GammaP + GIon2P + GIon + 2.0 * Muonium::GammaDecay)/2.0 * x[15] - Omega/sqrt(2.0) * x[8] - DCmixing * x[12] - sqrt(2.0) * DCmixing * (x[2]-x[3]);
  
  dxdt[16] =   GIon * x[2] + GIon2P * (x[1] + x[3]) - (Muonium::GammaDecay * x[16]);
  
  return;
}


int main(int argc, char **argv)
{

  //TFile* itf = new TFile( "../SeedFile/MuSeed.root" );
  TFile* itf = new TFile( "./Target_0417_Kapton_DG350.root" );
  TTree* seedTr = 0;
  itf->GetObject( "tree", seedTr );
  double X_sf=0, Y_sf=0, Z_sf=0, T_sf=0, VX_sf=0, VY_sf=0, VZ_sf=0;
  seedTr->SetBranchStatus( "*", 0 );
  seedTr->SetBranchStatus( "*_sf", 1 );
  seedTr->SetBranchAddress( "X_sf", &X_sf );
  seedTr->SetBranchAddress( "Y_sf", &Y_sf );
  seedTr->SetBranchAddress( "Z_sf", &Z_sf );
  seedTr->SetBranchAddress( "T_sf", &T_sf );
  seedTr->SetBranchAddress( "VX_sf", &VX_sf );
  seedTr->SetBranchAddress( "VY_sf", &VY_sf );
  seedTr->SetBranchAddress( "VZ_sf", &VZ_sf );
  
  cw.SetDetuning( 2*M_PI * atof(argv[2]) );
  cw.DumpSetting();
  pulse.DumpSetting();

  //constexpr double tStart = 0.28e-6;
  double tStart = pulse.GetPeakTime0() - 2.0e-8;
  constexpr int simCycle0 = 200;
  constexpr int simCycle  = 2000;
  double tPitch0 = tStart / simCycle0;
  constexpr double tPitch = 2e-11; // [s]
  
  bulirsch_stoer_dense_out< state_type > stepper( 1e-16, 1e-16, 1.0, 1.0 ); // looks OK?
  //bulirsch_stoer_dense_out< state_type > stepper( 1e-12, 1e-12, 1.0, 1.0 ); // sometimes cause unphysical results
  TreeManager tManager( argv[1], simCycle+1, tPitch, &cw, &pulse, false ); // ( filename, nPoints, tPitch, cw, pulse, saveFlag )
  //TreeManager tManager( argv[1], simCycle+1, tPitch, &cw, &pulse, true ); // ( filename, nPoints, tPitch, cw, pulse, saveFlag )

  const int nEntries = seedTr->GetEntries();
  for( int entry=0; entry<nEntries; entry++ ){
    if( entry % (nEntries/100) == 0 ){
      std::cout << entry << "/" << nEntries << "\r" << std::flush;
    }

    seedTr->GetEntry( entry );
    mu = new Muonium( {X_sf, Y_sf, Z_sf}, T_sf, {VX_sf, VY_sf, VZ_sf} );
    //mu->SetTemperature( 350. );
    tManager.Initialize( mu );
    state_type x = {};
    x[0] = 1.0;
    integrate_const( stepper, OpticalBloch, x , 0.0 , tPitch0*simCycle0, tPitch0 );
    integrate_const( stepper, OpticalBloch, x , tStart , tStart+tPitch*simCycle, tPitch , std::ref(tManager) );
    tManager.Fill();
    
    delete mu;
  }


  tManager.Write();

}



 

  
