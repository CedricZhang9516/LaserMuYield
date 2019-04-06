#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include "CWLaser.h"
#include "PulseLaser.h"
#include "Muonium.h"
//#include "SOA.h"
#include "TreeManager.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"


using namespace boost::numeric::odeint;

typedef std::array< double, 5 > state_type; // (rho_gg, Re(rho_ge_, Im(rho_ge), rho_ee, rho_ion)



Muonium* mu = 0;
CWLaser cw( 0, 0, 2, 0.3, 10 );
PulseLaser pulse( 0, 0, 3, 2e-6, 2e-9, 2.0, 0.1 );
//SOA soa( {0, 10}, {0, 0}, 6.1e-6, 1/25, 1e-7, 0 );

void OpticalBloch( const state_type &x , state_type &dxdt , double t )
{
  // current muonium position
  std::array<double,3> pos = mu->GetPosition( t );
  if( pos[2] < 0 ){
    // During the muonium is in the aerogel, the density matrix does not change.
    dxdt[0] = 0;
    dxdt[1] = 0;
    dxdt[2] = 0;
    dxdt[3] = 0;
    dxdt[4] = 0;
    return;
  }

  // laser intensity at the position
  double cIntensity = cw.GetIntensity( pos, t );
  double pIntensity = pulse.GetIntensity( pos, t );

  // Rabi frequency at the position
  double Omega = Muonium::GetOmega( cIntensity );

  // Ionization rate at the position due to 244 & 355
  double GIon = Muonium::GetIonizationRate( cIntensity + pIntensity );

  // actual detuning
  double dOmega = cw.GetDetuning();                                 // (2omega_L - omega_eg)
  dOmega -= Muonium::GetACStarkShift2S( cIntensity + pIntensity );  // 2S AC Stark shift
  dOmega += Muonium::GetACStarkShift1S( cIntensity + pIntensity );  // 1S AC Stark shift
  dOmega += Muonium::Get2ndDoppler( mu->GetAbsVelocity(t) );        // second doppler shift
  //dOmega += Muonium::GetDCStarkShift( soa.GetEField(pos.at(2),t) ); // DC Stark shift due to SOA


  // Optical-Bloch for the Muonium
  dxdt[0] = -  Omega * x[2] + (Muonium::GammaS * x[3]) - (Muonium::GammaDecay * x[0]);
  dxdt[1] =   dOmega * x[2] - (Muonium::GammaDecay + GIon/2 + Muonium::GammaS/2) * x[1];
  dxdt[2] = - dOmega * x[1] + Omega/2 * (x[0]-x[3]) - (Muonium::GammaDecay + GIon/2 + Muonium::GammaS/2) * x[2];
  dxdt[3] =    Omega * x[2] - (GIon + Muonium::GammaS + Muonium::GammaDecay) * x[3];
  dxdt[4] =     GIon * x[3] - (Muonium::GammaDecay * x[4]);

  return;
}


int main(int argc, char **argv)
{
  /*
    ./DecayEffect filename detune[Hz] (./DecayEffect +1kHz.root 1000)
   */


  //TFile* itf = new TFile( "/Users/taka/Documents/SPAN/Muonium/20181230Simulator/MuSeed.root" );
  TFile* itf = new TFile("/Users/zhangce/WorkArea/LaserMuYield/input.root");
  //TFile* itf = new TFile( "/Users/taka/Documents/SPAN/Muonium/20190106SlowMuEnhancement/SlowMuSeed1E7.root" ); // slower than 4e6 mm/s
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


  cw.SetDetuning( 2*M_PI * atoi(argv[2]) );
  cw.DumpSetting();
  pulse.DumpSetting();


  constexpr double tPitch = 1e-8; // [s]
  //constexpr int simCycle = 610;
  constexpr int simCycle = 210;


  bulirsch_stoer_dense_out< state_type > stepper( 1e-16, 1e-16, 1.0, 1.0 );
  TreeManager tManager( argv[1], simCycle+1, tPitch, &cw, &pulse, false );

  const int nEntries = seedTr->GetEntries();
  for( int entry=0; entry<nEntries; entry++ ){
    if( entry % (nEntries/100) == 0 ){
      std::cout << entry << "/" << nEntries << "\r" << std::flush;
    }

    seedTr->GetEntry( entry );
    mu = new Muonium( {X_sf, Y_sf, Z_sf}, T_sf, {VX_sf, VY_sf, VZ_sf} );
    //mu->SetTemperature( 350. );
    tManager.Initialize( mu );
    state_type x = { 1.0, 0, 0, 0, 0 };
    integrate_const( stepper, OpticalBloch, x , 0.0 , tPitch*simCycle, tPitch , std::ref(tManager) );
    tManager.Fill();

    delete mu;
  }


  tManager.Write();

}
