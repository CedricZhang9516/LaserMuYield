#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include "CWLaser.h"
#include "PulseLaser.h"
#include "Muonium.h"
#include "Muonium0.h"
//#include "SOA.h"
#include "TreeManager.h"
#include "TreeManager0.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"


using namespace boost::numeric::odeint;
using namespace std;

typedef std::array< double, 5 > state_type; // (rho_gg, Re(rho_ge_, Im(rho_ge), rho_ee, rho_ion)


Muonium* mu = 0;
Muonium0* mu0 = 0;

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
  


  double x_dec, y_dec, z_dec, t_dec;
  int Nentries;
  //TFile *Fdlinefile = new TFile("/Users/zhangce/WorkArea/LaserMuYield/Root/SimBeamStop_0406_365_tot7.1e8.root");
  TFile *Fdlinefile = new TFile("/Users/zhangce/WorkArea/LaserMuYield/Root/SimBeamStop_0406_DG275_SUS_tot.root");
  TTree * Tdlinefile = (TTree*) Fdlinefile->Get("position");
  Tdlinefile->SetBranchAddress("x", &x_dec);
  Tdlinefile->SetBranchAddress("y", &y_dec);
  Tdlinefile->SetBranchAddress("z", &z_dec);
  Tdlinefile->SetBranchAddress("glbt_gen", &t_dec);
  Nentries = Tdlinefile->GetEntries();

  TreeManager0 tManager0( argv[1] );
	
  //Nentries = 1000;

  for( int entry=0; entry<Nentries; entry++ ){
    if( entry % (Nentries/1000) == 0 ){
      std::cout << entry << "/" << Nentries << "\r" << std::flush;
    }
    Tdlinefile->GetEntry(entry);
    //t_dec from the stopping simulation glbt_gen, the mean in the exaple root file is about 0.35 us.
    mu0 = new Muonium0( x_dec, y_dec, z_dec, (t_dec-350)*1e-9 );
    mu0->Diffusion();
    //mu0->Emission();
    if(mu0->Get_flag_sf() == 1)tManager0.Fill(mu0);
    //if(mu0->Get_flag_laser_cw() == 1)tManager0.Fill(mu0);
  }

  tManager0.Write();
  //tManager0.Close();

  cout<<"filename Target_"<<argv[1]<<" closed."<<endl;



  /*
    ./DecayEffect filename detune[Hz] (./DecayEffect +1kHz.root 1000)
   */

  //TFile* itf = new TFile( "/Users/taka/Documents/SPAN/Muonium/20181230Simulator/MuSeed.root" );
  
  string filename2 = std::string(argv[1]);
  filename2 = "/Users/zhangce/WorkArea/LaserMuYield/Root/Target_" + filename2;
  TFile* itf = new TFile(filename2.c_str());

  //TFile* itf = new TFile("/Users/zhangce/WorkArea/LaserMuYield/Root/masuda_san_1224_tree_Type3_D87000_T322_Nrepeat3231566_Xfree0_Thick7.12_NewGeo0_.root");
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
  
  //TreeManager tManager( argv[1], simCycle+1, tPitch, &cw, &pulse, false );
  TreeManager tManager( argv[1], simCycle+1, tPitch, &cw, &pulse, true );

  const int nEntries = seedTr->GetEntries();
  for( int entry=0; entry<nEntries; entry++ ){
    if( entry % (nEntries/10) == 0 ){
      std::cout << entry << "/" << nEntries << "\r" << std::flush;
    }

    seedTr->GetEntry( entry );
    mu = new Muonium( {X_sf, Y_sf, Z_sf}, T_sf, {VX_sf, VY_sf, VZ_sf} );
    //mu = new Muonium( {X_sf, Y_sf, Z_sf}, 0 , {VX_sf, VY_sf, VZ_sf} );
    //mu->SetTemperature( 350. );
    tManager.Initialize( mu );
    state_type x = { 1.0, 0, 0, 0, 0 };
    integrate_const( stepper, OpticalBloch, x , 0.0 , tPitch*simCycle, tPitch , std::ref(tManager) );
    tManager.Fill();

    //testing

    double* Rho = tManager.GetRho();
    //cout<<"Rho "<<*(Rho+0)<<" "<<*(Rho+1)<<" "<<*(Rho+2)<<" "<<*(Rho+3)<<" "<<*(Rho+4)<<endl;    

    double test_Velocity[3]; // (vx, vy, vz)
    mu->GetVelocity( 0, test_Velocity[0], test_Velocity[1], test_Velocity[2] );
  	
    /*
    cout<<"test_Velocity "<<test_Velocity[0]<<" "<<test_Velocity[1]<<" "<<test_Velocity[2]<<endl;

  	double test_Temperature = mu->GetTemperature();
  	double test_SurfaceTime = mu->GetStartTime();
  	cout<<"test_Temperature "<<test_Temperature<<endl;
  	cout<<"test_SurfaceTime "<<test_SurfaceTime<<endl;
  	
  	double test_SurfacePos[3]; // (vx, vy, vz)
  	double test_SurfaceVelocity[3];
  	mu->GetPosition( test_SurfaceTime, test_SurfacePos[0], test_SurfacePos[1], test_SurfacePos[2] );
  	mu->GetVelocity( test_SurfaceTime, test_SurfaceVelocity[0], test_SurfaceVelocity[1], test_SurfaceVelocity[2] );
  	cout<<"test_SurfacePos "<<test_SurfacePos[0]<<" "<<test_SurfacePos[1]<<" "<<test_SurfacePos[2]<<endl;
  	cout<<"test_SurfaceVelocity "<<test_SurfaceVelocity[0]<<" "<<test_SurfaceVelocity[1]<<" "<<test_SurfaceVelocity[2]<<endl;
    */

    delete mu;
  }

  tManager.Write();


}
