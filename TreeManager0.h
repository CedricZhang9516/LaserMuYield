#ifndef __TREEMANAGER0_HH__
#define __TREEMANAGER0_HH__

#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "Muonium0.h"
#include "Muonium.h"
#include "CWLaser.h"
#include "PulseLaser.h"

typedef std::array< double , 5 > state_type;

using namespace std;

class TreeManager0
{

 public:
  TreeManager0(const char* filename);
  ~TreeManager0();
  void operator()( const state_type &x, const double t );
  void Initialize( Muonium0* mu );

  void Fill();// { fTree->Fill(); };
  void Write() { fFile->Write(); };
  void Close();

 private:
  Muonium0* fMu;
  CWLaser* fCW;
  PulseLaser* fPulse;
  bool fDetailSaveFlag;

  TFile* fFile;
  TTree* fTree;
  //TTree* tree;
  TTree* fSettingTree;

  // Calculation setting
  int fNPoints;// The number of points in an event
  double fTPitch; // [s]

  // CW setting
  double fDetuning;  // [Hz]
  double fCWPos[3];  // [mm]
  double fCWWaist0;
  double fCWPower;
  double fCWOnTime;
  double fCWOffTime;

  // Pulse setting
  double fPulsePos[3];
  double fPulseWaist0;
  double fPulseEnergy;

  int fPIndex;
  double fTemperature;
  double (*fPos)[3];     // (x , y, z) array
  double fVelocity[3]; // (vx, vy, vz)
  double* fTime;     // (t) array
  double (*fRho)[5];     // (rho_gg, Re(rho_ge), Im(rho_ge), rho_ee, rho_ion) array
  double (*fIntensity)[2];// Intensity array (244, 355)

  double fSurfacePos[3];
  double fSurfaceVelocity[3];
  double fSurfaceAbsVelocity;
  double fSurfaceTime;
  double fFinalTime;
  double fFinalPos[3];
  double fFinalVelocity[3];
  double fFinalAbsVelocity;
  double fFinalRho[5];

  double X_sf;// = -1; 
  double Y_sf;// = -1; 
  double Z_sf;// = -1; 
  double VX_sf;// = -1; 
  double VY_sf;// = -1; 
  double VZ_sf;// = -1;
  double T_sf;// = -1;
  double theta_sf;// = -1;
  double phi_sf;// = -1;

};


TreeManager0::TreeManager0(const char* filename)
{

  string filename2 = std::string(filename);
  filename2 = "/Users/zhangce/WorkArea/LaserMuYield/Root/Target_" + filename2;
  fFile = new TFile( filename2.c_str(), "RECREATE" );
  fSettingTree = new TTree( "Setting", "Setting" );
  fSettingTree->Fill();

  fTree = new TTree( "tree", "tree" );

  fTree->Branch("X_sf",&X_sf,"X_sf/D");
  fTree->Branch("Y_sf",&Y_sf,"Y_sf/D");
  fTree->Branch("Z_sf",&Z_sf,"Z_sf/D");
  fTree->Branch("VX_sf",&VX_sf,"VX_sf/D");
  fTree->Branch("VY_sf",&VY_sf,"VY_sf/D");
  fTree->Branch("VZ_sf",&VZ_sf,"VZ_sf/D");
  fTree->Branch("T_sf",&T_sf,"T_sf/D");
  fTree->Branch("theta_sf",&theta_sf,"theta_sf/D");
  fTree->Branch("phi_sf",&phi_sf,"phi_sf/D");

/*
  fPos = new double[fNPoints][3]();
  fTime = new double[fNPoints]();
  fRho = new double[fNPoints][5]();
  fIntensity = new double[fNPoints][2]();
*/
  /*
  fTree->Branch( "Temperature", &fTemperature, "Temperature/D" );
  if( fDetailSaveFlag ){
    fTree->Branch( "Position", fPos, Form( "Pos[%d][3]/D", fNPoints ) );
    fTree->Branch( "Velocity", fVelocity, "Velocity[3]/D" );
    fTree->Branch( "Time", fTime, Form( "Time[%d]/D", fNPoints) );
    fTree->Branch( "Rho", fRho, Form( "Rho[%d][5]/D", fNPoints) );
    fTree->Branch( "Intensity", fIntensity, Form( "Intensity[%d][2]/D", fNPoints ) );
  }
  fTree->Branch( "SurfacePosition", fSurfacePos, "SurfacePos[3]/D" );
  fTree->Branch( "SurfaceTime", &fSurfaceTime, "SurfaceTime/D" );
  fTree->Branch( "SurfaceVelocity", fSurfaceVelocity, "SurfaceVelocity[3]/D" );
  fTree->Branch( "SurfaceAbsVelocity", &fSurfaceAbsVelocity, "SurfaceAbsVelocity/D" );

  fTree->Branch( "FinalPosition", fFinalPos, "FinalPos[3]/D" );
  fTree->Branch( "FinalTime", &fFinalTime, "FinalTime/D" );
  fTree->Branch( "FinalVelocity", fFinalVelocity, "FinalVelocity[3]/D" );
  fTree->Branch( "FinalAbsVelocity", &fFinalAbsVelocity, "FinalAbsVelocity/D" );
  fTree->Branch( "FinalRho", fFinalRho, "FinalRho[5]/D" );
  */
  fMu = 0;
}


TreeManager0::~TreeManager0()
{
  fFile->Close();
  delete fFile;
  delete[] fPos;
  delete[] fTime;
  delete[] fRho;
  delete[] fIntensity;
}

void TreeManager0::Close(){
  fFile->Close();

}

void TreeManager0::Initialize( Muonium0* mu )
{
  fMu = mu;
  /*
  fPIndex = 0;

  fMu->GetVelocity( 0, fVelocity[0], fVelocity[1], fVelocity[2] );
  fTemperature = mu->GetTemperature();
  
  fSurfaceTime = fMu->GetStartTime();
  fMu->GetPosition( fSurfaceTime, fSurfacePos[0], fSurfacePos[1], fSurfacePos[2] );
  fMu->GetVelocity( fSurfaceTime, fSurfaceVelocity[0], fSurfaceVelocity[1], fSurfaceVelocity[2] );
  fSurfaceAbsVelocity = fMu->GetAbsVelocity( fSurfaceTime );
  */
}

void TreeManager0::Fill(){

  //fMu = mu;
  fFile->cd();

  X_sf = fMu->Get_X_sf();
  Y_sf = fMu->Get_Y_sf();
  Z_sf = fMu->Get_Z_sf();
  VX_sf = fMu->Get_VX_sf();
  VY_sf = fMu->Get_VY_sf();
  VZ_sf = fMu->Get_VZ_sf();
  T_sf = fMu->Get_T_sf();
  theta_sf = fMu->Get_theta_sf();
  phi_sf = fMu->Get_phi_sf();

  fTree->Fill();

}





void TreeManager0::operator()( const state_type &x, const double t )
{
  /*
  if( fDetailSaveFlag ){
    fTime[fPIndex] = t;
    fMu->GetPosition( fTime[fPIndex], fPos[fPIndex][0], fPos[fPIndex][1], fPos[fPIndex][2] );
    fIntensity[fPIndex][0] = fCW->GetIntensity( fPos[fPIndex][0], fPos[fPIndex][1], fPos[fPIndex][2] );
    fIntensity[fPIndex][1] = fPulse->GetIntensity( fPos[fPIndex][0], fPos[fPIndex][1], fPos[fPIndex][2], t );
    
    for( int i=0; i<5; i++ ){
      fRho[fPIndex][i] = x[i];
    }
  }

  fPIndex++;

  if( fPIndex == fNPoints ){
    // Last step
    fFinalTime = t;
    for( int i=0; i<3; i++){
      fMu->GetPosition( t, fFinalPos[0], fFinalPos[1], fFinalPos[2] );
      fMu->GetVelocity( t, fFinalVelocity[0], fFinalVelocity[1], fFinalVelocity[2] );
    }
    fFinalAbsVelocity = fMu->GetAbsVelocity( t );
    for( int i=0; i<5; i++ ){
      fFinalRho[i] = x[i];
    }
  }
  */
}
    
  
#endif //__TREEMANAGER0_HH__  


