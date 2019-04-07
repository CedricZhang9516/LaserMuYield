#ifndef __TREEMANAGER_HH__
#define __TREEMANAGER_HH__

#include <iostream>
#include <fstream>
#include <array>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "Muonium.h"
#include "CWLaser.h"
#include "PulseLaser.h"

typedef std::array< double , 5 > state_type;

using namespace std;

class TreeManager
{

 public:
  TreeManager( char* filename, int nPoints, double tPitch, CWLaser* cw, PulseLaser* pulse, bool saveFlag );
  ~TreeManager();
  void operator()( const state_type &x, const double t );
  void Initialize( Muonium* mu );
  void Fill() { fTree->Fill();};
  void Write() { fFile->Write();  };
  double* GetRho(){
    static double r[5]; 
    r[0] = fFinalRho[0];
    r[1] = fFinalRho[1];
    r[2] = fFinalRho[2];
    r[3] = fFinalRho[3];
    r[4] = fFinalRho[4];
    return r;
  };

 private:
  Muonium* fMu;
  CWLaser* fCW;
  PulseLaser* fPulse;
  bool fDetailSaveFlag;

  TFile* fFile;
  TTree* fTree;
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

};


TreeManager::TreeManager( char* filename, int nPoints, double tPitch, CWLaser* cw, PulseLaser* pulse, bool saveFlag=false ):
fCW(cw), fPulse(pulse), fDetailSaveFlag(saveFlag), fNPoints( nPoints ), fTPitch(tPitch), fPIndex(0)
{

  std::string filename2 = std::string(filename);
  filename2 = "/Users/zhangce/WorkArea/LaserMuYield/Root/" + filename2;

  fFile = new TFile( filename2.c_str(), "RECREATE" );
  fSettingTree = new TTree( "Setting", "Setting" );
  // Calculation setting
  fSettingTree->Branch( "NPoints", &fNPoints, "NPoints/I" );
  fSettingTree->Branch( "TPitch", &fTPitch, "TPitch/D" );

  // CW setting
  fDetuning = cw->GetDetuning() / (2*M_PI); // [/s]
  cw->GetCenterPosition( fCWPos[0], fCWPos[1], fCWPos[2] );
  fCWWaist0  = cw->GetWaist0();
  fCWPower   = cw->GetPower();
  fCWOnTime  = cw->GetOnTime();
  fCWOffTime = cw->GetOffTime();
  fSettingTree->Branch( "Detuning", &fDetuning, "Detuning/D" );
  fSettingTree->Branch( "CWPosition", fCWPos, "CWPosition[3]/D" );
  fSettingTree->Branch( "CWWaist0", &fCWWaist0, "CWWaist0/D" );
  fSettingTree->Branch( "CWPower", &fCWPower, "CWPower/D" );
  fSettingTree->Branch( "CWOnTime", &fCWOnTime, "CWOnTime/D" );
  fSettingTree->Branch( "CWOffTime", &fCWOffTime, "CWOffTime/D" );

  // Pulse setting
  pulse->GetCenterPosition( fPulsePos[0], fPulsePos[1], fPulsePos[2] );
  fPulseWaist0 = pulse->GetWaist0();
  fPulseEnergy = pulse->GetEnergy();
  fSettingTree->Branch( "PulsePosition", fPulsePos, "PulsePosition[3]/D" );
  fSettingTree->Branch( "PulseWaist0", &fPulseWaist0, "PulseWaist0/D" );
  fSettingTree->Branch( "PulseEnergy", &fPulseEnergy, "PulseEnergy/D" );

  fSettingTree->Fill();



  fTree = new TTree( "MuTree", "MuTree" );
  fPos = new double[fNPoints][3]();
  fTime = new double[fNPoints]();
  fRho = new double[fNPoints][5]();
  fIntensity = new double[fNPoints][2]();
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

  fMu = 0;
}


TreeManager::~TreeManager()
{
  fFile->Close();
  delete fFile;
  delete[] fPos;
  delete[] fTime;
  delete[] fRho;
  delete[] fIntensity;
}


void TreeManager::Initialize( Muonium* mu )
{
  fMu = mu;
  fPIndex = 0;

  fMu->GetVelocity( 0, fVelocity[0], fVelocity[1], fVelocity[2] );
  fTemperature = mu->GetTemperature();
  
  fSurfaceTime = fMu->GetStartTime();
  fMu->GetPosition( fSurfaceTime, fSurfacePos[0], fSurfacePos[1], fSurfacePos[2] );
  fMu->GetVelocity( fSurfaceTime, fSurfaceVelocity[0], fSurfaceVelocity[1], fSurfaceVelocity[2] );
  fSurfaceAbsVelocity = fMu->GetAbsVelocity( fSurfaceTime );
}



void TreeManager::operator()( const state_type &x, const double t )
{
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
}
    
  
#endif //__TREEMANAGER_HH__  


