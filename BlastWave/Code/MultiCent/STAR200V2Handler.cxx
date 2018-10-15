#include "STAR200V2Handler.h"
#include "BlastWaveFitter.h"
#include "DataForFit.h"
#include "TFile.h"
#include <cstdlib>
using namespace std;
STAR200V2Handler::STAR200V2Handler(BlastWaveFitter* aBWFitter, const char* aInputFileName)
  :IOHandler(aBWFitter,aInputFileName)
{
  mSpectra =0;
  mHBT=0;
  mV2=1;
}

void STAR200V2Handler::initLoad(){
  // set Files
  setOutputFile("data/v2/STARPidAuAu200.fit.root\0");
  if(!mInputFile){
    mInputFile = new TFile("data/v2/STAR200V2PRC.root");//data/v2/STARPidAuAu200Art.root");
  }
  // Set centralities
  mNCent = 9;
  mNParticipant = new double[mNCent];
  //double NParticipantSTARTOfdAuPP200Spectra[] = {
  //15.1211, 11.3394, 4.92598, 8.2, 
  //1.};
  mNParticipant[0] = 351.459;
  mNParticipant[1] = 299.039;
  mNParticipant[2] = 234.478;
  mNParticipant[3] = 167.192;
  mNParticipant[4] = 116.075;
  mNParticipant[5] = 76.779;
  mNParticipant[6] = 48.2101;
  mNParticipant[7] = 27.9088;
  mNParticipant[8] = 14.6283;
   
  // Set particle infos
  mNParticles = 3;
  mParticleName = new char*[mNParticles];
  char tParticleName[10][6] = {"Pi","K","P"};
  mParticleMass = new double[mNParticles];
  mParticleMass[0] = MPi;
  mParticleMass[1] = MK;
  mParticleMass[2] = MP;
  for(int ti=0; ti<mNParticles; ti++){
    mParticleName[ti] = new char[10];
    strcpy(mParticleName[ti],tParticleName[ti]);
  }
}

void STAR200V2Handler::load(){
  double tLowLimPt[3] = {0.4,0.,0.3};
  double tHighLimPt[3] = {0.85, 0.5, 1.1};
  for(int ti=0; ti<mNParticles; ti++){
    char tGraphName[100];
    sprintf(tGraphName,"%sC%i",mParticleName[ti],mCent);
    mBWFitter->addV2((TGraphErrors*) mInputFile->Get(tGraphName),
		     mParticleMass[ti],tLowLimPt[ti], tHighLimPt[ti]);
  }
}
			  
