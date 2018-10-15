#include "STAR200HBTHandler.h"
#include "BlastWaveFitter.h"
#include "DataForFit.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include <cstdlib>
using namespace std;
STAR200HBTHandler::STAR200HBTHandler(BlastWaveFitter* aBWFitter, const char* aInputFileName)
  //:STAR200asHBTHandler(aBWFitter,aInputFileName)
  :IOHandler(aBWFitter,aInputFileName)
{
  mV2=0;
  mHBT=1;
  mSpectra=0;
}

void STAR200HBTHandler::initLoad(){
  // set Files
  setOutputFile("data/STAR200HBT.fit.root\0");
  if(!mInputFile){
    mInputFile = new TFile("data/HBT/STARAuAu200.root");
  }
  // Set centralities
  mNCent = 6;
  mNParticipant = new double[mNCent];
  mNParticipant[0] = 351.5; // 0-5 C0
  mNParticipant[1] = 299.0;  // 5-10 C1
  mNParticipant[2] = 234.5;  // 10-20 C2
  mNParticipant[3] = 167.2; // 20-30 C3
  mNParticipant[4] = 96.5;  // 30-80 C4-5-7-8
  mNParticipant[5] = 30.2;  // 30-80 C4-5-7-8

  // Set particle infos
  mNParticles = 1;
  mParticleName = new char*[mNParticles];
  char tParticleName[10][6] = {"Pip","Pim","Kp","Km","Pp","Pm"};
  for(int ti=0; ti<mNParticles; ti++){
    mParticleName[ti] = new char[10];
    strcpy(mParticleName[ti],tParticleName[ti]);
  }
  mParticleMass = new double[mNParticles];
  mParticleMass[0] = MPi;
  // mParticleMass[6] = MPi;
  //mParticleMass[7] = MK;
  //mParticleMass[8] = MP;

}

void STAR200HBTHandler::load(){
  char tNameo[100];
  sprintf(tNameo,"pimROutC%i",mCent);
  char tNames[100];
  sprintf(tNames,"pimRSideC%i",mCent);
  char tNamel[100];
  sprintf(tNamel,"pimRLongC%i",mCent);
  mBWFitter->addHbtRadii((TGraphErrors*) mInputFile->Get(tNameo),
			 (TGraphErrors*) mInputFile->Get(tNames),
			 (TGraphErrors*) mInputFile->Get(tNamel),
			 MPi);
  mValidParticle[0]=1;
}



