#include "STAR200V2V4Handler.h"
#include "BlastWaveFitter.h"
#include "DataForFit.h"
#include "TFile.h"
#include <iostream>
#include <cstdlib>
using namespace std;

STAR200V2V4Handler::STAR200V2V4Handler(BlastWaveFitter* aBWFitter,
				       int aV4,
				       const char* aInputFileName)
  :IOHandler(aBWFitter,aInputFileName),mV4(aV4)
{
  mSpectra =1;
  mHBT=0;
  mV2=1;
}

void STAR200V2V4Handler::initLoad(){
  // set Files
  if(mV4){
    setOutputFile("data/STAR200ChV2V4Spectra.fit.root\0");
  }
  else{
    setOutputFile("data/STAR200ChV2Spectra.fit.root\0");
  }
  if(!mInputFile){
    mInputFile = new TFile("data/v2/STAR200V2V4.root");
  }
  if(!mSpectraInputFile){
    mSpectraInputFile = new TFile("data/Spectra/STARdEdx200Spectra.root");
  }
  // Set centralities
  mNCent = 9;
  mNParticipant = new double[mNCent];
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
  mNParticles = 6;
  mParticleName = new char*[mNParticles];
  char tParticleName[10][6] = {"Pip","Pim","Kp","Km","Pp","Pm"};
  for(int ti=0; ti<mNParticles; ti++){
    mParticleName[ti] = new char[10];
    strcpy(mParticleName[ti],tParticleName[ti]);
  }
  mParticleMass = new double[mNParticles];
  mParticleMass[0] = MPi;
  mParticleMass[1] = MPi;
  mParticleMass[2] = MK;
  mParticleMass[3] = MK;
  mParticleMass[4] = MP;
  mParticleMass[5] = MP;
}

void STAR200V2V4Handler::load(){
  char tGraphName[100];
  sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[0],mCent);
  mSpectras[0] = 
    mBWFitter->addSpectra((TGraphErrors*) mSpectraInputFile->Get(tGraphName),
			  MPi,0.5,1.2);
  sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[1],mCent);
  mSpectras[1] = 
    mBWFitter->addSpectra((TGraphErrors*) mSpectraInputFile->Get(tGraphName),
			  MPi,0.5,1.2);
  sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[2],mCent);
  mSpectras[2] = 
    mBWFitter->addSpectra((TGraphErrors*) mSpectraInputFile->Get(tGraphName),
			  MK,0.5,1.5);
  sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[3],mCent);
  mSpectras[3] = 
    mBWFitter->addSpectra((TGraphErrors*) mSpectraInputFile->Get(tGraphName),
			    MK,0.5,1.5);
  sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[4],mCent);
  mSpectras[4] = 
    mBWFitter->addSpectra((TGraphErrors*) mSpectraInputFile->Get(tGraphName),
			  MP,0.5,2.);
  sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[5],mCent);
  mSpectras[5] = 
    mBWFitter->addSpectra((TGraphErrors*) mSpectraInputFile->Get(tGraphName),
			  MP,0.5,2.);
  sprintf(tGraphName,"V2C%i",mCent);
  cout << tGraphName << " " << mInputFile << endl;
  cout <<  mInputFile->Get(tGraphName) << endl;
  mBWFitter->addV2((TGraphErrors*) mInputFile->Get(tGraphName),MPi,0.,1.2);
  sprintf(tGraphName,"V4C%i",mCent);
  if(mV4) mBWFitter->addV4((TGraphErrors*) mInputFile->Get(tGraphName),MPi,0.,1.2);
}


void STAR200V2V4Handler::fixParameters(){
  mBWFitter->freeAllParameters();
  if(!mV4) mBWFitter->fixRho4();
  mBWFitter->fixAs(); // need an extra switch
  mBWFitter->fixRy();
  mBWFitter->fixTau();
  mBWFitter->fixDeltat();
}
			  
