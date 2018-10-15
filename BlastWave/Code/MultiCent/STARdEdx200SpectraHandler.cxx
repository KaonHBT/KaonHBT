#include "STARdEdx200SpectraHandler.h"
#include "BlastWaveFitter.h"
#include "DataForFit.h"
#include "TFile.h"
#include <cstdlib>
using namespace std;

STARdEdx200SpectraHandler::STARdEdx200SpectraHandler(BlastWaveFitter* aBWFitter, const char* aInputFileName)
  :IOHandler(aBWFitter,aInputFileName)
{
  mSpectra =1;
  mHBT=0;
  mV2=0;
}

void STARdEdx200SpectraHandler::initLoad(){
  // set Files
  setOutputFile("data/Spectra/STARdEdx200Spectra.fit.root\0");
  if(!mInputFile){
    mInputFile = new TFile("data/Spectra/STARdEdx200Spectra.root");
  }
  // Set centralities
  mNCent = 10;
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
  mNParticipant[9] = 1.;
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

void STARdEdx200SpectraHandler::load(){
  if(mCent==9){
    char tGraphName[100];
    sprintf(tGraphName,"STARdEdx%sPP",mParticleName[0]);
    mSpectras[0] = 
      mBWFitter->addSpectra((TGraphErrors*) mInputFile->Get(tGraphName),
			    MPi,0.5,1.2);
    sprintf(tGraphName,"STARdEdx%sPP",mParticleName[1]);
    mSpectras[1] = 
      mBWFitter->addSpectra((TGraphErrors*) mInputFile->Get(tGraphName),
			    MPi,0.5,1.2);
    sprintf(tGraphName,"STARdEdx%sPP",mParticleName[2]);
    mSpectras[2] = 
    mBWFitter->addSpectra((TGraphErrors*) mInputFile->Get(tGraphName),
			  MK,0.5,1.2);
    sprintf(tGraphName,"STARdEdx%sPP",mParticleName[3]);
    mSpectras[3] = 
      mBWFitter->addSpectra((TGraphErrors*) mInputFile->Get(tGraphName),
			    MK,0.5,1.2);
    sprintf(tGraphName,"STARdEdx%sPP",mParticleName[4]);
    mSpectras[4] = 
      mBWFitter->addSpectra((TGraphErrors*) mInputFile->Get(tGraphName),
			    MP,0.5,1.2);
    sprintf(tGraphName,"STARdEdx%sPP",mParticleName[5]);
    mSpectras[5] = 
      mBWFitter->addSpectra((TGraphErrors*) mInputFile->Get(tGraphName),
			    MP,0.5,1.2);
  }
  else{
    char tGraphName[100];
    sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[0],mCent);
    mSpectras[0] = 
      mBWFitter->addSpectra((TGraphErrors*) mInputFile->Get(tGraphName),
			    MPi,0.5,1.2);
    sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[1],mCent);
    mSpectras[1] = 
      mBWFitter->addSpectra((TGraphErrors*) mInputFile->Get(tGraphName),
			    MPi,0.5,1.2);
    sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[2],mCent);
    mSpectras[2] = 
    mBWFitter->addSpectra((TGraphErrors*) mInputFile->Get(tGraphName),
			  MK,0.5,1.2);
    sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[3],mCent);
    mSpectras[3] = 
      mBWFitter->addSpectra((TGraphErrors*) mInputFile->Get(tGraphName),
			    MK,0.5,1.2);
    sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[4],mCent);
    mSpectras[4] = 
      mBWFitter->addSpectra((TGraphErrors*) mInputFile->Get(tGraphName),
			    MP,0.5,1.2);
    sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[5],mCent);
    mSpectras[5] = 
      mBWFitter->addSpectra((TGraphErrors*) mInputFile->Get(tGraphName),
			    MP,0.5,1.2);
  }
}
			  
