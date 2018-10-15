#include "DataTestHandler.h"
#include "BlastWaveFitter.h"
#include "DataForFit.h"

#include "TFile.h"
#include "TGraphErrors.h"

#include <cstdlib>
using namespace std;

DataTestHandler::DataTestHandler(BlastWaveFitter* aBWFitter, const char* aInputFileName)
  :IOHandler(aBWFitter,aInputFileName)
{
  mSpectra =1;
  mHBT=0;
  mV2=1;
}

void DataTestHandler::initLoad(){
  // set Files
  setOutputFile("data/DataTest.fit.root\0");
  if(!mSpectraInputFile){
    mSpectraInputFile = new TFile("data/Spectra/STARdEdx200Spectra.root");
  }
  if(!mInputFile){
    mInputFile = new TFile("data/v2/V2V4ForArt.root");
  }
  // Set centralities
  mNCent = 1;
  mNParticipant = new double[mNCent];
  mNParticipant[0] = 300.; // 0-5 C0
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

  // Set blast wave fitter parameters
  mBWFitter->fixAs();

}


void DataTestHandler::loadSpectra(int iPart, int iStat){
  int tCentMerged[] = {0,9};
  if(tCentMerged[mCent] == tCentMerged[mCent+1]-1){ // no merging
    char tGraphName[100];
    sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[iPart],mCent);
    mSpectras[iPart] = 
      mBWFitter->addSpectra((TGraphErrors*) mSpectraInputFile->Get(tGraphName),
			    mParticleMass[iPart],0.,1.2,iStat);
  }
  else{
    int tCent = tCentMerged[mCent];
    int tN(0);
    double* tX(0);
    double* tEX(0);
    double* tY(0);
    double* tEY(0);
    while(tCent<tCentMerged[mCent+1]){
      char tGraphName[100];
      sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[iPart],tCent);
      TGraphErrors* tG = (TGraphErrors*) mSpectraInputFile->Get(tGraphName);
      if(!tN){
	tN = tG->GetN();
	tX = tG->GetX();
	tEX = tG->GetEX();
	tY = tG->GetY();
	tEY = tG->GetEY();
	for(int ti=0; ti<tN; ti++){
	  tEY[ti] = tEY[ti]*tEY[ti];
	}
      }
      else{
	double* ttY = tG->GetY();
	double* ttEY = tG->GetEY();
	for(int ti=0; ti<tN; ti++){
	  tY[ti]+= ttY[ti];
	  tEY[ti] += (ttEY[ti]*ttEY[ti]);
	}
      }
      tCent++;
    }
    int tNCent = tCentMerged[mCent+1]-tCentMerged[mCent];
    for(int ti=0; ti<tN; ti++){
      tY[ti] /= tNCent;
      tEY[ti] = sqrt(tEY[ti])/tNCent;
    }
    mSpectras[iPart] = 
      mBWFitter->addSpectra(new TGraphErrors(tN,tX,tY,tEX,tEY),
			    mParticleMass[iPart],0.5,1.2);   
  }
}


void DataTestHandler::load(){
  int tStat[] = {1,1,0,0,0,0};
  for(int ti=0; ti<mNParticles; ti++){
    loadSpectra(ti,tStat[ti]); // ugly hack
  }

  mBWFitter->addV2((TGraphErrors*) mInputFile->Get("v2Pi"),MPi,0.2,1.,1);
  mBWFitter->addV4((TGraphErrors*) mInputFile->Get("v4Pi"),MPi,0.2,1.,1);
  mBWFitter->addV2((TGraphErrors*) mInputFile->Get("v2Pr"),MP,0.,2.);
  mBWFitter->addV4((TGraphErrors*) mInputFile->Get("v4Pr"),MP,0.,2.);

  //for(int ti=3; ti<mNParticles; ti++){
  mValidParticle[1]=0;
  mValidParticle[3]=0;  
  mValidParticle[5]=0;

}
			  

void DataTestHandler::close(){
  IOHandler::close();
  mSpectraInputFile->Close();
}

void DataTestHandler::fixParameters(){
  mBWFitter->freeAllParameters();
  mBWFitter->fixAs(); // need an extra switch
  mBWFitter->fixRy();
  mBWFitter->fixTau();
  mBWFitter->fixDeltat();
  //mBWFitter->fixf4();
  //mBWFitter->fixf2();
  //mBWFitter->fixRx();
}
