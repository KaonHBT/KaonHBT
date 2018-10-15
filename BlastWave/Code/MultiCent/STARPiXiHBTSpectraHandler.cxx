#include "STARPiXiHBTSpectraHandler.h"
#include "BlastWaveFitter.h"
#include "DataForFit.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include <cstdlib>
using namespace std;

STARPiXiHBTSpectraHandler::STARPiXiHBTSpectraHandler(BlastWaveFitter* aBWFitter, const char* aInputFileName)
  //:STARPiXiasHBTSpectraHandler(aBWFitter,aInputFileName)
  :IOHandler(aBWFitter,aInputFileName)
{
  mV2=0;
  mHBT=1;
  mSpectra=1;
}

void STARPiXiHBTSpectraHandler::initLoad(){
  // set Files
  setOutputFile("data/STARPiXiHBTSpectra.fit.root\0");
  //setOutputFile("data/STARPiXiHBTSpectra.fit.root\0");
  if(!mSpectraInputFile){
    mSpectraInputFile = new TFile("data/STARAuAuPP200GeVSpectra.root");
  }
  if(!mInputFile){
    mInputFile = new TFile("data/PiXiData.root");
    //mInputFile = new TFile("data/HBT/HBTparam_vsmt_5cent_TGraphs.root");
  }
  // Set centralities
  mNCent = 3;
  mNParticipant = new double[mNCent];
  mNParticipant[0] = 351.5; // 0-5 C0
  mNParticipant[1] = 299.0;  // 5-10 C1
  mNParticipant[2] = 234.5;  // 10-20 C2
  //mNParticipant[3] = 167.2; // 20-30 C3
  //mNParticipant[4] = 56.7;  // 30-80 C4-5-7-8
  //mNParticipant[4] = 96.5;  // 30-50 C4-5-7-8
  //mNParticipant[5] = 30.2;  // 50-80 C4-5-7-8

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
  // mParticleMass[6] = MPi;
  //mParticleMass[7] = MK;
  //mParticleMass[8] = MP;

}

void STARPiXiHBTSpectraHandler::loadSpectra(int iPart){
  //int tCentMerged[] = {0,1,2,3,4,6,9};
  int tCentMerged[] = {0,2,5,9};
  if(tCentMerged[mCent] == tCentMerged[mCent+1]-1){ // no merging
    char tGraphName[100];
    sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[iPart],mCent);
    mSpectras[iPart] = 
      mBWFitter->addSpectra((TGraphErrors*) mSpectraInputFile->Get(tGraphName),
			    mParticleMass[iPart],0.5,1.2);
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

void STARPiXiHBTSpectraHandler::load(){
  cout<<"kk check"<<endl;
  for(int ti=0; ti<6; ti++){
    loadSpectra(ti);
    mValidParticle[ti]=0;
  }
  char tNameo[100];
  // 0 is 0 to 5
  sprintf(tNameo,"GRoutC%i",mCent+1);
  //sprintf(tNameo,"pimROutC%i",mCent);
  char tNames[100];
  sprintf(tNames,"GRsideC%i",mCent+1);
  //sprintf(tNames,"pimRSideC%i",mCent);
  char tNamel[100];
  sprintf(tNamel,"GRlongC%i",mCent+1);
  //sprintf(tNamel,"pimRLongC%i",mCent);
  mBWFitter->addHbtRadii((TGraphErrors*) mInputFile->Get(tNameo),
			 (TGraphErrors*) mInputFile->Get(tNames),
			 (TGraphErrors*) mInputFile->Get(tNamel),
			 MPi);
  mValidParticle[0]=1;
}
			  

void STARPiXiHBTSpectraHandler::close(){
  //STARPiXiasHBTSpectraHandler::close();
  IOHandler::close();
  mSpectraInputFile->Close();
}

//void STARPiXiHBTSpectraHandler::fixParameters(){
//mBWFitter->freeAllParameters();
// mBWFitter->useR();
//mBWFitter->fixRhoA();
//}


