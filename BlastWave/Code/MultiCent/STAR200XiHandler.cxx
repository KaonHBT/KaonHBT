#include "STAR200XiHandler.h"
#include "BlastWaveFitter.h"
#include "DataForFit.h"
#include "TFile.h"
#include <cstdlib>
using namespace std;
STAR200XiHandler::STAR200XiHandler(BlastWaveFitter* aBWFitter, const char* aInputFileName)
  :IOHandler(aBWFitter,aInputFileName)
{
  mSpectra =1;
  mHBT=0;
  mV2=1;
}

void STAR200XiHandler::initLoad(){
  // set Files
  setOutputFile("data/STARAuAu200Xi.fit.root\0");
  if(!mInputFile){
    mInputFile = new TFile("data/STARAuAu200Xi.root");
  }
  // Set centralities
  mNCent = 2;
  mNParticipant = new double[mNCent];
  //double NParticipantSTARTOfdAuPP200Spectra[] = {
  //15.1211, 11.3394, 4.92598, 8.2, 
  //1.};
  mNParticipant[0] = 300.;
  mNParticipant[1] = 100.;
   
  // Set particle infos
  mNParticles = 1;
  mParticleName = new char*[mNParticles];
  char tParticleName[10][6] = {"Xi"};
  mParticleMass = new double[mNParticles];
  mParticleMass[0] = MXi;
  for(int ti=0; ti<mNParticles; ti++){
    mParticleName[ti] = new char[10];
    strcpy(mParticleName[ti],tParticleName[ti]);
  }
  // Set blast wave fitter parameters
  mBWFitter->fixAs();
  mBWFitter->fixRy();
  mBWFitter->fixTau();
  mBWFitter->fixDeltat();
}

void STAR200XiHandler::load(){
  char tCentName[2][10] = {"Cen","MB"};
  for(int ti=0; ti<mNParticles; ti++){
    char tGraphName1[100];
    sprintf(tGraphName1,"Spec%s%s",mParticleName[ti],tCentName[mCent]);
    //cout << tGraphName1 << endl;
    mSpectras[ti] = 
      mBWFitter->addSpectra((TGraphErrors*) mInputFile->Get(tGraphName1),
			    mParticleMass[ti], 0., 2.5);
    char tGraphName2[100];
    sprintf(tGraphName2,"V2%s%s",mParticleName[ti],tCentName[mCent]);
    //cout << tGraphName2 << endl;
    mBWFitter->addV2((TGraphErrors*) mInputFile->Get(tGraphName2),
		     mParticleMass[ti], 0.6, 2.5);

  }
}
			  
