#include "BlastWaveFitter.h"

#include "TROOT.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH2.h"

#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <cstring>
using namespace std;

int main(int argc, char **argv){
  if(argc<3){
    cout << "Need at least 2 arguments" << endl;
    return 1;
  }
  int iDataType = atoi(argv[1]); // 2: 
  int iMinimize = atoi(argv[2]); // 0: none, 1: migrad, 2: migrad+minos, 
  int aContour = argc==4 ? atoi(argv[3]) : 0;

  double MPi = 0.1396; 
  double MK = 0.4937;
  double MKStar = 0.892;
  double MP = 0.9383;
  double MLa = 1.1157;
  double MXi = 1.32132;
  double MOm = 1.67245;
  double  MPhi = 1.019413;
  double MD0 = 1.8646;
  double MLac = 2.2849;


  TROOT ("myROOT","ttt");
  // _____________________________________________________________________
  // --- Instantiate Blast wave
  BlastWaveFitter* BWFitter =  BlastWaveFitter::instance();
  SpectraForFit* SpecPim =0;
  SpectraForFit* SpecKm =0;
  SpectraForFit* SpecP =0;
  SpectraForFit* SpecLa =0;
  SpectraForFit* SpecPip =0;
  SpectraForFit* SpecKp =0;
  SpectraForFit* SpecPB =0;
  SpectraForFit* SpecLaB =0;
  SpectraForFit* SpecK0s =0;
  SpectraForFit* SpecXi =0;
  SpectraForFit* SpecXiB =0;
  SpectraForFit* SpecOm =0;
   SpectraForFit* SpecPhi =0; 
   SpectraForFit* SpecKStar =0; 
  // _____________________________________________________________________
  // --- Get data and plug them in

  
  double T;      // Temperature GeV
  double Rho0;    // Flow rapidity
  double Rho2;    // Flow modulation (phi dependence)
  double Rho4;
  double Rx;      // in-plane radius (fm) 
  double Ry;      // out-of-plane radius (fm)
  double As;      // Smoothness of the source, 0= box -> 0.3 = smooth
  double Tau;      // System life time (fm/c)
  double Deltat; // Emission duration
  char fOutName[200];

  switch(iDataType){
  case 0:  // central    
    {
      TFile fIn("forFabrice.root");
      SpecXi = BWFitter->addSpectra((TGraphErrors*) fIn.Get("specXiAxiPtCen"), MXi);
      BWFitter->addV2((TGraphErrors*) fIn.Get("v2XiAxiPtCen"), MXi, 0.7, 3.);
      fIn.Close();
      strcpy(fOutName,"XiCent200GeV.root");
      break;
    }
  case 1:  // min bias
    {
      TFile fIn1("forFabriceMB.root");
      BWFitter->addV2((TGraphErrors*) fIn1.Get("v2XiAxiPtMB"), MXi, 0.7, 3.);
      fIn1.Close();
      TFile fIn2("dndptXiMB.root");
      cout << fIn2.Get("mtSpectraTotSet") << endl;
      SpecXi = BWFitter->addSpectra((TH1*) fIn2.Get("mtSpectraTotSet"), MXi);
      strcpy(fOutName,"XiMB200GeV.root");
      fIn2.Close();
      break;
    }
  case 2: // published data. Min bias only
    {
      TFile fIn("data/XiSpectraV2Y2.root");
      BWFitter->addV2((TGraphErrors*) fIn.Get("XiV2"), MXi, 0.7, 3.);
      SpecXi = BWFitter->addSpectra((TGraphErrors*) fIn.Get("XiMBSpectra"), MXi);
      SpecXiB = BWFitter->addSpectra((TGraphErrors*) fIn.Get("AxiMBSpectra"), MXi);
      fIn.Close();
      strcpy(fOutName,"data/XiSpectraV2Y2.fit.root");
      break;
    }
  case 3: // Y4 data from Markus. Published spectra. MB sum.
    {
      TFile fIn("data/v2/XiV2Y4.root");
      BWFitter->addV2((TGraphErrors*) fIn.Get("v2_cent0_80"), MXi, 0.7, 3.);
      fIn.Close();
      TFile fIn2("data/XiSpectraV2Y2.root");
      SpecXi = BWFitter->addSpectra((TGraphErrors*) fIn2.Get("XiMBSpectra"), MXi);
      SpecXiB = BWFitter->addSpectra((TGraphErrors*) fIn2.Get("AxiMBSpectra"), MXi);
      fIn2.Close();
      strcpy(fOutName,"data/XiSpectraV2Y4.fit.root");
      break;
    }
  }


  T = 0.14;// 0.111;
  Rho0 = 0.8; //0.818;
  Rho2 = 0.04; //0.058;
  Rho4 = 0.;
  Rx = 10.0;//8.3;
  Ry = 10.0;//9.6;      
  As = 0.;
  Tau = 6.5;//5.9;
  Deltat = 0.9;//0.0002;

   

  // _____________________________________________________________________
  // --- Fitting



  BWFitter->setParameters(T,Rho0,Rho2,Rho4,Rx,Ry,As,Tau,Deltat);
  BWFitter->fixAs();
  //BWFitter->fixT();
  //BWFitter->fixRho0();
  //BWFitter->fixRho2();
  BWFitter->fixRho4();
  //BWFitter->useR(); // merge Rx and Ry
  BWFitter->fixRx();
  //BWFitter->fixRy();
  BWFitter->fixTau();
  BWFitter->fixDeltat();
  if(iMinimize>0) 
    cout << "Ready to minimize. NDOF = " << BWFitter->nDOF() << endl;
  switch(iMinimize){
  case 0: // nothing
    break;
  case 1: // minimize
    cout << "MINIMIZE" << endl;
    BWFitter->minimize();
    break;
  case 2: // minimize and minos
    cout << "MINIMIZE + MINOS" << endl;
    BWFitter->minimize();
    BWFitter->calculatePreciseErrors();
    break;
  }
  BWFitter->coutFitResults();
    
  // _____________________________________________________________________
  // --- Output 
  TFile* fOut = new TFile(fOutName,"RECREATE");  
  if(aContour){
    cout << "Calculate Beta vs T contout" << endl;
    BWFitter->betaTContour(1)->Write(); 
    BWFitter->betaTContour(2)->Write();   
    BWFitter->betaTContour(3)->Write(); 
    BWFitter->S2VsRho2Contour(1)->Write();
    BWFitter->S2VsRho2Contour(2)->Write();
    BWFitter->S2VsRho2Contour(3)->Write();
  }
  // _____________________________________________________________________
  // --- Get data againwithout including them in fit

  

  fOut->cd();
  // Do it on the spectra to get the scaling (dN/dY) right
  if(SpecPim) 
    SpecPim->histo("BestPimPt", 100, BWFitter->blastWave())->Write();
  if(SpecPip) 
    SpecPip->histo("BestPipPt", 100, BWFitter->blastWave())->Write();
  if(SpecKm) SpecKm->histo("BestKmPt", 100, BWFitter->blastWave())->Write();
  if(SpecKp) SpecKp->histo("BestKpPt", 100, BWFitter->blastWave())->Write();
  if(SpecK0s) 
    SpecK0s->histo("BestK0sPt", 100, BWFitter->blastWave())->Write();
  if(SpecP) SpecP->histo("BestPPt", 100, BWFitter->blastWave())->Write();
  if(SpecPB) SpecPB->histo("BestPBPt", 100, BWFitter->blastWave())->Write();
  if(SpecLa) SpecLa->histo("BestLaPt", 100, BWFitter->blastWave())->Write();
  if(SpecLaB) 
    SpecLaB->histo("BestLaBPt", 100, BWFitter->blastWave())->Write();
  if(SpecXi) 
    SpecXi->histo("BestXiPt", 100, BWFitter->blastWave())->Write();
  if(SpecXiB) 
    SpecXiB->histo("BestXiBPt", 100, BWFitter->blastWave())->Write();
  if(SpecOm) 
    SpecOm->histo("BestOmPt", 100, BWFitter->blastWave())->Write();
  if(SpecPhi) 
    SpecPhi->histo("BestPhiPt", 100, BWFitter->blastWave())->Write();
  if(SpecKStar) 
    SpecKStar->histo("BestKStarPt", 100, BWFitter->blastWave())->Write();
    
  BWFitter->blastWave()->v2VsPt("Pi",MPi)->Write();
  BWFitter->blastWave()->v2VsPt("P",MP)->Write();
  BWFitter->blastWave()->v2VsPt("La",MLa)->Write();
  BWFitter->blastWave()->v2VsPt("Xi",MXi)->Write();
  BWFitter->blastWave()->v2VsPt("D0",MD0)->Write();
  BWFitter->blastWave()->v2VsPt("Lac",MLac)->Write();

  fOut->Close();

  //BWFitter->coutFitResults();
  
  delete BWFitter;


}
