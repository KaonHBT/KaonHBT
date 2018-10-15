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
  double MPi = 0.1396; 
  double MK = 0.4937;
  double MKStar = 0.892;
  double MP = 0.9383;
  double MLa = 1.1157;
  double MXi = 1.32132;
  double MOm = 1.67245;
  double  MPhi = 1.019413;
  if(argc!=3){
    cout << "Need 2 arguments" << endl;
    return 1;
  }
  int aWhatData = atoi(argv[1]); // 0: central 200, 1: mid-periph 200, 2: periph 200
  int aWhatMinimize = atoi(argv[2]); // 0: none, 1: migrad, 2: migrad+minos, 3: migrad spectra, 4: 

  TROOT ("myROOT","ttt");
  // _____________________________________________________________________
  // --- Instantiate Blast wave
  BlastWaveFitter* BWFitter =  BlastWaveFitter::instance();

  // _____________________________________________________________________
  // --- Get data and plug them in
  TFile fIn("hbtRPdata.root");
  
  double T=0.1;      // Temperature GeV
  double Rho0=0.9;    // Flow rapidity
  double RhoA=0.05;    // Flow modulation (phi dependence)
  double Rx=13;      // in-plane radius (fm) 
  double Ry=13;      // out-of-plane radius (fm)
  double As=0.1;
  double Tau=10;      // System life time (fm/c)
  double Deltat=1; // Emission duration
  char fOutName[50];

  switch(aWhatData){
	  case 0: // cent
  		BWFitter->addHbtRadiiSquare((TGraphErrors*) fIn.Get("CentROut_0"),
  													 (TGraphErrors*) fIn.Get("CentRSide_0"),
  													 (TGraphErrors*) fIn.Get("CentROutSide_0"),
  													 (TGraphErrors*) fIn.Get("CentRLong_0"),
  													 MPi);
  		BWFitter->addHbtRadiiSquareOsc((TGraphErrors*) fIn.Get("CentROut_2"),
  													 	  (TGraphErrors*) fIn.Get("CentRSide_2"),
  													 	  (TGraphErrors*) fIn.Get("CentROutSide_2"),
  													 	  (TGraphErrors*) fIn.Get("CentRLong_2"),
  													 	  MPi);  		
        strcpy(fOutName,"BWHbtRP0.root");
      break;
	  case 1: // mid-periph
  		BWFitter->addHbtRadiiSquare((TGraphErrors*) fIn.Get("MidROut_0"),
  													 (TGraphErrors*) fIn.Get("MidRSide_0"),
  													 (TGraphErrors*) fIn.Get("MidROutSide_0"),
  													 (TGraphErrors*) fIn.Get("MidRLong_0"),
  													 MPi);
  		BWFitter->addHbtRadiiSquareOsc((TGraphErrors*) fIn.Get("MidROut_2"),
  													 	  (TGraphErrors*) fIn.Get("MidRSide_2"),
  													 	  (TGraphErrors*) fIn.Get("MidROutSide_2"),
  													 	  (TGraphErrors*) fIn.Get("MidRLong_2"),
  													 	  MPi);  		
        strcpy(fOutName,"BWHbtRP1.root");
      break;
	  case 2: // periph
  		BWFitter->addHbtRadiiSquare((TGraphErrors*) fIn.Get("PerROut_0"),
  													 (TGraphErrors*) fIn.Get("PerRSide_0"),
  													 (TGraphErrors*) fIn.Get("PerROutSide_0"),
  													 (TGraphErrors*) fIn.Get("PerRLong_0"),
  													 MPi);
  		BWFitter->addHbtRadiiSquareOsc((TGraphErrors*) fIn.Get("PerROut_2"),
  													 	  (TGraphErrors*) fIn.Get("PerRSide_2"),
  													 	  (TGraphErrors*) fIn.Get("PerROutSide_2"),
  													 	  (TGraphErrors*) fIn.Get("PerRLong_2"),
  													 	  MPi);  		
        strcpy(fOutName,"BWHbtRP2.root");
      break;      
	}    

  // _____________________________________________________________________
  // --- Fitting

  BWFitter->setParameters(T,Rho0,RhoA,Rx,Ry,As,Tau,Deltat);
  //BWFitter->fixAs();
  //BWFitter->fixRho0();
  //BWFitter->fixRhoA();
  //BWFitter->useR(); // merge Rx and Ry
  //BWFitter->fixRx();
  //BWFitter->fixRy();
  //BWFitter->fixTau();
  //BWFitter->fixDeltat();
  cout << "Ready to minimize. NDOF = " << BWFitter->nDOF() << endl;
  switch(aWhatMinimize){
  case 0: // nothing
    break;
  case 1: // minimize
    BWFitter->minimize();
    break;
  case 2: // minimize and minos
    BWFitter->minimize();
    BWFitter->calculatePreciseErrors();
    break;
  case 3: // mininize only spectra
    BWFitter->fixRhoA();
    BWFitter->fixRx();
    BWFitter->fixRy();
    BWFitter->fixTau();
    BWFitter->fixDeltat();
    BWFitter->minimize();
    break;
   case 4: // minimize only time
   BWFitter->fixT();
   BWFitter->fixRho0();
    BWFitter->fixRhoA();
    BWFitter->fixRx();
    BWFitter->fixRy();   
    BWFitter->minimize();    
    BWFitter->calculatePreciseErrors();    
    break;
    case 5: // minimize only Hbt
   BWFitter->fixT();
   BWFitter->fixRho0();
    BWFitter->fixRhoA();    
	BWFitter->useR();
    BWFitter->minimize();        
    break;
  }
  BWFitter->coutFitResults();
    
  // _____________________________________________________________________
  // --- Output 
  TFile* fOut = new TFile(fOutName.str(),"RECREATE"); 
   		   //BWFitter->S2VsRhoAContour(1)->Write();
   		   //BWFitter->betaTContour(1)->Write(); 

    
    
  BWFitter->blastWave()->rOut("Pi",MPi)->Write();
  BWFitter->blastWave()->rSide("Pi",MPi)->Write();
  BWFitter->blastWave()->rLong("Pi",MPi)->Write();
  BWFitter->blastWave()->rOutVsMt("Pi",MPi)->Write();
  BWFitter->blastWave()->rSideVsMt("Pi",MPi)->Write();
  BWFitter->blastWave()->rLongVsMt("Pi",MPi)->Write();
  BWFitter->blastWave()->rOutSquare("Pi",MPi)->Write();
  BWFitter->blastWave()->rSideSquare("Pi",MPi)->Write();
  BWFitter->blastWave()->rLongSquare("Pi",MPi)->Write();
  BWFitter->blastWave()->rOutSideSquare("Pi",MPi)->Write();
  BWFitter->blastWave()->rOutSquareOsc("Pi",MPi)->Write();
  BWFitter->blastWave()->rSideSquareOsc("Pi",MPi)->Write();
  BWFitter->blastWave()->rLongSquareOsc("Pi",MPi)->Write();
  BWFitter->blastWave()->rOutSideSquareOsc("Pi",MPi)->Write();
  BWFitter->blastWave()->v2VsPt("Pi",MPi)->Write();
  BWFitter->blastWave()->v2VsPt("P",MP)->Write();
  BWFitter->blastWave()->v2VsPt("La",MLa)->Write();
  


  BWFitter->blastWave()->avSepOut("PiK",MPi,MK)->Write();
  BWFitter->blastWave()->avSepOutSpace("PiK",MPi,MK)->Write();
  BWFitter->blastWave()->avSepOutTime("PiK",MPi,MK)->Write();
  BWFitter->blastWave()->avSepSide("PiK",MPi,MK)->Write();  
  BWFitter->blastWave()->avSepLong("PiK",MPi,MK)->Write();
  BWFitter->blastWave()->avSepOut("PiP",MPi,MP)->Write();
  BWFitter->blastWave()->avSepOutSpace("PiP",MPi,MP)->Write();
  BWFitter->blastWave()->avSepOutTime("PiP",MPi,MP)->Write();
  BWFitter->blastWave()->avSepSide("PiP",MPi,MP)->Write();  
  BWFitter->blastWave()->avSepLong("PiP",MPi,MP)->Write();
  BWFitter->blastWave()->avSepOut("KP",MK,MP)->Write();
  BWFitter->blastWave()->avSepOutSpace("KP",MK,MP)->Write();
  BWFitter->blastWave()->avSepOutTime("KP",MK,MP)->Write();
  BWFitter->blastWave()->avSepSide("KP",MK,MP)->Write();  
  BWFitter->blastWave()->avSepLong("KP",MK,MP)->Write();

  BWFitter->blastWave()->rOutSquareVsPhiP("Pi200",MPi,0.2)->Write();
  BWFitter->blastWave()->rSideSquareVsPhiP("Pi200",MPi,0.2)->Write();
  BWFitter->blastWave()->rLongSquareVsPhiP("Pi200",MPi,0.2)->Write();
  BWFitter->blastWave()->rOutSideSquareVsPhiP("Pi200",MPi,0.2)->Write();

  BWFitter->blastWave()->rOutSquareVsPhiP("Pi300",MPi,0.3)->Write();
  BWFitter->blastWave()->rSideSquareVsPhiP("Pi300",MPi,0.3)->Write();
  BWFitter->blastWave()->rLongSquareVsPhiP("Pi300",MPi,0.3)->Write();
  BWFitter->blastWave()->rOutSideSquareVsPhiP("Pi300",MPi,0.3)->Write();

  BWFitter->blastWave()->rOutSquareVsPhiP("Pi400",MPi,0.4)->Write();
  BWFitter->blastWave()->rSideSquareVsPhiP("Pi400",MPi,0.4)->Write();
  BWFitter->blastWave()->rLongSquareVsPhiP("Pi400",MPi,0.4)->Write();
  BWFitter->blastWave()->rOutSideSquareVsPhiP("Pi400",MPi,0.4)->Write();  
  
  BWFitter->blastWave()->rOutSquareVsPhiP("Pi525",MPi,0.525)->Write();
  BWFitter->blastWave()->rSideSquareVsPhiP("Pi525",MPi,0.525)->Write();
  BWFitter->blastWave()->rLongSquareVsPhiP("Pi525",MPi,0.525)->Write();
  BWFitter->blastWave()->rOutSideSquareVsPhiP("Piv",MPi,0.525)->Write();

  BWFitter->blastWave()->v2VsPt("K",MK)->Write();
  BWFitter->blastWave()->rOut("K",MK)->Write();
  BWFitter->blastWave()->rSide("K",MK)->Write();
  BWFitter->blastWave()->rLong("K",MK)->Write();
  BWFitter->blastWave()->rOutVsMt("K",MK)->Write();
  BWFitter->blastWave()->rSideVsMt("K",MK)->Write();
  BWFitter->blastWave()->rLongVsMt("K",MK)->Write();

  fOut->Close();

  BWFitter->coutFitResults();
  
  delete BWFitter;


}
