#include "BlastWaveFitter.h"

#include "TROOT.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH2.h"

#include <cmath>
#include <iostream>
#include <cstdlib>
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
  int iSetOfPar = atoi(argv[1]); //
  int iDataType = atoi(argv[2]); // 0: none, 1: central130, 2: mid-perif130, 3: perif130, 4: cent200, 5: mid200, 6: periRho200
  

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
  TFile fIn("DataForFullBW.root");

  if(iDataType!=0){
  	if(iDataType<100){
	  	iDataType-=1;
    	switch(iDataType){
    	case 0:
      		SpecPim  = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130PimCent0to5"), 
					     MPi, 0.4, 1.2);
      		SpecPip  = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130PipCent0to5"), 
			    	 MPi, 0.4, 1.2);
      		SpecKm   = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130KmCent0to5"), 
			    	 MK, 0.3, 1.5);
      		SpecKp   = 
					BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130KpCent0to5"), 
			     		MK, 0.3, 1.5);
      		SpecPB = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130pbarCent0to5"), 
			    	 MP, 0.3, 2.);
      		SpecP    = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130pCent0to5"), 
			    	 MP, 0.3, 2.);
      		SpecLa    = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130LaCent0to5"), 
			    	 MLa, 0.3, 2.);
      		SpecLaB    = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130LaBCent0to5"), 
					     MLa, 0.3, 2.);
      		break;
    	case 1:
      		SpecPim  = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130PimCent15to30"), 
			    	 MPi, 0.4, 1.2);
      		SpecPip  = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130PipCent15to30"), 
			     	MPi, 0.4, 1.2);
      		SpecKm   = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130KmCent15to30"), 
			    	 MK, 0.3, 1.5);
      		SpecKp   = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130KpCent15to30"), 
			    	 MK, 0.3, 1.5);
      		SpecPB = 
				BWFitter->addSpectra((TGraphErrors*)fIn.Get("Phenix130pbarCent15to30"), 
			    	 MP, 0.3, 2.);
      		SpecP = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130pCent15to30"), 
			    	 MP, 0.3, 2.);
      		SpecLa = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130LaCent20to35"), 
			    	 MLa, 0.3, 2.);
      		SpecLaB = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130LaBCent20to35"), 
			    	 MLa, 0.3, 2.);
      		break;
    	case 2:
      		SpecPim  = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130PimCent60to92"), 
			    	 MPi, 0.4, 1.2);
      		SpecPip  = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130PipCent60to92"), 
			    	 MPi, 0.4, 1.2);
      		SpecKm   = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130KmCent60to92"), 
			    	 MK, 0.3, 1.5);
      		SpecKp   = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130KpCent60to92"), 
			    	 MK, 0.3, 1.5);
      		SpecPB =
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130pbarCent60to92"), 
			    	 MP, 0.3, 2.);
      		SpecP = 
					BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130pCent60to92"), 
			     		MP, 0.3, 2.);
      		SpecLa = 
					BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130LaCent35to75"), 
			     		MLa, 0.3, 2.);
      		SpecLaB = 
					BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130LaBCent35to75"), 
			     		MLa, 0.3, 2.);
 
      		break;
    	case 3: // preliminary central AuAu 200
      		SpecK0s = BWFitter->addSpectra((TGraphErrors*) fIn.Get("PaulsSpecK0s0to5")
				     , MK, 0.5, 2.);
    		  SpecLa = BWFitter->addSpectra((TGraphErrors*) fIn.Get("HuisSpecLa0to5"),
				    MLa, 0., 2.5);
      		SpecPim = BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaisPimC0"), 
				     MPi);
      		SpecKm = BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaisKmC0"), MK);
      		SpecP = BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaispC0"), MP);
      		SpecPB = BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaispBC0"), MP);
      		break;
    	case 4: // preliminary MidPerif AuAu 200  
      		SpecPim = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaisPimC2"), MPi);
      		SpecKm = BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaisKmC2"), MK);
      		SpecP = BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaispC2"), MP);
      		SpecPB = BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaispBC2"), MP);  
      		break;
    	case 5: // preliminary Perif AuAu 200    
        	SpecPim = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaisPimC6"), MPi);
      		SpecKm = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaisKmC6"), MK);
      		SpecP = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaispC6"), MP);
      		SpecPB = 
      			BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaispBC6"), MP);
      		SpecK0s = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("PaulsSpecK0s40to60"), 
			    	 MK, 0.5, 2.);
      		SpecLa = 
				BWFitter->addSpectra((TGraphErrors*) fIn.Get("HuisSpecLa40to60"), 
			    	 MLa, 0., 2.5);
      		break;
    	case 6:
      		SpecPim = 
         		BWFitter->addSpectra((TGraphErrors*) fIn.Get("MasPim"), MPi,0.3,2.);
      		SpecKm = 
         		BWFitter->addSpectra((TGraphErrors*) fIn.Get("MasKm"), MK);
      		SpecKp = 
         		BWFitter->addSpectra((TGraphErrors*) fIn.Get("MasKp"), MK);
     		SpecP = 
       			BWFitter->addSpectra((TGraphErrors*) fIn.Get("MasP"), MP);
      		SpecPB = 
        		BWFitter->addSpectra((TGraphErrors*) fIn.Get("MasPB"), MP);
      		SpecLa = 
         		BWFitter->addSpectra((TGraphErrors*) fIn.Get("MasLa"), MLa);
      		SpecLaB = 
         		BWFitter->addSpectra((TGraphErrors*) fIn.Get("MasLaB"), MLa);         
      		break;    
    	case 7:
      		SpecXi = 
         		BWFitter->addSpectra((TGraphErrors*) fIn.Get("Xi130C0to10"), MXi);      
      		SpecXiB = 
         		BWFitter->addSpectra((TGraphErrors*) fIn.Get("AXi130C0to10"), MXi);                   
      		SpecOm = 
         		BWFitter->addSpectra((TGraphErrors*) fIn.Get("Om130C0to10"), MOm);      
      		SpecPhi = 
         		BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phi130Cent0to11"), MPhi);              
      		break;    
    	}
  	}
  	else{
	  int iCent = iDataType- ((int) floor(iDataType/100.))*100;
	  char tPimName[50];
	  sprintf(tPimName,"KaisPimC%i",iCent);
	  SpecPim = BWFitter->addSpectra((TGraphErrors*) fIn.Get(tPimName),
					 MPi);
	  char tKmName[50];
	  sprintf(tKmName,"KaisKmC%i",iCent);
	  SpecKm = BWFitter->addSpectra((TGraphErrors*) fIn.Get(tKmName), 
					MK);
	  char tpName[50];
	  sprintf(tpName,"KaispC%i",iCent);
	  SpecP = BWFitter->addSpectra((TGraphErrors*) fIn.Get(tpName), MP);
	  char tpBName[50];
	  sprintf(tpBName,"KaispBC%i",iCent);
	  SpecPB = BWFitter->addSpectra((TGraphErrors*) fIn.Get(tpBName), MP);
  		}
	}
	fIn.Close();
	
  // _____________________________________________________________________
  // --- Set Parameters
  double T;      // Temperature GeV
  double Rho0;    // Flow rapidity
  double Rho2;    // Flow modulation (phi dependence)
  double Rho4;
  double Rx;      // in-plane radius (fm) 
  double Ry;      // out-of-plane radius (fm)
  double As;
  double Tau;      // System life time (fm/c)
  double Deltat; // Emission duration
  char fOutName[200];
  sprintf(fOutName,"Paper/Calc/BWCalc%i.root",iSetOfPar);
  T = 0.1;
  Rho0 = 0.9; 
  Rho2 = 0.; 
  Rho4 = 0.;
  Rx = 13.;//12.04;
  Ry = 13.;//12.04;
  As = 0.;
  Tau = 9.;
  Deltat = 2.;
  int NStat=0;
  switch(iSetOfPar){
  case 0:
    break;
  case 1:
    T = 0.108;// 0.111;
    Rho0 = 0.88; //0.818;
    Rho2 = 0.06; //0.058;
    Rx = 12.9;//8.3;
    Ry = 12.8;//9.6;
    Tau = 8.9;//5.9;
    Deltat = 0.;//0.0002;
    break;
  case 2:
    T = 0.106;// 0.111;
    Rho0 = 0.87; //0.818;
    Rho2 = 0.052; //0.058;
    Rx = 10.2;//8.3;
    Ry = 11.8;//9.6;
    Tau = 7.4;//5.9;
    Deltat = 0.8;//0.0002;
    break;		
  case 10:		
    T =0.08;
    break;
  case 11:
    T=0.1;
    break;
  case 12:
    T=0.12;
    break;
  case 13:
    T=0.14;
    break;
  case 14:
    T = 0.01;
    break;
  case 15:
    T = 1.;
    break;	  	
  case 20:
    Rho0 = 0.;
    break;
  case 21:
    Rho0 = 0.6;
    break;
  case 22:
    Rho0 = 0.9;	  		  		  	
    break;
  case 23:
    Rho0 = 1.2;	  		  		  	
    break;
  case 24:
    Rho0 = 0.01;	  		  		  	
    break;	  	
  case 25:
    Rho0 = 2.;	  		  		  	
    break;	  		  	
  case 30:
    Rx = 7.;
    Ry = 7.;
    break;
  case 31:
    Rx = 10.;
    Ry = 10.;
    break;
  case 32:
    Rx = 13.;
    Ry = 13.;
    break;
  case 33:
    Rx = 16.;
    Ry = 16.;
    break;	    
  case 40:
    Tau = 7.;
    break;
  case 41:
    Tau = 9.;
    break;
  case 42:
    Tau = 11.;	  	
    break;	 
  case 43:
    Tau = 13.;	  	
    break;		  	
  case 50:
    Deltat=0.;
    break;
  case 51:
    Deltat=3.;
    break; 	
  case 52:
    Deltat=6.;
    break; 	
  case 60:
    Rho2=0.02;
    break;
  case 61:
    Rho2=0.05;
    break;
  case 62:
    Rho2=0.1;
    break;	  		  	
  case 70:
    Rx=11.;
    break;
  case 71:
    Rx=13.;
    break;
  case 72:
    Rx=16.;
    break;
  case 80:
    Ry = 13.;
    Rx = 11.;
    Rho2 = 0.05;
    break;
  case 81:
    Ry = 13.;
    Rx = 11.;
    break;	
  case 82:
    Rho2 = 0.05;	  	  	
    break;
  case 90:
    As = 0.; 
    break;
  case 91:
    As = 0.1;
    Rho0 = 0.845; 
    break;
  case 92:
    As = 0.2;
    Rho0 = 0.732; 
    break;
  case 93:
    As = 0.3;
    Rho0 = 0.612; 
    break;
  case 94:
    As = 0.01;
    Rho0 = 0.9; 
    break;
  case 1000: 
    T = 0.1;
    Rho0 = 1.; 
    Rho2 = 0.05; 
    Rx = 13.;
    Ry = 13.;
    Tau = 10.;
    Deltat = 2.;	  
    break;
  case 1001: 
    T = 0.1;
    Rho0 = 1.; 
    Rho2 = 0.05; 
    Rx = 13.;
    Ry = 13.;
    Tau = 10.;
    Deltat = 2.;	  
    NStat=1; // bose-einstein    
    break;
  case 2000:
    T = 0.099;
    Rho0 = 0.97;
    Rho2 = -0.018;
    Rx = 10.9;
    Ry = 11.8;
    Tau = 7.6;
    Deltat = 2.5;
    break;
  }

  BWFitter->setParameters(T,Rho0,Rho2,Rho4,Rx,Ry,As,Tau,Deltat);

  
  // _____________________________________________________________________
  // --- Output 
  TFile* fOut = new TFile(fOutName,"RECREATE");  

  cout << "Estruct " << endl;
  BWFitter->blastWave()->meanPtVsPhi("Pi",MPi)->Write();
  BWFitter->blastWave()->autoCorrVsPhi("Pi",MPi)->Write();

  /*  BWFitter->blastWave()->setBoseEinsteinStatistics();
  BWFitter->blastWave()->spectra("pi",MPi)->Write();
  // Do it on the spectra to get the scaling (dN/dY) right
  BWFitter->blastWave()->avSepOutVsPhiP("PiK300",MPi,MK, 0.3)->Write();  
  BWFitter->blastWave()->avSepSideVsPhiP("PiK300",MPi,MK, 0.3)->Write();   
  BWFitter->blastWave()->avSepOutVsPhiP("PiK600",MPi,MK, 0.6)->Write();  
  BWFitter->blastWave()->avSepSideVsPhiP("PiK600",MPi,MK, 0.6)->Write();    
  BWFitter->blastWave()->avSepOutVsPhiP("PiK900",MPi,MK, 0.9)->Write();  
  BWFitter->blastWave()->avSepSideVsPhiP("PiK900",MPi,MK, 0.9)->Write();    
  BWFitter->blastWave()->avSepOutVsPhiP("PiK1200",MPi,MK, 1.2)->Write();  
  BWFitter->blastWave()->avSepSideVsPhiP("PiK1200",MPi,MK, 1.2)->Write();     

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
  BWFitter->blastWave()->rOutSideSquareVsPhiP("Pi525",MPi,0.525)->Write();  
  
  //if(iSetOfPar==1000){
  //	fOut->Close();
  //	delete BWFitter;	    
  //    return 0;
  //}  


  	  
  if(SpecPim) 
    SpecPim->histo("BestPimPt", 100, BWFitter->blastWave())->Write();
  if(SpecPip) 
    SpecPip->histo("BestPipPt", 100, BWFitter->blastWave())->Write();
  if(SpecKm) SpecKm->histo("BestKmPt", 100, BWFitter->blastWave())->Write();
  if(SpecKp) SpecKp->histo("BestKpPt", 100, BWFitter->blastWave())->Write();
  if(SpecK0s) 
    SpecK0s->histo("BestK0sPt", 100, BWFitter->blastWave())->Write();
  if(SpecKStar) 
    SpecKStar->histo("BestKStarPt", 100, BWFitter->blastWave())->Write();
        
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
  BWFitter->blastWave()->v4VsPt("Pi",MPi)->Write();
  BWFitter->blastWave()->v6VsPt("Pi",MPi)->Write();
  BWFitter->blastWave()->v8VsPt("Pi",MPi)->Write();
  BWFitter->blastWave()->avOut("Pi",MPi)->Write();
  BWFitter->blastWave()->avTime("Pi",MPi)->Write();
  BWFitter->blastWave()->avSide("Pi",MPi)->Write();  
    
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
  
  BWFitter->blastWave()->rOutSquareVsPhiP("Pi600",MPi,0.6)->Write();
  BWFitter->blastWave()->rSideSquareVsPhiP("Pi600",MPi,0.6)->Write();
  BWFitter->blastWave()->rLongSquareVsPhiP("Pi600",MPi,0.6)->Write();
  BWFitter->blastWave()->rOutSideSquareVsPhiP("Pi600",MPi,0.6)->Write();

  BWFitter->blastWave()->v2VsPt("K",MK)->Write();
  BWFitter->blastWave()->rOut("K",MK)->Write();
  BWFitter->blastWave()->rSide("K",MK)->Write();
  BWFitter->blastWave()->rLong("K",MK)->Write();
  BWFitter->blastWave()->rOutVsMt("K",MK)->Write();
  BWFitter->blastWave()->rSideVsMt("K",MK)->Write();
  BWFitter->blastWave()->rLongVsMt("K",MK)->Write();
  BWFitter->blastWave()->avOut("K",MK)->Write();
  BWFitter->blastWave()->avTime("K",MK)->Write();
  BWFitter->blastWave()->avSide("K",MK)->Write();      
     
  cout << "emMomDist" << endl;
  BWFitter->blastWave()->emMomDist(MPi)->Write();
  cout << "emPosDist" << endl;
  BWFitter->blastWave()->emPosDist(MPi, 0.25, 0.)->Write();
  BWFitter->blastWave()->emPosDist(MK, 0.25*MK/MPi, 0.)->Write();   
  BWFitter->blastWave()->emPosDist(MPi, 0.5, 0.)->Write();
  BWFitter->blastWave()->emPosDist(MK, 0.5*MK/MPi, 0.)->Write();   
  BWFitter->blastWave()->emPosDist(MPi, 0.25, 0.)->Write();
  BWFitter->blastWave()->emPosDist(MPi, 0.25, TMath::Pi()/4.)->Write();   
  BWFitter->blastWave()->emPosDist(MPi, 0.25, TMath::Pi()/2.)->Write();
  BWFitter->blastWave()->emPosDist(MPi, 0.5, 0.)->Write();
  BWFitter->blastWave()->emPosDist(MPi, 0.5, TMath::Pi()/4.)->Write();   
  BWFitter->blastWave()->emPosDist(MPi, 0.5, TMath::Pi()/2.)->Write();
  cout << "ZDist" << endl;
  BWFitter->blastWave()->ZDist(MPi, 0.25)->Write();
  BWFitter->blastWave()->ZDist(MPi, 0.5)->Write(); 

  // ---
  BWFitter->blastWave()->setFermiDiracStatistics();
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
  BWFitter->blastWave()->v2VsPt("P",MP)->Write();
  BWFitter->blastWave()->v4VsPt("P",MP)->Write();
  BWFitter->blastWave()->v6VsPt("P",MP)->Write();
  BWFitter->blastWave()->v8VsPt("P",MP)->Write();
  BWFitter->blastWave()->v2VsPt("La",MLa)->Write();
  BWFitter->blastWave()->emPosDist(MP, 0.25*MP/MPi, 0.)->Write();
  BWFitter->blastWave()->emPosDist(MP, 0.5*MP/MPi, 0.)->Write();
  BWFitter->blastWave()->avOut("P",MP)->Write();
  BWFitter->blastWave()->avTime("P",MP)->Write();
  BWFitter->blastWave()->avSide("P",MP)->Write();  

  BWFitter->blastWave()->setBoltzmanStatistics();
  BWFitter->blastWave()->v2VsPt("Xi",MXi)->Write();
  BWFitter->blastWave()->rOut("Xi",MXi)->Write();
  BWFitter->blastWave()->rSide("Xi",MXi)->Write();
  BWFitter->blastWave()->rLong("Xi",MXi)->Write();
  BWFitter->blastWave()->rOutVsMt("Xi",MXi)->Write();
  BWFitter->blastWave()->rSideVsMt("Xi",MXi)->Write();
  BWFitter->blastWave()->rLongVsMt("Xi",MXi)->Write();
  BWFitter->blastWave()->avOut("Xi",MXi)->Write();
  BWFitter->blastWave()->avTime("Xi",MXi)->Write();
  BWFitter->blastWave()->avSide("Xi",MXi)->Write();   
  BWFitter->blastWave()->avSepOut("PiXi",MPi,MXi)->Write();
  BWFitter->blastWave()->avSepOutSpace("PiXi",MPi,MXi)->Write();
  BWFitter->blastWave()->avSepOutTime("PiXi",MPi,MXi)->Write();
  BWFitter->blastWave()->avSepSide("PiXi",MPi,MXi)->Write();  
  BWFitter->blastWave()->avSepLong("PiXi",MPi,MXi)->Write();
  */  
  fOut->Close();


  delete BWFitter;

	return 0;
}
