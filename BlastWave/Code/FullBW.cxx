// g++ StRoot/macros/FullBW.cxx -o FullBW.exe -I$ROOTSYS/include -L${ROOTSYS}/lib -lNew -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMinuit -lm -ldl -rdynamic 
#include "TROOT.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TH2.h"
#include "TMath.h"

#include <cmath>
#include <iomanip>
using namespace std;

double fBW(double aMass, 
	   double aX, double aY, double aEta,
	   double aPt, double aPhiP,
	   double aT,
	   double aRho0, double aRhoA,
	   double aRx, double aRy){
  double tR = sqrt(aX*aX/aRx/aRx+aY*aY/aRy/aRy);
  if(tR<1){
    //double tPhiS = (aX!=0)? TMath::ATan2(aY,aX) : 
    //  (aY>0? TMath::Pi()/2. : - TMath::Pi()/2.);
    double tPhiB = TMath::ATan2(aY*aRx*aRx/aRy/aRy,aX);//+
      //(tPhiS>0? (tPhiS>TMath::Pi()/2. ? TMath::Pi() : 0.) :
      //        (tPhiS<-TMath::Pi()/2. ? -TMath::Pi() : 0.));
    double tRho = tR*(aRho0+aRhoA*cos(2*tPhiB));
    double tCoshEta = TMath::CosH(aEta);
    return exp(-sqrt(aPt*aPt+aMass*aMass)/aT*TMath::CosH(tRho)*tCoshEta)*
      tCoshEta*exp(aPt/aT*TMath::SinH(tRho)*cos(aPhiP-tPhiB));    
  }
  else{
    return 0.;
  }
}

void fBWSpaceTime(double* aResult, 
		  double aMass, 
		  double aPt, double aPhiP,
		  double aT,
		  double aRho0, double aRhoA,
		  double aRx, double aRy, 
		  double aTau, double aDeltat){
  int NBinEta = 100;
  double tMinEta = -5.;
  double tMaxEta = 5.;
  double tStepEta = (tMaxEta-tMinEta)/NBinEta;
  tMinEta += tStepEta/2.;
  int NBinX = 20;
  double tMinX = -1.*aRx;
  double tMaxX = aRx;
  double tStepX = (tMaxX-tMinX)/NBinX;
  tMinX += tStepX/2.;
  int NBinY = 20;
  double tMinY = -1.*aRy;
  double tMaxY = aRy;
  double tStepY = (tMaxY-tMinY)/NBinY;
  tMinY += tStepY/2.;
 
  double tMt = sqrt(aPt*aPt+aMass*aMass);
  double tSum = 0.;
  double tX = 0.;
  double tX2 = 0.;
  double tY = 0.;
  double tY2 = 0.;
  double tXY = 0.;
  double tZ2 = 0.;
  double tT = 0.;
  double tT2 = 0.;
  double tXT = 0.;
  double tYT = 0.;
  double tR, tPhiB;
  double tRho, tBeta, tK0Beta, tK1Beta, tK2Beta, tBetaWeight, tAlphaWeight;
  //  for(double tEta = tMinEta; tEta < tMaxEta; tEta += tStepEta ){
  //double tCoshEta = TMath::CosH(tEta);
  //double tSinhEta = TMath::SinH(tEta);
    
  double tH1OverH0 = (aDeltat*aDeltat+aTau*aTau)/aTau;
  double tH2OverH0 = (3*aDeltat*aDeltat+aTau*aTau);

  for(double tLX =tMinX; tLX < tMaxX; tLX += tStepX){
    for(double tLY =tMinY; tLY < tMaxY; tLY += tStepY){
      tR = sqrt(tLX*tLX/aRx/aRx+tLY*tLY/aRy/aRy);

      if(tR<1.){
	tPhiB = TMath::ATan2(tLY*aRx*aRx/aRy/aRy,tLX);
	//  (tPhiS>0? (tPhiS>TMath::Pi()/2. ? TMath::Pi() : 0.) :
	//   (tPhiS<-TMath::Pi()/2. ? -TMath::Pi() : 0.));
	tRho = tR*(aRho0+aRhoA*cos(2*tPhiB));	
	tBeta = tMt/aT*TMath::CosH(tRho);
	tK0Beta = TMath::BesselK0(tBeta);
	tK1Beta = TMath::BesselK1(tBeta);
	tK2Beta = 2./tBeta*tK1Beta+tK0Beta;
	tBetaWeight = 2*tK1Beta;
	tAlphaWeight = exp(aPt/aT*TMath::SinH(tRho)*cos(aPhiP-tPhiB));
	  
	tSum += (tAlphaWeight*tBetaWeight);
	tX   += (tAlphaWeight*tBetaWeight * tLX);
	tX2  += (tAlphaWeight*tBetaWeight * tLX * tLX);
	tY   += (tAlphaWeight*tBetaWeight * tLY);
	tY2  += (tAlphaWeight*tBetaWeight * tLY * tLY);
	tXY  += (tAlphaWeight*tBetaWeight * tLX * tLY);
 	tBetaWeight = 2.*(tK1Beta/tBeta+tK0Beta);
	tT   += (tAlphaWeight*tBetaWeight * tH1OverH0);
	tXT  += (tAlphaWeight*tBetaWeight * tLX * tH1OverH0);
	tYT  += (tAlphaWeight*tBetaWeight * tLY * tH1OverH0);
	tBetaWeight = 2./tBeta*tK2Beta;
	tZ2  += (tAlphaWeight*tBetaWeight * tH2OverH0); 
	tBetaWeight += (2*tK1Beta);
	tT2  += (tAlphaWeight*tBetaWeight * tH2OverH0);	
      }
    }
  }
  tX  /= tSum;
  tX2 /= tSum; tX2 -= (tX*tX);
  tY  /= tSum;
  tY2 /= tSum; tY2 -= (tY*tY);
  tXY /= tSum; tXY -= (tX*tY);
  tZ2 /= tSum;
  tT  /= tSum;
  tT2 /= tSum; tT2 -= (tT*tT);
  tXT /= tSum; tXT -= (tX*tT);
  tYT /= tSum; tYT -= (tY*tT);

  double tBetaT = aPt/tMt;
  double tGammaT = tMt/aMass;
  double tCosPhiP = cos(aPhiP);
  double tSinPhiP = sin(aPhiP);
  // <Rout> No Time
  aResult[0] = tX * cos(aPhiP) + tY * sin(aPhiP);
  // <Rout> Time
  aResult[1] = aResult[0] - tBetaT * tT;
  // <Rout*>
  aResult[2] = aResult[1] * tGammaT; 
  // Rout^2
  aResult[3] = 0.5*(tX2+tY2)+0.5*(tX2-tY2)*cos(2*aPhiP)+tXY*sin(2*aPhiP)
    - 2*tBetaT*(tXT*cos(aPhiP)+tYT*sin(aPhiP))+tBetaT*tBetaT*tT2;
  // <Rside> No Time
  aResult[4] = tY * cos(aPhiP) - tX * sin(aPhiP);
  // Rside^2
  aResult[5] = 0.5*(tX2+tY2)-0.5*(tX2-tY2)*cos(2*aPhiP)-tXY*sin(2*aPhiP);
  // Routside^2
  aResult[6] = tXY*cos(2*aPhiP)-0.5*(tX2-tY2)*sin(2*aPhiP)
    + tBetaT*(tXT*sin(aPhiP)-tYT*cos(aPhiP));
  // <RLong>
  aResult[7] = 0;
  // RLong^2
  aResult[8] = tZ2;
  // Weight for v2 calculation or pt/mt spectra.
  aResult[9] = tSum;
}


	   
	   
int main(){
  double MPi = 0.139; 
  double MK = 0.495;
  double MP = 0.938;
  double T = 0.11; // GeV
  double Rho0 = 0.9; // Flow rapidity
  double RhoA = 0.; // Flow modulation (phi dependence)
  double Rx = 13.;//1.35; // Box Radius (fm) 
  double Ry = 13.; // ellipsoid box
  double Tau = 10.; // System life time (fm/c)
  double Deltat = 2.; // Emission duration

  TROOT ("myROOT","ttt");
  TFile outFile("FullBW.root","RECREATE");
  int N = 100;
  TH2D HTest("HTest","HTest",N,-Rx,Rx,N,-Ry,Ry);
  for(int ti=1; ti<= N; ti++){
    double tX = HTest.GetXaxis()->GetBinCenter(ti);
    for(int tj=1; tj<=N; tj++){
      double tY = HTest.GetYaxis()->GetBinCenter(tj);
      HTest.SetCellContent(ti,tj,fBW(MPi, tX, tY, 0., 0.15, 0.,
				  T, Rho0, RhoA, Rx, Ry));
      
    }
  }
  HTest.Write();

  double Result[10];
  fBWSpaceTime(Result,MPi, 0.15, 0., T, Rho0, RhoA, Rx,Ry, Tau, Deltat);
  for(int ti=0; ti<10; ti++) cout << Result[ti] << " ";
  cout << endl;

  int NBinKt = 20;
  int NBinPhi = 10;
  TH2D HRO2Pi("HRO2Pi","HROPi",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HRS2Pi("HRS2Pi","HRS2Pi",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HROS2Pi("HROS2Pi","HROS2Pi",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HRL2Pi("HRL2Pi","HRL2Pi",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HROPi("HROPi","HROPi",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HRSPi("HRSPi","HRSPi",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH1D HV2Pi("HV2Pi","HV2Pi",NBinKt,0.05,2.05);
  TH1D HPtPi("HPtPi","HPtPi",NBinKt,0.05,2.05);
  TH2D HRO2K("HRO2K","HROK",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HRS2K("HRS2K","HRS2K",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HROS2K("HROS2K","HROS2K",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HRL2K("HRL2K","HRL2K",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HROK("HROK","HROK",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HRSK("HRSK","HRSK",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH1D HV2K("HV2K","HV2K",NBinKt,0.05,2.05);
  TH1D HPtK("HPtK","HPtK",NBinKt,0.05,2.05);
  TH2D HRO2P("HRO2P","HROP",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HRS2P("HRS2P","HRS2P",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HROS2P("HROS2P","HROS2P",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HRL2P("HRL2P","HRL2P",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HROP("HROP","HROP",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HRSP("HRSP","HRSP",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH1D HV2P("HV2P","HV2P",NBinKt,0.05,2.05);
  TH1D HPtP("HPtP","HPtP",NBinKt,0.05,2.05);

  TH2D HROPiK("HROPiK","HROPiK",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HRSPiK("HRSPiK","HRSPiK",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HROPiP("HROPiP","HROPiP",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HRSPiP("HRSPiP","HRSPiP",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HROKP("HROKP","HROKP",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());
  TH2D HRSKP("HRSKP","HRSKP",NBinKt,0.05,2.05,NBinPhi,0.,TMath::Pi());

  for(int ti=1; ti <= NBinKt ; ti++){
   
    double v2Pi = 0.;
    double yieldPi = 0.;
    double v2K = 0.;
    double yieldK = 0.;
    double v2P = 0.;
    double yieldP = 0.;
    double tPt = HRO2Pi.GetXaxis()->GetBinCenter(ti);

    cout << ti << " " << tPt << endl;
    for(int tj=1; tj <= NBinPhi ; tj++){
      double tPhi = HRO2Pi.GetYaxis()->GetBinLowEdge(tj);
      fBWSpaceTime(Result,MPi, tPt, tPhi,
		   T, Rho0, RhoA, Rx, Ry, Tau, Deltat);
      v2Pi += Result[9] * cos(2*tPhi);
      yieldPi += Result[9];
      HRO2Pi.SetCellContent(ti,tj,Result[3]);
      HRS2Pi.SetCellContent(ti,tj,Result[5]);
      HROS2Pi.SetCellContent(ti,tj,Result[6]);
      HRL2Pi.SetCellContent(ti,tj,Result[8]);
      HROPi.SetCellContent(ti,tj,Result[2]);
      HRSPi.SetCellContent(ti,tj,Result[4]);
      double tROPi = Result[2];
      double tRSPi = Result[4];

      fBWSpaceTime(Result,MK, tPt, tPhi,
		   T, Rho0, RhoA, Rx, Ry, Tau, Deltat);
      v2K += Result[9] * cos(2*tPhi);
      yieldK += Result[9];
      HRO2K.SetCellContent(ti,tj,Result[3]);
      HRS2K.SetCellContent(ti,tj,Result[5]);
      HROS2K.SetCellContent(ti,tj,Result[6]);
      HRL2K.SetCellContent(ti,tj,Result[8]);
      HROK.SetCellContent(ti,tj,Result[2]);
      HRSK.SetCellContent(ti,tj,Result[4]);
      double tROK = Result[2];
      double tRSK = Result[4];

      fBWSpaceTime(Result,MP, tPt, tPhi,
		   T, Rho0, RhoA, Rx, Ry, Tau, Deltat);
      v2P += Result[9] * cos(2*tPhi);
      yieldP += Result[9];
      HRO2P.SetCellContent(ti,tj,Result[3]);
      HRS2P.SetCellContent(ti,tj,Result[5]);
      HROS2P.SetCellContent(ti,tj,Result[6]);
      HRL2P.SetCellContent(ti,tj,Result[8]);
      HROP.SetCellContent(ti,tj,Result[2]);
      HRSP.SetCellContent(ti,tj,Result[4]);


      fBWSpaceTime(Result,MK, MK/MPi*tPt, tPhi,
		   T, Rho0, RhoA, Rx, Ry, Tau, Deltat);
      HROPiK.SetCellContent(ti,tj,tROPi-Result[2]);
      HRSPiK.SetCellContent(ti,tj,tRSPi-Result[4]);
      fBWSpaceTime(Result,MP, MP/MPi*tPt, tPhi,
		   T, Rho0, RhoA, Rx, Ry, Tau, Deltat);
      HROPiP.SetCellContent(ti,tj,tROPi-Result[2]);
      HRSPiP.SetCellContent(ti,tj,tRSPi-Result[4]);
      fBWSpaceTime(Result,MP, MP/MK*tPt, tPhi,
		   T, Rho0, RhoA, Rx, Ry, Tau, Deltat);
      HROKP.SetCellContent(ti,tj,tROK-Result[2]);
      HRSKP.SetCellContent(ti,tj,tRSK-Result[4]);

    }
    HV2Pi.SetBinContent(ti,v2Pi/yieldPi);
    HPtPi.SetBinContent(ti,yieldPi);
    HV2K.SetBinContent(ti,v2K/yieldK);
    HPtK.SetBinContent(ti,yieldK);
    HV2P.SetBinContent(ti,v2P/yieldP);
    HPtP.SetBinContent(ti,yieldP);
  }
  HRO2Pi.Write();
  HRS2Pi.Write();
  HROS2Pi.Write();
  HRL2Pi.Write();
  HROPi.Write();
  HRSPi.Write();
  HV2Pi.Write();
  HPtPi.Write();
  HRO2K.Write();
  HRS2K.Write();
  HROS2K.Write();
  HRL2K.Write();
  HROK.Write();
  HRSK.Write();
  HV2K.Write();
  HPtK.Write();
  HRO2P.Write();
  HRS2P.Write();
  HROS2P.Write();
  HRL2P.Write();
  HROP.Write();
  HRSP.Write();
  HV2P.Write();
  HPtP.Write();
  HROPiK.Write();
  HRSPiK.Write();
  HROPiP.Write();
  HRSPiP.Write();
  HROKP.Write();
  HRSKP.Write();
  outFile.Close();
}
