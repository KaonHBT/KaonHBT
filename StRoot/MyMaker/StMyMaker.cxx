//  Include header files. 
#include "TFile.h"

#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StRoot/StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"


#include "StMyMaker.h"


ClassImp(StMyMaker) //ROOT related
  
//=========================================================
StMyMaker::StMyMaker(const Char_t *name) : StMaker(name) {
  mFileName = "";
}

//=========================================================
StMyMaker::~StMyMaker(){}

//=========================================================
Int_t StMyMaker::Init() {
   cout<<"volam StMyMaker::Init()"<<endl;
  
  mNEventsFailed=0;
  mNEventsPassed=0;

  //---- Event-wise histograms ----//
  hRefMult = new TH1F("hRefMult","Ref Mult",100,0,1000);
  hVz_passed = new TH1F("hVz_passed","Primary Vz of good events",100,-30.,30.);
  hVz_failed = new TH1F("hVz_failed","Primary Vz of bad events",100,-150.,150.);
  hVx_position = new TH1F("hVx_position","Primary Vx of all events", 100, -10.,10.);//edit 
  hVy_position = new TH1F("hVy_position","Primary Vy of all events", 100, -10.,10.);//edit
  hVxVy_position = new TH2F("hVxVy_position", "Primary VxVy position of all events", 100, -5., 5., 100, -5., 5.);//edit
  hdEdx = new TH2F("hdEdx","dEdx", 100, 0., 3., 100, 0., 0.0004 );//edit
          
          
  //---- Track-wise histograms ----//
  hPhi = new TH1F("hPhi","Phi",64,-3.2,3.2);
  hPt = new TH1F("hPt","Pt",100,2.,10.);
  
  return StMaker::Init();
 
}

//=========================================================
void StMyMaker::Clear(Option_t *option){
  cout<<"volam StMyMaker::Clear() - zatim prazdne"<<endl;
}

//=========================================================
Int_t StMyMaker::Finish() {
  cout<<"volam StMyMaker::Finish()"<<endl;
  // Output file
  TFile *mFile = new TFile(mFileName, "RECREATE");
  cout << "The output filename is " << mFileName.Data() << endl;
  cout << "Events passed: " << mNEventsPassed << endl;
  cout << "Events failed: " << mNEventsFailed << endl;

  hRefMult->Write();
  hVz_passed->Write();
  hVz_failed->Write();
  hVx_position->Write();//edit
  hVy_position->Write();//edit
  hVxVy_position->Write();//edit
  hdEdx->Write();//edit
  
  hPhi->Write();
  hPt->Write();
 
  mFile->Close();
  delete mFile;
  
  return kStOK;
}


//=========================================================
Int_t StMyMaker::Make() {
  cout<<"volam StMyMaker::Make()"<<endl;

  mCurrentMu= (StMuDst*)GetInputDS("MuDst");
  
  if(!mCurrentMu || !mCurrentMu->event()) {
    cout<<" POZOR: prazdny event - asi neni nacten!!!!"<<endl;
    mNEventsFailed++;
    return kStWarn;
  }


  //event cut
  if(!acceptEvent(mCurrentMu)) {
    mNEventsFailed++;
    //------ Fill event-wise histograms for failed events here -----//
    hVz_failed->Fill(mCurrentMu->event()->primaryVertexPosition().z());
    return kStOk;
  }
  mNEventsPassed++;
  
  //------ Fill event-wise histograms for passed event here -----//
  hVz_passed->Fill(mCurrentMu->event()->primaryVertexPosition().z());
  hRefMult->Fill(mCurrentMu->event()->refMult() );
  hVx_position->Fill(mCurrentMu->event()->primaryVertexPosition().x() ); //edit
  hVy_position->Fill(mCurrentMu->event()->primaryVertexPosition().y() ); //edit
  hVxVy_position->Fill(mCurrentMu->event()->primaryVertexPosition().x(), mCurrentMu->event()->primaryVertexPosition().y() );//edit
  //------loop over tracks
  TObjArray *tracks = mCurrentMu->primaryTracks();
  TIter next(tracks);
  
  StMuTrack *track=0;
  while ( (track = (StMuTrack*)next()) ) {
  
    if(!acceptTrack(track)) continue; // track cut
    //------- Fill track-wise histograms -------//
    hPhi->Fill(track->phi());
    hPt->Fill(track->pt());
    hdEdx->Fill(track->pt(), track->dEdx() );
  }

  return kStOK;
}

//=========================================================
Bool_t StMyMaker::acceptEvent(StMuDst* mu){
  float z= fabs(mu->event()->primaryVertexPosition().z());
  if (z<25.) return kTRUE; //return 1;
  else return kFALSE; //return 0;
}

/*
//=========================================================
Bool_t StMyMaker::acceptEvent(StMuDst* mu){
  float x= fabs(mu->event()->primaryVertexPosition().x());
  if (x<25.) return kTRUE; //return 1;
  else return kFALSE; //return 0;
}
 */
 
//=========================================================
Bool_t StMyMaker::acceptTrack(StMuTrack *trk){
  return kTRUE; //beru vsechno..
}
