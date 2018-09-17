#ifndef PionDstMaker_hh
#define PionDstMaker_hh
#include "StMaker.h"
#include "TString.h"
#include "TArrayI.h"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "StFlowAnalysisMaker/StFlowAnalysisMaker.h"
#include "StRefMultCorr/StRefMultCorr.h"

//Leave this in
#ifndef St_NO_NAMESPACES
using std::string;
#endif

class Event;
class Track;

class TFile;
class StMuDst;
class StMuDstMaker;
class StMuEvent;
class StMuTrack;

class StFlowAnalysisMaker;
//class StFlowPhiWgtMaker;
//class StFlowMaker;

class PionDstMaker : public StMaker {

public:
  
  PionDstMaker(StMuDstMaker* maker);
  ~PionDstMaker() {;}
  
  void Clear(Option_t *option="");
  Int_t Init();
  Int_t Make();
  Int_t Finish();
  
  StMuDst* muDst();
  void SetMuDstMaker(StMuDstMaker*);
  void SetRefmultCorrUtil(StRefMultCorr*);
  void SetNFiles(Int_t);
  void SetStFAM(StFlowAnalysisMaker*);
  void SetFileName(TString fileName) {mFileName = fileName;}
  
private:
  
  bool acceptPionTrack(StMuTrack*);
  
  StMuDstMaker *mMuDstMaker; //!
  StFlowAnalysisMaker* mStFAM;
  TFile* mFile[16]; //!
  TString mFileName; //!
  TString mLastMuFile; //!
  StMuEvent *mMuEvent;         //! pointer to Mu-DST Event array
  TObjArray* mMuTracks;        //! Mu-DST Primary Tracks
  TObjArray* mMuGlobalTracks;  //! Mu-DST Global Tracks
  
  TTree* mTree[16]; //!;
  Int_t nFiles;                 // Number of output files, generally for centrality dependence
  int nEventsPassed;
  int nEventsFailed;
  Event* mEvent;
  int nBytes;
  int bin;

  StRefMultCorr* mRefmultCorrUtil;
  Double_t mZdcCoincidenceRate;

  Float_t CalcDcaSigned(const StThreeVectorF vertex, const StPhysicalHelixD helix);

  ClassDef(PionDstMaker,1)
};

inline void PionDstMaker::SetMuDstMaker(StMuDstMaker* f) {mMuDstMaker=f;}

inline void PionDstMaker::SetRefmultCorrUtil(StRefMultCorr* f) {mRefmultCorrUtil=f;}

inline void PionDstMaker::SetStFAM(StFlowAnalysisMaker* f) {mStFAM = f;}

inline void PionDstMaker::SetNFiles(Int_t n) {nFiles = n;}

#endif
