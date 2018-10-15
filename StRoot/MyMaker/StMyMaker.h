#ifndef StMyMaker_hh     
#define StMyMaker_hh
//
//  Include files
#include "StMaker.h"
#include "StMuDSTMaker/COMMON/StMuTypes.hh"
#include <TH1.h>
#include <TH2.h>
//#include <TH3.h> //edit 
#include <TString.h> 

//  Forward declarations
class StMuTrack;
class TFile;
class TH1F;
class TH2F;
//class TH3F;//edit

class StMuDst;
class StMuTrack;

#ifndef ST_NO_NAMESPACES
using std::string;
#endif 

class StMyMaker : public StMaker {
 public:
  
  StMyMaker(const Char_t *name="myAnalysis");   // constructor
  virtual ~StMyMaker();                         // destructor - idelane visrtualni!
  
  //zakladni zdedene procedury z StMaker
  virtual Int_t    Init();                   // called once at the beginning of your job
  virtual void     Clear(Option_t *option=""); // called after every event to cleanup
  virtual Int_t    Make();                   // invoked for every event
  virtual Int_t    Finish();                 // called once at the end
    
  void SetOutputFileName(TString name);
  
 protected:
  
  virtual Bool_t acceptEvent(StMuDst* mu);
  virtual Bool_t acceptTrack(StMuTrack *trk);
  
 private:
    StMuDst*     mCurrentMu;
    TString      mFileName;       //! Name of output root file
   
    Int_t        mNEventsPassed;  //!
    Int_t        mNEventsFailed;  //!
   
    //histogramy
    TH1F* hRefMult;
    TH1F* hVz_passed;
    TH1F* hVz_failed;
    TH1F* hPhi;
    TH1F* hPt;
    TH1F* hVx_position;
    TH1F* hVy_position;
    TH2F* hVxVy_position;
    TH2F* hdEdx;
    //  TH3F* hVxVyVz;
   
    ClassDef(StMyMaker,1)     //ROOT related
};

inline void StMyMaker::SetOutputFileName(TString name){mFileName =name;}

#endif                
