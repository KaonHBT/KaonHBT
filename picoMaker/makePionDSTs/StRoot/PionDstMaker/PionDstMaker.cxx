#include "PionDstMaker.h"
#include "StEventTypes.h"

#include "TFile.h"
#include "TChain.h"
#include "StMessMgr.h"

#include "StMuDSTMaker/COMMON/StMuException.hh"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuDebug.h"
#include "StMuDSTMaker/COMMON/StMuCut.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"

#include "TClonesArray.h"
#include "TBranch.h"
#include "TTree.h"

#include "Event.h"
#include "Track.h"
#include "StFlowAnalysisMaker/StFlowAnalysisMaker.h"
#include "StRefMultCorr/StRefMultCorr.h"

ClassImp(PionDstMaker)

PionDstMaker::PionDstMaker(StMuDstMaker* maker) : StMaker("PionDstMaker") {
	mLastMuFile = "";
	nEventsPassed = nEventsFailed = 0;
    nFiles = 0;

	mMuDstMaker = maker;

}

Int_t PionDstMaker::Init() {
	
	mEvent = new Event();

    for(Int_t i = 0; i <= (nFiles - 1); i++)
    {
        TString* fName = new TString(mFileName.Data());
        if(nFiles > 1)
        {
            fName->Append("_Cent");
            *fName += i;
        }
        fName->Append(".root");
        mFile[i] = new TFile(fName->Data(), "RECREATE");

        delete fName;
    
        mTree[i] = new TTree("PionDst","PionDst");
        mTree[i]->AutoSave();
        mFile[i]->SaveSelf();
        mTree[i]->Branch("Event","Event",&mEvent);

    }

	nBytes = 0;

    if(!mRefmultCorrUtil){ mRefmultCorrUtil = new StRefMultCorr("refmult"); }
	
	cout << "end of PDMB::Init()\n";
	return StMaker::Init();
}

Int_t PionDstMaker::Make() {
	
	StMuDst* muDst = mMuDstMaker->muDst();
	mMuEvent = muDst->event();

	TString currentMuFile = mMuDstMaker->chain()->GetFile()->GetName();
	if(currentMuFile != mLastMuFile) {
        mRefmultCorrUtil->init(mMuEvent->runId());
		cout << "New File: " << currentMuFile.Data() << endl;
		mLastMuFile = currentMuFile;
	}
	
	
	//Get pointers for pMuEvent, pMuTracks
	mMuTracks = 0;
	mMuGlobalTracks = 0;
	mMuTracks = muDst->primaryTracks();
	mMuGlobalTracks = muDst->globalTracks();
    mRefmultCorrUtil->initEvent(mMuEvent->refMult(), mMuEvent->primaryVertexPosition().z(), mMuEvent->runInfo().zdcCoincidenceRate());
 
	//Set the Event quantities (mimic StFlowMaker:FillFromMuDst())
	mEvent->SetRunID(mMuEvent->runId());
	mEvent->SetEventID(mMuEvent->eventId());
	mEvent->SetRefMult(mMuEvent->refMult());
	mEvent->SetRefMultCorr(mRefmultCorrUtil->getRefMultCorr());
	mEvent->SetRefMultCorrWeight(mRefmultCorrUtil->getWeight());
	mEvent->SetCent9(mRefmultCorrUtil->getCentralityBin9());
	mEvent->SetCent16(mRefmultCorrUtil->getCentralityBin16());
	mEvent->SetRefMultPos(mMuEvent->refMultPos());
	mEvent->SetRefMultNeg(mMuEvent->refMultNeg());
	mEvent->SetTofMult(mMuEvent->btofTrayMultiplicity());
	mEvent->SetVertexPos(mMuEvent->primaryVertexPosition());
	mEvent->SetVzVpd(mMuEvent->vpdVz());
	mEvent->SetZDCe(mMuEvent->zdcTriggerDetector().adc(4));
	mEvent->SetZDCw(mMuEvent->zdcTriggerDetector().adc(0));
    mEvent->SetCtbMultiplicity(mMuEvent->ctbMultiplicity());
    mEvent->SetNumberOfTracks(mMuEvent->eventSummary().numberOfTracks());
    mEvent->SetNumberOfGoodTracks(mMuEvent->eventSummary().numberOfGoodTracks());
    mEvent->SetMagField(mMuEvent->magneticField());
	
 
	for(Int_t i = 0; i <= 1; ++i)
	{
		for(Int_t j = 0; j <= 2; ++j) 
        { mEvent->SetPsi(mStFAM->Psi2(i,j),i,j); 
        }
		mEvent->Setq2(mStFAM->q2(i),i);
	}
	 
	//Set the Track quantities for the current Event
	Track *track = 0;
	StMuTrack *muTrack = 0;
	
	int nMuTracks = mMuTracks->GetEntries();
	int goodTracks = 0;
	
	for(int j=0; j<nMuTracks; j++)
	{
		//Get quantities from the track
		muTrack = (StMuTrack*)mMuTracks->UncheckedAt(j);
		if(acceptPionTrack(muTrack)) {
			
			track = mEvent->AddPionTrack();
			track->SetPt(muTrack->pt());
			
			if(muTrack->index2Global()<0) {
				gMessMgr->Info() << "PionDstMaker: index2Global<0" << endl;
				continue;
			}

			StMuTrack* muGlobalTrack = (StMuTrack*)mMuGlobalTracks->At(muTrack->index2Global());
			if (!muGlobalTrack) {
				gMessMgr->Info() << "PionDstMaker: no global track" << endl;
				continue;
			}
			
			track->SetMsquared( muTrack->p().mag() * muTrack->p().mag() *
								(1./(muTrack->btofPidTraits().beta() * muTrack->btofPidTraits().beta()) - 1.) );
			track->SetBeta(muTrack->btofPidTraits().beta());
			track->SetPtGlobal(muGlobalTrack->pt());
			track->SetTrackType(muTrack->type());
			track->SetPhi(muTrack->phi());
			track->SetPhiGlobal(muGlobalTrack->phi());
			track->SetEta(muTrack->eta());
			track->SetEtaGlobal(muGlobalTrack->eta());
			track->SetDedx(muTrack->dEdx());
			track->SetCharge(muTrack->charge());
			track->SetDcaSigned(CalcDcaSigned(mMuEvent->primaryVertexPosition(),muTrack->helix()));
			track->SetDca(muTrack->dca().mag());
			track->SetDca3(muTrack->dca());
			track->SetDcaGlobal(muTrack->dcaGlobal().mag());
			track->SetChi2(muTrack->chi2xy()); 
			track->SetTopologyMap(muTrack->topologyMap());
			track->SetNhits(muTrack->nHits());
			track->SetFitPts(muTrack->nHitsFit() - muTrack->nHitsFit(kSvtId) -
					 muTrack->nHitsFit(kSsdId) - 1); //remove additional points
			track->SetMaxPts(muTrack->nHitsPoss() - muTrack->nHitsPoss(kSvtId) -
					 muTrack->nHitsPoss(kSsdId) - 1); //remove additional points
			track->SetNdedxPts(muTrack->nHitsDedx());
			track->SetDcaGlobal3(muTrack->dcaGlobal());
			track->SetTrackLength(muTrack->helix().pathLength(mMuEvent->primaryVertexPosition())); //???
			track->SetZFirstPoint(muTrack->firstPoint().z());
			track->SetZLastPoint(muTrack->lastPoint().z());
			track->SetFlag(muTrack->flag());
			track->SetZFirstPointX(muTrack->firstPoint().x());
			track->SetZFirstPointY(muTrack->firstPoint().y());
			track->SetIndex2Global(muTrack->index2Global());
			track->SetHelix(muTrack->helix());
			
			//Set the track momentum 3-vector (needed for Femtoscopy)
			track->SetP(muTrack->p().x(),muTrack->p().y(),muTrack->p().z());
			
			//Now the PID stuff
			if (muTrack->charge() < 0) {
				track->SetPidPiMinus(muTrack->nSigmaPion()); 
				track->SetPidAntiProton(muTrack->nSigmaProton());
				track->SetPidKaonMinus(muTrack->nSigmaKaon());
				track->SetPidAntiDeuteron( 999.0 );
				track->SetPidElectron(muTrack->nSigmaElectron());
			} else {
				track->SetPidPiPlus(muTrack->nSigmaPion()); 
				track->SetPidProton(muTrack->nSigmaProton()); 
				track->SetPidKaonPlus(muTrack->nSigmaKaon()); 
				track->SetPidDeuteron( 999.0 );
				track->SetPidPositron(muTrack->nSigmaElectron());
			}
			
			//Now dEdx stuff
			if ( muTrack->nSigmaKaon() > 2.0 ) {
				if (muTrack->charge() > 0 ) {
					track->SetMostLikelihoodPID(14); // proton
					track->SetMostLikelihoodProb( 0.99 ); // guaranteed
				}	else {
					track->SetMostLikelihoodPID(15); // anti-proton
					track->SetMostLikelihoodProb( 0.99 ); // guaranteed
				}
			}
			
			if ( muTrack->nSigmaPion() > 2.0 ) {
				if (muTrack->charge() > 0 ) {
					track->SetMostLikelihoodPID(11); // kaon
					track->SetMostLikelihoodProb( 0.99 ); // guaranteed
				}	else {
					track->SetMostLikelihoodPID(12); // anti-kaon
					track->SetMostLikelihoodProb( 0.99 ); // guaranteed
				}
			}
			
			if ( muTrack->nSigmaPion() < -2.0 ) {
				if (muTrack->charge() < 0 ) {
					track->SetMostLikelihoodPID(3); // electron
					track->SetMostLikelihoodProb( 0.99 ); // guaranteed
				}	else {
					track->SetMostLikelihoodPID(2); // positron
					track->SetMostLikelihoodProb( 0.99 ); // guaranteed
				}
			}
			
			track->SetExtrapTag(0); // none are in the PID merging area
			track->SetElectronPositronProb(muTrack->pidProbElectron());
			track->SetPionPlusMinusProb(muTrack->pidProbPion());
			track->SetKaonPlusMinusProb(muTrack->pidProbKaon());
			track->SetProtonPbarProb(muTrack->pidProbProton());
			
			goodTracks++;
		} //end if(acceptTrack)
	} // end loop over tracks

    if(nFiles == 1){
        nBytes += mTree[0]->Fill();
    } else if (nFiles == 9) {
        nBytes += mTree[mRefmultCorrUtil->getCentralityBin9()]->Fill();
    } else if (nFiles == 16) {
        nBytes += mTree[mRefmultCorrUtil->getCentralityBin16()]->Fill();
    }

	mEvent->Clear();
	mMuTracks->Clear();
	mMuGlobalTracks->Clear();
		
	return kStOK;
}

Int_t PionDstMaker::Finish() {

    for(Int_t i = 0; i <= (nFiles - 1); i++)
    {
        mFile[i]->cd();
        mTree[i]->Write();
        mFile[i]->Close();
    }
	
	cout << "PionDstMaker::Finish()\n\n";
	cout << "\t  nEventsPassed: " << nEventsPassed << " events.\n";
	cout << "\t  nEventsFailed: " << nEventsFailed << " events.\n";
	cout << "Finish() ended : PionDstMaker.cxx\n";
	
	return kStOK;
}
	
//__________________________________________________________________//
void PionDstMaker::Clear(Option_t *opt) {
	StMaker::Clear();
}

bool PionDstMaker::acceptPionTrack(StMuTrack* track) {
if(track->btofPidTraits().beta()==-999){
    bool trackPass = track && 
	( track->flag() > 0 ) &&
	( track->flag() < 1000) &&
	( track->pt() > 0.1 ) &&
	( track->pt() < 0.55 ) && 
	( track->topologyMap().numberOfHits(kTpcId) >= 10 ) &&
	( fabs(track->eta()) <= 1.5 ) && 
	( track->dca().mag() < 5.0 ) &&
	( fabs(track->nSigmaKaon()) <= 3.0) &&
	( fabs(track->nSigmaPion()) > 2.0) &&
	( fabs(track->nSigmaProton()) > 2.0) &&
	( fabs(track->nSigmaElectron()) > 2.0);
	return trackPass;
	}
	else{
	bool trackPass = track && 
	( track->flag() > 0 ) &&
	( track->flag() < 1000) &&
	( track->pt() > 0.1 ) &&
	( track->pt() < 2. ) && 
	( track->topologyMap().numberOfHits(kTpcId) >= 10 ) &&
	( fabs(track->eta()) <= 1.5 ) && 
	( track->dca().mag() < 5.0 ) &&
	( fabs(track->nSigmaKaon())<=3.0) &&
	( track->p().mag() * track->p().mag() *
	  (1./(track->btofPidTraits().beta() * track->btofPidTraits().beta()) - 1.)>0.16 )&&
	( track->p().mag() * track->p().mag() *
	  (1./(track->btofPidTraits().beta() * track->btofPidTraits().beta()) - 1.)<0.35 );

  return trackPass;
  }
}

Float_t PionDstMaker::CalcDcaSigned(const StThreeVectorF vertex, 
																	 const StPhysicalHelixD helix) {
	// find the distance between the center of the circle and vertex.
	// if the radius of curvature > distance, then call it positive
	// Bum Choi

	double xCenter = helix.xcenter();
	double yCenter = helix.ycenter();
	double radius = 1.0/helix.curvature();

	double dPosCenter = ::sqrt( (vertex.x() - xCenter) * (vertex.x() - xCenter) +
														(vertex.y() - yCenter) * (vertex.y() - yCenter));

	return (Float_t)(radius - dPosCenter);
}
