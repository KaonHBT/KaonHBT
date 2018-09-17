TEventList* makeEventList(TTree* tree, const Float_t zdcCut);
void makeBES19(const TString fileList = "test200.list"
				, const TString outFile = "testOut.root"
                , const Int_t nOutFiles = 1
				, Int_t nEvents = 999 
				, const Int_t nPass = 1 
                , const Int_t centLo = 0
                , const Int_t centHigh = 15
				)
{

    Bool_t  dailyZdc = kFALSE;        

    // UU/AuAu200
    Int_t  nTofMatchCut = -1;             // Number of Tof matched tracks 
    Float_t VzCut = 70;                 // Cut on z-vertex
    Float_t VrCut = 2.;                // Cut on r-vertex
    Float_t VrCenter[2] = {0.,0.};        // Beam center - (0, -0.89) for 14.6 GeV
    Float_t etaGap = 10;
    Float_t triggers[3] = {340001, 340011, 340021}; 
    Int_t nTriggers = 3;
	// Float_t triggers[6] = {400104, 400108, 400114, 400118, 400124, 400134}; // U+U 
    // Int_t nTriggers = 6; // U+U
    Int_t  zdcCut = 0;             // Zdc threshold cut - Au+Au 200 GeV
	

    Float_t ptRange[2] = {0.15,2.};      // pT Cut
    Int_t cent9Cut[2] = {0,8};
    Int_t cent16Cut[2] = {centLo,centHigh};

	TStopwatch* timer = new TStopwatch();

	//----------------- Load Libraries ---------------//
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	gSystem->Load("StFlowMaker");
	gSystem->Load("StFlowAnalysisMaker");
	gSystem->Load("PionDstMaker");
	gSystem->Load("StRefMultCorr");
	gSystem->Load("EventCut");

	//----------------- Instantiate chain and MuDst reader ---------------//
	StChain* chain = new StChain("StChain");
	chain->SetDebug(0);
	StMuDebug::setLevel(0);

	StMuDstMaker* muMaker = new StMuDstMaker(0,0,"",fileList.Data(),"st:MuDst.root", 100, "MuDst");
    StRefMultCorr* refmultCorrUtil = new StRefMultCorr("refmult");

	//----------------------- Define analysis makers -----------------------//

	StFlowMaker* flowMaker = new StFlowMaker();
	flowMaker->MuEventRead(kTRUE);
	flowMaker->SetMuDstMaker(muMaker);
	flowMaker->SetPhiWgtCalc(kFALSE);
	flowMaker->SetReCentCalc(kFALSE);

	StFlowAnalysisMaker* StFAM = new StFlowAnalysisMaker();
	StFAM->SetHistoRanges(kFALSE);
	StFAM->SetPtRange_for_vEta(ptRange[0],ptRange[1]);      // For filling v_n(eta) histograms
	StFAM->SetEtaRange_for_vPt(0.,1.3);                     // For filling v_n(pT) histograms
	StFAM->SetV1Ep1Ep2(kFALSE);

	Bool_t bFillPsi = kFALSE; //Should we calculate the psi-shift correction coefficients?
	Bool_t bDoPsi = kFALSE; //Should we *perform* the psi-shift correction?

	if(nPass == 1){ bFillPsi = kTRUE; }
	if(nPass == 2){ bDoPsi = kFALSE; }

	StFAM->SetFillPsiShift(bFillPsi);
	StFAM->SetDoPsiShift(bDoPsi);

	if(nPass == 2){
		PionDstMaker* dstMaker = new PionDstMaker(muMaker);
        dstMaker->SetRefmultCorrUtil(refmultCorrUtil);
		dstMaker->SetStFAM(StFAM);
		dstMaker->SetFileName(outFile.Data());
        dstMaker->SetNFiles(nOutFiles);
	}
	//----------------------- Init Chain and set cuts -----------------------//
	chain->Init();
//	StFlowCutEvent::SetVertexZ(-1*VzCut, VzCut);
	StFlowCutEvent::SetVertexZ(0,0);
	StFlowCutEvent::SetVertexR(0.,0);
	StFlowCutEvent::SetVertexRCenter(VrCenter[0],VrCenter[1]);
	
	StFlowCutTrack::IncludeTpcTracks(kFALSE);
	StFlowCutTrack::IncludeFtpcTracks(kFALSE);
	StFlowCutTrack::SetPtTpc(ptRange[0],ptRange[1]);
	StFlowCutTrack::SetEtaTpc(10.0,10.0);
	StFlowCutTrack::SetFitOverMaxPts(2.,2.05);
	//StFlowCutTrack::SetDcaGlobalTpc(0.,2.);
	
	StFlowEvent::SetEtaTpcCut(-1., -0.1, 0, 0);  // harmonic 1, selection 0, East TPC 
	StFlowEvent::SetEtaTpcCut(0.1, 1., 0, 1);  // harmonic 1, selection 1, West TPC 
	StFlowEvent::SetEtaTpcCut(-0.75, 0.75, 1, 0);  // harmonic 2, selection 0, mult analysis
	StFlowEvent::SetEtaTpcCut(-1., 0., 1, 1);  // harmonic 2, selection 1, q2 analysis
	StFlowEvent::SetPtWgt(kFALSE);
	StFlowEvent::SetEtaq2Cut(0.25,1);  // harmonic 2, selection 1, q2 analysis
	StFlowEvent::SetEtaGap(etaGap);  // harmonic 2, selection 1, q2 analysis
	StFlowEvent::SetEtaWgt(kFALSE);
	StFlowEvent::SetEtaSubs();

	//--------------------- Get the event list ---------------------//
    EventCut* eventCut = new EventCut(muMaker->chain());
    eventCut->SetVzCut(-1.0*VzCut, VzCut);
    eventCut->SetVrCut(VrCut);
    eventCut->SetTriggers(triggers, nTriggers);

    muMaker->SetEventList( eventCut->makeEventList() );



	//--------------------- The STAR chain Event loop ---------------------//
	// Start  bookeeping stuff. Timer, event counter, etc.
	Double_t overheadTime = timer->RealTime();
	cout << "***** Starting Event Loop. " << overheadTime  << " seconds elapsed so far .*****" << endl;
	timer->Start();
	Int_t iReturn = 0;
	Float_t percentCounter = 0.01;
	Int_t nEventsProcessed = 0;

	// Actual Event Loop
	for (Int_t iev=0;iev<nEvents; iev++) {

		Float_t progress = (Float_t)iev / (Float_t)nEvents;	
		nEventsProcessed++;

		if(progress >= percentCounter){
			cout << percentCounter * 100 << "% done" << endl;
			percentCounter += .01;
		}

		chain->Clear();
		iReturn = chain->Make(iev); 

		if (iReturn) {
			cout << "Ran out of events!" << endl;
			break;
		}
	} 

	Double_t eventLoopTime = timer->RealTime();
	cout << endl << "***** Finished Event Loop. " << eventLoopTime << " seconds to process " << nEventsProcessed << " events. "; 
	cout << (Double_t)nEventsProcessed / eventLoopTime << " events/s." << endl;

	chain->Finish();

	delete muMaker;
	delete flowMaker;
	delete StFAM;
	if(nPass == 2){ delete dstMaker;}
	delete chain;
}

TEventList* makeEventList(TTree* tree, const Float_t zdcCut) 
{

	//----------------- Read data using trees ---------------//

	TStopwatch* timer = new TStopwatch();
    TEventList* eventList = new TEventList();
    Long64_t nEntries = tree->GetEntries(); 
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("MuEvent.mZdcTriggerDetector.mAdc*",1);
	Double_t overheadTime = timer->RealTime();

	timer->Start();
    for(Int_t i = 0; i <= (nEntries-1); i++)
    {
        tree->GetEntry(i);
        TLeaf* muLeaf = tree->GetLeaf("MuEvent.mZdcTriggerDetector.mAdc");
        Float_t zdcWest = muLeaf->GetValue(0);
        Float_t zdcEast = muLeaf->GetValue(4);

        if( (zdcWest <= zdcCut) && (zdcEast <= zdcCut)) {eventList->Enter(i);}

    }
	Double_t eventLoopTime = timer->RealTime();

    cout << "\n\n***** makeEventList() Finished *****\n";
    cout << "Overhead time: " << overheadTime << "  Event loop time: " << eventLoopTime << endl;
    cout << "Processed " << nEntries << " entries. Added " << eventList->GetN() << " entries to the event list.\n\n";

    tree->SetBranchStatus("*",1);


    delete timer;

    return eventList;
}

