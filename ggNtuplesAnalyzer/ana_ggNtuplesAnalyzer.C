
#include <iostream>
#include <math.h>  

#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TROOT.h"
#include "TH1.h"

#include "plugins/ggNtuplesFilesReader.h"

using namespace std;

const auto muonMass = 105.6583745/1000.0; //GeV


void ana_ggNtuplesAnalyzer(vector<string> ggNtuplesFiles, int nFiles = -1, string outFileAppend = "")  
{
	// output file
	string outputFileName = "outputFiles/histos_ggNtuplesAnalyzer_"+ outFileAppend + ".root";
	TFile * outFile  = new TFile(outputFileName.c_str(),"RECREATE","ggNtuplesAnalyzer histograms");

	// loads the ggNtuples 
	if (nFiles > 0) {
		ggNtuplesFiles.resize(nFiles);	
	}
	// more info: https://root.cern.ch/7-using-ttreereader
	TTreeReader * dataReader = ggNtuplesFilesReader( ggNtuplesFiles, "ggNtuplizer/EventTree" );
	TTree * dataTree = dataReader->GetTree();

	

	////////////////////////////////////////////////////////////////////
	// activate branchs and define readers
	dataTree->SetBranchStatus("*",0);

	dataTree->SetBranchStatus("HLTEleMuX",1);
	TTreeReaderValue< ULong64_t > HLTEleMuX(*dataReader, "HLTEleMuX");

	// muons
	dataTree->SetBranchStatus("nMu",1);
	TTreeReaderValue< int > nMu(*dataReader, "nMu");
	dataTree->SetBranchStatus("muPt",1);
	TTreeReaderArray< float > muPt(*dataReader, "muPt");
	dataTree->SetBranchStatus("muEta",1);
	TTreeReaderArray< float > muEta(*dataReader, "muEta");
	dataTree->SetBranchStatus("muPhi",1);
	TTreeReaderArray< float > muPhi(*dataReader, "muPhi");
	dataTree->SetBranchStatus("muCharge",1);
	TTreeReaderArray< int > muCharge(*dataReader, "muCharge");
	dataTree->SetBranchStatus("muIDbit",1);
	TTreeReaderArray< unsigned short > muIDbit(*dataReader, "muIDbit");
 
	auto * h_LeadingMuon_Pt = new TH1D("h_LeadingMuon_Pt", ";p_{T}^{lead #mu} (GeV);", 150, 0.0, 150.0);
	auto * h_LeadingMuon_Eta = new TH1D("h_LeadingMuon_Eta", ";#eta^{lead #mu};", 30, -2.4, 2.4);
	auto * h_LeadingMuon_Phi = new TH1D("h_LeadingMuon_Phi", ";#phi^{lead #mu};", 70, -3.2, 3.2);

	auto * h_TrailingMuon_Pt = new TH1D("h_TrailingMuon_Pt", ";p_{T}^{trail #mu} (GeV);", 150, 0.0, 150.0);
	auto * h_TrailingMuon_Eta = new TH1D("h_TrailingMuon_Eta", ";#eta_{trail #mu};", 30, -2.4, 2.4);
	auto * h_TrailingMuon_Phi = new TH1D("h_TrailingMuon_Phi", ";#phi_{trail #mu};", 70, -3.2, 3.2);

	auto * h_DiMuon_Pt = new TH1D("h_DiMuon_Pt", ";p_{T}^{#mu#mu} (GeV);", 150, 0.0, 150.0);
	auto * h_DiMuon_Eta = new TH1D("h_DiMuon_Eta", ";#eta_{#mu#mu};", 30, -2.4, 2.4);
	auto * h_DiMuon_Phi = new TH1D("h_DiMuon_Phi", ";#phi_{#mu#mu};", 70, -3.2, 3.2);
	auto * h_DiMuon_Mass = new TH1D("h_DiMuon_Mass", ";M_{#mu#mu} (GeV);", 80, 50, 130);

	////////////////////////////////////////////////////////////////////
	// numer of entries
	auto totalEvts = dataTree->GetEntries();
	auto printEvery = 10000;
	cout << "\nN. Entries (" << outFileAppend <<  "): " << totalEvts << endl;
	cout << "\nPrinting every: " << printEvery << " events" << endl;
	cout << "\nLooping over events... \n" << endl;


	////////////////////////////////////////////////////////////////////
	// main loop
	auto iEvt = 0;
	while (dataReader->Next()) { // loop over events
		if (iEvt % printEvery == 0) cout << "----------------------------------------> Events read (" << outFileAppend <<  "): " << iEvt << " / " << totalEvts << " - ~"<< round(((float)iEvt/(float)totalEvts)*100) << "%"<< endl;
		iEvt++;


		////////////////////////////////////////////////////////////////////
		// trigger selection 
		// more info: https://github.com/cmkuo/ggAnalysis/blob/master/ggNtuplizer/plugins/ggNtuplizer_globalEvent.cc#L151
		auto goodTriggerEvt = true;
		// goodTriggerEvt = (((*HLTEleMuX >> 8) & 1) == 1) ? true : false; // HLT_Mu17_Photon30_CaloIdL_L1ISO_v*
		goodTriggerEvt = (((*HLTEleMuX >> 19) & 1) == 1) ? true : false; // HLT_IsoMu24_v*


		////////////////////////////////////////////////////////////////////
		// muons pre-selection  
		if (*nMu >= 2 && goodTriggerEvt) {
			int nGoodMuons = 0; //number of good muons

			int leadingMuonIndex = -99; //leading muon index
			bool leadingMuonIsTight = false; //leading muon id (tight)

			int trailingMuonIndex = -99; //leading muon index	
			bool trailingMuonIsTight = false; //leading muon id (tight)

			for (int iMuon = 0; iMuon < *nMu; iMuon++) { //loop over muons looking for the leading muon
				if (muPt[iMuon] > 27.0 && fabs(muEta[iMuon]) < 2.4) { // pt > 27 and Abs(Eta) < 2.4
					leadingMuonIndex = iMuon;
					// more info: https://github.com/cmkuo/ggAnalysis/blob/master/ggNtuplizer/plugins/ggNtuplizer_muons.cc#L170
					leadingMuonIsTight = (((muIDbit[iMuon] >> 2) & 1) == 1) ? true : false; // is tight muon
					nGoodMuons++;
					break;
				}
			}

			if (nGoodMuons > 0 ) {
				if (leadingMuonIndex + 1 < *nMu) {
					for (int iMuon = leadingMuonIndex + 1; iMuon < *nMu; iMuon++) { //loop over muons looking for the trailing muon
						if (muPt[iMuon] > 2.0 && fabs(muEta[iMuon]) < 2.4) { // pt > 27 and Abs(Eta) < 2.4
							if (muCharge[iMuon] * muCharge[leadingMuonIndex] < 0) { // the muons should have opposite charges
								trailingMuonIndex = iMuon;
								// more info: https://github.com/cmkuo/ggAnalysis/blob/master/ggNtuplizer/plugins/ggNtuplizer_muons.cc#L170
								trailingMuonIsTight = (((muIDbit[iMuon] >> 2) & 1) == 1) ? true : false; // is tight muon
								nGoodMuons++;
								break;
							}
						}
					}
				}
			}


			// pass trigger AND two good muons AND both muons are tight muons
			if (goodTriggerEvt == true && nGoodMuons > 1 && trailingMuonIsTight == true && leadingMuonIsTight == true) {
				//leading muon
				auto * leadingMuon = new TLorentzVector();
				leadingMuon->SetPtEtaPhiM(muPt[leadingMuonIndex], muEta[leadingMuonIndex], muPhi[leadingMuonIndex], muonMass);

				//trailing muon
				auto * trailingMuon = new TLorentzVector();
				trailingMuon->SetPtEtaPhiM(muPt[trailingMuonIndex], muEta[trailingMuonIndex], muPhi[trailingMuonIndex], muonMass);	

				// dimuon
				auto * dimuon = new TLorentzVector(); // dimuon = leading muon + trailing muon
				dimuon->SetPtEtaPhiM(
						(*leadingMuon+*trailingMuon).Pt(), 
						(*leadingMuon+*trailingMuon).Eta(),
						(*leadingMuon+*trailingMuon).Phi(),
						(*leadingMuon+*trailingMuon).M()
						);	

				// FILL HISTOS
				h_LeadingMuon_Pt->Fill(leadingMuon->Pt());
				h_LeadingMuon_Eta->Fill(leadingMuon->Eta());
				h_LeadingMuon_Phi->Fill(leadingMuon->Phi());

				h_TrailingMuon_Pt->Fill(trailingMuon->Pt());
				h_TrailingMuon_Eta->Fill(trailingMuon->Eta());
				h_TrailingMuon_Phi->Fill(trailingMuon->Phi());

				h_DiMuon_Pt->Fill(dimuon->Pt());
				h_DiMuon_Eta->Fill(dimuon->Eta());
				h_DiMuon_Phi->Fill(dimuon->Phi());
				h_DiMuon_Mass->Fill(dimuon->M());
			}

		}

    } // end loop over events

    // post-processing 
    cout << "\n\n\nWriting output file... (" << outFileAppend <<  ")"  << endl;
    outFile->Write();
    cout << "\nDone  (" << outFileAppend <<  ")" << "!\n\n\n\n\n" << endl;

} //end ana_ggNtuplesAnalyzer


