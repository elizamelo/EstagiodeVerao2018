
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
#include "plugins/deltaR_deltaPhi.h"

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
	auto * h_DiMuon_Mass = new TH1D("h_DiMuon_Mass", ";M_{#mu#mu} (GeV);", 80, 0, 200);
        
	//missing ET, neutrinos
	


        // photons
        dataTree->SetBranchStatus("nPho",1);
       	TTreeReaderValue< int > nPho(*dataReader, "nPho");
        dataTree->SetBranchStatus("phoEt",1);
	TTreeReaderArray< float > phoEt(*dataReader, "phoEt");
        dataTree->SetBranchStatus("phoEta",1);
	TTreeReaderArray< float > phoEta(*dataReader, "phoEta");
        dataTree->SetBranchStatus("phoPhi",1);
	TTreeReaderArray< float > phoPhi(*dataReader, "phoPhi");
        dataTree->SetBranchStatus("phoSCEta",1);
        TTreeReaderArray< float > phoSCEta(*dataReader, "phoSCEta");
        dataTree->SetBranchStatus("phoIDMVA",1);
	TTreeReaderArray< float > phoIDMVA(*dataReader, "phoIDMVA");
        dataTree->SetBranchStatus("phoR9",1);
	TTreeReaderArray< float > phoR9(*dataReader, "phoR9");
        dataTree->SetBranchStatus("phoIDbit",1);
	TTreeReaderArray< unsigned short > phoIDbit(*dataReader, "phoIDbit");
        dataTree->SetBranchStatus("phoEleVeto",1);
	TTreeReaderArray< int > phoEleVeto(*dataReader, "phoEleVeto");

	auto * h_LeadingPhoton_Pt = new TH1D("h_LeadingPhoton_Pt", ";p_{T}^{lead #gamma} (GeV);", 150, 0.0, 150.0);
	auto * h_LeadingPhoton_Eta = new TH1D("h_LeadingPhoton_Eta", ";#eta^{lead #gamma};", 30, -2.4, 2.4);
	auto * h_LeadingPhoton_Phi = new TH1D("h_LeadingPhoton_Phi", ";#phi^{lead #gamma};", 70, -3.2, 3.2);

	auto * h_DiMuonPhoton_Pt = new TH1D("h_DiMuonPhoton_Pt", ";p_{T}^{#mu#mu#gamma} (GeV);", 150, 0.0, 150.0);
	auto * h_DiMuonPhoton_Eta = new TH1D("h_DiMuonPhoton_Eta", ";#eta_{#mu#mu#gamma};", 30, -2.4, 2.4);
	auto * h_DiMuonPhoton_Phi = new TH1D("h_DiMuonPhoton_Phi", ";#phi_{#mu#mu#gamma};", 70, -3.2, 3.2);
	auto * h_DiMuonPhoton_Mass = new TH1D("h_DiMuonPhoton_Mass", ";M_{#mu#mu#gamma} (GeV);", 80, 50, 200);


	auto * h_deltaR_Leading_Trailing = new TH1D("h_deltaR_Leading_Trailing", ";#DeltaR(lead #mu, trail #mu);", 100, 0.0, 1.0);
	h_deltaR_Leading_Trailing->Sumw2();
	auto * h_deltaR_Leading_Photon = new TH1D("h_deltaR_Leading_Photon", ";#DeltaR(lead #mu, #gamma);", 100, 0.0, 5.0);
	h_deltaR_Leading_Photon->Sumw2();
	auto * h_deltaR_Trailing_Photon = new TH1D("h_deltaR_Trailing_Photon", ";#DeltaR(trail #mu, #gamma);", 100, 0.0, 5.0);
	h_deltaR_Trailing_Photon->Sumw2();
	auto * h_deltaR_Meson_Photon = new TH1D("h_deltaR_Meson_Photon", ";#DeltaR(#mu#mu, #gamma);", 100, 0.0, 5.0);
	h_deltaR_Meson_Photon->Sumw2();
	auto * h_deltaPhi_Meson_Photon = new TH1D("h_deltaPhi_Meson_Photon", ";|#Delta#phi(#mu#mu, #gamma)|;", 100, 0.0, 4.0);
	h_deltaPhi_Meson_Photon->Sumw2();


        auto * h_kin_deltaR_Leading_Trailing = new TH1D("h_kin_deltaR_Leading_Trailing", ";#DeltaR(lead #mu, trail #mu);", 100, 0.0, 1.0);
        h_kin_deltaR_Leading_Trailing->Sumw2();
        auto * h_kin_deltaR_Leading_Photon = new TH1D("h_kin_deltaR_Leading_Photon", ";#DeltaR(lead #mu, #gamma);", 100, 0.0, 5.0);
        h_kin_deltaR_Leading_Photon->Sumw2();
        auto * h_kin_deltaR_Trailing_Photon = new TH1D("h_kin_deltaR_Trailing_Photon", ";#DeltaR(trail #mu, #gamma);", 100, 0.0, 5.0);
        h_kin_deltaR_Trailing_Photon->Sumw2();
        auto * h_kin_deltaR_Meson_Photon = new TH1D("h_kin_deltaR_Meson_Photon", ";#DeltaR(#mu#mu, #gamma);", 100, 0.0, 5.0);
        h_kin_deltaR_Meson_Photon->Sumw2();
        auto * h_kin_deltaPhi_Meson_Photon = new TH1D("h_kin_deltaPhi_Meson_Photon", ";|#Delta#phi(#mu#mu, #gamma)|;", 100, 0.0, 4.0);
        h_kin_deltaPhi_Meson_Photon->Sumw2();
        auto * h_MesonPhoton_Pt = new TH1D("h_MesonPhoton_Pt", ";p_{T}^{#mu#mu#gamma} (GeV);", 150, 0.0, 150.0);
        auto * h_MesonPhoton_Eta = new TH1D("h_MesonPhoton_Eta", ";#eta_{#mu#mu#gamma};", 30, -2.4, 2.4);
        auto * h_MesonPhoton_Phi = new TH1D("h_MesonPhoton_Phi", ";#phi_{#mu#mu#gamma};", 70, -3.2, 3.2);
        auto * h_MesonPhoton_Mass = new TH1D("h_MesonPhoton_Mass", ";M_{#mu#mu#gamma} (GeV);", 80, 50, 200);

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
                        double phiMassMin = 1012.0/1000.0; //GeV
                        double phiMassMax = 1028.0/1000.0; //GeV
                        
			double miniMesonMass = 2800.0/1000.0; //GeV
                        double maxMesonMass = 3200.0/1000.0; //GeV

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
		//	}
   ////////////////////////////////////////////////////////////////////
                        // photons pre-selection  
                        if (*nPho >= 1) {
                            int nGoodPhotons = 0; //number of good photons
                            int leadingPhotonIndex = -99; //leading photons index
                            bool leadingPhotonIsSCEta = false; // 
                            bool leadingPhotonIsSCEtaG = false; //
  			    bool leadingPhotonEleVeto = false; // photon electron veto
 			    bool leadingPhotonIDMVA = false; //photon MVA ID
                            
			   for (int iPhoton = 0; iPhoton < *nPho; iPhoton++) { //loop over photon looking for the leading muon
				if (phoEt[iPhoton] > 20.0 && fabs(phoEta[iPhoton]) < 2.5) { // et > 20 and Abs(Eta) < 2.5
					leadingPhotonIndex = iPhoton;
					// more info: https://github.com/cmkuo/ggAnalysis/blob/master/ggNtuplizer/plugins/ggNtuplizer_muons.cc
					//leadingMuIsTight = (((muIDbit[iMuon] >> 2) & 1) == 1) ? true : false; // is tight muon
			                leadingPhotonIsSCEta = (fabs(phoSCEta[iPhoton]) < 2.5) ? true : false; 
  			                leadingPhotonIsSCEtaG = (!(fabs(phoSCEta[iPhoton]) > 1.4442 && fabs(phoSCEta[iPhoton]) < 1.566)) ? true : false;
			                leadingPhotonEleVeto = (phoEleVeto[iPhoton] != 0) ? true : false; 
			                leadingPhotonIDMVA = (phoIDMVA[iPhoton] > 0.2) ? true : false;
					nGoodPhotons++;
					break;
				}
			}
                          
                          // pass trigger AND two good muons AND both muons are tight muons && one photon
			  //if (goodTriggerEvt == true && nGoodMuons > 0 && trailingMuonIsTight == true && leadingMuonIsTight == true && leadingPhotonIsSCEta == true && leadingPhotonIsSCEtaG == true && leadingPhotonEleVeto == true && leadingPhotonIDMVA == true && (dimuon->M()>phiMassMin && dimuon->M()< phiMassMax ) ) {
			  if (goodTriggerEvt == true && nGoodMuons > 0 && trailingMuonIsTight == true && leadingMuonIsTight == true && leadingPhotonIsSCEta == true && leadingPhotonIsSCEtaG == true && leadingPhotonEleVeto == true && leadingPhotonIDMVA == true ) {
				//leading photon
				auto * leadingPhoton = new TLorentzVector();
				leadingPhoton->SetPtEtaPhiM(phoEt[leadingPhotonIndex], phoEta[leadingPhotonIndex], phoPhi[leadingPhotonIndex], 0);
				// dimuon+pho
				auto * dimuonpho = new TLorentzVector(); // dimuonpho = leading muon + trailing muon + leading photon
				dimuonpho->SetPtEtaPhiM(
						(*leadingMuon+*trailingMuon+*leadingPhoton).Pt(), 
						(*leadingMuon+*trailingMuon+*leadingPhoton).Eta(),
						(*leadingMuon+*trailingMuon+*leadingPhoton).Phi(),
						(*leadingMuon+*trailingMuon+*leadingPhoton).M()
						);


				// FILL HISTOS
				h_LeadingPhoton_Pt->Fill(leadingPhoton->Pt());
				h_LeadingPhoton_Eta->Fill(leadingPhoton->Eta());
				h_LeadingPhoton_Phi->Fill(leadingPhoton->Phi());	

				h_DiMuonPhoton_Pt->Fill(dimuonpho->Pt());
				h_DiMuonPhoton_Eta->Fill(dimuonpho->Eta());
				h_DiMuonPhoton_Phi->Fill(dimuonpho->Phi());
				h_DiMuonPhoton_Mass->Fill(dimuonpho->M());
				
				h_deltaR_Leading_Trailing->Fill(deltaR(leadingMuon->Eta(), leadingMuon->Phi(), trailingMuon->Eta(), trailingMuon->Phi()));
                                h_deltaR_Leading_Photon->Fill(deltaR(leadingMuon->Eta(), leadingMuon->Phi(), leadingPhoton->Eta(), leadingPhoton->Phi()));
                                h_deltaR_Trailing_Photon->Fill(deltaR(trailingMuon->Eta(), trailingMuon->Phi(), leadingPhoton->Eta(), leadingPhoton->Phi()));
                                h_deltaR_Meson_Photon->Fill(deltaR(dimuon->Eta(), dimuon->Phi(), leadingPhoton->Eta(), leadingPhoton->Phi()));
                                h_deltaPhi_Meson_Photon->Fill(fabs(deltaPhi(dimuon->Phi(), leadingPhoton->Phi())));
                              
                                auto goodDeltaRMuonsPhoton = true;
				auto goodDeltaRPhiMesonPhoton = true; 
				auto goodMesonMassCut = true;
                                goodDeltaRMuonsPhoton = ( (deltaR(leadingMuon->Eta(), leadingMuon->Phi(), leadingPhoton->Eta(), leadingPhoton->Phi()) > 1.0) && (deltaR(trailingMuon->Eta(), trailingMuon->Phi(), leadingPhoton->Eta(), leadingPhoton->Phi()) > 1.0)) ? true : false;
                                goodDeltaRPhiMesonPhoton = ( (deltaR(dimuon->Eta(), dimuon->Phi(), leadingPhoton->Eta(), leadingPhoton->Phi()) > 2.0) && (fabs(deltaPhi(dimuon->Phi(), leadingPhoton->Phi())) > 1.5) ) ? true : false;
                                goodMesonMassCut = ((dimuon->M() > miniMesonMass && dimuon->M() < maxMesonMass)) ? true : false;
				if(goodDeltaRMuonsPhoton==true && goodDeltaRPhiMesonPhoton==true && goodMesonMassCut==true){
				       h_MesonPhoton_Pt->Fill(dimuonpho->Pt());
	                               h_MesonPhoton_Eta->Fill(dimuonpho->Eta());
	                               h_MesonPhoton_Phi->Fill(dimuonpho->Phi());
         	                       h_MesonPhoton_Mass->Fill(dimuonpho->M());

                	               h_kin_deltaR_Leading_Trailing->Fill(deltaR(leadingMuon->Eta(), leadingMuon->Phi(), trailingMuon->Eta(), trailingMuon->Phi()));
                        	       h_kin_deltaR_Leading_Photon->Fill(deltaR(leadingMuon->Eta(), leadingMuon->Phi(), leadingPhoton->Eta(), leadingPhoton->Phi()));
                                       h_kin_deltaR_Trailing_Photon->Fill(deltaR(trailingMuon->Eta(), trailingMuon->Phi(), leadingPhoton->Eta(), leadingPhoton->Phi()));
                                       h_kin_deltaR_Meson_Photon->Fill(deltaR(dimuon->Eta(), dimuon->Phi(), leadingPhoton->Eta(), leadingPhoton->Phi()));
                                       h_kin_deltaPhi_Meson_Photon->Fill(fabs(deltaPhi(dimuon->Phi(), leadingPhoton->Phi())));


				
                                }//deltasRPhi
                           }//good pho
                              
                                  
			
                        }//pho sel
                    }//dimuon sel

		}//ja estava

    } // end loop over events

    // post-processing 
    cout << "\n\n\nWriting output file... (" << outFileAppend <<  ")"  << endl;
    outFile->Write();
    cout << "\nDone  (" << outFileAppend <<  ")" << "!\n\n\n\n\n" << endl;

} //end ana_ggNtuplesAnalyzer


