#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TMath.h"
#include "TArrow.h"
#include "TGraph.h"
#include "RooPolynomial.h"
#include "RooBinning.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooErrorVar.h"
#include "RooVoigtian.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooNLLVar.h"
#include "TLatex.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "TPaveLabel.h"
#include "RooNLLVar.h"
#include "RooProfileLL.h"
#include "RooProdPdf.h"
#include "RooSimultaneous.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooProduct.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "TCut.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TRandom3.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "TStyle.h"
#include "RooCBShape.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooHistPdf.h"
#include "RooHistPdf.h"
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/SPlot.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"
#include "RooChebychev.h"
#include "RooDLLSignificanceMCSModule.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/SPlot.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"
#include "RooStats/NumberCountingPdfFactory.h"

using namespace RooFit ;
using namespace RooStats ;

// Enables residuals on fits
#define DrawResiduals

// fix Upsi mass for the fit
#define fixmass

const float upsi_pdg_mass = 9.46030;
const float z_pdg_mass = 91.1876;

void FIT();
void fit() { FIT(); }
void FIT() {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  double Z_right_limit = 116.;
  double Z_left_limit  = 66.;
  double UpsiLeftCut = 8.5;
  double UpsiRightCut = 11.;

  TF1 *f_straighline = new TF1("f_straighline", "0", 0, 1000);
  f_straighline->SetLineColor(kBlack);
  TF1 *f_1 = new TF1("f_1", "1",  0, Z_right_limit); f_1->SetLineColor(kBlue); f_1->SetLineStyle(7);
  TF1 *f_2 = new TF1("f_2", "-1", 0, Z_right_limit); f_2->SetLineColor(kBlue); f_2->SetLineStyle(7);

  int ZBins = (int)(Z_right_limit - Z_left_limit)/2.5;
  int Upsi_bins = (int)(UpsiRightCut - UpsiLeftCut)/0.05;

  TFile *ntuple = new TFile("file:/lstore/cms/eliza/ntupleMass/outTreeToFit_HZtoUpsilonPhoton_Data_ZtoUpsilon_Cat0_ZZZZZZ.root");
  TTree* tree = (TTree*) ntuple->Get("default/outTree_HZtoUpsilonPhoton");

  RooRealVar mHZ("mHZ","#mu^{+}#mu^{-}#gamma invariant mass [GeV]", Z_left_limit, Z_right_limit);
  RooRealVar mMuMNU("mMuMNU", "m_{#mu^{+}#mu^{-}} [GeV]", UpsiLeftCut, UpsiRightCut);

  RooArgSet Variables(mHZ, mMuMNU);

  RooDataSet *newdata = new RooDataSet("newdata", "newdata", tree, Variables);

  // Upsi 
#ifdef fixmass
	RooRealVar mean_m("mean_m","mean of gaussian",upsi_pdg_mass);
#else
	RooRealVar mean_m("mean_m","mean of gaussian",upsi_pdg_mass,upsi_pdg_mass-.20,upsi_pdg_mass+.20);
#endif
	RooRealVar sigma_m("sigma_m","Scale Factor 1", 3.4e-2, 0.1e-2, 10e-1);
	RooGaussian mSig ("mSig","signal p.d.f.", mMuMNU, mean_m, sigma_m);

  RooRealVar       c0("c0", "c0", -.732850, -2., 0.);
  RooExponential mBkg("mBkg", "exponential", mMuMNU, c0);

  // Z 
#ifdef fixmass
  RooRealVar  mean  ("mean" , "Control mean of gaussian", z_pdg_mass);//, 89.0, 93.0);
#else
  RooRealVar  mean  ("mean" , "Control mean of gaussian", z_pdg_mass, 89.0, 93.0);
#endif
  RooRealVar  sigma ("sigma", "Control sigma", 2., .1, 4.); // 1.6
  RooRealVar  resolution_sigma("resolution_sigma","sigma of the resolution function", 1.8, 1.7, 2.6);
  RooVoigtian  gauss ("gauss", "first gaussian PDF",  mHZ, mean, sigma, resolution_sigma);

  RooRealVar   c1("c1","coefficient #1", -.0276343, -.2, 0.) ; // 1.1
  RooExponential expo("expo","background p.d.f.", mHZ, c1) ;

  RooProdPdf OniaBckgZBckg("OniaBckgZBckg", "OniaBckgZBckg", RooArgSet(mBkg, expo)); 
  RooProdPdf OniaBckgZSig ("OniaBckgZSig",  "OniaBckgZSig ", RooArgSet(mBkg, gauss)); 
  RooProdPdf OniaSigZBckg ("OniaSigZBckg",  "OniaSigZBckg ", RooArgSet(mSig, expo)); 
  RooProdPdf OniaSigZSig  ("OniaSigZSig",   "OniaSigZSig  ", RooArgSet(mSig, gauss));  

  RooRealVar NOniaBckgZBckg("NOniaBckgZBckg", "Upsi background - Z background yield", 10,  0, 10000); 
  RooRealVar NOniaBckgZSig ("NOniaBckgZSig",  "Upsi background - Z signal yield", 30,  0, 10000); 
  RooRealVar NOniaSigZBckg ("NOniaSigZBckg",  "Upsi signal - Z background yield", 10,  0, 10000); 
  RooRealVar NOniaSigZSig  ("NOniaSigZSig",   "Upsi signal - Z signal yield", 100, 0, 200); 

  RooExtendPdf eOniaBckgZBckg("eOniaBckgZBckg", "eOniaBckgZBckg", OniaBckgZBckg,  NOniaBckgZBckg); 
  RooExtendPdf eOniaBckgZSig ("eOniaBckgZSig",  "eOniaBckgZSig ", OniaBckgZSig ,  NOniaBckgZSig );
  RooExtendPdf eOniaSigZBckg ("eOniaSigZBckg",  "eOniaSigZBckg ", OniaSigZBckg ,  NOniaSigZBckg );
  RooExtendPdf eOniaSigZSig  ("eOniaSigZSig",   "eOniaSigZSig  ", OniaSigZSig  ,  NOniaSigZSig  );

  RooAddPdf simPdf ("simPdf", "simPdf", RooArgList(eOniaBckgZBckg, eOniaBckgZSig, eOniaSigZBckg, eOniaSigZSig));

  RooFitResult *fr = simPdf.fitTo(*newdata, NumCPU(4, kTRUE), Save(), Extended());

  RooPlot* frame = mHZ.frame(Title("Z Mass Fit"), Bins(ZBins), Range(Z_left_limit, Z_right_limit)) ;
  newdata->plotOn(frame, XErrorSize(0), Name("ZPlotData"));
  simPdf.plotOn(frame, Range(Z_left_limit, Z_right_limit), Name("ZPlotModel"));
  simPdf.plotOn(frame, Components(expo),LineStyle(kDashed),LineColor(kBlue), Name("Bckg"), Name("ZUpsiBackgroundLine"));

  double totalSignal = NOniaSigZSig.getVal();
  double totalBackground = NOniaBckgZSig.getVal() + NOniaSigZSig.getVal();
  double NormFactor  = totalSignal / totalBackground;
  double ZNormFactor = NOniaBckgZSig.getVal() / totalBackground; 
	simPdf.plotOn(frame, Components(gauss), LineColor(kGreen), DrawOption("F"), FillColor(kGreen), Name("SignalGaus2"),Normalization(NormFactor));
	simPdf.plotOn(frame, Components(gauss), LineColor(kOrange), Name("SignalGaus"),Normalization(ZNormFactor));
  newdata->plotOn(frame, XErrorSize(0), Name("ZPlotData"));
	frame->GetYaxis()->SetTitle("Events / 2.5 GeV");

	frame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}#gamma} [GeV]");

  RooPlot* dummy_frame_Z = mHZ.frame(Title("dummy frame to extract residuals"), Bins(ZBins), Range(Z_left_limit, Z_right_limit));
  newdata->plotOn(dummy_frame_Z,XErrorSize(0)); 
  simPdf.plotOn(dummy_frame_Z);

  RooHist* h_residuals_mass = dummy_frame_Z->pullHist();
  RooPlot* frame_residuals_mass_Z = mHZ.frame(Bins(ZBins), Range(Z_left_limit, Z_right_limit));
  frame_residuals_mass_Z->GetYaxis()->SetTitle("pulls"); 
  frame_residuals_mass_Z->GetYaxis()->SetTitleSize(0.34);
  frame_residuals_mass_Z->GetYaxis()->SetLabelSize(.30);
  frame_residuals_mass_Z->GetYaxis()->SetTitleOffset(.2);
  frame_residuals_mass_Z->GetYaxis()->SetNdivisions(5);
  frame_residuals_mass_Z->addPlotable(h_residuals_mass, "P");

  frame->GetXaxis()->SetTitleOffset(1.);

  TCanvas *mass_projection = new TCanvas("mass_projection", "Z boson", 900, 900); mass_projection->cd();
#ifdef DrawResiduals
  TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.05,1.0,1.0,21);                                       
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.1,22);                                       
  pad1->SetFillColor(0); pad2->SetFillColor(0);
  pad1->Draw(); 
	pad2->Draw(); 
  pad2->cd(); frame_residuals_mass_Z->Draw(); f_straighline->Draw("same"); f_1->Draw("same"); f_2->Draw("same");
  pad1->cd(); 
#endif
  frame->SetMinimum(1e-5);
  frame->Draw();

  //Plots data on to frame
  RooPlot* mFrame = mMuMNU.frame(Bins(Upsi_bins), Range(UpsiLeftCut,UpsiRightCut));
  newdata->plotOn(mFrame, Name("PlotData"),XErrorSize(0));

  //Plots full model, prompt and non-prompt models to frame
  simPdf.plotOn(mFrame, LineColor(kBlue),LineWidth(3),NumCPU(1,kTRUE), Name("PlotModel"));
  simPdf.plotOn(mFrame, Components("mBkg"), LineWidth(3),LineColor(kAzure),LineStyle(2),Name("UpsiBackgroundLine"));
  simPdf.plotOn(mFrame, Components("mSig"), LineWidth(3),DrawOption("F"),FillColor(kGreen),LineColor(kGreen),Normalization(NOniaSigZSig.getVal()/(NOniaSigZBckg.getVal()+NOniaSigZSig.getVal())), Name("UpsiSignalZSignalLine"));
  simPdf.plotOn(mFrame, Components("mSig"), LineWidth(3),LineColor(kOrange),Normalization(NOniaSigZBckg.getVal()/(NOniaSigZBckg.getVal()+NOniaSigZSig.getVal())),Name("UpsiSignalZBackgroundLine"));
  newdata->plotOn(mFrame,XErrorSize(0));

  RooPlot* dummy_frame_jpsi = mMuMNU.frame(Title("dummy frame to extract residuals"), Bins(Upsi_bins), Range(UpsiLeftCut,UpsiRightCut));
  newdata->plotOn(dummy_frame_jpsi,XErrorSize(0)); 
  simPdf.plotOn(dummy_frame_jpsi);

  RooHist* h_residuals_mass_jpsi = dummy_frame_jpsi->pullHist();
  RooPlot* frame_residuals_mass_jpsi = mMuMNU.frame(Title("Residual Distribution #mu^{+}#mu^{-} mass"), Range(UpsiLeftCut,UpsiRightCut));
	frame_residuals_mass_jpsi->GetYaxis()->SetTitle("pulls"); 
	frame_residuals_mass_jpsi->GetYaxis()->SetTitleSize(.34);
  frame_residuals_mass_jpsi->GetYaxis()->SetLabelSize(.30);
	frame_residuals_mass_jpsi->GetYaxis()->SetTitleOffset(.2);
	frame_residuals_mass_jpsi->GetYaxis()->SetNdivisions(5);
  frame_residuals_mass_jpsi->addPlotable(h_residuals_mass_jpsi, "P");

  mFrame->GetXaxis()->SetTitleOffset(1.);
  mFrame->GetYaxis()->SetTitle("Events / 50 MeV");

  // Creates and fills canvas with plot and info 
  TCanvas *canvas1 = new TCanvas("MassFit_", "Upsilon", 900, 900); canvas1->cd();
#ifdef DrawResiduals
  TPad *pad1_jpsi = new TPad("pad1_jpsi", "The pad 80% of the height",0.0,0.05,1.0,1.0,21);
  TPad *pad2_jpsi = new TPad("pad2_jpsi", "The pad 20% of the height",0.0,0.0,1.0,0.1,22);
  pad1_jpsi->Draw(); pad2_jpsi->Draw();
  pad1_jpsi->SetFillColor(0); pad2_jpsi->SetFillColor(0);
  pad2_jpsi->cd();
  frame_residuals_mass_jpsi->Draw(); f_straighline->Draw("same"); f_1->Draw("same"); f_2->Draw("same");
  pad1_jpsi->cd();
#endif

  mFrame->SetMinimum(1e-5);
  mFrame->Draw();
  mFrame->SetTitle("upsi Mass [GeV]");

}
