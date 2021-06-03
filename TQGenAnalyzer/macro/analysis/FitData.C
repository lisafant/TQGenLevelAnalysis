#include "TH1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
using namespace RooFit;

void fit_data(std::string CR_input){

//open trees
  TFile *fSR = new TFile("fout_mRun2018_A_B_C_D_trigger_TQ.root", "READ");
  TTree *treeSR = (TTree*)fSR->Get("tree_red");
  TFile *fCR = new TFile(("fout_mRun2018_A_B_C_D_trigger_"+CR_input+".root").c_str());  
  TTree *treeCR = (TTree*)fCR->Get("tree_red");
  TFile *fSPS = new TFile("fout_mSPS_TQ.root");
  TTree *treeSPS = (TTree*)fSPS->Get("tree_red");
  //TFile *f26 = new TFile("fout_m26_TQ.root");
  //TTree *tree26 = (TTree*)f26->Get("tree_red");

  TFile *f26 = new TFile("workspace_fit26.root", "READ");
  RooWorkspace *w;
  f26->GetObject("w", w);
  w->Print();
    

  Long64_t nentriesCR = (Long64_t)treeCR->GetEntries();
  Long64_t nentriesSPS = (Long64_t)treeSPS->GetEntries();
  
  RooRealVar TQ_mass_tilde("TQ_mass_tilde", "TQ_mass_tilde", 0, 60);
  //RooAbsReal TQ_mass_tilde_func("TQ_mass_tilde_func", "TQ_mass_tilde_func", "TQ_mass_tilde");

  //RooDataSet SR("SR", "SR", treeSR, RooArgSet(TQ_mass_tilde));
  RooDataSet CR("CR", "CR", treeCR, RooArgSet(TQ_mass_tilde));
  RooDataSet SPS("SPS", "SPS", treeSPS, RooArgSet(TQ_mass_tilde));
  
  RooDataSet* sig_data = w->pdf("pdf3")->generate(TQ_mass_tilde, 1200);
  RooKeysPdf signal_pdf("signal_pdf", "signal_pdf", TQ_mass_tilde, *sig_data);
  RooKeysPdf CR_pdf("CR_pdf", "CR_pdf", TQ_mass_tilde,CR);
  RooRealVar n_CR("n_CR", "n_CR", nentriesCR);
  RooRealVar n_sig("n_sig", "n_sig", 1200.); 
  RooAddPdf model("model", "", RooArgList(signal_pdf, CR_pdf), RooArgList(n_CR, n_sig));
 
  model.fitTo(CR); 
  RooPlot* CRegion = TQ_mass_tilde.frame();
  //CR.plotOn(CRegion);
  CR_pdf.plotOn(CRegion, LineColor(kRed));
  signal_pdf.plotOn(CRegion, LineColor(kGreen));
  //w->pdf("pdf3")->plotOn(CRegion, LineColor(kGreen));
  model.plotOn(CRegion);
  TCanvas* ContRegion = new TCanvas("ContRegion", "", 1500, 1000);
  CRegion->GetXaxis()->SetTitleOffset(0.9);
  CRegion->GetXaxis()->SetTitle("#tilde{m}_{TQ}");
  CRegion->GetYaxis()->SetTitle("Events / 0.6 GeV");
  CRegion->Draw();
  ContRegion->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/Fit.png").c_str());
  delete ContRegion;
  
}
