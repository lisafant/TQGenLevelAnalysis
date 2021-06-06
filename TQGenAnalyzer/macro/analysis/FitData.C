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
#include "RooFFTConvPdf.h"
#include "RooExtendPdf.h"
#include "RooBreitWigner.h"
using namespace RooFit;

void fit_data(std::string CR_input){

//open trees
  TFile *fSR = new TFile("fout_mRun2018_A_B_C_D_trigger_TQ.root", "READ");
  TTree *treeSR = (TTree*)fSR->Get("tree_red");
  TFile *fCR_tree = new TFile(("fout_mRun2018_A_B_C_D_trigger_"+CR_input+".root").c_str(), "READ");  
  TTree *treeCR = (TTree*)fCR_tree->Get("tree_red");
  TFile *fSPS = new TFile("fout_mSPS_TQ.root", "READ");
  TTree *treeSPS = (TTree*)fSPS->Get("tree_red");
  TFile *f26 = new TFile("fout_m26_TQ.root", "READ");
  TTree *tree26 = (TTree*)f26->Get("tree_red");

  TFile *fworksp_26 =  new TFile("workspace_fit26.root", "READ");
  RooWorkspace *w;
  fworksp_26->GetObject("w", w);
  w->Print();
 
  TFile *fworksp_CR = new TFile("workspace_fitCR.root", "READ");
  RooWorkspace *w_CR;
  fworksp_CR->GetObject("w_CR", w_CR);
  w_CR->Print();

  TFile *fworksp_SPS = new TFile("workspace_fitSPS.root", "READ");   
  RooWorkspace *w_SPS;
  fworksp_SPS->GetObject("w_SPS", w_SPS);
  w_SPS->Print();

  Long64_t nentriesCR = (Long64_t)treeCR->GetEntries();
  Long64_t nentriesSPS = (Long64_t)treeSPS->GetEntries();
  Long64_t nentries26 = (Long64_t)tree26->GetEntries();
  cout<<"N. entries CR: "<<nentriesCR<<endl;
  cout<<"N. entries SPS: "<<nentriesSPS<<endl;
  cout<<"N. entries signal 26 GeV: "<<nentries26<<endl;
  
  RooRealVar TQ_mass_tilde("TQ_mass_tilde", "TQ_mass_tilde", 0, 60);
 
  RooDataSet CR("CR", "CR", treeCR, RooArgSet(TQ_mass_tilde));
  cout<<CR.sumEntries()<<endl;
  cout<<CR.numEntries()<<endl; 
 
  RooAbsPdf* signal_pdf = w->pdf("pdf3");
  RooAbsPdf* CR_pdf = w_CR->pdf("model_CR");
  RooAbsPdf* SPS_pdf = w_SPS->pdf("sigmoid_pdf");

  // RooKeysPdf signal_pdf("signal_pdf", "signal_pdf", TQ_mass_tilde, *sig_data);
  //RooKeysPdf CR_pdf("CR_pdf", "CR_pdf", TQ_mass_tilde,CR);
  //RooKeysPdf SPS_pdf("SPS_pdf", "SPS_pdf", TQ_mass_tilde, SPS);
 
  RooRealVar n_CR("n_CR", "n_CR", nentriesCR, 0.5*nentriesCR, 2*nentriesCR);
  RooRealVar n_sig("n_sig", "n_sig", nentries26, 0., 4*nentries26);
  RooRealVar n_SPS("n_SPS", "n_SPS", 1500.);

  RooExtendPdf e_CR_pdf("e_CR_pdf", "", *CR_pdf, n_CR);
  RooExtendPdf e_signal_pdf("e_signal_pdf", "", *signal_pdf, n_sig);
  RooExtendPdf e_SPS_pdf("e_SPS_pdf", "", *SPS_pdf, n_SPS); 
  RooAddPdf e_model("e_model", "e_signal_pdf+e_CR_pdf+e_SPS_pdf", RooArgList(e_signal_pdf, e_CR_pdf, e_SPS_pdf));
 

  //TQ_mass_tilde.setRange("range", 9, 40);
  e_model.fitTo(CR, Save()); 

  RooPlot* CRegion = TQ_mass_tilde.frame();
  RooPlot* pdfs = TQ_mass_tilde.frame();
  CR.plotOn(CRegion);
  e_model.plotOn(CRegion);
  e_model.plotOn(CRegion, Components("e_CR_pdf"), LineColor(kRed));
  e_model.plotOn(CRegion, Components("e_signal_pdf"), LineColor(kGreen));
  e_model.plotOn(CRegion, Components("e_SPS_pdf"), LineColor(kMagenta));
  e_signal_pdf.plotOn(pdfs, LineColor(kGreen));
  e_SPS_pdf.plotOn(pdfs, LineColor(kMagenta));
  e_CR_pdf.plotOn(pdfs, LineColor(kRed));
  TCanvas* ContRegion = new TCanvas("ContRegion", "", 1500, 1000);
  ContRegion->Divide(1, 2);
  ContRegion->cd(1);
  CRegion->GetXaxis()->SetTitleOffset(0.9);
  CRegion->GetXaxis()->SetTitle("#tilde{m}_{TQ}");
  CRegion->GetYaxis()->SetTitle("Events / 0.6 GeV");
  CRegion->Draw();
  ContRegion->cd(2);
  pdfs->GetXaxis()->SetTitle("#tilde{m}_{TQ}");
  pdfs->GetYaxis()->SetTitle("p.d.f.");
  pdfs->Draw();
  ContRegion->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/Fit.png").c_str());
  delete ContRegion;
  
}
