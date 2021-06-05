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
  TFile *fCR_tree = new TFile(("fout_mRun2018_A_B_C_D_trigger_"+CR_input+".root").c_str());  
  TTree *treeCR = (TTree*)fCR_tree->Get("tree_red");
  //TFile *fSPS = new TFile("fout_mSPS_TQ.root");
  //TTree *treeSPS = (TTree*)fSPS->Get("tree_red");
  //TFile *f26 = new TFile("fout_m26_TQ.root");
  //TTree *tree26 = (TTree*)f26->Get("tree_red");

  TFile *f26 = new TFile("workspace_fit26.root", "READ");
  RooWorkspace *w;
  f26->GetObject("w", w);
  w->Print();
 
  TFile *fCR = new TFile("workspace_fitCR.root", "READ");
  RooWorkspace *w_CR;
  fCR->GetObject("w_CR", w_CR);
  w_CR->Print();

  TFile *fSPS = new TFile("workspace_fitSPS.root", "READ");   
  RooWorkspace *w_SPS;
  fSPS->GetObject("w_SPS", w_SPS);
  w_SPS->Print();

  Long64_t nentriesCR = (Long64_t)treeCR->GetEntries();
//  Long64_t nentriesSPS = (Long64_t)treeSPS->GetEntries();
  
  RooRealVar TQ_mass_tilde("TQ_mass_tilde", "TQ_mass_tilde", 10, 60);
  //RooAbsReal TQ_mass_tilde_func("TQ_mass_tilde_func", "TQ_mass_tilde_func", "TQ_mass_tilde");

  //RooDataSet SR("SR", "SR", treeSR, RooArgSet(TQ_mass_tilde));
  RooDataSet CR("CR", "CR", treeCR, RooArgSet(TQ_mass_tilde));
  //RooDataSet SPS("SPS", "SPS", treeSPS, RooArgSet(TQ_mass_tilde));
  //cout<<CR.sumEntries()<<endl;
  //cout<<CR.numEntries()<<endl; 
  
  RooAbsPdf* signal_pdf = w->pdf("pdf3");
  RooAbsPdf* CR_pdf = w_CR->pdf("cheb");
  RooAbsPdf* SPS_pdf = w_SPS->pdf("sigmoid_pdf");
 // RooKeysPdf signal_pdf("signal_pdf", "signal_pdf", TQ_mass_tilde, *sig_data);
  //RooKeysPdf CR_pdf("CR_pdf", "CR_pdf", TQ_mass_tilde,CR);
  //RooKeysPdf SPS_pdf("SPS_pdf", "SPS_pdf", TQ_mass_tilde, SPS);
  RooRealVar n_CR("n_CR", "n_CR", nentriesCR, 0.5*nentriesCR, 1.5*nentriesCR);
  RooRealVar n_sig("n_sig", "n_sig", 200., 0., 1200.);
  RooRealVar n_SPS("n_SPS", "n_SPS", 231.);
  RooExtendPdf e_CR_pdf("e_CR_pdf", "", *CR_pdf, n_CR);
  RooExtendPdf e_signal_pdf("e_signal_pdf", "", *signal_pdf, n_sig);
  RooExtendPdf e_SPS_pdf("e_SPS_pdf", "", *SPS_pdf, n_SPS); 
  RooAddPdf e_model("e_model", "e_signal_pdf+e_CR_pdf+e_SPS_pdf", RooArgList(e_signal_pdf, e_CR_pdf, e_SPS_pdf));
 
  e_model.fitTo(CR); 
  RooPlot* CRegion = TQ_mass_tilde.frame();
  RooPlot* prova = TQ_mass_tilde.frame();
  CR.plotOn(CRegion);
  //CR_pdf.plotOn(CRegion, LineColor(kRed));
  //signal_pdf->plotOn(CRegion, LineColor(kGreen));
  //w->pdf("pdf3")->plotOn(CRegion, LineColor(kMagenta));
  //#include "RooExtendPdf.h"On(CRegion, LineColor(kGreen));
  e_model.plotOn(CRegion);
  e_model.plotOn(prova, Components("e_signal_pdf"), LineColor(kGreen));
  e_model.plotOn(CRegion, Components("e_CR_pdf"), LineColor(kRed));
  e_model.plotOn(CRegion, Components("e_signal_pdf"), LineColor(kGreen));
  e_model.plotOn(CRegion, Components("e_SPS_pdf"), LineColor(kMagenta));
  TCanvas* ContRegion = new TCanvas("ContRegion", "", 1500, 1000);
  ContRegion->Divide(1, 2);
  ContRegion->cd(1);
  ContRegion->SetLogy();
  //ContRegion->DrawFrame(0, 0, 60, 0.0005);
  CRegion->GetXaxis()->SetTitleOffset(0.9);
  CRegion->GetXaxis()->SetTitle("#tilde{m}_{TQ}");
  CRegion->GetYaxis()->SetTitle("Events / 0.6 GeV");
  CRegion->Draw();
  ContRegion->cd(2);
  prova->Draw();
  ContRegion->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/Fit.png").c_str());
  delete ContRegion;
  
}
