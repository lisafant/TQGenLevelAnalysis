#include "TH1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"
#include "RooArgSet.h"
#include "RooChebychev.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooKeysPdf.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooAddition.h"
#include "RooCBShape.h"
#include "RooLandau.h"
using namespace RooFit;


RooFitResult* Fit_SPS(){

  TFile *fSPS = new TFile("fout_mSPS_TQ.root", "READ");
  TTree* treeSPS = (TTree*)fSPS->Get("tree_red");

  Long64_t nentries = (Long64_t)treeSPS->GetEntries();
  cout<<"N. entries SPS: "<<nentries<<endl;

  RooRealVar TQ_mass_tilde("TQ_mass_tilde", "TQ_mass_tilde", 13, 30);
//  RooRealVar weight("weight", "weight", 0.004);
  RooDataSet *data = new RooDataSet("data", "data", treeSPS, RooArgSet(TQ_mass_tilde));
  
//  RooFormulaVar wFunc("weig","event weight","0.004", TQ_mass_tilde);
//  RooRealVar* weig = (RooRealVar*) data->addColumn(wFunc);
//  RooDataSet wdata(data->GetName(),data->GetTitle(),data,*data->get(),0,weig->GetName());
 

  RooRealVar mean("mean", "mean", 21., 18., 24.);
  RooRealVar sigma("sigma", "sigma", 0.5, 0.05, 5.);
  RooRealVar alpha("alpha", "alpha", 0.8, 0.1, 5);
  RooRealVar n("n", "n", 1.8, 0.1, 20.);
  //RooRealVar mean_pos("mean_pos", "mean_pos", 21., 18., 24.);
  RooRealVar sigma_pos("sigma_pos", "sigma_pos", 0.5, 0.05, 5.);
  RooRealVar alpha_pos("alpha_pos", "alpha_pos", -0.2, -5., -0.01);
  RooRealVar n_pos("n_pos", "n_pos", 1.8, 0.1, 20.);
  RooRealVar fraction("fraction", "fraction", 0.8, 0., 1.);

  RooCBShape cryst("cryst", "cryst", TQ_mass_tilde, mean, sigma, alpha, n);
  RooCBShape cryst_pos("cryst_pos", "cryst_ps", TQ_mass_tilde, mean, sigma_pos, alpha_pos, n_pos);
  RooAddPdf pdfSPS("pdfSPS", "pdfSPS", RooArgList(cryst, cryst_pos), fraction);

  RooFitResult *r = pdfSPS.fitTo(*data, Save());



 

  RooPlot* plotSPS = TQ_mass_tilde.frame();
  plotSPS->SetTitle( "SPS" );
  plotSPS->GetXaxis()->SetTitle("#tilde{m}_{TQ} [GeV]");
  plotSPS->GetYaxis()->SetTitleOffset( 1.5 );
  plotSPS->GetYaxis()->SetTitle("Events / 0.3 GeV");
  data->plotOn(plotSPS, Name ("MC")); 
  pdfSPS.plotOn(plotSPS, Name ("fit"));
  pdfSPS.plotOn(plotSPS, Name ("cryst1"), Components("cryst"), LineColor(kBlack));
  pdfSPS.plotOn(plotSPS, Name ("cryst2"), Components("cryst_pos"), LineColor(kRed));
  TCanvas* canvSPS = new TCanvas("canvSPS", "canvSPS", 1500, 1000);
  plotSPS->Draw();
  TLegend *legend_SPS = new TLegend(0.20, 0.70, 0.40, 0.90);
  legend_SPS->AddEntry("fit", "CB1 + CB2 Fit", "l");
  legend_SPS->AddEntry("cryst1", "Crystal Ball 1", "l");
  legend_SPS->AddEntry("cryst2", "Crystal Ball 2", "l");
  legend_SPS->SetBorderSize(0);
  legend_SPS->Draw();
  canvSPS->SaveAs("/eos/home-l/lfantini/www/Background/CR1a/SPSFit.png");
  canvSPS->SaveAs("/eos/home-l/lfantini/www/Background/CR1b/SPSFit.png");
  canvSPS->SaveAs("/eos/home-l/lfantini/www/Background/CR2/SPSFit.png");  
  canvSPS->SaveAs("/eos/home-l/lfantini/www/Background/CR3/SPSFit.png");

  TQ_mass_tilde.setConstant(kTRUE);
  mean.setConstant(kTRUE);
  sigma.setConstant(kTRUE);
  alpha.setConstant(kTRUE);
  n.setConstant(kTRUE);
  sigma_pos.setConstant(kTRUE);
  alpha_pos.setConstant(kTRUE);
  n_pos.setConstant(kTRUE);
  

  RooWorkspace *w_SPS = new RooWorkspace("w_SPS", "workspace_SPS"); 
  w_SPS->import(pdfSPS);
  w_SPS->writeToFile("workspace_fitSPS.root");
  gDirectory->Add(w_SPS);
  

  delete canvSPS;
  delete data;
  delete legend_SPS;
  return r;
}





RooFitResult* Fit_CR(std::string CR_input){

  TFile *fCR = new TFile(("fout_mRun2018_A_B_C_D_trigger_"+CR_input+".root").c_str());
  TTree* treeCR = (TTree*)fCR->Get("tree_red");

  Long64_t nentries = (Long64_t)treeCR->GetEntries();
  cout<<"N. entries CR: "<<nentries<<endl;

  RooRealVar TQ_mass_tilde("TQ_mass_tilde", "TQ_mass_tilde", 13, 30);
 
  RooDataSet *data_CR = new RooDataSet("data_CR", "data_CR", treeCR, RooArgSet(TQ_mass_tilde));

  Double_t n_data_CR =data_CR->sumEntries("TQ_mass_tilde>13. && TQ_mass_tilde<30.");
  cout<<"N. Entries: "<<n_data_CR<<endl;


  RooRealVar a0("a0", "a0", 0.1, -0.5, 1.);
  RooRealVar a1("a1", "a1", -0.1, -0.5, 1.);
  RooRealVar a2("a2", "a2", -0.1, -0.5, 1.);
  RooRealVar a3("a3", "a3", -0.1, -0.5, 1.);
  RooRealVar a4("a4", "a4", -0.1, -0.5, 1.);

  RooChebychev cheb("cheb", "cheb", TQ_mass_tilde, RooArgSet(a0, a1, a2, a3, a4));
  
  RooFitResult *r = cheb.fitTo(*data_CR, Save());


  
  RooPlot* plotCR = TQ_mass_tilde.frame();
  plotCR->SetTitle( "CR" );
  plotCR->GetXaxis()->SetTitle("#tilde{m}_{TQ} [GeV]");
  plotCR->GetYaxis()->SetTitleOffset( 1.5 );
  plotCR->GetYaxis()->SetTitle("Events / 0.5 GeV");
  data_CR->plotOn(plotCR);
  cheb.plotOn(plotCR, Name ("fit"));
  TCanvas* canvCR = new TCanvas("canvCR", "canvCR", 1500, 1000);
  plotCR->Draw();
  TLegend *legend_CR = new TLegend(0.20, 0.70, 0.40, 0.90);
  legend_CR->AddEntry("fit", "Chebychev(5) Fit", "l");
  legend_CR->SetBorderSize(0);
  legend_CR->Draw();
  canvCR->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/CRFit.png").c_str());  

  TQ_mass_tilde.setConstant(kTRUE);
  a0.setConstant(kTRUE);
  a1.setConstant(kTRUE);
  a2.setConstant(kTRUE);
  a3.setConstant(kTRUE);
  a4.setConstant(kTRUE);
  

  RooWorkspace *w_CR = new RooWorkspace("w_CR", "workspace_CR");
  w_CR->import(cheb);
  w_CR->writeToFile("workspace_fitCR.root");
  gDirectory->Add(w_CR);

  delete canvCR;
  delete data_CR;
  return r;
}
























































