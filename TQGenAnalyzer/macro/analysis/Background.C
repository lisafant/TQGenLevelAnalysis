#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TAxis.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TMath.h"
#include "THStack.h"
#include <iostream>
using namespace std;



//Function for SPS weight
double WeightSPS(Long64_t N_gen, double xsec_tot, double lum_mis){
  double BR_e =0.0238;
  double BR_mu=0.0248;
  //double BR = BR_mu + BR_e;
  double BR = TMath::Power((BR_e + BR_mu), 2);
  double xsec = xsec_tot*BR;
  double lum_eq = (double)N_gen/xsec;
  double weight = lum_mis/lum_eq;
  return (weight);
}




void GraphBG(std::string CR_input){
  
///////////TREES//////////


//Signal region
  TFile *fTQ = new TFile("fout_mRun2018_A_B_C_D_trigger_TQ.root", "READ");
  TTree *treeTQ = (TTree*)fTQ->Get("tree_red");

//Control region
  TFile *fCR = new TFile(("fout_mRun2018_A_B_C_D_trigger_"+CR_input+".root").c_str(), "READ");
  TTree *treeCR = (TTree*)fCR->Get("tree_red");
 
//SPS
  TFile *fSPS = new TFile("fout_mSPS_TQ.root", "READ");
  TTree *treeSPS = (TTree*)fSPS->Get("tree_red");


  Long64_t nentriesSR = (Long64_t)treeTQ->GetEntries();
  Long64_t nentriesCR = (Long64_t)treeCR->GetEntries();
  Long64_t nentriesSPS = (Long64_t)treeSPS->GetEntries();
  cout<<"N entries SR: "<<nentriesSR<<endl;
  cout<<"N entries CR: "<<nentriesCR<<endl;
  cout<<"N entries SPS: "<<nentriesSPS<<endl;

  double lum_mis = 59.97; // fb^-1 
  TLeaf *xsec = treeSPS->GetLeaf("x_sec");
  treeSPS->GetEntry(2);
  Double_t cross_sec= (xsec->GetValue())*1000.; //fb 
  cout<<"Cross section SPS: "<<cross_sec<<endl;

  double weightSPS = WeightSPS(nentriesSPS, cross_sec, lum_mis)*0.15*nentriesSPS;
  cout<<"Weight lum for SPS: "<<weightSPS<<endl;

/////////STATISTICAL TESTS////////
  THStack *test = new THStack("test", "Mass tilde stat test; #tilde{m}_{TQ} [GeV]; Events / 0.21 GeV");
  
  TH1F *SR_test = new TH1F("SR_test", "", 500, 13, 28);
  TH1F *CR_test = new TH1F("CR_test", "", 500, 13, 28);
 
  Float_t mass_tilde_SR;
  Float_t mass_tilde_CR;
  treeTQ->SetBranchAddress("TQ_mass_tilde", &mass_tilde_SR);
  treeCR->SetBranchAddress("TQ_mass_tilde", &mass_tilde_CR);

  for(int i = 0; i<nentriesSR; i++){
    treeTQ->GetEntry(i);
    if (mass_tilde_SR>17. && mass_tilde_SR<19.) continue;
    SR_test->Fill(mass_tilde_SR);
  }

  for(int j = 0; j<nentriesCR; j++){
    treeCR->GetEntry(j);
    if (mass_tilde_CR>17. && mass_tilde_CR<19.) continue;
    CR_test->Fill(mass_tilde_CR);
  }
  
  SR_test->Sumw2();
  CR_test->Sumw2();

  double int_SR_test = SR_test->Integral();        
  CR_test->Scale(1./CR_test->Integral());
  SR_test->Scale(1./SR_test->Integral());

  SR_test->SetMarkerColor(1);
  CR_test->SetMarkerColor(2);

  test->Add(SR_test);
  test->Add(CR_test);

  Double_t KSTest = SR_test->KolmogorovTest(CR_test, "D");
  Double_t X2Test = SR_test->Chi2Test(CR_test, "NORM UU P");

  TCanvas* Test = new TCanvas("Test", "Hist for stat tests", 1500, 1000);
  test->Draw("E0 nostack");
  Test->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/test.png").c_str());
  delete test;
  delete Test;
  delete SR_test;
  delete CR_test;
   

//////////HISTOGRAMS//////////////

/////number of vertices
  THStack *n_vtx = new THStack("n_vtx", "Number of vertices; N vertices; Events / 1");

  TH1F *SR_n_vtx = new TH1F("SR_n_vtx", "", 60, 0, 60);
  treeTQ->Draw("n_vtx>>SR_n_vtx");
  TH1F *SPS_n_vtx = new TH1F("SPS_n_vtx", "", 60, 0, 60);
  treeSPS->Draw("n_vtx>>SPS_n_vtx");

  SR_n_vtx->Sumw2();
  SPS_n_vtx->Sumw2();

  SR_n_vtx->Scale(1./SR_n_vtx->Integral());
  SPS_n_vtx->Scale(1./SPS_n_vtx->Integral());


  SR_n_vtx->SetLineColor(1);
  SPS_n_vtx->SetLineColor(4);

  n_vtx->Add(SR_n_vtx);
  n_vtx->Add(SPS_n_vtx);

  TCanvas* n_Vtx = new TCanvas("n_Vtx", "Number of vertices", 1500, 1000);
  n_vtx->Draw("HIST nostack");
  TLegend* legend_nvtx = new TLegend(0.60, 0.78, 0.80, 0.88);
  legend_nvtx->AddEntry(SR_n_vtx, "Signal Region", "l");
  legend_nvtx->AddEntry(SPS_n_vtx, "SPS", "l");
  legend_nvtx->SetBorderSize(0);
  legend_nvtx->Draw();
  n_Vtx->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/n_vtx.png").c_str());
  delete n_vtx;
  delete n_Vtx;
  delete SR_n_vtx;
  delete SPS_n_vtx;
  delete legend_nvtx;
    


/////TQ mass
  THStack *TQ_mass = new THStack("TQ_mass", "Tetraquark mass; m_{TQ} [GeV]; Events / 0.6 GeV");

  TH1F *SR_TQ_mass = new TH1F("SR_TQ_mass", "TQ mass SR", 100, 0, 60);
  treeTQ->Draw("TQ_mass>>SR_TQ_mass");
  TH1F *CR_TQ_mass = new TH1F("CR_TQ_mass", "TQ mass CR", 100, 0, 60);
  treeCR->Draw("TQ_mass>>CR_TQ_mass");
  TH1F *SPS_TQ_mass = new TH1F("SPS_TQ_mass", "TQ mass SPS", 100, 0, 60);
  treeSPS->Draw("TQ_mass>>SPS_TQ_mass");
  //cout<<SPS_TQ_mass->Integral()<<endl;
  
 
  SR_TQ_mass->Sumw2();
  CR_TQ_mass->Sumw2();
  SPS_TQ_mass->Sumw2();

//cout<<SPS_TQ_mass->Integral()<<endl;
  double int_TQ_mass = SR_TQ_mass->Integral();
//cout<<int_TQ_mass<<endl;
  CR_TQ_mass->Scale(int_TQ_mass/CR_TQ_mass->Integral());  
   SPS_TQ_mass->Scale((weightSPS)/SPS_TQ_mass->Integral()); 
cout<<SPS_TQ_mass->Integral()<<endl;


  SR_TQ_mass->SetLineColor(1);
  CR_TQ_mass->SetLineColor(2);
  SPS_TQ_mass->SetLineColor(4);
  SR_TQ_mass->SetMarkerColor(1);
  CR_TQ_mass->SetMarkerColor(2);
  SPS_TQ_mass->SetMarkerColor(4);

  int nbinsx_mass = SR_TQ_mass->GetXaxis()->GetNbins();
  for(int j=0; j<nbinsx_mass; j++){
    double up = SR_TQ_mass->GetXaxis()->GetBinUpEdge(j);
    if(up>16. && up<20.){
      SR_TQ_mass-> SetBinContent(j, 0);
    }
  }
  

  TQ_mass->Add(SR_TQ_mass);
  TQ_mass->Add(CR_TQ_mass);
  TQ_mass->Add(SPS_TQ_mass);

  TCanvas *TQ_Mass = new TCanvas("TQ_Mass", "Tetraquark mass", 1500, 1000);
  TQ_mass->Draw("HIST nostack");
  TLegend* legend_TQ_mass = new TLegend(0.60, 0.78, 0.80, 0.88);
  legend_TQ_mass->AddEntry(SR_TQ_mass, "Signal Region", "l");
  legend_TQ_mass->AddEntry(CR_TQ_mass, "Control Region", "l");
  legend_TQ_mass->AddEntry(SPS_TQ_mass, "SPS", "l");
  legend_TQ_mass->SetBorderSize(0);
  legend_TQ_mass->Draw();
  TQ_Mass->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/massTQ.png").c_str());
  delete TQ_mass;
  delete TQ_Mass;
  delete SR_TQ_mass;
  delete CR_TQ_mass;
  delete SPS_TQ_mass;
  delete legend_TQ_mass;

  

//TQ mass tilde
  THStack *TQ_mass_tilde = new THStack("TQ_mass_tilde", "Tetraquark mass tilde; #tilde{m}_{TQ} [GeV]; Events / 0.6 GeV");

  TH1F *SR_TQ_mass_tilde = new TH1F("SR_TQ_mass_tilde", "TQ mass tilde SR", 100, 0, 60);
  treeTQ->Draw("TQ_mass_tilde>>SR_TQ_mass_tilde");
  TH1F *CR_TQ_mass_tilde = new TH1F("CR_TQ_mass_tilde", "TQ mass tilde CR", 100, 0, 60);
  treeCR->Draw("TQ_mass_tilde>>CR_TQ_mass_tilde");
  TH1F *SPS_TQ_mass_tilde = new TH1F("SPS_TQ_mass_tilde", "TQ mass tilde SPS", 100, 0, 60);
  treeSPS->Draw("TQ_mass_tilde>>SPS_TQ_mass_tilde");
  //cout<<SPS_TQ_mass_tilde->Integral()<<endl;
  //treeSPS->Draw("TQ_mass_tilde>>SPS_TQ_mass_tilde_prova(100, 0, 60)", "puw2018");
  //TH1F *SPS_TQ_mass_tilde_prova = (TH1F*)gDirectory->Get("SPS_TQ_mass_tilde_prova");
  //cout<<SPS_TQ_mass_tilde_prova->Integral()<<endl;
    
 
  SR_TQ_mass_tilde->Sumw2();
  CR_TQ_mass_tilde->Sumw2();
  SPS_TQ_mass_tilde->Sumw2();

  double int_TQ_mass_tilde = SR_TQ_mass_tilde->Integral();
  CR_TQ_mass_tilde->Scale(int_TQ_mass_tilde/CR_TQ_mass_tilde->Integral());  
  SPS_TQ_mass_tilde->Scale(weightSPS/SPS_TQ_mass_tilde->Integral()); 

  SR_TQ_mass_tilde->SetLineColor(1);
  CR_TQ_mass_tilde->SetLineColor(2);
  SPS_TQ_mass_tilde->SetLineColor(4);
  SR_TQ_mass_tilde->SetMarkerColor(1);
  CR_TQ_mass_tilde->SetMarkerColor(2);
  SPS_TQ_mass_tilde->SetMarkerColor(4);

  int nbinsx_mass_tilde = SR_TQ_mass_tilde->GetXaxis()->GetNbins();
  for(int j=0; j<nbinsx_mass_tilde; j++){
    double up_tilde = SR_TQ_mass_tilde->GetXaxis()->GetBinUpEdge(j);
    if(up_tilde>16. && up_tilde<20.){
      SR_TQ_mass_tilde-> SetBinContent(j, 0);
    }
  }
  

  TQ_mass_tilde->Add(SR_TQ_mass_tilde);
  TQ_mass_tilde->Add(CR_TQ_mass_tilde);
  TQ_mass_tilde->Add(SPS_TQ_mass_tilde);

  TCanvas *TQ_Mass_Tilde = new TCanvas("TQ_Mass_tilde", "Tetraquark mass tilde", 1500, 1000);
  TQ_mass_tilde->Draw("HIST nostack");
  //gPad->SetLogy();
  TLegend* legend_TQ_mass_tilde = new TLegend(0.60, 0.78, 0.80, 0.88);
  legend_TQ_mass_tilde->AddEntry(SR_TQ_mass_tilde, "Signal Region", "l");
  legend_TQ_mass_tilde->AddEntry(CR_TQ_mass_tilde, "Control Region", "l");
  legend_TQ_mass_tilde->AddEntry(SPS_TQ_mass_tilde, "SPS", "l");
  legend_TQ_mass_tilde->SetBorderSize(0);
  legend_TQ_mass_tilde->Draw();
  TQ_Mass_Tilde->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/massTQ_tilde.png").c_str());
  delete TQ_mass_tilde;
  delete TQ_Mass_Tilde;
  delete SR_TQ_mass_tilde;
  delete CR_TQ_mass_tilde;
  delete SPS_TQ_mass_tilde;
  delete legend_TQ_mass_tilde;

//TQ pt
  THStack *TQ_pt = new THStack("TQ_pt", "Tetraquark pt; p_{T}_{TQ} [GeV]; Events / 0.6 [GeV]");

  TH1F *SR_TQ_pt = new TH1F("SR_TQ_pt", "TQ pt SR", 100, 0, 60);
  treeTQ->Draw("TQ_pt>>SR_TQ_pt");
  TH1F *CR_TQ_pt = new TH1F("CR_TQ_pt", "TQ pt CR", 100, 0, 60);
  treeCR->Draw("TQ_pt>>CR_TQ_pt");  
  TH1F *SPS_TQ_pt = new TH1F("SPS_TQ_pt", "TQ pt SPS", 100, 0, 60);
  treeSPS->Draw("TQ_pt>>SPS_TQ_pt");

  SR_TQ_pt->Sumw2();
  CR_TQ_pt->Sumw2();
  SPS_TQ_pt->Sumw2();
                                                                                             
  double int_TQ_pt = SR_TQ_pt->Integral();
  CR_TQ_pt->Scale(int_TQ_pt/CR_TQ_pt->Integral());  
  SPS_TQ_pt->Scale(weightSPS/SPS_TQ_pt->Integral()); 
                                                                                             
  SR_TQ_pt->SetLineColor(1);
  CR_TQ_pt->SetLineColor(2);
  SPS_TQ_pt->SetLineColor(4);
  SR_TQ_pt->SetMarkerColor(1);
  CR_TQ_pt->SetMarkerColor(2);
  SPS_TQ_pt->SetMarkerColor(4);
                                                                                             
  TQ_pt->Add(SR_TQ_pt);
  TQ_pt->Add(CR_TQ_pt);
  TQ_pt->Add(SPS_TQ_pt);
                                                                                             
  TCanvas *TQ_Pt = new TCanvas("TQ_Pt", "Tetraquark pt", 1500, 1000);
  TQ_pt->Draw("HIST nostack");
  TLegend* legend_TQ_pt = new TLegend(0.60, 0.78, 0.80, 0.88);
  legend_TQ_pt->AddEntry(SR_TQ_pt, "Signal Region", "l");
  legend_TQ_pt->AddEntry(CR_TQ_pt, "Control Region", "l");
  legend_TQ_pt->AddEntry(SPS_TQ_pt, "SPS", "l");
  legend_TQ_pt->SetBorderSize(0);
  legend_TQ_pt->Draw();
  TQ_Pt->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/ptTQ.png").c_str());
  delete TQ_pt;
  delete TQ_Pt;
  delete SR_TQ_pt;
  delete CR_TQ_pt;
  delete SPS_TQ_pt;
  delete legend_TQ_pt;


//TQ eta
  THStack *TQ_eta = new THStack("TQ_eta", "Tetraquark eta; #eta_{TQ}; Events / 0.2");

  TH1F *SR_TQ_eta = new TH1F("SR_TQ_eta", "TQ eta SR", 100, -10, 10);
  treeTQ->Draw("TQ_eta>>SR_TQ_eta");  
  TH1F *CR_TQ_eta = new TH1F("CR_TQ_eta", "TQ eta CR", 100, -10, 10);
  treeCR->Draw("TQ_eta>>CR_TQ_eta");
  TH1F *SPS_TQ_eta = new TH1F("SPS_TQ_eta", "TQ eta SPS", 100, -10, 10);
  treeSPS->Draw("TQ_eta>>SPS_TQ_eta");           

  SR_TQ_eta->Sumw2();
  CR_TQ_eta->Sumw2();
  SPS_TQ_eta->Sumw2();
                                                                                             
  double int_TQ_eta = SR_TQ_eta->Integral();
  CR_TQ_eta->Scale(int_TQ_eta/CR_TQ_eta->Integral());  
  SPS_TQ_eta->Scale(weightSPS/SPS_TQ_eta->Integral()); 
                                                                                             
  SR_TQ_eta->SetLineColor(1);
  CR_TQ_eta->SetLineColor(2);
  SPS_TQ_eta->SetLineColor(4);
  SR_TQ_eta->SetMarkerColor(1);
  CR_TQ_eta->SetMarkerColor(2);
  SPS_TQ_eta->SetMarkerColor(4);
                                                                                             
  TQ_eta->Add(SR_TQ_eta);
  TQ_eta->Add(CR_TQ_eta);
  TQ_eta->Add(SPS_TQ_eta);
                                                                                             
  TCanvas *TQ_Eta = new TCanvas("TQ_Eta", "Tetraquark eta", 1500, 1000);
  TQ_eta->Draw("HIST nostack");
  TLegend* legend_TQ_eta = new TLegend(0.60, 0.78, 0.80, 0.88);  
  legend_TQ_eta->AddEntry(SR_TQ_eta, "Signal Region", "l");
  legend_TQ_eta->AddEntry(CR_TQ_eta, "Control Region", "l");
  legend_TQ_eta->AddEntry(SPS_TQ_eta, "SPS", "l");
  legend_TQ_eta->SetBorderSize(0);
  legend_TQ_eta->Draw();
  TQ_Eta->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/etaTQ.png").c_str());
  delete TQ_eta;
  delete TQ_Eta;
  delete SR_TQ_eta;
  delete CR_TQ_eta;
  delete SPS_TQ_eta;
  delete legend_TQ_eta;


//Ym mass
  THStack *Ym_mass = new THStack("Ym_mass", "Ym mass; m_{Ym} [GeV]; Events / 0.2 GeV");
  THStack *Ym_mass_small = new THStack("Ym_mass_small", "Ym mass; m_{Ym} [GeV]; Events / 0.05 GeV");

  TH1F *SR_Ym_mass = new TH1F("SR_Ym_mass", "Ym mass SR", 100, 0, 20);
  treeTQ->Draw("Ym_mass>>SR_Ym_mass");
  TH1F *CR_Ym_mass = new TH1F("CR_Ym_mass", "Ym mass CR", 100, 0, 20);
  treeCR->Draw("Ym_mass>>CR_Ym_mass");
  TH1F *SPS_Ym_mass = new TH1F("SPS_Ym_mass", "Ym mass CR", 100, 0, 20);
  treeSPS->Draw("Ym_mass>>SPS_Ym_mass");          
  TH1F *SR_Ym_mass_small = new TH1F("SR_Ym_mass_small", "Ym mass small", 100, 7, 12);
  treeTQ->Draw("Ym_mass>>SR_Ym_mass_small");
  TH1F *CR_Ym_mass_small = new TH1F("CR_Ym_mass_small", "Ym mass CR", 100, 7, 12);
  treeCR->Draw("Ym_mass>>CR_Ym_mass_small");
  TH1F *SPS_Ym_mass_small = new TH1F("SPS_Ym_mass_small", "Ym mass SPS", 100, 7, 12);
  treeSPS->Draw("Ym_mass>>SPS_Ym_mass_small");
 

  SR_Ym_mass->Sumw2();
  CR_Ym_mass->Sumw2();
  SPS_Ym_mass->Sumw2();
  SR_Ym_mass_small->Sumw2();
  CR_Ym_mass_small->Sumw2();
  SPS_Ym_mass_small->Sumw2();
                                                                                             
  double int_Ym_mass = SR_Ym_mass->Integral();
  CR_Ym_mass->Scale(int_Ym_mass/CR_Ym_mass->Integral());  
  SPS_Ym_mass->Scale(weightSPS/SPS_Ym_mass->Integral());
  double int_Ym_mass_small = SR_Ym_mass_small->Integral();
  CR_Ym_mass_small->Scale(int_Ym_mass_small/CR_Ym_mass_small->Integral());
  SPS_Ym_mass_small->Scale(weightSPS/SPS_Ym_mass_small->Integral());
                                                                                       
  SR_Ym_mass->SetLineColor(1);
  CR_Ym_mass->SetLineColor(2);
  SPS_Ym_mass->SetLineColor(4);
  SR_Ym_mass->SetMarkerColor(1);
  CR_Ym_mass->SetMarkerColor(2);
  SPS_Ym_mass->SetMarkerColor(4);
  SR_Ym_mass_small->SetLineColor(1);
  CR_Ym_mass_small->SetLineColor(2);
  SPS_Ym_mass_small->SetLineColor(4);
  SR_Ym_mass_small->SetMarkerColor(1);
  CR_Ym_mass_small->SetMarkerColor(2);
  SPS_Ym_mass_small->SetMarkerColor(4);
                                                                                             
  Ym_mass->Add(SR_Ym_mass);
  Ym_mass->Add(CR_Ym_mass);
  Ym_mass->Add(SPS_Ym_mass);
  Ym_mass_small->Add(SR_Ym_mass_small);
  Ym_mass_small->Add(CR_Ym_mass_small);
  Ym_mass_small->Add(SPS_Ym_mass_small);
                                                                                             
  TCanvas *Ym_Mass = new TCanvas("Ym_Mass", "Ym mass", 1500, 1500);
  Ym_Mass->Divide(1, 2);
  Ym_Mass->cd(1);
  Ym_mass->Draw("HIST nostack");
  TLegend* legend_Ym_mass = new TLegend(0.60, 0.78, 0.80, 0.88);   
  legend_Ym_mass->AddEntry(SR_Ym_mass, "Signal Region", "l");
  legend_Ym_mass->AddEntry(CR_Ym_mass, "Control Region", "l");
  legend_Ym_mass->AddEntry(SPS_Ym_mass, "SPS", "l");
  legend_Ym_mass->SetBorderSize(0);
  legend_Ym_mass->Draw();
  Ym_Mass->cd(2);
  Ym_mass_small->Draw("HIST nostack");    
  //Ym_mass_small->GetXaxis()->SetLimits(7., 12.);                                 
  TLegend* legend_Ym_mass_small = new TLegend(0.60, 0.78, 0.80, 0.88);   
  legend_Ym_mass_small->AddEntry(SR_Ym_mass_small, "Signal Region", "l");
  legend_Ym_mass_small->AddEntry(CR_Ym_mass_small, "Control Region", "l");
  legend_Ym_mass_small->AddEntry(SPS_Ym_mass_small, "SPS", "l");
  legend_Ym_mass_small->SetBorderSize(0);
  legend_Ym_mass_small->Draw();
  Ym_Mass->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/massYm.png").c_str());
  delete Ym_mass;
  delete Ym_mass_small;
  delete Ym_Mass;
  delete SR_Ym_mass;
  delete SR_Ym_mass_small;
  delete CR_Ym_mass;
  delete CR_Ym_mass_small;
  delete SPS_Ym_mass;
  delete SPS_Ym_mass_small;
  delete legend_Ym_mass;
  delete legend_Ym_mass_small;


//Ym pt
  THStack *Ym_pt = new THStack("Ym_pt", "Ym pt; p_{T}_{Ym} [GeV]; Events / 0.6 GeV");

  TH1F *SR_Ym_pt = new TH1F("SR_Ym_pt", "", 100, 0, 60);  
  treeTQ->Draw("Ym_pt>>SR_Ym_pt");
  TH1F *CR_Ym_pt = new TH1F("CR_Ym_pt", "", 100, 0, 60);
  treeCR->Draw("Ym_pt>>CR_Ym_pt");
  TH1F *SPS_Ym_pt = new TH1F("SPS_Ym_pt", "", 100, 0, 60);
  treeSPS->Draw("Ym_pt>>SPS_Ym_pt");

  SR_Ym_pt->Sumw2();
  CR_Ym_pt->Sumw2();
  SPS_Ym_pt->Sumw2();
                                                                                             
  double int_Ym_pt = SR_Ym_pt->Integral();
  CR_Ym_pt->Scale(int_Ym_pt/CR_Ym_pt->Integral());  
  SPS_Ym_pt->Scale(weightSPS/SPS_Ym_pt->Integral()); 
                                                                                             
  SR_Ym_pt->SetLineColor(1);
  CR_Ym_pt->SetLineColor(2);
  SPS_Ym_pt->SetLineColor(4);
  SR_Ym_pt->SetMarkerColor(1);
  CR_Ym_pt->SetMarkerColor(2);
  SPS_Ym_pt->SetMarkerColor(4);
                                                                                             
  Ym_pt->Add(SR_Ym_pt);
  Ym_pt->Add(CR_Ym_pt);
  Ym_pt->Add(SPS_Ym_pt);
                                                                                             
  TCanvas *Ym_Pt = new TCanvas("Ym_Pt", "Ym pt", 1500, 1000);
  Ym_pt->Draw("HIST nostack");
  TLegend* legend_Ym_pt = new TLegend(0.60, 0.78, 0.80, 0.88);     
  legend_Ym_pt->AddEntry(SR_Ym_pt, "Signal Region", "l");
  legend_Ym_pt->AddEntry(CR_Ym_pt, "Control Region", "l");
  legend_Ym_pt->AddEntry(SPS_Ym_pt, "SPS", "l");
  legend_Ym_pt->SetBorderSize(0);
  legend_Ym_pt->Draw();
  Ym_Pt->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/ptYm.png").c_str());
  delete Ym_pt;
  delete Ym_Pt;
  delete SR_Ym_pt;
  delete CR_Ym_pt;
  delete SPS_Ym_pt;
  delete legend_Ym_pt;

//Ym eta
  THStack *Ym_eta = new THStack("Ym_eta", "Ym eta; #eta_{Ym}; Events / 0.2");

  TH1F *SR_Ym_eta = new TH1F("SR_Ym_eta", "", 100, -10, 10);
  treeTQ->Draw("Ym_eta>>SR_Ym_eta");  
  TH1F *CR_Ym_eta = new TH1F("CR_Ym_eta", "", 100, -10, 10);
  treeCR->Draw("Ym_eta>>CR_Ym_eta");
  TH1F *SPS_Ym_eta = new TH1F("SPS_Ym_eta", "", 100, -10, 10);
  treeSPS->Draw("Ym_eta>>SPS_Ym_eta");           

  SR_Ym_eta->Sumw2();
  CR_Ym_eta->Sumw2();
  SPS_Ym_eta->Sumw2();
                                                                                             
  double int_Ym_eta = SR_Ym_eta->Integral();
  CR_Ym_eta->Scale(int_Ym_eta/CR_Ym_eta->Integral());  
  SPS_Ym_eta->Scale(weightSPS/SPS_Ym_eta->Integral()); 
                                                                                             
  SR_Ym_eta->SetLineColor(1);
  CR_Ym_eta->SetLineColor(2);
  SPS_Ym_eta->SetLineColor(4);
  SR_Ym_eta->SetMarkerColor(1);
  CR_Ym_eta->SetMarkerColor(2);
  SPS_Ym_eta->SetMarkerColor(4);
                                                                                             
  Ym_eta->Add(SR_Ym_eta);
  Ym_eta->Add(CR_Ym_eta);
  Ym_eta->Add(SPS_Ym_eta);
                                                                                             
  TCanvas *Ym_Eta = new TCanvas("Ym_Eta", "Ym eta", 1500, 1000);
  Ym_eta->Draw("HIST nostack");
  TLegend* legend_Ym_eta = new TLegend(0.60, 0.78, 0.80, 0.88);       
  legend_Ym_eta->AddEntry(SR_Ym_eta, "Signal Region", "l");
  legend_Ym_eta->AddEntry(CR_Ym_eta, "Control Region", "l");
  legend_Ym_eta->AddEntry(SPS_Ym_eta, "SPS", "l");
  legend_Ym_eta->SetBorderSize(0);
  legend_Ym_eta->Draw();
  Ym_Eta->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/etaYm.png").c_str());
  delete Ym_eta;
  delete Ym_Eta;
  delete SR_Ym_eta; 
  delete CR_Ym_eta;
  delete SPS_Ym_eta;
  delete legend_Ym_eta;

//Ye mass
  THStack *Ye_mass = new THStack("Ye_mass", "Ye mass; m_{Ye} [GeV]; Events / 0.3 GeV");

  TH1F *SR_Ye_mass = new TH1F("SR_Ye_mass", "", 100, 0, 30);
  treeTQ->Draw("Ye_mass>>SR_Ye_mass");  
  TH1F *CR_Ye_mass = new TH1F("CR_Ye_mass", "", 100, 0, 30);
  treeCR->Draw("Ye_mass>>CR_Ye_mass");
  TH1F *SPS_Ye_mass = new TH1F("SPS_Ye_mass", "", 100, 0, 30);
  treeSPS->Draw("Ye_mass>>SPS_Ye_mass");           

  SR_Ye_mass->Sumw2();
  CR_Ye_mass->Sumw2();
  SPS_Ye_mass->Sumw2();
                                                                                             
  double int_Ye_mass = SR_Ye_mass->Integral();
  CR_Ye_mass->Scale(int_Ye_mass/CR_Ye_mass->Integral());  
  SPS_Ye_mass->Scale(weightSPS/SPS_Ye_mass->Integral()); 
                                                                                             
  SR_Ye_mass->SetLineColor(1);
  CR_Ye_mass->SetLineColor(2);
  SPS_Ye_mass->SetLineColor(4);
  SR_Ye_mass->SetMarkerColor(1);
  CR_Ye_mass->SetMarkerColor(2);
  SPS_Ye_mass->SetMarkerColor(4);
                                                                                             
  Ye_mass->Add(SR_Ye_mass);
  Ye_mass->Add(CR_Ye_mass);
  Ye_mass->Add(SPS_Ye_mass);
                                                                                             
  TCanvas *Ye_Mass = new TCanvas("Ye_Mass", "Ye mass", 1500, 1000);
  Ye_mass->Draw("HIST nostack");
  TLegend* legend_Ye_mass = new TLegend(0.60, 0.78, 0.80, 0.88);   
  legend_Ye_mass->AddEntry(SR_Ye_mass, "Signal Region", "l");
  legend_Ye_mass->AddEntry(CR_Ye_mass, "Control Region", "l");
  legend_Ye_mass->AddEntry(SPS_Ye_mass, "SPS", "l");
  legend_Ye_mass->SetBorderSize(0);
  legend_Ye_mass->Draw();
  Ye_Mass->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/massYe.png").c_str());
  delete Ye_mass;
  delete Ye_Mass;
  delete SR_Ye_mass;
  delete CR_Ye_mass;
  delete SPS_Ye_mass;
  delete legend_Ye_mass;


//Ye pt
  THStack *Ye_pt = new THStack("Ye_pt", "Ye pt; p_{T}_{Ye} [GeV]; Events / 0.6 GeV");

  TH1F *SR_Ye_pt= new TH1F("SR_Ye_pt", "", 100, 0, 60);  
  treeTQ->Draw("Ye_pt>>SR_Ye_pt");
  TH1F *CR_Ye_pt = new TH1F("CR_Ye_pt", "", 100, 0, 60);
  treeCR->Draw("Ye_pt>>CR_Ye_pt");
  TH1F *SPS_Ye_pt = new TH1F("SPS_Ye_pt", "", 100, 0, 60);
  treeSPS->Draw("Ye_pt>>SPS_Ye_pt");
  

  SR_Ye_pt->Sumw2();
  CR_Ye_pt->Sumw2();
  SPS_Ye_pt->Sumw2();
                                                                                             
  double int_Ye_pt = SR_Ye_pt->Integral();
  CR_Ye_pt->Scale(int_Ye_pt/CR_Ye_pt->Integral());  
  SPS_Ye_pt->Scale(weightSPS/SPS_Ye_pt->Integral()); 
                                                                                             
  SR_Ye_pt->SetLineColor(1);
  CR_Ye_pt->SetLineColor(2);
  SPS_Ye_pt->SetLineColor(4);
  SR_Ye_pt->SetMarkerColor(1);
  CR_Ye_pt->SetMarkerColor(2);
  SPS_Ye_pt->SetMarkerColor(4);
                                                                                             
  Ye_pt->Add(SR_Ye_pt);
  Ye_pt->Add(CR_Ye_pt);
  Ye_pt->Add(SPS_Ye_pt);
                                                                                             
  TCanvas *Ye_Pt = new TCanvas("Ye_Pt", "Ye pt", 1500, 1000);
  Ye_pt->Draw("HIST nostack");
  TLegend* legend_Ye_pt = new TLegend(0.60, 0.78, 0.80, 0.88);     
  legend_Ye_pt->AddEntry(SR_Ye_pt, "Signal Region", "l");
  legend_Ye_pt->AddEntry(CR_Ye_pt, "Control Region", "l");
  legend_Ye_pt->AddEntry(SPS_Ye_pt, "SPS", "l");
  legend_Ye_pt->SetBorderSize(0);
  legend_Ye_pt->Draw();
  Ye_Pt->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/ptYe.png").c_str());
  delete Ye_pt;
  delete Ye_Pt;
  delete SR_Ye_pt;
  delete CR_Ye_pt;
  delete SPS_Ye_pt;
  delete legend_Ye_pt;

//Ye eta
  THStack *Ye_eta = new THStack("Ye_eta", "Ye eta; #eta_{Ye}; Events / 0.2");

  TH1F *SR_Ye_eta = new TH1F("SR_Ye_eta", "", 100, -10, 10);
  treeTQ->Draw("Ye_eta>>SR_Ye_eta");  
  TH1F *CR_Ye_eta = new TH1F("CR_Ye_eta", "", 100, -10, 10);
  treeCR->Draw("Ye_eta>>CR_Ye_eta");
  TH1F *SPS_Ye_eta = new TH1F("SPS_Ye_eta", "", 100, -10, 10);
  treeSPS->Draw("Ye_eta>>SPS_Ye_eta");           

  SR_Ye_eta->Sumw2();
  CR_Ye_eta->Sumw2();
  SPS_Ye_eta->Sumw2();
                                                                                             
  double int_Ye_eta = SR_Ye_eta->Integral();
  CR_Ye_eta->Scale(int_Ye_eta/CR_Ye_eta->Integral());  
  SPS_Ye_eta->Scale(weightSPS/SPS_Ye_eta->Integral()); 
                                                                                             
  SR_Ye_eta->SetLineColor(1);
  CR_Ye_eta->SetLineColor(2);
  SPS_Ye_eta->SetLineColor(4);
  SR_Ye_eta->SetMarkerColor(1);
  CR_Ye_eta->SetMarkerColor(2);
  SPS_Ye_eta->SetMarkerColor(4);
                                                                                             
  Ye_eta->Add(SR_Ye_eta);
  Ye_eta->Add(CR_Ye_eta);
  Ye_eta->Add(SPS_Ye_eta);
                                                                                             
  TCanvas *Ye_Eta = new TCanvas("Ye_Eta", "Ye eta", 1500, 1000);
  Ye_eta->Draw("HIST nostack");
  TLegend* legend_Ye_eta = new TLegend(0.60, 0.78, 0.80, 0.88);       
  legend_Ye_eta->AddEntry(SR_Ye_eta, "Signal Region", "l");
  legend_Ye_eta->AddEntry(CR_Ye_eta, "Control Region", "l");
  legend_Ye_eta->AddEntry(SPS_Ye_eta, "SPS", "l");
  legend_Ye_eta->SetBorderSize(0);
  legend_Ye_eta->Draw();
  Ye_Eta->SaveAs(("/eos/home-l/lfantini/www/Background/"+CR_input+"/etaYe.png").c_str());
  delete Ye_eta;
  delete Ye_Eta;
  delete SR_Ye_eta;
  delete CR_Ye_eta;
  delete SPS_Ye_eta;
  delete legend_Ye_eta;

}














 


  



