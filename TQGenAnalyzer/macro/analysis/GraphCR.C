#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TAxis.h"
#include "TTree.h"
#include <iostream>
using namespace std;


void Graph(std::string CR_input){
  ///////////Open the tree////////////////
  TFile *fTQ = new TFile("fout_mRun2018_B_C_trigger_TQ.root", "READ");
  TTree *treeTQ = (TTree*)fTQ->Get("tree_red");
  TFile *fCR = new TFile(("fout_mRun2018_B_C_trigger_"+CR_input+".root").c_str(), "READ");
  TTree *treeCR = (TTree*)fCR->Get("tree_red");
  
  ///////////Histograms////////////
  //SR
  TH1F *SR = new TH1F("SR", "TQ invariant mass CR", 100, 0, 60);
  SR->GetXaxis()->SetTitle("m_{TQ} [GeV]");
  SR->GetYaxis()->SetTitle("Events");

  TH1F *SR_pt = new TH1F("SR_pt", "TQ p_{T} CR", 100, 0, 60);
  SR_pt->GetXaxis()->SetTitle("p_{T} [GeV]");
  SR_pt->GetYaxis()->SetTitle("Events");

  TH1F *SR_Ye_mass = new TH1F("SR_Ye_mass", "Ye invariant mass CR", 100, 0, 30);
  SR_Ye_mass->GetXaxis()->SetTitle("m_{Y} [GeV]");
  SR_Ye_mass->GetYaxis()->SetTitle("Events");

  TH1F *SR_mass_tilde = new TH1F("SR_mass_tilde", "TQ tilde mass CR", 100, 0, 60);
  SR_mass_tilde->GetXaxis()->SetTitle("#tilde{m}_{TQ} [GeV]");
  SR_mass_tilde->GetYaxis()->SetTitle("Events");

  TH1F *SR_eta = new TH1F("SR_eta", "TQ eta CR", 100, -10, 10);
  SR_eta->GetXaxis()->SetTitle("#eta");
  SR_eta->GetYaxis()->SetTitle("Events");

  TH1F *SR_Ye_pt = new TH1F("SR_Ye_pt", "Ye pt CR", 100, 0, 60);
  SR_Ye_pt->GetXaxis()->SetTitle("p_{T} [GeV]");
  SR_Ye_pt->GetYaxis()->SetTitle("Events");

  TH1F *SR_Ye_eta = new TH1F("SR_Ye_eta", "Ye eta CR", 100, -10, 10);
  SR_Ye_eta->GetXaxis()->SetTitle("#eta");
  SR_Ye_eta->GetYaxis()->SetTitle("Events");

  TH1F *SR_Ym_pt = new TH1F("SR_Ym_pt", "Ym pt CR", 100, 0, 60);     
  SR_Ym_pt->GetXaxis()->SetTitle("p_{T} [GeV]");
  SR_Ym_pt->GetYaxis()->SetTitle("Events");
                                                                       
  TH1F *SR_Ym_eta = new TH1F("SR_Ym_eta", "Ym eta CR", 100, -10, 10);
  SR_Ym_eta->GetXaxis()->SetTitle("#eta");
  SR_Ym_eta->GetYaxis()->SetTitle("Events");


  TH1F *SR_Ym_mass = new TH1F("SR_Ym_mass", "Ym invariant mass CR", 100, 0, 25);
  SR_Ym_mass->GetXaxis()->SetTitle("m_{Y} [GeV]");
  SR_Ym_mass->GetYaxis()->SetTitle("Events");

  TH2F *SR_Y_mass_2d = new TH2F("SR_Y_mass_2d", "Y invariant mass SR", 100, -1, 25, 100, -1, 30);
  SR_Y_mass_2d->GetXaxis()->SetTitle("m_{#mu #mu} [GeV]");
  SR_Y_mass_2d->GetYaxis()->SetTitle("m_{ee} [GeV]");




  //CR
  TH1F *CR = new TH1F("CR", "TQ invariant mass CR", 100, 0, 60);
  CR->GetXaxis()->SetTitle("m_{TQ} [GeV]");
  CR->GetYaxis()->SetTitle("Events"); 
 
  TH1F *CR_pt = new TH1F("CR_pt", "TQ p_{T} CR", 100, 0, 60);
  CR_pt->GetXaxis()->SetTitle("p_{T} [GeV]");
  CR_pt->GetYaxis()->SetTitle("Events");
                                                                                  
  TH1F *CR_Ye_mass = new TH1F("CR_Ye_mass", "Ye invariant mass CR", 100, 0, 30);
  CR_Ye_mass->GetXaxis()->SetTitle("m_{Y} [GeV]");
  CR_Ye_mass->GetYaxis()->SetTitle("Events");
                                                                                  
  TH1F *CR_mass_tilde = new TH1F("CR_mass_tilde", "TQ tilde mass CR", 100, 0, 60);
  CR_mass_tilde->GetXaxis()->SetTitle("#tilde{m}_{TQ} [GeV]");
  CR_mass_tilde->GetYaxis()->SetTitle("Events");

  TH1F *CR_eta = new TH1F("CR_eta", "TQ eta CR", 100, -10, 10);
  CR_eta->GetXaxis()->SetTitle("#eta");
  CR_eta->GetYaxis()->SetTitle("Events");
 
  TH1F *CR_Ye_pt = new TH1F("CR_Ye_pt", "Ye pt CR", 100, 0, 60);
  CR_Ye_pt->GetXaxis()->SetTitle("p_{T} [GeV]");
  CR_Ye_pt->GetYaxis()->SetTitle("Events");
 
  TH1F *CR_Ye_eta = new TH1F("CR_Ye_eta", "Ye eta CR", 100, -10, 10);
  CR_Ye_eta->GetXaxis()->SetTitle("#eta");
  CR_Ye_eta->GetYaxis()->SetTitle("Events");
 
  TH1F *CR_Ym_pt = new TH1F("CR_Ym_pt", "Ym pt CR", 100, 0, 60);     
  CR_Ym_pt->GetXaxis()->SetTitle("p_{T} [GeV]");
  CR_Ym_pt->GetYaxis()->SetTitle("Events");
                                                                        
  TH1F *CR_Ym_eta = new TH1F("CR_Ym_eta", "Ym eta CR", 100, -10, 10);
  CR_Ym_eta->GetXaxis()->SetTitle("#eta");
  CR_Ym_eta->GetYaxis()->SetTitle("Events");
 
 
  TH1F *CR_Ym_mass = new TH1F("CR_Ym_mass", "Ym invariant mass CR", 100, 0, 25);
  CR_Ym_mass->GetXaxis()->SetTitle("m_{Y} [GeV]");
  CR_Ym_mass->GetYaxis()->SetTitle("Events");
 
  TH2F *CR_Y_mass_2d = new TH2F("CR_Y_mass_2d", "Y invariant mass CR", 100, -1, 25, 100, -1, 30);
  CR_Y_mass_2d->GetXaxis()->SetTitle("m_{#mu #mu} [GeV]");
  CR_Y_mass_2d->GetYaxis()->SetTitle("m_{ee} [GeV]");




                           
  ///////////Set branches/////////////
  //SR
  Float_t mass_TQ;
  Float_t pt_TQ;
  Float_t mass_Ye_TQ;
  Float_t mass_tilde_TQ;
  Float_t eta_TQ;
  Float_t pt_Ye_TQ;
  Float_t eta_Ye_TQ;
  Float_t pt_Ym_TQ;
  Float_t eta_Ym_TQ;
  Float_t mass_Ym_TQ;
  treeTQ->SetBranchAddress("TQ_mass", &mass_TQ);
  treeTQ->SetBranchAddress("TQ_pt", &pt_TQ);
  treeTQ->SetBranchAddress("Ye_mass", &mass_Ye_TQ);
  treeTQ->SetBranchAddress("TQ_mass_tilde", &mass_tilde_TQ);
  treeTQ->SetBranchAddress("TQ_eta", &eta_TQ);
  treeTQ->SetBranchAddress("Ye_pt", &pt_Ye_TQ);
  treeTQ->SetBranchAddress("Ye_eta", &eta_Ye_TQ);
  treeTQ->SetBranchAddress("Ym_pt", &pt_Ym_TQ);
  treeTQ->SetBranchAddress("Ym_eta", &eta_Ym_TQ);
  treeTQ->SetBranchAddress("Ym_mass", &mass_Ym_TQ);

  //CR
  Float_t mass_CR;
  Float_t pt_CR;
  Float_t mass_Ye_CR;
  Float_t mass_tilde_CR;
  Float_t eta_CR;
  Float_t pt_Ye_CR;
  Float_t eta_Ye_CR;
  Float_t pt_Ym_CR;
  Float_t eta_Ym_CR;
  Float_t mass_Ym_CR;
  treeCR->SetBranchAddress("TQ_mass", &mass_CR);
  treeCR->SetBranchAddress("TQ_pt", &pt_CR);
  treeCR->SetBranchAddress("Ye_mass", &mass_Ye_CR);
  treeCR->SetBranchAddress("TQ_mass_tilde", &mass_tilde_CR); 
  treeCR->SetBranchAddress("TQ_eta", &eta_CR);
  treeCR->SetBranchAddress("Ye_pt", &pt_Ye_CR);
  treeCR->SetBranchAddress("Ye_eta", &eta_Ye_CR);
  treeCR->SetBranchAddress("Ym_pt", &pt_Ym_CR);
  treeCR->SetBranchAddress("Ym_eta", &eta_Ym_CR);
  treeCR->SetBranchAddress("Ym_mass", &mass_Ym_CR);



  ///////loops////////

  Long64_t nentriesTQ = (Long64_t)treeTQ->GetEntries();
  Long64_t nentriesCR = (Long64_t)treeCR->GetEntries();

  for (int i = 0; i<nentriesTQ; i++){
    treeTQ->GetEntry(i);
   
    SR->Fill(mass_TQ);
    SR_pt->Fill(pt_TQ);
    SR_Ye_mass->Fill(mass_Ye_TQ);
    SR_eta->Fill(eta_TQ);
    SR_Ye_pt->Fill(pt_Ye_TQ);
    SR_Ye_eta->Fill(eta_Ye_TQ);
    SR_Ym_pt->Fill(pt_Ym_TQ);
    SR_Ym_eta->Fill(eta_Ym_TQ);
    SR_Ym_mass->Fill(mass_Ym_TQ);
    SR_Y_mass_2d->Fill(mass_Ym_TQ, mass_Ye_TQ);
    SR_mass_tilde->Fill(mass_tilde_TQ); 
  }

  for (int i = 0; i<nentriesCR; i++){
    treeCR->GetEntry(i); 

    CR->Fill(mass_CR);
    CR_pt->Fill(pt_CR);
    CR_Ye_mass->Fill(mass_Ye_CR);
    CR_mass_tilde->Fill(mass_tilde_CR);
    CR_eta->Fill(eta_CR);
    CR_Ye_pt->Fill(pt_Ye_CR);
    CR_Ye_eta->Fill(eta_Ye_CR);
    CR_Ym_pt->Fill(pt_Ym_CR);
    CR_Ym_eta->Fill(eta_Ym_CR);
    CR_Ym_mass->Fill(mass_Ym_CR);
    CR_Y_mass_2d->Fill(mass_Ym_CR, mass_Ye_CR);
  }  


  //////Canvas/////////

  TCanvas *Mass = new TCanvas("Mass", "CR_SR_mass", 1500, 1000);
  double int_SR = SR->Integral();
  CR->Scale(int_SR/CR->Integral());
  //cout<<SR->Integral()<<endl;
  //cout<<CR->Integral()<<endl;
  CR->SetLineColor(2);
  int nbinsx = SR->GetXaxis()->GetNbins();
  for(int j=0; j<nbinsx; j++){
    double up = SR->GetXaxis()->GetBinUpEdge(j);
    if(up>16. && up<20.){
     SR-> SetBinContent(j, 0);
    }
  }
  SR->Draw("HIST");
  CR->Draw("HIST SAME");
  TLegend* legend_mass = new TLegend(0.15, 0.78, 0.35, 0.88);
  legend_mass->AddEntry(SR, "Signal Region");
  legend_mass->AddEntry(CR, "Control Region");
  legend_mass->SetTextSize(0.03);
  legend_mass->Draw();
  Mass->SaveAs(("/eos/user/l/lfantini/www/"+CR_input+"/mass.png").c_str());
  delete Mass;
  delete legend_mass;

  for(int j=0; j<nbinsx; j++){
    double up = CR->GetXaxis()->GetBinUpEdge(j);
    if(up>16. && up<20.){
      CR-> SetBinContent(j, 0);
     }
  }

  int_SR = SR->Integral();   
  CR->Scale(int_SR/CR->Integral());
  Double_t test = SR->KolmogorovTest(CR);

  cout<<"Kolmogorov test: "<<test<<endl;

  //cout<<SR->Integral()<<endl;
  //cout<<CR->Integral()<<endl;

  delete SR;
  delete CR;


  TCanvas *Pt = new TCanvas("Pt", "CR_SR_pt", 1500, 1000);
  double int_SR_pt = SR_pt->Integral();
  CR_pt->Scale(int_SR_pt/CR_pt->Integral());
  CR_pt->SetLineColor(2);
  SR_pt->Draw("HIST");
  CR_pt->Draw("HIST SAME");
  TLegend* legend_pt = new TLegend(0.45, 0.78, 0.65, 0.88);
  legend_pt->AddEntry(SR_pt, "Signal Region");
  legend_pt->AddEntry(CR_pt, "Control Region");
  legend_pt->SetTextSize(0.03);
  legend_pt->Draw();

  //cout<<SR_pt->Integral()<<endl;
  //cout<<CR_pt->Integral()<<endl;
  Pt->SaveAs(("/eos/user/l/lfantini/www/"+CR_input+"/pt.png").c_str());
  delete Pt;
  delete SR_pt;
  delete CR_pt;
  delete legend_pt;


  TCanvas *Ye_Mass = new TCanvas("Ye_Mass", "CR_SR_Ye_mass", 1500, 1000);
  SR_Ye_mass->Draw("HIST");
  double int_SR_Ye_mass = SR_Ye_mass->Integral();
  CR_Ye_mass->Scale(int_SR_Ye_mass/CR_Ye_mass->Integral());
  CR_Ye_mass->SetLineColor(2);
  CR_Ye_mass->Draw("HIST SAME");
  TLegend* legend_Ye_mass = new TLegend(0.15, 0.78, 0.35, 0.88);
  legend_Ye_mass->AddEntry(SR_Ye_mass, "Signal Region");
  legend_Ye_mass->AddEntry(CR_Ye_mass, "Control Region");
  legend_Ye_mass->SetTextSize(0.03);
  legend_Ye_mass->Draw();
  Ye_Mass->SaveAs(("/eos/user/l/lfantini/www/"+CR_input+"/Ye_mass.png").c_str());
  delete Ye_Mass;
  delete SR_Ye_mass;
  delete CR_Ye_mass;
  delete legend_Ye_mass;


  TCanvas *Mass_Tilde = new TCanvas("Mass_Tilde", "CR_SR_mass_tilde", 1500, 1000);
  double int_SR_mass_tilde = SR_mass_tilde->Integral();
  CR_mass_tilde->Scale(int_SR_mass_tilde/CR_mass_tilde->Integral());
  CR_mass_tilde->SetLineColor(2);

  int nbinsx_tilde = SR_mass_tilde->GetXaxis()->GetNbins();
  for(int j=0; j<nbinsx_tilde; j++){
    double up_tilde = SR_mass_tilde->GetXaxis()->GetBinUpEdge(j);
    if(up_tilde>16. && up_tilde<20.){
     SR_mass_tilde-> SetBinContent(j, 0);
    }
  }

  SR_mass_tilde->Draw("HIST");
  CR_mass_tilde->Draw("HIST SAME");
  TLegend* legend_mass_tilde = new TLegend(0.15, 0.78, 0.35, 0.88);
  legend_mass_tilde->AddEntry(SR_mass_tilde, "Signal Region");
  legend_mass_tilde->AddEntry(CR_mass_tilde, "Control Region");
  legend_mass_tilde->SetTextSize(0.03);
  legend_mass_tilde->Draw();
  Mass_Tilde->SaveAs(("/eos/user/l/lfantini/www/"+CR_input+"/mass_tilde.png").c_str());
  delete Mass_Tilde;
  delete SR_mass_tilde;
  delete CR_mass_tilde;
  delete legend_mass_tilde;


  TCanvas *Eta = new TCanvas("Eta", "CR_SR_eta", 1500, 1000);
  double int_SR_eta = SR_eta->Integral();
  CR_eta->Scale(int_SR_eta/CR_eta->Integral());
  CR_eta->SetLineColor(2);
  CR_eta->Draw("HIST");
  SR_eta->Draw("HIST SAME");
  TLegend* legend_eta = new TLegend(0.15, 0.78, 0.35, 0.88);
  legend_eta->AddEntry(SR_eta, "Signal Region");
  legend_eta->AddEntry(CR_eta, "Control Region");
  legend_eta->SetTextSize(0.03);
  legend_eta->Draw();
  Eta->SaveAs(("/eos/user/l/lfantini/www/"+CR_input+"/eta.png").c_str());
  delete Eta;
  delete SR_eta;
  delete CR_eta;
  delete legend_eta;


  TCanvas *Ye_Pt = new TCanvas("Ye_Pt", "CR_SR_Ye_pt", 1500, 1000);                         
  SR_Ye_pt->Draw("HIST");
  double int_SR_Ye_pt = SR_Ye_pt->Integral();
  CR_Ye_pt->Scale(int_SR_Ye_pt/CR_Ye_pt->Integral());
  CR_Ye_pt->SetLineColor(2);
  CR_Ye_pt->Draw("HIST SAME");
  TLegend* legend_Ye_pt = new TLegend(0.35, 0.78, 0.55, 0.88);
  legend_Ye_pt->AddEntry(SR_Ye_pt, "Signal Region");
  legend_Ye_pt->AddEntry(CR_Ye_pt, "Control Region");
  legend_Ye_pt->SetTextSize(0.03);
  legend_Ye_pt->Draw();
  Ye_Pt->SaveAs(("/eos/user/l/lfantini/www/"+CR_input+"/Ye_pt.png").c_str());
  delete Ye_Pt;
  delete SR_Ye_pt;
  delete CR_Ye_pt;
  delete legend_Ye_pt;


  TCanvas *Ye_Eta = new TCanvas("Ye_Eta", "CR_SR_Ye_eta", 1500, 1000);                         
  double int_SR_Ye_eta = SR_Ye_eta->Integral();
  CR_Ye_eta->Scale(int_SR_Ye_eta/CR_Ye_eta->Integral());
  CR_Ye_eta->SetLineColor(2);
  CR_Ye_eta->Draw("HIST");
  SR_Ye_eta->Draw("HIST SAME");
  TLegend* legend_Ye_eta = new TLegend(0.15, 0.78, 0.35, 0.88);
  legend_Ye_eta->AddEntry(SR_Ye_eta, "Signal Region");
  legend_Ye_eta->AddEntry(CR_Ye_eta, "Control Region");
  legend_Ye_eta->SetTextSize(0.03);
  legend_Ye_eta->Draw();
  Ye_Eta->SaveAs(("/eos/user/l/lfantini/www/"+CR_input+"/Ye_eta.png").c_str());
  delete Ye_Eta;
  delete SR_Ye_eta;
  delete CR_Ye_eta;
  delete legend_Ye_eta;


  TCanvas *Ym_Eta = new TCanvas("Ym_Eta", "CR_SR_Ym_eta", 1500, 1000);                         
  double int_SR_Ym_eta = SR_Ym_eta->Integral();
  CR_Ym_eta->Scale(int_SR_Ym_eta/CR_Ym_eta->Integral());
  CR_Ym_eta->SetLineColor(2);
  CR_Ym_eta->Draw("HIST");
  SR_Ym_eta->Draw("HIST SAME");  
  TLegend* legend_Ym_eta = new TLegend(0.15, 0.78, 0.35, 0.88);
  legend_Ym_eta->AddEntry(SR_Ym_eta, "Signal Region");
  legend_Ym_eta->AddEntry(CR_Ym_eta, "Control Region");
  legend_Ym_eta->SetTextSize(0.03);
  legend_Ym_eta->Draw();
  Ym_Eta->SaveAs(("/eos/user/l/lfantini/www/"+CR_input+"/Ym_eta.png").c_str());
  delete Ym_Eta;
  delete SR_Ym_eta;
  delete CR_Ym_eta;
  delete legend_Ym_eta;


  TCanvas *Ym_Pt = new TCanvas("Ym_Pt", "CR_SR_Ym_pt", 1500, 1000);                        
  SR_Ym_pt->Draw("HIST");
  double int_SR_Ym_pt = SR_Ym_pt->Integral();
  CR_Ym_pt->Scale(int_SR_Ym_pt/CR_Ym_pt->Integral());
  CR_Ym_pt->SetLineColor(2);
  CR_Ym_pt->Draw("HIST SAME");
  TLegend* legend_Ym_pt = new TLegend(0.45, 0.78, 0.65, 0.88);
  legend_Ym_pt->AddEntry(SR_Ym_pt, "Signal Region");
  legend_Ym_pt->AddEntry(CR_Ym_pt, "Control Region");
  legend_Ym_pt->SetTextSize(0.03);
  legend_Ym_pt->Draw();
  Ym_Pt->SaveAs(("/eos/user/l/lfantini/www/"+CR_input+"/Ym_pt.png").c_str());
  delete Ym_Pt;
  delete SR_Ym_pt;
  delete CR_Ym_pt;
  delete legend_Ym_pt;


  TCanvas *Ym_Mass = new TCanvas("Ym_Mass", "CR_SR_Ym_mass", 1500, 1000);
  SR_Ym_mass->Draw("HIST");
  double int_SR_Ym_mass = SR_Ym_mass->Integral();
  CR_Ym_mass->Scale(int_SR_Ym_mass/CR_Ym_mass->Integral());
  CR_Ym_mass->SetLineColor(2);
  CR_Ym_mass->Draw("HIST SAME");
  TLegend* legend_Ym_mass = new TLegend(0.15, 0.78, 0.35, 0.88);
  legend_Ym_mass->AddEntry(SR_Ym_mass, "Signal Region");
  legend_Ym_mass->AddEntry(CR_Ym_mass, "Control Region");
  legend_Ym_mass->SetTextSize(0.03);
  legend_Ym_mass->Draw();
  Ym_Mass->SaveAs(("/eos/user/l/lfantini/www/"+CR_input+"/Ym_mass.png").c_str());
  delete Ym_Mass;
  delete SR_Ym_mass;
  delete CR_Ym_mass;
  delete legend_Ym_mass;


  TCanvas *SR_Y_Mass_2D = new TCanvas("SR_Y_Mass_2D", "SR_Y_mass_2d", 1500, 1000);
  SR_Y_mass_2d->SetStats(0);
  SR_Y_mass_2d->Draw("COLZ");
  SR_Y_Mass_2D->SaveAs(("/eos/user/l/lfantini/www/"+CR_input+"/Y_mass_2d_SR.png").c_str());
  delete SR_Y_Mass_2D;
  delete SR_Y_mass_2d;

  TCanvas *CR_Y_Mass_2D = new TCanvas("CR_Y_Mass_2D", "CR_Y_mass_2d", 1500, 1000);
  CR_Y_mass_2d->SetStats(0);
  CR_Y_mass_2d->Draw("COLZ");                                                      
  CR_Y_Mass_2D->SaveAs(("/eos/user/l/lfantini/www/"+CR_input+"/Y_mass_2d_CR.png").c_str());                 
  delete CR_Y_Mass_2D;
  delete CR_Y_mass_2d;
}
