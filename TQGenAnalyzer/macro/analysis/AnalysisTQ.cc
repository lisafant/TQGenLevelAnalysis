
#define AnalysisTQ_cxx
#include "AnalysisTQ.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
using namespace std;
#include <stdio.h>
#include <iostream>
void AnalysisTQ::Loop(std::string mass)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   
   int counter[9]={0};
   int counterBIS=0;
   int counterTRIS=0;

   double massnom;
   if(mass=="7" || mass=="9")massnom=3.09;
   else massnom=9.46;

   TTree tree_red("tree_red","");

   Float_t e_pt1; 
   Float_t e_pt2;
   Float_t m_pt1;
   Float_t m_pt2;
   Float_t e_eta1;
   Float_t e_eta2;
   Float_t e_phi1;
   Float_t e_phi2;
   Float_t m_eta1;
   Float_t m_eta2;
   Float_t m_phi1;
   Float_t m_phi2;
   Int_t e_isPF1;
   Int_t e_isPF2;
   Float_t TQ_vtxprob;
   Float_t Ym_vtxprob;
   Float_t Ye_vtxprob;
   Float_t TQ_mass; 
   Float_t Ym_mass; 
   Float_t Ye_mass; 
   Float_t TQ_pt; 
   Float_t Ym_pt; 
   Float_t Ye_pt; 
   Float_t TQ_eta; 
   Float_t Ym_eta; 
   Float_t Ye_eta; 


   
   tree_red.Branch("e_isPF1",&e_isPF1,"e_isPF1/I");
   tree_red.Branch("e_isPF2",&e_isPF2,"e_isPF2/I");
   tree_red.Branch("e_pt1",&e_pt1,"e_pt1/F");
   tree_red.Branch("e_pt2",&e_pt2,"e_pt2/F");
   tree_red.Branch("m_pt1",&m_pt1,"m_pt1/F");
   tree_red.Branch("m_pt2",&m_pt2,"m_pt2/F");
   tree_red.Branch("e_eta1",&e_eta1,"e_eta1/F");
   tree_red.Branch("e_eta2",&e_eta2,"e_eta2/F");
   tree_red.Branch("e_phi1",&e_phi1,"e_phi1/F");
   tree_red.Branch("e_phi2",&e_phi2,"e_phi2/F");
   tree_red.Branch("m_eta1",&m_eta1,"m_eta1/F");
   tree_red.Branch("m_eta2",&m_eta2,"m_eta2/F");
   tree_red.Branch("m_phi1",&m_phi1,"m_phi1/F");
   tree_red.Branch("m_phi2",&m_phi2,"m_phi2/F");
   tree_red.Branch("TQ_vtxprob",&TQ_vtxprob,"TQ_vtxprob/F");
   tree_red.Branch("Ym_vtxprob",&Ym_vtxprob,"Ym_vtxprob/F");
   tree_red.Branch("Ye_vtxprob",&Ye_vtxprob,"Ye_vtxprob/F");
   tree_red.Branch("TQ_mass",&TQ_mass,"TQ_mass/F");
   tree_red.Branch("Ym_mass",&Ym_mass,"Ym_mass/F");
   tree_red.Branch("Ye_mass",&Ye_mass,"Ye_mass/F");
   tree_red.Branch("TQ_pt",&TQ_pt,"TQ_pt/F");
   tree_red.Branch("Ym_pt",&Ym_pt,"Ym_pt/F");
   tree_red.Branch("Ye_pt",&Ye_pt,"Ye_pt/F");
   tree_red.Branch("TQ_eta",&TQ_eta,"TQ_eta/F");
   tree_red.Branch("Ym_eta",&Ym_eta,"Ym_eta/F");
   tree_red.Branch("Ye_eta",&Ye_eta,"Ye_eta/F");



   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      int nTQwithAtLeast2mu=0;
      int nTQwithAtLeast2muPtEtaID=0;
      int nTQwithAtLeast2muPtEtaIDVtxProb=0;
      int nTQwithAtLeast2el=0;
      int nTQwithAtLeast2elOver=0;
      int nTQwithAtLeast2elOverPtEtaID=0;
      int nTQwithAtLeast2elOverPtEtaIDVtxProb=0;
      int nTQwithVtxProb=0;
      //      std::cout<< "---------------------------------------------" <<std::endl;

      counter[0]++;

      int best_cand=999;
      float prob_cand=0;

      // ------------- TQ -------------- //
      for(int i=0; i< nTQReco; i++){

	//muons
	if(((*recoTQ_charge1)[i]*(*recoTQ_charge2)[i])>0)continue;

	nTQwithAtLeast2mu++;

	if((*recoTQ_pt1)[i]<1 || (*recoTQ_pt2)[i]<1 || fabs((*recoTQ_eta1)[i])>2.5 || fabs((*recoTQ_eta2)[i])>2.5 || (*recoTQ_softID1)[i]==0 || (*recoTQ_softID2)[i]==0) continue;

	nTQwithAtLeast2muPtEtaID++;

	if((*recoTQ_Y1vtxprob)[i]<0.1)continue;

	nTQwithAtLeast2muPtEtaIDVtxProb++;


	//electrons
	if(((*recoTQ_charge3)[i]*(*recoTQ_charge4)[i])>0 )continue;

        nTQwithAtLeast2el++;

	if((*recoTQ_isPFoverlap3)[i] || (*recoTQ_isPFoverlap4)[i])continue;

        nTQwithAtLeast2elOver++;

	bool isPt3wrong = ((*recoTQ_isPF3)[i]==1 && (*recoTQ_pt3)[i]<1) || ((*recoTQ_isPF3)[i]==0 && (*recoTQ_pt3mode)[i]<1);
	bool isPt4wrong = ((*recoTQ_isPF4)[i]==1 && (*recoTQ_pt4)[i]<1) || ((*recoTQ_isPF4)[i]==0 && (*recoTQ_pt4mode)[i]<1);
	bool isEta3wrong = fabs((*recoTQ_eta3)[i])>2.5;
	bool isEta4wrong = fabs((*recoTQ_eta4)[i])>2.5;
	bool isID3wrong = ((*recoTQ_isPF3)[i]==1 && (*recoTQ_mvaPFValue3)[i]<-2) || ((*recoTQ_isPF3)[i]==0 && (*recoTQ_mvaValue3)[i]<3);
	bool isID4wrong = ((*recoTQ_isPF4)[i]==1 && (*recoTQ_mvaPFValue4)[i]<-2) || ((*recoTQ_isPF4)[i]==0 && (*recoTQ_mvaValue4)[i]<3);

	if(isPt3wrong || isPt4wrong || isEta3wrong || isEta4wrong || isID3wrong || isID4wrong)continue;
	
        nTQwithAtLeast2elOverPtEtaID++;


	//select TQ with highest vtx prob
	if((*recoTQ_Y2vtxprob)[i]<0.1)continue;

        nTQwithAtLeast2elOverPtEtaIDVtxProb++;


	if((*recoTQ_vtxprob)[i]<0.1)continue;

        nTQwithVtxProb++;

	
	if((*recoTQ_Y2vtxprob)[i]> prob_cand){

	  prob_cand = (*recoTQ_Y2vtxprob)[i];
	  best_cand= i;
	}

	

      }
      

      //increase counters
      if(nTQwithAtLeast2mu>0)counter[1]++;
      if(nTQwithAtLeast2muPtEtaID>0)counter[2]++;
      if(nTQwithAtLeast2muPtEtaIDVtxProb>0)counter[3]++;
      if(nTQwithAtLeast2el>0)counter[4]++;
      if(nTQwithAtLeast2elOver>0)counter[5]++;
      if(nTQwithAtLeast2elOverPtEtaID>0)counter[6]++;
      if(nTQwithAtLeast2elOverPtEtaIDVtxProb>0)counter[7]++;
      if(nTQwithVtxProb>0)counter[8]++;
      
      //fillign reduced tree with candidate infos
      if(best_cand<999){

	e_pt1=(*recoTQ_isPF3)[best_cand]==1?  (*recoTQ_pt3)[best_cand] : (*recoTQ_pt3mode)[best_cand]; 
	e_pt2=(*recoTQ_isPF4)[best_cand]==1?  (*recoTQ_pt4)[best_cand] : (*recoTQ_pt4mode)[best_cand];
	e_eta1=(*recoTQ_eta3)[best_cand];
	e_eta2=(*recoTQ_eta4)[best_cand];
	e_phi1=(*recoTQ_phi3)[best_cand];
	e_phi2=(*recoTQ_phi4)[best_cand];
	e_isPF1=(*recoTQ_isPF3)[best_cand];
	e_isPF2=(*recoTQ_isPF4)[best_cand];


	m_pt1=(*recoTQ_pt1)[best_cand];
	m_pt2=(*recoTQ_pt2)[best_cand];
	m_eta1=(*recoTQ_eta1)[best_cand];
	m_eta2=(*recoTQ_eta2)[best_cand];
	m_phi1=(*recoTQ_phi1)[best_cand];
	m_phi2=(*recoTQ_phi2)[best_cand];

	Ym_vtxprob=(*recoTQ_Y1vtxprob)[best_cand];
	Ym_mass=(*recoTQ_Y1mass)[best_cand]; 
	Ym_pt=(*recoTQ_Y1pt)[best_cand]; 
	Ym_eta=(*recoTQ_Y1eta)[best_cand]; 

	Ye_vtxprob=(*recoTQ_Y2vtxprob)[best_cand];
	Ye_mass=(*recoTQ_Y2mass)[best_cand]; 
	Ye_pt=(*recoTQ_Y2pt)[best_cand]; 
	Ye_eta=(*recoTQ_Y2eta)[best_cand]; 

	TQ_vtxprob=(*recoTQ_vtxprob)[best_cand]; 
	TQ_mass=(*recoTQ_mass)[best_cand]; 
	TQ_pt=(*recoTQ_pt)[best_cand]; 
	TQ_eta=(*recoTQ_eta)[best_cand]; 
       
	tree_red.Fill();
      }
 }

   //printing results
   std::cout<<" tot: "<<counter[0]<<std::endl;
   std::cout<<" 2mu reco: "<<counter[1]<<std::endl;
   std::cout<<" 2mu reco pt eta id: "<<counter[2]<<std::endl;
   std::cout<<" 2mu reco pt eta id vtx: "<<counter[3]<<std::endl;
   std::cout<<" 2e reco: "<<counter[4]<<std::endl;
   std::cout<<" 2e reco no overlap: "<<counter[5]<<std::endl;
   std::cout<<" 2e reco no overlap pt eta id : "<<counter[6]<<std::endl;
   std::cout<<" 2e reco no overlap pt eta id vtx: "<<counter[7]<<std::endl;
   std::cout<<" TQ vtx: "<<counter[8]<<std::endl;
      

   //creating efficiency counters
   TH1F* h_counter = new TH1F("h_counter", "",9,0,9);
   TH1F* eff_counter = new TH1F("eff_counter", "",9,0,9);
   TH1F* effrel_counter = new TH1F("effrel_counter", "",9,0,9);
   for(int i=0;i<9;i++) h_counter->SetBinContent(i+1, counter[i]);
   for(int i=0;i<9;i++) {
     double eff= (float)counter[i]/counter[0];
     eff_counter->SetBinContent(i+1, eff);
     eff_counter->SetBinError(i+1, sqrt((eff*(1-eff))/counter[0]));
   }

   for(int i=0;i<9;i++) {
     double eff;
     if(i!=0) eff= (float)counter[i]/counter[i-1];
     else eff=1;
     double efferr;
     if(i!=0)efferr= sqrt((eff*(1-eff))/counter[i-1]);
     else efferr=0.;
     effrel_counter->SetBinContent(i+1, eff);
     effrel_counter->SetBinError(i+1, efferr);
   }
   
   h_counter->GetXaxis()->SetBinLabel(1 ," Total");
   h_counter->GetXaxis()->SetBinLabel(2 ," nDiMu >1");
   h_counter->GetXaxis()->SetBinLabel(3 ,"softID==1");
   h_counter->GetXaxis()->SetBinLabel(4 ,"#mu#mu vtx prob >0.1");
   h_counter->GetXaxis()->SetBinLabel(5,"nDiEle>1");
   h_counter->GetXaxis()->SetBinLabel(6,"no overlap");
   h_counter->GetXaxis()->SetBinLabel(7,"ID==1");
   h_counter->GetXaxis()->SetBinLabel(8 ,"ee vtx prob >0.1");
   h_counter->GetXaxis()->SetBinLabel(9 ,"TQ vtx prob >0.1");

   eff_counter->GetXaxis()->SetBinLabel(1 ," Total");
   eff_counter->GetXaxis()->SetBinLabel(2 ," nDiMu >1");
   eff_counter->GetXaxis()->SetBinLabel(3 ,"softID==1");
   eff_counter->GetXaxis()->SetBinLabel(4 ,"#mu#mu vtx prob >0.1");
   eff_counter->GetXaxis()->SetBinLabel(5,"nDiEle>0");
   eff_counter->GetXaxis()->SetBinLabel(6,"no overlap");
   eff_counter->GetXaxis()->SetBinLabel(7,"ID==1");
   eff_counter->GetXaxis()->SetBinLabel(8 ,"ee vtx prob >0.1");
   eff_counter->GetXaxis()->SetBinLabel(9 ,"TQ vtx prob >0.1");

   effrel_counter->GetXaxis()->SetBinLabel(1 ," Total");
   effrel_counter->GetXaxis()->SetBinLabel(2 ," nDiMu >0");
   effrel_counter->GetXaxis()->SetBinLabel(3 ,"softID==1");
   effrel_counter->GetXaxis()->SetBinLabel(4 ,"#mu#mu vtx prob >0.1");
   effrel_counter->GetXaxis()->SetBinLabel(5,"nDiEle>1");
   effrel_counter->GetXaxis()->SetBinLabel(6,"no overlap");
   effrel_counter->GetXaxis()->SetBinLabel(7,"ID==1");
   effrel_counter->GetXaxis()->SetBinLabel(8 ,"ee vtx prob >0.1");
   effrel_counter->GetXaxis()->SetBinLabel(9 ,"TQ vtx prob >0.1");

   //saving on file   
   TFile* fout = new TFile(("fout_m"+mass+"_TQ.root").c_str(),"RECREATE");
   fout->cd();
   tree_red.Write();
   h_counter->Write();
   eff_counter->Write();
   effrel_counter->Write();
   fout->Write();
   fout->Close();
}




