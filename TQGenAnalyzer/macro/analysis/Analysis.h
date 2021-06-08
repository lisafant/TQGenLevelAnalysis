//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May  5 16:37:50 2021 by ROOT version 6.14/09
// from TTree tree/tree
// found on file: ../../crab/ntuple_tetraquarks_wLowPt.root
//////////////////////////////////////////////////////////

#ifndef Analysis_h
#define Analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class Analysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           sampleID;
   Float_t         xsec;
   Int_t           run;
   Int_t           event;
   Int_t           lumi;
   Int_t           nEle;
   Int_t           nMu;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Int_t           nvtx;
   Float_t         rho;
   Int_t           npu;
   Float_t         puw_2016;
   Float_t         puw_2017;
   Float_t         puw_2018;
   Float_t         puw_ALL;
   Float_t         HLT2016_Dimuon0_Jpsi_Muon;
   Float_t         HLT2016_Dimuon0_Upsilon_Muon;
   Float_t         HLT2016_Mu8_DiEle12_CaloIdL_TrackIdL;
   Float_t         HLT2016_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;
   Float_t         HLT2017_Dimuon0_Jpsi3p5_Muon2;
   Float_t         HLT2017_Dimuon12_Upsilon_eta1p5;
   Float_t         HLT2017_Mu8_DiEle12_CaloIdL_TrackIdL;
   Float_t         HLT2017_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;
   Float_t         HLT2018_Dimuon0_Jpsi3p5_Muon2;
   Float_t         HLT2018_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon;
   Float_t         HLT2018_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL;
   Float_t         HLT2018_Mu8_DiEle12_CaloIdL_TrackIdL;
   Float_t         HLT2018_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;
   Float_t         TQ_genMass;
   vector<float>   *genLep_pt;
   vector<float>   *genLep_mass;
   vector<float>   *genLep_eta;
   vector<float>   *genLep_phi;
   vector<int>     *genLep_pdgId;
   vector<float>   *genMom_pt;
   vector<float>   *genMom_mass;
   vector<float>   *genMom_eta;
   vector<float>   *genMom_phi;
   vector<int>     *genMom_pdgId;
   Int_t           nMuReco;
   vector<float>   *recoMu_pt;
   vector<float>   *recoMu_mass;
   vector<float>   *recoMu_eta;
   vector<float>   *recoMu_phi;
   vector<int>     *recoMu_charge;
   vector<float>   *recoMu_dR1;
   vector<int>     *recoMu_matchid1;
   vector<float>   *recoMu_dR2;
   vector<int>     *recoMu_matchid2;
   vector<int>     *recoMu_looseID;
   vector<int>     *recoMu_softID;
   vector<int>     *recoMu_nTrkLayers;
   vector<int>     *recoMu_nPixLayers;
   vector<int>     *recoMu_isHQ;
   vector<float>   *recoMu_dxy;
   vector<float>   *recoMu_dz;
   Int_t           nDimuReco;
   vector<float>   *recoDimu_vx;
   vector<float>   *recoDimu_vy;
   vector<float>   *recoDimu_vz;
   vector<float>   *recoDimu_vtxchi2;
   vector<float>   *recoDimu_vtxndof;
   vector<float>   *recoDimu_vtxprob;
   vector<float>   *recoDimu_pt1;
   vector<float>   *recoDimu_eta1;
   vector<float>   *recoDimu_phi1;
   vector<float>   *recoDimu_charge1;
   vector<float>   *recoDimu_mass1;
   vector<int>     *recoDimu_softID1;
   vector<float>   *recoDimu_pt2;
   vector<float>   *recoDimu_eta2;
   vector<float>   *recoDimu_phi2;
   vector<float>   *recoDimu_charge2;
   vector<float>   *recoDimu_mass2;
   vector<int>     *recoDimu_softID2;
   vector<float>   *recoDimu_pt;
   vector<float>   *recoDimu_eta;
   vector<float>   *recoDimu_phi;
   vector<float>   *recoDimu_mass;
   vector<float>   *recoDimu_massErr;
   Int_t           nEleReco;
   Int_t           nPFEleReco;
   Int_t           nLowPtEleReco;
   vector<float>   *recoEle_pt;
   vector<float>   *recoEle_mass;
   vector<float>   *recoEle_eta;
   vector<float>   *recoEle_phi;
   vector<int>     *recoEle_charge;
   vector<float>   *recoEle_ptmode;
   vector<float>   *recoEle_etamode;
   vector<float>   *recoEle_phimode;
   vector<float>   *recoEle_p;
   vector<float>   *recoEle_vx;
   vector<float>   *recoEle_vy;
   vector<float>   *recoEle_vz;
   vector<float>   *recoEle_dR1;
   vector<int>     *recoEle_matchid1;
   vector<float>   *recoEle_dR2;
   vector<int>     *recoEle_matchid2;
   vector<float>   *recoEle_E_ecal_preReg;
   vector<float>   *recoEle_E_ecal_postReg;
   vector<float>   *recoEle_E_ecaltrk_preReg;
   vector<float>   *recoEle_E_ecaltrk_postReg;
   vector<float>   *recoEle_rawE;
   vector<float>   *recoEle_corrEcalE;
   vector<float>   *recoEle_gsfTrkChi2;
   vector<float>   *recoEle_passConvVeto;
   vector<float>   *recoEle_isPF;
   vector<float>   *recoEle_isLowPt;
   vector<float>   *recoEle_isPFoverlap;
   vector<float>   *recoEle_mvaPFValue;
   vector<float>   *recoEle_mvaValue;
   Int_t           nDieleReco;
   vector<float>   *recoDiele_vx;
   vector<float>   *recoDiele_vy;
   vector<float>   *recoDiele_vz;
   vector<float>   *recoDiele_vtxchi2;
   vector<float>   *recoDiele_vtxndof;
   vector<float>   *recoDiele_vtxprob;
   vector<float>   *recoDiele_pt1;
   vector<float>   *recoDiele_eta1;
   vector<float>   *recoDiele_phi1;
   vector<float>   *recoDiele_charge1;
   vector<float>   *recoDiele_mass1;
   vector<float>   *recoDiele_mvaValue1;
   vector<float>   *recoDiele_mvaPFValue1;
   vector<int>     *recoDiele_isPF1;
   vector<int>     *recoDiele_isLowPt1;
   vector<int>     *recoDiele_isPFoverlap1;
   vector<float>   *recoDiele_pt2;
   vector<float>   *recoDiele_eta2;
   vector<float>   *recoDiele_phi2;
   vector<float>   *recoDiele_charge2;
   vector<float>   *recoDiele_mass2;
   vector<float>   *recoDiele_mvaValue2;
   vector<float>   *recoDiele_mvaPFValue2;
   vector<int>     *recoDiele_isPF2;
   vector<int>     *recoDiele_isLowPt2;
   vector<int>     *recoDiele_isPFoverlap2;
   vector<float>   *recoDiele_pt;
   vector<float>   *recoDiele_eta;
   vector<float>   *recoDiele_phi;
   vector<float>   *recoDiele_mass;
   vector<float>   *recoDiele_massErr;
   vector<float>   *recoDiele_pt1mode;
   vector<float>   *recoDiele_eta1mode;
   vector<float>   *recoDiele_phi1mode;
   vector<float>   *recoDiele_pt2mode;
   vector<float>   *recoDiele_eta2mode;
   vector<float>   *recoDiele_phi2mode;

   // List of branches
   TBranch        *b_sampleID;   //!
   TBranch        *b_xsec;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_puw_2016;   //!
   TBranch        *b_puw_2017;   //!
   TBranch        *b_puw_2018;   //!
   TBranch        *b_puw_ALL;   //!
   TBranch        *b_HLT2016_Dimuon0_Jpsi_Muon;   //!
   TBranch        *b_HLT2016_Dimuon0_Upsilon_Muon;   //!
   TBranch        *b_HLT2016_Mu8_DiEle12_CaloIdL_TrackIdL;   //!
   TBranch        *b_HLT2016_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;   //!
   TBranch        *b_HLT2017_Dimuon0_Jpsi3p5_Muon2;   //!
   TBranch        *b_HLT2017_Dimuon12_Upsilon_eta1p5;   //!
   TBranch        *b_HLT2017_Mu8_DiEle12_CaloIdL_TrackIdL;   //!
   TBranch        *b_HLT2017_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;   //!
   TBranch        *b_HLT2018_Dimuon0_Jpsi3p5_Muon2;   //!
   TBranch        *b_HLT2018_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon;   //!
   TBranch        *b_HLT2018_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL;   //!
   TBranch        *b_HLT2018_Mu8_DiEle12_CaloIdL_TrackIdL;   //!
   TBranch        *b_HLT2018_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;   //!
   TBranch        *b_TQ_genMass;   //!
   TBranch        *b_genLep_pt;   //!
   TBranch        *b_genLep_mass;   //!
   TBranch        *b_genLep_eta;   //!
   TBranch        *b_genLep_phi;   //!
   TBranch        *b_genLep_pdgId;   //!
   TBranch        *b_genMom_pt;   //!
   TBranch        *b_genMom_mass;   //!
   TBranch        *b_genMom_eta;   //!
   TBranch        *b_genMom_phi;   //!
   TBranch        *b_genMom_pdgId;   //!
   TBranch        *b_nMuReco;   //!
   TBranch        *b_recoMu_pt;   //!
   TBranch        *b_recoMu_mass;   //!
   TBranch        *b_recoMu_eta;   //!
   TBranch        *b_recoMu_phi;   //!
   TBranch        *b_recoMu_charge;   //!
   TBranch        *b_recoMu_dR1;   //!
   TBranch        *b_recoMu_matchid1;   //!
   TBranch        *b_recoMu_dR2;   //!
   TBranch        *b_recoMu_matchid2;   //!
   TBranch        *b_recoMu_looseID;   //!
   TBranch        *b_recoMu_softID;   //!
   TBranch        *b_recoMu_nTrkLayers;   //!
   TBranch        *b_recoMu_nPixLayers;   //!
   TBranch        *b_recoMu_isHQ;   //!
   TBranch        *b_recoMu_dxy;   //!
   TBranch        *b_recoMu_dz;   //!
   TBranch        *b_nDimuReco;   //!
   TBranch        *b_recoDimu_vx;   //!
   TBranch        *b_recoDimu_vy;   //!
   TBranch        *b_recoDimu_vz;   //!
   TBranch        *b_recoDimu_vtxchi2;   //!
   TBranch        *b_recoDimu_vtxndof;   //!
   TBranch        *b_recoDimu_vtxprob;   //!
   TBranch        *b_recoDimu_pt1;   //!
   TBranch        *b_recoDimu_eta1;   //!
   TBranch        *b_recoDimu_phi1;   //!
   TBranch        *b_recoDimu_charge1;   //!
   TBranch        *b_recoDimu_mass1;   //!
   TBranch        *b_recoDimu_softID1;   //!
   TBranch        *b_recoDimu_pt2;   //!
   TBranch        *b_recoDimu_eta2;   //!
   TBranch        *b_recoDimu_phi2;   //!
   TBranch        *b_recoDimu_charge2;   //!
   TBranch        *b_recoDimu_mass2;   //!
   TBranch        *b_recoDimu_softID2;   //!
   TBranch        *b_recoDimu_pt;   //!
   TBranch        *b_recoDimu_eta;   //!
   TBranch        *b_recoDimu_phi;   //!
   TBranch        *b_recoDimu_mass;   //!
   TBranch        *b_recoDimu_massErr;   //!
   TBranch        *b_nEleReco;   //!
   TBranch        *b_nPFEleReco;   //!
   TBranch        *b_nLowPtEleReco;   //!
   TBranch        *b_recoEle_pt;   //!
   TBranch        *b_recoEle_mass;   //!
   TBranch        *b_recoEle_eta;   //!
   TBranch        *b_recoEle_phi;   //!
   TBranch        *b_recoEle_charge;   //!
   TBranch        *b_recoEle_ptmode;   //!
   TBranch        *b_recoEle_etamode;   //!
   TBranch        *b_recoEle_phimode;   //!
   TBranch        *b_recoEle_p;   //!
   TBranch        *b_recoEle_vx;   //!
   TBranch        *b_recoEle_vy;   //!
   TBranch        *b_recoEle_vz;   //!
   TBranch        *b_recoEle_dR1;   //!
   TBranch        *b_recoEle_matchid1;   //!
   TBranch        *b_recoEle_dR2;   //!
   TBranch        *b_recoEle_matchid2;   //!
   TBranch        *b_recoEle_E_ecal_preReg;   //!
   TBranch        *b_recoEle_E_ecal_postReg;   //!
   TBranch        *b_recoEle_E_ecaltrk_preReg;   //!
   TBranch        *b_recoEle_E_ecaltrk_postReg;   //!
   TBranch        *b_recoEle_rawE;   //!
   TBranch        *b_recoEle_corrEcalE;   //!
   TBranch        *b_recoEle_gsfTrkChi2;   //!
   TBranch        *b_recoEle_passConvVeto;   //!
   TBranch        *b_recoEle_isPF;   //!
   TBranch        *b_recoEle_isLowPt;   //!
   TBranch        *b_recoEle_isPFoverlap;   //!
   TBranch        *b_recoEle_mvaPFValue;   //!
   TBranch        *b_recoEle_mvaValue;   //!
   TBranch        *b_nDieleReco;   //!
   TBranch        *b_recoDiele_vx;   //!
   TBranch        *b_recoDiele_vy;   //!
   TBranch        *b_recoDiele_vz;   //!
   TBranch        *b_recoDiele_vtxchi2;   //!
   TBranch        *b_recoDiele_vtxndof;   //!
   TBranch        *b_recoDiele_vtxprob;   //!
   TBranch        *b_recoDiele_pt1;   //!
   TBranch        *b_recoDiele_eta1;   //!
   TBranch        *b_recoDiele_phi1;   //!
   TBranch        *b_recoDiele_charge1;   //!
   TBranch        *b_recoDiele_mass1;   //!
   TBranch        *b_recoDiele_mvaValue1;   //!
   TBranch        *b_recoDiele_mvaPFValue1;   //!
   TBranch        *b_recoDiele_isPF1;   //!
   TBranch        *b_recoDiele_isLowPt1;   //!
   TBranch        *b_recoDiele_isPFoverlap1;   //!
   TBranch        *b_recoDiele_pt2;   //!
   TBranch        *b_recoDiele_eta2;   //!
   TBranch        *b_recoDiele_phi2;   //!
   TBranch        *b_recoDiele_charge2;   //!
   TBranch        *b_recoDiele_mass2;   //!
   TBranch        *b_recoDiele_mvaValue2;   //!
   TBranch        *b_recoDiele_mvaPFValue2;   //!
   TBranch        *b_recoDiele_isPF2;   //!
   TBranch        *b_recoDiele_isLowPt2;   //!
   TBranch        *b_recoDiele_isPFoverlap2;   //!
   TBranch        *b_recoDiele_pt;   //!
   TBranch        *b_recoDiele_eta;   //!
   TBranch        *b_recoDiele_phi;   //!
   TBranch        *b_recoDiele_mass;   //!
   TBranch        *b_recoDiele_massErr;   //!
   TBranch        *b_recoDiele_pt1mode;   //!
   TBranch        *b_recoDiele_eta1mode;   //!
   TBranch        *b_recoDiele_phi1mode;   //!
   TBranch        *b_recoDiele_pt2mode;   //!
   TBranch        *b_recoDiele_eta2mode;   //!
   TBranch        *b_recoDiele_phi2mode;   //!

   Analysis(TTree *tree=0);
   virtual ~Analysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string mass);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Analysis_cxx
Analysis::Analysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../crab/ntuple_tetraquarks_wLowPt.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../crab/ntuple_tetraquarks_wLowPt.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("../../crab/ntuple_tetraquarks_wLowPt.root:/GenAnalysis");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

Analysis::~Analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Analysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   genLep_pt = 0;
   genLep_mass = 0;
   genLep_eta = 0;
   genLep_phi = 0;
   genLep_pdgId = 0;
   genMom_pt = 0;
   genMom_mass = 0;
   genMom_eta = 0;
   genMom_phi = 0;
   genMom_pdgId = 0;
   recoMu_pt = 0;
   recoMu_mass = 0;
   recoMu_eta = 0;
   recoMu_phi = 0;
   recoMu_charge = 0;
   recoMu_dR1 = 0;
   recoMu_matchid1 = 0;
   recoMu_dR2 = 0;
   recoMu_matchid2 = 0;
   recoMu_looseID = 0;
   recoMu_softID = 0;
   recoMu_nTrkLayers = 0;
   recoMu_nPixLayers = 0;
   recoMu_isHQ = 0;
   recoMu_dxy = 0;
   recoMu_dz = 0;
   recoDimu_vx = 0;
   recoDimu_vy = 0;
   recoDimu_vz = 0;
   recoDimu_vtxchi2 = 0;
   recoDimu_vtxndof = 0;
   recoDimu_vtxprob = 0;
   recoDimu_pt1 = 0;
   recoDimu_eta1 = 0;
   recoDimu_phi1 = 0;
   recoDimu_charge1 = 0;
   recoDimu_mass1 = 0;
   recoDimu_softID1 = 0;
   recoDimu_pt2 = 0;
   recoDimu_eta2 = 0;
   recoDimu_phi2 = 0;
   recoDimu_charge2 = 0;
   recoDimu_mass2 = 0;
   recoDimu_softID2 = 0;
   recoDimu_pt = 0;
   recoDimu_eta = 0;
   recoDimu_phi = 0;
   recoDimu_mass = 0;
   recoDimu_massErr = 0;
   recoEle_pt = 0;
   recoEle_mass = 0;
   recoEle_eta = 0;
   recoEle_phi = 0;
   recoEle_charge = 0;
   recoEle_ptmode = 0;
   recoEle_etamode = 0;
   recoEle_phimode = 0;
   recoEle_p = 0;
   recoEle_vx = 0;
   recoEle_vy = 0;
   recoEle_vz = 0;
   recoEle_dR1 = 0;
   recoEle_matchid1 = 0;
   recoEle_dR2 = 0;
   recoEle_matchid2 = 0;
   recoEle_E_ecal_preReg = 0;
   recoEle_E_ecal_postReg = 0;
   recoEle_E_ecaltrk_preReg = 0;
   recoEle_E_ecaltrk_postReg = 0;
   recoEle_rawE = 0;
   recoEle_corrEcalE = 0;
   recoEle_gsfTrkChi2 = 0;
   recoEle_passConvVeto = 0;
   recoEle_isPF = 0;
   recoEle_isLowPt = 0;
   recoEle_isPFoverlap = 0;
   recoEle_mvaPFValue = 0;
   recoEle_mvaValue = 0;
   recoDiele_vx = 0;
   recoDiele_vy = 0;
   recoDiele_vz = 0;
   recoDiele_vtxchi2 = 0;
   recoDiele_vtxndof = 0;
   recoDiele_vtxprob = 0;
   recoDiele_pt1 = 0;
   recoDiele_eta1 = 0;
   recoDiele_phi1 = 0;
   recoDiele_charge1 = 0;
   recoDiele_mass1 = 0;
   recoDiele_mvaValue1 = 0;
   recoDiele_mvaPFValue1 = 0;
   recoDiele_isPF1 = 0;
   recoDiele_isLowPt1 = 0;
   recoDiele_isPFoverlap1 = 0;
   recoDiele_pt2 = 0;
   recoDiele_eta2 = 0;
   recoDiele_phi2 = 0;
   recoDiele_charge2 = 0;
   recoDiele_mass2 = 0;
   recoDiele_mvaValue2 = 0;
   recoDiele_mvaPFValue2 = 0;
   recoDiele_isPF2 = 0;
   recoDiele_isLowPt2 = 0;
   recoDiele_isPFoverlap2 = 0;
   recoDiele_pt = 0;
   recoDiele_eta = 0;
   recoDiele_phi = 0;
   recoDiele_mass = 0;
   recoDiele_massErr = 0;
   recoDiele_pt1mode = 0;
   recoDiele_eta1mode = 0;
   recoDiele_phi1mode = 0;
   recoDiele_pt2mode = 0;
   recoDiele_eta2mode = 0;
   recoDiele_phi2mode = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("sampleID", &sampleID, &b_sampleID);
   fChain->SetBranchAddress("xsec", &xsec, &b_xsec);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("puw_2016", &puw_2016, &b_puw_2016);
   fChain->SetBranchAddress("puw_2017", &puw_2017, &b_puw_2017);
   fChain->SetBranchAddress("puw_2018", &puw_2018, &b_puw_2018);
   fChain->SetBranchAddress("puw_ALL", &puw_ALL, &b_puw_ALL);
   fChain->SetBranchAddress("HLT2016_Dimuon0_Jpsi_Muon", &HLT2016_Dimuon0_Jpsi_Muon, &b_HLT2016_Dimuon0_Jpsi_Muon);
   fChain->SetBranchAddress("HLT2016_Dimuon0_Upsilon_Muon", &HLT2016_Dimuon0_Upsilon_Muon, &b_HLT2016_Dimuon0_Upsilon_Muon);
   fChain->SetBranchAddress("HLT2016_Mu8_DiEle12_CaloIdL_TrackIdL", &HLT2016_Mu8_DiEle12_CaloIdL_TrackIdL, &b_HLT2016_Mu8_DiEle12_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT2016_DiMu9_Ele9_CaloIdL_TrackIdL_DZ", &HLT2016_DiMu9_Ele9_CaloIdL_TrackIdL_DZ, &b_HLT2016_DiMu9_Ele9_CaloIdL_TrackIdL_DZ);
   fChain->SetBranchAddress("HLT2017_Dimuon0_Jpsi3p5_Muon2", &HLT2017_Dimuon0_Jpsi3p5_Muon2, &b_HLT2017_Dimuon0_Jpsi3p5_Muon2);
   fChain->SetBranchAddress("HLT2017_Dimuon12_Upsilon_eta1p5", &HLT2017_Dimuon12_Upsilon_eta1p5, &b_HLT2017_Dimuon12_Upsilon_eta1p5);
   fChain->SetBranchAddress("HLT2017_Mu8_DiEle12_CaloIdL_TrackIdL", &HLT2017_Mu8_DiEle12_CaloIdL_TrackIdL, &b_HLT2017_Mu8_DiEle12_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT2017_DiMu9_Ele9_CaloIdL_TrackIdL_DZ", &HLT2017_DiMu9_Ele9_CaloIdL_TrackIdL_DZ, &b_HLT2017_DiMu9_Ele9_CaloIdL_TrackIdL_DZ);
   fChain->SetBranchAddress("HLT2018_Dimuon0_Jpsi3p5_Muon2", &HLT2018_Dimuon0_Jpsi3p5_Muon2, &b_HLT2018_Dimuon0_Jpsi3p5_Muon2);
   fChain->SetBranchAddress("HLT2018_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon", &HLT2018_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon, &b_HLT2018_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon);
   fChain->SetBranchAddress("HLT2018_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL", &HLT2018_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL, &b_HLT2018_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT2018_Mu8_DiEle12_CaloIdL_TrackIdL", &HLT2018_Mu8_DiEle12_CaloIdL_TrackIdL, &b_HLT2018_Mu8_DiEle12_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT2018_DiMu9_Ele9_CaloIdL_TrackIdL_DZ", &HLT2018_DiMu9_Ele9_CaloIdL_TrackIdL_DZ, &b_HLT2018_DiMu9_Ele9_CaloIdL_TrackIdL_DZ);
   fChain->SetBranchAddress("TQ_genMass", &TQ_genMass, &b_TQ_genMass);
   fChain->SetBranchAddress("genLep_pt", &genLep_pt, &b_genLep_pt);
   fChain->SetBranchAddress("genLep_mass", &genLep_mass, &b_genLep_mass);
   fChain->SetBranchAddress("genLep_eta", &genLep_eta, &b_genLep_eta);
   fChain->SetBranchAddress("genLep_phi", &genLep_phi, &b_genLep_phi);
   fChain->SetBranchAddress("genLep_pdgId", &genLep_pdgId, &b_genLep_pdgId);
   fChain->SetBranchAddress("genMom_pt", &genMom_pt, &b_genMom_pt);
   fChain->SetBranchAddress("genMom_mass", &genMom_mass, &b_genMom_mass);
   fChain->SetBranchAddress("genMom_eta", &genMom_eta, &b_genMom_eta);
   fChain->SetBranchAddress("genMom_phi", &genMom_phi, &b_genMom_phi);
   fChain->SetBranchAddress("genMom_pdgId", &genMom_pdgId, &b_genMom_pdgId);
   fChain->SetBranchAddress("nMuReco", &nMuReco, &b_nMuReco);
   fChain->SetBranchAddress("recoMu_pt", &recoMu_pt, &b_recoMu_pt);
   fChain->SetBranchAddress("recoMu_mass", &recoMu_mass, &b_recoMu_mass);
   fChain->SetBranchAddress("recoMu_eta", &recoMu_eta, &b_recoMu_eta);
   fChain->SetBranchAddress("recoMu_phi", &recoMu_phi, &b_recoMu_phi);
   fChain->SetBranchAddress("recoMu_charge", &recoMu_charge, &b_recoMu_charge);
   fChain->SetBranchAddress("recoMu_dR1", &recoMu_dR1, &b_recoMu_dR1);
   fChain->SetBranchAddress("recoMu_matchid1", &recoMu_matchid1, &b_recoMu_matchid1);
   fChain->SetBranchAddress("recoMu_dR2", &recoMu_dR2, &b_recoMu_dR2);
   fChain->SetBranchAddress("recoMu_matchid2", &recoMu_matchid2, &b_recoMu_matchid2);
   fChain->SetBranchAddress("recoMu_looseID", &recoMu_looseID, &b_recoMu_looseID);
   fChain->SetBranchAddress("recoMu_softID", &recoMu_softID, &b_recoMu_softID);
   fChain->SetBranchAddress("recoMu_nTrkLayers", &recoMu_nTrkLayers, &b_recoMu_nTrkLayers);
   fChain->SetBranchAddress("recoMu_nPixLayers", &recoMu_nPixLayers, &b_recoMu_nPixLayers);
   fChain->SetBranchAddress("recoMu_isHQ", &recoMu_isHQ, &b_recoMu_isHQ);
   fChain->SetBranchAddress("recoMu_dxy", &recoMu_dxy, &b_recoMu_dxy);
   fChain->SetBranchAddress("recoMu_dz", &recoMu_dz, &b_recoMu_dz);
   fChain->SetBranchAddress("nDimuReco", &nDimuReco, &b_nDimuReco);
   fChain->SetBranchAddress("recoDimu_vx", &recoDimu_vx, &b_recoDimu_vx);
   fChain->SetBranchAddress("recoDimu_vy", &recoDimu_vy, &b_recoDimu_vy);
   fChain->SetBranchAddress("recoDimu_vz", &recoDimu_vz, &b_recoDimu_vz);
   fChain->SetBranchAddress("recoDimu_vtxchi2", &recoDimu_vtxchi2, &b_recoDimu_vtxchi2);
   fChain->SetBranchAddress("recoDimu_vtxndof", &recoDimu_vtxndof, &b_recoDimu_vtxndof);
   fChain->SetBranchAddress("recoDimu_vtxprob", &recoDimu_vtxprob, &b_recoDimu_vtxprob);
   fChain->SetBranchAddress("recoDimu_pt1", &recoDimu_pt1, &b_recoDimu_pt1);
   fChain->SetBranchAddress("recoDimu_eta1", &recoDimu_eta1, &b_recoDimu_eta1);
   fChain->SetBranchAddress("recoDimu_phi1", &recoDimu_phi1, &b_recoDimu_phi1);
   fChain->SetBranchAddress("recoDimu_charge1", &recoDimu_charge1, &b_recoDimu_charge1);
   fChain->SetBranchAddress("recoDimu_mass1", &recoDimu_mass1, &b_recoDimu_mass1);
   fChain->SetBranchAddress("recoDimu_softID1", &recoDimu_softID1, &b_recoDimu_softID1);
   fChain->SetBranchAddress("recoDimu_pt2", &recoDimu_pt2, &b_recoDimu_pt2);
   fChain->SetBranchAddress("recoDimu_eta2", &recoDimu_eta2, &b_recoDimu_eta2);
   fChain->SetBranchAddress("recoDimu_phi2", &recoDimu_phi2, &b_recoDimu_phi2);
   fChain->SetBranchAddress("recoDimu_charge2", &recoDimu_charge2, &b_recoDimu_charge2);
   fChain->SetBranchAddress("recoDimu_mass2", &recoDimu_mass2, &b_recoDimu_mass2);
   fChain->SetBranchAddress("recoDimu_softID2", &recoDimu_softID2, &b_recoDimu_softID2);
   fChain->SetBranchAddress("recoDimu_pt", &recoDimu_pt, &b_recoDimu_pt);
   fChain->SetBranchAddress("recoDimu_eta", &recoDimu_eta, &b_recoDimu_eta);
   fChain->SetBranchAddress("recoDimu_phi", &recoDimu_phi, &b_recoDimu_phi);
   fChain->SetBranchAddress("recoDimu_mass", &recoDimu_mass, &b_recoDimu_mass);
   fChain->SetBranchAddress("recoDimu_massErr", &recoDimu_massErr, &b_recoDimu_massErr);
   fChain->SetBranchAddress("nEleReco", &nEleReco, &b_nEleReco);
   fChain->SetBranchAddress("nPFEleReco", &nPFEleReco, &b_nPFEleReco);
   fChain->SetBranchAddress("nLowPtEleReco", &nLowPtEleReco, &b_nLowPtEleReco);
   fChain->SetBranchAddress("recoEle_pt", &recoEle_pt, &b_recoEle_pt);
   fChain->SetBranchAddress("recoEle_mass", &recoEle_mass, &b_recoEle_mass);
   fChain->SetBranchAddress("recoEle_eta", &recoEle_eta, &b_recoEle_eta);
   fChain->SetBranchAddress("recoEle_phi", &recoEle_phi, &b_recoEle_phi);
   fChain->SetBranchAddress("recoEle_charge", &recoEle_charge, &b_recoEle_charge);
   fChain->SetBranchAddress("recoEle_ptmode", &recoEle_ptmode, &b_recoEle_ptmode);
   fChain->SetBranchAddress("recoEle_etamode", &recoEle_etamode, &b_recoEle_etamode);
   fChain->SetBranchAddress("recoEle_phimode", &recoEle_phimode, &b_recoEle_phimode);
   fChain->SetBranchAddress("recoEle_p", &recoEle_p, &b_recoEle_p);
   fChain->SetBranchAddress("recoEle_vx", &recoEle_vx, &b_recoEle_vx);
   fChain->SetBranchAddress("recoEle_vy", &recoEle_vy, &b_recoEle_vy);
   fChain->SetBranchAddress("recoEle_vz", &recoEle_vz, &b_recoEle_vz);
   fChain->SetBranchAddress("recoEle_dR1", &recoEle_dR1, &b_recoEle_dR1);
   fChain->SetBranchAddress("recoEle_matchid1", &recoEle_matchid1, &b_recoEle_matchid1);
   fChain->SetBranchAddress("recoEle_dR2", &recoEle_dR2, &b_recoEle_dR2);
   fChain->SetBranchAddress("recoEle_matchid2", &recoEle_matchid2, &b_recoEle_matchid2);
   fChain->SetBranchAddress("recoEle_E_ecal_preReg", &recoEle_E_ecal_preReg, &b_recoEle_E_ecal_preReg);
   fChain->SetBranchAddress("recoEle_E_ecal_postReg", &recoEle_E_ecal_postReg, &b_recoEle_E_ecal_postReg);
   fChain->SetBranchAddress("recoEle_E_ecaltrk_preReg", &recoEle_E_ecaltrk_preReg, &b_recoEle_E_ecaltrk_preReg);
   fChain->SetBranchAddress("recoEle_E_ecaltrk_postReg", &recoEle_E_ecaltrk_postReg, &b_recoEle_E_ecaltrk_postReg);
   fChain->SetBranchAddress("recoEle_rawE", &recoEle_rawE, &b_recoEle_rawE);
   fChain->SetBranchAddress("recoEle_corrEcalE", &recoEle_corrEcalE, &b_recoEle_corrEcalE);
   fChain->SetBranchAddress("recoEle_gsfTrkChi2", &recoEle_gsfTrkChi2, &b_recoEle_gsfTrkChi2);
   fChain->SetBranchAddress("recoEle_passConvVeto", &recoEle_passConvVeto, &b_recoEle_passConvVeto);
   fChain->SetBranchAddress("recoEle_isPF", &recoEle_isPF, &b_recoEle_isPF);
   fChain->SetBranchAddress("recoEle_isLowPt", &recoEle_isLowPt, &b_recoEle_isLowPt);
   fChain->SetBranchAddress("recoEle_isPFoverlap", &recoEle_isPFoverlap, &b_recoEle_isPFoverlap);
   fChain->SetBranchAddress("recoEle_mvaPFValue", &recoEle_mvaPFValue, &b_recoEle_mvaPFValue);
   fChain->SetBranchAddress("recoEle_mvaValue", &recoEle_mvaValue, &b_recoEle_mvaValue);
   fChain->SetBranchAddress("nDieleReco", &nDieleReco, &b_nDieleReco);
   fChain->SetBranchAddress("recoDiele_vx", &recoDiele_vx, &b_recoDiele_vx);
   fChain->SetBranchAddress("recoDiele_vy", &recoDiele_vy, &b_recoDiele_vy);
   fChain->SetBranchAddress("recoDiele_vz", &recoDiele_vz, &b_recoDiele_vz);
   fChain->SetBranchAddress("recoDiele_vtxchi2", &recoDiele_vtxchi2, &b_recoDiele_vtxchi2);
   fChain->SetBranchAddress("recoDiele_vtxndof", &recoDiele_vtxndof, &b_recoDiele_vtxndof);
   fChain->SetBranchAddress("recoDiele_vtxprob", &recoDiele_vtxprob, &b_recoDiele_vtxprob);
   fChain->SetBranchAddress("recoDiele_pt1", &recoDiele_pt1, &b_recoDiele_pt1);
   fChain->SetBranchAddress("recoDiele_eta1", &recoDiele_eta1, &b_recoDiele_eta1);
   fChain->SetBranchAddress("recoDiele_phi1", &recoDiele_phi1, &b_recoDiele_phi1);
   fChain->SetBranchAddress("recoDiele_charge1", &recoDiele_charge1, &b_recoDiele_charge1);
   fChain->SetBranchAddress("recoDiele_mass1", &recoDiele_mass1, &b_recoDiele_mass1);
   fChain->SetBranchAddress("recoDiele_mvaValue1", &recoDiele_mvaValue1, &b_recoDiele_mvaValue1);
   fChain->SetBranchAddress("recoDiele_mvaPFValue1", &recoDiele_mvaPFValue1, &b_recoDiele_mvaPFValue1);
   fChain->SetBranchAddress("recoDiele_isPF1", &recoDiele_isPF1, &b_recoDiele_isPF1);
   fChain->SetBranchAddress("recoDiele_isLowPt1", &recoDiele_isLowPt1, &b_recoDiele_isLowPt1);
   fChain->SetBranchAddress("recoDiele_isPFoverlap1", &recoDiele_isPFoverlap1, &b_recoDiele_isPFoverlap1);
   fChain->SetBranchAddress("recoDiele_pt2", &recoDiele_pt2, &b_recoDiele_pt2);
   fChain->SetBranchAddress("recoDiele_eta2", &recoDiele_eta2, &b_recoDiele_eta2);
   fChain->SetBranchAddress("recoDiele_phi2", &recoDiele_phi2, &b_recoDiele_phi2);
   fChain->SetBranchAddress("recoDiele_charge2", &recoDiele_charge2, &b_recoDiele_charge2);
   fChain->SetBranchAddress("recoDiele_mass2", &recoDiele_mass2, &b_recoDiele_mass2);
   fChain->SetBranchAddress("recoDiele_mvaValue2", &recoDiele_mvaValue2, &b_recoDiele_mvaValue2);
   fChain->SetBranchAddress("recoDiele_mvaPFValue2", &recoDiele_mvaPFValue2, &b_recoDiele_mvaPFValue2);
   fChain->SetBranchAddress("recoDiele_isPF2", &recoDiele_isPF2, &b_recoDiele_isPF2);
   fChain->SetBranchAddress("recoDiele_isLowPt2", &recoDiele_isLowPt2, &b_recoDiele_isLowPt2);
   fChain->SetBranchAddress("recoDiele_isPFoverlap2", &recoDiele_isPFoverlap2, &b_recoDiele_isPFoverlap2);
   fChain->SetBranchAddress("recoDiele_pt", &recoDiele_pt, &b_recoDiele_pt);
   fChain->SetBranchAddress("recoDiele_eta", &recoDiele_eta, &b_recoDiele_eta);
   fChain->SetBranchAddress("recoDiele_phi", &recoDiele_phi, &b_recoDiele_phi);
   fChain->SetBranchAddress("recoDiele_mass", &recoDiele_mass, &b_recoDiele_mass);
   fChain->SetBranchAddress("recoDiele_massErr", &recoDiele_massErr, &b_recoDiele_massErr);
   fChain->SetBranchAddress("recoDiele_pt1mode", &recoDiele_pt1mode, &b_recoDiele_pt1mode);
   fChain->SetBranchAddress("recoDiele_eta1mode", &recoDiele_eta1mode, &b_recoDiele_eta1mode);
   fChain->SetBranchAddress("recoDiele_phi1mode", &recoDiele_phi1mode, &b_recoDiele_phi1mode);
   fChain->SetBranchAddress("recoDiele_pt2mode", &recoDiele_pt2mode, &b_recoDiele_pt2mode);
   fChain->SetBranchAddress("recoDiele_eta2mode", &recoDiele_eta2mode, &b_recoDiele_eta2mode);
   fChain->SetBranchAddress("recoDiele_phi2mode", &recoDiele_phi2mode, &b_recoDiele_phi2mode);
   Notify();
}

Bool_t Analysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analysis_cxx
