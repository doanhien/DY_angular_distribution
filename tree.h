#ifndef tree_h
#define tree_h

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include <iostream>

float deltaPhi(float phi1, float phi2) {

  float dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
  if (dPhi <= -TMath::Pi()) dPhi += 2* TMath::Pi();

  return dPhi;

}

float deltaR(float eta1, float phi1, float eta2, float phi2) {

  float dEta = eta1 - eta2;
  float dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
  if (dPhi <= -TMath::Pi()) dPhi += 2* TMath::Pi();

  return sqrt(pow(dEta,2) + pow (dPhi,2));

}

const int maxPar = 50;

Int_t run_;
Int_t lumi_;
Long64_t event_;
Int_t nVtx_;
Float_t vx_;
Float_t vy_;
Float_t vz_;
Float_t genWeight_;
Float_t puweigj;
Float_t puweigj_69p2nb;
Float_t puweigj_65nb;
Float_t puweigj_63nb;
Float_t weight;
Float_t rho_;
Float_t pdfWeight_;
Float_t nnpdf_nlo_wei_[102];
Float_t nnpdf_nlo_scale_wei_[7];
Float_t PDF4LHC_wei_[103];
Float_t PDF4LHC_scale_wei_[6];

Int_t     trig_Ele27_WPTight;
Int_t     trig_Ele17_Ele12;
Int_t     trig_Ele23_Ele12;
Int_t     trig_Ele23_Ele12_nonDZ;
Int_t     trig_Ele23_Ele12_DZ;

Int_t     trig_Mu17_Mu8;
Int_t     trig_Mu17_Mu8_nonDZ;
Int_t     trig_Mu17_Mu8_DZ;


Int_t aMCNLO;
Int_t nPU;
Int_t leptType;
Float_t lept0_pt;
Float_t lept0_eta;
Float_t lept0_phi;
Float_t lept0_sceta;
Float_t lept0_miniRelIso;
Float_t lept0_trkIso;
Int_t lept0_pdgId;
Float_t lept0_normalizedChi2;
Int_t lept0_nValidMuonHits;
Int_t lept0_nMatchedStations;
Int_t lept0_nValidPixelHits;
Int_t lept0_nTrackerLayers;
Float_t lept0_muonBestTrack_dxyVTX;
Float_t lept0_muonBestTrack_dzVTX;
Float_t lept0_ptError;
Float_t lept0_sigmaIetaIeta;
Float_t lept0_dEtaIn;
Float_t lept0_dPhiIn;
Float_t lept0_hOverE;
Float_t lept0_ooEmooP;
Float_t lept0_d0;
Float_t lept0_dz;
Float_t lept0_mva;
Float_t lept0_chiso;
Float_t lept0_phoiso;
Float_t lept0_neuiso;
Float_t lept0_SIP;
Float_t lept0_sigEOverE;
Int_t lept0_expectedMissingInnerHits;
Float_t lept0_RecoSF;
Float_t lept0_SelSF;
Float_t lept0_trigSF;
Float_t lept0_RecoSF_Err;
Float_t lept0_SelSF_Err;
Float_t lept0_trigSF_Err;
Float_t lept0_EnErr;
Float_t lept0_En;

Float_t lept1_pt;
Float_t lept1_eta;
Float_t lept1_phi;
Float_t lept1_sceta;
Float_t lept1_miniRelIso;
Float_t lept1_trkIso;
Int_t lept1_pdgId;
Float_t lept1_normalizedChi2;
Int_t lept1_nValidMuonHits;
Int_t lept1_nMatchedStations;
Int_t lept1_nValidPixelHits;
Int_t lept1_nTrackerLayers;
Float_t lept1_muonBestTrack_dxyVTX;
Float_t lept1_muonBestTrack_dzVTX;
Float_t lept1_ptError;
Float_t lept1_sigmaIetaIeta;
Float_t lept1_dEtaIn;
Float_t lept1_dPhiIn;
Float_t lept1_hOverE;
Float_t lept1_ooEmooP;
Float_t lept1_d0;
Float_t lept1_dz;
Float_t lept1_mva;
Float_t lept1_chiso;
Float_t lept1_phoiso;
Float_t lept1_neuiso;
Float_t lept1_SIP;
Float_t lept1_sigEOverE;
Int_t lept1_expectedMissingInnerHits;
Float_t lept1_RecoSF;
Float_t lept1_SelSF;
Float_t lept1_trigSF;
Float_t lept1_RecoSF_Err;
Float_t lept1_SelSF_Err;
Float_t lept1_trigSF_Err;
Float_t lept1_EnErr;
Float_t lept1_En;
Float_t lept_dzSF;
Float_t lept_dzSF_Err;
Float_t pair_dPhi;

Float_t deltaPhi_lept;
Float_t deltaR_lept;

Float_t z_pt;
Float_t z_eta;
Float_t z_phi;
Float_t z_y;
Float_t z_mass;
Int_t   z_charge;
Float_t boss_pt;
Float_t boss_eta;
Float_t boss_phi;
Float_t boss_mass;
Float_t boss_massrel;
Float_t lep0_theta_Zrest;
Float_t lep0_costheta_Zrest;
Float_t lep0_phi_Zrest;
Float_t lep1_theta_Zrest;
Float_t lep1_costheta_Zrest;
Float_t lep1_phi_Zrest;
Float_t lep_theta_Zrest;
Float_t lep_costheta_Zrest;
Float_t lep_phi_Zrest;


Int_t nselJet;
Float_t jetPt_[40];
Float_t jetEta_[40];
Float_t jetJECUnc_[40];
Float_t jetP4SmearUp_[40];
Float_t jetP4SmearDo_[40];
Float_t jetP4Smear_[40];
Float_t dRLep1Jet[40];
Float_t dRLep2Jet[40];

//eg resolution and scale
Float_t eleScale_stat_up_1;
Float_t eleScale_stat_dn_1;
Float_t eleScale_syst_up_1;
Float_t eleScale_syst_dn_1;
Float_t eleScale_gain_up_1;
Float_t eleScale_gain_dn_1;
Float_t eleResol_rho_up_1;
Float_t eleResol_rho_dn_1;
Float_t eleResol_phi_up_1;
Float_t eleResol_phi_dn_1;

Float_t eleScale_stat_up_2;
Float_t eleScale_stat_dn_2;
Float_t eleScale_syst_up_2;
Float_t eleScale_syst_dn_2;
Float_t eleScale_gain_up_2;
Float_t eleScale_gain_dn_2;
Float_t eleResol_rho_up_2;
Float_t eleResol_rho_dn_2;
Float_t eleResol_phi_up_2;
Float_t eleResol_phi_dn_2;


Int_t ngenEle;
Int_t ngenMu;
Int_t ngenlep;
//Int_t genlepType;
Float_t   genMuPt[maxPar];
Float_t   genMuPhi[maxPar];
Float_t   genMuEta[maxPar];
Int_t     genMuCh[maxPar];
Float_t   genlepPt[maxPar];
Float_t   genlepPhi[maxPar];
Float_t   genlepEta[maxPar];
Int_t     genlepCh[maxPar];

Float_t   genZm;
Float_t   genZpt;
Float_t   genZy;
Int_t     genZch;
Int_t     genZtype;
Int_t     lep1fsr;
Int_t     lep2fsr;
Float_t   mcZm;

Float_t gen_lep0_theta_Zrest;
Float_t gen_lep0_costheta_Zrest;
Float_t gen_lep0_phi_Zrest;
Float_t gen_lep1_theta_Zrest;
Float_t gen_lep1_costheta_Zrest;
Float_t gen_lep1_phi_Zrest;
Float_t gen_lep_theta_Zrest;
Float_t gen_lep_costheta_Zrest;
Float_t gen_lep_phi_Zrest;


Float_t   lheEleEta1;
Float_t   lheEleEta2;
Float_t   lheMuEta1;
Float_t   lheMuEta2;

Float_t   lhe_Zm;
Int_t     lhe_type;

Int_t ngenjet_;
Float_t genJetPt[60];
Float_t genJetEta[60];

Int_t nlheJet;
Float_t lheJetPt[3];
Float_t lheJetEta[3];

void inittree(TTree* tree) {

  tree->Branch("run",         &run_);
  tree->Branch("lumi",        &lumi_);
  tree->Branch("event",       &event_);
  tree->Branch("nVtx",        &nVtx_);
  tree->Branch("vx",          &vx_);
  tree->Branch("vy",          &vy_);
  tree->Branch("vz",          &vz_);
  tree->Branch("genWeight",   &genWeight_);
  tree->Branch("weight",      &weight);
  tree->Branch("puweigj",     &puweigj);
  tree->Branch("puweigj_69p2nb",     &puweigj_69p2nb);
  tree->Branch("puweigj_65nb",       &puweigj_65nb);
  tree->Branch("puweigj_63nb",       &puweigj_63nb);
  tree->Branch("rho", &rho_);

  //tree->Branch("trig_Ele27_WPTight",      &trig_Ele27_WPTight);
  //tree->Branch("trig_Ele17_Ele12",        &trig_Ele17_Ele12);
  tree->Branch("trig_Ele23_Ele12",        &trig_Ele23_Ele12);
  tree->Branch("trig_Ele23_Ele12_nonDZ",  &trig_Ele23_Ele12_nonDZ);
  tree->Branch("trig_Ele23_Ele12_DZ",     &trig_Ele23_Ele12_DZ);
  tree->Branch("trig_Mu17_Mu8",           &trig_Mu17_Mu8);
  tree->Branch("trig_Mu17_Mu8_nonDZ",  &trig_Mu17_Mu8_nonDZ);
  tree->Branch("trig_Mu17_Mu8_DZ",     &trig_Mu17_Mu8_DZ);

  tree->Branch("leptType",                &leptType);
  tree->Branch("lept0_pt",                &lept0_pt);
  tree->Branch("lept0_eta",               &lept0_eta);
  tree->Branch("lept0_phi",               &lept0_phi);
  tree->Branch("lept0_sceta",             &lept0_sceta);

  tree->Branch("lept0_RecoSF",               &lept0_RecoSF);
  tree->Branch("lept0_SelSF",                &lept0_SelSF);
  tree->Branch("lept0_EnErr",                &lept0_EnErr);
  tree->Branch("lept0_En",                   &lept0_En);
  tree->Branch("lept0_trigSF",               &lept0_trigSF);
  tree->Branch("lept0_RecoSF_Err",           &lept0_RecoSF_Err);
  tree->Branch("lept0_SelSF_Err",            &lept0_SelSF_Err);
  tree->Branch("lept0_trigSF_Err",           &lept0_trigSF_Err);

  tree->Branch("lept1_pt",                &lept1_pt);
  tree->Branch("lept1_eta",               &lept1_eta);
  tree->Branch("lept1_phi",               &lept1_phi);
  tree->Branch("lept1_sceta",             &lept1_sceta);
  /*tree->Branch("lept1_miniRelIso",        &lept1_miniRelIso);
  tree->Branch("lept1_trkIso",            &lept1_trkIso);
  tree->Branch("lept1_pdgId",             &lept1_pdgId);
  tree->Branch("lept1_normalizedChi2",    &lept1_normalizedChi2);
  tree->Branch("lept1_nValidMuonHits",    &lept1_nValidMuonHits);
  tree->Branch("lept1_nMatchedStations",  &lept1_nMatchedStations);
  tree->Branch("lept1_nValidPixelHits",   &lept1_nValidPixelHits);
  tree->Branch("lept1_nTrackerLayers",    &lept1_nTrackerLayers);
  tree->Branch("lept1_muonBestTrack_dxyVTX", &lept1_muonBestTrack_dxyVTX);
  tree->Branch("lept1_muonBestTrack_dzVTX",  &lept1_muonBestTrack_dzVTX);
  tree->Branch("lept1_ptError",              &lept1_ptError);
  tree->Branch("lept1_sigmaIetaIeta",        &lept1_sigmaIetaIeta);
  tree->Branch("lept1_dEtaIn",               &lept1_dEtaIn);
  tree->Branch("lept1_dPhiIn",               &lept1_dPhiIn);
  tree->Branch("lept1_hOverE",               &lept1_hOverE);
  tree->Branch("lept1_ooEmooP",              &lept1_ooEmooP);
  tree->Branch("lept1_d0",                   &lept1_d0);
  tree->Branch("lept1_dz",                   &lept1_dz);
  tree->Branch("lept1_mva",                  &lept1_mva);
  tree->Branch("lept1_chiso",                &lept1_chiso);
  tree->Branch("lept1_phoiso",               &lept1_phoiso);
  tree->Branch("lept1_neuiso",               &lept1_neuiso);
  tree->Branch("lept1_SIP",                  &lept1_SIP);
  tree->Branch("lept1_sigEOverE",            &lept1_sigEOverE);
  tree->Branch("lept1_expectedMissingInnerHits", &lept1_expectedMissingInnerHits);*/
  tree->Branch("lept1_RecoSF",               &lept1_RecoSF);
  tree->Branch("lept1_SelSF",                &lept1_SelSF);
  tree->Branch("lept1_EnErr",                &lept1_EnErr);
  tree->Branch("lept1_En",                   &lept1_En);
  tree->Branch("lept1_trigSF",               &lept1_trigSF);
  tree->Branch("lept1_RecoSF_Err",           &lept1_RecoSF_Err);
  tree->Branch("lept1_SelSF_Err",            &lept1_SelSF_Err);
  tree->Branch("lept1_trigSF_Err",           &lept1_trigSF_Err);
  tree->Branch("lept_dzSF",                  &lept_dzSF);
  tree->Branch("lept_dzSF_Err",              &lept_dzSF_Err);
  tree->Branch("pair_dPhi",                  &pair_dPhi);

  tree->Branch("deltaR_lept",    &deltaR_lept);
  tree->Branch("deltaPhi_lept",  &deltaPhi_lept);

  tree->Branch("z_pt",                &z_pt);
  tree->Branch("z_eta",               &z_eta);
  //tree->Branch("z_phi",               &z_phi);
  tree->Branch("z_y",                 &z_y);
  tree->Branch("boss_pt",             &boss_pt);
  tree->Branch("boss_eta",            &boss_eta);
  //tree->Branch("boss_phi",            &boss_phi);
  //tree->Branch("boss_massrel",        &boss_massrel);
  tree->Branch("z_mass",              &z_mass);
  tree->Branch("z_charge",            &z_charge);
  tree->Branch("boss_mass",           &boss_mass);
  tree->Branch("lep0_theta_Zrest",    &lep0_theta_Zrest);
  tree->Branch("lep0_costheta_Zrest", &lep0_costheta_Zrest);
  tree->Branch("lep0_phi_Zrest",      &lep0_phi_Zrest);
  tree->Branch("lep1_theta_Zrest",    &lep1_theta_Zrest);
  tree->Branch("lep1_costheta_Zrest", &lep1_costheta_Zrest);
  tree->Branch("lep1_phi_Zrest",      &lep1_phi_Zrest);
  tree->Branch("lep_theta_Zrest",     &lep_theta_Zrest);
  tree->Branch("lep_costheta_Zrest",  &lep_costheta_Zrest);
  tree->Branch("lep_phi_Zrest",       &lep_phi_Zrest);
  

  tree->Branch("nselJet",             &nselJet);
  tree->Branch("jetPt",               &jetPt_,   "jetPt[nselJet]/F");
  tree->Branch("jetEta",              &jetEta_,  "jetEta[nselJet]/F");
  tree->Branch("jetJECUnc",           &jetJECUnc_,     "jetJECUnc[nselJet]/F");
  tree->Branch("jetP4SmearUp",        &jetP4SmearUp_,  "jetP4SmearUp[nselJet]/F");
  tree->Branch("jetP4SmearDo",        &jetP4SmearDo_,  "jetP4SmearDo[nselJet]/F");
  tree->Branch("jetP4Smear",          &jetP4Smear_,    "jetP4Smear[nselJet]/F");
  tree->Branch("dRLep1Jet",           &dRLep1Jet,      "dRLep1Jet[nselJet]/F");
  tree->Branch("dRLep2Jet",           &dRLep2Jet,      "dRLep2Jet[nselJet]/F");

  tree->Branch("eleScale_stat_up_1",    &eleScale_stat_up_1);
  tree->Branch("eleScale_stat_dn_1",    &eleScale_stat_dn_1);
  tree->Branch("eleScale_syst_up_1",    &eleScale_syst_up_1);
  tree->Branch("eleScale_syst_dn_1",    &eleScale_syst_dn_1);
  tree->Branch("eleScale_gain_up_1",    &eleScale_gain_up_1);
  tree->Branch("eleScale_gain_dn_1",    &eleScale_gain_dn_1);
  tree->Branch("eleResol_rho_up_1",     &eleResol_rho_up_1);
  tree->Branch("eleResol_rho_dn_1",     &eleResol_rho_dn_1);
  tree->Branch("eleResol_phi_up_1",     &eleResol_phi_up_1);
  tree->Branch("eleResol_phi_dn_1",     &eleResol_phi_dn_1);

  tree->Branch("eleScale_stat_up_2",    &eleScale_stat_up_2);
  tree->Branch("eleScale_stat_dn_2",    &eleScale_stat_dn_2);
  tree->Branch("eleScale_syst_up_2",    &eleScale_syst_up_2);
  tree->Branch("eleScale_syst_dn_2",    &eleScale_syst_dn_2);
  tree->Branch("eleScale_gain_up_2",    &eleScale_gain_up_2);
  tree->Branch("eleScale_gain_dn_2",    &eleScale_gain_dn_2);
  tree->Branch("eleResol_rho_up_2",     &eleResol_rho_up_2);
  tree->Branch("eleResol_rho_dn_2",     &eleResol_rho_dn_2);
  tree->Branch("eleResol_phi_up_2",     &eleResol_phi_up_2);
  tree->Branch("eleResol_phi_dn_2",     &eleResol_phi_dn_2);


  tree->Branch("ngenEle",     &ngenEle);
  tree->Branch("ngenMu",      &ngenMu);
  tree->Branch("ngenlep",      &ngenlep);
  tree->Branch("genlepPt",    &genlepPt,   "genlepPt[ngenlep]/F");
  tree->Branch("genlepPhi",   &genlepPhi,  "genlepPhi[ngenlep]/F");
  tree->Branch("genlepEta",   &genlepEta,  "genlepEta[ngenlep]/F");
  tree->Branch("genlepCh",    &genlepCh,   "genlepCh[ngenlep]/I");

  tree->Branch("genZm",       &genZm);
  tree->Branch("genZpt",      &genZpt);
  tree->Branch("genZy",       &genZy);
  tree->Branch("genZch",      &genZch);
  tree->Branch("genZtype",    &genZtype);

  tree->Branch("gen_lep0_theta_Zrest",    &gen_lep0_theta_Zrest);
  tree->Branch("gen_lep0_costheta_Zrest", &gen_lep0_costheta_Zrest);
  tree->Branch("gen_lep0_phi_Zrest",      &gen_lep0_phi_Zrest);
  tree->Branch("gen_lep1_theta_Zrest",    &gen_lep1_theta_Zrest);
  tree->Branch("gen_lep1_costheta_Zrest", &gen_lep1_costheta_Zrest);
  tree->Branch("gen_lep1_phi_Zrest",      &gen_lep1_phi_Zrest);
  tree->Branch("gen_lep_theta_Zrest",     &gen_lep_theta_Zrest);
  tree->Branch("gen_lep_costheta_Zrest",  &gen_lep_costheta_Zrest);
  tree->Branch("gen_lep_phi_Zrest",       &gen_lep_phi_Zrest);

  /*
  tree->Branch("lheEle1",     &lheEle1);
  tree->Branch("lheEle2",     &lheEle2);
  tree->Branch("lheMu1",      &lheMu1);
  tree->Branch("lheMu2",      &lheMu2);

  tree->Branch("lheEleEta1",    &lheEleEta1);
  tree->Branch("lheEleEta2",    &lheEleEta2);
  tree->Branch("lheMuEta1",     &lheMuEta1);
  tree->Branch("lheMuEta2",     &lheMuEta2);

  tree->Branch("lhe_Zm",        &lhe_Zm);
  tree->Branch("lhe_type",      &lhe_type);

  tree->Branch("ngenjet",      &ngenjet_);
  tree->Branch("genJetPt",     &genJetPt,    "genJetPt[ngenjet]/F");
  tree->Branch("genJetEta",    &genJetEta,   "genJetEta[ngenjet]/F");
  
  tree->Branch("nlheJet",      &nlheJet);
  tree->Branch("lheJetPt",     &lheJetPt,    "lheJetPt[nlheJet]/F");
  tree->Branch("lheJetEta",    &lheJetEta,   "lheJetEta[nlheJet]/F");
  */

}

void inittreeGen(TTree *tree){

  tree->Branch("run",         &run_);
  tree->Branch("lumi",        &lumi_);
  tree->Branch("event",       &event_);
  tree->Branch("genWeight",   &genWeight_);
  tree->Branch("puweigj",     &puweigj);
  tree->Branch("puweigj_69p2nb",     &puweigj_69p2nb);
  tree->Branch("puweigj_65nb",       &puweigj_65nb);
  tree->Branch("puweigj_63nb",       &puweigj_63nb);
  tree->Branch("pdfWeight",          &pdfWeight_);
  tree->Branch("nnpdf_nlo_wei",           &nnpdf_nlo_wei_,            "nnpdf_nlo_wei[102]/F");
  tree->Branch("nnpdf_nlo_scale_wei",     &nnpdf_nlo_scale_wei_,      "nnpdf_nlo_scale_wei[7]/F");
  tree->Branch("PDF4LHC_wei",             &PDF4LHC_wei_,              "PDF4LHC_wei[103]/F");
  tree->Branch("PDF4LHC_scale_wei",       &PDF4LHC_scale_wei_,        "PDF4LHC_scale_wei[6]/F");

  tree->Branch("ngenEle",     &ngenEle);
  tree->Branch("ngenMu",      &ngenMu);
  tree->Branch("ngenlep",      &ngenlep);
  //tree->Branch("genlepType",   &genlepType);
  //tree->Branch("genMuPt",     &genMuPt,    "genMuPt[ngenMu]/F");
  //tree->Branch("genMuPhi",    &genMuPhi,   "genMuPhi[ngenMu]/F");
  //tree->Branch("genMuEta",    &genMuEta,   "genMuEta[ngenMu]/F");
  //tree->Branch("genMuCh",     &genMuCh,    "genMuCh[ngenMu]/I");
  tree->Branch("genlepPt",    &genlepPt,   "genlepPt[ngenlep]/F");
  tree->Branch("genlepPhi",   &genlepPhi,  "genlepPhi[ngenlep]/F");
  tree->Branch("genlepEta",   &genlepEta,  "genlepEta[ngenlep]/F");
  tree->Branch("genlepCh",    &genlepCh,   "genlepCh[ngenlep]/I");
  tree->Branch("genZm",       &genZm);
  tree->Branch("genZpt",      &genZpt);
  tree->Branch("genZy",       &genZy);
  tree->Branch("genZtype",    &genZtype);
  tree->Branch("lep1fsr",     &lep1fsr);
  tree->Branch("lep2fsr",     &lep2fsr);
  tree->Branch("mcZm",        &mcZm);

  tree->Branch("gen_lep0_theta_Zrest",    &gen_lep0_theta_Zrest);
  tree->Branch("gen_lep0_costheta_Zrest", &gen_lep0_costheta_Zrest);
  tree->Branch("gen_lep0_phi_Zrest",      &gen_lep0_phi_Zrest);
  tree->Branch("gen_lep1_theta_Zrest",    &gen_lep1_theta_Zrest);
  tree->Branch("gen_lep1_costheta_Zrest", &gen_lep1_costheta_Zrest);
  tree->Branch("gen_lep1_phi_Zrest",      &gen_lep1_phi_Zrest);
  tree->Branch("gen_lep_theta_Zrest",     &gen_lep_theta_Zrest);
  tree->Branch("gen_lep_costheta_Zrest",  &gen_lep_costheta_Zrest);
  tree->Branch("gen_lep_phi_Zrest",       &gen_lep_phi_Zrest);

  /*
  tree->Branch("lheEle1",     &lheEle1);
  tree->Branch("lheEle2",     &lheEle2);
  tree->Branch("lheMu1",      &lheMu1);
  tree->Branch("lheMu2",      &lheMu2);

  tree->Branch("lheEleEta1",    &lheEleEta1);
  tree->Branch("lheEleEta2",    &lheEleEta2);
  tree->Branch("lheMuEta1",     &lheMuEta1);
  tree->Branch("lheMuEta2",     &lheMuEta2);

  tree->Branch("lhe_Zm",        &lhe_Zm);
  tree->Branch("lhe_type",      &lhe_type);
  */
  tree->Branch("ngenjet",      &ngenjet_);
  tree->Branch("genJetPt",     &genJetPt,    "genJetPt[ngenjet]/F");
  tree->Branch("genJetEta",    &genJetEta,   "genJetEta[ngenjet]/F");

  //tree->Branch("nlheJet",      &nlheJet);
  //tree->Branch("lheJetPt",     &lheJetPt,    "lheJetPt[nlheJet]/F");
  //tree->Branch("lheJetEta",    &lheJetEta,   "lheJetEta[nlheJet]/F");

  

}

Int_t     run;
Long64_t  event;
Int_t     lumis;
Bool_t    isData;
ULong64_t hlt;
Int_t     nVtx;
Float_t   vz;
Float_t   vx; 
Float_t   vy;
Float_t   genWeight;
Float_t   rho;
Bool_t    isPVGood;
Float_t   pdfWeight;
Float_t*  pdfSystWeight;

Int_t nEle ;
Int_t* eleCharge;
Int_t* eleChargeConsistent;
Float_t *eleEn;
Float_t* elePt ;
Float_t* eleEta ;
Float_t* elePhi ;
Float_t* eleHoverE ;
Float_t* eleEoverP ;
Float_t* eleEoverPInv ;
Float_t* eleSigmaIEtaIEta_Full5x5 ;
Float_t* eleSigmaIPhiIPhi;
Int_t* eleMissHits ;
Float_t* eleD0 ;
Float_t* eleDz ;
Float_t* eledEtaAtVtx ;
Float_t* eledPhiAtVtx ;
Float_t* eledEtaseedAtVtx ;
Int_t* eleConvVeto;
Int_t* eleEcalDrivenSeed;
Float_t* eleE1x5Full5x5;
Float_t* eleE2x5Full5x5;
Float_t* eleE5x5Full5x5;
Float_t* eleSCEta ;
Float_t* eleSCPhi;
Float_t* eleSCEtaWidth;
Float_t* eleSCPhiWidth;
Float_t* eleSCEn ;
Float_t* elePFChIso ;
Float_t* elePFPhoIso ;
Float_t* elePFNeuIso ;
Float_t* elePFPUIso ;
Float_t* elePFMiniIso ;
Short_t* eleID;
Float_t* eleIDMVA;
Float_t* eleIDMVAHZZ;
Float_t* elePFClusEcalIso;
Float_t* elePFClusHcalIso;
Float_t* eleDr03TkSumPt;
Float_t* eleGSFChi2;
Float_t* eleEcalEnErr;
Float_t* eleSIP;


Int_t    nMu;
Int_t*   muType;
Float_t* muPt;
Float_t* muEta;
Float_t* muPhi;
Int_t*   muCh;
Float_t*  muSIP;
Float_t* muChi2NDF;
Int_t*   muMuonHits;
Int_t*   muStations;
Int_t*   muTrkLayers;
Int_t*   muPixelHits;
Float_t* muInnerD0;
Float_t* muInnerDz;
Float_t* muD0;
Float_t* muDz;
Float_t* muBestTrkPtError;
Float_t* muBestTrkPt;
Int_t*   muBestTrkType;
Float_t* muPFChIso;
Float_t* muPFPhoIso;
Float_t* muPFNeuIso;
Float_t* muPFPUIso;
Float_t* muPFChIso03;
Float_t* muPFPhoIso03;
Float_t* muPFNeuIso03;
Float_t* muPFPUIso03;
Float_t* muPFMiniIso;
Float_t* muIsoTrk;
Int_t* muIDbit;
//Short_t* muIDbit;


//jet
Int_t nJet;
Float_t* jetPt;
Float_t* jetEta;
Float_t* jetPhi;
Float_t* jetEn;
Float_t* jetJECUnc;
Float_t* jetP4SmearUp;
Float_t* jetP4SmearDo;
Float_t* jetP4Smear;
Float_t* jetNHF;
Float_t* jetCHF;
Float_t* jetNEF;
Float_t* jetCEF;
Int_t* jetNCH;
Int_t* jetNNP;

//vector<bool> jetPFLooseId;

//gen level
Int_t nMC;
Int_t* mcPID;
Int_t* mcMomPID;
Int_t* mcGMomPID;
Float_t* mcPt;
Float_t* mcEt;
Float_t* mcPhi;
Float_t* mcEta;
Float_t* mcMomPt;
Float_t* mcMomEta;
Float_t* mcMomPhi;
Float_t* mcMomMass;
UShort_t* mcStatusFlag;

Float_t genEle1;
Float_t genEle2;
Float_t genMu1;
Float_t genMu2;
Float_t genEleEta1;
Float_t genEleEta2;
Float_t genMuEta1;
Float_t genMuEta2;
Float_t genElePhi1;
Float_t genElePhi2;
Float_t genMuPhi1;
Float_t genMuPhi2;

Int_t ngenjet;
Float_t* genjetPt_all;
Float_t* genjetEta_all;
Float_t* genjetPhi_all;
Float_t* genjetEn_all;

Int_t nLHEJet;
Float_t* lheJetPt_;
Float_t* lheJetEta_;
Float_t* lheJetPhi_;


//egamma resolution and scale
Float_t* eleScale_stat_up;
Float_t* eleScale_stat_dn;
Float_t* eleScale_syst_up;
Float_t* eleScale_syst_dn;
Float_t* eleScale_gain_up;
Float_t* eleScale_gain_dn;
Float_t* eleResol_rho_up;
Float_t* eleResol_rho_dn;
Float_t* eleResol_phi_up;
Float_t* eleResol_phi_dn;


void readggtree(TreeReader &data) {

  run = data.GetInt("run");
  event = data.GetLong64("event");
  lumis = data.GetInt("lumis");
  hlt = data.GetLong64("HLTEleMuX");
  nVtx = data.GetInt("nVtx");
  vx = data.GetFloat("vtx");
  vy = data.GetFloat("vty");
  vz = data.GetFloat("vtz");
  rho = data.GetFloat("rho");
  isPVGood = data.GetBool("isPVGood");
  //pdfWeight = data.GetFloat("pdfWeight");
  //pdfSystWeight = data.GetPtrFloat("pdfSystWeight");

  nEle = data.GetInt("nEle");
  eleCharge = data.GetPtrInt("eleCharge");
  elePt = data.GetPtrFloat("eleCalibPt");
  //elePt = data.GetPtrFloat("elePt");
  eleEn = data.GetPtrFloat("eleCalibEn");
  eleEta = data.GetPtrFloat("eleEta");
  elePhi = data.GetPtrFloat("elePhi");
  eleHoverE = data.GetPtrFloat("eleHoverE");
  eleEoverP = data.GetPtrFloat("eleEoverP");
  eleEoverPInv = data.GetPtrFloat("eleEoverPInv");
  eleSigmaIEtaIEta_Full5x5 = data.GetPtrFloat("eleSigmaIEtaIEtaFull5x5");
  eleMissHits = data.GetPtrInt("eleMissHits");
  eleD0 = data.GetPtrFloat("eleD0");
  eleDz = data.GetPtrFloat("eleDz");
  eledEtaAtVtx = data.GetPtrFloat("eledEtaAtVtx");
  eledPhiAtVtx = data.GetPtrFloat("eledPhiAtVtx");
  //eledEtaseedAtVtx = data.GetPtrFloat("eledEtaseedAtVtx");
  eleConvVeto = data.GetPtrInt("eleConvVeto");
  eleEcalDrivenSeed = data.GetPtrInt("eleEcalDrivenSeed");
  //eleE1x5Full5x5 = data.GetPtrFloat("eleE1x5Full5x5");
  //eleE2x5Full5x5 = data.GetPtrFloat("eleE2x5Full5x5");
  //eleE5x5Full5x5 = data.GetPtrFloat("eleE5x5Full5x5");
  eleSCEta = data.GetPtrFloat("eleSCEta");
  eleSCPhi = data.GetPtrFloat("eleSCPhi");
  eleSCEtaWidth = data.GetPtrFloat("eleSCEtaWidth");
  eleSCPhiWidth = data.GetPtrFloat("eleSCPhiWidth");
  eleSCEn = data.GetPtrFloat("eleSCEn");
  elePFChIso = data.GetPtrFloat("elePFChIso");
  elePFPhoIso = data.GetPtrFloat("elePFPhoIso");
  elePFNeuIso = data.GetPtrFloat("elePFNeuIso");
  //elePFMiniIso = data.GetPtrFloat("elePFMiniIso");
  eleID = data.GetPtrShort("eleIDbit");
  //eleIDMVA = data.GetPtrFloat("eleIDMVA");
  //eleIDMVAHZZ = data.GetPtrFloat("eleIDMVAHZZ");
  //eleIDMVANonTrg = data.GetPtrFloat("eleIDMVANonTrg");
  elePFClusEcalIso = data.GetPtrFloat("elePFClusEcalIso");
  elePFClusHcalIso = data.GetPtrFloat("elePFClusHcalIso");
  //eleDr03TkSumPt = data.GetPtrFloat("eleDr03TkSumPt");
  //eleGSFChi2 = data.GetPtrFloat("eleGSFChi2");
  //eleEcalEnErr = data.GetPtrFloat("eleEcalEnErr");
  eleSIP = data.GetPtrFloat("eleSIP");

  nMu = data.GetInt("nMu");
  muType = data.GetPtrInt("muType");
  //muPt = data.GetPtrFloat("muCalibPt");
  //muPt = data.GetPtrFloat("muCalibPtSyst2");
  muPt = data.GetPtrFloat("muPt");
  muEta = data.GetPtrFloat("muEta");
  muPhi = data.GetPtrFloat("muPhi");
  muCh = data.GetPtrInt("muCharge");
  muSIP = data.GetPtrFloat("muSIP");
  muChi2NDF = data.GetPtrFloat("muChi2NDF");
  muMuonHits = data.GetPtrInt("muMuonHits");
  muStations = data.GetPtrInt("muStations");
  muTrkLayers = data.GetPtrInt("muTrkLayers");
  muPixelHits = data.GetPtrInt("muPixelHits");
  muInnerD0 = data.GetPtrFloat("muInnerD0");
  muInnerDz = data.GetPtrFloat("muInnerDz");
  muD0 = data.GetPtrFloat("muD0");
  muDz = data.GetPtrFloat("muDz");
  muBestTrkPtError = data.GetPtrFloat("muBestTrkPtError");
  muBestTrkPt = data.GetPtrFloat("muBestTrkPt");
  muBestTrkType = data.GetPtrInt("muBestTrkType");
  muPFChIso = data.GetPtrFloat("muPFChIso");
  muPFNeuIso = data.GetPtrFloat("muPFNeuIso");
  muPFPhoIso = data.GetPtrFloat("muPFPhoIso");
  muPFPUIso = data.GetPtrFloat("muPFPUIso");
  muPFChIso03 = data.GetPtrFloat("muPFChIso03");
  muPFNeuIso03 = data.GetPtrFloat("muPFNeuIso03");
  muPFPhoIso03 = data.GetPtrFloat("muPFPhoIso03");
  muPFPUIso03 = data.GetPtrFloat("muPFPUIso03");
  //muPFMiniIso = data.GetPtrFloat("muPFMiniIso");
  muIsoTrk = data.GetPtrFloat("muIsoTrk");
  muIDbit = data.GetPtrInt("muIDbit");
  //muIDbit = data.GetPtrShort("muIDbit");


  nJet = data.GetInt("nJet");
  jetPt = data.GetPtrFloat("jetPt");
  jetEta = data.GetPtrFloat("jetEta");
  jetPhi = data.GetPtrFloat("jetPhi");
  jetEn  = data.GetPtrFloat("jetEn");
  if (data.HasMC()) {
    jetP4Smear = data.GetPtrFloat("jetP4Smear");
    jetP4SmearUp = data.GetPtrFloat("jetP4SmearUp");
    jetP4SmearDo = data.GetPtrFloat("jetP4SmearDo");
    jetJECUnc = data.GetPtrFloat("jetJECUnc");
  }
  jetNHF = data.GetPtrFloat("jetNHF");
  jetCHF = data.GetPtrFloat("jetCHF");
  jetNEF = data.GetPtrFloat("jetNEF");
  jetCEF = data.GetPtrFloat("jetCEF");
  jetNCH = data.GetPtrInt("jetNCH");
  jetNNP = data.GetPtrInt("jetNNP");

  //vector<bool> jetPFLooseId = ((vector<bool>) data.GetPtr("jetPFLooseId"));

  //get eg resolution scale   

  eleScale_stat_up = data.GetPtrFloat("eleScale_stat_up");
  eleScale_stat_dn = data.GetPtrFloat("eleScale_stat_dn");
  eleScale_syst_up = data.GetPtrFloat("eleScale_syst_up");
  eleScale_syst_dn = data.GetPtrFloat("eleScale_syst_dn");
  eleScale_gain_up = data.GetPtrFloat("eleScale_gain_up");
  eleScale_gain_dn = data.GetPtrFloat("eleScale_gain_dn");
  eleResol_rho_up = data.GetPtrFloat("eleResol_rho_up");
  eleResol_rho_dn = data.GetPtrFloat("eleResol_rho_dn");
  eleResol_phi_up = data.GetPtrFloat("eleResol_phi_up");
  eleResol_phi_dn = data.GetPtrFloat("eleResol_phi_dn");


  if (data.HasMC()) {
  nMC = data.GetInt("nMC");
  mcPID = data.GetPtrInt("mcPID");
  mcMomPID = data.GetPtrInt("mcMomPID");
  mcGMomPID = data.GetPtrInt("mcGMomPID");
  mcPt = data.GetPtrFloat("mcPt");
  mcEt = data.GetPtrFloat("mcEt");
  mcPhi = data.GetPtrFloat("mcPhi");
  mcEta = data.GetPtrFloat("mcEta");
  mcMomPt = data.GetPtrFloat("mcMomPt");
  mcMomEta = data.GetPtrFloat("mcMomEta");
  mcMomPhi = data.GetPtrFloat("mcMomPhi");
  mcMomMass = data.GetPtrFloat("mcMomMass");
  mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");
  genWeight = data.GetFloat("genWeight");

  /*
  ngenjet = data.GetInt("ngenjet");
  genjetPt_all = data.GetPtrFloat("genjetPt_all");
  genjetEta_all = data.GetPtrFloat("genjetEta_all");
  genjetPhi_all = data.GetPtrFloat("genjetPhi_all");
  genjetEn_all = data.GetPtrFloat("genjetEn_all");
  */
  /*
  genPho1 = data.GetFloat("lhePho1");
  genEle1 = data.GetFloat("lheEle1");
  genEle2 = data.GetFloat("lheEle2");
  genMu1  = data.GetFloat("lheMu1");
  genMu2  = data.GetFloat("lheMu2");
  
  genPhoEta1 = data.GetFloat("lhePhoEta1");
  genEleEta1 = data.GetFloat("lheEleEta1");
  genEleEta2 = data.GetFloat("lheEleEta2");
  genMuEta1  = data.GetFloat("lheMuEta1");
  genMuEta2  = data.GetFloat("lheMuEta2");

  genPhoPhi1 = data.GetFloat("lhePhoPhi1");
  genElePhi1 = data.GetFloat("lheElePhi1");
  genElePhi2 = data.GetFloat("lheElePhi2");
  genMuPhi1  = data.GetFloat("lheMuPhi1");
  genMuPhi2  = data.GetFloat("lheMuPhi2");

  nLHEJet = data.GetInt("nLHEJet");
  lheJetPt_ = data.GetPtrFloat("lheJetPt");
  lheJetEta_ = data.GetPtrFloat("lheJetEta");
  lheJetPhi_ = data.GetPtrFloat("lheJetPhi");
  */

  }

}


#endif
