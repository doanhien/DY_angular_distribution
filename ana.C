#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH1.h"
#include "TCanvas.h"

#include "tree.h"

#include <iostream>

using namespace std;

void ana(bool data = 1, bool ele =1) {

  TString fname;
  if (!data) fname = "minitrees/ZJets_aMCatNLO_correctZaxis.root";
  else {
    if (ele) fname = "minitrees/DoubleEG_NoRochesterCorr.root";
    else fname = "minitrees/DoubleMu_RochesterCorr_correctZaxis.root";
  }

  TFile *infile = new TFile(fname, "read");
  TTree *intree = (TTree*) infile->Get("outtree");

  TTree *ingentree;
  if (!data) ingentree = (TTree*) infile->Get("outtreeGen");

  if(!intree) {
    cout<<"could not read input tree\n";
  }



  Float_t puweigj_69p2nb;
  Float_t genWeight;
  
  Int_t     trig_Mu17_Mu8;
  Int_t     trig_Ele23_Ele12;
  
  Float_t lept0_pt;
  Float_t lept0_eta;
  Float_t lept0_phi;
  Float_t lept0_sceta;
  
  Float_t lept1_pt;
  Float_t lept1_eta;
  Float_t lept1_phi;
  Float_t lept1_sceta;
  
  Float_t z_pt;
  Float_t z_eta;
  Float_t z_phi;
  Float_t z_y;
  Float_t z_mass;
  Int_t   z_charge;
  Int_t   leptType;
  Int_t   lept0_pdgId;
  Int_t   lept1_pdgId;
  
  Float_t lept0_RecoSF;
  Float_t lept0_SelSF;
  Float_t lept0_trigSF;
  
  Float_t lept1_RecoSF;
  Float_t lept1_SelSF;
  Float_t lept1_trigSF;

  Float_t lep0_theta_Zrest;
  Float_t lep0_costheta_Zrest;
  Float_t lep0_phi_Zrest;
  Float_t lep1_theta_Zrest;
  Float_t lep1_costheta_Zrest;
  Float_t lep1_phi_Zrest;
  Float_t lep_theta_Zrest;
  Float_t lep_costheta_Zrest;
  Float_t lep_phi_Zrest;
  
  Int_t ngenEle;
  Int_t ngenMu;
  
  Float_t   genZm;
  Float_t   genZpt;
  Float_t   genZy;
  Int_t     genZch;
  Int_t     genZtype;
  
  Float_t gen_lep_theta_Zrest;
  Float_t gen_lep_costheta_Zrest;
  Float_t gen_lep_phi_Zrest;
  
  Float_t gen_lep0_theta_Zrest;
  Float_t gen_lep0_costheta_Zrest;
  Float_t gen_lep0_phi_Zrest;
  
  Float_t gen_lep1_theta_Zrest;
  Float_t gen_lep1_costheta_Zrest;
  Float_t gen_lep1_phi_Zrest;
  
  //const Int_t maxPar = 5;
  Float_t   genlepPt[5];
  Float_t   genlepPhi[5];
  Float_t   genlepEta[5];
  Int_t     genlepCh[5];

  intree->ResetBranchAddresses();
  intree->SetBranchAddress("genWeight",   &genWeight);
  intree->SetBranchAddress("puweigj_69p2nb",     &puweigj_69p2nb);
  intree->SetBranchAddress("trig_Ele23_Ele12",        &trig_Ele23_Ele12);
  intree->SetBranchAddress("trig_Mu17_Mu8",           &trig_Mu17_Mu8);

  intree->SetBranchAddress("leptType",                &leptType);
  intree->SetBranchAddress("lept0_pdgId",                &lept0_pdgId);
  intree->SetBranchAddress("lept1_pdgId",                &lept1_pdgId);
  intree->SetBranchAddress("lept0_pt",                &lept0_pt);
  intree->SetBranchAddress("lept0_eta",               &lept0_eta);
  intree->SetBranchAddress("lept0_phi",               &lept0_phi);
  intree->SetBranchAddress("lept0_sceta",             &lept0_sceta);

  intree->SetBranchAddress("lept0_RecoSF",               &lept0_RecoSF);
  intree->SetBranchAddress("lept0_SelSF",                &lept0_SelSF);
  intree->SetBranchAddress("lept0_trigSF",               &lept0_trigSF);

  intree->SetBranchAddress("lept1_pt",                &lept1_pt);
  intree->SetBranchAddress("lept1_eta",               &lept1_eta);
  intree->SetBranchAddress("lept1_phi",               &lept1_phi);
  intree->SetBranchAddress("lept1_sceta",             &lept1_sceta);

  intree->SetBranchAddress("lept1_RecoSF",               &lept1_RecoSF);
  intree->SetBranchAddress("lept1_SelSF",                &lept1_SelSF);
  intree->SetBranchAddress("lept1_trigSF",               &lept1_trigSF);

  intree->SetBranchAddress("z_pt",                &z_pt);
  intree->SetBranchAddress("z_eta",               &z_eta);
  //intree->SetBranchAddress("z_phi",               &z_phi);
  intree->SetBranchAddress("z_y",                 &z_y);
  intree->SetBranchAddress("z_mass",              &z_mass);
  intree->SetBranchAddress("z_charge",            &z_charge);
  intree->SetBranchAddress("leptType",            &leptType);

  intree->SetBranchAddress("lep0_theta_Zrest",    &lep0_theta_Zrest);
  intree->SetBranchAddress("lep0_costheta_Zrest", &lep0_costheta_Zrest);
  intree->SetBranchAddress("lep0_phi_Zrest",      &lep0_phi_Zrest);
  intree->SetBranchAddress("lep1_theta_Zrest",    &lep1_theta_Zrest);
  intree->SetBranchAddress("lep1_costheta_Zrest", &lep1_costheta_Zrest);
  intree->SetBranchAddress("lep1_phi_Zrest",      &lep1_phi_Zrest);
  intree->SetBranchAddress("lep_theta_Zrest",     &lep_theta_Zrest);
  intree->SetBranchAddress("lep_costheta_Zrest",  &lep_costheta_Zrest);
  intree->SetBranchAddress("lep_phi_Zrest",       &lep_phi_Zrest);

  if (!data) {
    ingentree->ResetBranchAddresses();
    ingentree->SetBranchAddress("genWeight",   &genWeight);
    ingentree->SetBranchAddress("puweigj_69p2nb",     &puweigj_69p2nb);
    
    ingentree->SetBranchAddress("ngenEle",     &ngenEle);
    ingentree->SetBranchAddress("ngenMu",      &ngenMu);
    ingentree->SetBranchAddress("genlepPt",    genlepPt);
    ingentree->SetBranchAddress("genlepPhi",   genlepPhi);
    ingentree->SetBranchAddress("genlepEta",   genlepEta);
    ingentree->SetBranchAddress("genlepCh",    genlepCh);
    ingentree->SetBranchAddress("genZm",       &genZm);
    ingentree->SetBranchAddress("genZpt",      &genZpt);
    ingentree->SetBranchAddress("genZy",       &genZy);


    ingentree->SetBranchAddress("gen_lep_theta_Zrest",      &gen_lep_theta_Zrest);
    ingentree->SetBranchAddress("gen_lep_costheta_Zrest",   &gen_lep_costheta_Zrest);
    ingentree->SetBranchAddress("gen_lep_phi_Zrest",        &gen_lep_phi_Zrest);

    ingentree->SetBranchAddress("gen_lep0_theta_Zrest",      &gen_lep0_theta_Zrest);
    ingentree->SetBranchAddress("gen_lep0_costheta_Zrest",   &gen_lep0_costheta_Zrest);
    ingentree->SetBranchAddress("gen_lep0_phi_Zrest",        &gen_lep0_phi_Zrest);
    
    ingentree->SetBranchAddress("gen_lep1_theta_Zrest",      &gen_lep1_theta_Zrest);
    ingentree->SetBranchAddress("gen_lep1_costheta_Zrest",   &gen_lep1_costheta_Zrest);
    ingentree->SetBranchAddress ("gen_lep1_phi_Zrest",        &gen_lep1_phi_Zrest);
  }
  
  //outfile
  TString outfilename;
  if (!data) outfilename = "minitrees/ZJets_aMCatNLO_Boost_test.root";
  else {
    if (ele) outfilename = "minitrees/DoubleEG_Boost.root";
    else outfilename = "minitrees/DoubleMu_RochesterCorr_Boost.root";
  }


  TFile *fout = new TFile(outfilename, "recreate");

  TTree *outtree = new TTree("outtree", "");
  inittree(outtree);

  TTree *outtreeGen = new TTree("outtreeGen", "output tree at gen level");
  inittreeGen(outtreeGen);

  TH1F *hphi = new TH1F("hphi", "hphi", 50, -3.15, 3.15);
  Long64_t totEv = intree->GetEntriesFast();

  //loop for reconstruction tree
  for (Long64_t iev = 0; iev < totEv; iev++) {
  //for (Long64_t iev = 0; iev < 1000000; iev++) {
    intree->GetEntry(iev);
    if(intree->GetEntry(iev)<=0) break;

    if (iev%20000 == 0) cout << "processing event: " << iev+1 << " th of " << totEv << " events" << endl;

    //if (z_mass>80) cout << z_mass << endl;
    if (z_mass<80. || z_mass>100.) continue;
    if (lept0_pt<25.) continue;
    if (lept1_pt<20.) continue;
    if (fabs(lept0_eta)>2.5) continue;
    if (fabs(lept1_eta)>2.5) continue;
    if (z_charge != 0) continue;

    //cout << "event passing selection" << endl;

    puweigj_69p2nb_ = puweigj_69p2nb;
    genWeight_ = genWeight;

    trig_Mu17_Mu8_ = trig_Mu17_Mu8;
    trig_Ele23_Ele12_ = trig_Ele23_Ele12;

    lept0_pt_ = lept0_pt;
    lept0_eta_ = lept0_eta;
    lept0_phi_ = lept0_phi;
    //lept0_sceta_ = lept0_sceta_;

    lept1_pt_ = lept1_pt;
    lept1_eta_ = lept1_eta;
    lept1_phi_ = lept1_phi;
    //lept1_sceta_ ;

    z_pt_ = z_pt;
    z_eta_ = z_eta;
    //z_phi_ = z_phi;
    z_y_ = z_y;
    z_mass_ = z_mass;
    z_charge_ = z_charge;
    leptType_ = leptType;

    lept0_RecoSF_ = lept0_RecoSF;
    lept0_SelSF_ = lept0_SelSF;
    lept0_trigSF_ = lept0_trigSF;

    lept1_RecoSF_ = lept1_RecoSF;
    lept1_SelSF_ = lept1_SelSF;
    lept1_trigSF_ = lept1_trigSF;

    lep0_theta_Zrest_ = -99.;
    lep0_costheta_Zrest_ = -99.;
    lep0_phi_Zrest_ = -99.;
    lep1_theta_Zrest_ = -99.;
    lep1_costheta_Zrest_ = -99.;
    lep1_phi_Zrest_ = -99.;

    lep_theta_Zrest_ = -99.;
    lep_costheta_Zrest_ = -99.;
    lep_phi_Zrest_ = -99.;


    TLorentzVector lep1, lep2, Zlab;

    if (leptType==11 && trig_Ele23_Ele12==1) {
      lep1.SetPtEtaPhiM(lept0_pt, lept0_eta, lept0_phi, 0.511*0.001);
      lep2.SetPtEtaPhiM(lept1_pt, lept1_eta, lept1_phi, 0.511*0.001);
    }
    else if (leptType==13 && trig_Mu17_Mu8==1) {
      lep1.SetPtEtaPhiM(lept0_pt, lept0_eta, lept0_phi, 104.*0.001);
      lep2.SetPtEtaPhiM(lept1_pt, lept1_eta, lept1_phi, 104.*0.001);
    }


    Zlab = lep1+lep2;
    z_phi_ = Zlab.Phi();
    float Zlab_pz = Zlab.Pz();

    //TLorentzVector vec4_Everything;
    //vec4_Everything . RotateZ (-Zlab.Phi());
    //cout << "Zlab phi before boost: " << Zlab.Phi() << endl;

    //cout << "Zlab phi after boost: " << Zlab.Phi() << endl;
    TVector3 vec3_BoostToRest;
    TVector3 vec3_BoostToPz, vec3_BoostToPt;
    //TVector3 vec3_BoostToPz(0., 0., Zlab.Pz()/Zlab.E());
    //TVector3 vec3_BoostToPt(Zlab.Px()/Zlab.E(), Zlab.Py()/Zlab.E(), 0.);

    vec3_BoostToRest = Zlab.BoostVector();
    vec3_BoostToPz.SetX(0.);
    vec3_BoostToPz.SetY(0.);
    vec3_BoostToPz.SetZ(vec3_BoostToRest.Z());

    //boost to Zpz first
    lep1.Boost(-vec3_BoostToPz);
    lep2.Boost(-vec3_BoostToPz);
    Zlab.Boost(-vec3_BoostToPz);
    //cout << "Pz of Z: " << Zlab.Pz() << endl;

    lep1.RotateZ (-z_phi_);
    lep2.RotateZ (-z_phi_);
    Zlab.RotateZ (-z_phi_);


    //then boost to Zpt
    vec3_BoostToPt = Zlab.BoostVector();
    Zlab.Boost(-vec3_BoostToPt);
    lep1.Boost(-vec3_BoostToPt);
    lep2.Boost(-vec3_BoostToPt);

    //cout << "Z px: " << Zlab.Px() << "  py : " << Zlab.Py() << endl;
    //TVector3 vec3_Boost;
    //vec3_Boost = Zlab.BoostVector();

    
    lep0_theta_Zrest_ = lep1.Theta();
    lep0_costheta_Zrest_ = lep1.CosTheta();
    lep0_phi_Zrest_ = lep1.Phi();
    lep1_theta_Zrest_ = lep2.Theta();
    lep1_costheta_Zrest_ = lep2.CosTheta();
    lep1_phi_Zrest_ = lep2.Phi();

    /*
    if (lept0_pdgId < 0) {
      lep_theta_Zrest = lep0_theta_Zrest;
      lep_costheta_Zrest = lep0_costheta_Zrest;
      lep_phi_Zrest = lep0_phi_Zrest;
    }
    else {
      lep_theta_Zrest = lep1_theta_Zrest;
      lep_costheta_Zrest = lep1_costheta_Zrest;
      lep_phi_Zrest = lep1_phi_Zrest;
    }
    */

    //lep1.Boost(-vec3_Boost);
    //lep2.Boost(-vec3_Boost);

    //cout << "angle of 2 leptons after boost: " << lep1.Angle(lep2.Vect()) << endl;
    //cout << "Px of lep1: " << lep1.Px() << " Py: " << lep1.Py() << " Pz: " << lep1.Pz() << endl;
    //cout << "Px of lep2: " << lep2.Px() << " Py: " << lep2.Py() << " Pz: " << lep2.Pz() << endl;
    //cout << "phi of lep1: " << lep1.Phi() << endl;
    //cout << "phi of lep2: " << lep2.Phi() << endl;

    //cout << endl;
    if (fabs(z_y) < 1. && z_pt<10.) hphi->Fill(lep1.Phi());

    outtree->Fill();

  }

  
  //loop for generated tree for mc only
  if (!data) {

    Long64_t totEv_gen = ingentree->GetEntriesFast();

    for (Long64_t iev = 0; iev < totEv_gen; iev++) {
      //for (Long64_t iev = 0; iev < 1000000; iev++) {
      ingentree->GetEntry(iev);

      if(ingentree->GetEntry(iev)<=0) break;
      
      if (iev%20000 == 0) cout << "processing event: " << iev+1 << " th of " << totEv_gen << " events" << endl;

      if (genZm<80. || genZm>100.) continue;
      if (fabs(genZy) >2.5) continue;

      puweigj_69p2nb_ = puweigj_69p2nb;
      genWeight_ = genWeight;

      ngenEle_ = ngenEle;
      ngenMu_ = ngenMu;

      genlepPt_[0] = genlepPt[0];
      genlepEta_[0] = genlepEta[0];
      genlepPhi_[0] = genlepPhi[0];
      genlepCh_[0] = genlepCh[0];

      genlepPt_[1] = genlepPt[1];
      genlepEta_[1] = genlepEta[1];
      genlepPhi_[1] = genlepPhi[1];
      genlepCh_[1] = genlepCh[1];


      genZm_ = genZm;
      genZpt_ = genZpt;
      genZy_ = genZy;
      genZch_ = genZch;
      genZtype_ = genZtype;

      gen_lep_theta_Zrest_ = -99.;
      gen_lep_costheta_Zrest_ = -99.;
      gen_lep_phi_Zrest_ = -99.;

      gen_lep0_theta_Zrest_ = -99.;
      gen_lep0_costheta_Zrest_ = -99.;
      gen_lep0_phi_Zrest_ = -99.;

      gen_lep1_theta_Zrest_ = -99.;
      gen_lep1_costheta_Zrest_ = -99.;
      gen_lep1_phi_Zrest_ = -99.;

      TLorentzVector lep1, lep2, Zlab;

      if (ngenEle>1) {
	lep1.SetPtEtaPhiM(genlepPt[0], genlepEta[0], genlepPhi[0], 0.511*0.001);
	lep2.SetPtEtaPhiM(genlepPt[1], genlepEta[1], genlepPhi[1], 0.511*0.001);
      }
      else if (ngenMu > 1) {
	lep1.SetPtEtaPhiM(genlepPt[0], genlepEta[0], genlepPhi[0], 104.*0.001);
        lep2.SetPtEtaPhiM(genlepPt[1], genlepEta[1], genlepPhi[1], 104.*0.001);
      }

      Zlab = lep1+lep2;
      float Zlab_pz = Zlab.Pz();
      float zphi = Zlab.Phi();
      lep1.RotateZ (-zphi);
      lep2.RotateZ (-zphi);
      Zlab.RotateZ (-zphi);

      TVector3 vec3_BoostToRest;
      TVector3 vec3_BoostToPz, vec3_BoostToPt;
      
      vec3_BoostToRest = Zlab.BoostVector();
      vec3_BoostToPz.SetX(0.);
      vec3_BoostToPz.SetY(0.);
      vec3_BoostToPz.SetZ(vec3_BoostToRest.Z());
      

      //boost to Zpz first
      Zlab.Boost(-vec3_BoostToPz);
      lep1.Boost(-vec3_BoostToPz);
      lep2.Boost(-vec3_BoostToPz);

      //cout << "Pz of Z: " << Zlab.Pz() << endl;
      
      vec3_BoostToPt = Zlab.BoostVector();
      Zlab.Boost(-vec3_BoostToPt);
      //cout << "Z px: " << Zlab.Px() << "  py : " << Zlab.Py() << endl;
      //TVector3 vec3_Boost;
      //vec3_Boost = Zlab.BoostVector();
            
      //then boost to Zpt
      lep1.Boost(-vec3_BoostToPt);
      lep2.Boost(-vec3_BoostToPt);
      
      gen_lep0_theta_Zrest_ = lep1.Theta();
      gen_lep0_costheta_Zrest_ = lep1.CosTheta();
      gen_lep0_phi_Zrest_ = lep1.Phi();
      gen_lep1_theta_Zrest_ = lep2.Theta();
      gen_lep1_costheta_Zrest_ = lep2.CosTheta();
      gen_lep1_phi_Zrest_ = lep2.Phi();

      if (genlepCh[0] < 0) {
	gen_lep_theta_Zrest_ = gen_lep0_theta_Zrest_;
	gen_lep_costheta_Zrest_ = gen_lep0_costheta_Zrest_;
	gen_lep_phi_Zrest_ = gen_lep0_phi_Zrest_;
      }
      else {
	gen_lep_theta_Zrest_ = gen_lep1_theta_Zrest_;
	gen_lep_costheta_Zrest_ = gen_lep1_costheta_Zrest_;
	gen_lep_phi_Zrest_ = gen_lep1_phi_Zrest_;
      }

      
      
      outtreeGen->Fill();
    }
  }
  
  /*
  TCanvas *c1 = new TCanvas("c1", "c1", 650, 650);
  c1->cd();
  hphi->SetMinimum(hphi->GetMinimum()/1.2);
  hphi->SetMaximum(hphi->GetMaximum()*1.2);
  hphi->Draw();
  */


  fout->cd();
  outtree->Write();
  hphi->Write();
  outtreeGen->Write();
  fout->Write();
  fout->Close();

  //hphi->Delete();
  

}



