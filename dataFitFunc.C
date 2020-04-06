#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "THStack.h"
#include "TPad.h"
#include "TF1.h"



Double_t fphi(double *x, Double_t *par) {

  Float_t xx = x[0];
  Double_t f = par[0] *(1+ 3.*TMath::Pi()/16*par[1]*TMath::Cos(xx) + 0.25*par[2]*TMath::Cos(2*xx));
  //Double_t f = par[0]*(1 - par[2] + par[1]*xx + 0.5*par[2]*xx*xx);
  return f;

}


Double_t fcostheta(double *x, Double_t *par) {

  Float_t xx = x[0];
  Double_t f = par[0]*(1+xx*xx + 0.5*par[1]*(1-3*xx*xx) + par[2]*xx);
  //Double_t f = par[0]*(xx*xx) + par[1]*xx + par[3];
  return f;

}


void SetHistStyle_sig(TH1F* h1) {

  //h1->SetFillColor(kRed-9);
  h1->SetLineColor(kRed-9);
  h1->SetLineWidth(2);
  //h1->SetFillColorAlpha(kRed-9, 0.5);

}


void dataFitFunc(bool ele = true, float min_Zpt = 0., float max_Zpt = 2.5, int Zy_bin = 1) {

  gStyle->SetOptStat(0);

  TString inf_da;
  //if (ele) inf_da = "../main/minitrees/DoubleEG_NoRochesterCorr.root";
  //else inf_da = "../main/minitrees/DoubleMu_RochesterCorr_correctZaxis.root";
  if (ele) inf_da = "../main/minitrees/DoubleEG_Boost.root";
  else inf_da = "../main/minitrees/DoubleMu_RochesterCorr_Boost.root";

  TFile *fda = new TFile(inf_da, "read");
  //TFile *fsig = new TFile("../main/minitrees/ZJets_aMCatNLO_correctZaxis.root", "read");
  TFile *fsig = new TFile("../main/minitrees/ZJets_aMCatNLO_Boost_test.root", "read");

  TTree *tda = (TTree*)  fda->Get("outtree");
  TTree *tsig = (TTree*) fsig->Get("outtree");
  TTree *tgen = (TTree*) fsig->Get("outtreeGen");

  TCut cutmc, cutgen;
  TCut cutda;

  //cut at reco level
  cutmc = "z_mass>80 && z_mass<100 && lept0_pt>25 && lept1_pt>20";
  if (ele) cutmc += "leptType==11 && trig_Ele23_Ele12==1 && fabs(lept0_eta)<2.5 && fabs(lept1_eta)<2.5";
  else cutmc += "leptType==13 && trig_Mu17_Mu8==1 && fabs(lept0_eta)<2.5 && fabs(lept1_eta)<2.5";

  if (Zy_bin == 1) cutmc += "fabs(z_y)<1.";
  else cutmc += "fabs(z_y)>=1. && fabs(z_y)<2.0";

  cutmc += Form("z_pt>%f && z_pt<%f", min_Zpt, max_Zpt);

  cutda = cutmc;
  //cut at gen level
  cutgen = "genZm>80 && genZm<100 && genlepPt[0]>0 && genlepPt[1]>0";
  if (ele) cutgen += "ngenEle==2";
  else cutgen += "ngenMu==2";

  if (Zy_bin == 1) cutgen += "fabs(genZy)<1.";
  else cutgen += "fabs(genZy)>=1. && fabs(genZy)<2.0";

  cutgen += Form("genZpt>%f && genZpt<%f", min_Zpt, max_Zpt);


  TCut mcwei = "genWeight*puweigj_69p2nb";
  //TCut mcwei = "genWeight*puweigj_69p2nb*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*lept0_trigSF*lept1_trigSF*lept_dzSF";
  //TCut mcwei = "genWeight*puweigj_69p2nb*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF";
  cutmc *= mcwei;
  cutgen *= mcwei;


  TH1F *hcosTheta_da = new TH1F("hcosTheta_da", "cosTheta in da", 20, -1., 1.);
  TH1F *hphi_da      = new TH1F("hphi_da", "phi in da", 15, 0., 3.15);

  TH1F *hcosTheta_da_l2 = new TH1F("hcosTheta_da_l2", "cosTheta in da", 20, -1., 1.);
  TH1F *hphi_da_l2      = new TH1F("hphi_da_l2", "phi in da", 15, 0., 3.15);
  //TH1F *hphi_da      = new TH1F("hphi_da", "phi in da", 40, -1., 1.);

  TH1F *hcosTheta_sig = new TH1F("hcosTheta_sig", "cosTheta in sig", 20, -1., 1.);
  TH1F *hphi_sig      = new TH1F("hphi_sig", "phi in sig", 15, 0., 3.15);

  TH1F *hcosTheta_sig_l2 = new TH1F("hcosTheta_sig_l2", "cosTheta in sig", 20, -1., 1.);
  TH1F *hphi_sig_l2      = new TH1F("hphi_sig_l2", "phi in sig", 15, 0., 3.15);
  //TH1F *hphi_sig      = new TH1F("hphi_sig", "phi in sig", 40, -1., 1.);

  TH1F *hcosTheta_gen = new TH1F("hcosTheta_gen", "cosTheta in gen", 20, -1., 1.);
  TH1F *hphi_gen      = new TH1F("hphi_gen", "phi in gen", 15, 0., 3.15);

  TH1F *hcosTheta_gen_l2 = new TH1F("hcosTheta_gen_l2", "cosTheta in gen", 20, -1., 1.);
  TH1F *hphi_gen_l2      = new TH1F("hphi_gen_l2", "phi in gen", 15, 0., 3.15);
  //TH1F *hphi_gen      = new TH1F("hphi_gen", "phi in gen", 40, -1., 1.);

  hcosTheta_da->Sumw2();
  hcosTheta_sig->Sumw2();
  hcosTheta_gen->Sumw2();

  hcosTheta_da_l2->Sumw2();
  hcosTheta_sig_l2->Sumw2();
  hcosTheta_gen_l2->Sumw2();

  hphi_da->Sumw2();
  hphi_sig->Sumw2();
  hphi_gen->Sumw2();

  hphi_da_l2->Sumw2();
  hphi_sig_l2->Sumw2();
  hphi_gen_l2->Sumw2();

  //filling histogram
  //tda->Draw("lep0_costheta_Zrest >> hcosTheta_da", cutda, "goff");
  //tda->Draw("lep1_costheta_Zrest >> hcosTheta_da_l2", cutda, "goff");
  tda->Draw("fabs(lep0_phi_Zrest) >> hphi_da", cutda, "goff");
  tda->Draw("fabs(lep1_phi_Zrest) >> hphi_da_l2", cutda, "goff");
  //tda->Draw("cos(lep0_phi_Zrest) >> hphi_da", cutda, "goff");

  //tsig->Draw("lep0_costheta_Zrest >> hcosTheta_sig", cutmc, "goff");
  //tsig->Draw("lep1_costheta_Zrest >> hcosTheta_sig_l2", cutmc, "goff");
  tsig->Draw("fabs(lep0_phi_Zrest) >> hphi_sig", cutmc, "goff");
  tsig->Draw("fabs(lep1_phi_Zrest) >> hphi_sig_l2", cutmc, "goff");
  //tsig->Draw("cos(lep0_phi_Zrest) >> hphi_sig", cutmc, "goff");

  //tgen->Draw("gen_lep0_costheta_Zrest >> hcosTheta_gen", cutgen, "goff");
  //tgen->Draw("gen_lep1_costheta_Zrest >> hcosTheta_gen_l2", cutgen, "goff");
  tgen->Draw("fabs(gen_lep0_phi_Zrest) >> hphi_gen", cutgen, "goff");
  tgen->Draw("fabs(gen_lep1_phi_Zrest) >> hphi_gen_l2", cutgen, "goff");
  //tgen->Draw("cos(gen_lep0_phi_Zrest) >> hphi_gen", cutgen, "goff");

  //hcosTheta_sig->Print();

  hcosTheta_da->Add(hcosTheta_da_l2);
  hcosTheta_sig->Add(hcosTheta_sig_l2);
  hcosTheta_gen->Add(hcosTheta_gen_l2);

  hphi_da->Add(hphi_da_l2);
  hphi_sig->Add(hphi_sig_l2);
  hphi_gen->Add(hphi_gen_l2);

  TH1F *hacc_theta = (TH1F*) hcosTheta_sig->Clone();
  TH1F *hacc_phi = (TH1F*) hphi_sig->Clone();

  //hacc_theta->Divide(hcosTheta_gen);
  hacc_phi->Divide(hphi_gen);

  //hcosTheta_sig->Divide(hacc_theta);
  //hcosTheta_da->Divide(hacc_theta);

  hphi_sig->Divide(hacc_phi);
  hphi_da->Divide(hacc_phi);

  SetHistStyle_sig(hcosTheta_sig);
  SetHistStyle_sig(hphi_sig);


  //call function for fit
  TF1 *f1 = new TF1("f1", fcostheta, -0.9, 0.9, 3);
  f1->SetParameters(1,0.01, 0.1);
  f1->SetParNames("Norm", "A0", "A4");
				    
  f1->SetLineColor(kBlue);
  f1->SetLineWidth(2);

  TF1 *f2 = new TF1("f2", fphi, 0., 2.8, 3);
  //TF1 *f2 = new TF1("f2", fphi, -0.8, 0.8, 3);
  f2->SetParameters(1,0.01, 0.1);
  f2->SetParNames("Norm", "A3", "A2");
				    
  f2->SetLineColor(kBlue);
  f2->SetLineWidth(2);
  
  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.05);
  tx.SetTextFont(42);

  /*
  TCanvas *cos = new TCanvas("cos", "cos", 650, 650);
  cos->cd();

  hcosTheta_da->GetYaxis()->SetTitle("Events");
  hcosTheta_da->GetXaxis()->SetTitle("cos#theta_{CS}");
  hcosTheta_da->SetTitleSize(0.05, "XYZ");
  hcosTheta_da->SetTitleFont(42, "XYZ");
  hcosTheta_da->SetTitleOffset(1.3, "XYZ");
  hcosTheta_da->SetTitleOffset(1.5, "Y");
  hcosTheta_da->SetMaximum(hcosTheta_da->GetMaximum()*1.2);
  hcosTheta_da->SetMinimum(0);
  hcosTheta_da->Draw("e");
  //f1->SetParameter(0,0.08);
  f1->SetParLimits(1, 0., 1.);
  //f1->SetParameter(2,0.001);
  hcosTheta_da->Fit(f1, "LM", "", -0.9, 0.9);
  //f1->Draw("same");

  cos->RedrawAxis();
  
  tx.DrawLatex(0.5, 0.3, Form("%.1f < p_{T}^{Z} < %.1f", min_Zpt, max_Zpt));
  if (Zy_bin == 1) tx.DrawLatex(0.5, 0.24, "0 < |Y^{Z}| < 1");
  else tx.DrawLatex(0.5, 0.24, "1 < |Y^{Z}| < 2.0");
  */
  
  TCanvas *cphi = new TCanvas("cphi", "cphi", 650, 650);
  cphi->cd();

  hphi_da->GetYaxis()->SetTitle("Events");
  hphi_da->GetXaxis()->SetTitle("|#phi_{CS}|");
  hphi_da->SetTitleSize(0.05, "XYZ");
  hphi_da->SetTitleFont(42, "XYZ");
  hphi_da->SetTitleOffset(1.3, "XYZ");
  hphi_da->SetTitleOffset(1.5, "Y");
  hphi_da->SetMaximum(hphi_da->GetMaximum()*1.2);
  hphi_da->SetMinimum(0);
  hphi_da->Draw("e");
  f2->SetParLimits(2, -1.,1.);
  f2->SetParameter(2, 0.);
  //f2->SetParameter(0, 10000.);
  //f1->SetParameter(2,0.001);
  hphi_da->Fit(f2, "LM", "", 0., 2.8);
  f2->Draw("same");

  cphi->RedrawAxis();

  tx.DrawLatex(0.5, 0.32, Form("%.1f < p_{T}^{Z} < %.1f", min_Zpt, max_Zpt));
  if (Zy_bin == 1) tx.DrawLatex(0.5, 0.26, "0 < |Y^{Z}| < 1");
  else tx.DrawLatex(0.5, 0.26, "1 < |Y^{Z}| < 2.5");


  //saving
  TString cname;
  cname = Form("_Zpt_%.1fto%.1f", min_Zpt, max_Zpt);
  if (Zy_bin==1) cname += Form("_ZY_%.1fto%.1f", 0., 1.);
  else cname += Form("_ZY_%.1fto%.1f", 1., 2.0);
  if (ele) cname += "_ele";
  else cname += "_mu";
  //cname += "_norm1";
  cname += "_minuslep"; 


  //cos->SaveAs("plots/inclusive/costheta_Zrest" + cname + "_func.pdf");
  //cos->SaveAs("plots/Reboost/costheta_Zrest" + cname + "_func_2lep.pdf");
  cphi->SaveAs("plots/Reboost/phi_Zrest" + cname + "_func_2lep.pdf");


}
