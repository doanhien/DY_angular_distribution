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



Double_t fcostheta(double *x, Double_t *par) {

  Float_t xx = x[0];
  //Double_t f = par[0]*(1+xx*xx) + 0.5*par[1]*(1-3*xx*xx) + par[2]*xx;
  Double_t f = par[0]*(xx*xx) + par[1]*xx + par[3];
  return f;

}


void SetHistStyle_gen(TH1F* h1) {

  //h1->SetFillColor(kRed-9);
  h1->SetLineColor(kRed-9);
  h1->SetLineWidth(2);
  //h1->SetFillColorAlpha(kRed-9, 0.5);

}

void SetHistStyle_rec(TH1F* h1) {

  //h1->SetFillColor(kGreen-9);
  h1->SetLineColor(kGreen-9);
  h1->SetLineWidth(2);
  //h1->SetFillColorAlpha(kGreen-9, 0.5);

}

void SetHistStyle_wj(TH1F* h1) {

  h1->SetFillColor(kCyan-9);
  h1->SetLineColor(kCyan-9);
  h1->SetFillColorAlpha(kCyan-9, 0.5);

}


void acceff(bool ele = true, float min_Zpt = 0., float max_Zpt = 2.5, int Zy_bin = 1) {

  gStyle->SetOptStat(0);

  TFile *fsig = new TFile("../main/minitrees/ZJets_aMCatNLO_correctZaxis.root", "read");

  TTree *tgen = (TTree*) fsig->Get("outtreeGen");
  TTree *trec = (TTree*) fsig->Get("outtree");


  TCut cutrec, cutgen;

  cutgen = "genZm>80 && genZm<100 && genlepPt[0]>0 && genlepPt[1]>0";
  if (ele) cutgen += "ngenEle==2 && fabs(genlepEta[0])<2.5 && fabs(genlepEta[1])<2.5";
  //else cutgen += "ngenMu==2 && fabs(genlepEta[0])<2.5 && fabs(genlepEta[1])<2.5";
  else cutgen += "ngenMu==2";

  if (Zy_bin == 1) cutgen += "fabs(genZy)<1.";
  else cutgen += "fabs(genZy)>=1. && fabs(genZy)<2.5";

  cutgen += Form("genZpt>%f && genZpt<%f", min_Zpt, max_Zpt);

  cutrec = "z_mass>80 && z_mass<100 && lept0_pt>25 && lept1_pt>20";
  if (ele) cutrec += "leptType==11 && trig_Ele23_Ele12==1 && fabs(lept0_eta)<2.5 && fabs(lept1_eta)<2.5";
  else cutrec += "leptType==13 && trig_Mu17_Mu8==1 && fabs(lept0_eta)<2.4 && fabs(lept1_eta)<2.4";

  if (Zy_bin == 1) cutrec+= "fabs(z_y)<1. ";
  else cutrec+= "fabs(z_y)>=1. && fabs(z_y)<2.5";

  cutrec += Form("z_pt>%f && z_pt<%f", min_Zpt, max_Zpt);

  //event weight
  TCut genwei = "genWeight*puweigj_69p2nb";
  //TCut mcwei = "genWeight*puweigj_69p2nb*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*lept0_trigSF*lept1_trigSF*lept_dzSF";
  TCut mcwei = "genWeight*puweigj_69p2nb*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF";

  cutgen *= genwei;
  cutrec *= mcwei;


  TH1F *hcosTheta_gen = new TH1F("hcosTheta_gen", "cosTheta in sig", 12, -1., 1.);
  TH1F *hphi_gen      = new TH1F("hphi_gen", "phi in sig", 12, 0., 3.15);

  TH1F *hcosTheta_rec = new TH1F("hcosTheta_rec", "cosTheta in sig", 12, -1., 1.);
  TH1F *hphi_rec      = new TH1F("hphi_rec", "phi in sig", 12, 0., 3.15);


  //filling histogram
  tgen->Draw("gen_lep_costheta_Zrest >> hcosTheta_gen", cutgen, "goff");
  tgen->Draw("fabs(gen_lep_phi_Zrest) >> hphi_gen", cutgen, "goff");

  trec->Draw("lep_costheta_Zrest >> hcosTheta_rec", cutrec, "goff");
  trec->Draw("fabs(lep_phi_Zrest) >> hphi_rec", cutrec, "goff");



  hcosTheta_gen->Print();


  SetHistStyle_gen(hcosTheta_gen);
  SetHistStyle_gen(hphi_gen);

  SetHistStyle_rec(hcosTheta_rec);
  SetHistStyle_rec(hphi_rec);


  TCanvas *cos = new TCanvas("cos", "cos", 650, 650);
  cos->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  pad1->SetFillColor(0);

  hcosTheta_gen->GetYaxis()->SetTitle("Events");
  //hcosTheta_gen->GetXaxis()->SetTitle("cos(#theta_{CS})");
  hcosTheta_gen->SetTitleSize(0.05, "XYZ");
  hcosTheta_gen->SetTitleFont(42, "XYZ");
  hcosTheta_gen->SetTitleOffset(1.3, "XYZ");
  hcosTheta_gen->SetTitleOffset(1.5, "Y");
  hcosTheta_gen->SetMaximum(hcosTheta_gen->GetMaximum()*1.2);
  hcosTheta_gen->SetMinimum(0);
  hcosTheta_gen->GetXaxis()->SetLabelSize(0);
  hcosTheta_gen->Draw("hist");
  hcosTheta_rec->Draw("histsame");
  TLegend *lg = new TLegend(0.4, 0.74, 0.55, 0.86);
  lg->AddEntry(hcosTheta_gen, "gen", "f");
  lg->AddEntry(hcosTheta_rec, "rec", "f");
  lg->SetBorderSize(0);
  lg->Draw();

  cos->RedrawAxis();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.05);
  tx.SetTextFont(42);
  tx.DrawLatex(0.5, 0.3, Form("%.1f < gen p_{T}^{Z} < %.1f", min_Zpt, max_Zpt));
  if (Zy_bin == 1) tx.DrawLatex(0.5, 0.24, "0 < gen |Y^{Z}| < 1");
  else tx.DrawLatex(0.5, 0.24, "1 < gen |Y^{Z}| < 2.5");

  cos->cd();
  TPad *pad2 = new TPad("pad2","pad2",0, 0, 1, 0.25);
  pad2->SetTopMargin(0.1);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad2->cd();

  TH1F *hratio_cos = (TH1F*) hcosTheta_rec->Clone();
  hratio_cos->Divide(hcosTheta_gen);

  //hratio_cos->SetMinimum(-6);
  hratio_cos->SetTitleSize(0.14, "XYZ");
  hratio_cos->SetLabelSize(0.14, "XYZ");
  hratio_cos->SetTitleFont(42, "XYZ");
  hratio_cos->GetYaxis()->SetTitle("efficiency");
  hratio_cos->GetXaxis()->SetTitle("cos#theta_{CS}");
  hratio_cos->GetYaxis()->SetTitleOffset(0.55);
  hratio_cos->GetXaxis()->SetTitleOffset(1.1);
  hratio_cos->GetYaxis()->SetTitleOffset(1.6);
  hratio_cos->GetYaxis()->SetRangeUser(0.0, 1.);
  hratio_cos->GetYaxis()->SetNdivisions(8);
  hratio_cos->Draw();



  TCanvas *cphi = new TCanvas("cphi", "cphi", 650, 650);
  cphi->cd();

  TPad *pad3 = new TPad("pad3", "pad3", 0, 0.25, 1, 1);
  pad3->SetBottomMargin(0.02);
  pad3->Draw();
  pad3->cd();
  pad3->SetFillColor(0);

  hphi_gen->GetYaxis()->SetTitle("Events");
  hphi_gen->GetXaxis()->SetTitle("#phi_{CS}");
  hphi_gen->SetTitleSize(0.05, "XYZ");
  hphi_gen->SetTitleFont(42, "XYZ");
  hphi_gen->SetTitleOffset(1.3, "XYZ");
  hphi_gen->GetYaxis()->SetRangeUser(hphi_gen->GetMinimum()/2, hphi_gen->GetMaximum()*1.2);
  hphi_gen->Draw("hist");
  hphi_rec->Draw("histsame");
  cphi->RedrawAxis();
  tx.DrawLatex(0.5, 0.38, Form("%.1f < gen p_{T}^{Z} < %.1f", min_Zpt, max_Zpt));
  if (Zy_bin == 1) tx.DrawLatex(0.5, 0.3, "0 < gen |Y^{Z}| < 1");
  else tx.DrawLatex(0.5, 0.3, "1 < gen |Y^{Z}| < 2.5");
  lg->Draw();

  cphi->cd();
  TPad *pad4 = new TPad("pad4","pad4",0, 0, 1, 0.25);
  pad4->SetTopMargin(0.1);
  pad4->SetBottomMargin(0.35);
  pad4->Draw();
  pad4->cd();

  TH1F *hratio_phi = (TH1F*) hphi_rec->Clone();
  hratio_phi->Divide(hphi_gen);

  //hratio_phi->SetMinimum(-6);
  hratio_phi->SetTitleSize(0.14, "XYZ");
  hratio_phi->SetLabelSize(0.14, "XYZ");
  hratio_phi->SetTitleFont(42, "XYZ");
  hratio_phi->GetYaxis()->SetTitle("efficiency");
  hratio_phi->GetXaxis()->SetTitle("|#phi_{CS}|");
  hratio_phi->GetYaxis()->SetTitleOffset(0.55);
  hratio_phi->GetXaxis()->SetTitleOffset(1.1);
  hratio_phi->GetYaxis()->SetTitleOffset(1.6);
  hratio_phi->GetYaxis()->SetRangeUser(0.4, 1.);
  hratio_phi->GetYaxis()->SetNdivisions(8);
  hratio_phi->Draw();


  //saving
  TString cname;
  cname = Form("_Zpt_%.1fto%.1f", min_Zpt, max_Zpt);
  if (Zy_bin==1) cname += Form("_ZY_%.1fto%.1f", 0., 1.);
  else cname += Form("_ZY_%.1fto%.1f", 1., 2.5);
  if (ele) cname += "_ele";
  else cname += "_mu";
  //cname += "_norm1";
  cname += "_minuslep"; 


  cos->SaveAs("plots/AccEff/Acc_costheta" + cname + ".pdf");
  cphi->SaveAs("plots/AccEff/Acc_phi" + cname + ".pdf");

  TString fileout;
  fileout = "AccEff/AccEff_costheta";
  fileout += cname;
  fileout += ".root";

  TFile *fout = new TFile(fileout, "recreate");
  fout->cd();

  hratio_cos->Write("AccEff_cosTheta");
  hratio_phi->Write("AccEff_phi");

  fout->Write();
  fout->Close();



}
