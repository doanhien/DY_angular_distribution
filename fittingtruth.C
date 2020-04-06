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
  Double_t f = par[0]*(1+xx*xx) + 0.5*par[1]*(1-3*xx*xx) + par[2]*xx;
  //Double_t f = par[0]*(xx*xx) + par[1]*xx + par[3];
  return f;

}


void SetHistStyle_sig(TH1F* h1) {

  //h1->SetFillColor(kRed-9);
  h1->SetLineColor(kRed-9);
  h1->SetLineWidth(2);
  //h1->SetFillColorAlpha(kRed-9, 0.5);

}

void SetHistStyle_ttb(TH1F* h1) {

  h1->SetFillColor(kGreen-9);
  h1->SetLineColor(kGreen-9);
  h1->SetFillColorAlpha(kGreen-9, 0.5);

}

void SetHistStyle_wj(TH1F* h1) {

  h1->SetFillColor(kCyan-9);
  h1->SetLineColor(kCyan-9);
  h1->SetFillColorAlpha(kCyan-9, 0.5);

}


void fittingtruth(bool ele = true, float min_Zpt = 0., float max_Zpt = 2.5, int Zy_bin = 1) {

  gStyle->SetOptStat(0);

  //use Zpt binning, less than ATLAS, more than CMS at 8 TeV
  int nZpt_bin = 16;
  float Zpt[] = {0., 2.5, 5.0, 8.0, 11., 15., 20., 25., 30., 35., 40., 45., 50., 80., 120., 200., 500.};

  //two bins of z rapidity, 1: 0 < |y| <1., 2: 1 < |y| < 2.5
  //int nZy_bin = 2;
  //float Zy[] = {0., 1., 2.5}

  //TFile *fsig = new TFile("minitrees/ZJets_aMCatNLO_NoRochesterCorr.root", "read");
  TFile *fsig = new TFile("../main/minitrees/ZJets_aMCatNLO_correctZaxis.root", "read");
  TFile *fttb = new TFile("../main/minitrees/TT_Powheg_NoRochesterCorr.root", "read");
  TFile *fwj  = new TFile("../main/minitrees/WJet_NoRochesterCorr.root", "read");


  TTree *tsig = (TTree*) fsig->Get("outtreeGen");
  TTree *ttb  = (TTree*) fttb->Get("outtreeGen");
  TTree *twj  = (TTree*) fwj->Get("outtreeGen");

  TH1D *hntotwei_sig = (TH1D*) fsig->Get("hntotweight");
  TH1D *hntotwei_ttb = (TH1D*) fttb->Get("hntotweight");
  TH1D *hntotwei_wj  = (TH1D*) fwj->Get("hntotweight");

  TCut cutmc;

  cutmc = "genZm>80 && genZm<100 && genlepPt[0]>0 && genlepPt[1]>0";
  if (ele) cutmc += "ngenEle==2 && fabs(genlepEta[0])<2.5 && fabs(genlepEta[1])<2.5";
  //else cutmc += "ngenMu==2 && fabs(genlepEta[0])<2.5 && fabs(genlepEta[1])<2.5";
  else cutmc += "ngenMu==2";

  if (Zy_bin == 1) cutmc += "fabs(genZy)<1.";
  else cutmc += "fabs(genZy)>=1. && fabs(genZy)<2.5";

  cutmc += Form("genZpt>%f && genZpt<%f", min_Zpt, max_Zpt);

  TCut mcwei = "genWeight*puweigj_69p2nb";
  //TCut mcwei = "genWeight*puweigj_69p2nb*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*lept0_trigSF*lept1_trigSF*lept_dzSF";
  //TCut mcwei = "genWeight*puweigj_69p2nb*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF";
  cutmc *= mcwei;


  TH1F *hcosTheta_sig = new TH1F("hcosTheta_sig", "cosTheta in sig", 12, -1., 1.);
  TH1F *hphi_sig      = new TH1F("hphi_sig", "phi in sig", 12, 0., 3.15);

  TH1F *hcosTheta_ttb = new TH1F("hcosTheta_ttb", "cosTheta in sig", 12, -1., 1.);
  TH1F *hphi_ttb      = new TH1F("hphi_ttb", "phi in sig", 12, 0., 3.15);

  TH1F *hcosTheta_wj = new TH1F("hcosTheta_wj", "cosTheta in sig", 12, -1., 1.);
  TH1F *hphi_wj      = new TH1F("hphi_wj", "phi in sig", 12, 0., 3.15);



  //filling histogram
  tsig->Draw("gen_lep_costheta_Zrest >> hcosTheta_sig", cutmc, "goff");
  tsig->Draw("fabs(gen_lep_phi_Zrest) >> hphi_sig", cutmc, "goff");


  //normalized to luminosity of data
  float lumi_2016 = 35900; //pb-1
  //float xs_sig = 6077.22; //pb
  float xs_sig = 6024.; //pb
  float xs_ttb = 831.76; //pb
  float xs_wj = 61526.7; //pb
  float xs_wz = 0.;
  float xs_ww = 0.;
  float xs_zz = 0.;

  float scale_sig = lumi_2016*xs_sig/hntotwei_sig->Integral();
  float scale_ttb = lumi_2016*xs_ttb/hntotwei_ttb->Integral();
  float scale_wj = lumi_2016*xs_wj/hntotwei_wj->Integral();
  //float scale_wz = lumi_2016*xs_wz/hntotwei_wz->Integral();

  //scaling to luminosity
  /*
  hcosTheta_sig->Scale(scale_sig);
  hphi_sig->Scale(scale_sig);

  hcosTheta_ttb->Scale(scale_ttb);
  hphi_ttb->Scale(scale_ttb);

  hcosTheta_wj->Scale(scale_wj);
  hphi_wj->Scale(scale_wj);
  */
  hcosTheta_sig->Print();

  hcosTheta_sig->Scale(1./hcosTheta_sig->Integral());
  hphi_sig->Scale(1./hphi_sig->Integral());


  SetHistStyle_sig(hcosTheta_sig);
  SetHistStyle_sig(hphi_sig);

  SetHistStyle_ttb(hcosTheta_ttb);
  SetHistStyle_ttb(hphi_ttb);

  SetHistStyle_wj(hcosTheta_wj);
  SetHistStyle_wj(hphi_wj);


  TF1 *f1 = new TF1("f1", fcostheta, -1, 1, 3);
  f1->SetParameters(1,0.01, 0.1);
  f1->SetParNames("Norm", "A0", "A4");
				    
  f1->SetLineColor(kBlue);
  f1->SetLineWidth(2);

  TCanvas *cos = new TCanvas("cos", "cos", 650, 650);
  cos->cd();

  hcosTheta_sig->GetYaxis()->SetTitle("Events");
  hcosTheta_sig->GetXaxis()->SetTitle("cos(#theta_{CS})");
  hcosTheta_sig->SetTitleSize(0.05, "XYZ");
  hcosTheta_sig->SetTitleFont(42, "XYZ");
  hcosTheta_sig->SetTitleOffset(1.3, "XYZ");
  hcosTheta_sig->SetTitleOffset(1.5, "Y");
  hcosTheta_sig->SetMaximum(hcosTheta_sig->GetMaximum()*1.2);
  hcosTheta_sig->SetMinimum(0);
  hcosTheta_sig->Draw("hist");
  //f1->SetParameter(0,0.08);
  //f1->SetParameter(1,0.1);
  //f1->SetParameter(2,0.001);
  hcosTheta_sig->Fit(f1, "", "");
  f1->Draw("same");


  cos->RedrawAxis();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.05);
  tx.SetTextFont(42);
  tx.DrawLatex(0.5, 0.3, Form("%.1f < gen p_{T}^{Z} < %.1f", min_Zpt, max_Zpt));
  if (Zy_bin == 1) tx.DrawLatex(0.5, 0.24, "0 < gen |Y^{Z}| < 1");
  else tx.DrawLatex(0.5, 0.24, "1 < gen |Y^{Z}| < 2.5");


  TCanvas *cphi = new TCanvas("cphi", "cphi", 650, 650);
  cphi->cd();

  hphi_sig->GetYaxis()->SetTitle("Events");
  hphi_sig->GetXaxis()->SetTitle("#phi_{CS}");
  hphi_sig->SetTitleSize(0.05, "XYZ");
  hphi_sig->SetTitleFont(42, "XYZ");
  hphi_sig->SetTitleOffset(1.3, "XYZ");
  hphi_sig->GetYaxis()->SetRangeUser(hphi_sig->GetMinimum()/2, hphi_sig->GetMaximum()*1.15);
  hphi_sig->Draw("hist");
  cphi->RedrawAxis();
  tx.DrawLatex(0.5, 0.38, Form("%.1f < gen p_{T}^{Z} < %.1f", min_Zpt, max_Zpt));
  if (Zy_bin == 1) tx.DrawLatex(0.5, 0.3, "0 < gen |Y^{Z}| < 1");
  else tx.DrawLatex(0.5, 0.3, "1 < gen |Y^{Z}| < 2.5");


  //saving
  TString cname;
  cname = Form("_Zpt_%.1fto%.1f", min_Zpt, max_Zpt);
  if (Zy_bin==1) cname += Form("_ZY_%.1fto%.1f", 0., 1.);
  else cname += Form("_ZY_%.1fto%.1f", 1., 2.5);
  if (ele) cname += "_ele";
  else cname += "_mu";
  //cname += "_norm1";
  cname += "_minuslep"; 


  //cos->SaveAs("plots/inclusive/costheta_Zrest" + cname + "_RochesterCorr.pdf");
  //cphi->SaveAs("plots/inclusive/phi_Zrest" + cname + "_RochesterCorr.pdf");


}
