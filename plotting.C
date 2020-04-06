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


void SetHistStyle_sig(TH1F* h1) {

  h1->SetFillColor(kRed-9);
  h1->SetLineColor(kRed-9);
  h1->SetFillColorAlpha(kRed-9, 0.5);

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


void plotting(bool ele = true, float min_Zpt = 0., float max_Zpt = 2.5, int Zy_bin = 1) {

  gStyle->SetOptStat(0);

  //use Zpt binning, less than ATLAS, more than CMS at 8 TeV
  int nZpt_bin = 16;
  float Zpt[] = {0.1, 2.5, 5.0, 8.0, 11., 15., 20., 25., 30., 35., 40., 45., 50., 80., 120., 200., 500.};

  //two bins of z rapidity, 1: 0 < |y| <1., 2: 1 < |y| < 2.5
  //int nZy_bin = 2;
  //float Zy[] = {0., 1., 2.5}

  TString fname;
  if (ele) fname = "minitrees/DoubleEG_Boost.root";
  else fname = "minitrees/DoubleMu_RochesterCorr_Boost.root";

  TFile *fda = new TFile(fname, "read");
  TFile *fsig = new TFile("minitrees/ZJets_aMCatNLO_Boost_test.root", "read");
  //TFile *fsig = new TFile("minitrees/ZJets_aMCatNLO_correctZaxis.root", "read");
  TFile *fttb = new TFile("minitrees/TT_Powheg_NoRochesterCorr.root", "read");
  TFile *fwj  = new TFile("minitrees/WJet_NoRochesterCorr.root", "read");


  TTree *tda  = (TTree*) fda->Get("outtree");
  TTree *tsig = (TTree*) fsig->Get("outtree");
  TTree *ttb  = (TTree*) fttb->Get("outtree");
  TTree *twj  = (TTree*) fwj->Get("outtree");

  //TH1D *hntotwei_sig = (TH1D*) fsig->Get("hntotweight");
  TH1D *hntotwei_ttb = (TH1D*) fttb->Get("hntotweight");
  TH1D *hntotwei_wj  = (TH1D*) fwj->Get("hntotweight");

  TCut cutda, cutmc;
  TCut cutpt_da, cutrap_da;
  TCut cutpt_mc, cutrap_mc;

  cutda = "z_mass>80 && z_mass<100 && lept0_pt>25 && lept1_pt>20";
  if (ele) cutda += "leptType==11 && fabs(lept0_eta)<2.5 && fabs(lept1_eta)<2.5 && trig_Ele23_Ele12==1";
  else cutda += "leptType==13 && fabs(lept0_eta)<2.4 && fabs(lept1_eta)<2.4 && trig_Mu17_Mu8==1";
  cutrap_da = cutda;

  if (Zy_bin == 1) cutda += "fabs(z_y)<1.";
  else cutda += "fabs(z_y)>=1. && fabs(z_y)<2.5";

  cutpt_da = cutda;

  cutda += Form("z_pt>%f && z_pt<%f", min_Zpt, max_Zpt);
  //cutrap_da += Form("z_pt>%f && z_pt<%f", min_Zpt, max_Zpt);

  TCut mcwei = "genWeight*puweigj_69p2nb";
  //TCut mcwei = "genWeight*puweigj_69p2nb*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*lept0_trigSF*lept1_trigSF*lept_dzSF";
  //TCut mcwei = "genWeight*puweigj_69p2nb*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF";
  cutmc = cutda * mcwei;

  cutpt_mc = cutpt_da * mcwei;
  cutrap_mc = cutrap_da * mcwei;

  //TCut cutpt = Form("z_pt>%f && z_pt<%f", min_Zpt, max_Zpt);
  //TCut cuty = Form("z_pt>%f && z_pt<%f", min_Zpt, max_Zpt);

  //TCut cutpt = cutda - Form("z_pt>%f && z_pt<%f", min_Zpt, max_Zpt);
  //TCut cuty = cutda - Form("z_pt>%f && z_pt<%f", min_Zpt, max_Zpt);

  TH1F *hcosTheta_da = new TH1F("hcosTheta_da", "cosTheta in da", 12, -1., 1.);
  TH1F *hphi_da      = new TH1F("hphi_da", "phi in da", 12, 0., 3.15);
  TH1F *hphi_da_l2      = new TH1F("hphi_da_l2", "phi in da", 12, 0., 3.15);
  TH1F *hzm_da       = new TH1F("hzm_da", "zmass in da", 30, 60, 120);
  //TH1F *hzpt_da      = new TH1F("hzpt_da", "zpt in da", nZpt_bin, Zpt);
  TH1F *hzpt_da      = new TH1F("hzpt_da", "zpt in da", 50, 0., 500);
  TH1F *hzy_da       = new TH1F("hzy_da", "z rapidity in da", 30, -2.5, 2.5);

  TH1F *hcosTheta_sig = new TH1F("hcosTheta_sig", "cosTheta in sig", 12, -1., 1.);
  TH1F *hphi_sig      = new TH1F("hphi_sig", "phi in sig", 12, 0., 3.15);
  TH1F *hphi_sig_l2      = new TH1F("hphi_sig_l2", "phi in sig", 12, 0., 3.15);
  TH1F *hzm_sig       = new TH1F("hzm_sig", "zmass in sig", 30, 60, 120);
  //TH1F *hzpt_sig      = new TH1F("hzpt_sig", "zpt in sig", nZpt_bin, Zpt);
  TH1F *hzpt_sig      = new TH1F("hzpt_sig", "zpt in sig", 50, 0., 500);
  TH1F *hzy_sig       = new TH1F("hzy_sig", "z rapidity in sig", 30, -2.5, 2.5);

  TH1F *hcosTheta_ttb = new TH1F("hcosTheta_ttb", "cosTheta in sig", 12, -1., 1.);
  TH1F *hphi_ttb      = new TH1F("hphi_ttb", "phi in sig", 12, 0., 3.15);
  TH1F *hzm_ttb       = new TH1F("hzm_ttb", "zmass in sig", 30, 60, 120);
  //TH1F *hzpt_ttb      = new TH1F("hzpt_ttb", "zpt in sig", nZpt_bin, Zpt);
  TH1F *hzpt_ttb      = new TH1F("hzpt_ttb", "zpt in sig", 50, 0., 500.);
  TH1F *hzy_ttb       = new TH1F("hzy_ttb", "z rapidity in sig", 30, -2.5, 2.5);

  TH1F *hcosTheta_wj = new TH1F("hcosTheta_wj", "cosTheta in sig", 12, -1., 1.);
  TH1F *hphi_wj      = new TH1F("hphi_wj", "phi in sig", 12, 0., 3.15);
  TH1F *hzm_wj       = new TH1F("hzm_wj", "zmass in sig", 30, 60, 120);
  //TH1F *hzpt_wj      = new TH1F("hzpt_wj", "zpt in sig", nZpt_bin, Zpt);
  TH1F *hzpt_wj      = new TH1F("hzpt_wj", "zpt in sig", 50, 0., 500);
  TH1F *hzy_wj       = new TH1F("hzy_wj", "z rapidity in sig", 30, -2.5, 2.5);


  //filling histogram
  tda->Draw("lep_costheta_Zrest >> hcosTheta_da", cutda, "goff");
  tda->Draw("fabs(lep0_phi_Zrest) >> hphi_da", cutda, "goff");
  tda->Draw("fabs(lep1_phi_Zrest) >> hphi_da_l2", cutda, "goff");
  tda->Draw("z_mass >> hzm_da", cutda, "goff");
  tda->Draw("z_pt >> hzpt_da", cutpt_da, "goff");
  tda->Draw("z_y >> hzy_da", cutrap_da, "goff");

  tsig->Draw("lep_costheta_Zrest >> hcosTheta_sig", cutmc, "goff");
  tsig->Draw("fabs(lep0_phi_Zrest) >> hphi_sig", cutmc, "goff");
  tsig->Draw("fabs(lep1_phi_Zrest) >> hphi_sig_l2", cutmc, "goff");
  tsig->Draw("z_mass >> hzm_sig", cutmc, "goff");
  tsig->Draw("z_pt >> hzpt_sig", cutpt_mc, "goff");
  tsig->Draw("z_y >> hzy_sig", cutrap_mc, "goff");

  ttb->Draw("lep_costheta_Zrest >> hcosTheta_ttb", cutmc, "goff");
  ttb->Draw("fabs(lep_phi_Zrest) >> hphi_ttb", cutmc, "goff");
  ttb->Draw("z_mass >> hzm_ttb", cutmc, "goff");
  ttb->Draw("z_pt >> hzpt_ttb", cutpt_mc, "goff");
  ttb->Draw("z_y >> hzy_ttb", cutrap_mc, "goff");

  twj->Draw("lep_costheta_Zrest >> hcosTheta_wj", cutmc, "goff");
  twj->Draw("fabs(lep_phi_Zrest) >> hphi_wj", cutmc, "goff");
  twj->Draw("z_mass >> hzm_wj", cutmc, "goff");
  twj->Draw("z_pt >> hzpt_wj", cutpt_mc, "goff");
  twj->Draw("z_y >> hzy_wj", cutrap_mc, "goff");

  //normalized to luminosity of data
  float lumi_2016 = 35900; //pb-1
  //float xs_sig = 6077.22; //pb
  float xs_sig = 6024.; //pb
  float xs_ttb = 831.76; //pb
  float xs_wj = 61526.7; //pb
  float xs_wz = 0.;
  float xs_ww = 0.;
  float xs_zz = 0.;

  //float scale_sig = lumi_2016*xs_sig/hntotwei_sig->Integral();
  float scale_sig = lumi_2016*xs_sig/80513408.0;
  float scale_ttb = lumi_2016*xs_ttb/hntotwei_ttb->Integral();
  float scale_wj = lumi_2016*xs_wj/hntotwei_wj->Integral();
  //float scale_wz = lumi_2016*xs_wz/hntotwei_wz->Integral();

  hcosTheta_da->Sumw2();
  hphi_da->Sumw2();
  hzm_da->Sumw2();
  hzpt_da->Sumw2();
  hzy_da->Sumw2();

  hphi_da->Add(hphi_da_l2);
  hphi_sig->Add(hphi_sig_l2);


  /*  //scaling to luminosity
  hcosTheta_sig->Scale(scale_sig);
  hphi_sig->Scale(scale_sig);
  hzm_sig->Scale(scale_sig);
  hzpt_sig->Scale(scale_sig);
  hzy_sig->Scale(scale_sig);

  hcosTheta_ttb->Scale(scale_ttb);
  hphi_ttb->Scale(scale_ttb);
  hzm_ttb->Scale(scale_ttb);
  hzpt_ttb->Scale(scale_ttb);
  hzy_ttb->Scale(scale_ttb);

  hcosTheta_wj->Scale(scale_wj);
  hphi_wj->Scale(scale_wj);
  hzm_wj->Scale(scale_wj);
  hzpt_wj->Scale(scale_wj);
  hzy_wj->Scale(scale_wj);
  */
  //hphi_da->Scale(1./hphi_da->Integral());
  //hphi_sig->Scale(1./hphi_sig->Integral());

  hcosTheta_da->Print();
  hcosTheta_sig->Print();
  hzpt_da->Print();
  hzpt_sig->Print();

  hphi_da->Print();
  hphi_sig->Print();


  cout << "entries of 1st bin of Zpt: " << hzpt_da->GetBinContent(1) << "\t" << hzpt_sig->GetBinContent(1) << endl;


  hcosTheta_da->Scale(1./hcosTheta_da->Integral());
  hphi_da->Scale(1./hphi_da->Integral());
  hzm_da->Scale(1./hzm_da->Integral());
  hzpt_da->Scale(1./hzpt_da->Integral());
  hzy_da->Scale(1./hzy_da->Integral());

  hcosTheta_sig->Scale(1./hcosTheta_sig->Integral());
  hphi_sig->Scale(1./hphi_sig->Integral());
  hzm_sig->Scale(1./hzm_sig->Integral());
  hzpt_sig->Scale(1./hzpt_sig->Integral());
  hzy_sig->Scale(1./hzy_sig->Integral());


  SetHistStyle_sig(hcosTheta_sig);
  SetHistStyle_sig(hphi_sig);
  SetHistStyle_sig(hzm_sig);
  SetHistStyle_sig(hzpt_sig);
  SetHistStyle_sig(hzy_sig);
  /*
  SetHistStyle_ttb(hcosTheta_ttb);
  SetHistStyle_ttb(hphi_ttb);
  SetHistStyle_ttb(hzm_ttb);
  SetHistStyle_ttb(hzpt_ttb);
  SetHistStyle_ttb(hzy_ttb);

  SetHistStyle_wj(hcosTheta_wj);
  SetHistStyle_wj(hphi_wj);
  SetHistStyle_wj(hzm_wj);
  SetHistStyle_wj(hzpt_wj);
  SetHistStyle_wj(hzy_wj);
  */
  /*
  THStack *hstack_cos = new THStack("hstack_cos", "");
  THStack *hstack_phi = new THStack("hstack_phi", "");
  THStack *hstack_zm  = new THStack("hstack_zm", "");
  THStack *hstack_zpt = new THStack("hstack_zpt", "");
  THStack *hstack_zy  = new THStack("hstack_zy", "");

  hstack_cos->Add(hcosTheta_sig);
  hstack_cos->Add(hcosTheta_ttb);
  hstack_cos->Add(hcosTheta_wj);

  hstack_phi->Add(hphi_sig);
  hstack_phi->Add(hphi_ttb);
  hstack_phi->Add(hphi_wj);

  hstack_zm->Add(hzm_sig);
  hstack_zm->Add(hzm_ttb);
  hstack_zm->Add(hzm_wj);

  hstack_zpt->Add(hzpt_sig);
  hstack_zpt->Add(hzpt_ttb);
  hstack_zpt->Add(hzpt_wj);

  hstack_zy->Add(hzy_sig);
  hstack_zy->Add(hzy_ttb);
  hstack_zy->Add(hzy_wj);
  */

  cout << "plotting cos theta" << endl;
  TCanvas *cos = new TCanvas("cos", "cos", 650, 650);
  cos->cd();

  hcosTheta_da->GetYaxis()->SetTitle("Events");
  hcosTheta_da->GetXaxis()->SetTitle("cos(#theta_{CS})");
  hcosTheta_da->SetTitleSize(0.05, "XYZ");
  hcosTheta_da->SetTitleFont(42, "XYZ");
  hcosTheta_da->SetTitleOffset(1.3, "XYZ");
  hcosTheta_da->SetTitleOffset(1.5, "Y");
  hcosTheta_da->SetMaximum(hcosTheta_da->GetMaximum()*1.2);
  hcosTheta_da->Draw();
  //hstack_cos->Draw("histsame");
  hcosTheta_sig->Draw("histsame");
  hcosTheta_da->Draw("same");
  cos->RedrawAxis();
  //TLegend *lg1 = new TLegend(0.5, 0.62, 0.68, 0.82);
  TLegend *lg1 = new TLegend(0.5, 0.35, 0.68, 0.57);
  //TLegend *lg1 = new TLegend(0.68, 0.58, 0.88, 0.78);
  lg1->AddEntry(hcosTheta_da, "data", "pe");
  lg1->AddEntry(hcosTheta_sig, "MC", "f");
  //lg1->AddEntry(hcosTheta_ttb, "t#bar{t}", "f");
  //lg1->AddEntry(hcosTheta_wj, "WJets", "f");
  lg1->SetBorderSize(0);
  lg1->SetTextFont(42);
  lg1->SetTextSize(0.05);
  lg1->Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.05);
  tx.SetTextFont(42);
  tx.DrawLatex(0.5, 0.3, Form("%.1f < p_{T}^{Z} < %.1f", min_Zpt, max_Zpt));
  if (Zy_bin == 1) tx.DrawLatex(0.5, 0.24, "0 < |Y^{Z}| < 1");
  else tx.DrawLatex(0.5, 0.24, "1 < |Y^{Z}| < 2.5");


  TCanvas *cphi = new TCanvas("cphi", "cphi", 650, 650);
  cphi->cd();

  hphi_da->GetYaxis()->SetTitle("Events");
  hphi_da->GetXaxis()->SetTitle("#phi_{CS}");
  hphi_da->SetTitleSize(0.05, "XYZ");
  hphi_da->SetTitleFont(42, "XYZ");
  hphi_da->SetTitleOffset(1.3, "XYZ");
  //hphi_da->SetMaximum(hphi_da->GetMaximum()*1.2);
  hphi_da->GetYaxis()->SetRangeUser(hphi_da->GetMinimum()/2, hphi_da->GetMaximum()*1.15);
  hphi_da->Draw();
  hphi_sig->Draw("histsame");
  //hstack_phi->Draw("histsame");
  hphi_da->Draw("same");
  cphi->RedrawAxis();
  TLegend *lg2 = new TLegend(0.25, 0.2, 0.48, 0.45);
  lg2->AddEntry(hphi_da, "data", "pe");
  lg2->AddEntry(hphi_sig, "MC", "f");
  //lg2->AddEntry(hphi_ttb, "t#bar{t}", "f");
  //lg2->AddEntry(hphi_wj, "WJets", "f");
  lg2->SetBorderSize(0);
  lg2->SetTextFont(42);
  lg2->SetTextSize(0.05);
  lg2->Draw();
  tx.DrawLatex(0.5, 0.38, Form("%.1f < p_{T}^{Z} < %.1f", min_Zpt, max_Zpt));
  if (Zy_bin == 1) tx.DrawLatex(0.5, 0.3, "0 < |Y^{Z}| < 1");
  else tx.DrawLatex(0.5, 0.3, "1 < |Y^{Z}| < 2.5");

  /*
  TCanvas *czm = new TCanvas("czm", "czm", 650, 650);
  czm->cd();

  hzm_da->GetYaxis()->SetTitle("Events/2 GeV");
  hzm_da->GetXaxis()->SetTitle("m_{ll} [GeV]");
  hzm_da->SetTitleSize(0.05, "XYZ");
  hzm_da->SetTitleFont(42, "XYZ");
  hzm_da->SetTitleOffset(1.4, "XYZ");
  hzm_da->Draw();
  hzm_sig->Draw("histsame");
  hzm_da->Draw("same");
  czm->RedrawAxis();
  TLegend *lg3 = new TLegend(0.6, 0.6, 0.78, 0.72);
  lg3->AddEntry(hcosTheta_da, "data", "pe");
  lg3->AddEntry(hcosTheta_sig, "DY", "f");
  lg3->SetBorderSize(0);
  lg3->SetTextFont(42);
  lg3->SetTextSize(0.05);
  lg3->Draw();

  tx.DrawLatex(0.6, 0.82, Form("%.1f < p_{T}^{Z} < %.1f", min_Zpt, max_Zpt));
  if (Zy_bin == 1) tx.DrawLatex(0.6, 0.75, "0 < |Y^{Z}| < 1");
  else tx.DrawLatex(0.6, 0.75, "1 < |Y^{Z}| < 2.5");
  */
  TCanvas *czpt = new TCanvas("czpt", "czpt", 650, 650);
  czpt->cd();
  //czpt->SetLogy();

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1., 1.);
  pad1->SetTopMargin(0.05);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();

  hzpt_da->GetYaxis()->SetTitle("Events");
  hzpt_da->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  hzpt_da->SetTitleSize(0.05, "XYZ");
  hzpt_da->SetTitleFont(42, "XYZ");
  hzpt_da->SetTitleOffset(1.3, "XYZ");
  hzpt_da->GetXaxis()->SetLabelSize(0);
  hzpt_da->Draw();
  hzpt_sig->Draw("histsame");
  hzpt_da->Draw("same");
  czpt->RedrawAxis();
  lg2->Draw();
  if (Zy_bin == 1) tx.DrawLatex(0.5, 0.36, "0 < |Y^{Z}| < 1");
  else tx.DrawLatex(0.5, 0.36, "1 < |Y^{Z}| < 2.5");

  czpt->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0., 1., 0.3);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.12);
  pad2->Draw();
  pad2->cd();


  TH1F *hratio = (TH1F*) hzpt_da->Clone();
  hratio->Divide(hzpt_sig);
  hratio->GetYaxis()->SetTitle("Data/MC");
  hratio->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");

  //hratio->GetXaxis()->SetLimits(xmin, xmax);
  hratio->GetXaxis()->SetLabelOffset(0.05);
  hratio->GetXaxis()->SetTitleOffset(1.5);
  hratio->GetYaxis()->SetLabelSize(0.1);
  hratio->GetXaxis()->SetLabelSize(0.1);
  hratio->GetXaxis()->SetTitleSize(0.1);
  hratio->Draw("e");

  TLine *l = new TLine(0.,1,500.,1);
  l->SetLineWidth(1);
  l->SetLineColor(1);
  l->SetLineStyle(kDashed);
  l->Draw("same");

  pad2->RedrawAxis();


  TCanvas *czy = new TCanvas("czy", "czy", 650, 650);
  czy->cd();

  hzy_da->GetYaxis()->SetTitle("Events");
  hzy_da->GetXaxis()->SetTitle("Y^{Z} [GeV]");
  hzy_da->SetTitleSize(0.05, "XYZ");
  hzy_da->SetTitleFont(42, "XYZ");
  hzy_da->SetTitleOffset(1.3, "XYZ");
  hzy_da->Draw();
  hzy_sig->Draw("histsame");
  hzy_da->Draw("same");
  czy->RedrawAxis();
  lg1->Draw();
  tx.DrawLatex(0.5, 0.38, Form("%.1f < p_{T}^{Z} < %.1f", min_Zpt, max_Zpt));
  

  //saving
  TString cname;
  cname = Form("_Zpt_%.1fto%.1f", min_Zpt, max_Zpt);
  if (Zy_bin==1) cname += Form("_ZY_%.1fto%.1f", 0., 1.);
  else cname += Form("_ZY_%.1fto%.1f", 1., 2.5);
  if (ele) cname += "_ele";
  else cname += "_mu";
  //cname += "_norm1";
  cname += "_minuslep"; 


  cos->SaveAs("plots/inclusive/costheta_Zrest" + cname + "_reboost.pdf");
  //cphi->SaveAs("plots/inclusive/phi_Zrest" + cname + "_reboost.pdf");
  //czm->SaveAs("plots/Zmass.pdf" + cname + ".pdf");
  //czpt->SaveAs("plots/Zpt.pdf" + cname + "_log.pdf");
  //czy->SaveAs("plots/Zrapidity.pdf" + cname + "_log.pdf");

}
