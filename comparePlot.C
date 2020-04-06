#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "THStack.h"


void SetHistStyle_ele(TH1F* h1) {

  //h1->SetFillColor(kRed-9);
  h1->SetLineColor(kRed-7);
  h1->SetLineWidth(2);
  //h1->SetFillColorAlpha(kRed-9, 0.5);

}

void SetHistStyle_mu(TH1F* h1) {

  //h1->SetFillColor(kGreen-9);
  h1->SetLineColor(kCyan-2);
  h1->SetLineWidth(2);
  //h1->SetFillColorAlpha(kGreen-9, 0.5);

}

void SetHistStyle_wj(TH1F* h1) {

  h1->SetFillColor(kCyan-9);
  h1->SetLineColor(kCyan-9);
  h1->SetFillColorAlpha(kCyan-9, 0.5);

}


void comparePlot() {

  gStyle->SetOptStat(0);

  TFile *fsig = new TFile("minitrees/ZJets_aMCatNLO_NoRochesterCorr.root", "read");

  TTree *tsig = (TTree*) fsig->Get("outtree");

  TH1D *hntotwei_sig = (TH1D*) fsig->Get("hntotweight");

  TCut cutele, cutmu;

  cutele = "leptType==11 && fabs(lept0_eta)<1.2 && fabs(lept1_eta)<1.2 && trig_Ele23_Ele12==1";
  cutmu = "leptType==13 && fabs(lept0_eta)<1.2 && fabs(lept1_eta)<1.2 && trig_Mu17_Mu8==1";

  cutele += "z_mass>80 && z_mass<100 && lept0_pt>25 && lept1_pt>20";
  cutmu  += "z_mass>80 && z_mass<100 && lept0_pt>25 && lept1_pt>20";

  TCut mcwei = "genWeight*puweigj_69p2nb";
  //TCut mcwei = "genWeight*puweigj_69p2nb*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*lept0_trigSF*lept1_trigSF*lept_dzSF";
  //TCut mcwei = "genWeight*puweigj_69p2nb*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF";
  cutele *= mcwei;
  cutmu *= mcwei;

  TH1F *hzm_ele       = new TH1F("hzm_ele", "zmass in ele", 30., 60, 120);
  TH1F *hzpt_ele      = new TH1F("hzpt_ele", "zpt in ele", 50, 0., 1000.);
  TH1F *hzy_ele       = new TH1F("hzy_ele", "z rapidity in ele", 30., -2.5, 2.5);
  TH1F *hlep1_pt_ele  = new TH1F("hlep1_pt_ele", "", 50., 25., 1025.);
  TH1F *hlep2_pt_ele  = new TH1F("hlep2_pt_ele", "", 50., 20., 520.);
  TH1F *hlep1_eta_ele  = new TH1F("hlep1_eta_ele", "", 20., -2.5, 2.5);
  TH1F *hlep2_eta_ele  = new TH1F("hlep2_eta_ele", "", 20., -2.5, 2.5);
  TH1F *hlep1_phi_ele  = new TH1F("hlep1_phi_ele", "", 20., -3.15, 3.15);
  TH1F *hlep2_phi_ele  = new TH1F("hlep2_phi_ele", "", 20., -3.15, 3.15);

  TH1F *hzm_mu       = new TH1F("hzm_mu", "zmass in mu", 30., 60, 120);
  TH1F *hzpt_mu      = new TH1F("hzpt_mu", "zpt in mu", 50, 0., 1000.);
  TH1F *hzy_mu       = new TH1F("hzy_mu", "z rapidity in mu", 30., -2.5, 2.5);
  TH1F *hlep1_pt_mu  = new TH1F("hlep1_pt_mu", "", 50., 25., 1025.);
  TH1F *hlep2_pt_mu  = new TH1F("hlep2_pt_mu", "", 50., 20., 520.);
  TH1F *hlep1_eta_mu  = new TH1F("hlep1_eta_mu", "", 20., -2.5, 2.5);
  TH1F *hlep2_eta_mu  = new TH1F("hlep2_eta_mu", "", 20., -2.5, 2.5);
  TH1F *hlep1_phi_mu  = new TH1F("hlep1_phi_mu", "", 20., -3.15, 3.15);
  TH1F *hlep2_phi_mu  = new TH1F("hlep2_phi_mu", "", 20., -3.15, 3.15);


  //filling histogram
  tsig->Draw("z_mass >> hzm_ele", cutele, "goff");
  tsig->Draw("z_pt >> hzpt_ele", cutele, "goff");
  tsig->Draw("z_y >> hzy_ele", cutele, "goff");
  tsig->Draw("lept0_pt >> hlep1_pt_ele", cutele, "goff");
  tsig->Draw("lept1_pt >> hlep2_pt_ele", cutele, "goff");
  tsig->Draw("lept0_eta >> hlep1_eta_ele", cutele, "goff");
  tsig->Draw("lept1_eta >> hlep2_eta_ele", cutele, "goff");
  tsig->Draw("lept0_phi >> hlep1_phi_ele", cutele, "goff");
  tsig->Draw("lept1_phi >> hlep2_phi_ele", cutele, "goff");

  tsig->Draw("z_mass >> hzm_mu", cutmu, "goff");
  tsig->Draw("z_pt >> hzpt_mu", cutmu, "goff");
  tsig->Draw("z_y >> hzy_mu", cutmu, "goff");
  tsig->Draw("lept0_pt >> hlep1_pt_mu", cutmu, "goff");
  tsig->Draw("lept1_pt >> hlep2_pt_mu", cutmu, "goff");
  tsig->Draw("lept0_eta >> hlep1_eta_mu", cutmu, "goff");
  tsig->Draw("lept1_eta >> hlep2_eta_mu", cutmu, "goff");
  tsig->Draw("lept0_phi >> hlep1_phi_mu", cutmu, "goff");
  tsig->Draw("lept1_phi >> hlep2_phi_mu", cutmu, "goff");

  hzm_ele->Scale(1./hzm_ele->Integral());
  hzpt_ele->Scale(1./hzpt_ele->Integral());
  hzy_ele->Scale(1./hzy_ele->Integral());
  hlep1_pt_ele->Scale(1./hlep1_pt_ele->Integral());
  hlep2_pt_ele->Scale(1./hlep2_pt_ele->Integral());
  hlep1_eta_ele->Scale(1./hlep1_eta_ele->Integral());
  hlep2_eta_ele->Scale(1./hlep2_eta_ele->Integral());
  hlep1_phi_ele->Scale(1./hlep1_phi_ele->Integral());
  hlep2_phi_ele->Scale(1./hlep2_phi_ele->Integral());

  hzm_mu->Scale(1./hzm_mu->Integral());
  hzpt_mu->Scale(1./hzpt_mu->Integral());
  hzy_mu->Scale(1./hzy_mu->Integral());
  hlep1_pt_mu->Scale(1./hlep1_pt_mu->Integral());
  hlep2_pt_mu->Scale(1./hlep2_pt_mu->Integral());
  hlep1_eta_mu->Scale(1./hlep1_eta_mu->Integral());
  hlep2_eta_mu->Scale(1./hlep2_eta_mu->Integral());
  hlep1_phi_mu->Scale(1./hlep1_phi_mu->Integral());
  hlep2_phi_mu->Scale(1./hlep2_phi_mu->Integral());


  SetHistStyle_ele(hzm_ele);
  SetHistStyle_ele(hzpt_ele);
  SetHistStyle_ele(hzy_ele);
  SetHistStyle_ele(hlep1_pt_ele);
  SetHistStyle_ele(hlep2_pt_ele);
  SetHistStyle_ele(hlep1_eta_ele);
  SetHistStyle_ele(hlep2_eta_ele);
  SetHistStyle_ele(hlep1_phi_ele);
  SetHistStyle_ele(hlep2_phi_ele);

  SetHistStyle_mu(hzm_mu);
  SetHistStyle_mu(hzpt_mu);
  SetHistStyle_mu(hzy_mu);
  SetHistStyle_mu(hlep1_pt_mu);
  SetHistStyle_mu(hlep2_pt_mu);
  SetHistStyle_mu(hlep1_eta_mu);
  SetHistStyle_mu(hlep2_eta_mu);
  SetHistStyle_mu(hlep1_phi_mu);
  SetHistStyle_mu(hlep2_phi_mu);


  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.05);
  tx.SetTextFont(42);

  TCanvas *czm = new TCanvas("czm", "czm", 650, 650);
  czm->cd();

  hzm_mu->GetYaxis()->SetTitle("Events/2 GeV");
  hzm_mu->GetXaxis()->SetTitle("m_{ll} [GeV/c^{2}]");
  hzm_mu->SetTitleSize(0.05, "XYZ");
  hzm_mu->SetTitleFont(42, "XYZ");
  hzm_mu->SetTitleOffset(1.4, "XYZ");
  hzm_mu->Draw("hist");
  hzm_ele->Draw("histsame");
  czm->RedrawAxis();
  TLegend *lg1 = new TLegend(0.65, 0.6, 0.8, 0.72);
  lg1->AddEntry(hzm_ele, "Z#rightarrow ee", "f");
  lg1->AddEntry(hzm_mu, "Z#rightarrow #mu#mu", "f");
  lg1->SetBorderSize(0);
  lg1->SetTextFont(42);
  lg1->SetTextSize(0.05);
  lg1->Draw();


  TCanvas *czpt = new TCanvas("czpt", "czpt", 650, 650);
  czpt->cd();
  czpt->SetLogy();

  hzpt_ele->GetYaxis()->SetTitle("Events");
  hzpt_ele->GetXaxis()->SetTitle("p_{T}^{Z} [GeV/c]");
  hzpt_ele->SetTitleSize(0.05, "XYZ");
  hzpt_ele->SetTitleFont(42, "XYZ");
  hzpt_ele->SetTitleOffset(1.3, "XYZ");
  hzpt_ele->Draw("hist");
  hzpt_mu->Draw("histsame");
  czpt->RedrawAxis();
  lg1->Draw();


  TCanvas *clep1_pt = new TCanvas("clep1_pt", "clep1_pt", 650, 650);
  clep1_pt->cd();
  clep1_pt->SetLogy();

  hlep1_pt_ele->GetYaxis()->SetTitle("Events");
  hlep1_pt_ele->GetXaxis()->SetTitle("leading p_{T} [GeV/c]");
  hlep1_pt_ele->SetTitleSize(0.05, "XYZ");
  hlep1_pt_ele->SetTitleFont(42, "XYZ");
  hlep1_pt_ele->SetTitleOffset(1.3, "XYZ");
  hlep1_pt_ele->Draw("hist");
  hlep1_pt_mu->Draw("histsame");
  clep1_pt->RedrawAxis();
  lg1->Draw();

  TCanvas *clep2_pt = new TCanvas("clep2_pt", "clep2_pt", 650, 650);
  clep2_pt->cd();
  //clep2_pt->SetLogy();

  hlep2_pt_ele->GetYaxis()->SetTitle("Events");
  hlep2_pt_ele->GetXaxis()->SetTitle("trailing p_{T} [GeV/c]");
  hlep2_pt_ele->SetTitleSize(0.05, "XYZ");
  hlep2_pt_ele->SetTitleFont(42, "XYZ");
  hlep2_pt_ele->SetTitleOffset(1.3, "XYZ");
  hlep2_pt_ele->Draw("hist");
  hlep2_pt_mu->Draw("histsame");
  clep2_pt->RedrawAxis();
  lg1->Draw();


  TCanvas *czy = new TCanvas("czy", "czy", 650, 650);
  czy->cd();

  hzy_ele->GetYaxis()->SetTitle("Events");
  hzy_ele->GetXaxis()->SetTitle("Y^{Z}");
  hzy_ele->SetTitleSize(0.05, "XYZ");
  hzy_ele->SetTitleFont(42, "XYZ");
  hzy_ele->SetTitleOffset(1.3, "XYZ");
  hzy_ele->Draw("hist");
  hzy_mu->Draw("histsame");
  czy->RedrawAxis();
  TLegend *lg2 = new TLegend(0.45, 0.5, 0.65, 0.62);
  lg2->AddEntry(hzm_ele, "Z#rightarrow ee", "f");
  lg2->AddEntry(hzm_mu, "Z#rightarrow #mu#mu", "f");
  lg2->SetBorderSize(0);
  lg2->SetTextFont(42);
  lg2->SetTextSize(0.05);
  lg2->Draw();



  TCanvas *clep1_eta = new TCanvas("clep1_eta", "clep1_eta", 650, 650);
  clep1_eta->cd();
  //clep1_eta->SetLogy();

  hlep1_eta_ele->GetYaxis()->SetTitle("Events");
  hlep1_eta_ele->GetXaxis()->SetTitle("leading #eta");
  hlep1_eta_ele->SetTitleSize(0.05, "XYZ");
  hlep1_eta_ele->SetTitleFont(42, "XYZ");
  hlep1_eta_ele->SetTitleOffset(1.3, "XYZ");
  hlep1_eta_ele->Draw("hist");
  hlep1_eta_mu->Draw("histsame");
  clep1_eta->RedrawAxis();
  lg2->Draw();

  TCanvas *clep2_eta = new TCanvas("clep2_eta", "clep2_eta", 650, 650);
  clep2_eta->cd();
  //clep2_eta->SetLogy();

  hlep2_eta_ele->GetYaxis()->SetTitle("Events");
  hlep2_eta_ele->GetXaxis()->SetTitle("trailing #eta");
  hlep2_eta_ele->SetTitleSize(0.05, "XYZ");
  hlep2_eta_ele->SetTitleFont(42, "XYZ");
  hlep2_eta_ele->SetTitleOffset(1.3, "XYZ");
  hlep2_eta_ele->Draw("hist");
  hlep2_eta_mu->Draw("histsame");
  clep2_eta->RedrawAxis();
  lg2->Draw();


  TCanvas *clep1_phi = new TCanvas("clep1_phi", "clep1_phi", 650, 650);
  clep1_phi->cd();
  //clep1_phi->SetLogy();

  hlep1_phi_ele->GetYaxis()->SetTitle("Events");
  hlep1_phi_ele->GetXaxis()->SetTitle("leading #phi [GeV/c]");
  hlep1_phi_ele->SetTitleSize(0.05, "XYZ");
  hlep1_phi_ele->SetTitleFont(42, "XYZ");
  hlep1_phi_ele->SetTitleOffset(1.3, "XYZ");
  hlep1_phi_ele->SetMaximum(hlep1_phi_ele->GetMaximum()*2);
  hlep1_phi_ele->SetMinimum(hlep1_phi_ele->GetMinimum()/2);
  hlep1_phi_ele->Draw("hist");
  hlep1_phi_mu->Draw("histsame");
  clep1_phi->RedrawAxis();
  lg2->Draw();

  TCanvas *clep2_phi = new TCanvas("clep2_phi", "clep2_phi", 650, 650);
  clep2_phi->cd();
  //clep2_phi->SetLogy();

  hlep2_phi_ele->GetYaxis()->SetTitle("Events");
  hlep2_phi_ele->GetXaxis()->SetTitle("trailing #phi");
  hlep2_phi_ele->SetTitleSize(0.05, "XYZ");
  hlep2_phi_ele->SetTitleFont(42, "XYZ");
  hlep2_phi_ele->SetTitleOffset(1.3, "XYZ");
  hlep2_phi_ele->SetMaximum(hlep2_phi_ele->GetMaximum()*2);
  hlep2_phi_ele->SetMinimum(hlep2_phi_ele->GetMinimum()/2);
  hlep2_phi_ele->Draw("hist");
  hlep2_phi_mu->Draw("histsame");
  clep2_phi->RedrawAxis();
  lg2->Draw();


  //saving
  czm->SaveAs("plots/Zmass_ele_mu_lepEta1p2.pdf");
  czpt->SaveAs("plots/Zpt_ele_mu_lepEta1p2pdf");
  czy->SaveAs("plots/Zrapidity_ele_mu_lepEta1p2pdf");
  clep1_pt->SaveAs("plots/leadingLepPt_ele_mu_lepEta1p2pdf");
  clep2_pt->SaveAs("plots/trailingLepPt_ele_mu_lepEta1p2pdf");
  clep1_eta->SaveAs("plots/leadingLepEta_ele_mu_lepEta1p2pdf");
  clep2_eta->SaveAs("plots/trailingLepEta_ele_mu_lepEta1p2pdf");
  clep1_phi->SaveAs("plots/leadingLepPhi_ele_mu_lepEta1p2pdf");
  clep2_phi->SaveAs("plots/trailingLepPhi_ele_mu_lepEta1p2pdf");


}
