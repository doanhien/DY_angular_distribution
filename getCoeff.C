#include "TH1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

void getCoeff() {

  const Int_t nbin = 16;
  Float_t Zpt_bins[nbin+1] = {0.2, 2.5, 5., 8., 11., 15., 20., 25., 30., 35., 40., 45., 50., 80., 120., 200., 3000};

  //for muon
  /*
  // boost 1 time
  Float_t A0_Y1_mu[nbin] = {0., 0.0149, 0.0308, 0.0684, 0.0845, 0.1196, 0.1784, 0.2309, 0.2942, 0.3338, 0.3787, 0.3926, 0.5318, 0.6767, 0.7983, 0.8256};
  Float_t ErrA0_Y1_mu[nbin] = {0.003, 0.0027, 0.0011, 0.0028, 0.0027, 0.0029, 0.0034, 0.0040, 0.0046, 0.0052, 0.0059, 0.0066, 0.0038, 0.0060, 0.0093, 0.0193};

  Float_t A0_Y2_mu[nbin] = {0.0631, 0.0785, 0.1236, 0.1268, 0.2056, 0.2389, 0.272, 0.3368, 0.3836, 0.4067, 0.4531, 0.4796, 0.5393, 0.5477, 0.5247, 0.3191};
  Float_t ErrA0_Y2_mu[nbin] = {0.0031, 0.0023, 0.0022, 0.0024, 0.0024, 0.0026, 0.003, 0.0036, 0.0041, 0.0047, 0.0053, 0.0059, 0.0034, 0.0054, 0.0085, 0.0194};


  Float_t A2_Y1_mu[nbin] = {0.0021, 0.004, 0.0144, 0.0084, 0.0041, -0.0046, -0.0064, 0.0029, 0.0065, -0.0153, 0.024, -0.059, -0.004, -0.0015, -0.01, -0.0166};
  Float_t ErrA2_Y1_mu[nbin] = {0.0112, 0.008, 0.0075, 0.0082, 0.008, 0.0081, 0.0095, 0.0109, 0.0125, 0.0143, 0.016, 0.0183, 0.0104, 0.0167, 0.026, 0.053};

  Float_t A2_Y2_mu[nbin] = {-0.0144, -0.0034, -0.0004, -0.002, 0.0137, -0.0055, -0.0127, -0.0193, -0.009, 0.0111, -0.026, -0.04, 0.0082, -0.022, 0.0375, -0.037};
  Float_t ErrA2_Y2_mu[nbin] = {0.0092, 0.0066, 0.0063, 0.0069, 0.0067, 0.007, 0.0082, 0.0094, 0.011, 0.0125, 0.0141, 0.0161, 0.0092, 0.015, 0.0237, 0.0523};


  //for electron
  Float_t A0_Y1_ele[nbin] = {0.0306, 0.0303, 0.0718, 0.0957, 0.126, 0.1724, 0.2019, 0.2379, 0.2777, 0.3441, 0.3782, 0.4231, 0.54, 0.7009, 0.8478, 0.8911};
  Float_t ErrA0_Y1_ele[nbin] = {0.0039, 0.0021, 0.0026, 0.0029, 0.0029, 0.003, 0.0036, 0.0042, 0.0048, 0.0055, 0.0061, 0.0069, 0.0001, 0.0063, 0.0097, 0.0194};

  Float_t A0_Y2_ele[nbin] = {0.0741, 0.0421, 0.0739, 0.1697, 0.171, 0.2281, 0.308, 0.3645, 0.3968, 0.4013, 0.4375, 0.4794, 0.5243, 0.4972, 0.4175, 0.2552};
  Float_t ErrA0_Y2_ele[nbin] = {0.0041,0.0029, 0.0028, 0.0025, 0.003, 0.0027, 0.003, 0.0037, 0.0043, 0.0049, 0.0055, 0.0061, 0.0035, 0.0055, 0.0086, 0.0192};

  Float_t A2_Y1_ele[nbin] = {};
  Float_t ErrA2_Y1_ele[nbin] = {};

  Float_t A2_Y2_ele[nbin] = {};
  Float_t ErrA2_Y2_ele[nbin] = {};
  */

  //boost twice and rotate
  //for muon
  Float_t A0_Y1_mu[nbin] = {};
  Float_t ErrA0_Y1_mu[nbin] = {};

  Float_t A0_Y2_mu[nbin] = {0.0891, 0.1067, 0.1442, 0.1494, 0.2284, 0.24, 0.308, 0.3883, 0.4217, 0.4677, 0.5168, 0.5079, 0.6315, 0.7781, 0.9138, 0.9901};
  Float_t ErrA0_Y2_mu[nbin] = {0.0027, 0.002, 0.0019, 0.0021, 0.0021, 0.0, 0.003, 0.003, 0.0035, 0.004, 0.0045, 0.005, 0.0028, 0.0045, 0.0071, 0.0159};

  Float_t A2_Y1_mu[nbin] = {};
  Float_t ErrA2_Y1_mu[nbin] = {};

  Float_t A2_Y2_mu[nbin] = {0.0112, 0.052, 0.0504, 0.0562, 0.096, 0.0784, 0.1335, 0.1566, 0.1823, 0.2302, 0.2671, 0.2834, 0.4267, 0.571, 0.6736, 0.7625};
  Float_t ErrA2_Y2_mu[nbin] = {0.0049, 0.0036, 0.0034, 0.0037, 0.0037, 0.0039, 0.0046, 0.0053, 0.0061, 0.0069, 0.0078, 0.0087, 0.0049, 0.008, 0.0121, 0.0273};

  //for ele
  Float_t A0_Y1_ele[nbin] = {0.03467, 0.0359, 0.0693, 0.0972, 0.1228, 0.1644, 0.1994, 0.2516, 0.28, 0.3312, 0.4095, 0.4269, 0.562, 0.7141, 0.9027, 0.9762};
  Float_t ErrA0_Y1_ele[nbin] = {0.0027, 0.0020, 0.0018, 0.002, 0.002, 0.0021, 0.0025, 0.0030, 0.0, 0.0039, 0.0044, 0.0049, 0.003, 0.0044, 0.0068, 0.0136};

  Float_t A0_Y2_ele[nbin] = {0.267, 0.1678, 0.1674, 0.2167, 0.282, 0.2714, 0.343, 0.4447, 0.4754, 0.5135, 0.5156, 0.581, 0.6232, 0.7537, 0.8172, 0.891};
  Float_t ErrA0_Y2_ele[nbin] = {0.0029,0.0021, 0.0019, 0.0022, 0.002, 0.0023, 0.003, 0.0032, 0.0037, 0.0042, 0.0046, 0.005, 0.003, 0.0047, 0.0073, 0.016};

  Float_t A2_Y1_ele[nbin] = {0.0427, 0.0581, 0.0893, 0.1254, 0.1446, 0.1538, 0.1613, 0.2453, 0.2483, 0.2411, 0.3364, 0.3632, 0.4699, 0.5668, 0.6652, 0.801};
  Float_t ErrA2_Y1_ele[nbin] = {0.0051, 0.0036,  0.0034, 0.0037, 0.0037, 0.0039, 0.0046, 0.0053, 0.006, 0.0069, 0.0077, 0.0086, 0.0048, 0.0076, 0.0116, 0.023};

  Float_t A2_Y2_ele[nbin] = {-0.0118, 0.061, 0.1121, 0.1778, 0.1633, 0.1223, 0.1449, 0.2293, 0.2957, 0.2689, 0.2443, 0.4499, 0.4014, 0.6067, 0.7109, 0.8685};
  Float_t ErrA2_Y2_ele[nbin] = {0.005, 0.0037, 0.0035, 0.0038, 0.0038, 0.004, 0.0047, 0.0055, 0.0063, 0.0072, 0.0081, 0.0089, 0.0051, 0.0079, 0.0123, 0.026};

  //convert to lambda, nu from A0, A2
  Float_t lamda_Y1_mu[nbin] = {0.};
  Float_t Errlamda_Y1_mu[nbin] = {0.};

  Float_t lamda_Y2_mu[nbin] = {0.};
  Float_t Errlamda_Y2_mu[nbin] = {0.};

  Float_t lamda_Y1_ele[nbin] = {0.};
  Float_t Errlamda_Y1_ele[nbin] = {0.};

  Float_t lamda_Y2_ele[nbin] = {0.};
  Float_t Errlamda_Y2_ele[nbin] = {0.};

  Float_t nu_Y1_mu[nbin] = {0.};
  Float_t Errnu_Y1_mu[nbin] = {0.};

  Float_t nu_Y2_mu[nbin] = {0.};
  Float_t Errnu_Y2_mu[nbin] = {0.};

  Float_t nu_Y1_ele[nbin] = {0.};
  Float_t Errnu_Y1_ele[nbin] = {0.};

  Float_t nu_Y2_ele[nbin] = {0.};
  Float_t Errnu_Y2_ele[nbin] = {0.};

  Float_t A0mA2_Y1_ele[nbin] = {0.};
  Float_t A0mA2_Y2_ele[nbin] = {0.};

  Float_t A0mA2_Y2_mu[nbin] = {0.};

  Float_t ErrA0mA2_Y1_ele[nbin] = {0.};
  Float_t ErrA0mA2_Y2_ele[nbin] = {0.};

  Float_t ErrA0mA2_Y2_mu[nbin] = {0.};



  for (int i = 0; i < nbin; i++) {
    lamda_Y1_mu[i] = (2.-3*A0_Y1_mu[i])/(2+A0_Y1_mu[i]);
    lamda_Y2_mu[i] = (2.-3*A0_Y2_mu[i])/(2+A0_Y2_mu[i]);
    
    Errlamda_Y1_mu[i] = 8.*ErrA0_Y1_mu[i]/TMath::Power(2+A0_Y1_mu[i],2);
    Errlamda_Y2_mu[i] = 8.*ErrA0_Y2_mu[i]/TMath::Power(2+A0_Y2_mu[i],2);

    //ele
    lamda_Y1_ele[i] = (2.-3*A0_Y1_ele[i])/(2+A0_Y1_ele[i]);
    lamda_Y2_ele[i] = (2.-3*A0_Y2_ele[i])/(2+A0_Y2_ele[i]);
    
    Errlamda_Y1_ele[i] = 8.*ErrA0_Y1_ele[i]/TMath::Power(2+A0_Y1_ele[i],2);
    Errlamda_Y2_ele[i] = 8.*ErrA0_Y2_ele[i]/TMath::Power(2+A0_Y2_ele[i],2);

    ///------- for nu -----//
    nu_Y1_mu[i] = (2.*A2_Y1_mu[i])/(2+A0_Y1_mu[i]);
    nu_Y2_mu[i] = (2.*A2_Y2_mu[i])/(2+A0_Y2_mu[i]);

    Errnu_Y1_mu[i]  = nu_Y1_mu[i]*(TMath::Power(ErrA2_Y1_mu[i]/A2_Y1_mu[i],2) + TMath::Power(ErrA0_Y1_mu[i]/(2+A0_Y1_mu[i]),2));
    Errnu_Y2_mu[i]  = nu_Y2_mu[i]*(TMath::Power(ErrA2_Y2_mu[i]/A2_Y2_mu[i],2) + TMath::Power(ErrA0_Y2_mu[i]/(2+A0_Y2_mu[i]),2));

    nu_Y1_ele[i] = (2.*A2_Y1_ele[i])/(2+A0_Y1_ele[i]);
    nu_Y2_ele[i] = (2.*A2_Y2_ele[i])/(2+A0_Y2_ele[i]);

    Errnu_Y1_ele[i]  = nu_Y1_ele[i]*(TMath::Power(ErrA2_Y1_ele[i]/A2_Y1_ele[i],2) + TMath::Power(ErrA0_Y1_ele[i]/(2+A0_Y1_ele[i]),2));
    Errnu_Y2_ele[i]  = nu_Y2_ele[i]*(TMath::Power(ErrA2_Y2_ele[i]/A2_Y2_ele[i],2) + TMath::Power(ErrA0_Y2_ele[i]/(2+A0_Y2_ele[i]),2));

    //lamtung relation
    A0mA2_Y1_ele[i] = A0_Y1_ele[i]-A2_Y1_ele[i];
    A0mA2_Y2_ele[i] = A0_Y2_ele[i]-A2_Y2_ele[i];
    A0mA2_Y2_mu[i] = A0_Y2_mu[i]-A2_Y2_mu[i];

    ErrA0mA2_Y1_ele[i] = ErrA0_Y1_ele[i]+ErrA2_Y1_ele[i];
    ErrA0mA2_Y2_ele[i] = ErrA0_Y2_ele[i]+ErrA2_Y2_ele[i];
    ErrA0mA2_Y2_mu[i] = ErrA0_Y2_mu[i]+ErrA2_Y2_mu[i];

  }


  //bins for graphs
  Float_t Zpt[nbin];
  Float_t Zpt_err[nbin];


  for (int i = 0; i < nbin; i++) {
    if (i == nbin - 1) {
      Zpt[i] = 350.;
      Zpt_err[i] = 350. - Zpt_bins[i];
    }
    else {
      Zpt[i] = ( Zpt_bins[i+1] + Zpt_bins[i] ) / 2;
      Zpt_err[i] = Zpt_bins[i+1] - Zpt[i];
    }

  }

  //graphic for angular distributions
  TGraphErrors *gr_A0_Y1_mu = new TGraphErrors(nbin, Zpt, A0_Y1_mu, Zpt_err, ErrA0_Y1_mu);
  TGraphErrors *gr_A0_Y2_mu = new TGraphErrors(nbin, Zpt, A0_Y2_mu, Zpt_err, ErrA0_Y2_mu);

  TGraphErrors *gr_A2_Y1_mu = new TGraphErrors(nbin, Zpt, A2_Y1_mu, Zpt_err, ErrA2_Y1_mu);
  TGraphErrors *gr_A2_Y2_mu = new TGraphErrors(nbin, Zpt, A2_Y2_mu, Zpt_err, ErrA2_Y2_mu);


  TGraphErrors *gr_lamda_Y1_mu = new TGraphErrors(nbin, Zpt, lamda_Y1_mu, Zpt_err, Errlamda_Y1_mu);
  TGraphErrors *gr_lamda_Y2_mu = new TGraphErrors(nbin, Zpt, lamda_Y2_mu, Zpt_err, Errlamda_Y2_mu);

  TGraphErrors *gr_nu_Y1_mu = new TGraphErrors(nbin, Zpt, nu_Y1_mu, Zpt_err, Errnu_Y1_mu);
  TGraphErrors *gr_nu_Y2_mu = new TGraphErrors(nbin, Zpt, nu_Y2_mu, Zpt_err, Errnu_Y2_mu);

  //------- ele ----------//
  TGraphErrors *gr_A0_Y1_ele = new TGraphErrors(nbin, Zpt, A0_Y1_ele, Zpt_err, ErrA0_Y1_ele);
  TGraphErrors *gr_A0_Y2_ele = new TGraphErrors(nbin, Zpt, A0_Y2_ele, Zpt_err, ErrA0_Y2_ele);

  TGraphErrors *gr_A2_Y1_ele = new TGraphErrors(nbin, Zpt, A2_Y1_ele, Zpt_err, ErrA2_Y1_ele);
  TGraphErrors *gr_A2_Y2_ele = new TGraphErrors(nbin, Zpt, A2_Y2_ele, Zpt_err, ErrA2_Y2_ele);

  TGraphErrors *gr_lamda_Y1_ele = new TGraphErrors(nbin, Zpt, lamda_Y1_ele, Zpt_err, Errlamda_Y1_ele);
  TGraphErrors *gr_lamda_Y2_ele = new TGraphErrors(nbin, Zpt, lamda_Y2_ele, Zpt_err, Errlamda_Y2_ele);

  TGraphErrors *gr_nu_Y1_ele = new TGraphErrors(nbin, Zpt, nu_Y1_ele, Zpt_err, Errnu_Y1_ele);
  TGraphErrors *gr_nu_Y2_ele = new TGraphErrors(nbin, Zpt, nu_Y2_ele, Zpt_err, Errnu_Y2_ele);

  TGraphErrors *gr_A0mA2_Y1_ele = new TGraphErrors(nbin, Zpt, A0mA2_Y1_ele, Zpt_err, ErrA0mA2_Y1_ele);
  TGraphErrors *gr_A0mA2_Y2_ele = new TGraphErrors(nbin, Zpt, A0mA2_Y2_ele, Zpt_err, ErrA0mA2_Y2_ele);
  TGraphErrors *gr_A0mA2_Y2_mu = new TGraphErrors(nbin, Zpt, A0mA2_Y2_mu, Zpt_err, ErrA0mA2_Y2_mu);

  gr_A0_Y1_mu->SetMarkerStyle(20);
  gr_A0_Y1_mu->SetMarkerColor(kGreen+2);
  gr_A0_Y1_mu->SetLineColor(kGreen+2);

  gr_A0_Y2_mu->SetMarkerStyle(21);
  gr_A0_Y2_mu->SetMarkerColor(kRed-7);
  gr_A0_Y2_mu->SetLineColor(kRed-7);

  gr_A2_Y1_mu->SetMarkerStyle(20);
  gr_A2_Y1_mu->SetMarkerColor(kGreen+2);
  gr_A2_Y1_mu->SetLineColor(kGreen+2);

  gr_A2_Y2_mu->SetMarkerStyle(21);
  gr_A2_Y2_mu->SetMarkerColor(kRed-7);
  gr_A2_Y2_mu->SetLineColor(kRed-7);

  gr_lamda_Y1_mu->SetMarkerStyle(20);
  gr_lamda_Y1_mu->SetMarkerColor(kGreen+2);
  gr_lamda_Y1_mu->SetLineColor(kGreen+2);

  gr_lamda_Y2_mu->SetMarkerStyle(21);
  gr_lamda_Y2_mu->SetMarkerColor(kRed-7);
  gr_lamda_Y2_mu->SetLineColor(kRed-7);

  gr_nu_Y1_mu->SetMarkerStyle(20);
  gr_nu_Y1_mu->SetMarkerColor(kGreen+2);
  gr_nu_Y1_mu->SetLineColor(kGreen+2);

  gr_nu_Y2_mu->SetMarkerStyle(21);
  gr_nu_Y2_mu->SetMarkerColor(kRed-7);
  gr_nu_Y2_mu->SetLineColor(kRed-7);

  gr_A0mA2_Y2_mu->SetMarkerStyle(21);
  gr_A0mA2_Y2_mu->SetMarkerColor(kRed-7);
  gr_A0mA2_Y2_mu->SetLineColor(kRed-7);

  //style for ele
  gr_A0_Y1_ele->SetMarkerStyle(20);
  gr_A0_Y1_ele->SetMarkerColor(kGreen+2);
  gr_A0_Y1_ele->SetLineColor(kGreen+2);

  gr_A0_Y2_ele->SetMarkerStyle(21);
  gr_A0_Y2_ele->SetMarkerColor(kRed-7);
  gr_A0_Y2_ele->SetLineColor(kRed-7);

  gr_A2_Y1_ele->SetMarkerStyle(20);
  gr_A2_Y1_ele->SetMarkerColor(kGreen+2);
  gr_A2_Y1_ele->SetLineColor(kGreen+2);

  gr_A2_Y2_ele->SetMarkerStyle(21);
  gr_A2_Y2_ele->SetMarkerColor(kRed-7);
  gr_A2_Y2_ele->SetLineColor(kRed-7);

  gr_lamda_Y1_ele->SetMarkerStyle(20);
  gr_lamda_Y1_ele->SetMarkerColor(kGreen+2);
  gr_lamda_Y1_ele->SetLineColor(kGreen+2);

  gr_lamda_Y2_ele->SetMarkerStyle(21);
  gr_lamda_Y2_ele->SetMarkerColor(kRed-7);
  gr_lamda_Y2_ele->SetLineColor(kRed-7);

  gr_nu_Y1_ele->SetMarkerStyle(20);
  gr_nu_Y1_ele->SetMarkerColor(kGreen+2);
  gr_nu_Y1_ele->SetLineColor(kGreen+2);

  gr_nu_Y2_ele->SetMarkerStyle(21);
  gr_nu_Y2_ele->SetMarkerColor(kRed-7);
  gr_nu_Y2_ele->SetLineColor(kRed-7);


  gr_A0mA2_Y1_ele->SetMarkerStyle(20);
  gr_A0mA2_Y1_ele->SetMarkerColor(kGreen+2);
  gr_A0mA2_Y1_ele->SetLineColor(kGreen+2);

  gr_A0mA2_Y2_ele->SetMarkerStyle(21);
  gr_A0mA2_Y2_ele->SetMarkerColor(kRed-7);
  gr_A0mA2_Y2_ele->SetLineColor(kRed-7);

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.05);
  tx.SetTextFont(42);

  TCanvas *cA0 = new TCanvas("cA0", "cA0", 650, 650);
  cA0->cd();
  cA0->SetLogx();

  gr_A0_Y2_mu->GetYaxis()->SetTitle("A_{0}");
  gr_A0_Y2_mu->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  gr_A0_Y2_mu->GetXaxis()->SetTitleOffset(1.3);
  gr_A0_Y2_mu->GetYaxis()->SetRangeUser(-0.1, 1.2);
  gr_A0_Y2_mu->Draw("APZ");
  //gr_A0_Y1_mu->Draw("PZ");

  TLegend *lg = new TLegend(0.2, 0.4, 0.5, 0.54);
  lg->AddEntry(gr_A0_Y1_mu, "|Y^{Z}| < 1", "pe");
  lg->AddEntry(gr_A0_Y2_mu, "1 < |Y^{Z}| < 2.0", "pe");
  lg->SetBorderSize(0);
  lg->SetTextSize(0.05);
  lg->Draw();
  tx.DrawLatex(0.3, 0.6, "Z#rightarrow #mu#mu");



  TCanvas *cA2 = new TCanvas("cA2", "cA2", 650, 650);
  cA2->cd();
  cA2->SetLogx();

  gr_A2_Y2_mu->GetYaxis()->SetTitle("A_{2}");
  gr_A2_Y2_mu->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  //gr_A2_Y2_mu->SetTitleFont(42, "XYZ");
  //gr_A2_Y2_mu->SetTitleSize(0.042, "XYZ");
  gr_A2_Y2_mu->GetXaxis()->SetTitleOffset(1.3);
  gr_A2_Y2_mu->GetYaxis()->SetRangeUser(-0.1, 0.5);
  gr_A2_Y2_mu->Draw("APZ");
  //gr_A2_Y1_mu->Draw("PZ");
  lg->Draw();
  tx.DrawLatex(0.3, 0.6, "Z#rightarrow #mu#mu");


  TCanvas *clamda = new TCanvas("clamda", "clamda", 650, 650);
  clamda->cd();
  clamda->SetLogx();

  gr_lamda_Y2_mu->GetYaxis()->SetTitle("#lambda");
  gr_lamda_Y2_mu->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  //gr_lamda_Y2_mu->SetTitleFont(42, "XYZ");
  //gr_lamda_Y2_mu->SetTitleSize(0.042, "XYZ");
  gr_lamda_Y2_mu->GetXaxis()->SetTitleOffset(1.3);
  gr_lamda_Y2_mu->GetYaxis()->SetRangeUser(-0.4, 1.2);
  gr_lamda_Y2_mu->Draw("APZ");
  //gr_lamda_Y1_mu->Draw("PZ");
  lg->Draw();
  tx.DrawLatex(0.3, 0.6, "Z#rightarrow #mu#mu");

  TCanvas *cA0mA2 = new TCanvas("cA0mA2", "cA0mA2", 650, 650);
  cA0mA2->cd();
  cA0mA2->SetLogx();

  gr_A0mA2_Y2_mu->GetYaxis()->SetTitle("A_{0}-A_{2}");
  gr_A0mA2_Y2_mu->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  //gr_A0mA2_Y2_mu->SetTitleFont(42, "XYZ");
  //gr_A0mA2_Y2_mu->SetTitleSize(0.042, "XYZ");
  gr_A0mA2_Y2_mu->GetXaxis()->SetTitleOffset(1.3);
  gr_A0mA2_Y2_mu->GetYaxis()->SetRangeUser(-0.1, 1.);
  gr_A0mA2_Y2_mu->Draw("APZ");
  //gr_A0mA2_Y1_mu->Draw("PZ");
  lg->Draw();
  tx.DrawLatex(0.3, 0.6, "Z#rightarrow #mu#mu");

  TCanvas *cnu = new TCanvas("cnu", "cnu", 650, 650);
  cnu->cd();
  cnu->SetLogx();

  gr_nu_Y2_mu->GetYaxis()->SetTitle("#nu");
  gr_nu_Y2_mu->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  //gr_nu_Y2_mu->SetTitleFont(42, "XYZ");
  //gr_nu_Y2_mu->SetTitleSize(0.042, "XYZ");
  gr_nu_Y2_mu->GetXaxis()->SetTitleOffset(1.3);
  gr_nu_Y2_mu->GetYaxis()->SetRangeUser(-0.2, 1.2);
  gr_nu_Y2_mu->Draw("APZ");
  //gr_nu_Y1_mu->Draw("PZ");
  lg->Draw();
  tx.DrawLatex(0.3, 0.6, "Z#rightarrow #mu#mu");


  TCanvas *cA0el = new TCanvas("cA0el", "cA0el", 650, 650);
  cA0el->cd();
  cA0el->SetLogx();

  gr_A0_Y1_ele->GetYaxis()->SetTitle("A_{0}");
  gr_A0_Y1_ele->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  gr_A0_Y1_ele->GetXaxis()->SetTitleOffset(1.3);
  gr_A0_Y1_ele->GetYaxis()->SetRangeUser(-0.1, 1.2);
  gr_A0_Y1_ele->Draw("APZ");
  gr_A0_Y2_ele->Draw("PZ");
  lg->Draw();
  tx.DrawLatex(0.3, 0.6, "Z#rightarrow ee");

  TCanvas *cA2el = new TCanvas("cA2el", "cA2el", 650, 650);
  cA2el->cd();
  cA2el->SetLogx();

  gr_A2_Y1_ele->GetYaxis()->SetTitle("A_{2}");
  gr_A2_Y1_ele->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  gr_A2_Y1_ele->GetXaxis()->SetTitleOffset(1.3);
  gr_A2_Y1_ele->GetYaxis()->SetRangeUser(-0.2, 1.2);
  gr_A2_Y1_ele->Draw("APZ");
  gr_A2_Y2_ele->Draw("PZ");
  lg->Draw();
  tx.DrawLatex(0.3, 0.6, "Z#rightarrow ee");

  TCanvas *clamdael = new TCanvas("clamdael", "clamdael", 650, 650);
  clamdael->cd();
  clamdael->SetLogx();

  gr_lamda_Y1_ele->GetYaxis()->SetTitle("#lambda");
  gr_lamda_Y1_ele->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  gr_lamda_Y1_ele->GetXaxis()->SetTitleOffset(1.3);
  gr_lamda_Y1_ele->GetYaxis()->SetRangeUser(-0.4, 1.2);
  gr_lamda_Y1_ele->Draw("APZ");
  gr_lamda_Y2_ele->Draw("PZ");
  lg->Draw();
  tx.DrawLatex(0.3, 0.6, "Z#rightarrow ee");


  TCanvas *cA0mA2el = new TCanvas("cA0mA2el", "cA0mA2el", 650, 650);
  cA0mA2el->cd();
  cA0mA2el->SetLogx();

  gr_A0mA2_Y1_ele->GetYaxis()->SetTitle("A_{0}-A_{2}");
  gr_A0mA2_Y1_ele->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  gr_A0mA2_Y1_ele->GetXaxis()->SetTitleOffset(1.3);
  gr_A0mA2_Y1_ele->GetYaxis()->SetRangeUser(-0.1, 0.8);
  gr_A0mA2_Y1_ele->Draw("APZ");
  gr_A0mA2_Y2_ele->Draw("PZ");
  lg->Draw();
  tx.DrawLatex(0.3, 0.6, "Z#rightarrow ee");

  TCanvas *cnuel = new TCanvas("cnuel", "cnuel", 650, 650);
  cnuel->cd();
  cnuel->SetLogx();

  gr_nu_Y1_ele->GetYaxis()->SetTitle("#nu");
  gr_nu_Y1_ele->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  gr_nu_Y1_ele->GetXaxis()->SetTitleOffset(1.3);
  gr_nu_Y1_ele->GetYaxis()->SetRangeUser(-0.2, 1.2);
  gr_nu_Y1_ele->Draw("APZ");
  gr_nu_Y2_ele->Draw("PZ");
  lg->Draw();
  tx.DrawLatex(0.3, 0.6, "Z#rightarrow ee");


  cA0->SaveAs("plots/angularCoeff/A0_mu_2016_reboot.pdf");
  cA2->SaveAs("plots/angularCoeff/A2_mu_2016_reboot.pdf");
  cnu->SaveAs("plots/angularCoeff/Nu_mu_2016_reboot.pdf");
  cA0el->SaveAs("plots/angularCoeff/A0_ele_2016_reboot.pdf");
  cA2el->SaveAs("plots/angularCoeff/A2_ele_2016_reboot.pdf");
  clamdael->SaveAs("plots/angularCoeff/lamda_ele_2016_reboot.pdf");
  cnuel->SaveAs("plots/angularCoeff/nu_ele_2016_reboot.pdf");
  cA0mA2->SaveAs("plots/angularCoeff/A0mA2_mu_2016_reboot.pdf");
  cA0mA2el->SaveAs("plots/angularCoeff/A0mA2_ele_2016_reboot.pdf");



}
