#include "xAna.C"
void runMC(bool, int);
void runData(bool, int);

void runxAna(){
  gSystem->SetBuildDir("tmpdir", kTRUE);
  gROOT->ProcessLine("xAna.C+");

  int iWP = 0;

  runData(0, iWP);
  //runMC(1, iWP);

}

void runData(bool mc = 0, int iWP = 0) {

  /*
  TString sample[7] = {
    "/data1/ggNtuples/V09_04_13_03/job_DoubleEG_Run2016B_Legacy/*.root",
    "/data1/ggNtuples/V09_04_13_03/job_DoubleEG_Run2016C_Legacy/*.root",
    "/data1/ggNtuples/V09_04_13_03/job_DoubleEG_Run2016D_Legacy/*.root",
    "/data1/ggNtuples/V09_04_13_03/job_DoubleEG_Run2016E_Legacy/*.root",
    "/data1/ggNtuples/V09_04_13_03/job_DoubleEG_Run2016F_Legacy/*.root",
    "/data1/ggNtuples/V09_04_13_03/job_DoubleEG_Run2016G_Legacy/*.root",
    "/data1/ggNtuples/V09_04_13_03/job_DoubleEG_Run2016H_Legacy/*.root",
  };
  */
  /*
  TString sample[7] = {
    "/data1/ggNtuples/V09_04_13_03/job_DoubleMu_Run2016B_Legacy/*.root",
    "/data1/ggNtuples/V09_04_13_03/job_DoubleMu_Run2016C_Legacy/*.root",
    "/data1/ggNtuples/V09_04_13_03/job_DoubleMu_Run2016D_Legacy/*.root",
    "/data1/ggNtuples/V09_04_13_03/job_DoubleMu_Run2016E_Legacy/*.root",
    "/data1/ggNtuples/V09_04_13_03/job_DoubleMu_Run2016F_Legacy/*.root",
    "/data1/ggNtuples/V09_04_13_03/job_DoubleMu_Run2016G_Legacy/*.root",
    "/data1/ggNtuples/V09_04_13_03/job_DoubleMu_Run2016H_Legacy/*.root"

  };
  */

  /*
  for (Int_t i=0; i<7; ++i) {
    TString inpaths = sample[i].Data();
    cout << "Process sample: " << inpaths << endl;
    xAna(inpaths, "out.root",  mc, iWP);
  }
  */

  //TString sample = "/data1/ggNtuples/V09_04_13_03/job_DoubleEG_*Legacy/*.root";
  TString sample = "/data1/ggNtuples/V09_04_13_03/job_DoubleMu_Run2016*Legacy/*.root";
  

  //TString inpaths = sample[i].Data();
  cout << "Process sample: " << sample << endl;
  xAna(sample, "out.root",  mc, iWP);

}

void runMC(bool mc = 1, int iWP = 0) {

  TString sample[11] = {
    "/data6/ggNtuples/V09_04_13_04/job_summer16_DYJetsToLL_m50_aMCatNLO_ext2/ggtree*.root",
    "/data6/ggNtuples/V09_04_13_04/job_summer16_WJets_madgraph_ext2/ggtree*.root",
    "/data6/ggNtuples/V09_04_13_04/job_summer16_TT_powheg/ggtree*.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_WWTo2L2Nu.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_WWToLNuQQ*.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_WZTo3LNu.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_WZTo2L2Q.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_ZZTo2L2Nu*.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_ZZTo2L2Q.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_ZZTo4L.root",
    "/data3/ggNtuples/Zg_summer16/TTGamma_Summer16_amcnlo_skim.root"
  };

  
  for (Int_t i=0; i<1; ++i) {
    TString inpaths = sample[i].Data();
    cout << "Process sample: " << inpaths << endl;
    xAna(inpaths, "out.root", mc, iWP);
  }

}

