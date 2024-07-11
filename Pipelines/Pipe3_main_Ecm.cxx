#ifndef Pipe3_Ecm_cxx
#define Pipe3_Ecm_cxx

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TROOT.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <string>

void Pipe3_main_Ecm(){

  // List of file names
  std::vector<std::string> fileNames = {
      "RootFiles/PID/alpha_1.root",
      //"RootFiles/Raw/run_test2.root",
  };
  const char *treeName = "DataTree"; // Tree name (common to all files)

  // Variables to hold leaf data
  double si;
  double si_cm;
  //parameters to compute ecm from emeas, obtained from polu_Eameas_Eacm.C in Simulation directory
  double p0 = 0.000772;  //0.000994
  double p1 = 0.277615;  //0.276885
  double p2 = -0.000192; //-0.000083

  TH1F *hist_E = new TH1F("hist_si_cal", "SI Calibrated Histogram", 500, 1, 10);
  TH1F *hist_Ecm = new TH1F("hist_si_cal_Ecm", "SI Calibrated Histogram", 500, 0, 6);

  // Loop over all files
  for (const auto &fileName : fileNames) {
    // Open the ROOT file
    TFile *file = TFile::Open(fileName.c_str());
    if (!file || file->IsZombie()) {
      std::cerr << "Error opening file: " << fileName << std::endl;
      continue;
    }

    // Get the tree
    TTree *tree = nullptr;
    file->GetObject(treeName, tree);
    if (!tree) {
      std::cerr << "Error: Tree " << treeName
                << " not found in file: " << fileName << std::endl;
      file->Close();
      continue;
    }

    // Set branch addresses
    tree->SetBranchAddress("Si_cal_1", &si);

    // Loop over the tree entries
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
      tree->GetEntry(i);
      hist_E->Fill(si);
      si_cm = p0+ p1* si+ p2 *si*si  ;
      hist_Ecm->Fill(si_cm);
    }

    // Close the file
    file->Close();
    delete file;
  }


  // Draw histograms
  TCanvas *c30 = new TCanvas("c30", "SI Histogram", 800, 600);
  hist_E->Draw();
  c30->SaveAs("hist_si_alpha.png");
  TCanvas *c31 = new TCanvas("c31", "SI Histogram Ecm", 800, 600);
  hist_Ecm->Draw();
  c31->SaveAs("hist_si_alpha_Ecm.png");
}
#endif
