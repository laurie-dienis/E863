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

void Pipe3_second_Ecm(){

  // List of file names
  std::vector<std::string> fileNames = {
      "RootFiles/PID/alpha_2.root",
      //"RootFiles/Raw/run_test2.root",
  };
  const char *treeName = "DataTree"; // Tree name (common to all files)

  // Variables to hold leaf data
  double si;
  double si_cm;

  TH1F *hist_E = new TH1F("hist_si_cal", "SI Calibrated Histogram", 2000, 1, 10);
  TH1F *hist_Ecm = new TH1F("hist_si_cal_Ecm", "SI Calibrated Histogram", 2000, 10, 20);

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
    tree->SetBranchAddress("Si_cal_2", &si);

    // Loop over the tree entries
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
      tree->GetEntry(i);
      hist_E->Fill(si);
      si_cm = si * 1. ;
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
