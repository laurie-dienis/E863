#ifndef Pipe2_PID_cxx
#define Pipe2_PID_cxx

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TROOT.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <string>

void Pipe2_second_PID() {

  // List of file names
  std::vector<std::string> fileNames = {
      "RootFiles/Cal/Cal.root",
      // Add more files as needed
  };
  const char *treeName = "DataTree"; // Tree name (common to all files)

  // Variables to hold leaf data
  double tac;
  double si_cal;

  // Create histograms
  TH2F *PID = new TH2F("PID", "PID", 3000, 0, 100, 3000, 0, 10);

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
      delete file;
      continue;
    }

    // Set branch addresses for reading
    tree->SetBranchAddress("Si_cal_b", &si_cal);
    tree->SetBranchAddress("TAC_cal_b", &tac);

    // Loop over the tree entries and apply the graphical cut if it exists
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
      tree->GetEntry(i);

      // Always fill the PID histogram
      PID->Fill(tac, si_cal);

      // Fill the new tree only if the cut exists and the point is inside the
      // cut
      if (cutG && cutG->IsInside(tac, si_cal)) {
        PIDTree->Fill();
      }
    }

    // Close the input file and clean up
    file->Close();
    delete file;
  }

  // Plot the PID histogram
  auto *c20 = new TCanvas("c20", "Pipe2 canvas 0");
  PID->GetXaxis()->SetTitle("TOF (s)");
  PID->GetYaxis()->SetTitle("E_Si (MeV)");
  PID->Draw("colz");
}
#endif
