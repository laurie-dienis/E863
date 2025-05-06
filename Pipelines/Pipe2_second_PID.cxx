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
  // Check if the cut file exists
  const char *cutFileName = "Cuts/cut_alpha_2.root";
  TFile *cutFile = TFile::Open(cutFileName);
  TCutG *cutG = nullptr;

  if (cutFile && !cutFile->IsZombie()) {
    // Try to read the graphical cut
    cutFile->GetObject("alpha_cut", cutG);
    cutFile->Close();
    delete cutFile;
  }

  // Test if the graphical cut is found
  if (!cutG) {
    std::cerr
        << "Error: Graphical cut 'cutG' not found in file 'cut_alpha_2.root'."
        << " Plotting PID without cut." << std::endl;
  }

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
  TH2F *PID = new TH2F("PID", "PID", 300, 0, 20, 300, 0, 109);

  // Create a new ROOT file only if the cut exists
  TFile *newFile = nullptr;
  TTree *PIDTree = nullptr;

  if (cutG) {
    newFile = TFile::Open("RootFiles/PID/alpha_2.root", "RECREATE");
    if (!newFile || newFile->IsZombie()) {
      std::cerr << "Error creating new file: RootFiles/PID/alpha_2.root"
                << std::endl;
      if (newFile)
        newFile->Close();
      delete newFile;
      delete PID;
      return;
    }

    // Create a new tree
    PIDTree = new TTree(treeName, "alpha Data Tree");

    // Set branch addresses for the new tree
    PIDTree->Branch("Si_cal_2", &si_cal, "Si_cal_2/D");
  }

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
    tree->SetBranchAddress("Si_cal_2", &si_cal);
    tree->SetBranchAddress("TAC_cal_2", &tac);

    // Loop over the tree entries and apply the graphical cut if it exists
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
      tree->GetEntry(i);

      // Always fill the PID histogram
      PID->Fill(si_cal, tac);

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
  PID->GetYaxis()->SetTitle("TOF (ns)");
  PID->GetXaxis()->SetTitle("E_Si (MeV)");
  PID->Draw("colz");

  // Plot the cut if it exists
  if (cutG) {
    cutG->SetLineColor(kRed); // Set cut line color
    cutG->SetLineWidth(2);    // Set cut line width
    cutG->Draw("L same");     // Draw cut on top of the histogram
  }
  c20->SaveAs("PID_2.png");

  // Write the histogram and the new tree (if the cut exists) to the new file
  if (newFile) {
    newFile->cd(); // Ensure the new file is the current directory
    PID->Write();
    PIDTree->Write();

    // Close the new file
    newFile->Close();
    delete newFile;
  }
}
#endif
