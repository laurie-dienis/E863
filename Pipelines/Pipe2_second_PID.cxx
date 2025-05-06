#ifndef Pipe2_PID_cxx
#define Pipe2_PID_cxx

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"

#include <fstream>
#include <iostream>
#include <string>

  // Create histograms
  TH2F* PID = new TH2F("PID", "PID", 800, 0, 20, 1100, 0, 110);
  TH2F* PID_faster =
  new TH2F("PID_faster", "PID_faster", 600, 200, 800, 2000, 0, 20);

  void Pipe2_second_PID(unsigned int runN) {
 TString runNumber = Form("%i",runN);
 
 // Check if the cut file exists
  const char *cutFileName = "Cuts/cut_alpha_2.root";
  TFile *cutFile = TFile::Open(cutFileName);
  TCutG *cutG = nullptr;

  if (cutFile && !cutFile->IsZombie()) {
    // Try to read the graphical cut
    cutFile->GetObject("CUTG", cutG);
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
 const char* treeName = "DataTree"; // Tree name (common to all files)
 TString fileName = "RootFiles/Cal/15O/Cal" + runNumber + ".root";
 TFile* file = TFile::Open(fileName);
 
 // Get the tree
 TTree* tree = nullptr;
 file->GetObject(treeName, tree);
 
 if (!tree) {
   std::cerr << "Error: Tree " << treeName
             << " not found in file: " << fileName << std::endl;
   file->Close();
   delete file;
   exit(1);
 }
 
 // Variables to hold leaf data
 double tac;
 double tac_diff;
 double si_cal;
 double faster_time;
 
 // Create a new ROOT file only if the cut exists
 TTree* PIDTree = nullptr;
 
 TFile* newFile = NULL;
 if (cutG) {
   // Create a new ROOT file to save the new tree and histograms
   newFile = TFile::Open("RootFiles/PID/15O/alpha_2_" + runNumber + ".root", "RECREATE");
   if (!newFile || newFile->IsZombie()) {
     std::cerr << "Error creating new file: RootFiles/Processed/output.root"
               << std::endl;
     return;
   }
  }
 
  // Create a new tree
  PIDTree = new TTree(treeName, "alpha Data Tree");
 
  // Set branch addresses for the new tree
  PIDTree->Branch("Si_cal_2", &si_cal, "Si_cal_2/D");
 
  // Set branch addresses for reading
  tree->SetBranchAddress("Si_cal_2", &si_cal);
  tree->SetBranchAddress("TAC_cal_2", &tac);
  tree->SetBranchAddress("faster_time_2", &faster_time);

  // Loop over the tree entries and apply the graphical cut if it exists
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    tac_diff = 100 - tac;
    // Always fill the PID histogram
    if(si_cal > 0 && tac > 0){
//     if (cutG && cutG->IsInside(si_cal, tac_diff)) {
       PID->Fill(si_cal, tac_diff);
      //  PIDTree->Fill();
//     } 
      // PID_faster->Fill(si_cal, faster_time);
    }
 
    // Fill the new tree only if the cut exists and the point is inside the
    // cut
    if (cutG && cutG->IsInside(si_cal, tac_diff)) {
      PIDTree->Fill();
    }
  }
 
  PIDTree->Write();
  // Plot the PID histogram
  auto *c20 = new TCanvas("c20", "Pipe2 canvas 0");
  PID->GetYaxis()->SetTitle("TOF (ns)");
  PID->GetXaxis()->SetTitle("E_Si (MeV)");
  PID->Draw("colz");
  std::cout << PID->GetEntries() << std::endl;

  // alpha
  TF1 *line2 = new TF1("line2", "(20.6504*x^(-0.499735))", 1, 7);
  line2->SetLineColor(kViolet);
  line2->SetLineWidth(4);
  line2->Draw("same");
 
  ////////Ti///////////////
  // //protons
  TF1 *line3 = new TF1("line3", "41.1506*x^(-0.499787)", 7, 13);
  line3->SetLineColor(kOrange);
  line3->SetLineWidth(4);
  line3->Draw("same");
 
  // Plot the cut if it exists
  if (cutG) {
    cutG->SetLineColor(kRed); // Set cut line color
    cutG->SetLineWidth(2);    // Set cut line width
    cutG->Draw("L same");     // Draw cut on top of the histogram
  }
  
  // auto *c21 = new TCanvas("c21", "Pipe2 canvas 1");
  // PID_faster->GetYaxis()->SetTitle("faster time (ns)");
  // PID_faster->GetXaxis()->SetTitle("E_Si (MeV)");
  // PID_faster->Draw("colz");

 // Close the input file and clean up
 file->Close();
 delete file;
}
#endif
