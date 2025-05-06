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

void Pipe2_main_PID_all() {
 
 TChain *chain = new TChain("DataTree");  
 chain->Add("RootFiles/Cal/15N/Cal*.root");
 
 // Variables to hold leaf data
 double tac;
 double tac_diff;
 double si_cal;
 double faster_time;
 
  // Set branch addresses for reading
  chain->SetBranchAddress("Si_cal_1", &si_cal);
  chain->SetBranchAddress("TAC_cal_1", &tac);
  chain->SetBranchAddress("faster_time_1", &faster_time);

  // Loop over the chain entries and apply the graphical cut if it exists
  Long64_t nEntries = chain->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    chain->GetEntry(i);
    tac_diff = 100 - tac;
    PID->Fill(si_cal, tac_diff);
      // PID_faster->Fill(si_cal, faster_time);
    }
 
  // Plot the PID histogram
  auto *c20 = new TCanvas("c20", "Pipe2 canvas 0");
  PID->GetYaxis()->SetTitle("TOF (ns)");
  PID->GetXaxis()->SetTitle("E_Si (MeV)");
  PID->Draw("colz");
  std::cout << PID->GetEntries() << std::endl;

  // alpha
  // TF1 *line2 = new TF1("line2", "(20.6504*x^(-0.499735))+60", 1, 7);
  // line2->SetLineColor(kViolet);
  // line2->SetLineWidth(4);
  // line2->Draw("same");
 
  // ////////Ti///////////////
  // // //protons
  // TF1 *line3 = new TF1("line3", "41.1506*x^(-0.499787)+60", 7, 13);
  // line3->SetLineColor(kOrange);
  // line3->SetLineWidth(4);
  // line3->Draw("same");
  
  // auto *c21 = new TCanvas("c21", "Pipe2 canvas 1");
  // PID_faster->GetYaxis()->SetTitle("faster time (ns)");
  // PID_faster->GetXaxis()->SetTitle("E_Si (MeV)");
  // PID_faster->Draw("colz");

}
#endif
