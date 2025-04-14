#ifndef Pipe1_Cal_cxx
#define Pipe1_Cal_cxx

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TROOT.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <string>

void Pipe1_Cal() {

  // List of file names
  std::vector<std::string> fileNames = {
      "faster-to-root/root/calib_Si_4.root",
  };

  //   std::vector<std::string> fileNames;
  //   std::ifstream runFile("runs.txt");
  //   std::string runNumber;

  //   if (!runFile) {
  //       std::cerr << "Error opening run.txt file." << std::endl;
  //   }
  //   while (std::getline(runFile, runNumber)) {
  //     // Ignorer les lignes vides
  //     if (runNumber.empty()) continue;

  //     std::string path = "faster-to-root/root/calib_Si_" + runNumber +
  //     ".root"; fileNames.push_back(path);
  // }

  const char *treeName = "DataTree"; // Tree name (common to all files)

  // Variables to hold leaf data
  double si_1;
  double tac_1;
  double si_cal_1;
  double tac_cal_1;
  double faster_time_1;
  double cal_parameter_Si_1 = 0.0920178;     // new parameter
  double cal_parameter_tac_1 = 0.0012941633; // new parameter

  double si_2;
  double tac_2;
  double si_cal_2;
  double tac_cal_2;
  double cal_parameter_Si_2 = 0.0736969;
  double cal_parameter_tac_2 = 0.0011900512; // new parameter

  double si_b;
  double tac_b;
  double si_cal_b;
  double tac_cal_b;
  double cal_parameter_Si_b = 0.389885;      // new parameter
  double cal_parameter_tac_b = 0.0011870845; // new parameter

  // Create histograms
  // TH1F *hist_si_1 = new TH1F("hist_si_1", "SI Histogram 1", 10000, 0,
  // 220000);
  TH1D *hist_si_1 = new TH1D("hist_si_1", "SI Histogram 1", 500, 53000, 65000);
  TH1F *hist_tac_1 = new TH1F("hist_tac_1", "TAC Histogram 1", 7000, 0, 850000);
  // TH1F *hist_si_cal_1 = new TH1F("hist_si_cal_1", "SI Calibrated Histogram
  // 1", 2000, 1, 22);
  TH1F *hist_si_cal_1 =
      new TH1F("hist_si_cal_1", "SI Calibrated Histogram 1", 500, 5, 6);
  TH1F *hist_tac_cal_1 =
      new TH1F("hist_tac_cal_1", "TAC calibrated Histogram 1", 5000, 0, 80);

  // TH1F *hist_si_2 = new TH1F("hist_si_2", "SI Histogram 2", 10000, 0,
  // 220000);
  TH1D *hist_si_2 = new TH1D("hist_si_2", "SI Histogram 2", 300, 63000, 82000);
  TH1F *hist_tac_2 = new TH1F("hist_tac_2", "TAC Histogram 2", 5000, 0, 950000);
  TH1F *hist_si_cal_2 =
      new TH1F("hist_si_cal_2", "SI Calibrated Histogram 2", 2000, 1, 22);
  TH1F *hist_tac_cal_2 =
      new TH1F("hist_tac_cal_2", "TAC Calibrated Histogram 2", 5000, 0, 80);
  TH2F *hist_faster_TAC = new TH2F("hist_faster_TAC", "faster time vs TAC", 500,
                                   0, 110, 100, 200, 800);

  TH1D *hist_si_b = new TH1D("hist_si_b", "SI Histogram b", 300, 12000, 16000);
  TH1F *hist_tac_b = new TH1F("hist_tac_b", "TAC Histogram b", 5000, 0, 950000);
  TH1F *hist_si_cal_b =
      new TH1F("hist_si_cal_b", "SI Calibrated Histogram b", 2000, 1, 22);
  TH1F *hist_tac_cal_b =
      new TH1F("hist_tac_cal_b", "TAC Calibrated Histogram b", 5000, 0, 80);

  // Create a new ROOT file to save the new tree and histograms
  TFile *newFile = TFile::Open("RootFiles/Cal/Cal.root", "RECREATE");
  if (!newFile || newFile->IsZombie()) {
    std::cerr << "Error creating new file: RootFiles/Processed/output.root"
              << std::endl;
    return;
  }

  // Create a new tree
  TTree *CalTree = new TTree(treeName, "Calibrated Data Tree");

  // Set branch addresses for the new tree
  CalTree->Branch("Si_1", &si_1, "Si_1/D");
  CalTree->Branch("TAC_1", &tac_1, "TAC_1/D");
  CalTree->Branch("Si_cal_1", &si_cal_1, "Si_cal_1/D");
  CalTree->Branch("TAC_cal_1", &tac_cal_1, "TAC_cal_1/D");
  CalTree->Branch("faster_time_1", &faster_time_1, "faster_time_1/D");
  CalTree->Branch("Si_2", &si_2, "Si_2/D");
  CalTree->Branch("TAC_2", &tac_2, "TAC_2/D");
  CalTree->Branch("Si_cal_2", &si_cal_2, "Si_cal_2/D");
  CalTree->Branch("TAC_cal_2", &tac_cal_2, "TAC_cal_2/D");
  CalTree->Branch("Si_B", &si_b, "Si_B/D");
  CalTree->Branch("TAC_B", &tac_b, "TAC_B/D");
  CalTree->Branch("Si_cal_B", &si_cal_b, "Si_cal_B/D");
  CalTree->Branch("TAC_cal_B", &tac_cal_b, "TAC_cal_B/D");

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
    tree->SetBranchAddress("Si_1", &si_1);
    tree->SetBranchAddress("TAC_1", &tac_1);
    tree->SetBranchAddress("faster_time_1", &faster_time_1);
    tree->SetBranchAddress("Si_2", &si_2);
    tree->SetBranchAddress("TAC_2", &tac_2);
    tree->SetBranchAddress("Si_B", &si_b);
    tree->SetBranchAddress("TAC_B", &tac_b);

    // Loop over the tree entries
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
      tree->GetEntry(i);
      if (tac_1 > 85) {
        hist_si_1->Fill(si_1);
      }
      hist_tac_1->Fill(tac_1);
      si_cal_1 = si_1 * cal_parameter_Si_1 * 1e-3;
      if (tac_1 > 85) {
        hist_si_cal_1->Fill(si_cal_1);
      }
      tac_cal_1 = tac_1 * cal_parameter_tac_1;
      hist_tac_cal_1->Fill(tac_cal_1);

      if (tac_2 > 85) {
        hist_si_2->Fill(si_2);
      }
      hist_tac_2->Fill(tac_2);
      si_cal_2 = si_2 * cal_parameter_Si_2 * 1e-3;
      if (tac_2 > 85) {
        hist_si_cal_2->Fill(si_cal_2);
      }
      tac_cal_2 = tac_2 * cal_parameter_tac_2;
      hist_tac_cal_2->Fill(tac_cal_2);
      // hist_si_cal_12->Fill(si_cal_1,si_cal_2);
      hist_faster_TAC->Fill(tac_cal_1, faster_time_1);

      // if (tac_b > 85) {
      hist_si_b->Fill(si_b);
      // }
      hist_tac_b->Fill(tac_b);
      si_cal_b = si_b * cal_parameter_Si_b * 1e-3;
      // std::cout << "sibcal = " << si_b << "\n";
      // if (tac_b > 85) {
      hist_si_cal_b->Fill(si_cal_b);
      // }
      tac_cal_b = tac_b * cal_parameter_tac_b;
      hist_tac_cal_b->Fill(tac_cal_b);

      // Fill the new tree
      CalTree->Fill();
    }

    // Close the file
    file->Close();
    delete file;
  }

  // Write the histograms and the new tree to the new file
  newFile->cd();

  hist_si_1->Write();
  hist_tac_1->Write();
  hist_si_cal_1->Write();
  hist_tac_cal_1->Write();

  hist_si_2->Write();
  hist_tac_2->Write();
  hist_si_cal_2->Write();
  hist_tac_cal_2->Write();

  hist_si_b->Write();
  hist_tac_b->Write();
  hist_si_cal_b->Write();
  hist_tac_cal_b->Write();
  CalTree->Write();

  // Close the new file
  newFile->Close();

  // Draw histograms
  TCanvas *c10 = new TCanvas("c10", "SI Histogram", 800, 600);
  c10->DivideSquare(2);
  c10->cd(1);
  hist_si_1->Draw();
  c10->cd(2);
  hist_si_2->Draw();
  c10->SaveAs("hist_si.png");

  TCanvas *c17 = new TCanvas("c17", "SI_B Histogram", 800, 600);
  hist_si_b->GetXaxis()->SetTitle("E (arb. units)");
  hist_si_b->Draw();

  TCanvas *c11 = new TCanvas("c11", "TAC Histogram", 800, 600);
  c11->DivideSquare(2);
  c11->cd(1);
  hist_tac_cal_1->Draw();
  c11->cd(2);
  hist_tac_cal_b->Draw();
  c11->SaveAs("hist_tac.png");

  TCanvas *c12 = new TCanvas("c12", "SI Calibrated Histogram n1", 800, 600);
  hist_si_cal_1->GetXaxis()->SetTitle("E (MeV)");
  hist_si_cal_1->Draw();
  c12->SaveAs("si_cal_1.png");

  TCanvas *c13 = new TCanvas("c13", "SI Calibrated Histogram n2", 800, 600);
  hist_si_cal_2->GetXaxis()->SetTitle("E (MeV)");
  hist_si_cal_2->Draw();
  c13->SaveAs("si_cal_2.png");

  TCanvas *c14 = new TCanvas("c14", "SI Calibrated Histogram b", 800, 600);
  hist_si_cal_b->GetXaxis()->SetTitle("E (MeV)");
  hist_si_cal_b->Draw();

  // TCanvas *c13 = new TCanvas("c13", "TAC unCalibrated Histogram n1", 800,
  // 600); hist_tac_1->GetXaxis()->SetTitle("time (arb. units)");
  // hist_tac_1->Draw();
  // c13->SaveAs("tac_cal_1.png");

  // TCanvas *c14 = new TCanvas("c14", "TAC Calibrated Histogram n1", 800, 600);
  // hist_tac_cal_1->GetXaxis()->SetTitle("time (ns)");
  // hist_tac_cal_1->Draw();
  // c14->SaveAs("tac_cal_1.png");

  // TCanvas *c15 = new TCanvas("c15", "TAC Calibrated Histogram n2", 800, 600);
  // hist_tac_cal_2->GetXaxis()->SetTitle("TOF (s)");
  // hist_tac_cal_2->Draw();
  // c14->SaveAs("tac_cal_2.png");

  // TCanvas *c16 = new TCanvas("c16", "Si Calibrated Histogram n12", 800, 600);
  // hist_faster_TAC->Draw("colz");
  // c16->SaveAs("tac_cal_12.png");
}
#endif
