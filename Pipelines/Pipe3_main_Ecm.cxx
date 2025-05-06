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

TH1F *extract_TH1F_from_TCanvas(const char *filename, const char *canvas_name) {
  std::cout << "Opening file: " << filename << std::endl;
  TFile *file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return nullptr;
  }

  TCanvas *canvas = dynamic_cast<TCanvas *>(file->Get(canvas_name));
  if (!canvas) {
    file->Close();
    return nullptr;
  }

  TList *primitives = canvas->GetListOfPrimitives();
  if (!primitives) {
    file->Close();
    return nullptr;
  }

  TH1F *hist = nullptr;
  for (TObject *obj : *primitives) {
    hist = dynamic_cast<TH1F *>(obj);
    if (hist)
      break;
  }

  if (!hist) {
    file->Close();
    return nullptr;
  }

  TH1F *hist_clone = dynamic_cast<TH1F *>(hist->Clone("hist_clone"));
  if (hist_clone)
    hist_clone->SetDirectory(0);

  file->Close();
  return hist_clone;
}

void Pipe3_main_Ecm() {
  std::vector<std::string> fileNames = {"RootFiles/PID/alpha_1.root"};
  const char *treeName = "DataTree";

  double si, si_cm;
  double p0 = 0.000772, p1 = 0.277615, p2 = -0.000192;

  TH1F *hist_E = new TH1F("hist_si_cal", "SI Calibrated Histogram", 500, 1, 10);
  TH1F *hist_Ecm =
      new TH1F("hist_si_cal_Ecm", "SI Calibrated Histogram", 500, 0, 6);

  for (const auto &fileName : fileNames) {
    TFile *file = TFile::Open(fileName.c_str());
    if (!file || file->IsZombie())
      continue;

    TTree *tree = nullptr;
    file->GetObject(treeName, tree);
    if (!tree) {
      file->Close();
      continue;
    }

    tree->SetBranchAddress("Si_cal_1", &si);
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
      tree->GetEntry(i);
      hist_E->Fill(si);
      si_cm = p0 + p1 * si + p2 * si * si;
      hist_Ecm->Fill(si_cm);
    }

    file->Close();
    delete file;
  }

  TH1F *hist_sim =
      extract_TH1F_from_TCanvas("Inputs/Sim/c_renorm_4.root", "c_renorm");
  if (hist_sim) {
    hist_sim->SetLineColor(kRed); // Ligne rouge
    hist_sim->SetLineWidth(2);    // Épaisseur de ligne
    hist_sim->SetMarkerStyle(1);  // Pas de marqueur (pas de points)
    hist_sim->Scale(0.007);       // Normalisation
    hist_sim->Rebin(4);           // Rebinning
  }

  // Plot hist_E
  TCanvas *c1 = new TCanvas("c1", "SI Histogram", 800, 600);
  hist_E->Draw();
  c1->SaveAs("hist_si_alpha.png");

  // Plot hist_Ecm
  TCanvas *c2 = new TCanvas("c2", "SI_CM Histogram", 800, 600);
  hist_Ecm->Draw();
  if (hist_sim)
    hist_sim->Draw("SAME C");
  c2->SaveAs("hist_si_alpha_Ecm.png");
}

#endif
