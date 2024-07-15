#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TUUID.h"

#include "CalibrationRunner.h"
#include "CalibrationSource.h"

#include <iostream>
#include <vector>

#include "CalibrationRunner.cxx"
#include "CalibrationSource.cxx"

void e863_Si2() {

  // Ouvrir le fichier ROOT
  TFile *file = TFile::Open("Inputs/E863/Si2_run008.root");
  if (!file || file->IsZombie()) {
    std::cerr << "Impossible d'ouvrir le fichier!" << std::endl;
  }

  // Extraire le TCanvas
  TCanvas *canvas = nullptr;
  file->GetObject("c17", canvas);
  if (!canvas) {
    std::cerr << "TCanvas non trouvé dans le fichier!" << std::endl;
    file->Close();
    return;
  }

  // Extraire l'histogramme TH1D du TCanvas
  TH1D *hist = dynamic_cast<TH1D *>(canvas->GetPrimitive("hist_si_2"));
  if (!hist) {
    std::cerr << "Histogramme TH1D non trouvé dans le TCanvas!" << std::endl;
    file->Close();
    return;
  }

  // Create source
  Calibration::Source source;
  source.Print();

  auto* histrebin {(TH1D*)hist->Clone("histRebin")};
  histrebin->Rebin(3);
  Calibration::Runner runner{&source, histrebin, hist, false};
  runner.SetRange(68000, 82000);
  runner.DoIt();
  runner.Draw(new TCanvas);
  cout << "SI_1" << endl;
  runner.PrintRes();

  // Plot
  auto *c0{new TCanvas{"c0", "Silicon canvas"}};
  hist->Draw();
}
