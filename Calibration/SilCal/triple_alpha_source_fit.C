#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TMath.h"
#include <iostream>

void triple_alpha_source_fit() {
  double i = 1;

  // Select the file
  TFile *f = new TFile("Inputs/E863/Si1_run020.root", "READ");

  // Check if the file is open successfully
  if (!f || f->IsZombie()) {
    std::cout << "Cannot open file" << std::endl;
    return;
  }

  // Extract the TCanvas
  TCanvas *canvas = nullptr;
  f->GetObject("c17", canvas);
  if (!canvas) {
    std::cerr << "TCanvas not found in the file!" << std::endl;
    f->Close();
    return;
  }

  // Extract the TH1D histogram from the TCanvas
  TH1D *histo = dynamic_cast<TH1D *>(canvas->GetPrimitive("hist_si_1"));
  if (!histo) {
    std::cerr << "TH1D histogram not found in the TCanvas!" << std::endl;
    f->Close();
    return;
  }

  std::cout << "Histogram found: " << histo->GetName() << std::endl;

  // Check histogram content
  int nBins = histo->GetNbinsX();
  double integral = histo->Integral();
  std::cout << "Number of bins: " << nBins << std::endl;
  std::cout << "Integral of histogram: " << integral << std::endl;

  if (integral == 0) {
    std::cerr << "Histogram is empty!" << std::endl;
    f->Close();
    return;
  }

  // Gaussian fit
  TF1 *fitFunction =
      new TF1("fitFunction", "gaus(0) + gaus(3) + gaus(6)",
              histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
  fitFunction->SetParNames("A1", "Mean1", "Sigma1", "A2", "Mean2", "Sigma2",
                           "A3", "Mean3", "Sigma3");

  // Initialize the values of the parameters
  fitFunction->SetParameters(700, 57232, 200, 800, 60898, 200, 500, 64442, 200);

  fitFunction->SetLineColor(kBlue);
  fitFunction->SetLineWidth(3);
  histo->Fit(fitFunction, "R");

  // Calibration
  double mean1 = fitFunction->GetParameter(1);
  double mean2 = fitFunction->GetParameter(4);
  double mean3 = fitFunction->GetParameter(7);

  // Energies of the peaks in keV
  double E1 = 5156;
  double E2 = 5485.6;
  double E3 = 5804.8;

  TF1 *linearFit =
      new TF1("linearFit", "[0]*x + [1]", mean1 - 1000, mean3 + 1000);
  linearFit->SetParameters((E3 - E1) / (mean3 - mean1),
                           E1 - (E3 - E1) / (mean3 - mean1) * mean1);

  // Create the graph
  TGraphErrors *graph = new TGraphErrors(3);
  graph->SetPoint(0, mean1, E1);
  graph->SetPoint(1, mean2, E2);
  graph->SetPoint(2, mean3, E3);
  graph->Fit(linearFit, "QR");

  double slope = linearFit->GetParameter(0);

  // Sigma values
  double sigma1 = fitFunction->GetParameter(2);
  double err_sigma1 = fitFunction->GetParError(2);
  double sigma2 = fitFunction->GetParameter(5);
  double err_sigma2 = fitFunction->GetParError(5);
  double sigma3 = fitFunction->GetParameter(8);
  double err_sigma3 = fitFunction->GetParError(8);

  double avgSigma =
      ((std::abs(sigma1) + std::abs(sigma2) + std::abs(sigma3)) / 3.0) *
      slope; // in keV
  double avg_err = sqrt(pow(err_sigma1,2)+pow(err_sigma2,2)+pow(err_sigma3,2)/3);

  std::cout << "slope = " << slope << std::endl;
  std::cout << "Sigma_" << i << " = " << avgSigma << "+-" << avg_err << " keV" << std::endl;
  std::cout << "FWHM_" << i << " = " << avgSigma*2.35 << "+-" << avg_err*2.35 << " keV" << std::endl;

  // Draw
  auto *c0{new TCanvas{"c0", "Silicon canvas"}};
  histo->Draw();

  // Close the file
  f->Close();
}