#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TMath.h"
#include <iostream>
#include <vector>

void triple_alpha_source_fit_Si1() {
  const int nPeaks = 9; // Number of peaks
  const int nIterations = 10; // Number of fit iterations
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

  // Energies of the peaks in keV
  std::vector<double> energies = {5156.59, 5114.43, 5105.50, 
  5485.56,5442.80, 5388.23, 5352.97, 
  5804.77, 5762.64};

  // Initial guess for parameters: amplitude, mean, sigma for each peak
  std::vector<double> initialParams = {
      3300, 57232, 100, 1000, 57000, 100, 500, 56670, 100,

      4000, 60898, 100, 1500, 60780, 100, 600, 60415, 100, 100, 59830, 100,

      2700, 64442, 100, 800,  63950, 100};

  // Gaussian fit with nPeaks components
  std::string formula = "";
  for (int j = 0; j < nPeaks; ++j) {
    formula += "gaus(" + std::to_string(j * 3) + ")";
    if (j < nPeaks - 1) {
      formula += " + ";
    }
  }
  
  TF1 *fitFunction =
      new TF1("fitFunction", formula.c_str(), histo->GetXaxis()->GetXmin(),
              histo->GetXaxis()->GetXmax());
  fitFunction->SetNpx(2000);

  for (int j = 0; j < nPeaks; ++j) {
    fitFunction->SetParameter(j * 3, initialParams[j * 3]);
    fitFunction->SetParameter(j * 3 + 1, initialParams[j * 3 + 1]);
    fitFunction->SetParameter(j * 3 + 2, initialParams[j * 3 + 2]);
    fitFunction->SetParLimits(j*3+2,70,73.2);
  }

  fitFunction->SetLineColor(kBlue);
  fitFunction->SetLineWidth(3);
  histo->Rebin(2);
  // Perform multiple fit iterations
  for (int iter = 0; iter < nIterations; ++iter) {
    std::cout << "Iteration " << iter + 1 << std::endl;
    histo->Fit(fitFunction, "R");

    // Extract the fit parameters and use them as initial guesses for the next
    // iteration
    for (int j = 0; j < nPeaks; ++j) {
      initialParams[j * 3] = fitFunction->GetParameter(j * 3);
      initialParams[j * 3 + 1] = fitFunction->GetParameter(j * 3 + 1);
      initialParams[j * 3 + 2] = fitFunction->GetParameter(j * 3 + 2);
    }
  }
  //histo->Fit(fitFunction, "R");

  std::vector<double> means(nPeaks);
  std::vector<double> sigmas(nPeaks);
  std::vector<double> err_sigmas(nPeaks);

  for (int j = 0; j < nPeaks; ++j) {
    means[j] = fitFunction->GetParameter(j * 3 + 1);
    sigmas[j] = fitFunction->GetParameter(j * 3 + 2);
    err_sigmas[j] = fitFunction->GetParError(j * 3 + 2);
  }

  TF1 *linearFit = new TF1("linearFit", "[0]*x + [1]", means[0] - 1000,
                           means[nPeaks - 1] + 1000);
  linearFit->SetParameters(
      (energies.back() - energies.front()) / (means.back() - means.front()),
      energies.front() - (energies.back() - energies.front()) /
                             (means.back() - means.front()) * means.front());

  // Create the graph
  TGraphErrors *graph = new TGraphErrors(nPeaks);
  for (int j = 0; j < nPeaks; ++j) {
    graph->SetPoint(j, means[j], energies[j]);
  }
  graph->Fit(linearFit, "QR");

  //double slope = linearFit->GetParameter(0);
  double slope = 0.0911739;

  double avgSigma = 0;
  double sum_errors_squared = 0;
  for (int j = 0; j < nPeaks; ++j) {
    avgSigma += std::abs(sigmas[j]);
    sum_errors_squared += pow(err_sigmas[j], 2);
    std::cout << "FWHM_" << j << " = " << sigmas[j] *slope * 2.35 << "keV" << std::endl;
  }
  avgSigma = (avgSigma / nPeaks) * slope; // in keV
  double avg_err = sqrt(sum_errors_squared / nPeaks);

  std::cout << "slope = " << slope << std::endl;
  std::cout << "Sigma_" << i << " = " << avgSigma << "+-" << avg_err << " keV"
            << std::endl;
  std::cout << "FWHM_" << i << " = " << avgSigma * 2.35 << "+-"
            << avg_err * 2.35 << " keV" << std::endl;

  // Draw
  auto *c0 = new TCanvas("c0", "Silicon canvas");
  histo->Draw();

  // Create a new canvas for the linear fit plot
  auto *c1 = new TCanvas("c1", "Energy vs Channel");
  graph->SetTitle("Energy vs Channel;Channel;Energy (keV)");
  graph->Draw("AP");
  linearFit->Draw("same");

  // Close the file
  f->Close();
}
