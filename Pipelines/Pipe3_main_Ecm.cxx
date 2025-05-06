#ifndef Pipe3_Ecm_cxx
#define Pipe3_Ecm_cxx

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TROOT.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"

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

void Pipe3_main_Ecm(TString runNumber) {
  std::string str_runN = std::string(runNumber.Data());

  bool comparison = true;
  
  TChain *chain = new TChain("DataTree");  

  //////////////////////////////////////////////
  ////////// all runs option : /////////////////
  //////////////////////////////////////////////
  if (runNumber == "NONE")
    chain->Add("RootFiles/PID/15N/alpha_1_*.root");
  else if (str_runN.find(":") != std::string::npos){
    int split_pos = str_runN.find(":");
    std::string str_runN1 = str_runN.substr(0, split_pos);
    unsigned int runN1 = stoi(str_runN1);
    std::string str_runN2 = str_runN.substr(split_pos + 1,str_runN.size());
    unsigned int runN2 = stoi(str_runN2);
    for (unsigned int irunN = runN1; irunN < runN2 + 1; irunN++){
     chain->Add(Form("RootFiles/PID/15N/alpha_1_%i.root", irunN));
    }
  }
  
  double si, si_cm;
//  double p0 = 0.368896, p1 = 0.309091, p2 = 0.000396; // 350mbar NOW WE RUN AT 570 mbar so it is wrong
//  double p0 = 0.595918, p1 = 0.269029, p2 = 0.00152183; // 428 Torr, equivalent to 570 mbar 
  // double p0 = 0.311254, p1 = 0.316234, p2 = 0.000258; // 428 Torr, equivalent to 570 mbar 15N
  // double p0 = 0.369155, p1 = 0.309036, p2 = 0.000399; // 428 Torr, equivalent to 570 mbar 15O
  double p0 = 0.485602, p1 = 0.292844, p2 = 0.001174; // 428 Torr, equivalent to 570 mbar 15O

  TH1F *hist_E = new TH1F("hist_si_cal", "SI Calibrated Histogram", 130, 4, 14);
  TH1F *hist_Ecm =
      new TH1F("hist_si_cal_Ecm", "SI Calibrated Histogram", 200, 1.9, 4.3);

  chain->SetBranchAddress("Si_cal_1", &si);
  Long64_t nEntries = chain->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
      chain->GetEntry(i);
      hist_E->Fill(si);
      si_cm = p0 + p1 * si + p2 * si * si;
      hist_Ecm->Fill(si_cm);
  }

  // TH1F *hist_sim =
      // extract_TH1F_from_TCanvas("Inputs/Sim/c_renorm_5.root", "c_renorm");
  TH1F *hist_sim =
      extract_TH1F_from_TCanvas("Inputs/Sim/c_renorm_5.root", "c2");
  if (hist_sim) {
    hist_sim->SetLineColor(kRed); // Ligne rouge
    hist_sim->SetLineWidth(2);    // Épaisseur de ligne
    hist_sim->SetMarkerStyle(1);  // Pas de marqueur (pas de points)
    // hist_sim->Scale(0.07);       // Normalisation
    // hist_sim->Scale(0.2);       // Normalisation
    hist_sim->Scale(6e-7*hist_E->Integral()/3.);
    hist_sim->Rebin(2);           // Rebinning
    hist_sim->Smooth(5);
  }

  std::vector<double> v_E;
  std::vector<double> v_N;
  unsigned int Nsize = 0;
  unsigned int hist_sim_size = hist_sim->GetNbinsX();
  for(unsigned int i = 0 ; i < hist_sim_size ; i++){
    double E = hist_sim->GetBinCenter(i)+0.18;
    double N = hist_sim->GetBinContent(i);
    if(N > 0){
     v_E.push_back(E);
     v_N.push_back(N);
     Nsize++;
    }
  }
  
  TGraph* gsim = new TGraph(Nsize, &v_E[0], &v_N[0]);
  gsim->SetLineColor(kRed); // Ligne rouge
  gsim->SetLineWidth(2);    // Épaisseur de ligne

  // Plot hist_E
  TCanvas *c1 = new TCanvas("c1", "SI Histogram", 800, 600);
  hist_E->Draw();
  // hist_E->Rebin(2);
  c1->SaveAs("hist_si_alpha.png");

  // Plot hist_Ecm
  TCanvas *c2 = new TCanvas("c2", "SI_CM Histogram", 800, 600);
  hist_Ecm->Draw();
  hist_Ecm->GetYaxis()->SetTitle("Counts/2keV");
  hist_Ecm->GetXaxis()->SetTitle("E_{CM} (MeV)");
  // hist_Ecm->Rebin(2);
  // gsim->Draw("SAME");
  c2->SaveAs("hist_si_alpha_Ecm.png");
}

#endif
