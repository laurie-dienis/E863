//***************************************
// Simulation resonant scattering reactions in gaseous target
//+ Si detection at 0 degree
//***************************************

//***************************************
// TO DO
// (ii) include contaminant reactions (p detected?)
// (iii) include SRIM expected 4He straggling
//***************************************

#include "../SRIM.cxx"
#include <iostream>
#include <stdlib.h>

//***************************************
// Extraction function for AZURE elastic CS
//***************************************
int extraction_data_CS(const Char_t *file, Double_t Ex_CS[], Double_t Ecm_CS[],
                       Double_t CS_diff[]) {
  std::ifstream in;
  in.open(file);
  Int_t nlines = 0;
  while (1) {
    in >> Ecm_CS[nlines] >> Ex_CS[nlines] >>
        CS_diff[nlines]; // energies in MeV, cs b/sr
    if (!in.good())
      break;
    nlines++;
  }
  in.close();
  return 1;
}

//***************************************
// Search function for cs@Ecm
//***************************************
double find_cs(double Ecm_reaction, Double_t Ecm_CS[], Double_t CS_diff[],
               int Nentry) {
  double cs_find = 0.;
  double diffEcm = 1000.;
  for (int i = 0; i < Nentry; i++) {
    if (TMath::Abs(Ecm_CS[i] - Ecm_reaction) < diffEcm) {
      diffEcm = TMath::Abs(Ecm_CS[i] - Ecm_reaction);
      cs_find = CS_diff[i];
    }
  }
  return cs_find;
}

//***************************************
// Poly2 function from MC simulations
// Ea_cm_reaction = f(Ea_lab_detector)
//***************************************
double conversion_Emeas_EcmReac(double Emeas) {
  // Ti    5 microns at exit
  // double p0 = 0.279022;
  // double p1 = 0.315438;
  // double p2 = 0.000164993;

  // Ti    6 microns at exit
  double p0 = 0.202280;
  double p1 = 0.339113;
  double p2 = -0.000262;

  double Ecm_reaction = p0 + p1 * Emeas + p2 * Emeas * Emeas;
  return Ecm_reaction;
}

//***************************************
// Derivation of target effective width
//***************************************
double derivation_eff_target_width(double Ecm, double deltaEcm,
                                   Physics::SRIM *srim) {
  double mu = 931.5;
  double m15N = 15 * mu + 0.101;
  double m4He = 4 * mu + 2.425;
  double mu4He = m4He / (m4He + m15N);
  double Eflab = (Ecm - deltaEcm) / mu4He;
  double elab_temp = Ecm / mu4He;
  double width = 0.;
  double deltaE = elab_temp - Eflab;
  double dx = 0.01; // mm
  double eloss = 0.;
  while (deltaE > 0) {
    width = width + dx;
    elab_temp = srim->Slow("light_gas", elab_temp, dx, 0.);
    deltaE = elab_temp - Eflab;
  }
  return width;
}

//***************************************
// Derivation scattered beam after entrance window
//***************************************
double scattered_theta() {
  TFile *file1 =
      new TFile("scattering/beamEntry/theta_distri_entrance.root", "READ");
  TH1F *srim_scattered = (TH1F *)file1->Get("theta_distri_entrance");
  gRandom->SetSeed(0);
  double theta = srim_scattered->GetRandom();
  file1->Close();
  return theta;
}

//***************************************
// renormalisation function to exp. stat
//***************************************
double renormalisation(TH1F *hMC, TH1F *hExp, double compt, int Nsim[3],
                       double bin_to_E, double Ibeam,
                       double concentration_target, double timeUT,
                       Physics::SRIM *srim) {
  int NbinX = hMC->GetNbinsX();
  TAxis *xaxis = hMC->GetXaxis();
  double acceptance_angular = 0.0012; // 1msr in deg <=> *pow(180./3.14,2)
  double conv_barn_mm2 = pow(10, -22);
  double renormStat = Ibeam * timeUT * concentration_target;
  // Poiss drawing from expected counts
  TRandom *rPois = new TRandom();
  double Ndraw;
  gRandom->SetSeed(0);
  double ecm, dx;
  double count;
  for (int b = 0; b < NbinX; b++) {
    ecm = xaxis->GetBinCenter(b);
    count = hMC->GetBinContent(b);
    if (count > 0) {
      dx = derivation_eff_target_width(ecm, bin_to_E, srim);
      Ndraw = rPois->Poisson(
          (renormStat * dx * conv_barn_mm2 * acceptance_angular * count) /
          (double(Nsim[2] * Nsim[1] * compt)));
      if (Ndraw > 0) {
        for (int kk = 0; kk < Ndraw; kk++) {
          hExp->Fill(ecm);
        }
      }
    }
  }
  return 1.0;
}

void printProgBar(float percent) {
  std::string bar;

  for (int i = 0; i < 50; i++) {
    if (i < (int(percent) / 2)) {
      bar.replace(i, 1, "=");
    } else if (i == int(percent / 2)) {
      bar.replace(i, 1, ">");
    } else {
      bar.replace(i, 1, " ");
    }
  }

  std::cout << "\r [" << bar << "]";
  std::cout.width(3);
  std::cout << std::fixed;
  std::cout << std::setprecision(2);

  std::cout << percent << "%" << std::flush;
}

//***************************************
// Monte Carlo function
// energy in MeV, pressure in Torr
// detector distance in mm, angle lab in deg
//***************************************
double MCsim(const Char_t *fileCS, int Nsim[3], double EbeamEntry,
             double bin_to_E, double Ibeam, int pressure_target,
             double target_width, double concentration_target, double timeUT,
             double entranceWindow_width, double exitWindow_width,
             double angle_max_Si, double dist_target_Si,
             TH2F *distriAngularAlpha, TH2F *distriToFAlpha,
             TH1F *distriEx_cm) {

  // SRIM stopping powers
  auto *srim{new Physics::SRIM};
  srim->ReadTable("light_window",
                  "/home/laurie/Analysis_e863/Inputs/SRIM/4He_Ti.txt");
  srim->ReadTable("beam_window",
                  "/home/laurie/Analysis_e863/Inputs/SRIM/15N_Ti.txt");
  srim->ReadTable(
      "beam_gas",
      TString::Format(
          "/home/laurie/Analysis_e863/Inputs/SRIM/15N_4He_%iTorr.txt",
          pressure_target)
          .Data());
  srim->ReadTable(
      "light_gas",
      TString::Format(
          "/home/laurie/Analysis_e863/Inputs/SRIM/4He_4He_%iTorr.txt",
          pressure_target)
          .Data());

  // Extraction cross_sections
  int N_entry = 25000;
  Double_t Ecm_CS[25000];
  Double_t Ex_CS[25000];
  Double_t CS_diff[25000];
  int extraction_data;
  extraction_data = extraction_data_CS(fileCS, Ex_CS, Ecm_CS, CS_diff);

  // statistical renormalisation
  int Nreac = Nsim[0];
  int Nbeam = Nsim[1];
  double acceptance_angular = 0.0012; // 1msr in deg <=> *pow(180./3.14,2)
  double conv_barn_mm2 = pow(10, -22);
  double renormStat = Ibeam * timeUT * concentration_target;
  // Poiss drawing from expected counts
  TRandom *rPois = new TRandom();
  int Ndraw;

  // Si detector response
  TRandom *rEloss = new TRandom(); // gaussian reponse of Si detector
  double sigma_e = 0.001 * 10.0;   // MeV
  srand(time(NULL));
  gRandom->SetSeed(0);

  // reaction features
  double c_light = 3.0 * pow(10, 2); // mm.ns-1
  double mu = 931.5;
  double m15N = 15 * mu + 0.101;
  double m4He = 4 * mu + 2.425;
  double mu4He = m4He / (m4He + m15N);
  double mu15N = m15N / (m4He + m15N);
  double m1H = mu + 7.289;
  gSystem->Load("libPhysics");
  double Pfais = sqrt(EbeamEntry * EbeamEntry + 2. * EbeamEntry * m15N);
  TLorentzVector target(0.0, 0.0, 0.0, m4He * 0.001);
  TLorentzVector beam(0.0, 0.0, Pfais * 0.001, (EbeamEntry + m15N) * 0.001);
  Double_t masses[2] = {0.001 * (m15N), 0.001 * m4He};
  TLorentzVector Tr = beam + target;
  TLorentzVector *pRecoil;
  TLorentzVector *pEjectil;
  TGenPhaseSpace event;
  event.SetDecay(Tr, 2, masses);
  // Variables
  double dx = 0.07; // mm
  double Ebeam, thetaBeam, beam_target_width, vbeamlab, ELossBeam, ElossRecoil;
  double thetaEjectil_lab, vEjectil_lab, eEjectil_lab, ElossEjectil,
      distEjectil, remainingDist, vEjectil_lab_meas, vProton_lab_meas,
      Ecm_reaction;
  double Emeas;
  double thetaExit;
  double positionBeam = 0;
  double ponderation_cs = 0.;
  double ToF = 0.0;
  double ToFproton = 0.0;
  double beam_ToF = 0.0;
  int compt;
  float step = Nbeam / 100.0;

  //********************************
  // Derivation reactions
  //********************************
  for (int b = 0; b < Nbeam; b++) {
    // std::cout << "Current iteration: " << b << std::endl << std::flush;
    beam_ToF = 0.;
    positionBeam = 0;
    Ebeam = EbeamEntry;
    // std::cout<<" beam entrance E "<<  Ebeam/15.001138  <<std::endl;
    thetaBeam = scattered_theta();
    beam_target_width = target_width / cos(3.14 * (thetaBeam / 180.));
    Ebeam = srim->Slow("beam_window", Ebeam, entranceWindow_width,
                       thetaBeam * TMath::DegToRad());
    // std::cout<<" beam entrance "<<
    // entranceWindow_width/cos(3.14*(thetaBeam/180.)) <<" loss E "<< ELossBeam
    // <<std::endl; std::cout<<" beam entrance after wind "<<  Ebeam/15.001138
    // <<std::endl;
    while (positionBeam < (beam_target_width - dx)) {
      Ecm_reaction = mu4He * Ebeam;
      // dx = derivation_eff_target_width(Ecm_reaction, bin_to_E, fileSRIM[1]);
      positionBeam = positionBeam + dx;
      // std::cout<< positionBeam<<" energy beam "<<Ebeam/15.001138<<" Ecm "<<
      // Ecm_reaction<<std::endl;
      Ebeam = srim->Slow("beam_gas", Ebeam, dx, thetaBeam * TMath::DegToRad());
      vbeamlab = sqrt(1 - 1 / pow((Ebeam / m15N + 1), 2));
      Pfais = sqrt(Ebeam * Ebeam + 2. * Ebeam * m15N);
      beam.SetPxPyPzE(0, 0, Pfais * 0.001, (Ebeam + m15N) * 0.001);
      Tr = beam + target;
      event.SetDecay(Tr, 2, masses);
      beam_ToF = beam_ToF + dx / (vbeamlab * c_light);
      ponderation_cs = find_cs((Ecm_reaction + Ebeam * mu4He) / 2.0, Ecm_CS,
                               CS_diff, N_entry);
      for (int j = 0; j < Nreac; j++) {
        event.Generate();
        pRecoil = event.GetDecay(0);
        pEjectil = event.GetDecay(1);
        thetaEjectil_lab = 180. / 3.14 * pEjectil->Theta(); // Rad
        vEjectil_lab = pEjectil->Beta();
        eEjectil_lab = (1 / sqrt(1 - pow(vEjectil_lab, 2)) - 1) * m4He;
        remainingDist = beam_target_width - positionBeam;
        if (eEjectil_lab > 0 && thetaEjectil_lab < angle_max_Si) {
          distEjectil = (remainingDist) / cos(3.14 * thetaEjectil_lab / 180.);
          eEjectil_lab = srim->Slow("light_gas", eEjectil_lab, distEjectil,
                                    thetaEjectil_lab * TMath::DegToRad());
          eEjectil_lab =
              srim->Slow("light_window", eEjectil_lab, exitWindow_width,
                         thetaEjectil_lab * TMath::DegToRad());
          Emeas = rEloss->Gaus(eEjectil_lab, sigma_e);
          vEjectil_lab_meas = sqrt(1 - 1 / pow((Emeas / m4He + 1), 2));
          vProton_lab_meas = sqrt(1 - 1 / pow((Emeas / m1H + 1), 2));
          ToF = beam_ToF + ((remainingDist + dist_target_Si) /
                            cos(3.14 * thetaEjectil_lab / 180.)) /
                               (vEjectil_lab_meas * c_light);
          ToFproton = beam_ToF + ((remainingDist + dist_target_Si) /
                                  cos(3.14 * thetaEjectil_lab / 180.)) /
                                     (vProton_lab_meas * c_light);
          // std::cout<<"meas. alphas E "<<Emeas<<std::endl;
          for (int kk = 0; kk < ponderation_cs * Nsim[2]; kk++) {
            distriAngularAlpha->Fill(Emeas, thetaEjectil_lab);
            distriToFAlpha->Fill(Emeas, ToF);
            distriEx_cm->Fill(conversion_Emeas_EcmReac(Emeas));
          }
          distriToFAlpha->Fill(Emeas, ToFproton);
        }
      }
      if (b % static_cast<int>(step) == 0)
        printProgBar(b * 100.0 / Nbeam);
    }
  }
  // std::cout<<"meas. alphas E "<<Emeas<<std::endl;
  // std::cout<<"end beam E "<<Ebeam/15.001138<<std::endl;
  return compt;
}

//***************************************
// Main Function
//***************************************
void MC15N4He() {

  std::cout << "Début de MC15N4He()" << std::endl;

  //***************************************
  // Experiment conditions
  //***************************************
  double Ibeam = 1.0 * pow(10, 6);     // pps
  double EbeamEntry = 1.7 * 15.001138; // MeV
  int choice_pressure_target = 3;
  int aivalable_pressure_target[4] = {92, 183, 274, 350};
  double target_width = 90.; // mm 85==Al   105 Au   90
  double concentration_target[4] = {
      3.01 * pow(10, 18) / 1000., 6.02 * pow(10, 18) / 1000.,
      9.03 * pow(10, 18) / 1000., 1.15 * pow(10, 19) / 1000.}; // at.mm^3
  double timeUT = 16.0 * 8. * 60. * 60.;                       // sec
  /*Ti windows*/
  double exitWindow_width = 0.006;     // mm
  double entranceWindow_width = 0.003; // mm

  double angle_max_Si = 2; // degree
  double frac_solid_angle = (1 - cos(angle_max_Si * 3.14 / 180.)) * 0.5;
  double acceptance_angular = 0.0012; // 1msr in deg <=> *pow(180./3.14,2)
  double dist_target_Si = 28.6;       // mm 343 57
  // std::cout << "Hey" << std::endl;

  //***************************************
  // MC statistics
  //***************************************
  int Nsim[3] = {50, 400, 500}; // reac(100000), beam, cs_azure
  // int Nsim[3] = {50, 400, 500}; // reac(100000), beam, cs_azure

  double compt_detec = frac_solid_angle * Nsim[0];

  //***************************************
  // histograms to plot
  //***************************************
  TCanvas *c = new TCanvas("c", "c", 800, 800);
  c->Divide(1, 2);
  TCanvas *c_renorm = new TCanvas("c_renorm ", "c_renorm ", 800, 800);
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
  auto legend2 = new TLegend(0.05, 0.2, 0.3, 0.3);
  auto legend2_renorm = new TLegend(0.05, 0.2, 0.3, 0.3);
  double Er[2] = {0., 20.};     // MeV
  double Ecmr[2] = {2.0, 10.0}; // MeV
  double thetar[2] = {0., 5.};  // deg
  double ToFr[2] = {0., 60.};   // ns
  double bin_to_theta = 0.05;
  double bin_to_E = 0.002; // kev bin
  double bin_to_ToF = 0.05;
  int binthetar = (thetar[1] - thetar[0]) / bin_to_theta;
  int binEr = (Er[1] - Er[0]) / bin_to_E;
  int binEcmr = (Ecmr[1] - Ecmr[0]) / bin_to_E;
  int binToF = (ToFr[1] - ToFr[0]) / bin_to_ToF;
  TH2F *distriAngularAlpha[2];
  TH1F *distriEx_cm[2];
  TH1F *distriEx_cm_renorm[2];
  TH2F *distriToFAlpha[2];

  //***************************************
  // Input cross-sections from AZURE2
  //***************************************
  // const Char_t *fileCS[1] = {"azure_results/cs180deg15N.dat",};
  const Char_t *fileCS[1] = {
      "azure_results/cs180deg_15N_2.dat",
  };
  // const Char_t *fileCS[1] = {"azure_results/AZUREOut_aa=1_R=1_
  // Laurie.extrap",};
  Double_t Ecm[25000];
  Double_t Ex[25000];
  Double_t CS_diff_theor[25000];
  Double_t N_theor_MC[25000];
  Double_t N_theor_renorm[25000];
  // int extraction_data = extraction_data_CS(fileCS[0],Ex, Ecm, CS_diff_theor);
  auto *srim{new Physics::SRIM};
  srim->ReadTable(
      "light_gas",
      TString::Format(
          "/home/laurie/Analysis_e863/Inputs/SRIM/4He_4He_%iTorr.txt",
          aivalable_pressure_target[choice_pressure_target])
          .Data());
  double mu = 931.5;
  double m15N = 15 * mu + 0.101;
  double m4He = 4 * mu + 2.425;
  double mu4He = m4He / (m4He + m15N);
  double effective_width = 0;
  double binAzure = 0.0002;
  for (int i = 0; i < 25000; i++) {
    // std::cout << "iteration"<< i << std::endl;
    // effective_width = derivation_eff_target_width(Ecm[i], bin_to_E, srim);
    N_theor_MC[i] = CS_diff_theor[i] * Nsim[2] * compt_detec * Nsim[1] *
                    (bin_to_E / binAzure);
    // std::cout << "iteration"<< i << std::endl;
    N_theor_renorm[i] = (CS_diff_theor[i] * pow(10, -22) * acceptance_angular) *
                        concentration_target[choice_pressure_target] * 0.08 *
                        0.1 * (Ibeam * timeUT) * 0.55 * (bin_to_E / binAzure);
    // std::cout << "CS_diff_theor = "<< CS_diff_theor[i] << std::endl;
    // std::cout << "Ntheorenorm = "<< N_theor_renorm[i] << std::endl;
  }
  TGraph *distriEx_cm_theor_MC = new TGraph(25000, Ecm, N_theor_MC);
  TGraph *distriEx_cm_theor_renorm = new TGraph(25000, Ecm, N_theor_renorm);

  //***************************************
  // Call of MC simulations
  //***************************************
  double simulation_reaction, applied_renorm;
  for (int t = 0; t < 1; t++) {
    distriAngularAlpha[t] =
        new TH2F(Form("distriAngularAlpha%i", t),
                 ";E_{#alpha} (MeV);#theta_{#alpha} (deg)", binEr, Er[0], Er[1],
                 binthetar, thetar[0], thetar[1]);
    distriEx_cm[t] =
        new TH1F(Form("distriEx_cm%i", t),
                 Form(";E_{cm} (MeV);Counts/(%.1fkeV)", 1000. * bin_to_E),
                 binEcmr, Ecmr[0], 10);
    distriEx_cm_renorm[t] =
        new TH1F(Form("distriEx_cm_renorm%i", t),
                 Form(";E_{cm} (MeV);Counts/(%.1fkeV)", 1000. * bin_to_E),
                 binEcmr, Ecmr[0], 10);
    distriToFAlpha[t] = new TH2F(Form("distriToFAlpha%i", t),
                                 ";E_{#alpha} (MeV);ToF_{#alpha} (ns)", binEr,
                                 Er[0], Er[1], binToF, ToFr[0], ToFr[1]);
    simulation_reaction = MCsim(
        fileCS[t], Nsim, EbeamEntry, bin_to_E, Ibeam,
        aivalable_pressure_target[choice_pressure_target], target_width,
        concentration_target[choice_pressure_target], timeUT,
        entranceWindow_width, exitWindow_width, angle_max_Si, dist_target_Si,
        distriAngularAlpha[t], distriToFAlpha[t], distriEx_cm[t]);
    applied_renorm = renormalisation(
        distriEx_cm[t], distriEx_cm_renorm[t], compt_detec, Nsim, bin_to_E,
        Ibeam, concentration_target[choice_pressure_target], timeUT, srim);
    std::cout << "simulation done for test " << t << std::endl;
  }

  //***************************************
  // Plotting
  //***************************************
  c->cd(1);
  distriToFAlpha[0]->Draw("colz");
  distriToFAlpha[0]->GetYaxis()->CenterTitle();
  distriToFAlpha[0]->GetYaxis()->SetTitleSize(0.05);
  distriToFAlpha[0]->GetYaxis()->SetLabelSize(0.04);
  distriToFAlpha[0]->GetXaxis()->CenterTitle();
  distriToFAlpha[0]->GetXaxis()->SetTitleSize(0.05);
  distriToFAlpha[0]->GetXaxis()->SetLabelSize(0.04);
  c->cd(2);
  distriAngularAlpha[0]->Draw("colz");
  distriAngularAlpha[0]->GetYaxis()->CenterTitle();
  distriAngularAlpha[0]->GetYaxis()->SetTitleSize(0.05);
  distriAngularAlpha[0]->GetYaxis()->SetLabelSize(0.04);
  distriAngularAlpha[0]->GetXaxis()->CenterTitle();
  distriAngularAlpha[0]->GetXaxis()->SetTitleSize(0.05);
  distriAngularAlpha[0]->GetXaxis()->SetLabelSize(0.04);
  //***************************************
  // Excitation spectra in E_cm
  //***************************************
  c2->cd();
  distriEx_cm[0]->Draw("her");
  legend2->AddEntry(distriEx_cm[0],
                    "Ex=6.441 MeV (J=#frac{3}{2}^{+}, #Gamma_{#alpha}=1.3 keV)",
                    "l");
  legend2->AddEntry(distriEx_cm_theor_MC, "Theoretical calculations", "l");
  distriEx_cm_theor_MC->Draw("sameC");
  distriEx_cm_theor_MC->SetLineColor(2);
  distriEx_cm[0]->GetYaxis()->CenterTitle();
  distriEx_cm[0]->GetYaxis()->SetTitleSize(0.05);
  distriEx_cm[0]->GetYaxis()->SetLabelSize(0.04);
  distriEx_cm[0]->GetXaxis()->CenterTitle();
  distriEx_cm[0]->GetXaxis()->SetTitleSize(0.05);
  distriEx_cm[0]->GetXaxis()->SetLabelSize(0.04);
  distriEx_cm[0]->GetXaxis()->SetRangeUser(2.0, 4.0);
  distriEx_cm[0]->GetYaxis()->SetRangeUser(0.0, 1000000.);
  c_renorm->cd();
  distriEx_cm_renorm[0]->Draw("her");
  legend2_renorm->AddEntry(
      distriEx_cm_renorm[0],
      "with level 5 (Ex=6.441 MeV, J=#frac{3}{2}^{+}, #Gamma_{#alpha}=1.3 keV)",
      "l");
  legend2_renorm->AddEntry(distriEx_cm_renorm[1], "without level 5", "l");
  legend2_renorm->AddEntry(distriEx_cm_theor_renorm, "Theory", "l");
  distriEx_cm_theor_renorm->Draw("sameC");
  distriEx_cm_theor_renorm->SetLineColor(2);
  distriEx_cm_renorm[0]->GetYaxis()->CenterTitle();
  distriEx_cm_renorm[0]->GetYaxis()->SetTitleSize(0.05);
  distriEx_cm_renorm[0]->GetYaxis()->SetLabelSize(0.04);
  distriEx_cm_renorm[0]->GetXaxis()->CenterTitle();
  distriEx_cm_renorm[0]->GetXaxis()->SetTitleSize(0.05);
  distriEx_cm_renorm[0]->GetXaxis()->SetLabelSize(0.04);
  distriEx_cm_renorm[0]->GetXaxis()->SetRangeUser(2.0, 4.0);
  // distriEx_cm_renorm[0]->GetYaxis()->SetRangeUser(0.0, 150.);
}
