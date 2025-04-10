//Simulation resonant scattering reactions in gaseous target
//+ Si detection at 0 degree
//***************************************

#include <iostream>
#include <stdlib.h>
#include "TH2F.h"
#include "../SRIM.cxx"


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

void FitPolynomial(TH2F* hist) {
    // Créer une fonction TF2 pour le fit
    TF2 *fitFunc = new TF2("fitFunc", "pol2", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), hist->GetYaxis()->GetXmin(), hist->GetYaxis()->GetXmax());

    // Faire l'ajustement
    hist->Fit(fitFunc, "Q");

    // Afficher les paramètres ajustés
    double par0 = fitFunc->GetParameter(0);
    double par1 = fitFunc->GetParameter(1);
    double par2 = fitFunc->GetParameter(2);

    printf("Fit Parameters:\n");
    printf("   par0: %f\n", par0);
    printf("   par1: %f\n", par1);
    printf("   par2: %f\n", par2);
}


//***************************************
//Derivation scattered beam after entrance window
//***************************************
double scattered_theta() {
    TFile* file1 = new TFile("scattering/beamEntry/theta_distri_entrance.root", "READ");
    TH1F* srim_scattered = (TH1F*)file1->Get("theta_distri_entrance");
    gRandom->SetSeed(0);
    double theta = srim_scattered->GetRandom();
    file1->Close();
    return theta;
}

//***************************************
// Main Function
//***************************************
void poly_Eameas_Eacm_15N(){
    
    //***************************************
    //Experiment conditions
    //***************************************
    int choice_pressure_target = 2;
    double EbeamEntry = 1.8*15.0030656 ;// MeV
    int aivalable_pressure_target[4]={100,150,350,410};
    double target_width = 90;//mm Al==85  Au105.
    /*Al windows*/
  /*  double exitWindow_width = 0.020; //mm
    double entranceWindow_width=0.004;//mm*/
    /*Au windows*/
    double exitWindow_width =0.006; //mm
    double entranceWindow_width=0.0030;//mm
    double angle_max_Si = 2;//degree
    //Si detector response
    TRandom* rEloss = new TRandom();//gaussian reponse of Si detector
    double sigma_e = 0.001*10.0;//MeV
    srand(time(NULL));
    gRandom->SetSeed(0);
    
    TCanvas* c = new TCanvas("c","c",800,800);
   

  auto *srim{new Physics::SRIM};
  srim->ReadTable(
      "light",
      TString::Format(
          "/home/laurie/Analysis_e863/Inputs/SRIM/4He_4He_%dTorr.txt", aivalable_pressure_target[choice_pressure_target])
          .Data());
  srim->ReadTable(
      "light_window","/home/laurie/Analysis_e863/Inputs/SRIM/4He_Ti.txt");
  srim->ReadTable(
      "beam",
      TString::Format(
          "/home/laurie/Analysis_e863/Inputs/SRIM/15N_4He_%dTorr.txt", aivalable_pressure_target[choice_pressure_target])
          .Data());
  srim->ReadTable(
      "beam_window","/home/laurie/Analysis_e863/Inputs/SRIM/15N_Ti.txt");
    
    //***************************************
    //Reaction initialization
    //***************************************
    double mu = 931.5;
    double m15N = 15 * mu +0.101 ;
    double m4He = 4 * mu + 2.425;
    double mu4He=m4He /(m4He + m15N);
    double mu15N=m15N /(m4He + m15N);
    gSystem->Load("libPhysics");
    double Pfais = sqrt(EbeamEntry * EbeamEntry + 2. * EbeamEntry * m15N);
    TLorentzVector target(0.0, 0.0, 0.0, m4He * 0.001);
    TLorentzVector beam(0.0, 0.0, Pfais * 0.001, (EbeamEntry + m15N) * 0.001);
    Double_t masses[2] ={0.001 * (m15N), 0.001 * m4He};
    TLorentzVector Tr = beam + target;
    TLorentzVector * pRecoil;
    TLorentzVector * pEjectil;
    TGenPhaseSpace event;
    event.SetDecay(Tr, 2, masses);
    //***************************************
    //Variables MC simu
    //***************************************
    int Nreac =50000;//500000
    int Nbeam =500;//5000
    float step = Nbeam / 100.0;
    double dx=0.07;//mm
    int Nentry = (target_width/dx)*double(Nbeam)-1;
    double Ebeam, thetaBeam, vbeamlab, ELossBeam, ElossRecoil;
    Ebeam = EbeamEntry;
    double thetaEjectil_lab, vEjectil_lab, eEjectil_lab, ElossEjectil, distEjectil, remainingDist, vEjectil_lab_meas, Ecm_reaction;
    double Emeas; int compteur=0; int compteur_reaction=0;
    double positionBeam=0;

    TH2F *h1 = new TH2F("h1", "Emeas_Ecm", 300, 6, 14, 300, 2, 5);
    h1->SetNameTitle(
        "1","""variation;Emeas;Ecm_reaction");
    
    //********************************
    //Derivation reactions
    //********************************
    double etemplossbeam=0;
    for(int b=0;b<Nbeam;b++){
        Ebeam = EbeamEntry;
        positionBeam= 0 ;   
        thetaBeam = 0.0;
        double beam_target_width =  target_width;
        std::cout<< b <<"   "<<std::endl;
        Ebeam = srim->Slow("beam_window", Ebeam, entranceWindow_width, thetaBeam);

        while(positionBeam<(beam_target_width-dx)){
        compteur_reaction=0;
        Ecm_reaction=mu4He*Ebeam;
        positionBeam = positionBeam + dx;
        std::cout<< Ebeam <<" Ebeam before  "<<std::endl;
        std::cout<< dx <<" dx  "<<std::endl;
        Ebeam = srim->Slow("beam", Ebeam, dx, thetaBeam);
        std::cout<< Ebeam <<" Ebeam after "<<std::endl;
        vbeamlab = sqrt(1 - 1 / pow((Ebeam / m15N + 1), 2));
        Pfais = sqrt(Ebeam * Ebeam + 2. * Ebeam * m15N);
        beam.SetPxPyPzE(0, 0, Pfais * 0.001, (Ebeam+ m15N) * 0.001);
        Tr = beam + target;
        event.SetDecay(Tr, 2, masses);
        Ecm_reaction = (Ecm_reaction+mu4He*Ebeam)/2.;
        // etemplossbeam=ELossBeam;
        for(int j=0;j<Nreac;j++){
            //std::cout<< j <<" j "<<std::endl;
            event.Generate();
            pRecoil = event.GetDecay(0); 
            pEjectil = event.GetDecay(1);
            thetaEjectil_lab = 180./3.14*pEjectil->Theta();//Rad
            vEjectil_lab =pEjectil->Beta();
            eEjectil_lab =(1 / sqrt(1 - pow(vEjectil_lab, 2)) - 1) * m4He;
            remainingDist = beam_target_width-positionBeam;
            if(eEjectil_lab>0 && thetaEjectil_lab<angle_max_Si && compteur_reaction==0){
                compteur_reaction=1;
                eEjectil_lab = srim->Slow("light", eEjectil_lab, remainingDist, thetaEjectil_lab*TMath::DegToRad());
                eEjectil_lab = srim->Slow("light_window", eEjectil_lab, exitWindow_width, thetaEjectil_lab*TMath::DegToRad());
                Emeas= rEloss->Gaus(eEjectil_lab,sigma_e);
                h1->Fill(Emeas, Ecm_reaction);
                j=Nreac;
            }
        }
    }
}
    //std::cout<<compteur<< " "<<Nentry<<std::endl;
    //std::cout<<"beam at the end "<<  Ebeam <<std::endl;
    
    //***************************************
    //Plotting
    //***************************************
    h1->Draw("col");
    //FitPolynomial(h1);

    // Fit function
    TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x + [2]*x*x", 0, 20);  // Adjust the range [0, 20] as needed

    // Set initial parameter values (optional)
    fitFunc->SetParameter(0, 0.3);
    fitFunc->SetParameter(1, 0.3);
    fitFunc->SetParameter(2, 0.0);

    // Fit the histogram
    h1->Fit(fitFunc, "Q");

    // Draw the fit function on top of the histogram
    fitFunc->Draw("same");
    
    // Afficher les paramètres ajustés
    double par0 = fitFunc->GetParameter(0);
    double par1 = fitFunc->GetParameter(1);
    double par2 = fitFunc->GetParameter(2);

    printf("Fit Parameters:\n");
    printf("   par0: %f\n", par0);
    printf("   par1: %f\n", par1);
    printf("   par2: %f\n", par2);
    
    //c->Close();
    //TGraph* function_EaMeas_EaCM = new TGraph(Nentry, Ea_meas, Ea_cm);
    //std::cout<<Ea_cm[Nentry-1]<< " "<<Ea_meas[Nentry-1]<<std::endl;
    //function_EaMeas_EaCM->SetTitle(";E_{#alpha}^{lab} at Si detector (MeV); E_{#alpha}^{cm} at elastic reaction (MeV)");
    //c->cd();function_EaMeas_EaCM->Draw("A*");
}
