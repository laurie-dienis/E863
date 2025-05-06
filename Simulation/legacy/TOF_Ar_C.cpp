//
// Created by De Oliveira francois
// version 05/10/2022
//
// see example   https://root.cern.ch/doc/master/PhaseSpace_8C.html
//
// Reaction  a(A,b)B
// Reaction p(6Li,p')6Li*
// Reaction B --> S + s
// Reaction 6Li* --> alpha + d
// to run :  root 6Li.cpp
//


#include "mass.hpp"


void TOF_Ar_C() {
    
    // *******************************************
    // variables
    // *******************************************
    // Elab = incident beam energy in MeV
    // iter = number of points simulated
	// Reaction  a(A,b)B
	// Reaction p(6Li,p')6Li*
	// Reaction B --> S + s
	// Reaction 6Li* --> alpha + d
    // *******************************************
    int iter=1000000;
    int Ghoix=1; // Ghoix = 1 : p(E) impulsion as a function of the energy 
				 // Ghoix = 2 : angular distribution for 6Li and p 
				 // Ghoix = 3 : E(6Li) as a function of Ep 
				 // Ghoix = 4 : E(Li) energy as a function of the angle (MeV)
				 // Ghoix = 5 : E(p) energy as a function of the angle (MeV)

    
    // variables
    double Elab,Ee,Qreac,Qreac2,target_thickness;

    // Reaction p(6Li,p')6Li*
        Elab=23;  //1.8 MeV/u 27MeV
		Ee=0; 
        ma=mc12;
        mA=mar40;
        mB=mv51+Ee;
        mb=mp;
		
	target_thickness=85; // Gold(Au) target, in micrometer 

    TCanvas *gr = new TCanvas;
    
	// Reaction p(6Li,p')6Li*

    //variables
    double x,betas,EB,Eb,ES,Es;
	

    //(Momentum, Energy units are in Gev/C, GeV)
    ma /= 1000.;
    mA /= 1000.;
    mB /= 1000.;
    mb /= 1000.;
    Elab /= 1000.;           
		
    //Random generator 
	std::mt19937_64 rng;
	std::uniform_real_distribution<double> unif(1,1.1);
    
    //random for E lab 6MeV/u until 0.35 MeV/u
	std::mt19937_64 rng2;
	std::uniform_real_distribution<double> unif2(1,239);
    
	// definition of the 2D and 1D spectra
    TH2F *h1 = new TH2F("h1","",800,1,20,800,0,109); // impulsion as a function of the energy

    double projB,projb,projs,projS;
    TLorentzVector ppx, difas;
    TVector3 imps,impbeam,impx;
    
// *******  LOOP ******************************************
    for (Int_t n=0;n<iter;n++) {
		
		
		//// Reaction p(6Li,p')6Li* ////
		
		TGenPhaseSpace event1;
    
		// final configuration 2 particles B+b
		Double_t masses1[2] = {mB, mb} ;
    

		// Lorentz vector for the target
		TLorentzVector target(0.0, 0.0, 0.0, ma);
    
		// beam = A
		// total energy of the beam  Etb = mc^2 + Ec
        Elab = unif2(rng2)/1000.;
		//cout << "\n Ebeam " << Elab << " MeV" ;
		double Etb1=mA+Elab;
		// impulsion beam
		double px1=sqrt(Etb1*Etb1-mA*mA);  // p of A
		// Lorentz vector for the beam
		TLorentzVector beam1;
		beam1.SetPxPyPzE(0,0,px1,Etb1); // Px, px, pz, E

    
		// entrance Lorentz vector
		TLorentzVector W1 = beam1 + target;
    
    
		// generate decay in 2 particles of mass mB and mb
		event1.SetDecay(W1, 2, masses1);
	
        Double_t weight1 = event1.Generate();
        TLorentzVector *pB = event1.GetDecay(0);
        TLorentzVector *pb = event1.GetDecay(1);
        TLorentzVector ppB = *pB ;
        TLorentzVector ppb = *pb ;
        // energies 
        EB=1000.*(ppB.E()-mB); // Ec 
		//cout << "\n E6Li is equal to " <33.12< EB << " MeV" ;
        Eb=1000.*(ppb.E()-mb);
		//cout << "\n Ealpha is equal to " << Eb << " MeV" ;
        // angles
        projB = ppB.Theta()*57.3; // angle
        projb = ppb.Theta()*57.3; // angle
		
		
		//taking into account the target thickness 
		
		//for the alpha
		//double rangeS_initial= 4.828*Eb-32.52;
		//double rangeS_final= rangeS_initial - target_thickness;
		//Eb= (rangeS_final+32.52)/4.828;
		//cout << "\n Ealpha(f) is equal to " << ES << " MeV" ;
        
        double vb = sqrt((2*Eb)/(mb*1000.))*30; //velocity in cm/ns
        //cout << "\n v alpha is equal to " << vb << "cm/ns" ;
        double vB = sqrt((2*EB)/(mB*1000.))*30; //velocity in cm/ns
        //cout << "\n v 15O is equal to " << vB << "cm/ns" ;

        double TOF_b = 177/vb + unif(rng)-4.4;  //28.6 is the distance for which we measure the TOF
		double TOF_B = 200/vB;
        //cout << "\n TOF alpha is equal to " << TOF_b << "ns" ;


// filling the spectra ************************************

        
		if (Ghoix == 1)
        {
            h1->Fill(Eb, TOF_b); // energy as a function of the angle (MeV)
            h1->SetMarkerColor(2);
        }
    }
// ***********  END LOOP ********************************
    
    
    
    

    
 //  Drawing the spectra ***********************************

	
	if (Ghoix == 1)
    {
    TF1 *fit = new TF1("fit","[0]*(x^[1])",10, -0.5);
    //protons
    cout << "protons"<<"\n" ;
    //h1->Fit("fit");
    cout << "\n" ;
    h1->GetXaxis()->SetTitle("E (MeV)");
    h1->GetYaxis()->SetTitle("TOF (ns)");
    h1->SetMarkerColor(kBlack);
    h1->Draw("colz");
    }
	
    
    // text and lines  *********************************
	
    TLatex *t1 = new TLatex();
    t1->SetTextFont(42);
    t1->SetTextSize(0.05);
    t1->SetTextAngle(0.0);
    
    TLatex *t2 = new TLatex();
    t2->SetTextFont(42);
    t2->SetTextSize(0.05);
    t2->SetTextAngle(90.0);
    
    TLatex *t3 = new TLatex();
    t3->SetTextFont(42);
    t3->SetTextSize(0.04);
    t3->SetTextAngle(0.0);
    
    TLine *line = new TLine();
    
    
    
    if (Ghoix == 1)  // impulsion as a function of the energy
    {
    }
    
    
    
    if (Ghoix ==2) // angular distribution
    {

    }
  
    
    
     if (Ghoix ==3)  // Ed as a function of Ea
     {
         t1->DrawLatex(22,-3,"E^{lab}_{t}(MeV)");
         t2->DrawLatex(-3.5,27,"E^{lab}_{#alpha}(MeV)");
//         t1->DrawLatex(5,20,"#Theta^{lab}_{#alpha_{1}}= 60#circ");
//         t1->DrawLatex(11,20,"#Theta^{lab}_{#alpha_{2}}= 73#circ");
     }

	 if (Ghoix == 4)  // energy #alpha as a function of the angle 
    {
        t1->DrawLatex(13.5,3,"#Theta_{#alpha}(deg)");
        t2->DrawLatex(-0.7,27,"E^{lab}_{#alpha}(MeV)");
    }
	
	if (Ghoix == 5)  // energy d as a function of the angle 
    {
        t1->DrawLatex(19,0.5,"#Theta_{d}(deg)");
        t2->DrawLatex(-1.5,15,"E^{lab}_{d}(MeV)");
    }
	
	
	
	
}
