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

//***************************************
//Extraction function for SRIM Stop. Pow.
//***************************************
int extraction_data_SRIM(const Char_t *file, Double_t EnergyPart[], Double_t RangePart_in_Au[], Double_t SPelPart_in_Au[], Double_t SPnuPart_in_Au[]) {
    std::ifstream in;
    in.open(file);
    Int_t nlines = 0;
    while (1) {
        in >> EnergyPart[nlines] >> SPelPart_in_Au[nlines] >> SPnuPart_in_Au[nlines] >> RangePart_in_Au[nlines]; // SP in keV/microns, E in keV, Range in Angstrom
    	//std::cout<< EnergyPart[nlines] <<" E  "<<std::endl;
        if (!in.good()) break;
        nlines++;
    }
    in.close();
    return 1;
}


//***************************************
//Interpolation on SRIM Stop. Pow.
//***************************************
double interpol(double x[], double y[], int const n, double p)
{
    int i = 1;
    double val, b1, c1, d1, c2, d2;
    // b1 false if something is wrong
    bool bb1 = true;
    // loop breaking
    bool bb2 = true;
    if (p <= x[1])
    {
        bb1 = false;
    }
    else
    {
        do
        {
            i++;
            if (p <= x[i])bb2 = false;
        } while (bb2 && (i < n - 2));
    }
    if (p > x[n - 2])bb1 = false;
    if (bb1 == true)
    {
        b1 = (y[i - 1] - y[i - 2]) / (x[i - 1] - x[i - 2]);
        b1 = b1 * p - b1 * x[i - 2] + y[i - 2];
        c1 = (y[i] - y[i - 2]) / (x[i] - x[i - 2]);
        c1 = c1 * p - c1 * x[i - 2] + y[i - 2];
        d1 = (y[i + 1] - y[i - 2]) / (x[i + 1] - x[i - 2]);
        d1 = d1 * p - d1 * x[i - 2] + y[i - 2];
        c2 = (c1 - b1) / (x[i] - x[i - 1]);
        c2 = c2 * p - c2 * x[i - 1] + b1;
        d2 = (d1 - b1) / (x[i + 1] - x[i - 1]);
        d2 = d2 * p - d2 * x[i - 1] + b1;
        val = (d2 - c2) / (x[i] - x[i - 1]);
        val = val * p - val * x[i] + c2;
    }
    else
    {
        val = 0.;
    }
    return val;
}


//***************************************
// Local derivation of energy losses
// energy in MeV, distance in mm,
// Sto. Pow. in keV/microns
//***************************************
double loss_E_srim(double Ei, double distance, Double_t EnergyPart[], Double_t SPel_Part[], Double_t SPnu_Part[], int N_energy) {
    //distance in mm
    double Etemp = Ei;
    double dpart = 0;
    double ddx = 1; //micron
    while (dpart <= distance*1000) {
        Etemp = Etemp - 0.001 * interpol(EnergyPart, SPel_Part, N_energy, Etemp * 1000.0) * ddx - 0.001 * interpol(EnergyPart, SPnu_Part, N_energy, Etemp * 1000.0) * ddx;
        //std::cout<< interpol(EnergyPart, SPel_Part, N_energy, Etemp * 1000.0)  <<" interpol  "<<std::endl;
        //std::cout<< Etemp * 1000.0  <<" interpol  "<<std::endl;
        dpart = dpart + ddx;
    }
    //std::cout <<" \n  "<<std::endl;
    //std::cout<< Etemp <<" Etemp  "<<std::endl;
    return Ei - Etemp;
}

void TOF_Ar_Fe() {
    
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
    int iter=10000;
    int Ghoix=1; // Ghoix = 1 : p(E) impulsion as a function of the energy 
				 // Ghoix = 2 : angular distribution for 6Li and p 
				 // Ghoix = 3 : E(6Li) as a function of Ep 
				 // Ghoix = 4 : E(Li) energy as a function of the angle (MeV)
				 // Ghoix = 5 : E(p) energy as a function of the angle (MeV)

    
    // variables
    double Elab,Ee,Qreac,Qreac2,target_thickness;
    
    int extraction_data;
    int N_energy = 132;
    const Char_t *fileSRIM ={"sp_srim/H_Si.dat"};
    Double_t EnergyRecoil[132];     
    Double_t RangeRecoil_in_Si[132];    
    Double_t SPelRecoil_in_Si[132];    
    Double_t SPnuRecoil_in_Si[132];
    extraction_data = extraction_data_SRIM(fileSRIM, EnergyRecoil, RangeRecoil_in_Si,SPelRecoil_in_Si, SPnuRecoil_in_Si);

    // Reaction p(6Li,p')6Li*
        Elab=23;  //1.8 MeV/u 27MeV
		Ee=0; 
        ma=mfe56;
        mA=mar40;
        mB=mtc95+Ee;
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
    TF1 *data = new TF1("data","[0]*x^[1]", 2, 18);
    data->SetParameter(0, 118.7339);
    data->SetParameter(1, -0.500555);		
    data->SetLineColor(kViolet);	
    TF1 *data2 = new TF1("data2","[0]+x*[1+x*x[2]", 2, 8.3);
    data2->SetParameter(0, -32.48);
    data2->SetParameter(1, 9.0437);		
    data2->SetParameter(2, -0.401009);		
    data2->SetLineColor(kViolet);	

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
        
        double ELoss = loss_E_srim(Eb, 0.5, EnergyRecoil, SPelRecoil_in_Si, SPnuRecoil_in_Si, N_energy);
        //cout << "\n ELoss " << ELoss*1000 << " MeV" ;
        double Eb = ELoss;


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
    h1->Fit("fit");
    cout << "\n" ;
    h1->GetXaxis()->SetTitle("E (MeV)");
    h1->GetYaxis()->SetTitle("TOF (ns)");
    h1->SetMarkerColor(kBlack);
    h1->Draw("colz");
    data->Draw("same");
    data2->Draw("same");
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
