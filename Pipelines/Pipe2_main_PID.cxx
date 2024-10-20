#ifndef Pipe2_PID_cxx
#define Pipe2_PID_cxx

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TROOT.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <string>

void Pipe2_main_PID() {
      // Check if the cut file exists
    const char *cutFileName = "Cuts/cut_alpha_1.root";
    TFile *cutFile = TFile::Open(cutFileName);
    TCutG *cutG = nullptr;

    if (cutFile && !cutFile->IsZombie()) {
        // Try to read the graphical cut
        cutFile->GetObject("alpha_cut", cutG);
        cutFile->Close();
        delete cutFile;
    }

    // Test if the graphical cut is found
    if (!cutG) {
        std::cerr << "Error: Graphical cut 'cutG' not found in file 'cut_alpha_1.root'."
                  << " Plotting PID without cut." << std::endl;
    }

    // List of file names
    std::vector<std::string> fileNames = {
        "RootFiles/Cal/Cal.root",
    };
    const char *treeName = "DataTree"; // Tree name (common to all files)

    // Variables to hold leaf data
    double tac;
    double si_cal;
    double faster_time;

    // Create histograms
    TH2F *PID = new TH2F("PID", "PID", 300, 0, 20, 300, 0, 109);
    TH2F *PID_faster = new TH2F("PID_faster", "PID_faster", 500, 200, 800, 3000, 0, 20);

    // Create a new ROOT file only if the cut exists
    TFile *newFile = nullptr;
    TTree *PIDTree = nullptr;
    
    if (cutG) {
        newFile = TFile::Open("RootFiles/PID/alpha_1.root", "RECREATE");
        if (!newFile || newFile->IsZombie()) {
            std::cerr << "Error creating new file: RootFiles/PID/alpha_1.root" << std::endl;
            if (newFile) newFile->Close();
            delete newFile;
            delete PID;
            return;
        }

        // Create a new tree
        PIDTree = new TTree(treeName, "alpha Data Tree");

        // Set branch addresses for the new tree
        PIDTree->Branch("Si_cal_1", &si_cal, "Si_cal_1/D");
    }

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
            delete file;
            continue;
        }

        // Set branch addresses for reading
        tree->SetBranchAddress("Si_cal_1", &si_cal);
        tree->SetBranchAddress("TAC_cal_1", &tac);
        tree->SetBranchAddress("faster_time_1", &faster_time);

        // Loop over the tree entries and apply the graphical cut if it exists
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            tree->GetEntry(i);

            // Always fill the PID histogram
            PID->Fill(si_cal, tac);
	    PID_faster->Fill(si_cal, faster_time);

            // Fill the new tree only if the cut exists and the point is inside the cut
            if (cutG && cutG->IsInside(tac, si_cal)) {
                PIDTree->Fill();
            }
        }

        // Close the input file and clean up
        file->Close();
        delete file;
    }

    // Plot the PID histogram
    auto *c20 = new TCanvas("c20", "Pipe2 canvas 0");
    PID->GetYaxis()->SetTitle("TOF (ns)");
    PID->GetXaxis()->SetTitle("E_Si (MeV)");
    PID->Draw("colz");

    // TF1 *line1 = new TF1("line1","261.7*x^(-0.4844)+10",0, 20);
    // line1->SetLineColor(kRed);
    // line1->Draw("same");
    
    // TF1 *line2 = new TF1("line2","534.226*x^(-0.493824)-200",0, 20);
    // line2->SetLineColor(kBlue);
    // line2->Draw("same");

    //////////////Al dist deg-si 575 /////////////
    // //protons
    // TF1 *line1 = new TF1("line1","(413.578*x^(-0.494904))-27",11.11, 20);
    // line1->SetLineColor(kCyan);
    // line1->SetLineWidth(4);
    // line1->Draw("same");
    
    // //alpha
    // TF1 *line2 = new TF1("line2","(2481.*x^(-1.22641))+5",11.3, 20);
    // line2->SetLineColor(kViolet);
    // line2->SetLineWidth(4);
    // line2->Draw("same");
    
    //////////////Al in col //////////////////
    // //protons
    // TF1 *line1 = new TF1("line1","(26.8864*x^(-0.439753))+62",13.3, 20);
    // line1->SetLineColor(kCyan);
    // line1->SetLineWidth(4);
    // line1->Draw("same");
    
    // //alpha
    // TF1 *line2 = new TF1("line2","(53.4505*x^(-0.463471))+56",11.3, 20);
    // line2->SetLineColor(kViolet);
    // line2->SetLineWidth(4);
    // line2->Draw("same");
    
    ///////Al dist deg-si 225 ///////////////
    //alphas
    // TF1 *line1 = new TF1("line1","(-3.79063+x*3.76577+x*x*(-0.0706469))+25",5.2, 20);
    // line1->SetLineColor(kBlue);
    // line1->SetLineWidth(4);
    // line1->Draw("same");
    
    //alpha
    TF1 *line2 = new TF1("line2","(329.41*x^(-0.521455))",11.9, 20);
    line2->SetLineColor(kViolet);
    line2->SetLineWidth(4);
    line2->Draw("same");
    
    ////////Ti///////////////
    // //protons
    TF1 *line3 = new TF1("line3","22.1702*x^(-0.722957)+70",2.3, 7.4);
    line3->SetLineColor(kOrange);
    line3->SetLineWidth(4);
    line3->Draw("same");
    
    // //alpha
    // TF1 *line4 = new TF1("line4","45.6722*x^(-0.592348)",0, 20);
    // line4->SetLineColor(kOrange+7);
    // line4->SetLineWidth(4);
    // line4->Draw("same");
    

    //////////////He gas////////////
    TF1 *line5 = new TF1("line5","18.34*x^(-0.374257)",0, 20);
    line5->SetLineColor(kRed);
    line5->SetLineWidth(4);
    line5->Draw("same");

    // //alpha
    TF1 *line6 = new TF1("line6","41.85433*x^(-0.47297)",0, 20);
    line6->SetLineColor(kOrange+1);
    line6->SetLineWidth(4);
    line6->Draw("same");

    /////////Carbon/////////////
    // TF1 *line3 = new TF1("line3","132.548*x^(-0.553591)",9.03, 20);
    // line3->SetLineColor(kOrange);
    // line3->SetLineWidth(4);
    // line3->Draw("same");

    /////////Fer//////////
    TF1 *line4 = new TF1("line4","126.18*x^(-0.533189)",2.3, 7.9);
    line4->SetLineColor(kCyan);
    line4->SetLineWidth(4);
    line4->Draw("same");
    // TF1 *line5 = new TF1("line5","13.52*x^(0.672484)",0, 7.9);
    // line5->SetLineColor(kCyan);
    // line5->SetLineWidth(4);
    // line5->Draw("same");

    ////////Al-Al///////////
    TF1 *line7 = new TF1("line6","18.3180953284087*(x)^(0.706276084966189)",0, 20);

    line7->SetLineColor(kGreen);
    line7->SetLineWidth(4);
    line7->Draw("same");

    // Plot the cut if it exists
    if (cutG) {
        cutG->SetLineColor(kRed); // Set cut line color
        cutG->SetLineWidth(2);    // Set cut line width
        cutG->Draw("L same");     // Draw cut on top of the histogram
    }
    c20->SaveAs("PID_Si1.png");
    
    // auto *c21 = new TCanvas("c21", "Pipe2 canvas 0, PID faster");
    // PID_faster->GetYaxis()->SetTitle("faster time (ns)");
    // PID_faster->GetXaxis()->SetTitle("E_Si (MeV)");
    // PID_faster->Draw("colz");

    // Write the histogram and the new tree (if the cut exists) to the new file
    if (newFile) {
        newFile->cd(); // Ensure the new file is the current directory
        PID->Write();
        PIDTree->Write();

        // Close the new file
        newFile->Close();
        delete newFile;
    }
}
#endif
