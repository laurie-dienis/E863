#include "TROOT.h"
#include "TString.h"

#include <iostream>
#include <string>

void Print(const std::string& beam, const std::string& target, const std::string& light, double e_beam_ini,
           int pressure, const std::string& what = "")
{
    std::cout << "···· Runner ····" << '\n';
    std::cout << "-> Beam   : " << beam << '\n';
    std::cout << "-> Target : " << target << '\n';
    std::cout << "-> Light  : " << light << '\n';
    std::cout << "-> EBeam  : " << e_beam_ini << " MeV / u" << '\n';
    std::cout << "-> P      : " << pressure << " mbar" << '\n';
    std::cout << "-> What   : " << what << '\n';
}

void Runner(TString what = "plot")
{
    std::string beam {"40Ar"};
    std::string target {"4He"};
    std::string light {"4He"};
    double ebeam_i {1.7}; //MeV/u
    int pressure {450}; // mbar

    auto args {TString::Format("(\"%s\", \"%s\", \"%s\")", beam.c_str(), target.c_str(), light.c_str())};
    TString path {"./Pipelines/"};
    TString func {};
    TString ext {".cxx"};
    if(what.Contains("1"))
    {
        func = "Pipe1_Cal";
        gROOT->LoadMacro(path + func + ext);
        gROOT->ProcessLine(func + "()");
    }

    if(what.Contains("2_main"))
    {
        func = "Pipe2_main_PID";
        gROOT->LoadMacro(path + func + ext);
        gROOT->ProcessLine(func + "()");
    }
    
    if(what.Contains("2_second"))
    {
        func = "Pipe2_second_PID";
        gROOT->LoadMacro(path + func + ext);
        gROOT->ProcessLine(func + "()");
    }

    if(what.Contains("3_main"))
    {
        func = "Pipe3_main_Ecm";
        gROOT->LoadMacro(path + func + ext);
        gROOT->ProcessLine(func + "()");
    }
    
    if(what.Contains("3_second"))
    {
        func = "Pipe3_second_Ecm";
        gROOT->LoadMacro(path + func + ext);
        gROOT->ProcessLine(func + "()");
    }

}
