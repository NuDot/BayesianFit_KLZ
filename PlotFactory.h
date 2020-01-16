// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__PLOT__H
#define __BAT__PLOT__H

#include <BAT/BCModel.h>
#include <BAT/BCGaussianPrior.h>
#include <string>
#include <vector>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>  
#include <iostream>
#include <sstream>
#include <map>

// This is a BLFitModel header file.
// Model source code is located in file BLFit/BLFitModel.cxx

using namespace std;

// ---------------------------------------------------------
class PlotFactory
{

public:

    // Constructor
    PlotFactory(std::map<std::string, TH3D*> map, TH3D* data, std::vector<double> param, int rbt, int tbt);

    void PlotROI(int rbin, int thetabin) {
        EnergyLow = 2.35;
        EnergyHigh = 2.7;
        logy = false;
        out_name = "hist_zoom";
        countLow = 0.01;
        countHigh = 10.0;
        SumHist(rbin, thetabin);
        plot(rbin, thetabin);
    }


    void PlotAll(int rbin, int thetabin) {
        EnergyLow = 0.5;
        EnergyHigh = 4.8;
        logy = true;
        out_name = "hist";
        countLow = 0.01;
        countHigh = 100000.0;
        SumHist(rbin, thetabin);
        plot(rbin, thetabin);
    }

    // Destructor
    ~PlotFactory();


private:

    void plot(int rbin, int thetabin);
    void SumHist(int rbin, int thetabin);


    int LINE_WIDTH=2;

    //Energy Binning
    int _NofEnergyBin = 90;
    double _EnergyMin = 0.5;
    double _EnergyMax = 5.0;
    double _EnergyBinWidth = (_EnergyMax - _EnergyMin) / double(_NofEnergyBin);

    TH3D* hist_data_all;
    std::map<std::string, TH3D*> MCMap_all;

    TH1D* hist_data;
    std::map<std::string, TH1D*> MCMap;
    std::vector<double> paramVector;

    double EnergyLow = 0.5;
    double EnergyHigh = 4.8;
    double countLow = 0.01;
    double countHigh = 100000.0;
    bool logy = true;
    string out_name = "hist";

    bool initialized = false;

    int bin_total;
    int thetabin_total;

};
// ---------------------------------------------------------

#endif
