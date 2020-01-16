// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__BLFITMODEL__H
#define __BAT__BLFITMODEL__H

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
#include <ParameterFactory.h>
#include <Parameter.h>
#include <DataFactory.h>
#include <EnergyResponse.h>

// This is a BLFitModel header file.
// Model source code is located in file BLFit/BLFitModel.cxx

using namespace std;

// ---------------------------------------------------------
class BLFitModel : public BCModel
{

public:

    // Constructor
    BLFitModel(const std::string& name);

    // Destructor
    ~BLFitModel();

    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double>& pars);


    // Overload LogAprioriProbability if not using built-in 1D priors
    // double LogAPrioriProbability(const std::vector<double> & pars);


    bool EditParameter(Parameter* param);

    void StackingDataPoints(double ** _E_Mean, double ** _Data, double ** _Data_e);
    
    double GetBackgroundCountsROI(const std::string& isotope_name, std::vector<double> pars);

    double GetDataCountsROI();

    std::vector<string> GetFittedRateParameter() {
        std::vector<string> fitted_param;
        for (int j=0; j < fitParameterVec.size(); j++) {
            if (fParameterInfoMap.find(fitParameterVec[j]) != fParameterInfoMap.end()) {
                fitted_param.push_back(fitParameterVec[j]);
            }
        }
        return fitted_param;
    }

    int GetFitParameterIndex(string histName) {
        return fParameterInfoMap[histName]->GetI(histName);
    }

    std::vector<double> GetFitParameterScale() {
        double low_rate_scale = 1E-6;
        double high_rate_scale = 1E-3;
        double energy_scale = 1.0;
        std::vector<double> scaleVector;
        for (int i=0; i<fittingVector.size(); i++) {
            if (fittingVector[i].find("Rate")!=string::npos) {
                if ((fittingVector[i] == "Rate_Xe136_2nu_XeLS") ||
                    (fittingVector[i] == "Rate_Kr85_XeLS")      ||
                    (fittingVector[i] == "Rate_U238_S2_XeLS")   ||
                    (fittingVector[i] == "Rate_Th232_S1_XeLS")  ||
                    (fittingVector[i] == "Rate_Bi210_XeLS")     ||
                    (fittingVector[i] == "Rate_Kr85_KamLS")     ||
                    (fittingVector[i] == "Rate_U238_S2_KamLS")
                    ) {
                    scaleVector.push_back(high_rate_scale);
                } else {
                    scaleVector.push_back(low_rate_scale);
                }
            } else {
                scaleVector.push_back(energy_scale);
            }
        }
        // std::vector<double> fitted_param;
        // for (int j=0; j < fitParameterVec.size(); j++) {//Rate Parameter
        //     if (fParameterInfoMap.find(fitParameterVec[j]) != fParameterInfoMap.end()) {
        //         fitted_param.push_back(rate_scale);

        //     }
        // }
        // for (int j=0; j < NofEResParam; j++) fitted_param.push_back(energy_scale);//Energy Scale
        // fitted_param.push_back(energy_scale);//Scaling
        if (scaleVector.size() != GetNParameters()) {
            cerr<<"parameter size ("<< scaleVector.size()<<" , "<<GetNParameters()<<") does not agree!"<<endl;
            abort();
        }
        return scaleVector;

    }

    int GetNParameter() {
        return param_id;
    }

    ifstream readParameterFile();

    double Spectrum();

    TH1D* GetDataHist(int rbin=30, int thetabin=-1);

    double Get0vbbCL() {
        if (fParameterInfoMap.find("Rate_Xe136_0nu_XeLS") == fParameterInfoMap.end()) return 0;
        if (GetParameter("Rate_Xe136_0nu_XeLS").Fixed()) return GetParameter("Rate_Xe136_0nu_XeLS").GetFixedValue();
        return GetMarginalized("Rate_Xe136_0nu_XeLS").GetLimit(0.9);
    }

    void UpdateParameterLimits() {
        const std::vector<double> pars = GetBestFitParameters();
        const std::vector<double> pars_err = GetBestFitParameterErrors();
        for(int i=0; i<GetNParameters(); i++) {
            double mean, variance;
            if (GetParameter(i).Fixed()) {
                mean = GetParameter(i).GetFixedValue();
                variance = 1.0;
            } else {
                mean = pars[i];
                variance = pars_err[i];
            }
            if (isnan(variance) || (!isfinite(variance)) || variance == 0.0) {
                variance == 1.0;
            }
            cout<<"mean: "<<mean<<"   "<<variance<<endl;
            GetParameter(i).SetLimits(std::max(mean - 5*variance, 0.0), mean + 5*variance);
        }
    }

    std::map<std::string, TH1D*> GetMCHist(std::vector<double>& pars, int rbin=99, int thetabin=-1);

    std::map<std::string, TH1D*> GetMCHist_Isotope(std::vector<double>& pars, int rbin=99, int thetabin=-1);

    std::vector<double> GetInitialPosition() { return initialPosition; }

    void BuildParameterGraph();


private:
    std::vector<string> fitParameterVec;
    std::vector<string> fittingVector;
    std::vector<string> fittingVectorRate;
    std::map<std::string,std::vector<string> > observableMap;
    std::vector<string> observableVec;
    std::vector<double> initialPosition;

    // std::map<std::string,std::vector<double> > fitXMap;
    // std::map<std::string,std::vector<double> > fitYMap;
    // std::map<std::string,std::vector<double> > fitZMap;
    // std::map<std::string,TH1D*> fit_map;
    // std::map<std::string,TH1D*> fit_map;
    int data_size;
    Int_t histBins;
    double histMin;
    double histMax;
    int param_id;

    int _nbin = 400;
    double _min = 0;
    double _max = 20.0;
    double _bin_width = (_max - _min) / double(_nbin);

    //Energy Binning
    const int _NofEnergyBin = 90;
    double _EnergyMin = 0.5;
    double _EnergyMax = 5.0;
    double _EnergyBinWidth = (_EnergyMax - _EnergyMin) / double(_NofEnergyBin);

    double _EnergyThreshold;
    double _UpperEnergyThreshold;


    ParameterFactory* Param;
    DataFactory* Data;
    EnergyResponse* ERes;
    Parameter* currentFitParameter;//dummy variable
    Parameter* buildParameter;//dummy variable

    std::map<std::string, Parameter*> fParameterInfoMap;
    TGraph** fParameterGraph[100];

    static const double _Min_kB = 0.00;
    static const double _Max_kB = 0.70;
    static const int _Bin_kB = 70;

    static const double _Min_R = 0.00;
    static const double _Max_R = 0.25;
    static const int _Bin_R = 50;

    static const int NBIN_kB = _Bin_kB + 1;
    static const int NBIN_R = _Bin_R + 1;
    static const int NBIN = NBIN_kB * NBIN_R;

    //Energy Response Param:
    static const int NofEResParam = 6;
    string EResParameterName[NofEResParam] = {"Alpha_Internal", "kB_Internal", "R_Internal", "Alpha_External", "kB_External", "R_External"};

    Parameter* currentEnergyResponse[20];

    bool fakeData = false;

    double ** _E_Mean;
    double ** _Data;
    double ** _Data_e;

    double * _E_Mean_fake[1000];
    double * _Data_fake[1000];
    double * _Data_e_fake[1000];

};
// ---------------------------------------------------------

#endif
