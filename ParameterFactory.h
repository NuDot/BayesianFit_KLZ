// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef _BFIT_PARAMETER_FACTORY
#define _BFIT_PARAMETER_FACTORY

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <cmath>
#include "TSystem.h"
#include <ParameterFactory.h>
#include <Parameter.h>

using namespace std;

// ---------------------------------------------------------
class ParameterFactory
{

public:

    // Constructor
    ParameterFactory();

    // Destructor
    ~ParameterFactory();

    bool ReadModelParameters();

    bool ReadGPParameters();

    //Start assigning fit parameters
    //fit parameters include two types of parameters:
    //   GPParameter: Parameter with a Gaussian Prior, it means that we have some previous knowledge about this parameter
    //   UPParameter: Parameter with an Uniform Prior, it means that we do not have any prior knowledge about this parameter,
    //                therefore we let it vary with uniform possibility throughout the given range.
    bool ReadFitParameters();

    bool ReadLivetimeParameters();

    bool ReadPhysicsParameters();

    Parameter* GetModelParameter(const std::string& name) {
        return gModelParameter[name];
    }

    Parameter* GetFitParameter(const std::string& name) {
        return gFitParameter[name];
    }

    std::vector<string> GetOrderedFitList() { return orderedRateName; };

    bool isFixedAll() {
        return _FixParameters_All;
    }

    bool CheckEnergyResponse(const std::string& type, const std::string& position) {
        if (type == "scale" && position == "internal") {
            return _FixEnergyScale_Internal;
        } else if (type == "scale" && position == "external") {
            return _FixEnergyScale_External;
        } else if (type == "nonlinearity" && position == "internal") {
            return _FixEnergyNonlinearity_Internal;
        } else if (type == "nonlinearity" && position == "external") {
            return _FixEnergyNonlinearity_External;
        } else {
            cerr << "ERROR : Unknown argument is given to checkEnergyResponse(type or position)" << endl;
            abort();
        }
    }

    Parameter* GetEnergyResponse(const std::string& name) {
        return gEnergyResposeParameter[name];
    }

    Parameter* setupEnergyResponse(Parameter* infoParameter, const std::vector<double>& params);

    std::vector<double>& GetEnergyResponseSpectrum(const std::string& name) { return gEnergyResponseMap[name]; };

    // bool isFixedEnergyScale_Internal() {
    //     return _FixEnergyScale_Internal;
    // }

    // bool isFixedEnergyNonlinearity_Internal() {
    //     return _FixEnergyNonlinearity_Internal;
    // }

    // bool isFixedEnergyScale_External() {
    //     return _FixEnergyScale_External;
    // }

    // bool isFixedEnergyNonlinearity_External() {
    //     return _FixEnergyNonlinearity_External;
    // }

    // vector<string> GetOrderedNameArray() {
    //     return orderedRateName;
    // }

    // //Read in required parameters from input .dat file or Monte Carlo tree.
    // bool readFVFile();

    // double Spectrum();

    // std::map<std::string,std::vector<TH1D*> >* GetEnergyMap() { return &fitEnergyMap;}

    // std::map<std::string,TH1D*>* GetDataMap() { return &fitDataMap; }

    // std::map<std::string, int>* GetParamMap() { return &parameterIndexMap; }


private:

    std::string Parameters_name;
    std::string FitParameters_name;
    std::string GPParameters_name;
    std::string EnergySpectrumDir;
    std::string Livetime_name;
    std::string RunList_name;
    std::string BiPo_name;
    std::string ChiSquare_BirksCherenkov_name;

    //dummy variables:
    double x[100000], y[100000];
    int N;
    Parameter * buildParameter;

    char buffer[256];
    double val, val_e;
    double dummy, obs_energy;
    string buf1, buf2, buf3, buf4, buf5, buf6, buf7;

    //List of Model Parameters:
    int _StartRun;
    int _LastRun;
    double _LiveTime;
    double _NumberOfBiPo214;
    double _R_Lower;
    double _R_Upper;
    double _Z_Lower;
    double _Z_Upper;
    double _DetectionEfficiency;
    double _DeltaNormalize;
    double _DeltaEnergy;

    //List of Gaussian Prior Parameters(GPParameter):
    static const int _NofGPParameter_Max = 100;
    int _NofGPParameter = 0;
    string _GPParameter[_NofGPParameter_Max];
    double _RateGPParameter[_NofGPParameter_Max];
    double _RateGPParameter_e[_NofGPParameter_Max];
    bool _IsGPParameter[_NofGPParameter_Max];

    //List of Fit Parameters:
    int _NofFitParameters = 0;

    bool _FixParameters_All = true;
    bool _FixEnergyScale_Internal = false;
    bool _FixEnergyNonlinearity_Internal = false;
    bool _FixEnergyScale_External = false;
    bool _FixEnergyNonlinearity_External = false;

    double _InitialAlpha_Internal; //overall normalization
    double _InitialkB_Internal; //birk's constant
    double _InitialR_Internal; //cherenkov/scintillation rate

    //detector response for external background
    double _InitialAlpha_External;
    double _InitialkB_External;
    double _InitialR_External;

    std::map<std::string,Parameter*> gModelParameter;
    std::map<std::string,Parameter*> gFitParameter;
    std::map<std::string,Parameter*> gEnergyResposeParameter;
    std::map<std::string,std::vector<double> > gEnergyResponseMap;

    vector<string> orderedRateName;

    // static const double _Min_kB = 0.00;
    // static const double _Max_kB = 0.70;
    // static const int _Bin_kB = 70;

    // static const double _Min_R = 0.00;
    // static const double _Max_R = 0.25;
    // static const int _Bin_R = 50;

    // static const int NBIN_kB = _Bin_kB + 1;
    // static const int NBIN_R = _Bin_R + 1;
    // static const int NBIN = NBIN_kB * NBIN_R;

    // double _kB_Nonlinearity[NBIN];
    // double _R_Nonlinearity[NBIN];

    // int _N_kB_Internal;
    // int _N_R_Internal;
    // double _r_kB_Internal;
    // double _r_R_Internal;
    // int _N_kB_External;
    // int _N_R_External;
    // double _r_kB_External;
    // double _r_R_External;

    // double _Alpha_Tagged_Bi214[NBIN];
    // double _ChiSquare_Tagged_Bi214[NBIN];




};
#endif
