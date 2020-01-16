#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <ParameterFactory.h>

using namespace std;

ParameterFactory::ParameterFactory() {
    if(
    getenv("DOUBLEBETA_ANALYSIS_DATA_DIR")==NULL
    ||
    getenv("DOUBLEBETA_ANALYSIS_ENERGY_SPECTRUM_DIR")==NULL
    ||
    getenv("DOUBLEBETA_ANALYSIS_LIVETIME_DIR")==NULL
    ||
    getenv("DOUBLEBETA_ANALYSIS_BIPO_DIR")==NULL
    ||
    getenv("DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_FILE_XELS")==NULL
    ||
    getenv("DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_FILE_KAMLS")==NULL
    ||
    getenv("DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_FILE_FILM")==NULL
    ||
    getenv("DOUBLEBETA_ANALYSIS_PARAMETER_DIR")==NULL
    ||
    getenv("DOUBLEBETA_ANALYSIS_RUNLIST_DIR")==NULL
    //||
    //getenv("DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_BY_PHOTON_FILE_FILM")==NULL
    ){
    cerr << "ERROR: Cannot set analysis parameters" << endl;
    cerr << "Please check environmental variable" << endl;
    cerr << "DOUBLEBETA_ANALYSIS_DATA_DIR: " << getenv("DOUBLEBETA_ANALYSIS_DATA_DIR") << endl;
    cerr << "DOUBLEBETA_ANALYSIS_ENERGY_SPECTRUM_DIR: " << getenv("DOUBLEBETA_ANALYSIS_ENERGY_SPECTRUM_DIR") << endl;
    cerr << "DOUBLEBETA_ANALYSIS_LIVETIME_DIR: " << getenv("DOUBLEBETA_ANALYSIS_LIVETIME_DIR") << endl;
    cerr << "DOUBLEBETA_ANALYSIS_BIPO_DIR: " << getenv("DOUBLEBETA_ANALYSIS_BIPO_DIR") << endl;
    cerr << "DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_FILE_XELS: " << getenv("DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_FILE_XELS") << endl;
    cerr << "DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_FILE_KAMLS: " << getenv("DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_FILE_KAMLS") << endl;
    cerr << "DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_FILE_FILM: " << getenv("DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_FILE_FILM") << endl;
    cerr << "DOUBLEBETA_ANALYSIS_PARAMETER_DIR: " << getenv("DOUBLEBETA_ANALYSIS_PARAMETER_DIR") << endl;
    cerr << "DOUBLEBETA_ANALYSIS_RUNLIST_DIR: " << getenv("DOUBLEBETA_ANALYSIS_RUNLIST_DIR") << endl;
    //cerr << "DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_BY_PHOTON_FILE_FILM: " << getenv("DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_BY_PHOTON_FILE_FILM") << endl;
    cerr << endl;

    abort();
    }
    Parameters_name = ((string) getenv("DOUBLEBETA_ANALYSIS_PARAMETER_DIR")) + "/Parameters.dat";
    FitParameters_name = ((string) getenv("DOUBLEBETA_ANALYSIS_PARAMETER_DIR")) + "/FitParameters_converge.dat";
    GPParameters_name = ((string) getenv("DOUBLEBETA_ANALYSIS_PARAMETER_DIR")) + "/PenaltyParameters.dat";
    Livetime_name = ((string) getenv("DOUBLEBETA_ANALYSIS_LIVETIME_DIR")) + "/livetime.dat";
    RunList_name = ((string) getenv("DOUBLEBETA_ANALYSIS_DATA_DIR")) + "/Single-DoubleBeta.list";
    BiPo_name = ((string) getenv("DOUBLEBETA_ANALYSIS_BIPO_DIR"));
    if (!ReadModelParameters()) {
        cerr << "ERROR : Failed reading model parameters!" << endl;
        abort();
    }
    if (!ReadLivetimeParameters()) {
        cerr << "ERROR : Failed reading livetime parameters!" << endl;
        //abort();
    }
    if (!ReadPhysicsParameters()) {
        cerr << "ERROR : Failed reading physics parameters!" << endl;
        //abort();
    }
    if (!ReadGPParameters()) {
        cerr << "ERROR : Failed reading GP parameters!" << endl;
        abort();
    }
    if (!ReadFitParameters()) {
        cerr << "ERROR : Failed reading fit parameters!" << endl;
        abort();
    }
}

ParameterFactory::~ParameterFactory() {
    //Destructor
}

bool ParameterFactory::ReadModelParameters() {
    //Setting model parameters
    ifstream fParameters(Parameters_name.c_str());
    if (!fParameters) {
        cerr << "ERROR : cannot open input file : Parameters" << endl;
        return false;
    }

    cerr << "=================================================================" << endl;
    cerr << "Set Model Parameters:" << endl;
    cerr << "=================================================================" << endl;
    N = 0;
    while (fParameters >> buffer >> val) {
        switch(N) {
            case 0: _StartRun = int(val + 0.1);
                    cerr << "StartRun = " << _StartRun << endl;
                    break;
            case 1: _LastRun = int(val + 0.1);
                    cerr << "LastRun = " << _LastRun << endl;
                    break;
            case 2: _R_Lower = val;
                    break;
            case 3: _R_Upper = val;
                    if (!(_R_Lower == 0 && _R_Upper == 0)) {
                        cerr << "R_Lower (m) = " << _R_Lower << endl;
                        cerr << "R_Upper (m) = " << _R_Upper << endl;
                    }
                    break;
            case 4: _Z_Lower = val;
                    cerr << "Z_Lower (m) = " << _Z_Lower << endl;
                    break;
            case 5: _Z_Upper = val;
                    cerr << "Z_Upper (m) = " << _Z_Upper << endl;
                    break;
            case 6: _DetectionEfficiency = val;
                    cerr << "DetectionEfficiency  = " << _DetectionEfficiency << endl;
                    break;
            case 7: _DeltaNormalize = val;
                    cerr << "DeltaNormalize = " << _DeltaNormalize << endl;
                    break;
            case 8: _DeltaEnergy = val;
                    cerr << "DeltaEnergy = " << _DeltaEnergy << endl;
                    break;
            default: {
                cerr << "ERROR : Parameters format error" << endl;
                return false;
            }
        }
        N++;
    }
    fParameters.close();
    cerr << "=================================================================" << endl;
    cerr << "Build Model Parameters:" << endl;

    buildParameter = new Parameter("run_range", "data_parameter", true);
    buildParameter->SetIValue("low", _StartRun);
    buildParameter->SetIValue("high", _LastRun);
    gModelParameter["run_range"] = buildParameter;

    buildParameter = new Parameter("z_fv", "data_parameter", true);
    buildParameter->SetDValue("low", _Z_Lower);
    buildParameter->SetDValue("high",_Z_Upper);
    gModelParameter["z_fv"] = buildParameter;

    buildParameter = new Parameter("r_fv", "data_parameter", true);
    buildParameter->SetDValue("low", _R_Lower);
    buildParameter->SetDValue("high",_R_Upper);
    gModelParameter["r_fv"] = buildParameter;

    buildParameter = new Parameter("DetectionEfficiency", "data_parameter", true);
    buildParameter->SetDValue("initial", _DetectionEfficiency);
    gModelParameter["DetectionEfficiency"] = buildParameter;

    buildParameter = new Parameter("DeltaNormalize", "data_parameter", true);
    buildParameter->SetDValue("initial", _DeltaNormalize);
    gModelParameter["DeltaNormalize"] = buildParameter;

    buildParameter = new Parameter("DeltaEnergy", "data_parameter", true);
    buildParameter->SetDValue("initial", _DeltaEnergy);
    gModelParameter["DeltaEnergy"] = buildParameter;
    cerr << "Finished Building Parameters!" << endl;
    

}

bool ParameterFactory::ReadGPParameters() {
    Parameter* buildParameter;
    ifstream fGPParameter(GPParameters_name.c_str());
    if (!fGPParameter) {
        cerr << "ERROR : cannot open input file : GPParameter" << endl;
        return false;
    }

    cerr << "=================================================================" << endl;
    cerr << "Set Gaussian Prior Parameters:" << endl;
    cerr << "=================================================================" << endl;

    _NofGPParameter = 0;

    while (fGPParameter.getline(buffer, sizeof(buffer))) {

        istringstream strout(buffer);

        if (!(strout >> buf1 >> buf2 >> buf3 >> buf4)) {
            cerr << "ERROR : GPParameter format error" << endl;
            return false;
        }

        if (_NofGPParameter >= _NofGPParameter_Max) {
            cerr << "ERROR : GPParameter format error (out of range)" << endl;
            return false;
        }

        cerr << buffer << endl;

        //Reading 4 fields from GPParameter.dat
        string parname = buf1;
        double rate = 0;
        double rate_e = 0;
        int isGP = 0;
        Parameter* physicsParam = gModelParameter["physics_constant"];

        //2 types of GPParameter are evaluated separately, these includes:
        //Parameter that can be estimated from past publication
        //Parameter that can be estimated from event tagging
        if (buf2 == "past") {
            if (parname == "Rate_Xe136_0nu_XeLS") {
                // estimation from past publication (initial value for fit)
                rate = physicsParam->GetD("Rate_Xe136_0nu_XeLS_Expected");
                rate_e = 0;

            } else if (parname == "Rate_Xe136_2nu_XeLS") {
                // estimation from past publication (initial value for fit)
                rate = physicsParam->GetD("Rate_Xe136_2nu_XeLS_Expected");
                rate_e = 0;
            } else {
                cerr << "ERROR : GPParameter format error (rate Xe136_2nu in XeLS)" << endl;
                return false;
            }
        } else if (buf2 == "tag") {
            if (parname == "Rate_U238_S2_XeLS") {
                // estimation from BiPo214
                rate = physicsParam->GetD("Rate_U238_S2_XeLS_Expected");
                rate_e = physicsParam->GetD("Rate_U238_S2_XeLS_Expected_e");
            } else if (parname == "Rate_Th232_S2_XeLS") {
                // estimation from Tl208 (tentative) --> "0" (no penalty)
                rate = physicsParam->GetD("Rate_Th232_S2_XeLS_Expected");
                rate_e = physicsParam->GetD("Rate_Th232_S2_XeLS_Expected_e");
            } else {
                cerr << "ERROR : GPParameter format error (rate U238_S2 in XeLS)" << endl;
                return false;
            }
        } else {
            if (!(istringstream(buf2) >> rate)) {
                cerr << "ERROR : GPParameter format error (rate)" << endl;
                return false;
            }

            if (!(istringstream(buf3) >> rate_e)) {
                cerr << "ERROR : GPParameter format error (rate error)" << endl;
                return false;
            }
        }

        if (!(istringstream(buf4) >> isGP)) {
            cerr << "ERROR : GPParameter format error (isGP flag)" << endl;
            return false;
        }

        if (!(isGP == 0 || isGP == 1)) {
            cerr << "ERROR : FitParameters format error (isGP)" << endl;
            return false;
        }

        _GPParameter[_NofGPParameter] = parname;
        _RateGPParameter[_NofGPParameter] = rate;
        _RateGPParameter_e[_NofGPParameter] = rate_e;


        //If the penalty parameter has 0 penalty, then it's not a penalty parameter
        if (isGP == 0) {
            _IsGPParameter[_NofGPParameter] = false;
        } else {
            _IsGPParameter[_NofGPParameter] = true;
        }
        _NofGPParameter++;
    }
    cerr << endl;

    fGPParameter.close();
    cerr << "=================================================================" << endl;
}

//Parameter Structure:
//  Name: name
//  Type: "fit_parameter"
//  PriorType: "GP" or "UP"
//  Unit: unit string
//  int value: color
//  double initial: initial value
//  double expected: expected rate
//  double error: expected error
//  double low: fit range lower bound
//  double high: fit range upper bound 

bool ParameterFactory::ReadFitParameters() {
    Parameter* buildParameter;
    ifstream fFitParameters(FitParameters_name.c_str());
    if (!fFitParameters) {
        cerr << "ERROR : cannot open input file : FitParameters" << endl;
        return false;
    }

    cerr << "=================================================================" << endl;
    cerr << "Set Fit Parameters:" << endl;
    cerr << "=================================================================" << endl;

    bool is_fixed_all = true;
    bool is_fixed_Alpha_Internal = false;
    bool is_fixed_kB_Internal = false;
    bool is_fixed_R_Internal = false;
    bool is_fixed_Alpha_External = false;
    bool is_fixed_kB_External = false;
    bool is_fixed_R_External = false;

    _NofFitParameters = 0;

    while (fFitParameters.getline(buffer, sizeof(buffer))) {
        istringstream strout(buffer);


        if (!(strout >> buf1 >> buf2 >> buf3 >> buf4 >> buf5 >> buf6 >> buf7)) {
            cerr << "ERROR : FitParameters format error" << endl;
            return false;
        }

        cerr << buffer << endl;

        //Reading these 7 types of information
        string parname = buf1;
        bool is_fixed = true;
        double initial_val = 0;
        double lower_bound_val = 0;
        double upper_bound_val = 0;
        int unit_type = 0;
        int color_val = 0;

        if (buf2 == "fix") {
            is_fixed = true;

            //Internal Detector Energy Response parameters
            if (parname == "Alpha_Internal") is_fixed_Alpha_Internal = true;
            if (parname == "kB_Internal") is_fixed_kB_Internal = true;
            if (parname == "R_Internal") is_fixed_R_Internal = true;

            //External Detector Energy Response parameters
            if (parname == "Alpha_External") is_fixed_Alpha_External = true;
            if (parname == "kB_External") is_fixed_kB_External = true;
            if (parname == "R_External") is_fixed_R_External = true;
        } else if (buf2 == "free") {
            is_fixed = false;
            is_fixed_all = false;
        } else {
            cerr << "ERROR : FitParameters format error (fix or free)" << endl;
            return false;
        }

        double expected_rate = 0;
        double expected_rate_e = 0;
        bool gaussian_prior = false;

        //If the initial value of fit parameter is "expected", the parameter is a GP parameter.
        //Otherwise the parameter is an UPParameter.
        if (buf3 == "expected") {
            bool is_expected = false;
            for (int j = 0; j < _NofGPParameter; j++) {
                if (parname == _GPParameter[j]) {
                    initial_val = _RateGPParameter[j];

                    expected_rate = _RateGPParameter[j];
                    expected_rate_e = _RateGPParameter_e[j];
                    gaussian_prior = _IsGPParameter[j];
                    is_expected = true;
                    break;
                }
            }


            if (is_expected == false) {
                cerr << "ERROR : FitParameters format error (initial value w/ expected)" << endl;
                return false;
            }
        } else {
            if (!(istringstream(buf3) >> initial_val)) {
                cerr << "ERROR : FitParameters format error (initial value)" << endl;
                return false;
            }
        }

        if (!(istringstream(buf4) >> lower_bound_val)) {
            cerr << "ERROR : FitParameters format error (lower bound value)" << endl;
            return false;
        }

        if (!(istringstream(buf5) >> upper_bound_val)) {
            cerr << "ERROR : FitParameters format error (upper bound value)" << endl;
            return false;
        }


        if (!(istringstream(buf7) >> color_val)) {
            cerr << "ERROR : FitParameters format error (color value)" << endl;
            return false;
        }

        buildParameter = new Parameter(parname, "fit_parameter", is_fixed);
        buildParameter->SetDValue("initial", initial_val);
        buildParameter->SetDValue("low", lower_bound_val);
        buildParameter->SetDValue("high", upper_bound_val);
        buildParameter->SetSValue("unit", buf6);


        //These are the parameters that's not event rate, a.k.a detector response etc.
        if (parname == "Mean_Monochromatic" || parname == "Sigma_Monochromatic" || parname == "Scaling" || parname == "Alpha_Internal" || parname == "kB_Internal" || parname == "R_Internal" || parname == "Alpha_External" || parname == "kB_External" || parname == "R_External") {
            _NofFitParameters++;
            gEnergyResposeParameter[parname] = buildParameter;
        } else {
            //These values will be 0(False) for non-penalized parameter.
            //For parameter with a penalty parameter, the corresponding rate can be set.
            buildParameter->SetDValue("expected", expected_rate);
            buildParameter->SetDValue("error", expected_rate_e);
            if (gaussian_prior) {
                buildParameter->SetSValue("prior", "GP");
            } else {
                buildParameter->SetSValue("prior", "UP");
            }
            _NofFitParameters++;
            gFitParameter[parname] = buildParameter;
            orderedRateName.push_back(parname);
            // _NofFitParameters_Rate++;

        }
    }
    cerr << endl;
    cerr << "=================================================================" << endl;

    if (gEnergyResposeParameter.size() + gFitParameter.size() != _NofFitParameters) {
        cerr << gFitParameter.size() << endl;
        cerr << gEnergyResposeParameter.size() << endl;
        cerr << _NofFitParameters << endl;
        cerr << "ERROR : FitParameters format error" << endl;
        return false;
    }

    //********************************************************************//
    //Internal Detector Response
    _FixParameters_All = is_fixed_all;

    if (_FixParameters_All == true) {
        cerr << "(fix all fit parameters)" << endl;
        cerr << endl;
    }

    if (is_fixed_Alpha_Internal == true) _FixEnergyScale_Internal = true;
    if (is_fixed_kB_Internal == true && is_fixed_R_Internal == true) _FixEnergyNonlinearity_Internal = true;

    if (_FixEnergyScale_Internal == true) {
        cerr << "(fix energy scale at internal)" << endl;
        cerr << endl;
    }

    if (_FixEnergyNonlinearity_Internal == true) {
        cerr << "(fix energy nonlinearity at internal)" << endl;
        cerr << endl;
    }

    //********************************************************************//
    //External Detector Response

    if (is_fixed_Alpha_External == true) _FixEnergyScale_External = true;
    if (is_fixed_kB_External == true && is_fixed_R_External == true) _FixEnergyNonlinearity_External = true;

    if (_FixEnergyScale_External == true) {
        cerr << "(fix energy scale at external)" << endl;
        cerr << endl;
    }

    if (_FixEnergyNonlinearity_External == true) {
        cerr << "(fix energy nonlinearity at external)" << endl;
        cerr << endl;
    }
    //********************************************************************//
    buildParameter = new Parameter("isFixedEnergyResponseParam", "status_parameter", true);
    buildParameter->SetBValue("all", _FixParameters_All);
    buildParameter->SetBValue("energy_scale_internal", _FixEnergyScale_Internal);
    buildParameter->SetBValue("nonlinearity_internal", _FixEnergyNonlinearity_Internal);
    buildParameter->SetBValue("energy_scale_external", _FixEnergyScale_External);
    buildParameter->SetBValue("nonlinearity_external", _FixEnergyNonlinearity_External);
    gEnergyResposeParameter["isFixedESP"] = buildParameter;


    fFitParameters.close();
    //This part of the code check the GP rate.
    //GP rate is a boolean parameter telling us if this parameter is a GP parameter.
    //If _GPRate[i] == 1, then this parameter is a GP parameter, otherwise it's a UPParameter.
    cerr << "=================================================================" << endl;
    cerr << "Check Rate GPParameter:" << endl;
    cerr << "=================================================================" << endl;
    map<string, Parameter*>::iterator it;
    for ( it = gFitParameter.begin(); it != gFitParameter.end(); it++ )
    {   
        if (it->second->GetS("prior") == "GP") {
            fprintf(stderr, "%-35s  %15f  %15f\n", it->second->GetS("name").c_str(), it->second->GetD("expected"), it->second->GetD("error"));
        }
    }
    cerr << endl;
    cerr << "=================================================================" << endl;
}

//Assign:
//"livetime"       : livetime_parameter, per-run livetime
//"livetime_eff"   : livetime_parameter, per-run livetime ratio
//"total_livetime" : data_parameter, stores "livetime" and "NofBiPo214"
bool ParameterFactory::ReadLivetimeParameters() {
    //Livetime Block
    ifstream fLiveTime(Livetime_name.c_str());

    if (!fLiveTime) {
        cerr << "ERROR : cannot open input file : LiveTime" << endl;
        return false;
    }

    map < int, double > LiveTime;
    map < int, double > LiveTimeEff;

    Parameter* livetime_param = new Parameter("livetime", "livetime_parameter", true);
    Parameter* livetime_eff_param =  new Parameter("livetime_eff", "livetime_parameter", true);

    int run;
    double runtime, vetotime, livetime, ratio, deadtime; //Reading in different component of livetime


    while (fLiveTime >> run >> runtime >> vetotime >> livetime >> ratio >> deadtime) {

        std::stringstream ss;
        ss << run;
        string run_s = ss.str();

        if (LiveTime.find(run) != LiveTime.end() || LiveTimeEff.find(run) != LiveTimeEff.end()) {
            cerr << "ERROR : invalid run number : LiveTime" << endl;
            return false;
        }

        livetime_param->SetDValue(run_s, livetime);
        //LiveTimeEff[run] = ratio;
        livetime_eff_param->SetDValue(run_s, ratio);
    }

    fLiveTime.close();

    ifstream fRunList(RunList_name.c_str());
    if (!fRunList) {
        cerr << "ERROR : cannot open input file : RunList" << endl;
        return false;
    }

    //   int run;
    double _LiveTime = 0;
    double _NumberOfBiPo214 = 0.0;

    //Summing over all livetime
    while (fRunList >> run) {

        std::stringstream ss;
        ss << run;
        string run_s = ss.str();

        if (run < _StartRun || run > _LastRun) continue;

        _LiveTime += livetime_param->GetD(run_s);

        char filename[256];

        sprintf(filename, "%s/run%06d.dat", BiPo_name.c_str(), run);

        ifstream fBiPo(filename);

        if (!fBiPo) {
            cerr << "ERROR : cannot open input file : BiPo -- > Run " << run << endl;
            return false;
        }

        double n_BiPo214;

        if (!(fBiPo >> dummy >> n_BiPo214 >> dummy)) {
            cerr << "ERROR : cannot read input file : BiPo -- > Run " << run << endl;
            return false;
        }

        fBiPo.close();

        _NumberOfBiPo214 += n_BiPo214;
    }

    fRunList.close();

    //_LiveTime = _LiveTime / 60.0 / 60.0 / 24.0 / 1000000.0; // microsec -> day
    _LiveTime = _LiveTime / 60.0 / 60.0 / 24.0; // sec -> day

    cerr << "LiveTime (days) = " << _LiveTime << endl;
    cerr << "NumberOfBiPo214 = " << _NumberOfBiPo214 << endl;
    cerr << endl;

    buildParameter = new Parameter("total_livetime", "data_parameter", true);
    buildParameter->SetDValue("livetime", _LiveTime);
    buildParameter->SetDValue("NofBiPo214", _NumberOfBiPo214);
    gModelParameter["total_livetime"] = buildParameter;
    gModelParameter["livetime"] = livetime_param;
    gModelParameter["livetime_eff"] = livetime_eff_param;
    cerr << "Finished Building Parameters!" << endl;
    //End of Livetime Block
}

//List of physics constant being dumped into physics constant parameter factory:
//FV_XeLS_All      : 16.2806
//FMass_XeLS_All   : buildParameter->GetD("FV_XeLS_All") * * 0.78013 * 1.0e-3
//Exposure_XeLS_All: 
bool ParameterFactory::ReadPhysicsParameters() {

    buildParameter = new Parameter("physics_constant", "physics_parameter", true);
    /////////////////////////////////////////////////////////////////
    /// DS-1 ///
    //  const double WeightFraction = 0.0252; // Ueshima-san's mail (12 1/6)
    //  const double WeightFraction = 0.0244; // Kengo-san's mail (12 May)
    //  const double EnrichmentFactor = 0.9093; // enriched Xe from Sasha's measurement
    //  const double MassNumber = 135.722; // enriched Xe from Sasha's measurement

    ///////////////////////////////////////////////////////////////////
    //Several different physical constants, including halflifes
    /// DS-5 ///
    //const double WeightFraction = 0.0296; // (= 0.3833 Xe-ton / 12.93 ton) from Ueshima's slide (20 Jan, 2014)
    //const double WeightFraction = 0.0293; // (= 0.3833 Xe-ton / 13.07 ton) from Ueshima's e-mail (6 Feb, 2016)
    const double WeightFraction = 0.031285;   // (= 0.745 Xe-ton / (30.53*0.78) ton) 
    const double EnrichmentFactor = 0.9086; // enriched Xe from Sasha's measurement (5 Nov, 2014)
    const double MassNumber = 135.80; // enriched Xe from Sasha's measurement (5 Nov, 2014)

    //const double NA = 6.022045e+23;
    const double NA = 6.0221409e+23;//from PDG2015
    const double mnu = 0.150; // eV (effective mass) : Klapdor et al., Phys. Lett. B 586 (2004) 198
    double HalfLife_0nu_mnu2 = 1.14e+24; // y eV2 : Nucl. Phys. A793 (2007) 213 (QRPA model)

    const double HalfLife_Xe136_2nu = 2.15e+21; // year : (2016) (zen 400 mean)
    const double HalfLife_Xe136_0nu = HalfLife_0nu_mnu2 / pow(mnu, 2);

    const double MeanLife_Xe136_0nu = HalfLife_Xe136_0nu / log(2.0);
    const double MeanLife_Xe136_2nu = HalfLife_Xe136_2nu / log(2.0);

    double R_scale=1.006;
    buildParameter->SetDValue("R_scale", 1.006);

    // Xe-LS all volume
    cerr << "set Xe-LS all volume [defined in MC]:" << endl;

    //_FiducialVolume_XeLS_All = 16.2806; // m3 (volume from MC)
    //Converting FV to Fiducial Mass
    // _FiducialMass_XeLS_All = _FiducialVolume_XeLS_All * 0.78013 * 1.0e-3; // kton
    // _Exposure_XeLS_All = _FiducialMass_XeLS_All * _LiveTime; // kton-day

    buildParameter->SetDValue("FiducialVolume_XeLS_All", 30.53);
    buildParameter->SetDValue("FiducialMass_XeLS_All", buildParameter->GetD("FiducialVolume_XeLS_All") * 0.78013 * 1.0e-3);
    buildParameter->SetDValue("Exposure_XeLS_All", buildParameter->GetD("FiducialMass_XeLS_All") * gModelParameter["total_livetime"]->GetD("livetime"));

    // //    _FiducialVolume_XeLS_1m = 4.188790205; // m3 
    // _FiducialVolume_XeLS_150cm = 14.13716694; // m3 
    // _FiducialMass_XeLS_150cm = _FiducialVolume_XeLS_150cm * 0.78013 * 1.0e-3; // kton
    // _Exposure_XeLS_150cm = _FiducialMass_XeLS_150cm * _LiveTime; // kton-day
    buildParameter->SetDValue("FiducialVolume_XeLS_150cm", 14.13716694);
    buildParameter->SetDValue("FiducialMass_XeLS_150cm", buildParameter->GetD("FiducialVolume_XeLS_150cm") * 0.78013 * 1.0e-3);
    buildParameter->SetDValue("Exposure_XeLS_150cm", buildParameter->GetD("FiducialMass_XeLS_150cm") * gModelParameter["total_livetime"]->GetD("livetime"));

    cerr << "FiducialVolume_XeLS_All (m3) [nominal] = " << buildParameter->GetD("FiducialVolume_XeLS_All") << endl;
    cerr << "FiducialMass_XeLS_All (kton) [nominal] = " << buildParameter->GetD("FiducialMass_XeLS_All") << endl;
    cerr << "Exposure_XeLS_All (kton-day) [nominal] = " << buildParameter->GetD("Exposure_XeLS_All") << endl;
    cerr << endl;

    // _Mass_Xe = (_FiducialMass_XeLS_All * 1.0e+6) * WeightFraction; // kg
    // _Mass_Xe136 = _Mass_Xe * EnrichmentFactor; // kg
    // _NTarget_Xe136 = (_Mass_Xe136 * 1.0e+3) / MassNumber * NA;
    // _DecayRate_Xe136_0nu_XeLS = _NTarget_Xe136 / MeanLife_Xe136_0nu / 365.25; // events / day
    // _DecayRate_Xe136_2nu_XeLS = _NTarget_Xe136 / MeanLife_Xe136_2nu / 365.25; // events / day
    // _Rate_Xe136_0nu_XeLS_Expected = _DecayRate_Xe136_0nu_XeLS / _FiducialMass_XeLS_All; // events / day / kton
    // _Rate_Xe136_2nu_XeLS_Expected = _DecayRate_Xe136_2nu_XeLS / _FiducialMass_XeLS_All; // events / day / kton

    buildParameter->SetDValue("Mass_Xe", buildParameter->GetD("FiducialMass_XeLS_All") * 1.0e+6 * WeightFraction);
    buildParameter->SetDValue("Mass_Xe136", buildParameter->GetD("Mass_Xe") * EnrichmentFactor);
    buildParameter->SetDValue("NTarget_Xe136", (buildParameter->GetD("Mass_Xe136")* 1.0e+3) / MassNumber * NA);
    buildParameter->SetDValue("DecayRate_Xe136_0nu_XeLS", buildParameter->GetD("NTarget_Xe136") / MeanLife_Xe136_0nu / 365.25 );
    buildParameter->SetDValue("DecayRate_Xe136_2nu_XeLS", buildParameter->GetD("NTarget_Xe136") / MeanLife_Xe136_2nu / 365.25 );
    buildParameter->SetDValue("Rate_Xe136_0nu_XeLS_Expected", buildParameter->GetD("DecayRate_Xe136_0nu_XeLS") / buildParameter->GetD("FiducialMass_XeLS_All"));
    buildParameter->SetDValue("Rate_Xe136_2nu_XeLS_Expected", buildParameter->GetD("DecayRate_Xe136_2nu_XeLS") / buildParameter->GetD("FiducialMass_XeLS_All"));


    cerr << "Mass_Xe (kg) [nominal]    = " << buildParameter->GetD("Mass_Xe") << endl;
    cerr << "Mass_Xe136 (kg) [nominal] = " << buildParameter->GetD("Mass_Xe136") << endl;
    cerr << "NTarget_Xe136 [nominal]   = " << buildParameter->GetD("NTarget_Xe136") << endl;
    cerr << "DecayRate_Xe136_0nu_XeLS (events/day) = " << buildParameter->GetD("DecayRate_Xe136_0nu_XeLS") << endl;
    cerr << "DecayRate_Xe136_2nu_XeLS (events/day) = " << buildParameter->GetD("DecayRate_Xe136_2nu_XeLS") << endl;
    cerr << "Rate_Xe136_0nu_XeLS_Expected (events/day/kton) = " << buildParameter->GetD("Rate_Xe136_0nu_XeLS_Expected") << endl;
    cerr << "Rate_Xe136_2nu_XeLS_Expected (events/day/kton) = " << buildParameter->GetD("Rate_Xe136_2nu_XeLS_Expected") << endl;
    cerr << endl;

    double Concentration_U238_S2_XeLS = gModelParameter["total_livetime"]->GetD("NofBiPo214") / gModelParameter["total_livetime"]->GetD("livetime") / (14.137*R_scale*R_scale*R_scale); // events / day / m3 --- from number of BiPo214 event in R 1.2 m
    //double Concentration_U238_S2_XeLS = 50000.0 / 600000.0 / 7.238; // events / day / m3 --- from number of BiPo214 event in R 1.2 m
    double Concentration_Th232_S2_XeLS =  0.1395; // events / day / m3 --- from Tl208 event in all volume (R distribution fit) --- to be updated?? 

    // _DecayRate_U238_S2_XeLS = Concentration_U238_S2_XeLS * _FiducialVolume_XeLS_All; // events / day
    // _DecayRate_Th232_S2_XeLS = Concentration_Th232_S2_XeLS * _FiducialVolume_XeLS_All; // events / day

    buildParameter->SetDValue("DecayRate_U238_S2_XeLS", Concentration_U238_S2_XeLS * buildParameter->GetD("FiducialVolume_XeLS_All"));
    buildParameter->SetDValue("DecayRate_Th232_S2_XeLS", Concentration_Th232_S2_XeLS * buildParameter->GetD("FiducialVolume_XeLS_All"));

    cerr << "DecayRate_U238_S2_XeLS (events/day)  = " << buildParameter->GetD("DecayRate_U238_S2_XeLS") << endl;
    cerr << "DecayRate_Th232_S2_XeLS (events/day) = " << buildParameter->GetD("DecayRate_Th232_S2_XeLS") << endl;

    // _Rate_U238_S2_XeLS_Expected = _DecayRate_U238_S2_XeLS / _FiducialMass_XeLS_All;
    // _Rate_Th232_S2_XeLS_Expected = _DecayRate_Th232_S2_XeLS / _FiducialMass_XeLS_All;

    buildParameter->SetDValue("Rate_U238_S2_XeLS_Expected", buildParameter->GetD("DecayRate_U238_S2_XeLS") / buildParameter->GetD("FiducialMass_XeLS_All"));
    buildParameter->SetDValue("Rate_Th232_S2_XeLS_Expected", buildParameter->GetD("DecayRate_Th232_S2_XeLS") / buildParameter->GetD("FiducialMass_XeLS_All"));

    // _Rate_U238_S2_XeLS_Expected_e  = _Rate_U238_S2_XeLS_Expected * 0.05; // 5% error (used until 2016.04.27)
    // _Rate_U238_S2_XeLS_Expected_e = _Rate_U238_S2_XeLS_Expected * 0.17; // 17% error (= error/tag ineff. = 0.00008 / 0.00046)
    // _Rate_Th232_S2_XeLS_Expected_e = _Rate_Th232_S2_XeLS_Expected * 0.17; // 17% error (10.3 +- 1.7 events / day) --- to be updated

    buildParameter->SetDValue("Rate_U238_S2_XeLS_Expected_e", buildParameter->GetD("Rate_U238_S2_XeLS_Expected")  * 0.17);//0.17 originally
    buildParameter->SetDValue("Rate_Th232_S2_XeLS_Expected_e", buildParameter->GetD("Rate_Th232_S2_XeLS_Expected") * 0.23);

    cerr << "Rate_U238_S2_XeLS_Expected (events/day/kton)   = " << buildParameter->GetD("Rate_U238_S2_XeLS_Expected") << " +- " << buildParameter->GetD("Rate_U238_S2_XeLS_Expected_e") << endl;
    cerr << "Rate_Th232_S2_XeLS_Expected (events/day/kton)  = " << buildParameter->GetD("Rate_Th232_S2_XeLS_Expected") << " +- " << buildParameter->GetD("Rate_Th232_S2_XeLS_Expected_e") << endl;

    // Kam-LS all volume
    cerr << "set Kam-LS all volume [defined in MC]:" << endl;

    // _FiducialVolume_KamLS_All = 1134.35; // m3 (volume from MC)
    // _FiducialMass_KamLS_All = _FiducialVolume_KamLS_All * 0.78013 * 1.0e-3; // kton
    // _Exposure_KamLS_All = _FiducialMass_KamLS_All * _LiveTime; // kton-day

    buildParameter->SetDValue("FiducialVolume_KamLS_All", 1135.35);
    buildParameter->SetDValue("FiducialMass_KamLS_All", buildParameter->GetD("FiducialVolume_KamLS_All") * 0.78013 * 1.0e-3);
    buildParameter->SetDValue("Exposure_KamLS_All", buildParameter->GetD("FiducialMass_KamLS_All") * gModelParameter["total_livetime"]->GetD("livetime"));

        //    _FiducialVolume_KamLS_norm = 593.7610115; // m3 ( 3 m - 6 m ,theta < 1/2 )
    // _FiducialVolume_KamLS_norm = 437.8594761; // m3 ( 3 m - 5.5 m ,theta < 1/2 )
    // _FiducialMass_KamLS_norm = _FiducialVolume_KamLS_norm * 0.78013 * 1.0e-3; // kton
    // _Exposure_KamLS_norm = _FiducialMass_KamLS_norm * _LiveTime; // kton-day

    buildParameter->SetDValue("FiducialVolume_KamLS_norm", 437.8594761);
    buildParameter->SetDValue("FiducialMass_KamLS_norm", buildParameter->GetD("FiducialVolume_KamLS_norm") * 0.78013 * 1.0e-3);
    buildParameter->SetDValue("Exposure_KamLS_norm", buildParameter->GetD("FiducialMass_KamLS_norm") * gModelParameter["total_livetime"]->GetD("livetime"));

    cerr << "FiducialVolume_KamLS_All (m3) [nominal] = " << buildParameter->GetD("FiducialVolume_KamLS_All") << endl;
    cerr << "FiducialMass_KamLS_All (kton) [nominal] = " << buildParameter->GetD("FiducialMass_KamLS_All") << endl;
    cerr << "Exposure_KamLS_All (kton-day) [nominal] = " << buildParameter->GetD("Exposure_KamLS_All") << endl;
    cerr << endl;
    gModelParameter["physics_constant"] = buildParameter;
    return true;
     //End of Physical Constant block
    //////////////////////////////////////////////////////////////////////////////////
}







