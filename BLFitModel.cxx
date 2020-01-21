// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************
#include "BLFitModel.h"
#include <BAT/BCLog.h>
#include "TSystem.h"
#include <assert.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1D.h>
#include <TCanvas.h>  
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>

// #include <BAT/BCMath.h>

using namespace std;
using namespace TMath;

void BLFitModel::StackingDataPoints(double ** _E_Mean, double ** _Data, double ** _Data_e) {
    std::vector<double > pars = GetInitialPosition();
    Parameter* currentEnergyResponse = ERes->setupEnergyResponse(fParameterInfoMap["energy_response_info"], pars);
    double scaling = pars[fParameterInfoMap["energy_response_info"]->GetI("Scaling")];
    double Rate=1.2;
    gRandom = new TRandom3(1000);
    for(int n=0; n < Data->GetFVParameter()->GetI("total_bin"); n++) {
        int n1 = n / Data->GetFVParameter()->GetI("theta_bin");
        int n2 = n % Data->GetFVParameter()->GetI("theta_bin");
        _E_Mean[n] = new double[_NofEnergyBin];//Average Energy
        _Data[n] = new double[_NofEnergyBin];//Data Counts
        _Data_e[n] = new double[_NofEnergyBin];//Error of Data
        for(int i=0; i<_NofEnergyBin; i++) {
            double expected_value = 0.0;
            double E_Mean  = _EnergyMin + _EnergyBinWidth * double(i+0.5);//Bin Center
            double E_Lower = _EnergyMin + _EnergyBinWidth * double(i+0.0);//Bin Lower Edge
            double E_Upper = _EnergyMin + _EnergyBinWidth * double(i+1.0);//Bin Upper Edge
            for (int j=0; j < fitParameterVec.size(); j++) {
                if (fParameterInfoMap.find(fitParameterVec[j]) != fParameterInfoMap.end()) {
                    expected_value += ERes->GetStackedHist(fParameterInfoMap[fitParameterVec[j]], currentEnergyResponse, fitParameterVec[j], E_Mean, n, pars);
                    //expected_value += Data->GetHistogramInterpolation(fParameterInfoMap[fitParameterVec[j]], E_Mean, n);
                }
            }
            _E_Mean[n][i] = E_Mean;
            _Data[n][i] = gRandom->Poisson(Rate*expected_value);
            _Data_e[n][i] = sqrt(_Data[n][i]);
        }
    }
}

//This group up isotopes
std::map<std::string, TH3D*> BLFitModel::GetMCHist(std::vector<double>& pars) {
    Parameter* currentEnergyResponse = ERes->setupEnergyResponse(fParameterInfoMap["energy_response_info"], pars);
    std::map<std::string, TH3D*> MCMap;
    // Int_t nbinx = 0;
    // Double_t xbins[1000];
    // for(int n=0; n<_NofEnergyBin; n++) {
    //     if(n+1==1000) {
    //         cerr << "ERROR : xbins is out of range" << endl;
    //         abort();
    //     }
    //     xbins[nbinx]   = _EnergyMin + _EnergyBinWidth * double(n+0.0);
    //     xbins[nbinx+1] = _EnergyMin + _EnergyBinWidth * double(n+1.0);
    //     nbinx++;
    // }
    for (int j=0; j < fitParameterVec.size(); j++) {
        if (fParameterInfoMap.find(fitParameterVec[j]) != fParameterInfoMap.end()) {
            string histName = fitParameterVec[j];
            string mapName;
            if (histName=="Rate_Xe136_0nu_XeLS") {
                mapName = "^{136}Xe 0#nu#beta#beta";
                pars[fParameterInfoMap["Rate_Xe136_0nu_XeLS"]->GetI("Rate_Xe136_0nu_XeLS")] = Get0vbbCL();
            } else if (histName=="Rate_Xe136_2nu_XeLS") {
                mapName = "^{136}Xe 2#nu#beta#beta";
            }else if (histName=="Rate_Ag110_XeLS") {
                mapName = "^{110m}Ag";
            } else if ((histName.find("XeLS")!=string::npos) && (
                        (histName.find("U238")!=string::npos) ||
                        (histName.find("Th232")!=string::npos) ||
                        (histName=="Rate_Bi210_XeLS") ||
                        (histName=="Rate_Po210_XeLS") ||
                        (histName=="Rate_Kr85_XeLS") ||
                        (histName=="Rate_K40_XeLS")))
            {
                mapName = "^{238}U+^{232}Th+^{210}Bi+^{210}Po+^{85}Kr+^{40}K";
            } else if ((histName.find("film")!=string::npos) || (histName.find("KamLS")!=string::npos)) {
                mapName = "IB/External";
            } else if ((histName.find("Spallation")!=string::npos) ||
                        (histName=="Rate_C11_XeLS") ||
                        (histName=="Rate_C10_XeLS"))

             {
                mapName = "Spallation";
            } else {
                mapName = "Other";
            }
            if (MCMap.find(mapName) == MCMap.end()) {
                MCMap[mapName] = Data->GetEmptyHistogram(mapName);
            }
            for(int i=0; i<_NofEnergyBin; i++) {
                double E_Mean  = _EnergyMin + _EnergyBinWidth * double(i+0.5);//Bin Center
                double E_Lower = _EnergyMin + _EnergyBinWidth * double(i+0.0);//Bin Lower Edge
                double E_Upper = _EnergyMin + _EnergyBinWidth * double(i+1.0);//Bin Upper Edge
                for(int n=0; n < Data->GetFVParameter()->GetI("total_bin"); n++) {
                    int n1 = n / Data->GetFVParameter()->GetI("theta_bin");
                    int n2 = n % Data->GetFVParameter()->GetI("theta_bin");
                    MCMap[mapName]->SetBinContent(i+1, n1+1, n2+1, MCMap[mapName]->GetBinContent(i+1, n1+1, n2+1) + ERes->GetStackedHist(fParameterInfoMap[histName], currentEnergyResponse, histName, E_Mean, n, pars));
                }
            }
        }
    }
    delete currentEnergyResponse;
    return MCMap;
}

//This displays every single isotope instead of grouping them up.
std::map<std::string, TH3D*> BLFitModel::GetMCHist_Isotope(std::vector<double>& pars, std::vector<string>& invesitgate) {


    Parameter* currentEnergyResponse = ERes->setupEnergyResponse(fParameterInfoMap["energy_response_info"], pars);
    std::map<std::string, TH3D*> MCMap;

    for (int j=0; j < fitParameterVec.size(); j++) {
        if (fParameterInfoMap.find(fitParameterVec[j]) != fParameterInfoMap.end()) {
            string histName = fitParameterVec[j];
            currentFitParameter = Param->GetFitParameter(histName);
            // if ((currentFitParameter->isFixed()) && currentFitParameter->GetD("initial") == 0) continue;

            string mapName = "Other";
            if (invesitgate.size() == 0) mapName=histName;
            for (int ivindex=0; ivindex<invesitgate.size(); ivindex++) {
                if (invesitgate[ivindex] == histName) mapName=histName;

            }

            //Initialize mapname TH1D
            if (MCMap.find(mapName) == MCMap.end()) {
                MCMap[mapName] = Data->GetEmptyHistogram(mapName);
            }

            for(int i=0; i<_NofEnergyBin; i++) {
                double E_Mean  = _EnergyMin + _EnergyBinWidth * double(i+0.5);//Bin Center
                double E_Lower = _EnergyMin + _EnergyBinWidth * double(i+0.0);//Bin Lower Edge
                double E_Upper = _EnergyMin + _EnergyBinWidth * double(i+1.0);//Bin Upper Edge
                for(int n=0; n < Data->GetFVParameter()->GetI("total_bin"); n++) {
                    int n1 = n / Data->GetFVParameter()->GetI("theta_bin");
                    int n2 = n % Data->GetFVParameter()->GetI("theta_bin");
                    MCMap[mapName]->SetBinContent(i+1, n1+1, n2+1, MCMap[mapName]->GetBinContent(i+1, n1+1, n2+1) + ERes->GetStackedHist(fParameterInfoMap[histName], currentEnergyResponse, histName, E_Mean, n, pars));
                }
            }
        }
    }
    delete currentEnergyResponse;
    return MCMap;
}

TH3D* BLFitModel::GetDataHist() {
    Int_t nbinx = 0;
    Double_t xbins[1000];
    for(int n=0; n<_NofEnergyBin; n++) {
        if(n+1==1000) {
            cerr << "ERROR : xbins is out of range" << endl;
            abort();
        }
        xbins[nbinx]   = _EnergyMin + _EnergyBinWidth * double(n+0.0);
        xbins[nbinx+1] = _EnergyMin + _EnergyBinWidth * double(n+1.0);
        nbinx++;
    }
    TH3D* dataHist = Data->GetEmptyHistogram("h_data");
    for(int i=0; i<_NofEnergyBin; i++) {
        for(int n=0; n < Data->GetFVParameter()->GetI("total_bin"); n++) {
            int n1 = n / Data->GetFVParameter()->GetI("theta_bin");
            int n2 = n % Data->GetFVParameter()->GetI("theta_bin");

            double binvalue = 0.0;
            if (fakeData) {
                binvalue += _Data_fake[n][i];
            } else{
                binvalue += _Data[n][i];
            }
            dataHist->SetBinContent(i+1, n1+1, n2+1, binvalue);
        }
    }
    return dataHist;
}

double BLFitModel::GetBackgroundCountsROI(const std::string& isotope_name, std::vector<double> pars) {
    
    Parameter* currentEnergyResponse = ERes->setupEnergyResponse(fParameterInfoMap["energy_response_info"], pars);
    double elow=2.35;
    double ehi=2.70;
    double r_threshold=160.0;

    double outputs = 0.0;
    if (fParameterInfoMap.find(isotope_name) != fParameterInfoMap.end()) {
        for(double E=0.0;E<20.0;E+=0.001){
            for(int n=0; n < Data->GetFVParameter()->GetI("total_bin"); n++) {

                int n1 = n / Data->GetFVParameter()->GetI("theta_bin");
                int n2 = n % Data->GetFVParameter()->GetI("theta_bin");

                double volume_multiplier = 0.0;
                if ((E>=elow) && (E<ehi) &&n<10){
                    volume_multiplier = 0.001/0.05;
                    // if (Data->GetRBinMax(n1)<=r_threshold) {
                    //     volume_multiplier = 0.001/0.05;
                    // } 
                    // // else if ((Data->GetRBinMin(n1)<r_threshold) and (Data->GetRBinMax(n1)>=r_threshold)) {
                    // //     double volume = Data->getFiducialVolume(Data->GetRBinMin(n1), r_threshold, Data->GetThetaBinMin(n2), Data->GetThetaBinMax(n2));
                    // //     double binvolume = Data->getFiducialVolume(Data->GetRBinMin(n1), Data->GetRBinMax(n1), Data->GetThetaBinMin(n2), Data->GetThetaBinMax(n2));
                    // //     if ((volume != 0.0) && (binvolume != 0.0)) {
                    // //         volume_multiplier = volume/binvolume;
                    // //     }
                    // // }
                }

                outputs += ERes->GetStackedHist(fParameterInfoMap[isotope_name], currentEnergyResponse, isotope_name, E, n, pars)* volume_multiplier;
            }
        }
    }
    return outputs;
}

double BLFitModel::GetDataCountsROI() {
    double elow=2.35;
    double ehi=2.70;

    _E_Mean = Data->GetData("Emean");
    _Data = Data->GetData("data");

    double outputs = 0.0;
    for(int n=0; n < Data->GetFVParameter()->GetI("total_bin"); n++) {
        for(int i=0; i<_NofEnergyBin; i++) {

            double volume_multiplier = 0.0;
            if ((_E_Mean[n][i]>=elow) && (_E_Mean[n][i]<ehi) && n<10){
                    volume_multiplier = 1.0;
            }

            cout<<_Data[n][i]<<"   "<<_E_Mean[n][i]<<"   "<<volume_multiplier<<endl;
            outputs+=_Data[n][i]*volume_multiplier;
        }
    }
    return outputs;
}

void BLFitModel::BuildParameterGraph() {
    cerr<<"Start building graph!"<<endl;
    double x[10000], y[10000];
    for (int j=0; j < fitParameterVec.size(); j++) {
        string histName = fitParameterVec[j];
        if (fParameterInfoMap.find(histName) != fParameterInfoMap.end()) {
            for(int n=0; n < Data->GetFVParameter()->GetI("total_bin"); n++) {
                int N = 0;
                for(int i=0; i<_nbin; i++) {
                    x[i] = _min + _bin_width * double(i+0.5);
                    y[i] = Data->GetHistogramInterpolation(fParameterInfoMap[histName], i, n);
                    N++;
                }
                fParameterInfoMap[histName]->SetGraphAddress(N, x, y, n);
            }

        }
    }
    cerr<<"Successfully Built Graph!"<<endl;

}

// ---------------------------------------------------------
BLFitModel::BLFitModel(const std::string& name)
    : BCModel(name)
{

    Param = new ParameterFactory();
    Data = new DataFactory(Param);
    ERes = new EnergyResponse(Param, Data);
    cerr<< "Run Succeed!"<<endl;
    // Comment this if running with fake data
    if (!fakeData) {
        _E_Mean = Data->GetData("Emean");
        _Data = Data->GetData("data");
        _Data_e = Data->GetData("error");
    }


    // if( getenv("DOUBLEBETA_ANALYSIS_ETH")==NULL
    // ||
    // getenv("DOUBLEBETA_ANALYSIS_UPP_ETH")==NULL
    // ){
    // cerr << "ERROR: Cannot set analysis parameters" << endl;
    // cerr << "Please check environmental variable" << endl;
    // cerr << "DOUBLEBETA_ANALYSIS_ETH: " << getenv("DOUBLEBETA_ANALYSIS_ETH") << endl;
    // cerr << "DOUBLEBETA_ANALYSIS_UPP_ETH: " << getenv("DOUBLEBETA_ANALYSIS_UPP_ETH") << endl;
    // cerr << endl;

    // abort();
    // }

    // _EnergyThreshold = atof(getenv("DOUBLEBETA_ANALYSIS_ETH"));
    // _UpperEnergyThreshold = atof(getenv("DOUBLEBETA_ANALYSIS_UPP_ETH"));

    _EnergyThreshold=0.5;
    _UpperEnergyThreshold=4.8;


    param_id = 0;
    fitParameterVec = Param->GetOrderedFitList();
    // for (int j=0; j < fitParameterVec.size(); j++) { cout<<fitParameterVec[j]<<endl;}
    // abort();
    for (int j=0; j < fitParameterVec.size(); j++) {
        currentFitParameter = Param->GetFitParameter(fitParameterVec[j]);
        //n==0 to make sure only do this once
        if (fitParameterVec[j] == "Rate_Monochromatic"){
            //Adding Rate_Monochromatic
            buildParameter = new Parameter(fitParameterVec[j], "info_parameter", true);
            AddParameter(fitParameterVec[j], currentFitParameter->GetD("low"), currentFitParameter->GetD("high"), fitParameterVec[j], currentFitParameter->GetS("unit"));
            GetParameters().Back().SetPriorConstant();
            fittingVector.push_back(fitParameterVec[j]);
            buildParameter->SetIValue(fitParameterVec[j], param_id);
            buildParameter->SetSValue("unit", currentFitParameter->GetS("unit"));
            initialPosition.push_back(currentFitParameter->GetD("initial"));
            if (currentFitParameter->GetB("is_fixed")) {
                GetParameters().Back().Fix(currentFitParameter->GetD("initial"));
            }
            param_id++;

            //Adding Mean_Monochromatic
            currentFitParameter = Param->GetEnergyResponse("Mean_Monochromatic");
            AddParameter("Mean_Monochromatic", currentFitParameter->GetD("low"), currentFitParameter->GetD("high"), "Mean_Monochromatic", currentFitParameter->GetS("unit"));
            GetParameters().Back().SetPriorConstant();
            fittingVector.push_back("Mean_Monochromatic");
            buildParameter->SetIValue("Mean_Monochromatic", param_id);
            initialPosition.push_back(currentFitParameter->GetD("initial"));
            if (currentFitParameter->GetB("is_fixed")) {
                if (fitParameterVec[j] == "Rate_Xe136_0nu_XeLS") {
                    GetParameters().Back().Fix((atof(getenv("ZERONU_RATE"))));
                } else {
                    GetParameters().Back().Fix(currentFitParameter->GetD("initial"));
                }
            }
            //cout<<GetParameters().Back().<<currentFitParameter->GetD("initial")<<"  "<<param_id<<"    "<<currentFitParameter->GetB("is_fixed")<<endl;
            param_id++;

            //Adding Sigma_Monochromatic
            currentFitParameter = Param->GetEnergyResponse("Sigma_Monochromatic");
            AddParameter("Sigma_Monochromatic", currentFitParameter->GetD("low"), currentFitParameter->GetD("high"), "Sigma_Monochromatic", currentFitParameter->GetS("unit"));
            GetParameters().Back().SetPriorConstant();
            fittingVector.push_back("Sigma_Monochromatic");
            buildParameter->SetIValue("Sigma_Monochromatic", param_id);
            initialPosition.push_back(currentFitParameter->GetD("initial"));
            if (currentFitParameter->GetB("is_fixed")) {
                GetParameters().Back().Fix(currentFitParameter->GetD("initial"));
            }
            param_id++;

            fParameterInfoMap[fitParameterVec[j]] = buildParameter;
            continue;
        }
        string param_name = fitParameterVec[j];
        if (fParameterInfoMap.find(fitParameterVec[j]) == fParameterInfoMap.end()) {
            buildParameter = new Parameter(fitParameterVec[j], "info_parameter", true);
            buildParameter->SetSValue("unit", currentFitParameter->GetS("unit"));
            buildParameter->SetIValue("i_ES", j);
            if (!buildParameter->readyToFit()) {
                delete buildParameter;
                continue;
            }
            fParameterInfoMap[fitParameterVec[j]] = buildParameter;
        }
        AddParameter(param_name, currentFitParameter->GetD("low"), currentFitParameter->GetD("high"), param_name, currentFitParameter->GetS("unit"));
        fittingVector.push_back(param_name);
        initialPosition.push_back(currentFitParameter->GetD("initial"));
        if (currentFitParameter->GetS("prior") == "GP") {
            GetParameters().Back().SetPrior(new BCGaussianPrior(currentFitParameter->GetD("expected"), currentFitParameter->GetD("error")));
        } else {
            GetParameters().Back().SetPriorConstant();
        }
        if (currentFitParameter->GetB("is_fixed")) {
            GetParameters().Back().Fix(currentFitParameter->GetD("initial"));
        }
        buildParameter = fParameterInfoMap[fitParameterVec[j]];
        //buildParameter->SetSValue(param_name, param_name);
        buildParameter->SetIValue(param_name, param_id);
        param_id++;
    }
    buildParameter = new Parameter("energy_response_info", "info_parameter", true);
    for( unsigned int i = 0; i < 6; i++ ) {
        currentFitParameter = Param->GetEnergyResponse(EResParameterName[i]);
        if (EResParameterName[i] == "Alpha_Internal") {
            // for(double Alpha=0.975;Alpha<=1.025+0.0001;Alpha+=0.005){ // default
            // for(double Alpha=0.986; Alpha<=1.006+0.0001; Alpha+=0.001) { // 21 loop for 0.8-4.8MeV fit #use this
            // for(double Alpha=0.975;Alpha<=1.000+0.0001;Alpha+=0.005){ // 6 loop for 0.5-4.8MeV fit
            // for(double Alpha=0.980;Alpha<=0.995+0.0001;Alpha+=0.003){ // 6 loop for 0.5-4.8MeV fit
            AddParameter(EResParameterName[i], 0.975, (1.025+0.0001), EResParameterName[i], currentFitParameter->GetS("unit"));
            fittingVector.push_back(EResParameterName[i]);
            GetParameters().Back().SetPrior(new BCGaussianPrior(1.0, Param->GetModelParameter("DeltaEnergy")->GetD("initial")));
        } else {
            AddParameter(EResParameterName[i], currentFitParameter->GetD("low"), currentFitParameter->GetD("high"), EResParameterName[i], currentFitParameter->GetS("unit"));
            fittingVector.push_back(EResParameterName[i]);
            if (EResParameterName[i] == "Alpha_External") {
                GetParameters().Back().SetPrior(new BCGaussianPrior(1.0, Param->GetModelParameter("DeltaEnergy")->GetD("initial")));
            } else {
                GetParameters().Back().SetPriorConstant();
            }
        }
        if (currentFitParameter->isFixed()) {
            GetParameters().Back().Fix(currentFitParameter->GetD("initial"));
        }
        buildParameter->SetIValue(EResParameterName[i], param_id);
        initialPosition.push_back(currentFitParameter->GetD("initial"));
        param_id++;
    }

    double deltaN = Param->GetModelParameter("DeltaNormalize")->GetD("initial");
    AddParameter("Scaling", 0.9, 1.1, "Scaling", "Unitless");
    fittingVector.push_back("Scaling");
    GetParameters().Back().SetPrior(new BCGaussianPrior(1.0, deltaN));
    buildParameter->SetIValue("Scaling", param_id);
    initialPosition.push_back(1.0);
    GetParameters().Back().Fix(1.0);
    // if (currentFitParameter->GetB("is_fixed")) {
    //     GetParameters().Back().Fix(1.0);
    // }
    param_id++;
    fParameterInfoMap["energy_response_info"] = buildParameter;
    BuildParameterGraph();
    if (fittingVector.size() != GetNParameters()) {
        cerr<<"Fitting vector size and parameter size do not match!"<<endl;
        abort();
    }
    for (int j=0; j < fittingVector.size(); j++) {
        if (fittingVector[j].find("Rate") != string::npos) {
            fittingVectorRate.push_back(fittingVector[j]);
        }
    }
    if (fakeData) {
        // Using fake data instead of real data
        StackingDataPoints(_E_Mean_fake, _Data_fake, _Data_e_fake);
        for (int i=0; i<initialPosition.size(); i++) { cout<<i<<"="<<initialPosition[i]<<"|";}
    }
}

// ---------------------------------------------------------
BLFitModel::~BLFitModel()
{
    // destructor
}


double BLFitModel::LogLikelihood(const std::vector<double>& pars)
{   
    //cout<<pars[fParameterInfoMap["energy_response_info"]->GetI("kB_Internal")]<<"|"<<pars[fParameterInfoMap["energy_response_info"]->GetI("R_Internal")]<<endl; ;
    cout<<"Chain "<<GetCurrentChain()<<": ";
    double logprob = 0.0;
    currentEnergyResponse[GetCurrentChain()] = ERes->setupEnergyResponse(fParameterInfoMap["energy_response_info"], pars);
    logprob -= 0.5*currentEnergyResponse[GetCurrentChain()]->GetD("ChiSquare_EnergyNonlinearity_External");
    logprob -= 0.5*currentEnergyResponse[GetCurrentChain()]->GetD("ChiSquare_EnergyNonlinearity_Internal");
    double scaling = pars[fParameterInfoMap["energy_response_info"]->GetI("Scaling")];

    for(int n=0; n < Data->GetFVParameter()->GetI("total_bin"); n++) {
        for(int energy=0; energy<_NofEnergyBin; energy++) {
            double observed_events;
            double accurateEnergy;
            if (fakeData) {
                observed_events = _Data_fake[n][energy];
                accurateEnergy = _E_Mean_fake[n][energy];
            } else {
                observed_events = _Data[n][energy];
                accurateEnergy = _E_Mean[n][energy];
            }
            if (accurateEnergy < _EnergyThreshold || accurateEnergy >= _UpperEnergyThreshold) continue;
            double expected_output = 0.0;
            for (int j=0; j < fittingVectorRate.size(); j++) {
                currentFitParameter = Param->GetFitParameter(fittingVectorRate[j]);
                if ((currentFitParameter->isFixed()) && currentFitParameter->GetD("initial") == 0) continue;
                //expected_output += Data->GetHistogramInterpolation(fParameterInfoMap[fitParameterVec[j]], accurateEnergy, n) * pars[fParameterInfoMap[fitParameterVec[j]]->GetI(fitParameterVec[j])];
                expected_output += ERes->GetStackedHist(fParameterInfoMap[fittingVectorRate[j]], currentEnergyResponse[GetCurrentChain()], fittingVectorRate[j], accurateEnergy, n, pars) * scaling;
                // if (fitParameterVec[j] == "Rate_Th232_S2_film_0") {
                //     cout<<accurateEnergy<<"   "<<ERes->GetStackedHist(fParameterInfoMap[fitParameterVec[j]], currentEnergyResponse[GetCurrentChain()], fitParameterVec[j], accurateEnergy, n, pars)<<endl;
                // }
                // if ((accurateEnergy > 2.2) && (accurateEnergy <2.7)) {
                // cout<<fitParameterVec[j]<<":"<<observed_events <<"    "<<expected_output<<"|";
                // }
            }
            if ((expected_output>0.0) && (observed_events >= 0.0)) {
                double factorial = 0.0;
                if (observed_events==0.0) {
                    factorial = 0.0;
                } else {
                    //factorial = 0.5 * log( 2 * TMath::Pi() * observed_events) + observed_events * log( observed_events / TMath::E());
                    factorial = observed_events * log(observed_events) - observed_events;
                    //factorial = observed_events * log(observed_events) - observed_events; //Stirling
                    //factorial = log(TMath::Factorial(observed_events));
                }
                logprob += observed_events * log(expected_output) - expected_output - factorial;
            } else if (observed_events == 0.0) {
                logprob -= (expected_output - observed_events);
            }
        }
    }
    cout<<" Logprob:"<<logprob<<endl;
    return logprob;
}


