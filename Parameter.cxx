#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <Parameter.h>
#include "TSystem.h"

using namespace std;

Parameter::Parameter(const std::string& n, const std::string& t, bool fixed) {
    stringParamMap["name"] = n;
    stringParamMap["type"] = t;
    boolParamMap["is_fixed"] = fixed;
    if(getenv("DOUBLEBETA_ANALYSIS_PARAMETER_DIR")==NULL
    ){
    cerr << "DOUBLEBETA_ANALYSIS_PARAMETER_DIR: " << getenv("DOUBLEBETA_ANALYSIS_PARAMETER_DIR") << endl;
    cerr << endl;

    abort();
    }
    RadiusBin_name = ((string) getenv("DOUBLEBETA_ANALYSIS_PARAMETER_DIR")) + "/RadiusBin.dat";
    ThetaBin_name = ((string) getenv("DOUBLEBETA_ANALYSIS_PARAMETER_DIR")) + "/ThetaBin.dat";
    if (t == "info_parameter") {
        if (BuildFiducialVolume()) {
    	   fitParameterBuilt = EditParameter(n);
        } else {
            cerr << "Error building fiducial volume for information parameter!" << endl;
        }
    }
}




bool Parameter::EditParameter(string name) {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Start of Fit Parameter Assignment
    if(name=="Rate_Xe136_0nu_XeLS") {
        basename.push_back("Xe136_0nu_XeLS");
    }else if(name=="Rate_Xe136_2nu_XeLS") {
        basename.push_back("Xe136L_2nu_XeLS");
        basename.push_back("Xe136H_2nu_XeLS");
    }else if(name=="Rate_U238_S1_XeLS") {
        basename.push_back("Pa234m_XeLS");
        // if(_FitParameters[i]=="Rate_U238_S1_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Pa234m_XeLS * h_Pa234m_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Pa234m_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name=="Rate_U238_S2_XeLS") {
        basename.push_back("Pb214m_XeLS");
        basename.push_back("Bi214m_XeLS");
        tagging.push_back(false);
        tagging.push_back(true);
        // if(_FitParameters[i]=="Rate_U238_S2_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Pb214m_XeLS * h_Pb214m_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Pb214m_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width) + (TaggingEfficiency_Bi214m_XeLS * Ratio_Bi214m_XeLS * h_Bi214m_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Bi214m_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name=="Rate_Th232_S1_XeLS") {
        basename.push_back("Ac228m_XeLS");
        //    if(_FitParameters[i]=="Rate_Th232_S1_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Ac228m_XeLS * h_Ac228m_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Ac228m_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        //    }
    }else if(name=="Rate_Th232_S2_XeLS") {
        basename.push_back("Bi212m_XeLS");
        basename.push_back("Pileup_XeLS");
        basename.push_back("Tl208m_XeLS");
        tagging.push_back(true);
        tagging.push_back(true);//Ratio:Bi214m
        tagging.push_back(false);
        // if(_FitParameters[i]=="Rate_Th232_S2_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N]  =(TaggingEfficiency_Bi212m_XeLS * Ratio_Bi212m_XeLS * h_Bi212m_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Bi212m_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width)
        //          + (TaggingEfficiency_Pileup_XeLS * Ratio_Bi212m_XeLS * h_Pileup_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Pileup_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width)
        //          + (Ratio_Tl208m_XeLS * h_Tl208m_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Tl208m_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name=="Rate_K40_XeLS") {
        basename.push_back("K40p_XeLS");
        basename.push_back("K40m_XeLS");
        // if(_FitParameters[i]=="Rate_K40_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_K40p_XeLS * h_K40p_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_K40p_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width) 
        //          + (Ratio_K40m_XeLS * h_K40m_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_K40m_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name=="Rate_Bi210_XeLS") {
        basename.push_back("Bi210m_XeLS");
        // if(_FitParameters[i]=="Rate_Bi210_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Bi210m_XeLS * h_Bi210m_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Bi210m_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name=="Rate_Kr85_XeLS") {
        basename.push_back("Kr85m_XeLS");
        // if(_FitParameters[i]=="Rate_Kr85_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Kr85m_XeLS * h_Kr85m_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Kr85m_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name=="Rate_Cs137_XeLS") {
        basename.push_back("Cs137m_gnd_XeLS");
        basename.push_back("Cs137m_1st_XeLS");
        basename.push_back("Cs137g_XeLS");
        // if(_FitParameters[i]=="Rate_Cs137_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Cs137m_gnd_XeLS * h_Cs137m_gnd_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Cs137m_gnd_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width) 
        //  + (Ratio_Cs137m_1st_XeLS * h_Cs137m_1st_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Cs137m_1st_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width)
        //  + (Ratio_Cs137g_XeLS * h_Cs137g_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Cs137g_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name=="Rate_Cs134_XeLS") {
        basename.push_back("Cs134m_XeLS");
        // if(_FitParameters[i]=="Rate_Cs134_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Cs134m_XeLS * h_Cs134m_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Cs134m_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name=="Rate_Bi208_XeLS") {
        basename.push_back("Bi208p_XeLS");
        // if(_FitParameters[i]=="Rate_Bi208_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Bi208p_XeLS * h_Bi208p_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Bi208p_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name=="Rate_Co60_XeLS") {
        basename.push_back("Co60m_XeLS");
        // if(_FitParameters[i]=="Rate_Co60_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Co60m_XeLS * h_Co60m_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Co60m_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name=="Rate_Y88_XeLS") {
        basename.push_back("Y88p_XeLS");
        // if(_FitParameters[i]=="Rate_Y88_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Y88p_XeLS * h_Y88p_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Y88p_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name=="Rate_Ag110_XeLS") {
        basename.push_back("Ag110m_XeLS");
        // if(_FitParameters[i]=="Rate_Ag110_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Ag110m_XeLS * h_Ag110m_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Ag110m_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name=="Rate_C10_XeLS") {
        basename.push_back("C10p_XeLS");
        // if(_FitParameters[i]=="Rate_C10_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_C10p_XeLS * h_C10p_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_C10p_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    } else if(name=="Rate_C11_XeLS") {
        basename.push_back("C11p_XeLS");
        // if(_FitParameters[i]=="Rate_C11_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_C11p_XeLS * h_C11p_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_C11p_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    } else if(name=="Rate_C11_XeLS") {
        basename.push_back("C11p_XeLS");
        // if(_FitParameters[i]=="Rate_C11_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_C11p_XeLS * h_C11p_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_C11p_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    } else if(name=="Rate_SolarNu_XeLS") {
        basename.push_back("SolarNu_XeLS");
        // if(_FitParameters[i]=="Rate_SolarNu_XeLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_SolarNu_XeLS * h_SolarNu_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_SolarNu_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    } else if(name=="Rate_Xe137_XeLS"){
        basename.push_back("Xe137m_XeLS");
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Xe137m_XeLS * h_Xe137m_XeLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Xe137m_XeLS) / _Exposure_XeLS_150cm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    } else if(name=="Rate_U238_S1_film"){
        basename.push_back("Pa234m_film");
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Pa234m_film * h_Pa234m_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Pa234m_film) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    } else if(name=="Rate_U238_S2_film") {
        basename.push_back("Pb214m_film");
        basename.push_back("Bi214m_film");
        tagging.push_back(false);
        tagging.push_back(true);
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Pb214m_film * h_Pb214m_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Pb214m_film) / _bin_width) 
        //     + (TaggingEfficiency_Bi214m_film * Ratio_Bi214m_film * h_Bi214m_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Bi214m_film) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    }else if(name.find("Rate_U238_S2_film_")!=string::npos) {
        string par_str = name;
        string target_str = "Rate_U238_S2_film_";
        par_str.erase(0, target_str.length());
        int k = atoi(par_str.c_str());

        if(k<0 || k>=thetaBin) {
            cerr << "ERROR : FitParameters format error (out of range) --> " << k << endl;
            return false;
        }

        basename.push_back("Pb214m_film");
        basename.push_back("Bi214m_film");
        tagging.push_back(false);
        tagging.push_back(true);
        theta.push_back(k);
        theta.push_back(k);

        // if(_FitParameters[i].find("Rate_U238_S2_film_")!=string::npos){
        // string par_str = _FitParameters[i];
        // string target_str = "Rate_U238_S2_film_";
        // par_str.erase(0, target_str.length());
        // int k = atoi(par_str.c_str());

        // if(k<0 || k>=thetaBin){
        //     cerr << "ERROR : FitParameters format error (out of range) --> " << k << endl;
        //     return false;
        // }

        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Pb214m_film * h_Pb214m_film_k[k]->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Pb214m_film) / _bin_width) 
        //     + (TaggingEfficiency_Bi214m_film * Ratio_Bi214m_film * h_Bi214m_film_k[k]->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Bi214m_film) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name=="Rate_Th232_S1_film") {
        basename.push_back("Ac228m_film");
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Ac228m_film * h_Ac228m_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Ac228m_film) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    }else if(name=="Rate_Th232_S2_film") {
        basename.push_back("Bi212m_film");
        basename.push_back("Pileup_film");
        basename.push_back("Tl208m_film");
        tagging.push_back(true);
        tagging.push_back(true);
        tagging.push_back(false);
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (TaggingEfficiency_Bi212m_film * Ratio_Bi212m_film * h_Bi212m_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Bi212m_film) / _bin_width)  
        //          + (TaggingEfficiency_Pileup_film * Ratio_Bi212m_film * h_Pileup_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Pileup_film) / _bin_width) 
        //          + (Ratio_Tl208m_film * h_Tl208m_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Tl208m_film) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    }else if(name.find("Rate_Th232_S2_film_")!=string::npos) {
        string par_str = name;
        string target_str = "Rate_Th232_S2_film_";
        par_str.erase(0, target_str.length());
        int k = atoi(par_str.c_str());

        if(k<0 || k>=thetaBin) {
            cerr << "ERROR : FitParameters format error (out of range) --> " << k << endl;
            return false;
        }

        basename.push_back("Bi212m_film");
        basename.push_back("Pileup_film");
        basename.push_back("Tl208m_film");
        tagging.push_back(true);
        tagging.push_back(true);
        tagging.push_back(false);
        theta.push_back(k);
        theta.push_back(k);
        theta.push_back(k);
    }else if(name=="Rate_K40_film") {
        basename.push_back("K40p_film");
        basename.push_back("K40m_film");
    }else if(name.find("Rate_K40_film_")!=string::npos) {
        string par_str = name;
        string target_str = "Rate_K40_film_";
        par_str.erase(0, target_str.length());
        int k = atoi(par_str.c_str());

        if(k<0 || k>=thetaBin) {
            cerr << "ERROR : FitParameters format error (out of range) --> " << k << endl;
            return false;
        }
        basename.push_back("K40p_film");
        basename.push_back("K40m_film");
        theta.push_back(k);
        theta.push_back(k);
    }else if(name=="Rate_Bi210_film") {
        basename.push_back("Bi210m_film");
        // if(_FitParameters[i]=="Rate_Bi210_film"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Bi210m_film * h_Bi210m_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Bi210m_film) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name.find("Rate_Bi210_film_")!=string::npos) {
        string par_str = name;
        string target_str = "Rate_Bi210_film_";
        par_str.erase(0, target_str.length());
        int k = atoi(par_str.c_str());

        if(k<0 || k>=thetaBin) {
            cerr << "ERROR : FitParameters format error (out of range) --> " << k << endl;
            return false;
        }

        basename.push_back("Bi210m_film");
        theta.push_back(k);
        // if(_FitParameters[i].find("Rate_Bi210_film_")!=string::npos){
        // string par_str = _FitParameters[i];
        // string target_str = "Rate_Bi210_film_";
        // par_str.erase(0, target_str.length());
        // int k = atoi(par_str.c_str());

        // if(k<0 || k>=thetaBin){
        //     cerr << "ERROR : FitParameters format error (out of range) --> " << k << endl;
        //     return false;
        // }

        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Bi210m_film * h_Bi210m_film_k[k]->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Bi210m_film) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name=="Rate_Cs137_film") {
        basename.push_back("Cs137m_gnd_film");
        basename.push_back("Cs137m_1st_film");
        basename.push_back("Cs137g_film");
        // for(int j=0; j<_nbin; j++) {
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Cs137m_gnd_film * h_Cs137m_gnd_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Cs137m_gnd_film) / _bin_width)
        //            + (Ratio_Cs137m_1st_film * h_Cs137m_1st_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Cs137m_1st_film) / _bin_width)
        //            + (Ratio_Cs137g_film * h_Cs137g_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Cs137g_film) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    }else if(name.find("Rate_Cs137_film_")!=string::npos) {
        string par_str = name;
        string target_str = "Rate_Cs137_film_";
        par_str.erase(0, target_str.length());
        int k = atoi(par_str.c_str());

        if(k<0 || k>=thetaBin) {
            cerr << "ERROR : FitParameters format error (out of range) --> " << k << endl;
            return false;
        }

        basename.push_back("Cs137m_gnd_film");
        basename.push_back("Cs137m_1st_film");
        basename.push_back("Cs137g_film");
        theta.push_back(k);
        theta.push_back(k);
        theta.push_back(k);

        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Cs137m_gnd_film * h_Cs137m_gnd_film_k[k]->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Cs137m_gnd_film) / _bin_width) 
        //  + (Ratio_Cs137m_1st_film * h_Cs137m_1st_film_k[k]->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Cs137m_1st_film) / _bin_width)
        //  + (Ratio_Cs137g_film * h_Cs137g_film_k[k]->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Cs137g_film) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name=="Rate_Cs134_film") {
        basename.push_back("Cs134m_film");
        // if(_FitParameters[i]=="Rate_Cs134_film"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Cs134m_film * h_Cs134m_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Cs134m_film) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name.find("Rate_Cs134_film_")!=string::npos) {
        string par_str = name;
        string target_str = "Rate_Cs134_film_";
        par_str.erase(0, target_str.length());
        int k = atoi(par_str.c_str());

        if(k<0 || k>=thetaBin) {
            cerr << "ERROR : FitParameters format error (out of range) --> " << k << endl;
            return false;
        }

        basename.push_back("Cs134m_film");
        theta.push_back(k);
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Cs134m_film * h_Cs134m_film_k[k]->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Cs134m_film) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    }else if(name=="Rate_Bi208_film") {
        basename.push_back("Bi208p_film");
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Bi208p_film * h_Bi208p_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Bi208p_film) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    }else if(name=="Rate_Co60_film") {
        basename.push_back("Co60m_film");
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Co60m_film * h_Co60m_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Co60m_film) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    }else if(name=="Rate_Y88_film") {
        basename.push_back("Y88p_film");
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Y88p_film * h_Y88p_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Y88p_film) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    }else if(name=="Rate_Ag110_film") {
        basename.push_back("Ag110m_film");
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Ag110m_film * h_Ag110m_film->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Ag110m_film) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    }else if(name.find("Rate_Ag110_film_")!=string::npos) {
        string par_str = name;
        string target_str = "Rate_Ag110_film_";
        par_str.erase(0, target_str.length());
        int k = atoi(par_str.c_str());

        if(k<0 || k>=thetaBin) {
            cerr << "ERROR : FitParameters format error (out of range) --> " << k << endl;
            return false;
        }

        basename.push_back("Ag110m_film");
        tagging.push_back(false);
        theta.push_back(k);
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Ag110m_film * h_Ag110m_film_k[k]->GetBinContent(j+1, n1+1, n2+1) / double(GenerateEvent_Ag110m_film) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    }else if(name=="Rate_U238_S2_KamLS") {
        basename.push_back("Pb214m_KamLS");
        basename.push_back("Bi214m_KamLS");
        tagging.push_back(false);
        tagging.push_back(true);
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Pb214m_KamLS * h_Pb214m_KamLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Pb214m_KamLS) / _Exposure_KamLS_norm * _ExposureBin[n]) / _bin_width) 
        //     + (TaggingEfficiency_Bi214m_KamLS * Ratio_Bi214m_KamLS * h_Bi214m_KamLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Bi214m_KamLS) / _Exposure_KamLS_norm * _ExposureBin[n]) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    }else if(name=="Rate_Th232_S2_KamLS") {
        basename.push_back("Bi212m_KamLS");
        basename.push_back("Tl208m_KamLS");
        tagging.push_back(true);
        tagging.push_back(false);
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (TaggingEfficiency_Bi212m_KamLS * Ratio_Bi212m_KamLS * h_Bi212m_KamLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Bi212m_KamLS) / _Exposure_KamLS_norm * _ExposureBin[n]) / _bin_width)
        //     + (Ratio_Tl208m_KamLS * h_Tl208m_KamLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Tl208m_KamLS) / _Exposure_KamLS_norm * _ExposureBin[n]) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    }else if(name=="Rate_K40_KamLS") {
        basename.push_back("K40p_KamLS");
        basename.push_back("K40m_KamLS");
        // if(_FitParameters[i]=="Rate_K40_KamLS"){
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_K40p_KamLS * h_K40p_KamLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_K40p_KamLS) / _Exposure_KamLS_norm * _ExposureBin[n]) / _bin_width) 
        //     + (Ratio_K40m_KamLS * h_K40m_KamLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_K40m_KamLS) / _Exposure_KamLS_norm * _ExposureBin[n]) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
        // }
    }else if(name.find("Rate_K40_KamLS_")!=string::npos) {
        string par_str = name;
        string target_str = "Rate_K40_KamLS_";
        par_str.erase(0, target_str.length());
        int k = atoi(par_str.c_str());

        if(k<0 || k>=thetaBin) {
            cerr << "ERROR : FitParameters format error (out of range) --> " << k << endl;
            return false;
        }

        basename.push_back("K40p_KamLS");
        basename.push_back("K40m_KamLS");
        theta.push_back(k);
        theta.push_back(k);
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_K40p_KamLS * h_K40p_KamLS_k[k]->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_K40p_KamLS) / _Exposure_KamLS_norm * _ExposureBin[n]) / _bin_width) 
        //     + (Ratio_K40m_KamLS * h_K40m_KamLS_k[k]->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_K40m_KamLS) / _Exposure_KamLS_norm * _ExposureBin[n]) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    }else if(name=="Rate_Bi210_KamLS") {
        basename.push_back("Bi210m_KamLS");
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Bi210m_KamLS * h_Bi210m_KamLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Bi210m_KamLS) / _Exposure_KamLS_norm * _ExposureBin[n]) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    }else if(name.find("Rate_Bi210_KamLS_")!=string::npos) {
        string par_str = name;
        string target_str = "Rate_Bi210_KamLS_";
        par_str.erase(0, target_str.length());
        int k = atoi(par_str.c_str());

        if(k<0 || k>=thetaBin) {
            cerr << "ERROR : FitParameters format error (out of range) --> " << k << endl;
            return false;
        }

        basename.push_back("Bi210m_KamLS");
        theta.push_back(k);
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Bi210m_KamLS * h_Bi210m_KamLS_k[k]->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Bi210m_KamLS) / _Exposure_KamLS_norm * _ExposureBin[n]) / _bin_width);
        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    }else if(name=="Rate_Kr85_KamLS") {
        basename.push_back("Kr85m_KamLS");
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_Kr85m_KamLS * h_Kr85m_KamLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_Kr85m_KamLS) / _Exposure_KamLS_norm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    } else if(name=="Rate_SolarNu_KamLS"){
        basename.push_back("SolarNu_KamLS");
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_SolarNu_KamLS * h_SolarNu_KamLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_SolarNu_KamLS) / _Exposure_KamLS_norm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    } else if(name=="Rate_C10_KamLS"){
        basename.push_back("C10p_KamLS");
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_C10p_KamLS * h_C10p_KamLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_C10p_KamLS) / _Exposure_KamLS_norm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    } else if(name=="Rate_C11_KamLS"){
        basename.push_back("C11p_KamLS");
        // N=0;
        // for(int j=0;j<_nbin;j++){
        //     x[N] = _min + _bin_width * double(j+0.5);
        //     y[N] = (Ratio_C11p_KamLS * h_C11p_KamLS->GetBinContent(j+1, n1+1, n2+1) / (double(GenerateEvent_C11p_KamLS) / _Exposure_KamLS_norm * _ExposureBin[n]) / _bin_width);

        //     N++;
        // }
        // _Spectrum_MC[n][i] = new kuFUNC(N,x,y);
    } else {
        //If the parameter is not found in this list
        return true;
    }
    return true;
} 

bool Parameter::BuildFiducialVolume() {

    char buffer[256];
    // multi Radius bin
    ifstream fRadiusBin(RadiusBin_name.c_str());
    if (!fRadiusBin) {
        cerr << "ERROR : cannot open input file" << endl;
        return false;
    }

    radiusBin = 0;
    while (fRadiusBin.getline(buffer, sizeof(buffer))) {
        istringstream strout(buffer);

        double radius_bin_min, radius_bin_max;
        int radius_bin_internal;

        if (!(strout >> radius_bin_min >> radius_bin_max >> radius_bin_internal)) {
            cerr << "ERROR : RadiusBin format error" << endl;
            return false;
        }

        // if radius_bin_internal == 1, then the bin is within the FV
        // if radius_bin_internal == 0, then the bin is outside the FV
        if (!(radius_bin_internal == 0 || radius_bin_internal == 1)) {
            cerr << "ERROR : RadiusBin format error" << endl;
            return false;
        }

        if (radiusBin >= 100) {
            cerr << "ERROR : RadiusBin format error (out of range)" << endl;
            return false;
        }

        radiusBin++;
    }

    fRadiusBin.close();

    // multi Theta bin
    ifstream fThetaBin(ThetaBin_name.c_str());
    if (!fThetaBin) {
        cerr << "ERROR : cannot open input file" << endl;
        return false;
    }

    thetaBin = 0;
    while (fThetaBin.getline(buffer, sizeof(buffer))) {
        istringstream strout(buffer);

        double theta_bin_min, theta_bin_max;
        int theta_bin_internal;

        if (!(strout >> theta_bin_min >> theta_bin_max >> theta_bin_internal)) {
            cerr << "ERROR : ThetaBin format error" << endl;
            return false;
        }

        if (!(theta_bin_internal == 0 || theta_bin_internal == 1)) {
            cerr << "ERROR : ThetaBin format error" << endl;
            return false;
        }

        if (thetaBin >= 10) {
            cerr << "ERROR : ThetaBin format error (out of range)" << endl;
            return false;
        }

        cerr << thetaBin << " :" << endl;
        cerr << endl;

        thetaBin++;
    }

    fThetaBin.close();

    return true;

}

Parameter::~Parameter() {
	//Destructor
}


double Parameter::GetD(const std::string& valType) {
    if ( doubleParamMap.find(valType) != doubleParamMap.end() ) {
	  return doubleParamMap[valType];
	} else {
		cerr << "Given " << valType << " value not found from double parameter map!" << endl;
		cerr<< "Avaliable Candidates are:";
		std::map<std::string, double>::iterator it;
		for (it = doubleParamMap.begin(); it != doubleParamMap.end(); it++)
		{	
			cerr<<it->first<<"|";
		}
		cerr<<endl;
	  abort();
	}
}

int Parameter::GetI(const std::string& valType) {
    if ( intParamMap.find(valType) != intParamMap.end() ) {
	  return intParamMap[valType];
	} else {
		cerr << "Given " << valType << " value not found from integer parameter map!" << endl;
		cerr<< "Avaliable Candidates are:";
		std::map<std::string, int>::iterator it;
		for (it = intParamMap.begin(); it != intParamMap.end(); it++)
		{	
			cerr<<it->first<<"|";
		}
		cerr<<endl;
	  abort();
	}
}

string Parameter::GetS(const std::string& valType) {
    if ( stringParamMap.find(valType) != stringParamMap.end() ) {
	  return stringParamMap[valType];
	} else if (valType=="unit") {
		return "";
	} else {
	  cerr << "Given " << valType << " value not found from string parameter map!" << endl;
	  	cerr<< "Avaliable Candidates are:";
		std::map<std::string, string>::iterator it;
		for (it = stringParamMap.begin(); it != stringParamMap.end(); it++)
		{	
			cerr<<it->first<<"|";
		}
		cerr<<endl;
	  abort();
	}
}

bool Parameter::GetB(const std::string& valType) {
    if ( boolParamMap.find(valType) != boolParamMap.end() ) {
	  return boolParamMap[valType];
	} else {
	  cerr << "Given " << valType << " value not found from boolean parameter map!" << endl;
	  abort();
	}
}

std::vector<string> Parameter::GetSArray(const std::string& valType) {
	std::map<std::string, std::string>::iterator it;
	std::vector<string> output;
	for (it = stringParamMap.begin(); it != stringParamMap.end(); it++)
    {
    	if (it->first.find(valType) != std::string::npos) {
    		output.push_back(it->second);
    	}

    }
    return output;
}

std::vector<int> Parameter::GetIArray(const std::string& valType) {
	std::map<std::string, int>::iterator it;
	std::vector<int> output;
	for (it = intParamMap.begin(); it != intParamMap.end(); it++)
    {
    	if (it->first.find(valType) != std::string::npos) {
    		output.push_back(it->second);
    	}

    }
    return output;
}

std::vector<bool> Parameter::GetBArray(const std::string& valType) {
	std::map<std::string, bool>::iterator it;
	std::vector<bool> output;
	for (it = boolParamMap.begin(); it != boolParamMap.end(); it++)
    {
    	if (it->first.find(valType) != std::string::npos) {
    		output.push_back(it->second);
    	}

    }
    return output;
}

std::vector<double> Parameter::GetDArray(const std::string& valType) {
	std::map<std::string, double>::iterator it;
	std::vector<double> output;
	for (it = doubleParamMap.begin(); it != doubleParamMap.end(); it++)
    {
    	if (it->first.find(valType) != std::string::npos) {
    		output.push_back(it->second);
    	}

    }
    return output;
}
