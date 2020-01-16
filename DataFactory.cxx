#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <DataFactory.h>
#include <stdio.h>
#include <stdlib.h>
#include "TSystem.h"
#include "TApplication.h"

#include "TCanvas.h"

#include "TNtuple.h"

#include "TH1.h"

#include "TH2.h"

#include "TH3.h"

#include "TMath.h"

#include "TF1.h"

#include "TStyle.h"

#include "TFunction.h"

#include "TObject.h"

#include "TPad.h"

#include "TText.h"

#include "TPostScript.h"

#include "TGraphErrors.h"

#include "TColor.h"

#include "TCut.h"

using namespace std;

DataFactory::DataFactory(ParameterFactory* Param) {
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

    SimulationSpectrumFile_XeLS = ((string) getenv("DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_FILE_XELS"));
    SimulationSpectrumFile_film = ((string) getenv("DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_FILE_FILM"));
    SimulationSpectrumFile_KamLS = ((string) getenv("DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_FILE_KAMLS"));
    RadiusBin_name = ((string) getenv("DOUBLEBETA_ANALYSIS_PARAMETER_DIR")) + "/RadiusBin.dat";
    ThetaBin_name = ((string) getenv("DOUBLEBETA_ANALYSIS_PARAMETER_DIR")) + "/ThetaBin.dat";
    Root_name = ((string) getenv("DOUBLEBETA_ANALYSIS_DATA_DIR")) + "/Single-DoubleBeta.root";

    rootf_mc_XeLS = TFile::Open(SimulationSpectrumFile_XeLS.c_str());
    rootf_mc_KamLS = TFile::Open(SimulationSpectrumFile_KamLS.c_str());
    rootf_mc_film = TFile::Open(SimulationSpectrumFile_film.c_str());

    fParameters = Param;
    _StartRun = fParameters->GetModelParameter("run_range")->GetI("low");
    _LastRun = fParameters->GetModelParameter("run_range")->GetI("high");
    _LiveTime = fParameters->GetModelParameter("total_livetime")->GetD("livetime");
    // double _NumberOfBiPo214;
    _R_Lower = fParameters->GetModelParameter("r_fv")->GetD("low");
    _R_Upper = fParameters->GetModelParameter("r_fv")->GetD("high");
    _Z_Lower = fParameters->GetModelParameter("z_fv")->GetD("low");
    _Z_Upper = fParameters->GetModelParameter("z_fv")->GetD("high");
    // double _DetectionEfficiency;
    // double _DeltaNormalize; //??
    // double _DeltaEnergy; //?

    Exposure_XeLS_All = fParameters->GetModelParameter("physics_constant")->GetD("Exposure_XeLS_All");
    Exposure_KamLS_All = fParameters->GetModelParameter("physics_constant")->GetD("Exposure_KamLS_All");
    Exposure_KamLS_norm = fParameters->GetModelParameter("physics_constant")->GetD("Exposure_KamLS_norm");
    Exposure_XeLS_150cm = fParameters->GetModelParameter("physics_constant")->GetD("Exposure_XeLS_150cm");

    if (!BuildFiducialVolume()) {
        cerr << "ERROR : Failed building Fiducial Volume!" << endl;
        abort();
    }

    if (!DataReader()) {
        cerr << "ERROR : Failed reading data!" << endl;
        abort();
    }

    gMCHistogram = MCReader();
}   

DataFactory::~DataFactory() {
    //Destructor

}

bool DataFactory::DataReader() {
    char run_cut[256];
    char run_cut2[256];
    char run_cut3[256];
    char r_cut[256];
    char r_cut2[256];
    char z_cut[256];
    char z_cut2[256];
    char r0_cut[256];
    char r0_cut2[256];
    char z0_cut[256];
    char z0_cut2[256];

    double alpha = 1.0;

    // For John's data:
    // sprintf(run_cut, "run_number>=%d && run_number<=%d", _StartRun, _LastRun);
    // sprintf(r_cut, "r>=%f && r<%f", _R_Lower * 100.0, _R_Upper * 100.0);
    // sprintf(r_cut2, "(r*%f)>=%f && (r*%f)<%f", alpha, _R_Lower * 100.0, alpha, _R_Upper * 100.0);
    // sprintf(z_cut, "z>=%f && z<%f", _Z_Lower * 100.0, _Z_Upper * 100.0);
    // sprintf(z_cut2, "(z*%f)>=%f && (z*%f)<%f", alpha, _Z_Lower * 100.0, alpha, _Z_Upper * 100.0);
    // sprintf(r0_cut, "r0>=%f && r0<%f", _R_Lower * 100.0, _R_Upper * 100.0);
    // sprintf(r0_cut2, "(r0*%f)>=%f && (r0*%f)<%f", alpha, _R_Lower * 100.0, alpha, _R_Upper * 100.0);
    // sprintf(z0_cut, "z0>=%f && z0<%f", _Z_Lower * 100.0, _Z_Upper * 100.0);
    // sprintf(z0_cut2, "(z0*%f)>=%f && (z0*%f)<%f", alpha, _Z_Lower * 100.0, alpha, _Z_Upper * 100.0);

    // For Touhoku Data:
    sprintf(run_cut, "run>=%d && run<%d", _StartRun, 15597);
    sprintf(run_cut2, "run>%d && run<%d", 15603, 15667);
    sprintf(run_cut3, "run>%d && run<=%d", 15688, _LastRun);
    sprintf(r_cut, "r>=%f && r<%f", _R_Lower * 100.0, _R_Upper * 100.0);
    sprintf(r_cut2, "(r*%f)>=%f && (r*%f)<%f", alpha, _R_Lower * 100.0, alpha, _R_Upper * 100.0);
    sprintf(z_cut, "z>=%f && z<%f", _Z_Lower * 100.0, _Z_Upper * 100.0);
    sprintf(z_cut2, "(z*%f)>=%f && (z*%f)<%f", alpha, _Z_Lower * 100.0, alpha, _Z_Upper * 100.0);

    sprintf(r0_cut, "r0>=%f && r0<%f", _R_Lower * 100.0, _R_Upper * 100.0);
    sprintf(r0_cut2, "(r0*%f)>=%f && (r0*%f)<%f", alpha, _R_Lower * 100.0, alpha, _R_Upper * 100.0);
    sprintf(z0_cut, "z0>=%f && z0<%f", _Z_Lower * 100.0, _Z_Upper * 100.0);
    sprintf(z0_cut2, "(z0*%f)>=%f && (z0*%f)<%f", alpha, _Z_Lower * 100.0, alpha, _Z_Upper * 100.0);

    TCut RunCut = run_cut;
    TCut RunCut2 = run_cut2;
    TCut RunCut3 = run_cut3;
    TCut RCut = r_cut;
    TCut RCut2 = r_cut2;
    TCut ZCut = z_cut;
    TCut ZCut2 = z_cut2;
    TCut R0Cut = r0_cut;
    TCut R0Cut2 = r0_cut2;
    TCut Z0Cut = z0_cut;
    TCut Z0Cut2 = z0_cut2;
    TCut DVeto = "Dveto==0";
    TCut C10Veto = "C10veto==0&&Trackveto==0";
    // TCut C10Veto = "C10veto==0";
    TCut RVeto = "Rveto==0";
    TCut RnVeto = "Rnveto==0";
    char hspot_data[256];
    double R_scale = fParameters->GetModelParameter("physics_constant")->GetD("R_scale");
    sprintf(hspot_data,"sqrt((z*%f+190.1)*(z*%f+190.1)+(x*x+y*y)*%f*%f)>70",R_scale,R_scale,R_scale,R_scale);

    TCut Hspot = hspot_data;
    //    TCut BadnessCut = "Badness<41.1*exp(-9.7*Evis)+2.31";
    //    TCut BadnessCut = "Badness<41.1*exp(-9.7*Evis)+2.31";
    TCut BadnessCut = "Badness<41.1*exp(-7.1*Evis)+2.7";  //

    string Epr="prompt_energy_multi*energy_corrected/(prompt_energy_multi+delayed_energy_multi)";
    string Ede="delayed_energy_multi*energy_corrected/(prompt_energy_multi+delayed_energy_multi)";
    string delT="prompt_to_delayed_time_multi";
    string delChi="chi2_multi-chi2_DP_multi";

    TCut EpMin ="prompt_energy_multi*energy_corrected/(prompt_energy_multi+delayed_energy_multi)>0.1";
    TCut EdMax ="delayed_energy_multi*energy_corrected/(prompt_energy_multi+delayed_energy_multi)<1.2";//11C
    TCut EdExp = (Ede + ">" + "0.7*exp(-1/5.5*(" + delT + "-13.5))+0.3*exp(-1/50.*(" + delT + "-13.5)) + 0.15").c_str();
    // TCut c_dp=EpMin&&EdMax&&EdExp;
    // TCut s_dp=!c_dp;
    TCut s_dp="Pileupveto==0";//new

    TFile *rootf = TFile::Open(Root_name.c_str());

    //The ntuple file of input events
    //TNtuple* nt = (TNtuple*) rootf->Get("klz_tree");
    TNtuple* nt = (TNtuple*) rootf->Get("tree");

    cerr << "set histogram:" << endl;

    for(int n=0; n<_NofEnergyBin; n++) {
        if(n+1==1000) {
            cerr << "ERROR : xbins is out of range" << endl;
            abort();
        }
        xbins[nbinx]   = _EnergyMin + _EnergyBinWidth * double(n+0.0);
        xbins[nbinx+1] = _EnergyMin + _EnergyBinWidth * double(n+1.0);
        nbinx++;
    }

    for(int n=0; n<_NofRadiusBin; n++) {
        if(n+1==100) {
            cerr << "ERROR : ybins is out of range" << endl;
            abort();
        }
        ybins[nbiny]   = _RadiusBinMin[n] * 100.0;
        ybins[nbiny+1] = _RadiusBinMax[n] * 100.0;
        nbiny++;
    }

    for(int n=0; n<_NofThetaBin; n++) {
        if(n+1==10) {
            cerr << "ERROR : zbins is out of range" << endl;
            abort();
        }
        zbins[nbinz]   = _ThetaBinMin[n];
        zbins[nbinz+1] = _ThetaBinMax[n];
        nbinz++;
    }

    //Constructing 3D histogram with energy, radius and theta bin.
    TH3D* h_Energy0 = new TH3D("h_Energy0","",nbinx, xbins, nbiny, ybins, nbinz, zbins);
    char hfill_data[256];
    sprintf(hfill_data,"z/r:r*%f:Evis>>h_Energy0",R_scale);

    //Draw TNtuple
    //nt->Draw("z/r:r:energy>>h_Energy0", RunCut, "goff");
    //nt->Draw("z/r:r:Evis>>h_Energy0", RunCut && ZCut && DVeto && MVeto && C10Veto && RVeto && RnVeto && BadnessCut, "goff");
    nt->Draw(hfill_data, (RunCut || RunCut2 || RunCut3) && DVeto  && C10Veto && RVeto && RnVeto && BadnessCut &&s_dp && Hspot, "goff");

    for(int n=0; n<_NofRadiusThetaBin; n++) {
        int n1 = n / _NofThetaBin;//radius bin
        int n2 = n % _NofThetaBin;//theta bin

        _E_Mean[n] = new double[_NofEnergyBin];//Average Energy
        _Data[n] = new double[_NofEnergyBin];//Data Counts
        _Data_e[n] = new double[_NofEnergyBin];//Error of Data

        _NofData = 0;
        for(int i=0; i<_NofEnergyBin; i++) {
            double E_Mean  = _EnergyMin + _EnergyBinWidth * double(i+0.5);//Bin Center
            double E_Lower = _EnergyMin + _EnergyBinWidth * double(i+0.0);//Bin Lower Edge
            double E_Upper = _EnergyMin + _EnergyBinWidth * double(i+1.0);//Bin Upper Edge

            double Obs = h_Energy0->GetBinContent(i+1, n1+1, n2+1);

            _E_Mean[n][_NofData] = E_Mean;
            _Data[n][_NofData] = Obs;
            _Data_e[n][_NofData] = sqrt(Obs);//Possion Error

            _NofData++;
        }
    }
    cerr << endl;

    rootf->Close();
    return true;
}

std::map<std::string,TH3D* > DataFactory::MCReader() {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Start Monte Carlo file input process

    std::map<std::string,TH3D* > map;
    //TFile *rootf_mc_by_photon_film = TFile::Open(SimulationSpectrumByPhotonFile_film.c_str());

    TTree* nt0_XeLS = (TNtuple*) rootf_mc_XeLS->Get("nt0");
    TTree* nt0_KamLS = (TNtuple*) rootf_mc_KamLS->Get("nt0");
    TTree* nt0_film = (TNtuple*) rootf_mc_film->Get("nt0");
    //TTree* nt0_by_photon_film = (TNtuple*) rootf_mc_by_photon_film->Get("nt0");

    // Xe136_0nu
    TTree* nt_Xe136_0nu_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_0nu_Xe136");

    // Xe136_2nu
    TTree* nt_Xe136_2nu_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_2nu_Xe136");

    // Xe136_2nu >2.0 MeV
    TTree* nt_Xe136H_2nu_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_2nu_Xe136H");

    // Xe136_2nu <2.0 MeV
    TTree* nt_Xe136L_2nu_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_2nu_Xe136L");

    // Xe134_2nu
    TTree* nt_Xe134_2nu_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_2nu_Xe134");

    // U238_S1
    TTree* nt_Pa234m_XeLS = (TNtuple*) rootf_mc_XeLS->Get("nt_Pa234m");

    // U238_S2
    TTree* nt_Pb214m_XeLS = (TNtuple*) rootf_mc_XeLS->Get("nt_Pb214m");
    TTree* nt_Bi214m_XeLS = (TNtuple*) rootf_mc_XeLS->Get("nt_Bi214m");

    // Th232_S1
    TTree* nt_Ac228m_XeLS = (TNtuple*) rootf_mc_XeLS->Get("nt_Ac228m");

    // Th232_S2
    TTree* nt_Bi212m_XeLS = (TNtuple*) rootf_mc_XeLS->Get("nt_Bi212m");
    TTree* nt_Pileup_XeLS = (TNtuple*) rootf_mc_XeLS->Get("nt_Pileup");
    TTree* nt_Tl208m_XeLS = (TNtuple*) rootf_mc_XeLS->Get("nt_Tl208m");

    // K40
    TTree* nt_K40p_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_K40p");
    TTree* nt_K40m_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_K40m");

    // Bi210
    TTree* nt_Bi210m_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_Bi210m");

    // Kr85
    TTree* nt_Kr85m_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_Kr85m");

    // Cs137
    TTree* nt_Cs137m_gnd_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_Cs137m_gnd");
    TTree* nt_Cs137m_1st_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_Cs137m_1st");
    TTree* nt_Cs137g_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_Cs137g");

    // Cs134
    TTree* nt_Cs134m_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_Cs134m");

    // Cs136
    TTree* nt_Cs136m_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_Cs136m");

    // Bi208
    TTree* nt_Bi208p_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_Bi208p");

    // Co60
    TTree* nt_Co60m_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_Co60m");

    // Y88
    TTree* nt_Y88p_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_Y88p");

    // Ag110
    TTree* nt_Ag110m_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_Ag110m");

    // C10
    TTree* nt_C10p_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_C10p");
    // C11
    TTree* nt_C11p_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_C11p");

    // B8Solar-nu
    TTree* nt_SolarNu_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_SolarB8ES");

    // Xe137m
    TTree* nt_Xe137m_XeLS  = (TNtuple*) rootf_mc_XeLS->Get("nt_Xe137m");

    // U238_S1_film
    TTree* nt_Pa234m_film = (TNtuple*) rootf_mc_film->Get("nt_Pa234m");

    // U238_S2_film
    TTree* nt_Pb214m_film = (TNtuple*) rootf_mc_film->Get("nt_Pb214m");
    //    TTree* nt_Bi214m_film = (TNtuple*) rootf_mc_by_photon_film->Get("nt_Bi214m");
    TTree* nt_Bi214m_film = (TNtuple*) rootf_mc_film->Get("nt_Bi214m");

    // Th232_S1_film
    TTree* nt_Ac228m_film = (TNtuple*) rootf_mc_film->Get("nt_Ac228m");

    // Th232_S2_film
    TTree* nt_Bi212m_film = (TNtuple*) rootf_mc_film->Get("nt_Bi212m");
    TTree* nt_Pileup_film = (TNtuple*) rootf_mc_film->Get("nt_Pileup");
    //    TTree* nt_Tl208m_film = (TNtuple*) rootf_mc_by_photon_film->Get("nt_Tl208m");
    TTree* nt_Tl208m_film = (TNtuple*) rootf_mc_film->Get("nt_Tl208m");

    // K40_film
    //TTree* nt_K40p_film  = (TNtuple*) rootf_mc_by_photon_film->Get("nt_K40p");
    //TTree* nt_K40m_film  = (TNtuple*) rootf_mc_by_photon_film->Get("nt_K40m");

    TTree* nt_K40p_film  = (TNtuple*) rootf_mc_film->Get("nt_K40p");
    TTree* nt_K40m_film  = (TNtuple*) rootf_mc_film->Get("nt_K40m");

    // Bi210_film
    TTree* nt_Bi210m_film  = (TNtuple*) rootf_mc_film->Get("nt_Bi210m");

    // Cs137_film
    TTree* nt_Cs137m_gnd_film  = (TNtuple*) rootf_mc_film->Get("nt_Cs137m_gnd");
    TTree* nt_Cs137m_1st_film  = (TNtuple*) rootf_mc_film->Get("nt_Cs137m_1st");
    TTree* nt_Cs137g_film  = (TNtuple*) rootf_mc_film->Get("nt_Cs137g");

    // Cs134_film
    //TTree* nt_Cs134m_film  = (TNtuple*) rootf_mc_by_photon_film->Get("nt_Cs134m");
    TTree* nt_Cs134m_film  = (TNtuple*) rootf_mc_film->Get("nt_Cs134m");

    // Bi208_film
    TTree* nt_Bi208p_film  = (TNtuple*) rootf_mc_film->Get("nt_Bi208p");

    // Co60_film
    TTree* nt_Co60m_film  = (TNtuple*) rootf_mc_film->Get("nt_Co60m");

    // Y88_film
    TTree* nt_Y88p_film  = (TNtuple*) rootf_mc_film->Get("nt_Y88p");

    // Ag110_film
    //TTree* nt_Ag110m_film  = (TNtuple*) rootf_mc_by_photon_film->Get("nt_Ag110m");
    TTree* nt_Ag110m_film  = (TNtuple*) rootf_mc_film->Get("nt_Ag110m");

    // U238_S2_KamLS
    TTree* nt_Pb214m_KamLS = (TNtuple*) rootf_mc_KamLS->Get("nt_Pb214m");
    TTree* nt_Bi214m_KamLS = (TNtuple*) rootf_mc_KamLS->Get("nt_Bi214m");

    // Th232_S2_KamLS
    TTree* nt_Bi212m_KamLS = (TNtuple*) rootf_mc_KamLS->Get("nt_Bi212m");
    TTree* nt_Tl208m_KamLS = (TNtuple*) rootf_mc_KamLS->Get("nt_Tl208m");

    // K40_KamLS
    TTree* nt_K40p_KamLS  = (TNtuple*) rootf_mc_KamLS->Get("nt_K40p");
    TTree* nt_K40m_KamLS  = (TNtuple*) rootf_mc_KamLS->Get("nt_K40m");

    // Bi210_KamLS
    TTree* nt_Bi210m_KamLS  = (TNtuple*) rootf_mc_KamLS->Get("nt_Bi210m");

    // Kr85_KamLS
    TTree* nt_Kr85m_KamLS  = (TNtuple*) rootf_mc_KamLS->Get("nt_Kr85m");

    // Solar-Nu_KamLS
    TTree* nt_SolarNu_KamLS  = (TNtuple*) rootf_mc_KamLS->Get("nt_SolarB8ES");

    // C10
    TTree* nt_C10p_KamLS  = (TNtuple*) rootf_mc_KamLS->Get("nt_C10p");

    // C11
    TTree* nt_C11p_KamLS  = (TNtuple*) rootf_mc_KamLS->Get("nt_C11p");

    double Ratio_Xe136_0nu_XeLS = 1.0;
    double Ratio_Xe136_2nu_XeLS = 1.0;
    double Ratio_Xe136L_2nu_XeLS = 0.99824952;
    double Ratio_Xe136H_2nu_XeLS = 0.00175048;
    double Ratio_Xe134_2nu_XeLS = 1.0;
    double Ratio_Pa234m_XeLS = 1.0;
    double Ratio_Pb214m_XeLS = 1.0;
    double Ratio_Bi214m_XeLS = 1.0;
    double Ratio_Ac228m_XeLS = 1.0;
    double Ratio_Bi212m_XeLS = 0.6406; // beta only
    double Ratio_Tl208m_XeLS = 0.3594;
    double Ratio_K40p_XeLS   = 0.1072;
    double Ratio_K40m_XeLS   = 0.8928;
    double Ratio_Bi210m_XeLS = 1.0;
    double Ratio_Kr85m_XeLS  = 1.0;
    double Ratio_Cs137m_gnd_XeLS = 0.056;
    double Ratio_Cs137m_1st_XeLS = 0.944;
    double Ratio_Cs137g_XeLS = 0.944;
    double Ratio_Cs134m_XeLS = 1.0;
    double Ratio_Bi208p_XeLS = 1.0;
    double Ratio_Co60m_XeLS = 1.0;
    double Ratio_Y88p_XeLS = 1.0;
    double Ratio_Ag110m_XeLS = 1.0;
    double Ratio_C10p_XeLS   = 1.0;
    double Ratio_C11p_XeLS   = 1.0;
    double Ratio_SolarNu_XeLS = 1.0;
    double Ratio_Xe137m_XeLS = 1.0;
    double Ratio_Cs136m_XeLS= 1.0;
    double Ratio_Pa234m_film = 1.0;
    double Ratio_Pb214m_film = 1.0;
    double Ratio_Bi214m_film = 1.0;
    double Ratio_Ac228m_film = 1.0;
    double Ratio_Bi212m_film = 0.6406; // beta only
    double Ratio_Tl208m_film = 0.3594;
    double Ratio_K40p_film   = 0.1072;
    double Ratio_K40m_film   = 0.8928;
    double Ratio_Bi210m_film = 1.0;
    double Ratio_Cs137m_gnd_film = 0.056;
    double Ratio_Cs137m_1st_film = 0.944;
    double Ratio_Cs137g_film = 0.944;
    double Ratio_Cs134m_film = 1.0;
    double Ratio_Bi208p_film = 1.0;
    double Ratio_Co60m_film = 1.0;
    double Ratio_Y88p_film = 1.0;
    double Ratio_Ag110m_film = 1.0;
    double Ratio_Pb214m_KamLS = 1.0;
    double Ratio_Bi214m_KamLS = 1.0;
    double Ratio_Bi212m_KamLS = 0.6406; // beta only
    double Ratio_Tl208m_KamLS = 0.3594;
    double Ratio_K40p_KamLS   = 0.1072;
    double Ratio_K40m_KamLS   = 0.8928;
    double Ratio_Bi210m_KamLS = 1.0;
    double Ratio_Kr85m_KamLS  = 1.0;
    double Ratio_SolarNu_KamLS = 1.0;
    double Ratio_C10p_KamLS   = 1.0;
    double Ratio_C11p_KamLS   = 1.0;

    gFitConstant["Ratio_Xe136_0nu_XeLS"] = 1.0;
    gFitConstant["Ratio_Xe136_2nu_XeLS"] = 1.0;
    gFitConstant["Ratio_Xe136L_2nu_XeLS"] = 0.99824952;
    gFitConstant["Ratio_Xe136H_2nu_XeLS"] = 0.00175048;
    gFitConstant["Ratio_Xe134_2nu_XeLS"] = 1.0;
    gFitConstant["Ratio_Pa234m_XeLS"] = 1.0;
    gFitConstant["Ratio_Pb214m_XeLS"] = 1.0;
    gFitConstant["Ratio_Bi214m_XeLS"] = 1.0;
    gFitConstant["Ratio_Ac228m_XeLS"] = 1.0;
    gFitConstant["Ratio_Bi212m_XeLS"] = 0.6406; // beta only
    gFitConstant["Ratio_Pileup_XeLS"] = 0.6406; // beta only
    gFitConstant["Ratio_Tl208m_XeLS"] = 0.3594;
    gFitConstant["Ratio_K40p_XeLS"] = 0.1072;
    gFitConstant["Ratio_K40m_XeLS"] = 0.8928;
    gFitConstant["Ratio_Bi210m_XeLS"] = 1.0;
    gFitConstant["Ratio_Kr85m_XeLS"] = 1.0;
    gFitConstant["Ratio_Cs137m_gnd_XeLS"] = 0.056;
    gFitConstant["Ratio_Cs137m_1st_XeLS"] = 0.944;
    gFitConstant["Ratio_Cs137g_XeLS"] = 0.944;
    gFitConstant["Ratio_Cs134m_XeLS"] = 1.0;
    gFitConstant["Ratio_Bi208p_XeLS"] = 1.0;
    gFitConstant["Ratio_Co60m_XeLS"] = 1.0;
    gFitConstant["Ratio_Y88p_XeLS"] = 1.0;
    gFitConstant["Ratio_Ag110m_XeLS"] = 1.0;
    gFitConstant["Ratio_C10p_XeLS"] = 1.0;
    gFitConstant["Ratio_C11p_XeLS"] = 1.0;
    gFitConstant["Ratio_SolarNu_XeLS"] = 1.0;
    gFitConstant["Ratio_Xe137m_XeLS"] = 1.0;
    gFitConstant["Ratio_Cs136m_XeLS"] = 1.0;
    gFitConstant["Ratio_Pa234m_film"] = 1.0;
    gFitConstant["Ratio_Pb214m_film"] = 1.0;
    gFitConstant["Ratio_Bi214m_film"] = 1.0;
    gFitConstant["Ratio_Ac228m_film"] = 1.0;
    gFitConstant["Ratio_Bi212m_film"] = 0.6406; // beta only
    gFitConstant["Ratio_Pileup_film"] = 0.6406; // beta only
    gFitConstant["Ratio_Tl208m_film"] = 0.3594;
    gFitConstant["Ratio_K40p_film"] = 0.1072;
    gFitConstant["Ratio_K40m_film"] = 0.8928;
    gFitConstant["Ratio_Bi210m_film"] = 1.0;
    gFitConstant["Ratio_Cs137m_gnd_film"] = 0.056;
    gFitConstant["Ratio_Cs137m_1st_film"] = 0.944;
    gFitConstant["Ratio_Cs137g_film"] = 0.944;
    gFitConstant["Ratio_Cs134m_film"] = 1.0;
    gFitConstant["Ratio_Bi208p_film"] = 1.0;
    gFitConstant["Ratio_Co60m_film"] = 1.0;
    gFitConstant["Ratio_Y88p_film"] = 1.0;
    gFitConstant["Ratio_Ag110m_film"] = 1.0;
    gFitConstant["Ratio_Pb214m_KamLS"] = 1.0;
    gFitConstant["Ratio_Bi214m_KamLS"] = 1.0;
    gFitConstant["Ratio_Bi212m_KamLS"] = 0.6406; // beta only
    gFitConstant["Ratio_Tl208m_KamLS"] = 0.3594;
    gFitConstant["Ratio_K40p_KamLS"] = 0.1072;
    gFitConstant["Ratio_K40m_KamLS"] = 0.8928;
    gFitConstant["Ratio_Bi210m_KamLS"] = 1.0;
    gFitConstant["Ratio_Kr85m_KamLS"] = 1.0;
    gFitConstant["Ratio_SolarNu_KamLS"] = 1.0;
    gFitConstant["Ratio_C10p_KamLS"] = 1.0;
    gFitConstant["Ratio_C11p_KamLS"] = 1.0;

    // double TaggingEfficiency_Bi214m_XeLS = 0.00028; // used until 2016.04.27
    double TaggingEfficiency_Bi214m_XeLS = 0.00046; // <- tag ineff.  tag eff: 0.99943 +- 0.00008
    //    double TaggingEfficiency_Pileup_XeLS = 0.0456; //( dT<20nsec )
    double TaggingEfficiency_Pileup_XeLS = 0.0214; //( dT<20nsec )

    //    double TaggingEfficiency_Pileup_XeLS = 0.0000; //( dT<20nsec )
    double TaggingEfficiency_Bi212m_XeLS = 0.0001;

    double TaggingEfficiency_Bi214m_KamLS = 0.00028; // temporary
    double TaggingEfficiency_Bi212m_KamLS = 0.0456*0.6; // temporary
    double TaggingEfficiency_Bi214m_film = 2.4 / (2.4+2.2);
    double TaggingEfficiency_Bi212m_film = 2.4 / (2.4+2.2); // assuming same with Bi214 // (249.5-145.9)/249.5 /new
    //    double TaggingEfficiency_Pileup_film = 0.0456; // ( dT<20nsec )
    double TaggingEfficiency_Pileup_film = 0.0211; // ( dT<20nsec )

    //  gFitConstant["TaggingEfficiency_Bi214m_XeLS"] = 0.00028; // used until 2016.04.27
     gFitConstant["TaggingEfficiency_Bi214m_XeLS"] = TaggingEfficiency_Bi214m_XeLS; // <- tag ineff.  tag eff: 0.99943 +- 0.00008
     gFitConstant["TaggingEfficiency_Pileup_XeLS"] = TaggingEfficiency_Pileup_XeLS; //( dT<20nsec )

    //     gFitConstant["TaggingEfficiency_Pileup_XeLS"] = 0.0000; //( dT<20nsec )
     gFitConstant["TaggingEfficiency_Bi212m_XeLS"] = TaggingEfficiency_Bi212m_XeLS;

     gFitConstant["TaggingEfficiency_Bi214m_KamLS"] = TaggingEfficiency_Bi214m_KamLS; // temporary
     gFitConstant["TaggingEfficiency_Bi212m_KamLS"] = TaggingEfficiency_Bi212m_KamLS; // temporary

     gFitConstant["TaggingEfficiency_Bi214m_film"] = TaggingEfficiency_Bi214m_film;
     gFitConstant["TaggingEfficiency_Pileup_film"] = TaggingEfficiency_Pileup_film; // ( dT<20nsec )
     gFitConstant["TaggingEfficiency_Bi212m_film"] = TaggingEfficiency_Bi212m_film; // assuming same with Bi214

    //    char r0cut[256];
    //    sprintf(r0cut,"r0*%f<150",R_scale);

    TCut R0150cm_cut="r0<150";
    TCut R0_cut_KamLS="r0>300&&r0<550&&z0/r0<0.5";

    //    int GenerateEvent_Xe136_0nu_XeLS_all;
    //    int GenerateEvent_Bi214m_XeLS_all;
    //    int GenerateEvent_Xe136_2nu_XeLS_all;
    //    int GenerateEvent_Bi214m_KamLS_all;
    //    nt0_KamLS->SetBranchAddress("GenerateEvent_Bi214m", &GenerateEvent_Bi214m_KamLS_all);

    //    nt0_KamLS->SetBranchAddress("GenerateEvent_Bi214m", &GenerateEvent_Bi214m_KamLS_all);
    //    nt0_KamLS->GetEntry(0);
    //    cout<< "Normilized : "<<((double)nt_Bi214m_KamLS->Draw("", R0_cut_KamLS))/_Exposure_KamLS_norm<<" all:  "<<(double)GenerateEvent_Bi214m_KamLS_all / _Exposure_KamLS_All<<endl;

    /*
    nt0_XeLS->SetBranchAddress("GenerateEvent_0nu_Xe136", &GenerateEvent_Xe136_0nu_XeLS_all);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Bi214m", &GenerateEvent_Bi214m_XeLS_all);
    nt0_XeLS->SetBranchAddress("GenerateEvent_2nu_Xe136", &GenerateEvent_Xe136_2nu_XeLS_all);
    nt0_XeLS->GetEntry(0);
    cout<< "150cm: "<<((double)nt_Xe136_0nu_XeLS->Draw("", R0150cm_cut))/_Exposure_XeLS_150cm<<" all:  "<<(double)GenerateEvent_Xe136_0nu_XeLS_all / _Exposure_XeLS_All<<endl;
    cout<< "150cm: "<<((double)nt_Xe136_2nu_XeLS->Draw("", R0150cm_cut))/_Exposure_XeLS_150cm<<" all:  "<<(double)GenerateEvent_Xe136_2nu_XeLS_all / _Exposure_XeLS_All<<endl;
    cout<< "150cm: "<<((double)nt_Bi214m_XeLS->Draw("", R0150cm_cut))/_Exposure_XeLS_150cm<<" all:  "<<(double)GenerateEvent_Bi214m_XeLS_all / _Exposure_XeLS_All<<endl;
    */

    int GenerateEvent_Xe136_0nu_XeLS  = nt_Xe136_0nu_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Xe136_2nu_XeLS  = nt_Xe136_2nu_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Xe136H_2nu_XeLS = nt_Xe136H_2nu_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Xe136L_2nu_XeLS = nt_Xe136L_2nu_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Xe134_2nu_XeLS = nt_Xe134_2nu_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Pa234m_XeLS    = nt_Pa234m_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Pb214m_XeLS    = nt_Pb214m_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Bi214m_XeLS    = nt_Bi214m_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Ac228m_XeLS    = nt_Ac228m_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Bi212m_XeLS    = nt_Bi212m_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Pileup_XeLS    = nt_Pileup_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Tl208m_XeLS    = nt_Tl208m_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_K40p_XeLS      = nt_K40p_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_K40m_XeLS      = nt_K40m_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Bi210m_XeLS    = nt_Bi210m_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Kr85m_XeLS     = nt_Kr85m_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Cs137m_gnd_XeLS    = nt_Cs137m_gnd_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Cs136m_XeLS        = nt_Cs136m_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Cs137m_1st_XeLS    = nt_Cs137m_1st_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Cs137g_XeLS    = nt_Cs137g_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Cs134m_XeLS    = nt_Cs134m_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Bi208p_XeLS    = nt_Bi208p_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Co60m_XeLS     = nt_Co60m_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Y88p_XeLS      = nt_Y88p_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Ag110m_XeLS    = nt_Ag110m_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_C10p_XeLS      = nt_C10p_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_C11p_XeLS      = nt_C11p_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_SolarNu_XeLS   = nt_SolarNu_XeLS->Draw("", R0150cm_cut);
    int GenerateEvent_Xe137m_XeLS   = nt_Xe137m_XeLS->Draw("", R0150cm_cut);
    /*    
    int GenerateEvent_Xe136_0nu_XeLS;
    int GenerateEvent_Xe136_2nu_XeLS;
    int GenerateEvent_Pa234m_XeLS;
    int GenerateEvent_Pb214m_XeLS;
    int GenerateEvent_Bi214m_XeLS;
    int GenerateEvent_Ac228m_XeLS;
    int GenerateEvent_Bi212m_XeLS;
    int GenerateEvent_Tl208m_XeLS;
    int GenerateEvent_K40p_XeLS;
    int GenerateEvent_K40m_XeLS;
    int GenerateEvent_Bi210m_XeLS;
    int GenerateEvent_Kr85m_XeLS;
    int GenerateEvent_Cs137m_gnd_XeLS;
    int GenerateEvent_Cs137m_1st_XeLS;
    int GenerateEvent_Cs137g_XeLS;
    int GenerateEvent_Cs134m_XeLS;
    int GenerateEvent_Bi208p_XeLS;
    int GenerateEvent_Co60m_XeLS;
    int GenerateEvent_Y88p_XeLS;
    int GenerateEvent_Ag110m_XeLS;
    int GenerateEvent_C10p_XeLS;
    int GenerateEvent_C11p_XeLS;
    int GenerateEvent_SolarNu_XeLS;
    */

    int GenerateEvent_Pa234m_film;
    int GenerateEvent_Pb214m_film;
    int GenerateEvent_Bi214m_film;
    int GenerateEvent_Ac228m_film;
    int GenerateEvent_Bi212m_film;
    int GenerateEvent_Pileup_film;
    int GenerateEvent_Tl208m_film;
    int GenerateEvent_K40p_film;
    int GenerateEvent_K40m_film;
    int GenerateEvent_Bi210m_film;
    int GenerateEvent_Cs137m_gnd_film;
    int GenerateEvent_Cs137m_1st_film;
    int GenerateEvent_Cs137g_film;
    int GenerateEvent_Cs134m_film;
    int GenerateEvent_Bi208p_film;
    int GenerateEvent_Co60m_film;
    int GenerateEvent_Y88p_film;
    int GenerateEvent_Ag110m_film;

    int GenerateEvent_Pb214m_KamLS   = nt_Pb214m_KamLS->Draw("", R0_cut_KamLS);
    int GenerateEvent_Bi214m_KamLS   = nt_Bi214m_KamLS->Draw("", R0_cut_KamLS);
    int GenerateEvent_Bi212m_KamLS   = nt_Bi212m_KamLS->Draw("", R0_cut_KamLS);
    int GenerateEvent_Tl208m_KamLS   = nt_Tl208m_KamLS->Draw("", R0_cut_KamLS);
    int GenerateEvent_K40p_KamLS     = nt_K40p_KamLS->Draw("", R0_cut_KamLS);
    int GenerateEvent_K40m_KamLS     = nt_K40m_KamLS->Draw("", R0_cut_KamLS);
    int GenerateEvent_Bi210m_KamLS   = nt_Bi210m_KamLS->Draw("", R0_cut_KamLS);
    int GenerateEvent_Kr85m_KamLS    = nt_Kr85m_KamLS->Draw("", R0_cut_KamLS);
    int GenerateEvent_SolarNu_KamLS  = nt_SolarNu_KamLS->Draw("", R0_cut_KamLS);
    int GenerateEvent_C10p_KamLS     = nt_C10p_KamLS->Draw("", R0_cut_KamLS);
    int GenerateEvent_C11p_KamLS     = nt_C11p_KamLS->Draw("", R0_cut_KamLS);


    /*
    nt0_XeLS->SetBranchAddress("GenerateEvent_0nu_Xe136", &GenerateEvent_Xe136_0nu_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_2nu_Xe136", &GenerateEvent_Xe136_2nu_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Pa234m", &GenerateEvent_Pa234m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Pb214m", &GenerateEvent_Pb214m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Bi214m", &GenerateEvent_Bi214m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Ac228m", &GenerateEvent_Ac228m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Bi212m", &GenerateEvent_Bi212m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Pileup", &GenerateEvent_Pileup_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Tl208m", &GenerateEvent_Tl208m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_K40p", &GenerateEvent_K40p_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_K40m", &GenerateEvent_K40m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Bi210m", &GenerateEvent_Bi210m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Kr85m", &GenerateEvent_Kr85m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Cs137m_gnd", &GenerateEvent_Cs137m_gnd_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Cs137m_1st", &GenerateEvent_Cs137m_1st_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Cs137g", &GenerateEvent_Cs137g_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Cs134m", &GenerateEvent_Cs134m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Bi208p", &GenerateEvent_Bi208p_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Co60m", &GenerateEvent_Co60m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Y88p", &GenerateEvent_Y88p_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Ag110m", &GenerateEvent_Ag110m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_C10p", &GenerateEvent_C10p_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_C11p", &GenerateEvent_C11p_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_SolarB8ES", &GenerateEvent_SolarNu_XeLS);*/

    //    cout<<GenerateEvent_Pb214m_KamLS<<" "<< GenerateEvent_Bi214m_KamLS 


    
    /*
    nt0_XeLS->SetBranchAddress("GenerateEvent_0nu_Xe136", &GenerateEvent_Xe136_0nu_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_2nu_Xe136", &GenerateEvent_Xe136_2nu_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Pa234m", &GenerateEvent_Pa234m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Pb214m", &GenerateEvent_Pb214m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Bi214m", &GenerateEvent_Bi214m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Ac228m", &GenerateEvent_Ac228m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Bi212m", &GenerateEvent_Bi212m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Pileup", &GenerateEvent_Pileup_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Tl208m", &GenerateEvent_Tl208m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_K40p", &GenerateEvent_K40p_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_K40m", &GenerateEvent_K40m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Bi210m", &GenerateEvent_Bi210m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Kr85m", &GenerateEvent_Kr85m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Cs137m_gnd", &GenerateEvent_Cs137m_gnd_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Cs137m_1st", &GenerateEvent_Cs137m_1st_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Cs137g", &GenerateEvent_Cs137g_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Cs134m", &GenerateEvent_Cs134m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Bi208p", &GenerateEvent_Bi208p_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Co60m", &GenerateEvent_Co60m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Y88p", &GenerateEvent_Y88p_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_Ag110m", &GenerateEvent_Ag110m_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_C10p", &GenerateEvent_C10p_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_C11p", &GenerateEvent_C11p_XeLS);
    nt0_XeLS->SetBranchAddress("GenerateEvent_SolarB8ES", &GenerateEvent_SolarNu_XeLS);*/

    nt0_film->SetBranchAddress("GenerateEvent_Pa234m", &GenerateEvent_Pa234m_film);
    nt0_film->SetBranchAddress("GenerateEvent_Pb214m", &GenerateEvent_Pb214m_film);
    //nt0_by_photon_film->SetBranchAddress("GenerateEvent_Bi214m", &GenerateEvent_Bi214m_film);
    nt0_film->SetBranchAddress("GenerateEvent_Bi214m", &GenerateEvent_Bi214m_film);
    nt0_film->SetBranchAddress("GenerateEvent_Ac228m", &GenerateEvent_Ac228m_film);
    nt0_film->SetBranchAddress("GenerateEvent_Bi212m", &GenerateEvent_Bi212m_film);
    nt0_film->SetBranchAddress("GenerateEvent_Pileup", &GenerateEvent_Pileup_film);
    //nt0_by_photon_film->SetBranchAddress("GenerateEvent_Tl208m", &GenerateEvent_Tl208m_film);
    nt0_film->SetBranchAddress("GenerateEvent_Tl208m", &GenerateEvent_Tl208m_film);
    //nt0_by_photon_film->SetBranchAddress("GenerateEvent_K40p", &GenerateEvent_K40p_film);
    //nt0_by_photon_film->SetBranchAddress("GenerateEvent_K40m", &GenerateEvent_K40m_film);
    nt0_film->SetBranchAddress("GenerateEvent_K40p", &GenerateEvent_K40p_film);
    nt0_film->SetBranchAddress("GenerateEvent_K40m", &GenerateEvent_K40m_film);
    nt0_film->SetBranchAddress("GenerateEvent_Bi210m", &GenerateEvent_Bi210m_film);
    nt0_film->SetBranchAddress("GenerateEvent_Cs137m_gnd", &GenerateEvent_Cs137m_gnd_film);
    nt0_film->SetBranchAddress("GenerateEvent_Cs137m_1st", &GenerateEvent_Cs137m_1st_film);
    nt0_film->SetBranchAddress("GenerateEvent_Cs137g", &GenerateEvent_Cs137g_film);
    //nt0_by_photon_film->SetBranchAddress("GenerateEvent_Cs134m", &GenerateEvent_Cs134m_film);
    nt0_film->SetBranchAddress("GenerateEvent_Cs134m", &GenerateEvent_Cs134m_film);
    nt0_film->SetBranchAddress("GenerateEvent_Bi208p", &GenerateEvent_Bi208p_film);
    nt0_film->SetBranchAddress("GenerateEvent_Co60m", &GenerateEvent_Co60m_film);
    nt0_film->SetBranchAddress("GenerateEvent_Y88p", &GenerateEvent_Y88p_film);
    //nt0_by_photon_film->SetBranchAddress("GenerateEvent_Ag110m", &GenerateEvent_Ag110m_film);
    nt0_film->SetBranchAddress("GenerateEvent_Ag110m", &GenerateEvent_Ag110m_film);
    /*    nt0_KamLS->SetBranchAddress("GenerateEvent_Pb214m", &GenerateEvent_Pb214m_KamLS);
    nt0_KamLS->SetBranchAddress("GenerateEvent_Bi214m", &GenerateEvent_Bi214m_KamLS);
    nt0_KamLS->SetBranchAddress("GenerateEvent_Bi212m", &GenerateEvent_Bi212m_KamLS);
    nt0_KamLS->SetBranchAddress("GenerateEvent_Tl208m", &GenerateEvent_Tl208m_KamLS);
    nt0_KamLS->SetBranchAddress("GenerateEvent_K40p", &GenerateEvent_K40p_KamLS);
    nt0_KamLS->SetBranchAddress("GenerateEvent_K40m", &GenerateEvent_K40m_KamLS);
    nt0_KamLS->SetBranchAddress("GenerateEvent_Bi210m", &GenerateEvent_Bi210m_KamLS);
    nt0_KamLS->SetBranchAddress("GenerateEvent_Kr85m", &GenerateEvent_Kr85m_KamLS);
    nt0_KamLS->SetBranchAddress("GenerateEvent_SolarB8ES", &GenerateEvent_SolarNu_KamLS);
    nt0_KamLS->SetBranchAddress("GenerateEvent_C10p", &GenerateEvent_C10p_KamLS);
    nt0_KamLS->SetBranchAddress("GenerateEvent_C11p", &GenerateEvent_C11p_KamLS);*/

    nt0_film->GetEntry(0);

    gFitConstant["GenerateEvent_Xe136_0nu_XeLS"] =   (double)GenerateEvent_Xe136_0nu_XeLS;
    gFitConstant["GenerateEvent_Xe136_2nu_XeLS"] =   (double)GenerateEvent_Xe136_2nu_XeLS;
    gFitConstant["GenerateEvent_Xe136H_2nu_XeLS"] =   (double)GenerateEvent_Xe136H_2nu_XeLS;//
    gFitConstant["GenerateEvent_Xe136L_2nu_XeLS"] =   (double)GenerateEvent_Xe136L_2nu_XeLS;//
    gFitConstant["GenerateEvent_Xe134_2nu_XeLS"] =   (double)GenerateEvent_Xe134_2nu_XeLS;
    gFitConstant["GenerateEvent_Pa234m_XeLS"] =   (double)GenerateEvent_Pa234m_XeLS;
    gFitConstant["GenerateEvent_Pb214m_XeLS"] =   (double)GenerateEvent_Pb214m_XeLS;
    gFitConstant["GenerateEvent_Bi214m_XeLS"] =   (double)GenerateEvent_Bi214m_XeLS;
    gFitConstant["GenerateEvent_Ac228m_XeLS"] =   (double)GenerateEvent_Ac228m_XeLS;
    gFitConstant["GenerateEvent_Bi212m_XeLS"] =   (double)GenerateEvent_Bi212m_XeLS;
    gFitConstant["GenerateEvent_Pileup_XeLS"] =   (double)GenerateEvent_Pileup_XeLS;
    gFitConstant["GenerateEvent_Tl208m_XeLS"] =   (double)GenerateEvent_Tl208m_XeLS;
    gFitConstant["GenerateEvent_K40p_XeLS"] =   (double)GenerateEvent_K40p_XeLS;
    gFitConstant["GenerateEvent_K40m_XeLS"] =   (double)GenerateEvent_K40m_XeLS;
    gFitConstant["GenerateEvent_Bi210m_XeLS"] =   (double)GenerateEvent_Bi210m_XeLS;
    gFitConstant["GenerateEvent_Kr85m_XeLS"] =   (double)GenerateEvent_Kr85m_XeLS;
    gFitConstant["GenerateEvent_Cs137m_gnd_XeLS"] =   (double)GenerateEvent_Cs137m_gnd_XeLS;
    gFitConstant["GenerateEvent_Cs137m_1st_XeLS"] =   (double)GenerateEvent_Cs137m_1st_XeLS;
    gFitConstant["GenerateEvent_Cs137g_XeLS"] =   (double)GenerateEvent_Cs137g_XeLS;
    gFitConstant["GenerateEvent_Cs134m_XeLS"] =   (double)GenerateEvent_Cs134m_XeLS;
    gFitConstant["GenerateEvent_Bi208p_XeLS"] =   (double)GenerateEvent_Bi208p_XeLS;
    gFitConstant["GenerateEvent_Co60m_XeLS"] =   (double)GenerateEvent_Co60m_XeLS;
    gFitConstant["GenerateEvent_Y88p_XeLS"] =   (double)GenerateEvent_Y88p_XeLS;
    gFitConstant["GenerateEvent_Ag110m_XeLS"] =   (double)GenerateEvent_Ag110m_XeLS;
    gFitConstant["GenerateEvent_C10p_XeLS"] =   (double)GenerateEvent_C10p_XeLS;
    gFitConstant["GenerateEvent_C11p_XeLS"] =   (double)GenerateEvent_C11p_XeLS;
    gFitConstant["GenerateEvent_SolarNu_XeLS"] =   (double)GenerateEvent_SolarNu_XeLS;
    gFitConstant["GenerateEvent_Xe137m_XeLS"] =   (double)GenerateEvent_Xe137m_XeLS;

    gFitConstant["GenerateEvent_Xe136_0nu_XeLS"] =   (double)GenerateEvent_Xe136_0nu_XeLS;
    gFitConstant["GenerateEvent_Xe136_2nu_XeLS"] =   (double)GenerateEvent_Xe136_2nu_XeLS;
    gFitConstant["GenerateEvent_Pa234m_XeLS"] =   (double)GenerateEvent_Pa234m_XeLS;
    gFitConstant["GenerateEvent_Pb214m_XeLS"] =   (double)GenerateEvent_Pb214m_XeLS;
    gFitConstant["GenerateEvent_Bi214m_XeLS"] =   (double)GenerateEvent_Pb214m_XeLS;
    gFitConstant["GenerateEvent_Ac228m_XeLS"] =   (double)GenerateEvent_Ac228m_XeLS;
    gFitConstant["GenerateEvent_Ac228m_XeLS"] =   (double)GenerateEvent_Ac228m_XeLS;
    gFitConstant["GenerateEvent_Tl208m_XeLS"] =   (double)GenerateEvent_Tl208m_XeLS;
    gFitConstant["GenerateEvent_K40p_XeLS"] =   (double)GenerateEvent_K40p_XeLS;
    gFitConstant["GenerateEvent_K40m_XeLS"] =   (double)GenerateEvent_K40m_XeLS;
    gFitConstant["GenerateEvent_Bi210m_XeLS"] =   (double)GenerateEvent_Bi210m_XeLS;
    gFitConstant["GenerateEvent_Kr85m_XeLS"] =   (double)GenerateEvent_Kr85m_XeLS;
    gFitConstant["GenerateEvent_Cs137m_gnd_XeLS"] =   (double)GenerateEvent_Cs137m_gnd_XeLS;
    gFitConstant["GenerateEvent_Cs137m_1st_XeLS"] =   (double)GenerateEvent_Cs137m_1st_XeLS;
    gFitConstant["GenerateEvent_Cs137g_XeLS"] =   (double)GenerateEvent_Cs137g_XeLS;
    gFitConstant["GenerateEvent_Cs134m_XeLS"] =   (double)GenerateEvent_Cs134m_XeLS;
    gFitConstant["GenerateEvent_Bi208p_XeLS"] =   (double)GenerateEvent_Bi208p_XeLS;
    gFitConstant["GenerateEvent_Co60m_XeLS"] =   (double)GenerateEvent_Co60m_XeLS;
    gFitConstant["GenerateEvent_Y88p_XeLS"] =   (double)GenerateEvent_Y88p_XeLS;
    gFitConstant["GenerateEvent_Ag110m_XeLS"] =   (double)GenerateEvent_Ag110m_XeLS;
    gFitConstant["GenerateEvent_C10p_XeLS"] =   (double)GenerateEvent_C10p_XeLS;
    gFitConstant["GenerateEvent_C11p_XeLS"] =   (double)GenerateEvent_C11p_XeLS;
    gFitConstant["GenerateEvent_SolarNu_XeLS"] =   (double)GenerateEvent_SolarNu_XeLS;

    gFitConstant["GenerateEvent_Pa234m_film"] =   (double)GenerateEvent_Pa234m_film;
    gFitConstant["GenerateEvent_Pb214m_film"] =   (double)GenerateEvent_Pb214m_film;
    gFitConstant["GenerateEvent_Bi214m_film"] =   (double)GenerateEvent_Bi214m_film;
    gFitConstant["GenerateEvent_Ac228m_film"] =   (double)GenerateEvent_Ac228m_film;
    gFitConstant["GenerateEvent_Bi212m_film"] =   (double)GenerateEvent_Bi212m_film;
    gFitConstant["GenerateEvent_Pileup_film"] =   (double)GenerateEvent_Pileup_film;
    gFitConstant["GenerateEvent_Tl208m_film"] =   (double)GenerateEvent_Tl208m_film;
    gFitConstant["GenerateEvent_K40p_film"] =   (double)GenerateEvent_K40p_film;
    gFitConstant["GenerateEvent_K40m_film"] =   (double)GenerateEvent_K40m_film;
    gFitConstant["GenerateEvent_Bi210m_film"] =   (double)GenerateEvent_Bi210m_film;
    gFitConstant["GenerateEvent_Cs137m_gnd_film"] =   (double)GenerateEvent_Cs137m_gnd_film;
    gFitConstant["GenerateEvent_Cs137m_1st_film"] =   (double)GenerateEvent_Cs137m_1st_film;
    gFitConstant["GenerateEvent_Cs137g_film"] =   (double)GenerateEvent_Cs137g_film;
    gFitConstant["GenerateEvent_Cs134m_film"] =   (double)GenerateEvent_Cs134m_film;
    gFitConstant["GenerateEvent_Bi208p_film"] =   (double)GenerateEvent_Bi208p_film;
    gFitConstant["GenerateEvent_Co60m_film"] =   (double)GenerateEvent_Co60m_film;
    gFitConstant["GenerateEvent_Y88p_film"] =   (double)GenerateEvent_Y88p_film;
    gFitConstant["GenerateEvent_Ag110m_film"] =   (double)GenerateEvent_Ag110m_film;

    gFitConstant["GenerateEvent_Pb214m_KamLS"] =   (double)GenerateEvent_Pb214m_KamLS;
    gFitConstant["GenerateEvent_Bi214m_KamLS"] =   (double)GenerateEvent_Bi214m_KamLS;
    gFitConstant["GenerateEvent_Bi212m_KamLS"] =   (double)GenerateEvent_Bi212m_KamLS;
    gFitConstant["GenerateEvent_Tl208m_KamLS"] =   (double)GenerateEvent_Tl208m_KamLS;
    gFitConstant["GenerateEvent_K40p_KamLS"] =   (double)GenerateEvent_K40p_KamLS;
    gFitConstant["GenerateEvent_K40m_KamLS"] =   (double)GenerateEvent_K40m_KamLS;
    gFitConstant["GenerateEvent_Bi210m_KamLS"] =   (double)GenerateEvent_Bi210m_KamLS;
    gFitConstant["GenerateEvent_Kr85m_KamLS"] =   (double)GenerateEvent_Kr85m_KamLS;
    gFitConstant["GenerateEvent_SolarNu_KamLS"] =   (double)GenerateEvent_SolarNu_KamLS;
    gFitConstant["GenerateEvent_C10p_KamLS"] =   (double)GenerateEvent_C10p_KamLS;
    gFitConstant["GenerateEvent_C11p_KamLS"] =   (double)GenerateEvent_C11p_KamLS;

    //Making xbins
    for(int n=0; n<_nbin; n++) {
        if(n+1==1000) {
            cerr << "ERROR : x1bins is out of range" << endl;
            abort();
        }
        x1bins[nbinx1]   = _min + _bin_width * double(n+0.0);
        x1bins[nbinx1+1] = _min + _bin_width * double(n+1.0);
        nbinx1++;
    }


    char hspot[256];
    sprintf(hspot,"sqrt((z+190.1)*(z+190.1)+(x*x+y*y))>70");
    Hspot = hspot;

    // (r, energy) --- reconstructed vertex
    const bool fill_MC_Spectrum = true;

    //Filling in reconstructed vertex position and energy
    // (r, energy) --- reconstructed vertex
    TH3D* h_Xe136_0nu_XeLS = new TH3D("h_Xe136_0nu_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Xe136_0nu_XeLS->Draw("z/r:r:EnergyA2>>h_Xe136_0nu_XeLS", Hspot,"goff");
    map["h_Xe136_0nu_XeLS"] = h_Xe136_0nu_XeLS;

    TH3D* h_Xe136_2nu_XeLS = new TH3D("h_Xe136_2nu_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Xe136_2nu_XeLS->Draw("z/r:r:EnergyA2>>h_Xe136_2nu_XeLS", Hspot,"goff");
    map["h_Xe136_2nu_XeLS"] = h_Xe136_2nu_XeLS;

    TH3D* h_Xe136H_2nu_XeLS = new TH3D("h_Xe136H_2nu_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Xe136H_2nu_XeLS->Draw("z/r:r:EnergyA2>>h_Xe136H_2nu_XeLS", Hspot,"goff");
    map["h_Xe136H_2nu_XeLS"] = h_Xe136H_2nu_XeLS;

    TH3D* h_Xe136L_2nu_XeLS = new TH3D("h_Xe136L_2nu_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Xe136L_2nu_XeLS->Draw("z/r:r:EnergyA2>>h_Xe136L_2nu_XeLS", Hspot,"goff");
    map["h_Xe136L_2nu_XeLS"] = h_Xe136L_2nu_XeLS;


    TH3D* h_Xe134_2nu_XeLS = new TH3D("h_Xe134_2nu_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Xe134_2nu_XeLS->Draw("z/r:r:EnergyA2>>h_Xe134_2nu_XeLS", Hspot,"goff");
    map["h_Xe134_2nu_XeLS"] = h_Xe134_2nu_XeLS;

    TH3D* h_Pa234m_XeLS = new TH3D("h_Pa234m_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Pa234m_XeLS->Draw("z/r:r:EnergyA2>>h_Pa234m_XeLS", Hspot,"goff");
    map["h_Pa234m_XeLS"] = h_Pa234m_XeLS;

    TH3D* h_Pb214m_XeLS = new TH3D("h_Pb214m_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Pb214m_XeLS->Draw("z/r:r:EnergyA2>>h_Pb214m_XeLS", Hspot,"goff");
    map["h_Pb214m_XeLS"] = h_Pb214m_XeLS;

    TH3D* h_Bi214m_XeLS = new TH3D("h_Bi214m_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Bi214m_XeLS->Draw("z/r:r:EnergyA2>>h_Bi214m_XeLS", Hspot,"goff");
    map["h_Bi214m_XeLS"] = h_Bi214m_XeLS;

    TH3D* h_Ac228m_XeLS = new TH3D("h_Ac228m_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Ac228m_XeLS->Draw("z/r:r:EnergyA2>>h_Ac228m_XeLS", Hspot,"goff");
    map["h_Ac228m_XeLS"] = h_Ac228m_XeLS;

    TH3D* h_Bi212m_XeLS = new TH3D("h_Bi212m_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Bi212m_XeLS->Draw("z/r:r:EnergyA2>>h_Bi212m_XeLS", Hspot,"goff");
    map["h_Bi212m_XeLS"] = h_Bi212m_XeLS;

    TH3D* h_Tl208m_XeLS = new TH3D("h_Tl208m_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Tl208m_XeLS->Draw("z/r:r:EnergyA2>>h_Tl208m_XeLS", Hspot,"goff");
    map["h_Tl208m_XeLS"] = h_Tl208m_XeLS;

    TH3D* h_Pileup_XeLS = new TH3D("h_Pileup_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Pileup_XeLS->Draw("z/r:r:(EnergyA2)>>h_Pileup_XeLS", Hspot,"goff");
    map["h_Pileup_XeLS"] = h_Pileup_XeLS;

    TH3D* h_K40p_XeLS = new TH3D("h_K40p_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_K40p_XeLS->Draw("z/r:r:EnergyA2>>h_K40p_XeLS", Hspot,"goff");
    map["h_K40p_XeLS"] = h_K40p_XeLS;

    TH3D* h_K40m_XeLS = new TH3D("h_K40m_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_K40m_XeLS->Draw("z/r:r:EnergyA2>>h_K40m_XeLS", Hspot,"goff");
    map["h_K40m_XeLS"] = h_K40m_XeLS;

    TH3D* h_Bi210m_XeLS = new TH3D("h_Bi210m_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Bi210m_XeLS->Draw("z/r:r:EnergyA2>>h_Bi210m_XeLS", Hspot,"goff");
    map["h_Bi210m_XeLS"] = h_Bi210m_XeLS;

    TH3D* h_Kr85m_XeLS = new TH3D("h_Kr85m_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Kr85m_XeLS->Draw("z/r:r:EnergyA2>>h_Kr85m_XeLS", Hspot,"goff");
    map["h_Kr85m_XeLS"] = h_Kr85m_XeLS;

    TH3D* h_Cs137m_gnd_XeLS = new TH3D("h_Cs137m_gnd_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Cs137m_gnd_XeLS->Draw("z/r:r:EnergyA2>>h_Cs137m_gnd_XeLS", Hspot,"goff");
    map["h_Cs137m_gnd_XeLS"] = h_Cs137m_gnd_XeLS;

    TH3D* h_Cs137m_1st_XeLS = new TH3D("h_Cs137m_1st_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Cs137m_1st_XeLS->Draw("z/r:r:EnergyA2>>h_Cs137m_1st_XeLS", Hspot,"goff");
    map["h_Cs137m_1st_XeLS"] = h_Cs137m_1st_XeLS;

    TH3D* h_Cs137g_XeLS = new TH3D("h_Cs137g_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Cs137g_XeLS->Draw("z/r:r:EnergyA2>>h_Cs137g_XeLS", Hspot,"goff");
    map["h_Cs137g_XeLS"] = h_Cs137g_XeLS;

    TH3D* h_Cs134m_XeLS = new TH3D("h_Cs134m_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Cs134m_XeLS->Draw("z/r:r:EnergyA2>>h_Cs134m_XeLS", Hspot,"goff");
    map["h_Cs134m_XeLS"] = h_Cs134m_XeLS;

    TH3D* h_Bi208p_XeLS = new TH3D("h_Bi208p_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Bi208p_XeLS->Draw("z/r:r:EnergyA2>>h_Bi208p_XeLS", Hspot,"goff");
    map["h_Bi208p_XeLS"] = h_Bi208p_XeLS;

    TH3D* h_Co60m_XeLS = new TH3D("h_Co60m_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Co60m_XeLS->Draw("z/r:r:EnergyA2>>h_Co60m_XeLS", Hspot,"goff");
    map["h_Co60m_XeLS"] = h_Co60m_XeLS;

    TH3D* h_Y88p_XeLS = new TH3D("h_Y88p_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Y88p_XeLS->Draw("z/r:r:EnergyA2>>h_Y88p_XeLS", Hspot,"goff");
    map["h_Y88p_XeLS"] = h_Y88p_XeLS;

    TH3D* h_Ag110m_XeLS = new TH3D("h_Ag110m_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Ag110m_XeLS->Draw("z/r:r:EnergyA2>>h_Ag110m_XeLS", Hspot,"goff");
    map["h_Ag110m_XeLS"] = h_Ag110m_XeLS;

    TH3D* h_C10p_XeLS = new TH3D("h_C10p_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_C10p_XeLS->Draw("z/r:r:EnergyA2>>h_C10p_XeLS", Hspot,"goff");
    map["h_C10p_XeLS"] = h_C10p_XeLS;


    TH3D* h_C11p_XeLS = new TH3D("h_C11p_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_C11p_XeLS->Draw("z/r:r:EnergyA2>>h_C11p_XeLS", Hspot,"goff");
    map["h_C11p_XeLS"] = h_C11p_XeLS;

    TH3D* h_SolarNu_XeLS = new TH3D("h_SolarNu_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_SolarNu_XeLS->Draw("z/r:r:EnergyA2>>h_SolarNu_XeLS", Hspot,"goff");
    map["h_SolarNu_XeLS"] = h_SolarNu_XeLS;


    TH3D* h_Xe137m_XeLS = new TH3D("h_Xe137m_XeLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Xe137m_XeLS->Draw("z/r:r:EnergyA2>>h_Xe137m_XeLS", Hspot,"goff");
    map["h_Xe137m_XeLS"] = h_Xe137m_XeLS;

    TH3D* h_Pa234m_film = new TH3D("h_Pa234m_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Pa234m_film->Draw("z/r:r:EnergyA2>>h_Pa234m_film", Hspot,"goff");
    map["h_Pa234m_film"] = h_Pa234m_film;

    TH3D* h_Pb214m_film = new TH3D("h_Pb214m_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Pb214m_film->Draw("z/r:r:EnergyA2>>h_Pb214m_film", Hspot,"goff");
    map["h_Pb214m_film"] = h_Pb214m_film;

    //These seems to be dealing with background that has a theta variation.
    //Instead of directly dumping everything into one histogram,
    //It applies a cut on z/r with respect to the current Theta Bin and then dump events.
    TH3D* h_Pb214m_film_k[10];
    for(int k=0; k<_NofThetaBin; k++) {
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_Pb214m_film_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_Pb214m_film_%d", k);
        //This seems to cut on the current theta bin. There are 10 theta bins in total.
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;
        h_Pb214m_film_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_Pb214m_film->Draw(varexp, Theta0Cut && Hspot,"goff");
        map[title] = h_Pb214m_film_k[k];
    }

    TH3D* h_Bi214m_film = new TH3D("h_Bi214m_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Bi214m_film->Draw("z/r:r:EnergyA2>>h_Bi214m_film", Hspot,"goff");
    map["h_Bi214m_film"] = h_Bi214m_film;

    TH3D* h_Bi214m_film_k[10];
    for(int k=0; k<_NofThetaBin; k++) {
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_Bi214m_film_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_Bi214m_film_%d", k);
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;

        h_Bi214m_film_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_Bi214m_film->Draw(varexp, Theta0Cut && Hspot,"goff");
        map[title] = h_Bi214m_film_k[k];
    }

    TH3D* h_Ac228m_film = new TH3D("h_Ac228m_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Ac228m_film->Draw("z/r:r:EnergyA2>>h_Ac228m_film", Hspot,"goff");
    map["h_Ac228m_film"] = h_Ac228m_film;

    TH3D* h_Bi212m_film = new TH3D("h_Bi212m_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Bi212m_film->Draw("z/r:r:EnergyA2>>h_Bi212m_film", Hspot,"goff");
    map["h_Bi212m_film"] = h_Bi212m_film;

    TH3D* h_Pileup_film = new TH3D("h_Pileup_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Pileup_film->Draw("z/r:r:EnergyA2>>h_Pileup_film", Hspot,"goff");
    map["h_Pileup_film"] = h_Pileup_film;

    TH3D* h_Bi212m_film_k[10];
    for(int k=0; k<_NofThetaBin; k++) {
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_Bi212m_film_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_Bi212m_film_%d", k);
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;

        h_Bi212m_film_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_Bi212m_film->Draw(varexp, Theta0Cut && Hspot,"goff");
        map[title] = h_Bi212m_film_k[k];
    }

    TH3D* h_Pileup_film_k[10];
    for(int k=0;k<_NofThetaBin;k++){
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_Pileup_film_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_Pileup_film_%d", k);
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;

        h_Pileup_film_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_Pileup_film->Draw(varexp, Theta0Cut&&Hspot,"goff");
        map[title] = h_Pileup_film_k[k];
    }

    TH3D* h_Tl208m_film = new TH3D("h_Tl208m_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Tl208m_film->Draw("z/r:r:EnergyA2>>h_Tl208m_film", Hspot,"goff");
    map["h_Tl208m_film"] = h_Tl208m_film;

    TH3D* h_Tl208m_film_k[10];
    for(int k=0; k<_NofThetaBin; k++) {
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_Tl208m_film_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_Tl208m_film_%d", k);
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;

        h_Tl208m_film_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_Tl208m_film->Draw(varexp, Theta0Cut && Hspot,"goff");
        map[title] = h_Tl208m_film_k[k];
    }

    TH3D* h_K40p_film = new TH3D("h_K40p_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_K40p_film->Draw("z/r:r:EnergyA2>>h_K40p_film", Hspot,"goff");
    map["h_K40p_film"] = h_K40p_film;

    TH3D* h_K40p_film_k[10];
    for(int k=0; k<_NofThetaBin; k++) {
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_K40p_film_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_K40p_film_%d", k);
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;

        h_K40p_film_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_K40p_film->Draw(varexp, Theta0Cut && Hspot,"goff");
        map[title] = h_K40p_film_k[k];
    }

    TH3D* h_K40m_film = new TH3D("h_K40m_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_K40m_film->Draw("z/r:r:EnergyA2>>h_K40m_film", Hspot,"goff");
    map["h_K40m_film"] = h_K40m_film;

    TH3D* h_K40m_film_k[10];
    for(int k=0; k<_NofThetaBin; k++) {
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_K40m_film_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_K40m_film_%d", k);
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;

        h_K40m_film_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_K40m_film->Draw(varexp, Theta0Cut && Hspot,"goff");
        map[title] = h_K40m_film_k[k];
    }

    TH3D* h_Bi210m_film = new TH3D("h_Bi210m_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Bi210m_film->Draw("z/r:r:EnergyA2>>h_Bi210m_film", Hspot,"goff");
    map["h_Bi210m_film"] = h_Bi210m_film;

    TH3D* h_Bi210m_film_k[10];
    for(int k=0; k<_NofThetaBin; k++) {
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_Bi210m_film_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_Bi210m_film_%d", k);
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;

        h_Bi210m_film_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_Bi210m_film->Draw(varexp, Theta0Cut && Hspot,"goff");
        map[title] = h_Bi210m_film_k[k];
    }

    TH3D* h_Cs137m_gnd_film = new TH3D("h_Cs137m_gnd_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Cs137m_gnd_film->Draw("z/r:r:EnergyA2>>h_Cs137m_gnd_film", Hspot,"goff");
    map["h_Cs137m_gnd_film"] = h_Cs137m_gnd_film;

    TH3D* h_Cs137m_gnd_film_k[10];
    for(int k=0; k<_NofThetaBin; k++) {
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_Cs137m_gnd_film_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_Cs137m_gnd_film_%d", k);
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;

        h_Cs137m_gnd_film_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_Cs137m_gnd_film->Draw(varexp, Theta0Cut && Hspot,"goff");
        map[title] = h_Cs137m_gnd_film_k[k];
    }

    TH3D* h_Cs137m_1st_film = new TH3D("h_Cs137m_1st_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Cs137m_1st_film->Draw("z/r:r:EnergyA2>>h_Cs137m_1st_film", Hspot,"goff");
    map["h_Cs137m_1st_film"] = h_Cs137m_1st_film;

    TH3D* h_Cs137m_1st_film_k[10];
    for(int k=0; k<_NofThetaBin; k++) {
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_Cs137m_1st_film_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_Cs137m_1st_film_%d", k);
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;

        h_Cs137m_1st_film_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_Cs137m_1st_film->Draw(varexp, Theta0Cut && Hspot,"goff");
        map[title] = h_Cs137m_1st_film_k[k];
    }

    TH3D* h_Cs137g_film = new TH3D("h_Cs137g_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Cs137g_film->Draw("z/r:r:EnergyA2>>h_Cs137g_film", Hspot,"goff");
    map["h_Cs137g_film"] = h_Cs137g_film;

    TH3D* h_Cs137g_film_k[10];
    for(int k=0; k<_NofThetaBin; k++) {
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_Cs137g_film_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_Cs137g_film_%d", k);
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;

        h_Cs137g_film_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_Cs137g_film->Draw(varexp, Theta0Cut && Hspot,"goff");
        map[title] = h_Cs137g_film_k[k];
    }

    TH3D* h_Cs134m_film = new TH3D("h_Cs134m_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Cs134m_film->Draw("z/r:r:EnergyA2>>h_Cs134m_film", Hspot,"goff");
    map["h_Cs134m_film"] = h_Cs134m_film;

    TH3D* h_Cs134m_film_k[10];
    for(int k=0; k<_NofThetaBin; k++) {
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_Cs134m_film_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_Cs134m_film_%d", k);
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;

        h_Cs134m_film_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_Cs134m_film->Draw(varexp, Theta0Cut && Hspot,"goff");
        map[title] = h_Cs134m_film_k[k];
    }

    TH3D* h_Bi208p_film = new TH3D("h_Bi208p_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Bi208p_film->Draw("z/r:r:EnergyA2>>h_Bi208p_film", Hspot,"goff");
    map["h_Bi208p_film"] = h_Bi208p_film;

    TH3D* h_Co60m_film = new TH3D("h_Co60m_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Co60m_film->Draw("z/r:r:EnergyA2>>h_Co60m_film", Hspot,"goff");
    map["h_Co60m_film"] = h_Co60m_film;

    TH3D* h_Y88p_film = new TH3D("h_Y88p_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Y88p_film->Draw("z/r:r:EnergyA2>>h_Y88p_film", Hspot,"goff");
    map["h_Y88p_film"] = h_Y88p_film;

    TH3D* h_Ag110m_film = new TH3D("h_Ag110m_film", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Ag110m_film->Draw("z/r:r:EnergyA2>>h_Ag110m_film", Hspot,"goff");
    map["h_Ag110m_film"] = h_Ag110m_film;

    TH3D* h_Ag110m_film_k[10];
    for(int k=0; k<_NofThetaBin; k++) {
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_Ag110m_film_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_Ag110m_film_%d", k);
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;

        h_Ag110m_film_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_Ag110m_film->Draw(varexp, Theta0Cut,"goff");
        map[title] = h_Ag110m_film_k[k];
    }

    TH3D* h_Pb214m_KamLS = new TH3D("h_Pb214m_KamLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Pb214m_KamLS->Draw("z/r:r:EnergyA2>>h_Pb214m_KamLS", Hspot,"goff");
    map["h_Pb214m_KamLS"] = h_Pb214m_KamLS;

    TH3D* h_Bi214m_KamLS = new TH3D("h_Bi214m_KamLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Bi214m_KamLS->Draw("z/r:r:EnergyA2>>h_Bi214m_KamLS", Hspot,"goff");
    map["h_Bi214m_KamLS"] = h_Bi214m_KamLS;

    TH3D* h_Bi212m_KamLS = new TH3D("h_Bi212m_KamLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Bi212m_KamLS->Draw("z/r:r:EnergyA2>>h_Bi212m_KamLS", Hspot,"goff");
    map["h_Bi212m_KamLS"] = h_Bi212m_KamLS;

    TH3D* h_Tl208m_KamLS = new TH3D("h_Tl208m_KamLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Tl208m_KamLS->Draw("z/r:r:EnergyA2>>h_Tl208m_KamLS", Hspot,"goff");
    map["h_Tl208m_KamLS"] = h_Tl208m_KamLS;

    TH3D* h_K40p_KamLS = new TH3D("h_K40p_KamLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_K40p_KamLS->Draw("z/r:r:EnergyA2>>h_K40p_KamLS", Hspot,"goff");
    map["h_K40p_KamLS"] = h_K40p_KamLS;

    TH3D* h_K40p_KamLS_k[10];
    for(int k=0; k<_NofThetaBin; k++) {
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_K40p_KamLS_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_K40p_KamLS_%d", k);
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;

        h_K40p_KamLS_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_K40p_KamLS->Draw(varexp, Theta0Cut && Hspot,"goff");
        map[title] = h_K40p_KamLS_k[k];
    }

    TH3D* h_K40m_KamLS = new TH3D("h_K40m_KamLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_K40m_KamLS->Draw("z/r:r:EnergyA2>>h_K40m_KamLS", Hspot,"goff");
    map["h_K40m_KamLS"] = h_K40m_KamLS;

    TH3D* h_K40m_KamLS_k[10];
    for(int k=0; k<_NofThetaBin; k++) {
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_K40m_KamLS_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_K40m_KamLS_%d", k);
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;

        h_K40m_KamLS_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_K40m_KamLS->Draw(varexp, Theta0Cut && Hspot,"goff");
        map[title] = h_K40m_KamLS_k[k];
    }

    TH3D* h_Bi210m_KamLS = new TH3D("h_Bi210m_KamLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Bi210m_KamLS->Draw("z/r:r:EnergyA2>>h_Bi210m_KamLS", Hspot,"goff");
    map["h_Bi210m_KamLS"] = h_Bi210m_KamLS;

    TH3D* h_Bi210m_KamLS_k[10];
    for(int k=0; k<_NofThetaBin; k++) {
        char title[256];
        char varexp[256];
        char theta0_cut[256];

        sprintf(title, "h_Bi210m_KamLS_%d", k);
        sprintf(varexp, "z/r:r:EnergyA2>>h_Bi210m_KamLS_%d", k);
        sprintf(theta0_cut, "z0/r0>=%f && z0/r0<%f", _ThetaBinMin[k], _ThetaBinMax[k]);

        TCut Theta0Cut = theta0_cut;

        h_Bi210m_KamLS_k[k] = new TH3D(title, "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
        if(fill_MC_Spectrum==true) nt_Bi210m_KamLS->Draw(varexp, Theta0Cut && Hspot,"goff");\
        map[title] = h_Bi210m_KamLS_k[k];
    }

    TH3D* h_Kr85m_KamLS = new TH3D("h_Kr85m_KamLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Kr85m_KamLS->Draw("z/r:r:EnergyA2>>h_Kr85m_KamLS", Hspot,"goff");
    map["h_Kr85m_KamLS"] = h_Kr85m_KamLS;

    TH3D* h_SolarNu_KamLS = new TH3D("h_SolarNu_KamLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_SolarNu_KamLS->Draw("z/r:r:EnergyA2>>h_SolarNu_KamLS", Hspot,"goff");
    map["h_SolarNu_KamLS"] = h_SolarNu_KamLS;

    TH3D* h_C10p_KamLS = new TH3D("h_C10p_KamLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_C10p_KamLS->Draw("z/r:r:EnergyA2>>h_C10p_KamLS", Hspot,"goff");
    map["h_C10p_KamLS"] = h_C10p_KamLS;

    TH3D* h_C11p_KamLS = new TH3D("h_C11p_KamLS", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_C11p_KamLS->Draw("z/r:r:EnergyA2>>h_C11p_KamLS", Hspot,"goff");
    map["h_C11p_KamLS"] = h_C11p_KamLS;

    //This hsitogram fills in the value r0 and z0 instead of r and z:
    //What does r0 and z0 mean?
    // (r0, energy) --- generated vertex
    TH3D* h_Xe136_0nu_XeLS_Generate = new TH3D("h_Xe136_0nu_XeLS_Generate", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Xe136_0nu_XeLS->Draw("z0/r0:r0:EnergyA2>>h_Xe136_0nu_XeLS_Generate", "","goff");
    map["h_Xe136_0nu_XeLS_Generate"] = h_Xe136_0nu_XeLS_Generate;

    TH3D* h_Xe136_2nu_XeLS_Generate = new TH3D("h_Xe136_2nu_XeLS_Generate", "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);
    if(fill_MC_Spectrum==true) nt_Xe136_2nu_XeLS->Draw("z0/r0:r0:EnergyA2>>h_Xe136_2nu_XeLS_Generate", "","goff");
    map["h_Xe136_2nu_XeLS_Generate"] = h_Xe136_2nu_XeLS_Generate;
    //End of Setting up MC file
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Start of checking the uniformity of events
    cerr << "check MC event uniformity in Xe-LS (Xe136_0nu, Xe136_2nu):" << endl;
    for(int n=0; n<_NofRadiusThetaBin; n++) {
        int n1 = n / _NofThetaBin;//radius bin number
        int n2 = n % _NofThetaBin;//theta bin number

        int EventBin_Xe136_0nu_XeLS = 0;
        int GenerateEventBin_Xe136_0nu_XeLS = 0;
        int EventBin_Xe136_2nu_XeLS = 0;
        int GenerateEventBin_Xe136_2nu_XeLS = 0;
        for(int j=-1; j<=_nbin; j++) { // include data in underflow and overflow bin
            EventBin_Xe136_0nu_XeLS += int(h_Xe136_0nu_XeLS->GetBinContent(j+1, n1+1, n2+1));
            GenerateEventBin_Xe136_0nu_XeLS += int(h_Xe136_0nu_XeLS_Generate->GetBinContent(j+1, n1+1, n2+1));
            EventBin_Xe136_2nu_XeLS += int(h_Xe136_2nu_XeLS->GetBinContent(j+1, n1+1, n2+1));
            GenerateEventBin_Xe136_2nu_XeLS += int(h_Xe136_2nu_XeLS_Generate->GetBinContent(j+1, n1+1, n2+1));
        }

        cerr<<GenerateEvent_Xe136_0nu_XeLS<<"    "<<GenerateEvent_Xe136_2nu_XeLS<<endl;

        //This calculate the ratio between number of events and number of generated events
        //Still need to understand what is "generated"
        double EventRatio_Xe136_0nu_XeLS = double(EventBin_Xe136_0nu_XeLS) / double(GenerateEvent_Xe136_0nu_XeLS);
        double EventRatio_Xe136_0nu_XeLS_e = sqrt(double(EventBin_Xe136_0nu_XeLS)) / double(GenerateEvent_Xe136_0nu_XeLS);
        double EventRatio_Xe136_2nu_XeLS = double(EventBin_Xe136_2nu_XeLS) / double(GenerateEvent_Xe136_2nu_XeLS);
        double EventRatio_Xe136_2nu_XeLS_e = sqrt(double(EventBin_Xe136_2nu_XeLS)) / double(GenerateEvent_Xe136_2nu_XeLS);
        double GenerateEventRatio_Xe136_0nu_XeLS = double(GenerateEventBin_Xe136_0nu_XeLS) / double(GenerateEvent_Xe136_0nu_XeLS);
        double GenerateEventRatio_Xe136_0nu_XeLS_e = sqrt(double(GenerateEventBin_Xe136_0nu_XeLS)) / double(GenerateEvent_Xe136_0nu_XeLS);
        double GenerateEventRatio_Xe136_2nu_XeLS = double(GenerateEventBin_Xe136_2nu_XeLS) / double(GenerateEvent_Xe136_2nu_XeLS);
        double GenerateEventRatio_Xe136_2nu_XeLS_e = sqrt(double(GenerateEventBin_Xe136_2nu_XeLS)) / double(GenerateEvent_Xe136_2nu_XeLS);
        double VolumeRatio = _ExposureBin[n] / fParameters->GetModelParameter("physics_constant")->GetD("Exposure_XeLS_All");//This measures the exposure in nth bin vs. total exposure

        //This measure the percentage deviation of GenerateEventRatio/VolumeRatio, it should be 1
        double GenerateEventDeviation_Xe136_0nu_XeLS = (GenerateEventRatio_Xe136_0nu_XeLS / VolumeRatio - 1.0) * 100.0;
        double GenerateEventDeviation_Xe136_0nu_XeLS_e = (GenerateEventRatio_Xe136_0nu_XeLS_e / VolumeRatio) * 100.0;
        double GenerateEventDeviation_Xe136_2nu_XeLS = (GenerateEventRatio_Xe136_2nu_XeLS / VolumeRatio - 1.0) * 100.0;
        double GenerateEventDeviation_Xe136_2nu_XeLS_e = (GenerateEventRatio_Xe136_2nu_XeLS_e / VolumeRatio) * 100.0;

        //SpaceEfficiency is defined by EventBin/GenerateEventBin
        double SpaceEfficiency_Xe136_0nu_XeLS = 0;
        double SpaceEfficiency_Xe136_0nu_XeLS_e = 0;
        double SpaceEfficiency_Xe136_2nu_XeLS = 0;
        double SpaceEfficiency_Xe136_2nu_XeLS_e = 0;

        if(GenerateEventBin_Xe136_0nu_XeLS>0) {
            SpaceEfficiency_Xe136_0nu_XeLS = double(EventBin_Xe136_0nu_XeLS) / double(GenerateEventBin_Xe136_0nu_XeLS);
            SpaceEfficiency_Xe136_0nu_XeLS_e = sqrt(double(EventBin_Xe136_0nu_XeLS)) / double(GenerateEventBin_Xe136_0nu_XeLS);
        }

        if(GenerateEventBin_Xe136_2nu_XeLS>0) {
            SpaceEfficiency_Xe136_2nu_XeLS = double(EventBin_Xe136_2nu_XeLS) / double(GenerateEventBin_Xe136_2nu_XeLS);
            SpaceEfficiency_Xe136_2nu_XeLS_e = sqrt(double(EventBin_Xe136_2nu_XeLS)) / double(GenerateEventBin_Xe136_2nu_XeLS);
        }

        cerr << n << " :" << endl;
        cerr << "R (m) = " << _RadiusBinMin[n1] << " - " << _RadiusBinMax[n1] << endl;
        cerr << "CosTheta = " << _ThetaBinMin[n2] << " - " << _ThetaBinMax[n2] << endl;
        cerr << "VolumeRatio = " << VolumeRatio << endl;
        cerr << endl;

        cerr << "EventRatio_Xe136_0nu_XeLS = " << EventRatio_Xe136_0nu_XeLS << " +- " << EventRatio_Xe136_0nu_XeLS_e << endl;
        cerr << "EventRatio_Xe136_2nu_XeLS = " << EventRatio_Xe136_2nu_XeLS << " +- " << EventRatio_Xe136_2nu_XeLS_e << endl;
        cerr << "GenerateEventRatio_Xe136_0nu_XeLS = " << GenerateEventRatio_Xe136_0nu_XeLS << " +- " << GenerateEventRatio_Xe136_0nu_XeLS_e << endl;
        cerr << "GenerateEventRatio_Xe136_2nu_XeLS = " << GenerateEventRatio_Xe136_2nu_XeLS << " +- " << GenerateEventRatio_Xe136_2nu_XeLS_e << endl;
        cerr << endl;

        cerr << "GenerateEventDeviation_Xe136_0nu_XeLS [%] = " << GenerateEventDeviation_Xe136_0nu_XeLS << " +- " << GenerateEventDeviation_Xe136_0nu_XeLS_e << endl;
        cerr << "GenerateEventDeviation_Xe136_2nu_XeLS [%] = " << GenerateEventDeviation_Xe136_2nu_XeLS << " +- " << GenerateEventDeviation_Xe136_2nu_XeLS_e << endl;
        cerr << endl;

        cerr << "SpaceEfficiency_Xe136_0nu_XeLS [nominal] = " << SpaceEfficiency_Xe136_0nu_XeLS << " +- " << SpaceEfficiency_Xe136_0nu_XeLS_e << endl;
        cerr << "SpaceEfficiency_Xe136_2nu_XeLS [nominal] = " << SpaceEfficiency_Xe136_2nu_XeLS << " +- " << SpaceEfficiency_Xe136_2nu_XeLS_e << endl;
        cerr << endl;
    }

    cerr << "note:" << endl;
    cerr << "- deviation (event_ratio/volume_ratio) should be small in R < 1.54 m, and if deviation errors are too large, statistics of MC event should be increased!!!" << endl;
    cerr << "- nominal space efficiency is close to 1.0, but efficiency loss of 2nu can be smaller due to worse vertex resolution in low energy events (below energy threshold)." << endl;
    cerr << endl;

    //cerr << "for film rate check:" << endl;
    //cerr << "GenerateEvent_Bi214m_film = " << GenerateEvent_Bi214m_film << endl;
    //cerr << "GenerateEvent_K40p_film   = " << GenerateEvent_K40p_film   << endl;
    //cerr << "GenerateEvent_K40m_film   = " << GenerateEvent_K40m_film   << endl;
    //cerr << "GenerateEvent_Cs134m      = " << GenerateEvent_Cs134m      << endl;
    //cerr << "GenerateEvent_Ag110m      = " << GenerateEvent_Ag110m      << enld;
    //cerr << endl;
    //End of Uniformity Check
    //End of fit parameter setting
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    return map;
}
double DataFactory::getFiducialVolume(double R_Lower, double R_Upper, double Theta_Lower, double Theta_Upper) {
    double FiducialVolume = 4.0 / 3.0 * M_PI * (pow(R_Upper, 3) - pow(R_Lower, 3)) * (Theta_Upper - Theta_Lower) / 2.0;
    
    if(R_Upper>1.901-0.70 && Theta_Lower<-sqrt(190.1*190.1-70.*70.)/190.1 ){
    
      if( Theta_Upper<-sqrt(190.1*190.1-70.*70.)/190.1 || Theta_Lower!=-1.0){
        cerr <<" Wrong theta cut: theta_upper "<<Theta_Upper<<" should be > "<<-sqrt(190.1*190.1-70.*70.)/190.1<<endl;
        cerr <<" Or"<<Theta_Lower<<" should be > "<<-1<<endl;
        abort();
      }

      double cos=(1.901*1.901+R_Upper*R_Upper-0.70*0.70)/(2.*1.901*R_Upper);
      double z=R_Upper*cos;
      double V1=M_PI*(R_Upper*R_Upper*z-1./3.*pow(z, 3)+2./3.*pow(R_Upper, 3));
      double zd=z-1.901;
      double V2=M_PI*(0.70*0.70*zd-1./3.*pow(zd, 3)+2./3.*pow(0.70, 3));
      double V_Upper=V1-V2;

      if(R_Upper>1.901) V_Upper = 4.0 / 3.0 * M_PI * (pow(R_Upper, 3) - pow(0.70, 3));

      double V_Lower=4.0/3.0 *M_PI * pow(R_Lower, 3);

      if(R_Lower>1.901-0.70 ){
        double cos_l=(1.901*1.901+R_Lower*R_Lower-0.70*0.70)/(2.*1.901*R_Lower);
        double z_l=R_Lower*cos;
        double V1_l=M_PI*(R_Lower*R_Lower*z_l-1./3.*pow(z_l, 3)+2./3.*pow(R_Lower, 3));
        double zd_l=z_l-1.901;
        double V2_l=M_PI*(0.70*0.70*zd-1./3.*pow(zd, 3)+2./3.*pow(0.70, 3));
        V_Lower=V1_l-V2_l;  
      }
      if(R_Lower>1.901) V_Lower = 4.0 / 3.0 * M_PI * (pow(R_Lower, 3) - pow(0.70, 3));

      FiducialVolume = (V_Upper - V_Lower) - (4.0 / 3.0 * M_PI * (pow(R_Upper, 3) - pow(R_Lower, 3)) * (1. - Theta_Upper) / 2.0);
    
      //      cout<<" V1 "<<V_Upper<<" V2 "<<V_Lower<< " V1-V2 " << V_Upper-V_Lower<<endl;
    }
    
    return FiducialVolume;
}

bool DataFactory::BuildFiducialVolume() {
    cerr << "=================================================================" << endl;
    cerr << "Set Radius Bin:" << endl;
    cerr << "=================================================================" << endl;

    double FiducialVolume_Total = 0;
    double FiducialMass_Total = 0;
    double Exposure_Total = 0;

    double FiducialVolume_Internal = 0;
    double FiducialMass_Internal = 0;
    double Exposure_Internal = 0;

    double FiducialVolume_External = 0;
    double FiducialMass_External = 0;
    double Exposure_External = 0;

    if (!(_R_Lower == 0 && _R_Upper == 0)) {
        cerr << "not ready ..." << endl;
        return false;

        //     // single Radius bin
        //     _NofRadiusBin = 0;

        //     double radius_bin_min = _R_Lower;
        //     double radius_bin_max = _R_Upper;
        //     int radius_bin_internal = 1;

        //     _RadiusBinMin[_NofRadiusBin] = radius_bin_min;
        //     _RadiusBinMax[_NofRadiusBin] = radius_bin_max;
        //     _RadiusBinInternal[_NofRadiusBin] = radius_bin_internal;

        //     _FiducialVolumeBin[_NofRadiusBin] = getFiducialVolume(radius_bin_min, radius_bin_max, _Z_Lower, _Z_Upper); // m3;
        //     _FiducialMassBin[_NofRadiusBin] = _FiducialVolumeBin[_NofRadiusBin] * 0.78013 * 1.0e-3; // kton
        //     _ExposureBin[_NofRadiusBin] = _FiducialMassBin[_NofRadiusBin] * _LiveTime; // kton-day

        //     _NofRadiusBin++;
    } else {
        // multi Radius bin
        ifstream fRadiusBin(RadiusBin_name.c_str());
        if (!fRadiusBin) {
            cerr << "ERROR : cannot open input file" << endl;
            return false;
        }

        _NofRadiusBin = 0;
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

            if (_NofRadiusBin >= 100) {
                cerr << "ERROR : RadiusBin format error (out of range)" << endl;
                return false;
            }

            _RadiusBinMin[_NofRadiusBin] = radius_bin_min;
            _RadiusBinMax[_NofRadiusBin] = radius_bin_max;
            _RadiusBinInternal[_NofRadiusBin] = radius_bin_internal;

            cerr << _NofRadiusBin << " :" << endl;
            cerr << "R (m) = " << _RadiusBinMin[_NofRadiusBin] << " - " << _RadiusBinMax[_NofRadiusBin] << endl;
            cerr << endl;

            _NofRadiusBin++;
        }

        fRadiusBin.close();
        cerr << "=================================================================" << endl;

        // multi Theta bin
        ifstream fThetaBin(ThetaBin_name.c_str());
        if (!fThetaBin) {
            cerr << "ERROR : cannot open input file" << endl;
            return false;
        }

        _NofThetaBin = 0;
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

            if (_NofThetaBin >= 10) {
                cerr << "ERROR : ThetaBin format error (out of range)" << endl;
                return false;
            }

            _ThetaBinMin[_NofThetaBin] = theta_bin_min;
            _ThetaBinMax[_NofThetaBin] = theta_bin_max;
            _ThetaBinInternal[_NofThetaBin] = theta_bin_internal;

            cerr << _NofThetaBin << " :" << endl;
            cerr << "CosTheta = " << _ThetaBinMin[_NofThetaBin] << " - " << _ThetaBinMax[_NofThetaBin] << endl;
            cerr << endl;

            _NofThetaBin++;
        }

        fThetaBin.close();
        cerr << "=================================================================" << endl;

        //This variable counts the total number of position bin
        _NofRadiusThetaBin = 0;

        //Loop through all the position bin, calculating the Fiducal Volume/Mass/Exposure per bin.
        for (int n1 = 0; n1 < _NofRadiusBin; n1++) {
            for (int n2 = 0; n2 < _NofThetaBin; n2++) {
                double radius_bin_min, radius_bin_max;
                int radius_bin_internal;

                radius_bin_min = _RadiusBinMin[n1];
                radius_bin_max = _RadiusBinMax[n1];
                radius_bin_internal = _RadiusBinInternal[n1];

                double theta_bin_min, theta_bin_max;
                int theta_bin_internal;

                theta_bin_min = _ThetaBinMin[n2];
                theta_bin_max = _ThetaBinMax[n2];
                theta_bin_internal = _ThetaBinInternal[n2];

                _FiducialVolumeBin[_NofRadiusThetaBin] = getFiducialVolume(radius_bin_min, radius_bin_max, theta_bin_min, theta_bin_max); // m3;
                _FiducialMassBin[_NofRadiusThetaBin] = _FiducialVolumeBin[_NofRadiusThetaBin] * 0.78013 * 1.0e-3; // kton
                _ExposureBin[_NofRadiusThetaBin] = _FiducialMassBin[_NofRadiusThetaBin] * _LiveTime; // kton-day

                cerr << _NofRadiusThetaBin << " :" << endl;
                cerr << "R (m) = " << _RadiusBinMin[n1] << " - " << _RadiusBinMax[n1] << endl;
                cerr << "CosTheta = " << _ThetaBinMin[n2] << " - " << _ThetaBinMax[n2] << endl;
                cerr << "FiducialVolumeBin (m3) = " << _FiducialVolumeBin[_NofRadiusThetaBin] << endl;
                cerr << "FiducialMassBin (kton) = " << _FiducialMassBin[_NofRadiusThetaBin] << endl;
                cerr << "ExposureBin (kton-day) = " << _ExposureBin[_NofRadiusThetaBin] << endl;
                cerr << endl;

                _NofRadiusThetaBin++;
            }
        }
    }
    cerr << "=================================================================" << endl;

    for (int n = 0; n < _NofRadiusThetaBin; n++) {
        FiducialVolume_Total += _FiducialVolumeBin[n];
        FiducialMass_Total += _FiducialMassBin[n];
        Exposure_Total += _ExposureBin[n];

        if (_RadiusBinInternal[n] == 1) {
            FiducialVolume_Internal += _FiducialVolumeBin[n];
            FiducialMass_Internal += _FiducialMassBin[n];
            Exposure_Internal += _ExposureBin[n];
        } else {
            FiducialVolume_External += _FiducialVolumeBin[n];
            FiducialMass_External += _FiducialMassBin[n];
            Exposure_External += _ExposureBin[n];
        }
    }

    fiducialVolume = new Parameter("fv", "data_parameter", true);
    fiducialVolume->SetIValue("total_bin", _NofRadiusThetaBin);
    fiducialVolume->SetIValue("radius_bin", _NofRadiusBin);
    fiducialVolume->SetIValue("theta_bin", _NofThetaBin);


    cerr << "FiducialVolume_Total (m3) = " << FiducialVolume_Total << endl;
    cerr << "FiducialMass_Total (kton) = " << FiducialMass_Total << endl;
    cerr << "Exposure_Total (kton-day) = " << Exposure_Total << endl;
    cerr << endl;

    cerr << "FiducialVolume_Internal (m3) = " << FiducialVolume_Internal << endl;
    cerr << "FiducialMass_Internal (kton) = " << FiducialMass_Internal << endl;
    cerr << "FiducialVolume_Internal (m3) = " << FiducialVolume_Internal << endl;
    cerr << endl;

    cerr << "FiducialVolume_External (m3) = " << FiducialVolume_External << endl;
    cerr << "FiducialMass_External (kton) = " << FiducialMass_External << endl;
    cerr << "Exposure_External (kton-day) = " << Exposure_External << endl;
    cerr << endl;
}

double DataFactory::GetHistogramInterpolation(Parameter* infoParam, int j, int n) {
    // cout<<infoParam->GetS("name")<<endl;
    // std::vector<string> basename = infoParam->GetSArray("basename");
    // std::vector<bool> TaggingEfficiency = infoParam->GetBArray("basename");
    // std::vector<int> thetaVariance = infoParam->GetIArray("basename");
    int n1 = n / _NofThetaBin;
    int n2 = n % _NofThetaBin;
    double y = 0.0;
    for (int base_index=0; base_index < infoParam->GetNFitParameter(); base_index++) {
        string basename = infoParam->GetBasename(base_index);
        int theta = infoParam->GetTheta(base_index);
        bool tagging = infoParam->GetTagging(base_index);
        double expectation = 0.0;
        double binContent = 0.0;
        if (theta != -1) {
            char title[256];
            sprintf(title, "h_%s_%d",basename.c_str(), theta);
            //binContent = gMCHistogram[title]->ProjectionX("", n1+1, n1+1, n2+1, n2+1)->Interpolate(x);
            binContent = gMCHistogram[title]->GetBinContent(j+1, n1+1, n2+1);
        } else {
            //binContent = gMCHistogram["h_" + basename]->ProjectionX("", n1+1, n1+1, n2+1, n2+1)->Interpolate(x);
            binContent = gMCHistogram["h_" + basename]->GetBinContent(j+1, n1+1, n2+1);
        }
        double normalization = gFitConstant["GenerateEvent_" + basename];
        if (normalization == 0.0) {
            //cerr<<basename<<" is not simulated! Skip this background component."<<endl;
            continue;
        }
        if (basename.find("XeLS") != std::string::npos) {
            // normalization = normalization / Exposure_XeLS_All * _ExposureBin[n];
            normalization = normalization / Exposure_XeLS_150cm * _ExposureBin[n];
        } else if (basename.find("KamLS") != std::string::npos) {
            // normalization = normalization / Exposure_KamLS_All * _ExposureBin[n];
            normalization = normalization / Exposure_KamLS_norm * _ExposureBin[n];
        }
        expectation = gFitConstant["Ratio_" + basename] * binContent / normalization / _bin_width;
        //cout<<Exposure_XeLS_150cm<<"     "<<_ExposureBin[n]<<"   "<<binContent<<"   "<<_bin_width<<endl;
        if (tagging) {
            expectation *= gFitConstant["TaggingEfficiency_" + basename];
        }
        y += expectation;
    }
    return y;
}
