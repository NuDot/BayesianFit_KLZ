#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <EnergyResponse.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

EnergyResponse::EnergyResponse(ParameterFactory* Param, DataFactory* Data) {
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

    EnergySpectrumDir = ((string) getenv("DOUBLEBETA_ANALYSIS_ENERGY_SPECTRUM_DIR"));;
    ChiSquare_BirksCherenkov_name = EnergySpectrumDir + "/Chi2-from-EnergySpectrum-Bi214-EnergyNonlinearity.dat";
    fParameters = Param;
    fData = Data;
    fiducialVolume = fData->GetFVParameter();

    _LiveTime = fParameters->GetModelParameter("total_livetime")->GetD("livetime");

    // Exposure_XeLS_All = fParameters->GetModelParameter("physics_constant")->GetD("Exposure_XeLS_All");
    // Exposure_KamLS_All = fParameters->GetModelParameter("physics_constant")->GetD("Exposure_KamLS_All");

    if (!ReadEnergyResponse()) {
        cerr << "ERROR : Failed reading energy response file!" << endl;
        abort();
    }

}   

EnergyResponse::~EnergyResponse() {
    //Destructor

}

bool EnergyResponse::ReadEnergyResponse() {
    //Start of Detector Response on Energy Spectrum
    //TFile* f = TFile::Open("eres.root", "RECREATE");
    std::vector<string> orderedRateName = fParameters->GetOrderedFitList();

    //Loop through all fit parameters
    for(int i=0; i<orderedRateName.size(); i++) {
        cout<<orderedRateName[i]<<endl;

        buildParameter = new Parameter(orderedRateName[i], "energy_response", true);
        buildParameter->SetSValue("file","");

        //Reads in all the energy spectrum file
        if(orderedRateName[i]=="Rate_Xe136_0nu_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/DoubleBeta_0nu_Xe136_Sum_Visible_XeLS-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Xe136_2nu_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/DoubleBeta_2nu_Xe136_Sum_Visible_XeLS-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_U238_S1_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-Series-Uranium-U238-Ra226-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_U238_S2_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-Series-Uranium-Ra226-Pb210-Veto-BadnessCut-from-BirksCherenkov.dat");
        //if(orderedRateName[i]=="Rate_Th232_S1_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-Series-Thorium-Pb212-Pb208-Veto-BadnessCut-from-BirksCherenkov.dat"); // assuming Pb212 remain
        if(orderedRateName[i]=="Rate_Th232_S1_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-Series-Thorium-Th232-Ra224-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Th232_S2_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-Series-Thorium-Ra224-Pb208-Veto-BadnessCut-from-BirksCherenkov.dat");
        //  if(_FitParameters[i]=="Rate_Pileup_XeLS") Spectrum_name[i] = EnergySpectrumDir + "/Spectrum-Series-Thorium-Ra224-Pb208-Veto-BadnessCut-from-BirksCherenkov.dat";
        if(orderedRateName[i]=="Rate_K40_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-K40-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Bi210_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Bi210m-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_C11_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-C11p-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_C10_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-C10p-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Kr85_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Kr85m-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Cs137_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Cs137m-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Cs136_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Cs136m-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Cs134_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Cs134m-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Bi208_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Bi208p-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Bi207_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Bi207p-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Co60_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Co60m-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Y90_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Y90m-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Y88_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Y88p-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Ag110_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Ag110m-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Te129_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Te129m_0keV-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_SolarNu_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-SolarB8ES-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Xe137_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Xe137m-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Nb95_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Nb95m-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Sr89_XeLS") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Sr89m-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_Xe137_Spallation") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Xe137m-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_C11_Spallation") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-C11p-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_C10_Spallation") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-C10p-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_He6_Spallation") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-He6m-from-BirksCherenkov.dat");
        if(orderedRateName[i]=="Rate_B12_Spallation") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-B12m-from-BirksCherenkov.dat");
        //if(orderedRateName[i]=="Rate_Li8_Spallation") buildParameter->SetSValue("file", EnergySpectrumDir + "/Spectrum-BetaPlusGamma-Li8p-from-BirksCherenkov.dat");

        //Skip some of the background events
        // temporary
        if(orderedRateName[i]=="Rate_Xe134_2nu_XeLS") continue;
        if(orderedRateName[i]=="Rate_Po210_XeLS") continue;
        if(orderedRateName[i]=="Rate_Ar39_XeLS") continue;
        if(orderedRateName[i]=="Rate_C14_XeLS") continue;
        if(orderedRateName[i]=="Rate_Pileup_XeLS") continue;
        if(orderedRateName[i]=="Rate_Po210_film") continue;
        if(orderedRateName[i]=="Rate_Li8_Spallation") continue; // tentative
        if(orderedRateName[i]=="Rate_ShortLived_Spallation") continue;
        if(orderedRateName[i]=="Rate_Monochromatic") continue;
        if(orderedRateName[i]=="Rate_U238_S1_film") continue;
        if(orderedRateName[i]=="Rate_U238_S2_film") continue;
        if(orderedRateName[i].find("Rate_U238_S2_film_")!=string::npos) continue;    
        if(orderedRateName[i]=="Rate_Th232_S1_film") continue;
        if(orderedRateName[i]=="Rate_Th232_S2_film") continue;
        if(orderedRateName[i].find("Rate_Th232_S2_film_")!=string::npos) continue;
        if(orderedRateName[i]=="Rate_Pileup_film") continue;
        if(orderedRateName[i].find("Rate_Pileup_film_")!=string::npos) continue;     
        if(orderedRateName[i]=="Rate_K40_film") continue;
        if(orderedRateName[i].find("Rate_K40_film_")!=string::npos) continue;    
        if(orderedRateName[i]=="Rate_Bi210_film") continue;
        if(orderedRateName[i].find("Rate_Bi210_film_")!=string::npos) continue;    
        if(orderedRateName[i]=="Rate_Cs137_film") continue;
        if(orderedRateName[i].find("Rate_Cs137_film_")!=string::npos) continue;    
        if(orderedRateName[i]=="Rate_Cs134_film") continue;
        if(orderedRateName[i].find("Rate_Cs134_film_")!=string::npos) continue;    
        if(orderedRateName[i]=="Rate_Bi208_film") continue;
        if(orderedRateName[i]=="Rate_Co60_film") continue;
        if(orderedRateName[i]=="Rate_Y88_film") continue;
        if(orderedRateName[i]=="Rate_Ag110_film") continue;
        //  if(orderedRateName[i]=="Rate_SolarNu_XeLS") continue;
        if(orderedRateName[i]=="Rate_SolarNu_KamLS") continue;
        if(orderedRateName[i].find("Rate_Ag110_film_")!=string::npos) continue;    
        if(orderedRateName[i]=="Rate_U238_S2_KamLS") continue;
        if(orderedRateName[i]=="Rate_Th232_S2_KamLS") continue;
        if(orderedRateName[i]=="Rate_K40_KamLS") continue;
        if(orderedRateName[i].find("Rate_K40_KamLS_")!=string::npos) continue;    
        if(orderedRateName[i]=="Rate_Bi210_KamLS") continue;
        if(orderedRateName[i].find("Rate_Bi210_KamLS_")!=string::npos) continue;    
        if(orderedRateName[i]=="Rate_Kr85_KamLS") continue;
        if(orderedRateName[i]=="Rate_C10_KamLS") continue;
        if(orderedRateName[i]=="Rate_C11_KamLS") continue;

        //Assign i_ES, which is used later to set up energy response
        //The purpose of i_ES_isotope is to group isotopes together
        //e.g. we have Bi214_0 and Bi214_1 on the film with different theta bin
        // both of them will use the same energy spectrum, therefore they will have
        // the same i_ES
        if(orderedRateName[i]=="Rate_Xe136_0nu_XeLS" || 
           orderedRateName[i]=="Rate_Xe136_2nu_XeLS" || 
           orderedRateName[i]=="Rate_U238_S1_XeLS" || 
           orderedRateName[i]=="Rate_U238_S2_XeLS" || 
           orderedRateName[i]=="Rate_Th232_S1_XeLS" || 
           orderedRateName[i]=="Rate_Th232_S2_XeLS" || 
           orderedRateName[i]=="Rate_K40_XeLS" || 
           orderedRateName[i]=="Rate_Bi210_XeLS" || 
           orderedRateName[i]=="Rate_Kr85_XeLS" || 
           orderedRateName[i]=="Rate_Cs137_XeLS" || 
           orderedRateName[i]=="Rate_Cs136_XeLS" || 
           orderedRateName[i]=="Rate_Cs134_XeLS" || 
           orderedRateName[i]=="Rate_Bi208_XeLS" || 
           orderedRateName[i]=="Rate_Bi207_XeLS" || 
           orderedRateName[i]=="Rate_Co60_XeLS" || 
           orderedRateName[i]=="Rate_Y90_XeLS" || 
           orderedRateName[i]=="Rate_Y88_XeLS" || 
           orderedRateName[i]=="Rate_Ag110_XeLS" || 
           orderedRateName[i]=="Rate_C10_XeLS" || 
           orderedRateName[i]=="Rate_C11_XeLS" || 
           orderedRateName[i]=="Rate_SolarNu_XeLS" || 
           orderedRateName[i]=="Rate_Xe137_XeLS" || 
           orderedRateName[i]=="Rate_Te129_XeLS" || 
           orderedRateName[i]=="Rate_Nb95_XeLS" || 
           orderedRateName[i]=="Rate_Sr89_XeLS" || 
           orderedRateName[i]=="Rate_Xe137_Spallation" || 
           orderedRateName[i]=="Rate_C11_Spallation" || 
           orderedRateName[i]=="Rate_C10_Spallation" ||
           orderedRateName[i]=="Rate_He6_Spallation" ||
           orderedRateName[i]=="Rate_B12_Spallation"){
            if(orderedRateName[i]=="Rate_Xe136_0nu_XeLS") _i_ES_Xe136_0nu = i;
            if(orderedRateName[i]=="Rate_Xe136_2nu_XeLS") _i_ES_Xe136_2nu = i;
            if(orderedRateName[i]=="Rate_U238_S1_XeLS")   _i_ES_U238_S1 = i;
            if(orderedRateName[i]=="Rate_U238_S2_XeLS")   _i_ES_U238_S2 = i;
            if(orderedRateName[i]=="Rate_Th232_S1_XeLS")  _i_ES_Th232_S1 = i;
            if(orderedRateName[i]=="Rate_Th232_S2_XeLS")  _i_ES_Th232_S2 = i;
            if(orderedRateName[i]=="Rate_K40_XeLS")       _i_ES_K40 = i;
            if(orderedRateName[i]=="Rate_Bi210_XeLS")     _i_ES_Bi210 = i;
            if(orderedRateName[i]=="Rate_Po210_XeLS")     _i_ES_Po210 = i;
            if(orderedRateName[i]=="Rate_Kr85_XeLS")      _i_ES_Kr85 = i;
            if(orderedRateName[i]=="Rate_Cs137_XeLS")     _i_ES_Cs137 = i;
            if(orderedRateName[i]=="Rate_Cs134_XeLS")     _i_ES_Cs134 = i;
            if(orderedRateName[i]=="Rate_Bi208_XeLS")     _i_ES_Bi208 = i;
            if(orderedRateName[i]=="Rate_Co60_XeLS")      _i_ES_Co60 = i;
            if(orderedRateName[i]=="Rate_Y88_XeLS")       _i_ES_Y88 = i;
            if(orderedRateName[i]=="Rate_Ag110_XeLS")     _i_ES_Ag110 = i;
            if(orderedRateName[i]=="Rate_C10_XeLS")       _i_ES_C10 = i;
            if(orderedRateName[i]=="Rate_C11_XeLS")       _i_ES_C11 = i;
            if(orderedRateName[i]=="Rate_SolarNu_XeLS")   _i_ES_SolarNu = i;
            if(orderedRateName[i]=="Rate_Xe137_XeLS")     _i_ES_Xe137 = i;

            //Reads in spectrum parameter into io buffer
            ifstream fSpectrum(buildParameter->GetS("file").c_str(), ios::in | ios::binary);
            if(!fSpectrum) {
                cerr << "ERROR : cannot open input file --> " << orderedRateName[i] << " : " << buildParameter->GetS("file") << endl;
                //return false;
                continue;
            }

            //Note: string ChiSquare_BirksCherenkov_name = EnergySpectrumDir + "/Chi2-from-EnergySpectrum-Bi214-EnergyNonlinearity.dat");
            //Reads in Bi214 Calibration information for energy nonlinearity
            ifstream fChiSquare_BirksCherenkov(ChiSquare_BirksCherenkov_name.c_str());
            if(!fChiSquare_BirksCherenkov) {
                cerr << "ERROR : cannot open input file : " << ChiSquare_BirksCherenkov_name << endl;
                return false;
            }
            //Detector Response parameters:
            //kB: Birks constant
            //R: Cherenkov/Scintillation ratio
            //Alpha: energy scale normalization
            double kB, R, Alpha, chi2;
            double kB_0, R_0, Alpha_0, chi2_0;
            double dummy;
            double chi2_Min = 1.0e+36;
            int nbin = 0;//This nbin counts how many detector energy response curve is inside fChiSquare_BirksCherenkov

            // get chi2 from neutron
            //while(fChiSquare_BirksCherenkov >> kB >> R >> Alpha >> chi2){

            // get chi2 from Bi214
            while(fChiSquare_BirksCherenkov >> kB >> R >> dummy >> dummy >> Alpha >> chi2) {
                // set Alpha and chi2 from Bi214
                //This reads in the Energy and Spectrum
                std::ostringstream basestring;
                basestring << orderedRateName[i] << nbin;
                string basestr = basestring.str();
                _Alpha_Tagged_Bi214[nbin]     = Alpha;
                _ChiSquare_Tagged_Bi214[nbin] = chi2;
                // buildParameter->SetDValue(basestr + "Alpha_Tagged_Bi214", Alpha);
                // buildParameter->SetDValue(basestr + "ChiSquare_Tagged_Bi214", chi2);

                if(chi2_Min>chi2) chi2_Min = chi2;

                //This reads in energy spectrum information from different isotopes
                // 4 x 8 byte
                //TECHNICAL DETAIL
                fSpectrum.read((char*) &kB_0,    sizeof(double));
                fSpectrum.read((char*) &R_0,     sizeof(double));
                fSpectrum.read((char*) &Alpha_0, sizeof(double));
                fSpectrum.read((char*) &chi2_0,  sizeof(double));

                //The detector energy response parameter has to agree with the Bi214
                if(!(kB==kB_0 && R==R_0 && Alpha==Alpha_0 && chi2==chi2_0)) {
                    cerr << "ERROR : invalid kB or R or Alpha or chi2" << endl;
                    abort();
                }

                N=0;
                int xmin = 0.0;
                int xmax = 0.0;
                for(int n=0; n<5000; n++) {

                    double E;
                    double Spectrum;

                    fSpectrum.read((char*) &E,        sizeof(double));
                    fSpectrum.read((char*) &Spectrum, sizeof(double));

                    //Save a point every 10 entry
                    if(n%10==0) {
                        //This store the energy response curve.
                        //x[n] corresponds to the energy span, y n corresponds to normalized number of events.
                        x[N] = E;
                        y[N] = Spectrum;

                        N++;
                    }
                }

                //This stores the best fit spectrum
                if(kB==0.34 && R==0.000) { // best-fit from Bi214 with calibration
                    //if(kB==0.28 && R==0.015){ // best-fit from Bi214 and Tl208
                    //if(kB==0.38 && R==0.010){ // best-fit from Bi214 and Tl208
                    //if(kB==0.24 && R==0.040){ // best-fit from all sources in Kam-LS
                    fEnergyResponseCurve[orderedRateName[i]] = new TGraph(N,x,y);
                    fEnergyResponseCurve[orderedRateName[i]]->SetName(orderedRateName[i].c_str());
                    // fEnergyResponseCurve[orderedRateName[i]]->Write();
                }

                _kB_Nonlinearity[nbin] = kB;
                _R_Nonlinearity[nbin] = R;
                fNonlinearityMap[buildParameter->GetS("name")][nbin] = new TGraph(N,x,y);

                nbin++;
            }
            fSpectrum.close();
            fChiSquare_BirksCherenkov.close();

            //This loop through the previously stored Chi2 for tagged Bi214, making sure the minimum chi2
            //in this array is 0
            for(int j=0; j<nbin; j++) {
                _ChiSquare_Tagged_Bi214[j] = _ChiSquare_Tagged_Bi214[j] - chi2_Min;
            }

            continue;
        }

        ifstream fSpectrum(buildParameter->GetS("file").c_str());
        if(!fSpectrum) {
            cerr << "ERROR : cannot open input file --> " << i << " : " << buildParameter->GetS("file") << endl;
            return false;
        }
        N=0;
        while(fSpectrum >> x[N] >> y[N]) {
            N++;
        }
        fSpectrum.close();
    }
    // f->Close();
    //abort();
    return true;
}

Parameter* EnergyResponse::setupEnergyResponse(Parameter* infoParameter, const std::vector<double>& params) {

    buildParameter = new Parameter("general", "energy_response", true);
    {   
        double kB = params[infoParameter->GetI("kB_Internal")];
        double R = params[infoParameter->GetI("R_Internal")];

        if(kB<_Min_kB) kB = _Min_kB;
        if(R<_Min_R) R = _Min_R;

        //This 1.0e-10 might be machine epsilon?
        if(kB>=_Max_kB) kB = _Max_kB - 1.0e-10;
        if(R>=_Max_R) R = _Max_R - 1.0e-10;

        double Ratio_kB = (kB - _Min_kB) / (_Max_kB - _Min_kB) * _Bin_kB;
        double Ratio_R = (R - _Min_R) / (_Max_R - _Min_R) * _Bin_R;

        buildParameter->SetIValue("N_kB_Internal", int(Ratio_kB));
        buildParameter->SetIValue("N_R_Internal", int(Ratio_R));
        int _N_kB_Internal = int(Ratio_kB);
        int _N_R_Internal = int(Ratio_R);

        double _r_kB_Internal = Ratio_kB - double(_N_kB_Internal);
        double _r_R_Internal = Ratio_R - double(_N_R_Internal);
        buildParameter->SetDValue("r_kB_Internal", _r_kB_Internal);
        buildParameter->SetDValue("r_R_Internal", _r_R_Internal);

        int N_kB = _N_kB_Internal;
        int N_R = _N_R_Internal;

        double r_kB = _r_kB_Internal;
        double r_R = _r_R_Internal;
        double Alpha = params[infoParameter->GetI("Alpha_Internal")];
        buildParameter->SetDValue("Alpha_Internal", Alpha);

        int N_Nonlinearity_00 = ((N_kB + 0) * NBIN_R) + (N_R + 0);
        int N_Nonlinearity_10 = ((N_kB + 1) * NBIN_R) + (N_R + 0);
        int N_Nonlinearity_01 = ((N_kB + 0) * NBIN_R) + (N_R + 1);
        int N_Nonlinearity_11 = ((N_kB + 1) * NBIN_R) + (N_R + 1);

        double Alpha_Tagged_Bi214 = ((1.0 - r_kB) * (1.0 - r_R) * _Alpha_Tagged_Bi214[N_Nonlinearity_00]
                                     + r_kB  * (1.0 - r_R) * _Alpha_Tagged_Bi214[N_Nonlinearity_10]
                                     + (1.0 - r_kB) *        r_R  * _Alpha_Tagged_Bi214[N_Nonlinearity_01]
                                     + r_kB  *        r_R  * _Alpha_Tagged_Bi214[N_Nonlinearity_11]);

        double ChiSquare_Tagged_Bi214 = ((1.0 - r_kB) * (1.0 - r_R) * _ChiSquare_Tagged_Bi214[N_Nonlinearity_00]
                                         + r_kB  * (1.0 - r_R) * _ChiSquare_Tagged_Bi214[N_Nonlinearity_10]
                                         + (1.0 - r_kB) *        r_R  * _ChiSquare_Tagged_Bi214[N_Nonlinearity_01]
                                         + r_kB  *        r_R  * _ChiSquare_Tagged_Bi214[N_Nonlinearity_11]);
        // if ((GetTaggedBi214(N_Nonlinearity_00)>0.3*9.0) &&
        //     (GetTaggedBi214(N_Nonlinearity_01)>0.3*9.0) &&
        //     (GetTaggedBi214(N_Nonlinearity_10)>0.3*9.0) &&
        //     (GetTaggedBi214(N_Nonlinearity_11)>0.3*9.0)) {
        //     ChiSquare_Tagged_Bi214 += 1E7;
        // }
        buildParameter->SetDValue("Alpha_EnergyNonlinearity_Internal", Alpha_Tagged_Bi214);
        buildParameter->SetDValue("ChiSquare_EnergyNonlinearity_Internal", ChiSquare_Tagged_Bi214);
        // _Alpha_EnergyNonlinearity_Internal     = Alpha_Tagged_Bi214;
        // _ChiSquare_EnergyNonlinearity_Internal = ChiSquare_Tagged_Bi214;
    }

    {
        double kB = params[infoParameter->GetI("kB_External")];
        double R = params[infoParameter->GetI("R_External")];

        if(kB<_Min_kB) kB = _Min_kB;
        if(R<_Min_R)   R  = _Min_R;

        if(kB>=_Max_kB) kB = _Max_kB - 1.0e-10;
        if(R>=_Max_R)   R  = _Max_R - 1.0e-10;

        double Ratio_kB = (kB - _Min_kB) / (_Max_kB - _Min_kB) * _Bin_kB;
        double Ratio_R  = (R - _Min_R) / (_Max_R - _Min_R) * _Bin_R;

        int _N_kB_External = int(Ratio_kB);
        int _N_R_External  = int(Ratio_R);
        buildParameter->SetIValue("N_kB_External", _N_kB_External);
        buildParameter->SetIValue("N_R_External", _N_R_External);

        double _r_kB_External = Ratio_kB - double(_N_kB_External);
        double _r_R_External  = Ratio_R  - double(_N_R_External);
        buildParameter->SetDValue("r_kB_External", _r_kB_External);
        buildParameter->SetDValue("r_R_External", _r_R_External);

        int N_kB = _N_kB_External;
        int N_R  = _N_R_External;

        double r_kB  = _r_kB_External;
        double r_R   = _r_R_External;
        double Alpha = params[infoParameter->GetI("Alpha_External")];
        buildParameter->SetDValue("Alpha_External", Alpha);

        int N_Nonlinearity_00 = ((N_kB + 0) * NBIN_R) + (N_R + 0);
        int N_Nonlinearity_10 = ((N_kB + 1) * NBIN_R) + (N_R + 0);
        int N_Nonlinearity_01 = ((N_kB + 0) * NBIN_R) + (N_R + 1);
        int N_Nonlinearity_11 = ((N_kB + 1) * NBIN_R) + (N_R + 1);

        double Alpha_Tagged_Bi214 = ((1.0 - r_kB) * (1.0 - r_R) * _Alpha_Tagged_Bi214[N_Nonlinearity_00]
                                     + r_kB  * (1.0 - r_R) * _Alpha_Tagged_Bi214[N_Nonlinearity_10]
                                     + (1.0 - r_kB) *        r_R  * _Alpha_Tagged_Bi214[N_Nonlinearity_01]
                                     + r_kB  *        r_R  * _Alpha_Tagged_Bi214[N_Nonlinearity_11]);

        double ChiSquare_Tagged_Bi214 = ((1.0 - r_kB) * (1.0 - r_R) * _ChiSquare_Tagged_Bi214[N_Nonlinearity_00]
                                         +        r_kB  * (1.0 - r_R) * _ChiSquare_Tagged_Bi214[N_Nonlinearity_10]
                                         + (1.0 - r_kB) *        r_R  * _ChiSquare_Tagged_Bi214[N_Nonlinearity_01]
                                         +        r_kB  *        r_R  * _ChiSquare_Tagged_Bi214[N_Nonlinearity_11]);

        // if ((GetTaggedBi214(N_Nonlinearity_00)>0.3*9.0) &&
        //     (GetTaggedBi214(N_Nonlinearity_01)>0.3*9.0) &&
        //     (GetTaggedBi214(N_Nonlinearity_10)>0.3*9.0) &&
        //     (GetTaggedBi214(N_Nonlinearity_11)>0.3*9.0)) {
        //     ChiSquare_Tagged_Bi214 += 1E7;
        // }
        buildParameter->SetDValue("Alpha_EnergyNonlinearity_External", Alpha_Tagged_Bi214);
        buildParameter->SetDValue("ChiSquare_EnergyNonlinearity_External", ChiSquare_Tagged_Bi214);
    }
    return buildParameter;
}

double EnergyResponse::GetStackedHist(Parameter* infoParam, Parameter* EResParam, string parname, double E, int n, const std::vector<double>& params) {
    if (params[infoParam->GetI(parname)] == 0.0) return 0.0;

    E = _min + ((int((E - _min) / _bin_width) + 0.5) * _bin_width);
     // temporary
    if(infoParam->GetS("name")=="Rate_Xe134_2nu_XeLS") return 0;
    if(infoParam->GetS("name")=="Rate_Po210_XeLS") return 0;
    if(infoParam->GetS("name")=="Rate_Ar39_XeLS") return 0;
    if(infoParam->GetS("name")=="Rate_C14_XeLS") return 0;
    if(infoParam->GetS("name")=="Rate_Po210_film") return 0;
    if(infoParam->GetS("name")=="Rate_Li8_Spallation") return 0; // tentative
    if(infoParam->GetS("name")=="Rate_ShortLived_Spallation") return 0;
    if(infoParam->GetS("name")=="Rate_SolarNu") return 0;

    double Unit;

    if(infoParam->GetS("unit")=="/day/kton") {
        Unit = fData->GetExposureBin(n); // kton-day
    }
    else if(infoParam->GetS("unit")=="/day") {
        Unit = _LiveTime; // day
    }
    else {
        cerr << "ERROR : invalid unit type for "<< infoParam->GetS("name") <<":"  << infoParam->GetS("unit") << endl;
        abort();
    }

    if (infoParam->GetS("name")=="Rate_Monochromatic") {
        //double Spectrum_Monochromatic = (TMath::Gaus(E, _Mean_Monochromatic, _Sigma_Monochromatic, true) * _Rate[i] * Unit * _EnergyBinWidth);
        //cout<<E<<"Mean"<<params[infoParam->GetI("Mean_Monochromatic")]<<"Sigma"<<params[infoParam->GetI("Sigma_Monochromatic")]<<"Gauss:"<<TMath::Gaus(E, params[infoParam->GetI("Mean_Monochromatic")], params[infoParam->GetI("Sigma_Monochromatic")], true)<<"Unit:"<<Unit<<"EBin:"<<_EnergyBinWidth<<endl;
        double Spectrum_Monochromatic = (TMath::Gaus(E, params[infoParam->GetI("Mean_Monochromatic")], params[infoParam->GetI("Sigma_Monochromatic")], true) * params[infoParam->GetI("Rate_Monochromatic")] * Unit * _EnergyBinWidth);

        if(Spectrum_Monochromatic < 0) Spectrum_Monochromatic = 0;

        if(E < _PromptThreshold) Spectrum_Monochromatic = Spectrum_Monochromatic * _PrescaleFactor;
        //Directly return the monochromatic spectrum
        return Spectrum_Monochromatic;
    }
    ////////////////////////////////////////////////////////////////////////////////////
    int i_ES = -1;
    double Rate = 0;
    bool NumericalSpectrum = false;

    // spectrum with energy non-linearity
    if(infoParam->GetS("name")=="Rate_Bi207_XeLS" || 
       infoParam->GetS("name")=="Rate_Y90_XeLS" || 
       infoParam->GetS("name")=="Rate_Xe134_2nu_XeLS" || 
       infoParam->GetS("name")=="Rate_Te129_XeLS" || 
       infoParam->GetS("name")=="Rate_Nb95_XeLS" || 
       infoParam->GetS("name")=="Rate_Sr89_XeLS" || 
       infoParam->GetS("name")=="Rate_Xe137_Spallation" ||
       infoParam->GetS("name")=="Rate_C11_Spallation" || 
       infoParam->GetS("name")=="Rate_C10_Spallation" || 
       infoParam->GetS("name")=="Rate_He6_Spallation" || 
       infoParam->GetS("name")=="Rate_B12_Spallation") {
        i_ES = infoParam->GetI("i_ES");
        Rate = params[infoParam->GetI(parname)];

        NumericalSpectrum = true;
    } else {
        if(infoParam->GetS("name")=="Rate_Xe136_0nu_XeLS") i_ES = _i_ES_Xe136_0nu;
        if(infoParam->GetS("name")=="Rate_Xe136_2nu_XeLS") i_ES = _i_ES_Xe136_2nu;
        if(infoParam->GetS("name")=="Rate_U238_S1_XeLS")   i_ES = _i_ES_U238_S1;
        if(infoParam->GetS("name")=="Rate_U238_S2_XeLS")   i_ES = _i_ES_U238_S2;
        if(infoParam->GetS("name")=="Rate_Th232_S1_XeLS")  i_ES = _i_ES_Th232_S1;
        if(infoParam->GetS("name")=="Rate_Th232_S2_XeLS")  i_ES = _i_ES_Th232_S2;
        if(infoParam->GetS("name")=="Rate_Pileup_XeLS")    i_ES = _i_ES_Th232_S2;
        if(infoParam->GetS("name")=="Rate_K40_XeLS")       i_ES = _i_ES_K40;
        if(infoParam->GetS("name")=="Rate_Bi210_XeLS")     i_ES = _i_ES_Bi210;
        if(infoParam->GetS("name")=="Rate_Po210_XeLS")     i_ES = _i_ES_Po210;
        if(infoParam->GetS("name")=="Rate_Kr85_XeLS")      i_ES = _i_ES_Kr85;
        if(infoParam->GetS("name")=="Rate_Cs137_XeLS")     i_ES = _i_ES_Cs137;
        if(infoParam->GetS("name")=="Rate_Cs134_XeLS")     i_ES = _i_ES_Cs134;
        if(infoParam->GetS("name")=="Rate_Bi208_XeLS")     i_ES = _i_ES_Bi208;
        if(infoParam->GetS("name")=="Rate_Co60_XeLS")      i_ES = _i_ES_Co60;
        if(infoParam->GetS("name")=="Rate_Y88_XeLS")       i_ES = _i_ES_Y88;
        if(infoParam->GetS("name")=="Rate_Ag110_XeLS")     i_ES = _i_ES_Ag110;
        if(infoParam->GetS("name")=="Rate_C10_XeLS")       i_ES = _i_ES_C10;
        if(infoParam->GetS("name")=="Rate_C11_XeLS")       i_ES = _i_ES_C11;
        if(infoParam->GetS("name")=="Rate_SolarNu_XeLS")   i_ES = _i_ES_SolarNu;
        if(infoParam->GetS("name")=="Rate_Xe137_XeLS")     i_ES = _i_ES_Xe137;
        if(infoParam->GetS("name")=="Rate_Cs136_XeLS")     i_ES = _i_ES_Cs136;
        if(infoParam->GetS("name")=="Rate_U238_S1_film")   i_ES = _i_ES_U238_S1;
        if(infoParam->GetS("name")=="Rate_U238_S2_film")   i_ES = _i_ES_U238_S2;
        if(infoParam->GetS("name").find("Rate_U238_S2_film_")!=string::npos) i_ES = _i_ES_U238_S2;
        if(infoParam->GetS("name")=="Rate_Th232_S1_film")  i_ES = _i_ES_Th232_S1;
        if(infoParam->GetS("name")=="Rate_Th232_S2_film")  i_ES = _i_ES_Th232_S2;
        if(infoParam->GetS("name").find("Rate_Th232_S2_film_")!=string::npos) i_ES = _i_ES_Th232_S2;
        if(infoParam->GetS("name")=="Rate_Pileup_film")    i_ES = _i_ES_Th232_S2;
        if(infoParam->GetS("name").find("Rate_Th232_S2_film_")!=string::npos) i_ES = _i_ES_Th232_S2;
        if(infoParam->GetS("name").find("Rate_Pileup_film_")!=string::npos) i_ES = _i_ES_Th232_S2;
        if(infoParam->GetS("name")=="Rate_K40_film")       i_ES = _i_ES_K40;
        if(infoParam->GetS("name").find("Rate_K40_film_")!=string::npos) i_ES = _i_ES_K40;
        if(infoParam->GetS("name")=="Rate_Bi210_film")     i_ES = _i_ES_Bi210;
        if(infoParam->GetS("name").find("Rate_Bi210_film_")!=string::npos) i_ES = _i_ES_Bi210;
        if(infoParam->GetS("name")=="Rate_Po210_film")     i_ES = _i_ES_Po210;
        if(infoParam->GetS("name")=="Rate_Kr85_film")      i_ES = _i_ES_Kr85;
        if(infoParam->GetS("name")=="Rate_Cs137_film")     i_ES = _i_ES_Cs137;
        if(infoParam->GetS("name").find("Rate_Cs137_film_")!=string::npos) i_ES = _i_ES_Cs137;
        if(infoParam->GetS("name")=="Rate_Cs134_film")     i_ES = _i_ES_Cs134;
        if(infoParam->GetS("name").find("Rate_Cs134_film_")!=string::npos) i_ES = _i_ES_Cs134;
        if(infoParam->GetS("name")=="Rate_Bi208_film")     i_ES = _i_ES_Bi208;
        if(infoParam->GetS("name")=="Rate_Co60_film")      i_ES = _i_ES_Co60;
        if(infoParam->GetS("name")=="Rate_Y88_film")       i_ES = _i_ES_Y88;
        if(infoParam->GetS("name")=="Rate_Ag110_film")     i_ES = _i_ES_Ag110;
        if(infoParam->GetS("name").find("Rate_Ag110_film_")!=string::npos) i_ES = _i_ES_Ag110;
        if(infoParam->GetS("name")=="Rate_U238_S2_KamLS")  i_ES = _i_ES_U238_S2;
        if(infoParam->GetS("name")=="Rate_Th232_S2_KamLS") i_ES = _i_ES_Th232_S2;
        if(infoParam->GetS("name")=="Rate_K40_KamLS")      i_ES = _i_ES_K40;
        if(infoParam->GetS("name").find("Rate_K40_KamLS_")!=string::npos) i_ES = _i_ES_K40;
        if(infoParam->GetS("name")=="Rate_Bi210_KamLS")    i_ES = _i_ES_Bi210;
        if(infoParam->GetS("name").find("Rate_Bi210_KamLS_")!=string::npos) i_ES = _i_ES_Bi210;
        if(infoParam->GetS("name")=="Rate_Kr85_KamLS")     i_ES = _i_ES_Kr85;
        if(infoParam->GetS("name").find("Rate_Kr85_KamLS_")!=string::npos) i_ES = _i_ES_Kr85;
        if(infoParam->GetS("name")=="Rate_SolarNu_KamLS")  i_ES = _i_ES_SolarNu;
        if(infoParam->GetS("name")=="Rate_C10_KamLS")      i_ES = _i_ES_C10;
        if(infoParam->GetS("name")=="Rate_C11_KamLS")      i_ES = _i_ES_C11;
        Rate = 1.0;
    }

    if(i_ES==-1) {
        cerr << "ERROR: invalid i_ES " <<infoParam->GetS("name")<< endl;
        abort();
    }

    int N_kB, N_R;
    int N_Nonlinearity_00, N_Nonlinearity_10, N_Nonlinearity_01, N_Nonlinearity_11;
    double r_kB, r_R, Alpha;
    int n1 = n / fiducialVolume->GetI("theta_bin");

    //if _RadiusBinInternal reads 1(True), it will use the internal parameter;
    //Otherwise it use external parameter.
    if(fData->GetRadiusBinInternal(n1)==1) {
        N_kB = EResParam->GetI("N_kB_Internal");
        N_R = EResParam->GetI("N_R_Internal");

        N_Nonlinearity_00 = ((N_kB + 0) * NBIN_R) + (N_R + 0);
        N_Nonlinearity_10 = ((N_kB + 1) * NBIN_R) + (N_R + 0);
        N_Nonlinearity_01 = ((N_kB + 0) * NBIN_R) + (N_R + 1);
        N_Nonlinearity_11 = ((N_kB + 1) * NBIN_R) + (N_R + 1);

        r_kB = EResParam->GetD("r_kB_Internal");
        r_R =  EResParam->GetD("r_R_Internal");
        Alpha = EResParam->GetD("Alpha_Internal");
    } else {
        N_kB = EResParam->GetI("N_kB_External");
        N_R = EResParam->GetI("N_R_External");

        N_Nonlinearity_00 = ((N_kB + 0) * NBIN_R) + (N_R + 0);
        N_Nonlinearity_10 = ((N_kB + 1) * NBIN_R) + (N_R + 0);
        N_Nonlinearity_01 = ((N_kB + 0) * NBIN_R) + (N_R + 1);
        N_Nonlinearity_11 = ((N_kB + 1) * NBIN_R) + (N_R + 1);

        r_kB = EResParam->GetD("r_kB_External");
        r_R =  EResParam->GetD("r_R_External");
        Alpha = EResParam->GetD("Alpha_External");
    }
    /////////////////////////////////////////////////////////////////////////////
    //Spectrum Nonlinearity(probably detector energy response?)
    //What dose linint do?
    //Alpha * Rate * Unit * _EnergyBinWidth: multiply the normalization back
    string basename = fParameters->GetOrderedFitList()[i_ES];
    double Spectrum_Nonlinearity = ((1.0 - r_kB) * (1.0 - r_R) * fNonlinearityMap[basename][N_Nonlinearity_00]->Eval(E / Alpha)
                                    + r_kB * (1.0 - r_R) * fNonlinearityMap[basename][N_Nonlinearity_10]->Eval(E / Alpha)
                                    + (1.0 - r_kB) * r_R * fNonlinearityMap[basename][N_Nonlinearity_01]->Eval(E / Alpha)
                                + r_kB * r_R * fNonlinearityMap[basename][N_Nonlinearity_11]->Eval(E / Alpha)) / Alpha * Rate * Unit * _EnergyBinWidth;

    if(Spectrum_Nonlinearity < 0) Spectrum_Nonlinearity = 0;

    if(E < _PromptThreshold) Spectrum_Nonlinearity = Spectrum_Nonlinearity * _PrescaleFactor;

    if(NumericalSpectrum==true) {
        return Spectrum_Nonlinearity;
    }

    //cout<<infoParam->GetS("name")<<" "<<fParameters->GetEnergyResponse(infoParam->GetS("name"))<<endl;
    double Spectrum_Numerical = (fEnergyResponseCurve[basename]->Eval(E / Alpha) / Alpha * Rate * Unit * _EnergyBinWidth);
    double Spectrum_Nonlinearity_Numerical = Spectrum_Nonlinearity;

    // if (infoParam->GetS("name").find("XeLS")==string::npos) {
    //     cout<<infoParam->Get("name")<<"  "<<Spectrum_Nonlinearity<<endl;
    // }

    double Spectrum = 0;

    if(Spectrum_Numerical>0) {
        Spectrum =  (infoParam->GetGraph(n)->Eval(E / Alpha) / Alpha * params[infoParam->GetI(parname)] * Unit * _EnergyBinWidth) / Spectrum_Numerical * Spectrum_Nonlinearity_Numerical;
        // if (infoParam->GetS("name").find("Th232")!=string::npos) {
        //     if ((E>=2.0) && (E<=3.0)) cout<<E<<"    "<<Spectrum<<"  "<<fData->GetHistogramInterpolation(infoParam, E / Alpha, n)<<"   "<<Spectrum_Numerical<<"     "<<Spectrum_Nonlinearity_Numerical<<endl;
        // }
    }

    if(Spectrum < 0) Spectrum = 0;

    if(E < _PromptThreshold) Spectrum = Spectrum * _PrescaleFactor;

    return Spectrum;
}
