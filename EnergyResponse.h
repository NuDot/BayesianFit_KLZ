#ifndef _BFIT_EnergyResponse
#define _BFIT_EnergyResponse

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <cmath>
#include <ParameterFactory.h>
#include <DataFactory.h>
#include "TH3.h"
#include "TFile.h"
#include "TSystem.h"
#include "TApplication.h"

#include "TCanvas.h"

#include "TGraph.h"

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

// ---------------------------------------------------------
class EnergyResponse
{

public:

    // Constructor
    EnergyResponse(ParameterFactory* Param, DataFactory* Data);

    // Destructor
    ~EnergyResponse();

    bool ReadEnergyResponse();

    double GetTaggedBi214(int n) { return _ChiSquare_Tagged_Bi214[n]; };

    double GetStackedHist(Parameter* infoParam, Parameter* EResParam, string parname, double E, int n, const std::vector<double>& params);

    Parameter* setupEnergyResponse(Parameter* infoParameter, const std::vector<double>& params);

    // double interpolate(string basename, int nbin, double x);

    // TH3D* GetEmptyHistogram(const std::string& name) {return new TH3D(name.c_str(), "", nbinx1, x1bins, nbiny, ybins, nbinz, zbins);}

private:
    ParameterFactory* fParameters;
    DataFactory* fData;

    static const int _NofFitParameters_Rate_Max = 100;
    static const int _NofFitParameters_NoRate = 9;

    double _PromptThreshold = 0.40; //???
    double _PrescaleFactor = 0.01024;

    double _LiveTime;
    // double _NumberOfBiPo214;

    //Energy Binning
    const int _NofEnergyBin = 90;
    double _EnergyMin = 0.5;
    double _EnergyMax = 5.0;
    double _EnergyBinWidth = (_EnergyMax - _EnergyMin) / double(_NofEnergyBin);

    static const int _nbin = 400;
    double _min = 0;
    double _max = 20.0;
    double _bin_width = (_max - _min) / double(_nbin);

    std::string EnergySpectrumDir;
    std::string ChiSquare_BirksCherenkov_name;

    //dummy variables:
    Parameter* buildParameter;
    double x[100000], y[100000];
    int N;

    //1000 energy bin, 100 radius bin, 10 theta bin
    double xbins[1000], ybins[100], zbins[10];
    int nbinx = 0;
    int nbiny = 0;
    int nbinz = 0;
    double x1bins[1000];
    int nbinx1 = 0;

    // //Binning information on Radius
    // int _NofRadiusBin = 0;
    // double _RadiusBinMin[100];
    // double _RadiusBinMax[100];
    // int _RadiusBinInternal[100];

    // //Binning information on Theta
    // int _NofThetaBin = 0;
    // double _ThetaBinMin[10];
    // double _ThetaBinMax[10];
    // int _ThetaBinInternal[10];

    // //Fiducial Mass/Exposure
    // int _NofRadiusThetaBin = 0;
    // double _FiducialVolumeBin[1000];
    // double _FiducialMassBin[1000];
    // double _ExposureBin[1000];

    //Overall data storage
    double * _E_Mean[1000];
    double * _Data[1000];
    double * _Data_e[1000];
    int _NofData;

    Parameter* fiducialVolume;


    static const double _Min_kB = 0.00;
    static const double _Max_kB = 0.70;
    static const int _Bin_kB = 70;

    static const double _Min_R = 0.00;
    static const double _Max_R = 0.25;
    static const int _Bin_R = 50;

    static const int NBIN_kB = _Bin_kB + 1;
    static const int NBIN_R = _Bin_R + 1;
    static const int NBIN = NBIN_kB * NBIN_R;

    double _kB_Nonlinearity[NBIN];
    double _R_Nonlinearity[NBIN];
    int _N_kB_Internal;
    int _N_R_Internal;
    double _r_kB_Internal;
    double _r_R_Internal;
    int _N_kB_External;
    int _N_R_External;
    double _r_kB_External;
    double _r_R_External;

    double _Alpha_Tagged_Bi214[NBIN];
    double _ChiSquare_Tagged_Bi214[NBIN];

    std::map<std::string, TGraph* > fEnergyResponseCurve;
    std::map<std::string, std::map<int, TGraph* > > fNonlinearityMap;
    //TGraph ** _Spectrum_Nonlinearity[100];


    int _i_ES_Xe136_0nu;
    int _i_ES_Xe136_2nu;
    int _i_ES_U238_S1;
    int _i_ES_U238_S2;
    int _i_ES_Th232_S1;
    int _i_ES_Th232_S2;
    int _i_ES_K40;
    int _i_ES_Bi210;
    int _i_ES_Po210;
    int _i_ES_Kr85;
    int _i_ES_Cs137;
    int _i_ES_Cs134;
    int _i_ES_Cs136;
    int _i_ES_Bi208;
    int _i_ES_Co60;
    int _i_ES_Y88;
    int _i_ES_Ag110;
    int _i_ES_C11;
    int _i_ES_C10;
    int _i_ES_SolarNu;
    int _i_ES_Xe137;



    



};

#endif