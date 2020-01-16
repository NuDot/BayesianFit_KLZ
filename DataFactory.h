#ifndef _BFIT_DATAFACTORY
#define _BFIT_DATAFACTORY

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <cmath>
#include <ParameterFactory.h>
#include "TH3.h"
#include "TFile.h"

#include "TSystem.h"
#include "TApplication.h"

#include "TCanvas.h"

#include "TNtuple.h"

#include "TH1.h"

#include "TH2.h"

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
class DataFactory
{
    
public:
    
    // Constructor
    DataFactory(ParameterFactory* Param);
    
    // Destructor
    ~DataFactory();
    
    bool DataReader();
    
    std::map<std::string,TH3D* > MCReader();
    
    double getFiducialVolume(double R_Lower, double R_Upper, double Theta_Lower, double Theta_Upper);
    
    bool BuildFiducialVolume();
    
    bool CheckUniformity();
    
    Parameter* GetFVParameter() { return fiducialVolume; };
    
    int GetRadiusBinInternal(int iint) { return _RadiusBinInternal[iint]; };
    
    double GetExposureBin(int ibin) { return _ExposureBin[ibin]; };
    
    TH3D* GetEmptyHistogram(const std::string& name) {
        TH3D* output = new TH3D(name.c_str(), "", nbinx, xbins, nbiny, ybins, nbinz, zbins);
        for (int i=0;i<=nbinx;i++) {
            for (int j=0;j<=nbiny;j++) {
                for (int k=0;k<=nbinz;k++) {
                    output->SetBinContent(i,j,k,0.0);
                }
            }
        }
        return output;
    }
    
    double GetHistogramInterpolation(Parameter* infoParam, int j, int n);
    
    double** GetData(const std::string& name) {
        if (name == "Emean") {
            return _E_Mean;
        } else if (name == "data") {
            return _Data;
        } else if (name == "error") {
            return _Data_e;
        } else {
            cerr << "Wrong input name" << name << "given!" << endl;
            abort();
        }
    }
    
    double GetRBinMin(int bin) {return ybins[bin];}
    double GetRBinMax(int bin) {return ybins[bin+1];}
    
    double GetThetaBinMin(int bin) {return zbins[bin];}
    double GetThetaBinMax(int bin) {return zbins[bin+1];}
    
    
    
    
private:
    ParameterFactory* fParameters;
    
    double _PromptThreshold = 0.40; //???
    double _PrescaleFactor = 0.01024;
    
    //Energy Binning
    const int _NofEnergyBin = 90;
    double _EnergyMin = 0.5;
    double _EnergyMax = 5.0;
    double _EnergyBinWidth = (_EnergyMax - _EnergyMin) / double(_NofEnergyBin);
    
    static const int _nbin = 400;
    double _min = 0;
    double _max = 20.0;
    double _bin_width = (_max - _min) / double(_nbin);
    
    string Root_name;
    string SimulationSpectrumFile_XeLS;
    string SimulationSpectrumFile_KamLS;
    string SimulationSpectrumFile_film;
    string RadiusBin_name;
    string ThetaBin_name;
    
    int _StartRun;
    int _LastRun;
    double _LiveTime;
    // double _NumberOfBiPo214;
    double _R_Lower;
    double _R_Upper;
    double _Z_Lower;
    double _Z_Upper;
    // double _DetectionEfficiency;
    // double _DeltaNormalize; //??
    // double _DeltaEnergy; //?
    double Exposure_XeLS_All;
    double Exposure_KamLS_All;
    double Exposure_KamLS_norm;
    double Exposure_XeLS_150cm;
    
    //dummy variables:
    double x[100000], y[100000];
    int N;
    
    char buffer[256];
    double val, val_e;
    double dummy, obs_energy;
    string buf1, buf2, buf3, buf4, buf5, buf6, buf7;
    
    //1000 energy bin, 100 radius bin, 10 theta bin
    double xbins[1000], ybins[100], zbins[10];
    int nbinx = 0;
    int nbiny = 0;
    int nbinz = 0;
    double x1bins[1000];
    int nbinx1 = 0;
    
    //Binning information on Radius
    int _NofRadiusBin = 0;
    double _RadiusBinMin[100];
    double _RadiusBinMax[100];
    int _RadiusBinInternal[100];
    
    //Binning information on Theta
    int _NofThetaBin = 0;
    double _ThetaBinMin[10];
    double _ThetaBinMax[10];
    int _ThetaBinInternal[10];
    
    //Fiducial Mass/Exposure
    int _NofRadiusThetaBin = 0;
    double _FiducialVolumeBin[1000];
    double _FiducialMassBin[1000];
    double _ExposureBin[1000];
    
    //Overall data storage
    double * _E_Mean[1000];
    double * _Data[1000];
    double * _Data_e[1000];
    int _NofData;
    
    std::map<std::string,TH3D* > gMCHistogram;
    std::map<std::string, double> gFitConstant;
    std::map<std::string, int> gGenerateEvent;
    
    Parameter* fiducialVolume;
    
    TFile *rootf_mc_XeLS;
    TFile *rootf_mc_KamLS;
    TFile *rootf_mc_film;
    
    TCut Hspot = "sqrt((z+190.1)*(z+190.1)+x*x+y*y)>70";
};

#endif

