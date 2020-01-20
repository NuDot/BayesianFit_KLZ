// ***************************************************************
// This file was created using the CreateProject.sh script
// for project BLFit.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCMTF.h>

#include "BLFitModel.h"
#include "THStack.h"
#include "PlotFactory.h"

int LINE_WIDTH=2;

//Energy Binning
int _NofEnergyBin = 90;
double _EnergyMin = 0.5;
double _EnergyMax = 5.0;
double _EnergyBinWidth = (_EnergyMax - _EnergyMin) / double(_NofEnergyBin);

 
int main()
{
    // open log file
    BCLog::OpenLog("log_plot.txt", BCLog::detail, BCLog::detail);
    BCLog::SetLogLevel(BCLog::detail);

    // create new BLFitModel object
    BLFitModel m("BLFit");

    m.PrepareToContinueMarginalization("chain.root", "BLFit_mcmc", "BLFit_parameters", false, true);

    double src[] = {2.1,109159,2413.51,1.99999,137.215,150.245,1242.57,24566.9,59474.2,1029.57,0.629182,0.8,0.722637,4.5,3.14471e-07,4.4455,0.000491212,7.21532,24.4762,15.7143,24.9963,15.8299,174.829,2.81805e-09,2650.04,2507.85,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0128184,2.9045,0.139891,0,4.5,188.924,11.3052,1366.44,2806.68,9663.37,12917.1,841.181,0.768671,3.4e+07,0.323,0.037,1.01568,0.289931,0.024997,1.01476,0.405408,0.0135369,1};

    std::vector<double> dest;

    for (int i=0; i<sizeof(src); i++) {
        dest.push_back(src[i]);
    }

    std::vector<double> paramVector = dest;
    std::vector<string> investigate_bkg;
    investigate_bkg.push_back("Rate_He6_Spallation");
    investigate_bkg.push_back("Rate_B12_Spallation");
    investigate_bkg.push_back("Rate_Li8_Spallation");


    PlotFactory * pf = new PlotFactory(m.GetMCHist(paramVector), m.GetDataHist(), paramVector, m.GetTotalBin(), m.GetThetaBin());

    for (int i=4;i<=11;i++) {
        for(int j=-1;j<=1;j++) {
            pf->PlotAll(i, j);
        }
    }

    pf->PlotAll(40, -1);
    pf->PlotROI(5, -1);

    std::vector<string> fitParamName = m.GetFittedRateParameter();
    for (int i=0;i<fitParamName.size();i++) {
        if (fitParamName[i]=="Rate_Xe136_2nu_XeLS" ||
            fitParamName[i]=="Rate_U238_S2_XeLS" ||
            fitParamName[i]=="Rate_Pileup_XeLS" ||
            fitParamName[i]=="Rate_C10_XeLS" ||
            fitParamName[i]=="Rate_Xe137_XeLS" ||
            fitParamName[i]=="Rate_U238_S2_film_0" ||
            fitParamName[i]=="Rate_U238_S2_film_1" ||
            fitParamName[i]=="Rate_He6_Spallation" ||
            fitParamName[i]=="Rate_B12_Spallation" ||
            fitParamName[i]=="Rate_Li8_Spallation"
            ) continue;
        m.GetParameter(fitParamName[i]).Fix(0.0);
    }

    // m.GetParameter("Alpha_Internal").Fix(0.0);
    // m.GetParameter("kB_Internal").Fix(0.0);
    // m.GetParameter("R_Internal").Fix(0.0);

    // for(int i=0;i<m.GetInitialPosition().size();i++) cout<<m.GetBestFitParameters()[i]<<",";

    // m.PrintCorrelationPlot(m.GetSafeName() + "_correlation.pdf");
    // m.PrintCorrelationMatrix(m.GetSafeName() + "_correlationMatrix.pdf");




    return 0;
}