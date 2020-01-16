// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef _BAT_PARAMETER
#define _BAT_PARAMETER

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <cmath>
#include "TGraph.h"

using namespace std;

// ---------------------------------------------------------
class Parameter
{

public:

    // Constructor
    Parameter(const std::string& n, const std::string& t, bool fixed);

    string GetName() {return stringParamMap["name"];};
    string GetType() {return stringParamMap["type"];};
    bool isFixed() {return boolParamMap["is_fixed"];};

    bool BuildFiducialVolume();

    void SetDValue(const std::string& valType, double val) { doubleParamMap[valType] = val; };
    void SetIValue(const std::string& valType, int val) { intParamMap[valType] = val; };
    void SetBValue(const std::string& valType, bool val) { boolParamMap[valType] = val; };
    void SetSValue(const std::string& valType, string val) { stringParamMap[valType] = val; };

    void SetGraphAddress(int N, double* x, double* y, int n) { graph[n] = new TGraph(N,x,y); };
    TGraph* GetGraph(int n) { return graph[n]; };

    int GetI(const std::string& valType);
    double GetD(const std::string& valType);
    string GetS(const std::string& valType);
    bool GetB(const std::string& valType);

    bool EditParameter(string name);

    int GetNFitParameter() {return basename.size();};
    string GetBasename(int i ) { return basename[i];};
    int GetTheta(int i ) { 
        if (theta.size() > i) {
            return theta[i];
        } else {
            return -1;
        }
    };
    bool GetTagging(int i ) {
        if (tagging.size() > i) {
            return tagging[i];
        } else {
            return false;
        }
    };

    std::vector<string> GetSArray(const std::string& valType);
    std::vector<int> GetIArray(const std::string& valType);
    std::vector<bool> GetBArray(const std::string& valType);
    std::vector<double> GetDArray(const std::string& valType);

    bool readyToFit() { return fitParameterBuilt;};

    // Destructor
    ~Parameter();



private:

    std::map<std::string,bool> boolParamMap;
    std::map<std::string,int> intParamMap;
    std::map<std::string,double> doubleParamMap;
    std::map<std::string,std::string> stringParamMap;

    std::vector<string> basename;
    std::vector<int> theta;
    std::vector<bool> tagging;

    string RadiusBin_name;
    string ThetaBin_name;

    bool fitParameterBuilt = false;

    TGraph* graph[100];

    int thetaBin = 0;
    int radiusBin = 0;


};

#endif