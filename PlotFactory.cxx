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


void PlotFactory::plot(int rbin, int thetabin) {

    ///////////////////////////////////////////////////////////////

    int n = 0;

    //THStack* stack = new THStack("", "");

    map<string, TH1D*>::iterator it;

    //Set Data Histogram
    int nbins = hist_data->GetNbinsX();
    hist_data->SetMarkerStyle(20);

    //Set Canvas
    TCanvas* c1 = new TCanvas();
    TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.2,1.0,1.0,0);
    TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.2,0);

    //Hist Sum
    TH1D* hist_sum = new TH1D(*hist_data);
    hist_sum->SetLineColor(kBlack);
    hist_sum->SetLineStyle(1);
    hist_sum->SetLineWidth(LINE_WIDTH);
    for (int i = 1; i <= nbins; ++i) hist_sum->SetBinContent(i, 0);

    // create a container of temporary histograms
    std::vector<TH1D*> histcontainer;
    
    TLegend *legend = new TLegend(0.4,0.5,0.88,0.88);
    legend->SetBorderSize(0);
    legend->SetColumnSeparation(0.001);
    legend->SetNColumns(2);
    legend->SetTextSize(0.025);
    int color_array[10] = {2,3,4,6,7,38,33,28,17,5};
    int style_array[10] = {2,3,4,5,6,7,8,9,10,7};
    int N=0;
    TH1D* zerovbbhist;
    bool foundzerovbb = false;
    for ( it = MCMap.begin(); it != MCMap.end(); it++ )
    {
        //binValue += it->second->GetBinContent(ibin) * pars[m.GetFitParameterIndex(fittedParameters[i])];
        // if (it->first == "Other") continue;
        TH1D* temphist = it->second;
        TH1D* hist(0);
        hist = new TH1D( *(temphist) );
        hist->SetLineColor(color_array[N]);
        hist->SetLineStyle(style_array[N]);
        if (it->first != "^{136}Xe 0#nu#beta#beta") {
            if (it->first == "^{136}Xe 2#nu#beta#beta"){
                hist->SetLineStyle(1);
                hist->SetLineWidth(LINE_WIDTH);
            }
            legend->AddEntry(hist, it->first.c_str(),"l");
        } else {
            for (int ibin = 1; ibin <= nbins; ++ibin) {

                // add to expectation
                double expectation = temphist->GetBinContent(ibin);

                // set bin content
                hist->SetBinContent(ibin, expectation);
            }
            zerovbbhist = hist;
            zerovbbhist->SetLineColor(42);
            zerovbbhist->SetLineStyle(7);
            foundzerovbb = true;
            N++;
            continue;
        }

        // scale histogram
        for (int ibin = 1; ibin <= nbins; ++ibin) {

            // add to expectation
            double expectation = temphist->GetBinContent(ibin);

            // set bin content
            hist->SetBinContent(ibin, expectation);

            // add bin content
            hist_sum->SetBinContent(ibin, hist_sum->GetBinContent(ibin) + expectation);
        }

        // add histogram to container (for memory management)
        histcontainer.push_back(hist);

        // // add histogram to stack
        // stack->Add(hist);
        N++;
    }

    if (logy) pad1->SetLogy();
    pad1->Draw();
    pad1->cd();
    TH1D* hist_sum_zerovbbhist = new TH1D(*hist_data);
    hist_sum_zerovbbhist->SetLineColor(kBlack);
    hist_sum_zerovbbhist->SetLineStyle(7);
    hist_sum_zerovbbhist->SetLineWidth(LINE_WIDTH);
    histcontainer[0]->SetAxisRange(countLow, countHigh,"Y");
    histcontainer[0]->GetXaxis()->SetRangeUser(EnergyLow,EnergyHigh);
    histcontainer[0]->Draw("HIST");


    for (int i=1;i<histcontainer.size();i++) histcontainer[i]->Draw("SAME HIST");

    hist_sum->Draw("SAME HIST");
    histcontainer.push_back(hist_sum);

    if (foundzerovbb) {
        zerovbbhist->SetLineStyle(1);
        zerovbbhist->SetLineWidth(LINE_WIDTH);
        zerovbbhist->Draw("SAME HIST");
        histcontainer.push_back(zerovbbhist);
        legend->AddEntry(zerovbbhist, "0#nu#beta#beta","l");
        for (int i = 1; i <= nbins; ++i) {
            hist_sum_zerovbbhist->SetBinContent(i, hist_sum->GetBinContent(i) + zerovbbhist->GetBinContent(i));
        }
        hist_sum_zerovbbhist->Draw("SAME HIST");
        histcontainer.push_back(hist_sum_zerovbbhist);
        legend->AddEntry(hist_sum_zerovbbhist, "Total(0#nu#beta#beta 90%CL)","l");
    }
    legend->AddEntry(hist_sum, "Total(No 0#nu#beta#beta)","l");
    hist_data->Draw("SAME E1");
    histcontainer.push_back(hist_data);
    legend->AddEntry(hist_data, "Data","p");
    legend->Draw();
    gStyle->SetOptStat(0);
    c1->cd();

    static const int NBIN_ARRAY = nbins;
    Double_t x[NBIN_ARRAY];
    Double_t y[NBIN_ARRAY];
    Double_t y_data[NBIN_ARRAY];
    Double_t ex[NBIN_ARRAY];
    Double_t ey1[NBIN_ARRAY];
    Double_t ey2[NBIN_ARRAY];
    Double_t ey3[NBIN_ARRAY];

    x[0] = _EnergyMin + _EnergyBinWidth * double(0.5);
    y[0] = 1.0;
    y_data[0] = 1.0;
    ex[0] = 0.0;
    ey1[0] = 0.0;
    ey2[0] = 0.0;
    ey3[0] = 0.0;
    for (int ibin = 1; ibin <= nbins; ++ibin) {

        // add to expectation
        Double_t expect = hist_sum->GetBinContent(ibin);
        Double_t data = hist_data->GetBinContent(ibin);
        Double_t data_e = TMath::Sqrt(data);
        Double_t expect_e = TMath::Sqrt(expect);
        Double_t E_Mean_val  = _EnergyMin + _EnergyBinWidth * double(ibin+0.5);//Bin Center
        Double_t E_Lower_val = _EnergyMin + _EnergyBinWidth * double(ibin+0.0);//Bin Lower Edge
        Double_t E_Upper_val = _EnergyMin + _EnergyBinWidth * double(ibin+1.0);//Bin Upper Edge

        x[ibin] = E_Mean_val;
        y[ibin] = 1.0;
        ex[ibin] = 0.0;
        if ((expect != 0.0) && (data != 0.0))  {
            y_data[ibin] = data / expect;
            ey1[ibin] = TMath::Sqrt(TMath::Power(data_e/data, 2) + TMath::Power(expect_e/expect, 2));
            ey2[ibin] =  2.0*ey1[ibin];
            ey3[ibin] = 3.0*ey1[ibin];
        } else {
            y_data[ibin] = 1.0;
            ey1[ibin] = 0.0;
            ey2[ibin] = 0.0;
            ey3[ibin] = 0.0;
        }
        // cout<<x[ibin]<<"    "<<y[ibin]<<"    "<<y_data[ibin]<<"    "<<ey1[ibin]<<"    "<<ey2[ibin]<<"    "<<ey3[ibin]<<endl;
    }
    // //Drawing error graph with respect to sigma
    pad2->Draw();
    pad2->cd();
    TLegend *legend2 = new TLegend(0.15,0.55,0.2,0.85);
    legend2->SetBorderSize(0);
    TGraphErrors* gr3 = new TGraphErrors(NBIN_ARRAY,x,y,ex,ey3);
    gr3->GetYaxis()->SetRangeUser(-1.0,4.0);
    gr3->GetXaxis()->SetLimits(EnergyLow,EnergyHigh);
    gr3->GetYaxis()->SetLabelSize(0.1);
    gr3->GetXaxis()->SetLabelSize(0.1);
    gr3->SetTitle("");
    gr3->SetFillColor(29);
    gr3->Draw("a3");
    legend2->AddEntry(gr3, "3#sigma","f");

    TGraphErrors* gr2 = new TGraphErrors(NBIN_ARRAY,x,y,ex,ey2);
    gr2->SetFillColor(8);
    gr2->SetTitle("");
    gr2->Draw("SAME e3");
    legend2->AddEntry(gr2, "2#sigma","f");

    TGraphErrors* gr = new TGraphErrors(NBIN_ARRAY,x,y,ex,ey1);
    gr->SetFillColor(7);
    gr->SetTitle("");
    gr->Draw("SAME e3");
    legend2->AddEntry(gr, "1#sigma","f");

    TGraph* datagr = new TGraph(NBIN_ARRAY,x,y_data);
    // gr->SetFillColor(4);
    datagr->SetMarkerStyle(7);
    datagr->Draw("SAME P");
    datagr->SetTitle("");
    legend2->Draw();
    gStyle->SetOptStat(0);



    char title[256];
    sprintf(title, "%s_%d_%d.pdf", out_name.c_str(), rbin, thetabin);
    // print
    c1->Print(title);

    for (unsigned int i = 0; i < histcontainer.size(); ++i) {
        TH1D* hist = histcontainer.at(i);
        delete hist;
    }

}

void PlotFactory::SumHist(int rbin, int thetabin) {
    map<string, TH3D*>::iterator it;
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

    for ( it = MCMap_all.begin(); it != MCMap_all.end(); it++ )
    {
        string mapName = it->first;

        cout<<mapName<<endl;


        // if (MCMap.find(mapName) != MCMap.end()) {
        //     delete MCMap[mapName];
        // }
        MCMap[mapName] = new TH1D(mapName.c_str(),"",nbinx, xbins);
        for (int i = 0; i <= nbinx; ++i) MCMap[mapName]->SetBinContent(i, 0);

        for(int i=0; i<_NofEnergyBin; i++) {
            double binvalue = 0.0;
            for(int n=0; n < bin_total; n++) {
                int n1 = n / thetabin_total;
                int n2 = n % thetabin_total;
                if(((thetabin==-1)||(n2==thetabin)) && (n1<rbin)) {
                    binvalue += it->second->GetBinContent(i+1, n1+1, n2+1);
                }
            }
            MCMap[mapName]->SetBinContent(i+1, binvalue);
        }
    }

    hist_data = new TH1D("h_data","",nbinx, xbins);
    for (int i = 0; i <= nbinx; ++i) hist_data->SetBinContent(i, 0);
    for(int i=0; i<_NofEnergyBin; i++) {
        double binvalue = 0.0;
        for(int n=0; n < bin_total; n++) {
            int n1 = n / thetabin_total;
            int n2 = n % thetabin_total;
            if(((thetabin==-1)||(n2==thetabin)) && (n1<rbin)) {
                binvalue += hist_data_all->GetBinContent(i+1, n1+1, n2+1);
            }
        }
        hist_data->SetBinContent(i+1, binvalue);
    }

}


PlotFactory::PlotFactory(std::map<std::string, TH3D*> map, TH3D* data, std::vector<double> param, int rbt, int tbt) {

    MCMap_all = map;
    hist_data_all = data;
    paramVector = param;
    hist_data;
    bin_total = rbt;
    thetabin_total = tbt;


}