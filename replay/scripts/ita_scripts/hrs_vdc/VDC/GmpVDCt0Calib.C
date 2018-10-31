#include "GmpVDCt0Calib.h"
#include "TTree.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TEntryList.h"
#include "TLine.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

const TString GmpVDCt0Calib::planeName[] = {"u1", "v1", "u2", "v2"};

GmpVDCt0Calib::GmpVDCt0Calib(const char* name, const char* description) : TNamed(name, description)
{
    memset(fOffset, 0, sizeof(Double_t)*nVDCPlanes*nWiresPerPlane);
    //nWiresPerGroup = 16;
    SetGroupWireNumber(16);
    arm = "";
}

GmpVDCt0Calib::~GmpVDCt0Calib()
{
}

void GmpVDCt0Calib::Calibrate(TTree* T, Int_t hrs, TCut cut, TString plotname)
{
    if (hrs==0) arm = "R";
    if (hrs==1) arm = "L";

    if (arm.IsNull()) {
        Error("Calibrate","Branch of event type bits is not found. Quit now.");
        return;
    }
   
    cout << "Calibrating " << arm.Data() << "HRS now." << endl;
    const Double_t MIN = 1500;
    const Double_t MAX = 3000;
    const Int_t NBins =400 ;

    TEntryList* elist = new TEntryList("elist","");
    T->Draw(">>elist",cut,"entrylist");
    T->SetEntryList(elist);

    Long64_t nevents = elist->GetN();
    Info("Calibrate", "Total number of events: %lld", nevents);

    //TCanvas* c1 = 0;
    TCanvas* c1 = new TCanvas("c1","VDC t0 calibration",1000,800);
    if (!plotname.IsNull()) {
        TString file = plotname + "[";
        //c1 = new TCanvas("c1","",1000,800);
        c1->Print(file.Data());
        //delete c1;
    }
    

    UInt_t nGroupsPerPlane = nWiresPerPlane/nWiresPerGroup;
    if (nWiresPerPlane%nWiresPerGroup!=0) nGroupsPerPlane++;

    for (UInt_t i=0; i<nVDCPlanes; i++) {

        TH2D* tdcvw = new TH2D(Form("tdcvw_%s",planeName[i].Data()),Form("Raw TDC vs. wire number (%s plane)",planeName[i].Data()),nWiresPerPlane,0,nWiresPerPlane,NBins,MIN,MAX);
        tdcvw->GetXaxis()->SetTitle("Wire number");
        tdcvw->GetYaxis()->SetTitle("Raw TDC");

        T->Project(tdcvw->GetName(),Form("%s.vdc.%s.rawtime:%s.vdc.%s.wire",arm.Data(),planeName[i].Data(),arm.Data(),planeName[i].Data()));

	//cout<< "check the arm: "<< arm.Data()<<endl;
        UInt_t l=0;

        for (UInt_t j=0; j<nGroupsPerPlane; j++) {
            if (j%6==0) {
                //c1 = new TCanvas(Form("%s.vdc.%s_group%d",arm.Data(),planeName[i].Data(),j+1),"VDC t0 calibration",1000,800);
                c1->Clear();
                c1->Divide(2,3);
            }

            c1->cd((j%6)+1);
            Int_t lower = j*nWiresPerGroup;
            Int_t upper = (j+1)*nWiresPerGroup;
            //TH1F* hist = new TH1F(Form("hist_%d_%d",i,j),Form("TDC spectrum of wires %d -- %d in VDC %s plane",lower+1, TMath::Min(upper,(Int_t)nWiresPerPlane), planeName[i].Data()),NBins,MIN,MAX);
            //TString wirename = arm + ".vdc." + planeName[i] + ".wire";
            //TString cutgroup = wirename + Form(">=%d", lower) + " && " + wirename + Form("<%d", upper);
            //T->Draw(Form("%s.vdc.%s.rawtime>>%s",arm.Data(),planeName[i].Data(),hist->GetName()),cutgroup.Data());

            TH1D* hist = tdcvw->ProjectionY(Form("%s_%d",tdcvw->GetName(),j),lower+1,upper,"eo");
            hist->SetTitle(Form("Raw TDC spectrum of vdc %s plane, wire %d -- wire %d",planeName[i].Data(),lower+1,TMath::Min(upper,(Int_t)nWiresPerPlane)));
            hist->Draw("hist");
            hist->Smooth();

            Double_t t0 = Findt0(hist);
            for (UInt_t k=0; k<nWiresPerGroup && l<nWiresPerPlane; k++)
                fOffset[i][l++] = t0;
            if (t0>0) {
                TLine* line = new TLine(t0,0,t0,hist->GetMaximum());
                line->SetLineColor(kRed);
                line->SetLineWidth(2);
                line->Draw();

            }

            gPad->Update();

            if ((j%6==5 || j==nGroupsPerPlane-1) && !plotname.IsNull()) c1->Print(plotname.Data());
        }
    }

    if (!plotname.IsNull()) {
        TString file = plotname + "]";
        c1->Print(file.Data());
    }

    T->SetEntryList(0);
    delete elist;
}

Double_t GmpVDCt0Calib::Findt0(const TH1* hnew) const
{
    //TH1* hnew = h1->Clone("hnew");
    //hnew->Smooth();
    //hnew->Draw();

    if (!hnew) {
        Error("Findt0","Empty time spectrum to calibrate???");
        return 0;
    }

    Double_t dx = hnew->GetBinWidth(1);
    Double_t slope_min = 0, slope = 0;
    Double_t sdbin = -1; // Bin number with steepest descent

    Int_t nbins = hnew->GetNbinsX();
    Int_t maxbin = hnew->GetMaximumBin();
    Double_t maxcont = hnew->GetBinContent(maxbin);

    for (Int_t i=maxbin+1; i<=nbins; i++) {
        if (hnew->GetBinContent(i)>0.5*maxcont) continue;
        if (i>1 && i<nbins) {
            slope = (hnew->GetBinContent(i+1)-hnew->GetBinContent(i-1))/(2.*dx);
        } else if (i==1) {
            slope = (hnew->GetBinContent(i+1)-hnew->GetBinContent(i))/dx;
        } else if (i==nbins) {
            slope = (hnew->GetBinContent(i)-hnew->GetBinContent(i-1))/dx;
        }

        if (TMath::Abs(slope)<1.e-4) slope = 0.;

        if (slope<slope_min) {
            slope_min = slope;
            sdbin = i;
        }

    }

    Double_t t0 = 0.;
    if (slope_min<-1 && hnew->GetEntries()>1000.) {
        t0 = hnew->GetBinCenter(sdbin)-hnew->GetBinContent(sdbin)/slope_min;
    }

    return t0;
}

void GmpVDCt0Calib::SaveDatabase(const char* database) const
{
    std::ofstream* outfile = new std::ofstream;
    outfile->open(database);

    if (!outfile->is_open()) {
        Error("SaveDatabase","Can not open the specified database file: %s", database);
        outfile->close();
        delete outfile;
        return;
    }

    *outfile << fixed << setprecision(1);
    for (UInt_t i=0; i<nVDCPlanes; i++) {
      //L.vdc.u1.tdc.offsets =
        *outfile<<arm<<".vdc."<<planeName[i]<<".tdc.offsets="<<endl;
      //   *outfile << "[ " << arm << ".vdc." << planeName[i] << " ]" << endl;
       //  *outfile << "TDC Offset Values" << endl; 
        for (UInt_t j=0; j<nWiresPerPlane; j++) {
           
            *outfile  << setw(6) << fOffset[i][j] << " ";
        //   *outfile << setw(3) << j+1 << " " << setw(6) << fOffset[i][j] << " ";
            if (j%8==7) *outfile << endl;
        }
        *outfile << endl;
    }
    outfile->close();
    delete outfile;

    Info("SaveDatabase","Calibration coefficients saved to %s",database);
}

ClassImp(GmpVDCt0Calib)
