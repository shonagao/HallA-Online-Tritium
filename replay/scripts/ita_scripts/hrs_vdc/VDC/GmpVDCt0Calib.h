#ifndef GMPVDCT0CALIB_H_
#define GMPVDCT0CALIB_H_

#include "TNamed.h"
#include "TCut.h"
#include "TString.h"

class TTree;
class TH1;

class GmpVDCt0Calib : public TNamed
{
    public:
        GmpVDCt0Calib(const char* name = "GmpVDCt0Calib", const char* description = "Calss for VDC t0 calibration");
        virtual ~GmpVDCt0Calib();

        void Calibrate(TTree*, Int_t hrs, TCut = "1", TString plotname = "");

        void SaveDatabase(const char*) const;
        void SetGroupWireNumber(UInt_t nwire) { nWiresPerGroup = nwire; }

        Double_t GetOffset(UInt_t plane, UInt_t wire) const { assert(plane<nVDCPlanes && wire<nWiresPerPlane); return fOffset[plane][wire]; }

        static UInt_t GetNumberOfPlanes() { return nVDCPlanes; }
        static UInt_t GetNumberOfWiresPerPlane() { return nWiresPerPlane; }
        static const char* GetPlaneName(UInt_t plane) { assert(plane<nVDCPlanes); return planeName[plane].Data(); }
    private:
        static const UInt_t nVDCPlanes = 4;
        static const UInt_t nWiresPerPlane = 368;
        //static const UInt_t nGroupsPerPlane = 23;
        static const TString planeName[nVDCPlanes];
        UInt_t nWiresPerGroup;
        TString arm;

        Double_t fOffset[nVDCPlanes][nWiresPerPlane];

        Double_t Findt0(const TH1*) const;

        ClassDef(GmpVDCt0Calib,0)

};

#endif
