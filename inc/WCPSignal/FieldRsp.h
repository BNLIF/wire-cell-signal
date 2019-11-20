#ifndef WIRECELLSIGNAL_FIELDRSP_H
#define WIRECELLSIGNAL_FIELDRSP_H

#include <string>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TVirtualFFT.h"
#include "TFile.h"

#include "WCPSignal/ElectronicsConfig.h"

namespace WCPSignal {

  class FieldRsp {

  public:
    FieldRsp(ElectronicsConfig *con, const float overall_time_offset = 0, const float uv_time_offset = 0, const float uw_time_offset = 0);
    FieldRsp(ElectronicsConfig *con, const std::string infilename, TFile *f = 0, const float overall_time_offset = -50/*us*/, const float uv_time_offset = 0, const float uw_time_offset = 0); 
    ~FieldRsp();

    TH1F* GetFieldResponseHist(int i) { return fieldResp[i]; }
    TH1* GetFieldResponseInFreq(int i) { return hFieldRespInFreq[i]; }
    TH1* GetFieldResponseInPhase(int i) { return hFieldRespInPhase[i]; }
    TH1F* GetFieldResponseHist(int i, int j, int k) { return hFieldResp2D[i][j][k]; }
    TH1* GetFieldResponseInFreq(int i, int j, int k) { return hFieldRespInFreq2D[i][j][k]; }
    TH1* GetFieldResponseInPhase(int i, int j, int k) { return hFieldRespInPhase2D[i][j][k]; }

    void DrawFieldResponse(const std::string canvas_name = "field_response");
    void DrawFieldResponseFFT(const std::string canvas_name = "field_response_fft");
    void SetOutfile(TFile *o) { outfile = o; }    

  private:
    void InitFieldResponse(const float overall_time_offset, const float uv_time_offset, const float uw_time_offset);
    void InitFieldResponse(const std::string infilename, const float overall_time_offset, const float uv_time_offset, const float uw_time_offset);
    void FFTFieldResponse();

    TFile *outfile;
    ElectronicsConfig *config;
    // single-wire response functions
    TH1F *fieldResp[3]; // [0]=U, [1]=V, [2]=Y (or W)
    TGraph *gFieldResp[3];
    TH1 *hFieldRespInFreq[3];
    TH1 *hFieldRespInPhase[3];

    // adjacent wires included response functions
    TGraph *gFieldResp2D[3][6][7];
    TH1F *hFieldResp2D[3][6][7];
    TH1 *hFieldRespInFreq2D[3][6][7];
    TH1 *hFieldRespInPhase2D[3][6][7];

    bool simu2D;
  };
  
}

#endif
