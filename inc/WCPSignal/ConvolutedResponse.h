#ifndef WIRECELLSIGNAL_CONVOLUTEDRESPONSE_H
#define WIRECELLSIGNAL_CONVOLUTEDRESPONSE_H

#include <string>
#include <stdio.h>
#include <iostream>
#include "TString.h"
#include "TH1.h"
#include "TH1F.h"
#include "TFile.h"

#include "WCPSignal/ElectronicsConfig.h"
#include "WCPSignal/GenElecRsp.h"
#include "WCPSignal/FieldRsp.h"

namespace WCPSignal {

  class ConvolutedResponse {

  public:
    ConvolutedResponse(const std::string infilename = "convoluted_response.root");
    ConvolutedResponse(ElectronicsConfig *conf, const std::string garfield_file = "./signal/dune.root", const float overall_time_offset = -50, const float uv_time_offset = 0, const float uw_time_offset = 0, const std::string outfilename = "convoluted_response.root");
    ConvolutedResponse(FieldRsp *fR, GenElecRsp *eR, const std::string outfilename = "convoluted_response.root");
    ~ConvolutedResponse();
    
    TH1F* GetResponseFunction(int plane, int index1, int index2, int Tshift = 0, double charge = 1);
    //TH1* GetResponseMag(int plane, int index1, int index2);
    //TH1* GetResponsePh(int plane, int index1, int index2);
    
    void BuildConvolutedResponse();
    void OutputConvolutedResponse();    
  private:
    std::string outfile;
    TFile *infile;
    TH1 *hRsp[3][6][7];
    TH1 *hRspMag[3][6][7];
    TH1 *hRspPh[3][6][7];
    TH1F *h;
    GenElecRsp *eRsp;
    FieldRsp *fRsp;
  };
  
}

#endif
