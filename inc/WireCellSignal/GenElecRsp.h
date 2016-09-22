#ifndef WIRECELLSIGNAL_GENELECRSP_H
#define WIRECELLSIGNAL_GENELECRSP_H

#include "TF1.h"
#include "TH1F.h"
#include "TVirtualFFT.h"
#include "TMath.h"

#include "WireCellSignal/ElectronicsConfig.h"

namespace WireCellSignal {

  class GenElecRsp {

  public:
    GenElecRsp(ElectronicsConfig*);
    GenElecRsp(ElectronicsConfig *c, double g, double t);
    ~GenElecRsp();

    void SetShapingFunction(TF1 *f);
    TF1 *GetShapingFunction() { return shapingFunction; }
    TH1F *GetShapingFunctionHist() { return hShapingF; }
    TH1 *GetShapingFInFreq() { return hShapingInFreq; }
    TH1 *GetShapingFInPhase() { return hShapingInPhase; }
    void DrawShapingFunction(const std::string canvas_name = "electronics_shaping_function.eps");
    
  private:
    void FFTShapingFunction();
    void Init();
    
    // Transfer function is given by electronics experts.
    // It is the Laplace transformation of shaping function in time domain.
    // To obtain the shaping function in time domain, inverse Laplace tranformation should be performed on the transfer function.
    // The transfer function depends on gain and shaping time.    
    ElectronicsConfig *config;
    TF1 *shapingFunction;
    TH1F *hShapingF;
    TH1 *hShapingInFreq;
    TH1 *hShapingInPhase;
    double shapingTime;
    double gain;
  };
  
}

#endif
