#ifndef WIRECELLSIGNAL_GENNOISE_H
#define WIRECELLSIGNAL_GENNOISE_H

#include <iostream>
#include <math.h>
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TVirtualFFT.h"

#include "ElectronicsConfig.h"
#include "GenElecRsp.h"

namespace WireCellSignal {

  class GenNoise {

  public:
    GenNoise(ElectronicsConfig& ocon, double rand = 1); // a random number between 0 and 1
    ~GenNoise();
    double NoiseRMS();
    //double NoiseENC();
    void PrintNoiseInTime(TCanvas *c = 0);
    void PrintNoiseInFrequency(TCanvas *c = 0);  
    TH1 *NoiseInTime();// { return fb; }
    void SetBaseline(double b) { baseline = b; }
    double GetBaseline();
    
  private:
    void PinkNoise(double rand = 0);
    void Init(double shapingtime);
    ElectronicsConfig *config;
    GenElecRsp *eRsp;
    GenElecRsp *eRspD;
    TH1F *hERspD;
    TH1F *hERsp;
    TH1 *eRspMag;
    TH1 *eRspPh;
    TH1 *eRspMagD;
    TH1 *eRspPhD;
    double MaxPoisson;
    TF1 *MyPoisson;
    TF1 *f1;    
    double *re;
    double *im;
    TH1 *fb;
    TH1F *hFreq;
    double digitNoise;
    double baseline;
  };
  
}

#endif
