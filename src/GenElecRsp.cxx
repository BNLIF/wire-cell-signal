#include "WireCellSignal/GenElecRsp.h"

double InvTransferFunction(double *x, double *par)
{
  double f = 4.31054*exp(-2.94809*x[0]/par[1])*par[0]-2.6202*exp(-2.82833*x[0]/par[1])*cos(1.19361*x[0]/par[1])*par[0]
    -2.6202*exp(-2.82833*x[0]/par[1])*cos(1.19361*x[0]/par[1])*cos(2.38722*x[0]/par[1])*par[0]
    +0.464924*exp(-2.40318*x[0]/par[1])*cos(2.5928*x[0]/par[1])*par[0]
    +0.464924*exp(-2.40318*x[0]/par[1])*cos(2.5928*x[0]/par[1])*cos(5.18561*x[0]/par[1])*par[0]
    +0.762456*exp(-2.82833*x[0]/par[1])*sin(1.19361*x[0]/par[1])*par[0]
    -0.762456*exp(-2.82833*x[0]/par[1])*cos(2.38722*x[0]/par[1])*sin(1.19361*x[0]/par[1])*par[0]
    +0.762456*exp(-2.82833*x[0]/par[1])*cos(1.19361*x[0]/par[1])*sin(2.38722*x[0]/par[1])*par[0]
    -2.6202*exp(-2.82833*x[0]/par[1])*sin(1.19361*x[0]/par[1])*sin(2.38722*x[0]/par[1])*par[0] 
    -0.327684*exp(-2.40318*x[0]/par[1])*sin(2.5928*x[0]/par[1])*par[0] + 
    +0.327684*exp(-2.40318*x[0]/par[1])*cos(5.18561*x[0]/par[1])*sin(2.5928*x[0]/par[1])*par[0]
    -0.327684*exp(-2.40318*x[0]/par[1])*cos(2.5928*x[0]/par[1])*sin(5.18561*x[0]/par[1])*par[0]
    +0.464924*exp(-2.40318*x[0]/par[1])*sin(2.5928*x[0]/par[1])*sin(5.18561*x[0]/par[1])*par[0];
  
  if (x[0] >0&&x[0] < 10){
    return f*10.12;
  }else{
    return 0;
  }

}

WireCellSignal::GenElecRsp::GenElecRsp(ElectronicsConfig *c)
  : config(c)
{
  shapingTime = config->ShapingTime();
  gain = config->Gain();
  Init();
  FFTShapingFunction();
}

WireCellSignal::GenElecRsp::GenElecRsp(ElectronicsConfig *c, double g, double t)
  : shapingTime(t)
  , gain(g)
  , config(c)
{
  Init();
  FFTShapingFunction();
}

WireCellSignal::GenElecRsp::~GenElecRsp()
{
  if(shapingFunction) delete shapingFunction;
  if(hShapingF) delete hShapingF;
  if(hShapingInFreq) delete hShapingInFreq;
  if(hShapingInPhase) delete hShapingInPhase;
}

void WireCellSignal::GenElecRsp::Init()
{
  shapingFunction = new TF1("transFunction", InvTransferFunction, 0, 10, 2);
  shapingFunction->SetParameters(gain, shapingTime);
  hShapingF = new TH1F("shaping_function", "shaping_function", config->NTDC()*5, 0, config->NTDC()/config->DigitFreq());
  for(int i = 0; i < config->NTDC()*5; ++i) {
    double time = hShapingF->GetBinCenter(i+1)/5.;
    hShapingF->SetBinContent(i+1, shapingFunction->Eval(time));
  }
}

void WireCellSignal::GenElecRsp::SetShapingFunction(TF1 *f)
{
  shapingFunction = f;
  hShapingF->Reset();
  for(int i = 0; i < config->NTDC()*5; ++i) {
    double time = hShapingF->GetBinCenter(i+1)/5.;
    hShapingF->SetBinContent(i+1, shapingFunction->Eval(time));
  }
}

void WireCellSignal::GenElecRsp::FFTShapingFunction()
{
  TVirtualFFT::SetTransform(0);
  hShapingInFreq = hShapingF->FFT(0, "MAG");
  hShapingInPhase = hShapingF->FFT(0, "PH");
  hShapingInFreq->SetName("hShapingInFreq");
  hShapingInPhase->SetName("hShapingInPhase");
  //hShapingInFreq->Scale(1./TMath::Sqrt(config->NTDC()));
}
