#ifndef WIRECELLSIGNAL_ELECTRONICSCONFIG_H
#define WIRECELLSIGNAL_ELECTRONICSCONFIG_H

#include <iostream>

namespace WireCellSignal {

  struct ElectronicsConfig {

    double gain; // gain in unit of mV/fC
    double shpTime; // shaping time in unit of micro second
    int nTdc; // number of TDC per frame
    double smplTime; // length per frame in unit of micro second    
    double adc; // ADC/mV
    double digitFreq; // digitization frequency in MHz
    
    ElectronicsConfig();
    ~ElectronicsConfig();    

    void SetGain(double g) { gain = g; }
    double Gain() { return gain; }  
    void SetShapingTime(double s) { shpTime = s; }
    double ShapingTime() { return shpTime; }
    void SetNTDC(int t) { nTdc = t; }
    int NTDC() { return nTdc; }
    void SetSamplingTime(double s) { smplTime = s; }
    double SamplingTime() { return smplTime; }
    void SetADC(double a) { adc = a; }
    double ADC() { return adc; }    
    void SetDigitFreq(double f) { digitFreq = f; }
    double DigitFreq() { return digitFreq; }
  };
  
}

#endif
