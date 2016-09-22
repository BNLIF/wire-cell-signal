#include "WireCellSignal/ElectronicsConfig.h"

WireCellSignal::ElectronicsConfig::ElectronicsConfig()
{
  gain = 14; // mV/fC
  shpTime = 2; // micro second
  nTdc = 9600;
  digitFreq = 2; // MHz
  smplTime = (double)nTdc/digitFreq; // micro second
  adc = 4096/1800; // ADC/mV 4096 bits, 1800 mV
  std::cout << "Please set electronics configurations:\n"
	    << "this->SetGain(double); default is 14 mV/fC\n"
	    << "this->SetShapingTime(double); default is 2 us\n"
	    << "this->SetNTDC(int); default is 9600 TDC per frame\n"
	    << "this->SetDigitFreq(double); default is 2 MHz\n"
	    << "this->SetSamplingTime(double); default is 1.6 ms\n"
	    << "Pleaset note that samplingTime = NTDC/digitizationFrequency\n"
	    << "this->SetADC(double); default is 4096/1800 ADC/mV"
	    << std::endl;
}

WireCellSignal::ElectronicsConfig::~ElectronicsConfig()
{
}
