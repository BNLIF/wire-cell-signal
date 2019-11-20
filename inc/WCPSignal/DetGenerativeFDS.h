#ifndef WIRECELLSIGNAL_DETGENERATIVEFDS
#define WIRECELLSIGNAL_DETGENERATIVEFDS

#include "WCPNav/FrameDataSource.h"
#include "WCPNav/DetectorGDS.h"
#include "WCPNav/WrappedGDS.h"
#include "WCPNav/SimDataSource.h"
#include "WCPNav/Depositor.h"
#include "WCPData/Units.h"
#include "WCPData/Vector.h"
#include "WCPSignal/ConvolutedResponse.h"

#include "TStopwatch.h"
#include "TH2F.h"
#include <math.h>
#include <fstream>

namespace WCPSignal {

    /** A FrameDataSource which generates frames of traces based on
     * hits from a Depositor object.
     *
     * This FDS models an infinite volume with wire planes embedded at
     * the locations given by the GeomDataSource.  Hits from provided
     * by the Depositor object are assumed to drift in the negative X
     * direction towards the wire planes at a constant
     * bin_drift_distance per digitized time bin.  A Frame of the
     * given bins_per_frame is digitized and hits are collected onto
     * wires.  Any hits outside the frame are ignored.
     */
  class DetGenerativeFDS : public WCP::FrameDataSource , public WCP::SimDataSource
    {
    public:
	/** Create a GenerativeFDS.
	 * 
	 * The dep is a Depositor.  See that class for details.
	 *
	 * The gds is a GeomDataSource that determines which wires may
	 * be "hit" by what is produced by the Depositor.
	 *
	 * The bins_per_frame determines how deep in the drift
	 * distance a single frame is.
	 *
	 * The nframes set how many frames to assume is in the data
	 * stream or negative to run forever.
	 *
	 * The bin_drift_distance is how far an electron will drift in
	 * the time collected in one time bin.  Ie, digitization
	 * period * drift speed.
	 */	
      DetGenerativeFDS(const WCP::Depositor& dep, const WCP::DetectorGDS& gds,
		       int bins_per_frame1 = 9600, int nframes_total = -1,
		       float bin_drift_distance = 0.5*1.6* units::millimeter,
		       float unit_dis=1.6,
		       float longest_drift = 3200 * units::millimeter);

	virtual ~DetGenerativeFDS();

	virtual int size() const;

	/// Explicitly set the "frame" (event) to process.  Frame number returned or -1 on error.
	virtual int jump(int frame_number); 

	virtual WCP::SimTruthSelection truth() const; 

	void clear(){frame.clear();simtruth.clear();};
	void SetResponseFunctions(ConvolutedResponse *r) { fRsp = r; }
	
	double *sideband_charge[2];
	double tot_true_charge() {return tot_charge;}
	double tot_collected_charge() {return collected_charge;}	

    private:    
	const WCP::Depositor& dep;
	const WCP::DetectorGDS& det_gds;
	TH2F *hraw[4][100][3][2][10]; // hard coded, bad. only works for no more than four cryostats and no more than 100 APAs per cryostat.
	//TH2D *hconv[4][100][3][2][10];	
	int  max_frames;
	float bin_drift_distance;
	float unit_dis;
	float max_drift;
	double tot_charge;
	double collected_charge;
	ConvolutedResponse *fRsp;
	mutable WCP::SimTruthSet simtruth;
    };

}

#endif
