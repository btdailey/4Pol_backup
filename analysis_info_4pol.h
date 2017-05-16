// making a struct to pass to tree in MyCorrelator.cxx 
// analysis info in arrays of length 4 
// for things that change with pol 

struct analysis_info_4pol { 
  int eventNumber;
  double anitaLat;
  double anitaLon;
  double anitaAlt;
  double heading;

  double deltaTheta[4]; //changes with map, needs 4 pols 
  double deltaPhi[4]; // 4
  double deltamcmTheta[4]; // 4
  double deltamcmPhi[4]; // 4
  double mapSNR[4]; // 4
  double peakVal[4]; // 4
  double ratioFirstToSecondPeak[4]; //4
  double snrCoherent[4]; // 4, changes with antennas you use which is based on map
  double snrPeakAnt[4]; // set really early, only based on vpol right now, but could go in getgraphsthisevent and change it
  double maxSignalPeak[4]; //only based on vpol right now, probably need for all 4 


  double peakHilbertCoherent[4]; // 4

  double snrPeakAfterFilter[4]; // 4
  double snrPeakAfterFilter2[4]; // 4
  int didIFilter[4]; // 4
  int triggerOrPhiMaskFlag[4]; // 4

  double thetaMap[4];
  double phiMap[4];


  double secondThetaMap[4]; // 4 
  double secondPhiMap[4]; // 4 
  double headingOfThisEvent[4]; // 4, outputs heading - peakphi of map, so changes with map 

  double tertiaryThetaMap[4]; // 4,  CHANGED NAME SO ALSO CHANGE NAME IN POINTTHISEVENT
  double tertiaryPhiMap[4]; // 4,    CHANGED NAME SO ALSO CHANGE NAME IN POINTTHISEVENT 
  int varnerFlag[4];  // 4, SEE IF CAN THROW OUT    
  int varnerFlag2[4]; // 4
 
  int phiMaskFlag[4]; // 4, dependent on peak phi 
  int hwTriggerFlag[4]; // 4
 
  double polAngleCoherent[4]; // 4 
  double polFractionCoherent[4]; // 4
  int didIFilterAboveSatellite[4]; // 4
  //int didIFilterHoriz[4]; // do not need 
  //int didIFilterAboveSatelliteHoriz[4]; // do not need this anymore 

  double meanFreqVert[4]; // 4, do we need all 4? 
  //double meanFreqHoriz[4]; // do not need this anymore 


  double peakThetaFinal[4]; // not needed 
  double peakPhiFinal[4]; // not needed 
  //double finaltheta; // not needed 

  int eventPointedFlag[4]; // 4 
  int eventTracedFlag[4]; // 4
  double sourceLat[4]; // 4 
  double sourceLon[4]; // 4
  double sourceAlt[4]; // 4


  double rmsNoiseCoherent[4]; // 4 
  double noiseBeforeFilter[4]; // 4  //we need to code this more
  double CWheight[4]; // 4
  double SNR_ant[4]; // 4, add in index 
  double SNR_ant_triggered[4]; // 4, add in index 
  double SNR_ant_closest[4]; // 4, add in index 
  double SNR_ant_coherent[4]; // 4, add in index 
  double PowerCut[4]; // 4, add in index, might not keep as array of 40 antennas
  int  CoherentAnts[4]; // 4, list of antennas 
  double distance_from_source[4]; // 4 
  double peak2peak_signal[4]; // 4, depend on after filtering which will change for each pol 

  double peakVoltage_2[4];
  double peakSNR_2[4];

  double peakVal_box[4];//peakVal in box surrounding peak in Vpol

  double max_correlation[4][40][40];
  double RMS_ants[4][40];
  double time_delay[4][40][40];
  int num_bins_filtered;
  /////flags////////
  int mainrfcmflag;   
  int bigenoughpeakflag;
  int dcoffsetflag;
  int shorttraceflag;
  
  int nadirrfcmflag;
  int payloadblastflag;
	  
  double SNR_coherent[4];
  double SNR_coherent2[4];//for powerSNR min value==1

  int SNR_noise_bin[4];
  int noiseFlag[4];

  double preFilter_power;
  double postFilter_power;
};  
