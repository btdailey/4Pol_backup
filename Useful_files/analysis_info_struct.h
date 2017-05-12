/////////////// OINDREE -- FOR MY REFERENCE ONLY ////////////////

double deltaTheta[NPOL_ANALYSIS]; //changes with map, needs 4 pols 
double deltaPhi[NPOL_ANALYSIS]; // 4
double deltamcmTheta[NPOL_ANALYSIS]; // 4
double deltamcmPhi[NPOL_ANALYSIS]; // 4
double mapSNR[NPOL_ANALYSIS]; // 4
double peakVal[NPOL_ANALYSIS]; // 4
double ratioFirstToSecondPeak[NPOL_ANALYSIS]; //4
double snrCoherent[NPOL_ANALYSIS]; // 4, changes with antennas you use which is based on map
double snrPeakAnt[NPOL_ANALYSIS]; // set really early, only based on vpol right now, but could go in getgraphsthisevent and change it
double maxSignalPeak[NPOL_ANALYSIS]; //only based on vpol right now, probably need for all 4 
double distanceTD[NPOL_ANALYSIS]; // only once 
double peakHilbertCoherent[NPOL_ANALYSIS]; // 4
double deltaTTD[NPOL_ANALYSIS]; // only once
double snrPeakAfterFilter[NPOL_ANALYSIS]; // 4
int didIFilter[NPOL_ANALYSIS]; // 4
int triggerOrPhiMaskFlag[NPOL_ANALYSIS]; // 4
double thetaTD[NPOL_ANALYSIS]; // only once, direction in payload coords to Taylor Dome 
double phiTD[NPOL_ANALYSIS]; // only once 
double thetaWilly[NPOL_ANALYSIS]; // only once 
double phiWilly[NPOL_ANALYSIS]; // only once 
double hwTriggerAngle[NPOL_ANALYSIS]; // only once, what angle did hardware go off? change for ANITA 3!!! 
int thisPhiMask[NPOL_ANALYSIS]; // only once 
double distanceMcM[NPOL_ANALYSIS]; // only once 
double secondThetaMap[NPOL_ANALYSIS]; // 4 
double secondPhiMap[NPOL_ANALYSIS]; // 4 
double headingOfThisEvent[NPOL_ANALYSIS]; // 4, outputs heading - peakphi of map, so changes with map 
int nadirFlag[NPOL_ANALYSIS]; // only once 
double tertiaryThetaMap[NPOL_ANALYSIS]; // 4,  CHANGED NAME SO ALSO CHANGE NAME IN POINTTHISEVENT
double tertiaryPhiMap[NPOL_ANALYSIS]; // 4,    CHANGED NAME SO ALSO CHANGE NAME IN POINTTHISEVENT 
int varnerFlag[NPOL_ANALYSIS];  // 4, SEE IF CAN THROW OUT    
int varnerFlag2[NPOL_ANALYSIS]; // 4
double pitch[NPOL_ANALYSIS]; // only once 
double roll[NPOL_ANALYSIS]; // only once 
double heading[NPOL_ANALYSIS]; // only once 
int phiMaskFlag[NPOL_ANALYSIS]; // 4, dependent on peak phi 
int hwTriggerFlag[NPOL_ANALYSIS]; // 4
int payloadBlastFlag[NPOL_ANALYSIS]; // only once 
double polAngleCoherent[NPOL_ANALYSIS]; // 4 
double polFractionCoherent[NPOL_ANALYSIS]; // 4
int didIFilterAboveSatellite[NPOL_ANALYSIS]; // 4
int didIFilterHoriz[NPOL_ANALYSIS]; // do not need this anymore  
int didIFilterAboveSatelliteHoriz[NPOL_ANALYSIS]; // do not need this anymore 
double meanFreqVert[NPOL_ANALYSIS]; // 4, do we need all 4? 
double meanFreqHoriz[NPOL_ANALYSIS]; // do not need this anymore 


double peakThetaFinal; // not needed 
double peakPhiFinal; // not needed 
double finaltheta; // not needed 
int eventNumber; // only once 
int eventPointedFlag; // 4 
int eventTracedFlag; // 4
double sourceLat=0; // 4 
double sourceLon=0; // 4
double sourceAlt=0; // 4

 
double anitaLatitude; // only once 
double anitaLongitude; // only once 
double anitaAltitude; // only once 
double anitaHeading; // only once, same as heading, can get rid of 
int realTime; //only once 
int mcmflag; // only once 
int tdflag; //once 
int shorttraceflag; //once 
int rfonlyflag; //once
int calpulserflag; //once
int channelsaturatedflag; //once 
int syncslipflag; //once 
int mainrfcmflag;//onc 
int nadirrfcmflag; //once 
int tdreflectionflag; // once 
int dcoffsetflag; //once 
int bigenoughpeakflag;  //once

int xCorPassFlag=0; // can get rid of and change input for pointhisevent 
int phase_info_flag = phase_flag; // only once 


double rmsNoiseCoherent; // 4 
double noiseBeforeFilter; // 4  //we need to code this more
double CWheight; // 4
double SNR_ant; // 4, add in index 
double SNR_ant_triggered; // 4, add in index 
double SNR_ant_closest; // 4, add in index 
double SNR_ant_coherent; // 4, add in index 
double PowerCut; // 4, add in index, might not keep as array of 40 antennas
int  CoherentAnts; // 4, list of antennas 
double distance_from_source; // 4 
double peak2peak_signal; // 4, depend on after filtering which will change for each pol 



//to do:: 
// look up how to make struct -- did this
// put in makefile -- did this
// make sure it can compile -- did this 
// set up pointer -- did this
// pointhisevent only needs 2 inputs (event number, drawmaps) 
// create tree fill it with pointer inside loopoverevents -- did this
//after tdataTracing->Fill() -- did this
// go to pointhisevent, set the stuff here .. struct-->thing = thing, might be in a separate function 



