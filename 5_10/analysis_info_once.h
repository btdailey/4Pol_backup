// making a struct to pass to tree in MyCorrelator.cxx 
// analysis info that is needed only once not for every pol 



//double snrPeakAnt[NPOL_ANALYSIS]; // set really early, only based on vpol right now, but could go in getgraphsthisevent and change it
//double maxSignalPeak[NPOL_ANALYSIS]; //only based on vpol right now, probably need for all 4 

struct analysis_info_once { 

  double distanceTD; // only once 

  double deltaTTD; // only once


  double thetaTD; // only once, direction in payload coords to Taylor Dome 
  double phiTD; // only once 
  double thetaWilly; // only once 
  double phiWilly; // only once 
  double hwTriggerAngle; // only once, what angle did hardware go off? change for ANITA 3!!! 
  int thisPhiMask; // only once 
  double distanceMcM; // only once 

  int nadirFlag; // only once 

  double pitch; // only once 
  double roll; // only once 
  double heading; // only once 

  int payloadBlastFlag; // only once 

  //int didIFilterHoriz[NPOL_ANALYSIS]; // do not need this anymore  
  //int didIFilterAboveSatelliteHoriz[NPOL_ANALYSIS]; // do not need this anymore 

  int eventNumber; // only once 


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
  int mainrfcmflag;//once
  int nadirrfcmflag; //once 
  int tdreflectionflag; // once 
  int dcoffsetflag; //once 
  int bigenoughpeakflag;  //once


  //int phase_info_flag = phase_flag; // only once 
  int phase_info_flag; 


}; 





