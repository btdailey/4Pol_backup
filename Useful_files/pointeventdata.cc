//tausgSystem->Reset();
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include "TTreeIndex.h"
#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TText.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include <unistd.h>
#include "TVector3.h"
#include "TRotation.h"
#include "TSpline.h"
#include "Math/InterpolationTypes.h"
#include "Math/Interpolator.h"
#include "TNtuple.h"
#include "MyCorrelator.h"
#include "analysis_info_4pol.h"

using namespace std;

class MyCorrelator;

//void Pointsimdata() {

int main() {
  TStyle* RootStyle();
  // TStyle *color=RootStyle();

  

  /*  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");
  gSystem->AddIncludePath("-I${EVENT_SIMULATION_DIR}");
  gSystem->AddIncludePath("-I${ANITA_ANALYSIS}/classes");
  gSystem->Load("libfftw3.so");
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");  
  gSystem->Load("libRootFftwWrapper.so");     	  
  gSystem->Load("libAnitaEvent.so");   	  
  gSystem->Load("libMyCorrelator.so");
  gSystem->Load("libAnitaCorrelator.so");
  gSystem->Load("libStephenAnalysisTools.so");
  */
  //int ctr_brian=0;
  int xCorPassFlag=0;
  int ratioOfPeaksPassFlag=0;
  int elevationAnglePassFlag=0;
  int peakCrossCorrFlag=0;
  int polFractionFlag=0;
  int peakHilbertFlag=0;
  int triggerflag=0;
  
  //int firstEvent;
  //int lastEvent;
  vector<int> vtimes_regular;
  vector<int> vtimes_regular1;
  vector<double> time_full;
  vector<double> time_spectrum;

  vector<int> eventNumber_vector;
  vector<int> eventNumber_tracker;

  /*double freq;
  double spec;
  double freq_cos;
  double freq_sin;
  double time;
  double time_spec;

  int times;
  int timens;
  //int eventNumber;
  int eventNumber_1;

  int startrun=12;//start of run 
  int endrun=224;//end of run 
  int basenumber=274;//base number 
  */
  // for (int run=startrun;run<=endrun;run++){
    int run;
    int eventNumber;

    //int eventNumber_array[18] = {647776,681482,1679977,2774152,6984814,6986208,9666408,11207106,12943715,15999070,16396728,17036354,17092799,17176190,17180868,17249263,17253694,26680038};

    //int run_array[18] = {13,13,18,25,61,61,89,110,127,155,157,160,160,161,161,161,161,257};


    eventNumber = 588404;//CW0
    run=12;
    
    //eventNumber = 604434;//CW1
    //run =13;

    //eventNumber =7563488;//CW2
    //run=63;

    //eventNumber = 2774152;
    //run=25;
    

    //eventNumber= 6805308;//run60
    //run=60;

    //eventNumber = 3182406;
    //run=28;

    //eventNumber =  4571978;
    //run=42;

    // eventNumber = 9986489;
    //run = 93;

    //eventNumber = 6984814;
    //run=61;

    // eventNumber = 471301;
    //run=12;

    
    //eventNumber = 681482;//one CW peak
    //run = 13;

    //eventNumber = 1679977;//one nice CW peak
    //run = 18;
    //eventNumber = 17036354;
    //run = 160;

    // eventNumber = 17180868;
    //run = 161;

    // eventNumber = 17253694;
    //run = 161;

    // eventNumber = 17176190;
    //run =161;

    //eventNumber= 17253694;
    //run = 161;
    
    //eventNumber = 687109;//TD event
    //run=13;
    //eventNumber = 731940;
    //eventNumber =  817731;//snr 20
    //run=14;
   
    // run=13;
    //eventNumber = 758742;
    //run=13;
    //eventNumber = 758986;
    //eventNumber = 759103;
    //eventNumber = 759248;
    // eventNumber = 759511;
    //run=14;
    //eventNumber = 17697490;
    //run =164;
    //eventNumber =17628458;
    //run =163;
    //SNR DISCREPANCY!
    // eventNumber= 1312785;
    // run=17;
    // eventNumber = 1389605;
    //run=17;
    //eventNumber = 1450347;
    // run=17;
    //eventNumber = 1504276;
    // run=18;
    //eventNumber = 1552864;
    //run=18;
      

    //eventNumber = 737676;//BAD RECONSTRUCTION FOR SINE SUBTRACTION
    //run=13;

    //eventNumber = 7672130;//BAD LOCATION for new reconstruction
    //run=64;
    
    //eventNumber=1193915;
    //run = 16;

    // eventNumber=8015760;
    //run = 67;
    //eventNumber = 710450;
    //run=13;
    //eventNumber = 710519;
    //eventNumber = 710696;
    //eventNumber = 711588;
    //run=13;
    //eventNumber =721107;
    //run=13;
    //eventNumber = 743421;
    //run=13;
    //eventNumber=1167901;
    
    //run=16;
    // eventNumber = 2530209;
    //eventNumber = 2531024;
    //run= 24;
    //eventNumber=17480645;
    //run=162;
    //eventNumber=2869411;
    //eventNumber=2872357;
    //run=26;
    //  int eventNumber_array[12] ={687109,708810,710450,710519,710639,710696,710745,710878, 710915, 711428,711531,711588}; 

    // run=13;

    // eventNumber = 17539931;
    //run = 163;

    //base209
    //eventNumber =3696611;
    //run=32;
    //eventNumber=3768621;
    //run = 33;
    // eventNumber=18361950;
    //run =167;


    //taylorDome bad reconstruction Abby
    //eventNumber = 2925035;
    //run = 26;
    // eventNumber = 6901374;
    //run = 60;
    //eventNumber = 7580550;
    //run = 64;
    //eventNumber = 17501145;
    //run = 163;
    //eventNumber = 17514713;
    //run = 163;
    //eventNumber = 18087843;
    //run = 166;
    ////
    //eventNumber = 6517007;
    //run = 58;


    //eventNumber = 1067605;//taylor dome event
    //run=16;
    //eventNumber = 710519;
   
    //eventNumber = 723300; //calpulser
    //eventNumber = 730902;
    //eventNumber = 711766;
    
    //run=13;
    int whichPol=0;

    //eventNumber = 768815;
    //eventNumber = 768231;
    //run = 14;
    //eventNumber = 909262;
    //run = 15;
    //eventNumber = 1622560;
    //run = 18;
    //eventNumber = 17494227;
    //run = 162;

    //eventNumber = 735221; //stripe TD event
    //eventNumber = 735326; //stripe TD event
    //eventNumber=742614;
    //run=13;
    
    //eventNumber = 17482461; //stripe TD event
    //run = 162;
    
    //eventNumber = 761555;
    //run = 14;

    //eventNumber= 17496901;
    //run = 162;
    //////////Coherent SNR low, high correlation
    
    //  eventNumber=792872;
    //eventNumber=794943;
    //run = 14;

    //eventNumber=17943984;
    //run=165;

    //eventNumber=1918483;
    //run = 20;
    

    //////////////SNR 6///////////
    //eventNumber = 739743;
    //run = 13;
    //eventNumber = 754720;
    //run = 13;
    //eventNumber = 1085169;
    //run = 16;
    /////////////SNR 30 ///////////

    //eventNumber = 1088841;
    //run = 16;
    //eventNumber = 1777256;
    //run = 19;
    //eventNumber = 2058739;
    //run = 20;

    // eventNumber = 589860;
    //run=12;

    //eventNumber = 1148994;
    //run=16;


    //HIGH COHERENT SNR//////
    //eventNumber = 21361067;
    //eventNumber = 21364081;
    //eventNumber = 21359908;
    //eventNumber = 21367635;
    
    //run = 200;
    //eventNumber = 21390221;//SNR 51
    //eventNumber= 21411486;//SNR465
    //eventNumber = 21413742;//SNR 77
    //eventNumber = 21413683;//SNR27
    //run = 201;

    //eventNumber = 21509061;//SNR40
    //eventNumber = 21553408;//SNR22
    //run=206;
    
    //eventNumber = 21574023;//SNR98
    //run = 207;

    //eventNumber = 21660369;//SNR25
    //run=208;
    //eventNumber = 21699696;
    //run = 209;
    //eventNumber = 14636821;
    //run = 143;
    //eventNumber = 16088312;
    //eventNumber = 16088304;
    //run = 155;
    //eventNumber = 16545910;
    //run = 158;
    //eventNumber = 14579114;
    //run = 142;
    //eventNumber = 8521625;
    //run = 76;
    //eventNumber = 19103714;
    //run = 174;
    //eventNumber = 1268913;
    //run = 16;

    //insertedEvents
    //eventNumber= 12943715;
    //run=127;

    //eventNumber = 25887362;
    //run = 250;

     //inserted event
    //eventNumber = 15636066;
    //run = 153;
    
    // eventNumber = 16014510;
    //run = 155;

    eventNumber = 14496361;//Hol Abby/MAtt
    run = 142;

    //eventNumber = 1457785;//Hol Abby/MAtt
    //run = 142;

    //eventNumber = 21684774;//Hol Abby/MAtt
    //run = 208;

    //eventNumber = 27146983;//Hol Abby/MAtt
    //run = 261;
    //eventNumber = 10345208;
    //run = 99;
    //eventNumber = 20564174;
    //run = 190;
    
    //10% sample/////
    //eventNumber = 8441658;
    //run = 75;
    // for (int i =0;i<5;i++){
         int i =0;
	 // eventNumber=eventNumber_array[i];
      //run = run_array[i];
    //   if(eventNumber==9666408 || eventNumber==15999070){
    //  	whichPol=1;
    //  }

      cout<<"eventNumber, run and pol are "<<eventNumber<<","<<run<<","<<whichPol<<"\n";

  //MyCorrelator *magicPtr = new MyCorrelator("/u/home/agoodhue/icemcEventMaker/runs",run, WaveCalType::kDefault);
      //MyCorrelator *magicPtr = new MyCorrelator("/data/anita/anitaplus/AnitaFlight0809RootData",run, WaveCalType::kDefault);
      MyCorrelator *magicPtr = new MyCorrelator("/data/anita/btdailey/anita_data/sample90",run, WaveCalType::kDefault);
    // MyCorrelator *magicPtr = new MyCorrelator("/home/btdailey/icemc_svn/outputs",run, WaveCalType::kDefault);
  //MyCorrelator *magicPtr = new MyCorrelator("/u/osgstorage/anita/data/mcm08",run, WaveCalType::kDefault);
  cout<<"RUN NUMBER: "<<run<<endl;
  analysis_info_4pol* pol4_Ptr = 0;

  pol4_Ptr = new analysis_info_4pol();
  //double deltaTheta, deltaPhi;
  double thetaMap,phiMap;
  
  // double sourceLon, sourceLat, sourceAlt;
  double finaltheta;
  
     magicPtr->pointThisEvent(eventNumber, 1, pol4_Ptr);
    
    //inputs are: event number, a flag that is 1 if you want to draw maps and 0 if you don't, first ntuple, 
   //                                       second ntuple, 
   //                                       third ntuple, fourth ntuple, where it pointed in theta, where it pointed in phi, 
   //                                       and polarization (0 is vertical pol, 1 h)
     cout<<"passed point this event \n\n\n";

  
     cout<<"deltaThetas are "<<pol4_Ptr->deltaTheta[0]<<" "<<pol4_Ptr->deltaTheta[1]<<" "<<pol4_Ptr->deltaTheta[2]<<" "<<pol4_Ptr->deltaTheta[3]<<"\n";
      cout<<"deltaPhis are "<<pol4_Ptr->deltaPhi[0]<<" "<<pol4_Ptr->deltaPhi[1]<<" "<<pol4_Ptr->deltaPhi[2]<<" "<<pol4_Ptr->deltaPhi[3]<<"\n";
      cout<<"mapSNR are "<<pol4_Ptr->mapSNR[0]<<" "<<pol4_Ptr->mapSNR[1]<<" "<<pol4_Ptr->mapSNR[2]<<" "<<pol4_Ptr->mapSNR[3]<<"\n";
       cout<<"mapSNR are "<<pol4_Ptr->peakVal[0]<<" "<<pol4_Ptr->peakVal[1]<<" "<<pol4_Ptr->peakVal[2]<<" "<<pol4_Ptr->peakVal[3]<<"\n";
    
       
   delete magicPtr;
   //  }//for loop
   
  
}//main
//}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TStyle* RootStyle() {

  // const char* modified = "Borrowed and adapted from paus et al";

  TStyle *RootStyle = new TStyle("Root-Style","The Perfect Style for Plots ;-)");

#ifdef __CINT__
  TStyle *GloStyle;
  GloStyle = gStyle;                          // save the global style reference

  gStyle = RootStyle;
#endif
  // otherwise you need to call TROOT::SetStyle("Root-Style")

  // Paper size

  RootStyle->SetPaperSize(TStyle::kUSLetter);

  // Canvas

  RootStyle->SetCanvasColor     (0);
  RootStyle->SetCanvasBorderSize(10);
  RootStyle->SetCanvasBorderMode(0);
  RootStyle->SetCanvasDefH      (600);
  RootStyle->SetCanvasDefW      (600);
  RootStyle->SetCanvasDefX      (10);
  RootStyle->SetCanvasDefY      (10);

  // Pads

  RootStyle->SetPadColor       (0);
  RootStyle->SetPadBorderSize  (10);
  RootStyle->SetPadBorderMode  (0);
  //  RootStyle->SetPadBottomMargin(0.13);
  RootStyle->SetPadBottomMargin(0.16);
  RootStyle->SetPadTopMargin   (0.08);
  RootStyle->SetPadLeftMargin  (0.13);
  RootStyle->SetPadRightMargin (.13);
  RootStyle->SetPadGridX       (0);
  RootStyle->SetPadGridY       (0);
  RootStyle->SetPadTickX       (1);
  RootStyle->SetPadTickY       (1);

  // Frames

  RootStyle->SetFrameFillStyle ( 0);
  RootStyle->SetFrameFillColor ( 0);
  RootStyle->SetFrameLineColor ( 1);
  RootStyle->SetFrameLineStyle ( 0);
  RootStyle->SetFrameLineWidth ( 2);
  RootStyle->SetFrameBorderSize(10);
  RootStyle->SetFrameBorderMode( 0);


  // Histograms

  RootStyle->SetHistFillColor(0);
  RootStyle->SetHistFillStyle(1);
  RootStyle->SetHistLineColor(1);
  RootStyle->SetHistLineStyle(0);
  RootStyle->SetHistLineWidth(2);

  // Functions

  RootStyle->SetFuncColor(1);
  RootStyle->SetFuncStyle(0);
  RootStyle->SetFuncWidth(1);

  //Legends 

  RootStyle->SetStatBorderSize(0);
  RootStyle->SetStatFont      (42);
  RootStyle->SetOptStat       (111111);
  //RootStyle->SetOptStat       ("ne");
  RootStyle->SetStatColor     (0);
  RootStyle->SetStatX         (0);
  RootStyle->SetStatY         (0);
  RootStyle->SetStatFontSize  (0.06);
  RootStyle->SetStatW         (0.2);
  RootStyle->SetStatH         (0.15);
  
  // Labels, Ticks, and Titles

  RootStyle->SetTickLength ( 0.015,"X");
  RootStyle->SetTitleSize  ( 0.055,"X");
  RootStyle->SetTitleOffset( 1.000,"X");
  RootStyle->SetTitleBorderSize(0);
  //  RootStyle->SetTitleFontSize((float)3.);
  RootStyle->SetLabelOffset( 0.015,"X");
  RootStyle->SetLabelSize  ( 0.040,"X");
  RootStyle->SetLabelFont  ( 42   ,"X");

  RootStyle->SetTickLength ( 0.015,"Y");
  RootStyle->SetTitleSize  ( 0.055,"Y");
  RootStyle->SetTitleOffset( 1.00,"Y");
  RootStyle->SetLabelOffset( 0.015,"Y");
  RootStyle->SetLabelSize  ( 0.040,"Y");
  RootStyle->SetLabelSize  ( 0.040,"Z");
  RootStyle->SetLabelFont  ( 42   ,"Y");
  RootStyle->SetLabelFont  ( 42   ,"Z");
  RootStyle->SetTitleOffset(1.00  ,"Z");
  RootStyle->SetTitleSize  ( 0.055,"Z");
  RootStyle->SetTitleFont  (42,"XYZ");
  RootStyle->SetTitleColor  (1);

  // Options

  RootStyle->SetOptFit     (1);

  RootStyle->SetMarkerStyle(20);
  RootStyle->SetMarkerSize(.8);

  //  cout << ">> Style initialized with the Root Style!" << endl;
  //  cout << ">> " << modified << endl << endl;
  return RootStyle;
}


