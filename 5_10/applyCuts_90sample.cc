#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cmath>
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
#include <iomanip>
#include "MyCorrelator.h"

//#include <fftw3.h>

//#ifdef ANITA_UTIL_EXISTS
#include "UsefulAnitaEvent.h"
#include "AnitaGeomTool.h"
#include "AnitaConventions.h"
#include "RawAnitaHeader.h"
#include "analysis_info_4pol.h"
RawAnitaHeader *fHeadPtr = 0;
UsefulAnitaEvent *fUsefulEventPtr;
UsefulAnitaEvent *fUsefulEventPtr_CW;
analysis_info_4pol *pol4_Ptr = 0;
analysis_info_4pol *simpol4_Ptr = 0;
//#endif
class MyCorrelator;
using namespace std;
TStyle* RootStyle();
TStyle *color=RootStyle();
/*
#include "Tools.h"
#include "Constants.h"
#include "Vector.h"
#include "Position.h"
#include "EarthModel.h"
#include "IceModel.h"
#include "Efficiencies.h"
#include "Spectra.h"
#include "Event.h"
#include "Trigger.h"
#include "Detector.h"
#include "Settings.h"
#include "counting.hh"
#include "Primaries.h"
#include "Report.h"

#include "Ray.h"

class EarthModel; //class
*/

 

int main() { 
  MyCorrelator *magicPtr = new MyCorrelator("/data/anita/anitaplus/AnitaFlight0809RootData",12, WaveCalType::kDefault);
   gStyle=color;
   gStyle->SetPalette(1);
   gStyle -> SetOptStat(0000000000000);

   int index;
   int eventNumber_cut=0;
   double peakHilbertCoherent=0.;
   double peakVal=0.;
   double ratioFirstToSecondPeak=0.;
   double thetaMap=0.;
   double phiMap=0.;
   double polFractionCoherent=0.;
   double snrcoherent=0.;
   int hwTriggerFlag=0;
   int varnerFlag=0;
   int eventTracedFlag2=0;

   int varnerFlag2;
  
   float nadirFlag;
   
   double AnitaLat;
   double AnitaLon;
   double AnitaAlt;
   double heading;

   double simAnitaLat;
   double simAnitaLon;
   double simAnitaAlt;
   double simheading;

   int badNoiseFlag;

   
   float didIFilter;
   float didIFilterHoriz;
   float didIFilterAboveSatelliteHoriz;
   float didIFilterAboveSatellite;
    
   int payloadBlastctr=0;
   int hilbertctr=0;
   int ratiopeaksctr=0;
   int crosscorrctr=0;
   int polfractionctr=0;
   int rotatedctr=0;
   int elevationctr=0;
   int triggeredctr=0;
   int tracedctr=0;
   int badNoisectr=0;
   
   int varnerevents=0;
   int datactr=0;

     int ratio_last=0;
   int peakVal_last=0;
   int hilbert_last=0;
   int polfrac_last=0;
   int rotated_last=0;
   int elevation_last=0;
   int traced_last=0;
   int triggered_last=0;

   double snrPeak;

   double sourceLon;
   double sourceLat;
   double sourceAlt;
   int false_neg_ctr=0;
   TH1D *hpeakVal;
   TH1D *hratiopeaks;
   TH1D *hhilbert;
   TH1D *hpolfrac;
   TH1D *helevation;
   
  int numbins = 50;
   int index_mult = 1;
   // numbins=50;
   // index_mult=1;
   //create histograms that are filled at the end of the code
   hratiopeaks = new TH1D("ratiopeaks","; Ratio 2nd/1st Peak;Number of Events",100,0,1);
   hpeakVal = new TH1D("peakVal",";Peak Val of Cross Correlation;Number of Events",100,0,1);
   hhilbert = new TH1D("hilbert",";Peak of Coherently Summed Waveform",250,0,500);
   hpolfrac = new TH1D("polfrac","; Polarization Fraction;Number of Events",100,0,1);
   helevation = new TH1D("elevation",";Elevation Angle;Number of Events",180,-90,90);
  
  

   TH1D *hpeakVal_sim;
   TH1D *hratiopeaks_sim;
   TH1D *hhilbert_sim;
   TH1D *hpolfrac_sim;
   TH1D *helevation_sim;
  
   hratiopeaks_sim = new TH1D("ratiopeaks sim","; Ratio 2nd/1st Peak ;Number of Events",100,0,1);
   hpeakVal_sim = new TH1D("peakVal sim",";Peak Val of Cross Correlation;Number of Events",100,0,1);
   
   hhilbert_sim = new TH1D("hilbert sim",";Peak of Coherently Summed Waveform",250,0,500);
   
   hpolfrac_sim = new TH1D("polfrac sim","; Polarization Fraction;Number of Events",100,0,1);
   
   
   helevation_sim = new TH1D("elevation sim",";Elevation Angle;Number of Events",180,-90,90);

   hratiopeaks_sim->SetLineColor(kBlue);
   hpeakVal_sim->SetLineColor(kBlue);
   hhilbert_sim->SetLineColor(kBlue);
   hpolfrac_sim->SetLineColor(kBlue);
   helevation_sim->SetLineColor(kBlue);
    
   TH2D *hrotated = new TH2D("rotated",";Peak Val;SNR of Coherent Waveform",100,0,1,5000,0,500);

   double weight;
   UInt_t eventNumber_icemc;
   int n_icemc=0;
   double lat_icemc;
   double lon_icemc;
  
   double arrival_times[48];
   double peak2peak_max=0.;
  
   
  //read in root file
   
    int  run=12;
  string temp = "/data/anita/btdailey/final_filter/90sample/geom_4pol_partial_0301/output%d_0.root";
  
  
  char *rootfile = Form(temp.c_str(),run);
  
 
  TChain *pol4_Tree = new TChain("analysis_info_4pol");
  TChain *Pointed_Tree =  new TChain("tdataTracing");
  // pol4_Tree->Add(rootfile);
  //Pointed_Tree->Add(rootfile);
  //extras
    for(int extra=13;extra<263;extra++){
	 if(extra==218 || extra==261){
       for(int delta=0;delta<1000000;delta+=1000){
	
	sprintf(rootfile,"/data/anita/btdailey/final_filter/90sample/geom_4pol_partial_kingbee_0301/output%d_%d.root",extra,delta);
	
	TFile *fpTest = TFile::Open(rootfile);
	if(!fpTest){ 
	  //cout<<"broke at extra "<<extra<<"_"<<delta<<"\n";
	  extra++;
	  break;
	  //break;
	}
	else {
	  delete fpTest;
	  
	  
	  
	  pol4_Tree->Add(rootfile);
	  Pointed_Tree->Add(rootfile);
	}
	
       }
    }
	 for(int delta=0;delta<2000000;delta+=9000){
	   
	   sprintf(rootfile,"/data/anita/btdailey/final_filter/90sample/geom_4pol_partial_0301/output%d_%d.root",extra,delta);
	   
	 TFile *fpTest = TFile::Open(rootfile);
	 if(!fpTest){ 
	   //cout<<"broke at extra "<<extra<<"_"<<delta<<"\n";
	   break;
	   //break;
	 }
	 else {
	   delete fpTest;
	   
	   
	   
	   pol4_Tree->Add(rootfile);
	   Pointed_Tree->Add(rootfile);
	 }
	  }//delta
       }//extra
  
    /*for(int extra=200;extra<261;extra++){
      
      for(int delta=0;delta<15000;delta+=1000){
	
	sprintf(rootfile,"/data/anita/btdailey/final_filter/10sample/geom_4pol_partial_1213/output%d_%d.root",extra,delta);
	
	TFile *fpTest = TFile::Open(rootfile);
	if(!fpTest){ 
	  cout<<"broke at extra "<<extra<<"_"<<delta<<"\n";
	  break;
	  //break;
	}
	else {
	  delete fpTest;
	   
	   
	   
	   pol4_Tree->Add(rootfile);
	   Pointed_Tree->Add(rootfile);
	 }
	  }//delta
       }//extra
    
      for(int extra=261;extra<263;extra++){
	 
	 for(int delta=0;delta<20000;delta+=3000){
	   
	   sprintf(rootfile,"/data/anita/btdailey/final_filter/10sample/geom_4pol_partial_1213/output%d_%d.root",extra,delta);
	   
	 TFile *fpTest = TFile::Open(rootfile);
	 if(!fpTest){ 
	   cout<<"broke at extra "<<extra<<"_"<<delta<<"\n";
	   break;
	   //break;
	 }
	 else {
	   delete fpTest;
	   
	   
	   
	   pol4_Tree->Add(rootfile);
	   Pointed_Tree->Add(rootfile);
	 }
	  }//delta
       }//extra
    */
      //////THESE TREES ARE FILLED FOR EVENTS PASSING ALL CUTS!////////////
      char filename[150];
 
      sprintf(filename,"/data/anita/btdailey/passingCuts/90sample_partial_passedcuts_0301_HPol.root");
  
      cout<<"outputting to file: "<<filename<<endl;
      TFile *rootfile_out = new TFile(filename,"RECREATE");
      TTree *eventTree_out = pol4_Tree->CloneTree(0);

      TTree *extraTree = new TTree("extraTree","extraTree");
     
      extraTree->Branch("anitaLat",&AnitaLat,"anitaLat/D");
      extraTree->Branch("anitaLon",&AnitaLon,"anitaLon/D");
      extraTree->Branch("anitaAlt",&AnitaAlt,"anitaAlt/D");
      extraTree->Branch("heading",&heading,"heading/D");

     

     int nevents0 = pol4_Tree->GetEntries();//dataTree->GetEntries();
     
     int nevents_pointed = Pointed_Tree->GetEntries();
     
     pol4_Tree->SetBranchAddress("pol4_Ptr",&pol4_Ptr);
     
     int eventNumber_pointed;
    
     Pointed_Tree->SetBranchAddress("eventNumber",&eventNumber_pointed);
     Pointed_Tree->SetBranchAddress("anitaLatitude",&AnitaLat);
     Pointed_Tree->SetBranchAddress("anitaLongitude",&AnitaLon);
     Pointed_Tree->SetBranchAddress("anitaAltitude",&AnitaAlt);
     Pointed_Tree->SetBranchAddress("heading",&heading);

    
     int mainrfcmflag;
     int bigenoughpeakflag;
     int dcoffsetflag;
     int shorttraceflag;
     int nadirrfcmflag;
     
     
     int n_old=0;


     int polarization =1;
     for(int m=0;m<nevents0;m++){
     //for(int m=1369290;m<nevents0;m++){
       // for(int m=0;m<0;m++){
       	if (m % (nevents0 / 100) == 0){
	  cout << m << " events. " <<(double(m)/double(nevents0)) * 100 << "% complete.\n";
	}
	//get events
      
      
       pol4_Tree->GetEvent(m);
       Pointed_Tree->GetEvent(m);
       
       //cout<<"Anita is at "<<AnitaLat<<" "<<AnitaLon<<" "<<AnitaAlt<<"\n";
       eventNumber_cut = pol4_Ptr->eventNumber;
       peakHilbertCoherent=pol4_Ptr->peakHilbertCoherent[polarization];
       peakVal=pol4_Ptr->peakVal[polarization];
       ratioFirstToSecondPeak=pol4_Ptr->ratioFirstToSecondPeak[polarization];
       thetaMap=pol4_Ptr->thetaMap[polarization];
       snrcoherent = pol4_Ptr->SNR_coherent[polarization];
       polFractionCoherent=pol4_Ptr->polFractionCoherent[polarization];
      
       //if(peakHilbertCoherent != peakHilbertCoherent) cout<<"eventnumber is "<<eventNumber_cut<<" snr is "<<snrcoherent<<"\n";
       hwTriggerFlag=pol4_Ptr->hwTriggerFlag[polarization];
       snrPeak = pol4_Ptr->snrPeakAfterFilter[polarization];
       
       varnerFlag = pol4_Ptr->varnerFlag[polarization];
       didIFilter = pol4_Ptr->didIFilter[polarization];
       eventTracedFlag2 = pol4_Ptr->eventTracedFlag[polarization];
       varnerFlag2 = pol4_Ptr->varnerFlag2[polarization];
      
      
       mainrfcmflag=  pol4_Ptr->mainrfcmflag;
       bigenoughpeakflag=pol4_Ptr->bigenoughpeakflag;
       dcoffsetflag=pol4_Ptr->dcoffsetflag;
       shorttraceflag=pol4_Ptr->shorttraceflag;
       nadirrfcmflag=pol4_Ptr->nadirrfcmflag;
       badNoiseFlag = pol4_Ptr->noiseFlag[polarization];
       if ( !mainrfcmflag || shorttraceflag !=0 || dcoffsetflag != 0 || bigenoughpeakflag !=1){
	 cout<<"quality cut event \n";
	 continue;
       }
	

       //PAYLOAD BLASTS
       if(pol4_Ptr->payloadblastflag ==1){
	 payloadBlastctr++;
	 continue;
       }
       
       index = (int) snrPeak;

      
       if(index >= numbins){
	 index=numbins;
	 
       }
       
       hpeakVal->Fill(peakVal);
       hratiopeaks->Fill(1./ratioFirstToSecondPeak);
       hhilbert->Fill(peakHilbertCoherent);
       hpolfrac->Fill(polFractionCoherent);
       helevation->Fill(thetaMap);
       
   

       float limit=0.9;

       
       if (didIFilterHoriz>1 || didIFilter>1){
	 limit=.85;
       }
       
       if(eventTracedFlag2!=1){
	 eventTracedFlag2=magicPtr->traceBackToContinent_Brian2(eventNumber_cut,pol4_Ptr->peakThetaFinal[polarization],pol4_Ptr->peakPhiFinal[polarization], sourceLon, sourceLat, sourceAlt,AnitaLat,AnitaLon,AnitaAlt,heading);
	 pol4_Ptr->sourceLon[polarization]=sourceLon;
	 pol4_Ptr->sourceLat[polarization]=sourceLat;
	 pol4_Ptr->sourceAlt[polarization]=sourceAlt;

	 if(eventTracedFlag2==1){
  
	   false_neg_ctr++;
	 }
       }
       
       //do ctr for  last cuts////

       if(peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap <0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	 if(ratioFirstToSecondPeak<=(1/limit)){
	     ratio_last++;//ratio last
	   }
       }
       if(ratioFirstToSecondPeak > (1/limit) && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap <0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	  if(peakVal<0.075){
	    peakVal_last++;//peakVal last
	   }
       }
       if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075  && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	  if(peakHilbertCoherent <=15){
	    hilbert_last++;//hilbert last
	   }
       }
     if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && polFractionCoherent>0.3){
	 if(peakHilbertCoherent<-350*peakVal+57.14){
	   rotated_last++; //rotated last
	 }
       }
       if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14)){
	 if(polFractionCoherent<=.3){
	   polfrac_last++;//rotated last	
       }
      
       }
       if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && eventTracedFlag2==1 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	  if(thetaMap<=-35 || thetaMap>=0){
	    elevation_last++;//elevation last
       }
      
       }

       if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap <0 && hwTriggerFlag !=0 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	 if(eventTracedFlag2!=1){
	   traced_last++;//traced Last
	 }
       }

       if(ratioFirstToSecondPeak > (1/limit) && peakVal >= 0.075 && peakHilbertCoherent >15 && thetaMap >-35 && thetaMap < 0 && eventTracedFlag2==1 && varnerFlag !=1 && peakHilbertCoherent>=(-350*peakVal+57.14) && polFractionCoherent>0.3){
	 if(hwTriggerFlag==0){
	   //cout<<"eventnumber being cut by hwtrigger is "<<eventNumber<<" hwTriggerFlag is "<<hwTriggerFlag<<"\n";
	   triggered_last++;//triggered_last
	 }
       }
      
       
       if(eventNumber_cut == 27146983) cout<<"inserted event \n";


       /////////////////////////////////////Stuff for plotting and in order cuts
      

       /* if(ratioFirstToSecondPeak <=(1/limit)){
	 ratiopeaksctr++;
	 continue;
       }
       */
       /*
       if(peakVal<.075){
	 //	 if(peakVal<.03){
	 crosscorrctr++;
	 
	 continue;
       }
       */
       if(peakHilbertCoherent <=15){
	 hilbertctr++;
	 continue;
       }
       	     
       /*
        if(polFractionCoherent<=.3){
	 polfractionctr++;
	 continue;
       }
       */
  
       
       // if(peakHilbertCoherent<-350*peakVal+57.14){
       // if(peakHilbertCoherent<-100*peakVal+10){
       // if(peakHilbertCoherent<-117*peakVal+14){
       // rotatedctr++;
       //	 continue;
       // }

      

       if(thetaMap<=-35 || thetaMap>=0){
	 elevationctr++;
	 continue;
       }
      
  
        if(eventTracedFlag2!=1){
	 tracedctr++;
	 continue;
       }
      
  
       if(hwTriggerFlag==0){
	 triggeredctr++;
	 continue;
       }
      
   
       if(varnerFlag==1){// || varnerFlag2==1){
	 varnerevents++;
	 continue;
       }

        if(badNoiseFlag==1){
	 badNoisectr++;
	 continue;
       }
	cout<<"eventnumber is "<<eventNumber_cut<<"\n";
	if(snrcoherent >20) cout<<"eventnumber is "<<eventNumber_cut<<" snr is "<<snrcoherent<<" noiseflag is "<<badNoiseFlag<<"\n";
       hrotated->Fill(peakVal,snrcoherent);
       /*  
       for(int n= n_old;n<nevents_pointed;n++){
	 Pointed_Tree->GetEvent(n);
	 //cout<<"eventNumber_pointed is "<<eventNumber_pointed<<" eventNumber is "<<pol4_Ptr->eventNumber<<"\n";
	 if(eventNumber_pointed == eventNumber_cut) {
	   n_old=n;
	   break;
	 }
	 
       }
       */
       
       datactr++;
       eventTree_out->Fill();
       extraTree->Fill();
       // cout<<"eventnumber cut is "<<eventnumber2<<"\n";
	
     }//nevents0

   cout<<"10sample cuts! /////////////// \n";
   cout<<"payload blast cut is: "<<payloadBlastctr<<"\n";
   cout<<"false_neg is "<<false_neg_ctr<<"\n";
   cout<<"badnoisectr is "<<badNoisectr<<"\n";
   cout<<"ratio of peaks cut is: "<<ratiopeaksctr<<" :"<<ratio_last<<" \n";
   cout<<"cross corr cut is: "<< crosscorrctr<<" :"<<peakVal_last<<"\n";
   cout<<"number of events with hilbert <15 is: "<<hilbertctr<<" :"<<hilbert_last<<"\n";
   cout<<"polfraction cut is: "<< polfractionctr<<" :"<<polfrac_last<<"\n";
   cout<<"rotated cross cut is: "<<rotatedctr<<" :"<<rotated_last<<" \n";
   cout<<"elevation cut is: "<<elevationctr<<" :"<<elevation_last<<"\n";
   cout<<"traced cut is: "<<tracedctr<<" :"<<traced_last<<"\n";
   cout<<"triggered cut is: "<<triggeredctr<<" :"<<triggered_last<<"\n";
   cout<<"varner events is: "<<varnerevents<<"\n";
   cout<<"number that pass all cuts is: "<<datactr<<"\n";

   cout<<"peak2peak_max is "<<peak2peak_max<<"\n";

  

   rootfile_out = eventTree_out->GetCurrentFile();
   rootfile_out->Write();
  


 
 delete hpeakVal;
 delete hratiopeaks;
 delete hhilbert;
 delete hpolfrac;

 delete helevation;
 delete hrotated;

  
}//main
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
  RootStyle->SetPadLeftMargin  (0.16);
  RootStyle->SetPadRightMargin (.2);
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
  RootStyle->SetLabelSize  ( 0.030,"X");
  RootStyle->SetLabelFont  ( 42   ,"X");

  RootStyle->SetTickLength ( 0.015,"Y");
  RootStyle->SetTitleSize  ( 0.050,"Y");
  RootStyle->SetTitleOffset( 1.4,"Y");
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

