
#include "TChain.h"
#include "TFile.h"
#include "Constants.h"
#include "TRandom3.h"
#include "Settings.h"
#include "earthmodel.hh"
#include "icemodel.hh"
#include "vector.hh"
#include "TF1.h"
#include "TH1F.h"

#include "position.hh"
#include "anita.hh"
#include "ray.hh"

#include "balloon.hh"
#include "Settings.h"
#include "Primaries.h"

using std::ifstream;

Balloon::Balloon() {
    MAXHORIZON=800000.;                // pick the interaction within this distance from the balloon so that it is within the horizon
    ibnposition=0;
    igps=0;
    horizcoord_bn=0; // x component of balloon position
    vertcoord_bn=0; // y component of balloon position
    BN_LONGITUDE=999; //balloon longitude for fixed balloon location
    BN_LATITUDE=999; //balloon latitude for fixed balloon location
}

void  Balloon::setObservationLocation(Interaction *interaction1,int inu,IceModel *antarctica,Settings *settings1) {
    interaction1->banana_volts = 0; //Zero the variable
    interaction1->banana_obs = Vector(0,0,Interaction::banana_observation_distance);
    
    // First pick theta of the obs vector from nu vector
    interaction1->banana_theta_obs = -0.5*PI/settings1->vertical_banana_points * ((int)(inu/settings1->horizontal_banana_points));
    interaction1->banana_phi_obs = 2*PI/settings1->horizontal_banana_points * (inu%settings1->horizontal_banana_points);
    interaction1->banana_weight = (double)fabs(sin(interaction1->banana_theta_obs)); //Set "weight" - purely a phase factor
    interaction1->banana_obs = interaction1->banana_obs.RotateY(interaction1->banana_theta_obs);
    
    // Now the phi of the obs vector relative to the nu vector
    interaction1->banana_obs = interaction1->banana_obs.RotateZ(interaction1->banana_phi_obs);
    
    // Finally rotate the obs vector to the nu vector's orientation (change coordinate systems, in effect)
    interaction1->banana_obs = interaction1->banana_obs.RotateY(PI-Interaction::theta_nu_banana);
    interaction1->banana_obs = interaction1->banana_obs.RotateZ(Interaction::phi_nu_banana);
    
    //Finally, set the balloon position to the location which we would like to observe
    r_bn = interaction1->nu_banana_surface + interaction1->banana_obs;
    
    theta_bn = r_bn.Theta();
    phi_bn = r_bn.Phi();
    surface_under_balloon = antarctica->Surface(r_bn);
    // position of balloon at earth's surface (the shadow it would cast with the sun overhead)
    r_bn_shadow = surface_under_balloon * r_bn.Unit(); 
    altitude_bn = r_bn.Mag() - surface_under_balloon;
    
    //Finished setting observation location
}
void Balloon::SetDefaultBalloonPosition(IceModel *antarctica1) { // position of surface of earth under balloon
    
    // set the default balloon position
    // if you are using real Anita-lite path, these get overwritten for each event
    //cout << "BN_LATITUDE is " << BN_LATITUDE << "\n";
    
    int BN_LATITUDE_SETTING = BN_LATITUDE;
    int BN_LONGITUDE_SETTING = BN_LONGITUDE;
    
    if(BN_LATITUDE_SETTING==999)
	theta_bn=10*RADDEG; // wrt south pole
    else 
	theta_bn=(90-BN_LATITUDE_SETTING)*RADDEG; 
    
    if(BN_LONGITUDE_SETTING==999)
	phi_bn=PI/4; //wrt 90E longitude
    else 
	phi_bn=EarthModel::LongtoPhi_0isPrimeMeridian(BN_LONGITUDE_SETTING); //remember input of LongtoPhi is between -180 and 180
    
    
    r_bn = Position(theta_bn,phi_bn); // direction of balloon- right now this is a unit vector
    
    if (BN_ALTITUDE==0) // if the altitude isn't set in the input file
	altitude_bn=120000*12.*CMINCH/100.; // 120000 ft.=36.6 m
    else
	altitude_bn=BN_ALTITUDE*12.*CMINCH/100.; // converts the altitude in the input file to meters
    
    surface_under_balloon = antarctica1->Surface(r_bn); // distance between center of earth and surface under balloon
    
    r_bn_shadow = surface_under_balloon * r_bn.Unit(); // position of surface under balloon
    
    r_bn = (antarctica1->Geoid(r_bn)+altitude_bn) * r_bn; // position of balloon
    
} // set default balloon position


void Balloon::ReadAnitaliteFlight() {
    
    
    
    
    ifstream flightfile(anitaliteflight.c_str()); //set file to the right one
    if (!flightfile) {
	cout << "Flight file not found.\n";
	exit(1);
    }
    
    
    
    string slatitude; // latitude in string format
    string slongitude; //longitude in string format
    string saltitude; // altitude in string format
    string sheading; // heading in string format
    string srealtime; // real life unix time in string format
    string junk; // string not used for anything
    
    string line; // for reading a whole line
    double latitude; // same quantities as above but in double format
    double longitude;
    double altitude;
    double heading;
    double realtime;
    
    
    getline(flightfile,junk); // first few lines are junk
    getline(flightfile,junk);
    getline(flightfile,junk);
    getline(flightfile,junk);
    
    NPOINTS=0; // number of points we have for the flight path
    
    
    
    while (!flightfile.eof()) {
	
	
	flightfile >> srealtime >> junk >> slatitude >> slongitude >> saltitude >> sheading >> junk >> junk >> junk;
	
	
	//cout << "they are " << NPOINTS << " " << srealtime << " " << junk << " " << slatitude << " " << slongitude << " " << saltitude << " " << sheading << "\n";
	//cout << "slatitude is " << slatitude << "\n";
	latitude=(double)atof(slatitude.c_str());
	longitude=(double)atof(slongitude.c_str());
	altitude=(double)atof(saltitude.c_str());
	heading=(double)atof(sheading.c_str());
	realtime=(double)atof(srealtime.c_str());
	
	
	// latitude and longitude of each balloon position
	latitude_bn_anitalite[NPOINTS]=latitude;
	longitude_bn_anitalite[NPOINTS]=longitude;
	altitude_bn_anitalite[NPOINTS]=altitude;
	heading_bn_anitalite[NPOINTS]=heading;
	realtime_bn_anitalite[NPOINTS]=realtime;
	
	
	NPOINTS++;
	
	
	//    cout << "NPOINTS, altitude are " << NPOINTS << " " << altitude << "\n";
	
	getline(flightfile,junk);
	
    }//while
    
    
    
    NPOINTS_MAX=NPOINTS-140; // exclude the fall 
    NPOINTS_MIN=0;
    
    
    
}//ReadAnitaliteFlight
void Balloon::InitializeBalloon() {
    
  WHICHPATH=7;//HACKED!!!!!!!
    // GPS positions of Anita or Anita-lite balloon flight
    if (WHICHPATH==2)
	ReadAnitaliteFlight();
    
    
    cout << "Inside balloon constructor.  whichpath is " << WHICHPATH << " which was hacked by Brian. \n";
    MINALTITUDE=30000; // balloon has to be 30 km altitude at least for us to read the event from the flight data file
    phiTrigMask=0;
    
    // initialisation of igps_previous
    if (WHICHPATH==6 || WHICHPATH==7)
	igps_previous=0; // which entry from the flight data file the previous event was
    if (WHICHPATH==2)
	igps_previous=NPOINTS_MIN; // initialise here to avoid times during launch
    
    // This for Anita 1 flight
    if (WHICHPATH==6) {
	
	flightdatachain = new TChain("foricemc");
	flightdatachain->Add("data/anita1flightdata.root");
	
	flightdatachain->SetBranchAddress("surfTrigBandMask",surfTrigBandMask);
	flightdatachain->SetBranchAddress("powerthresh",powerthresh);
	flightdatachain->SetBranchAddress("meanp",meanp);  
	flightdatachain->SetBranchAddress("longitude",&flongitude);
	flightdatachain->SetBranchAddress("latitude",&flatitude);
	flightdatachain->SetBranchAddress("altitude",&faltitude);
	flightdatachain->SetBranchAddress("realTime_surfhk",&realTime_flightdata);
	flightdatachain->SetBranchAddress("heading",&fheading);
	
    }
    
    // This for Anita 2 flight
    if (WHICHPATH==7) {
	
	flightdatachain = new TChain("adu5PatTree");
	flightdatachain->SetMakeClass(1);
	flightdatachain->Add("data/gpsFileAnita2.root");
	//    flightdatachain->Add("data/anita2gpsdata_2.root");
	
	flightdatachain->SetBranchAddress("longitude",&flongitude);
	flightdatachain->SetBranchAddress("latitude",&flatitude);
	flightdatachain->SetBranchAddress("altitude",&faltitude);
	flightdatachain->SetBranchAddress("heading",&fheading);
	flightdatachain->SetBranchAddress("realTime",&realTime_flightdata);
	
	//cout<<"realTime is "<<realTime_flightdata<<"\n";
	//cout << "Loading file.  n events is " << flightdatachain->GetEntries() << "\n";
	
	
	fturf=new TFile("data/turfrate_icemc.root");
	//turfratechain = new TChain("turfrate_icemc");
	//turfratechain->Add("data/turfrate_icemc.root");
	
	turfratechain=(TTree*)fturf->Get("turfrate_icemc");
	//    turfratechain->Add("data/anita2turfrate_2.root");
	turfratechain->SetMakeClass(1);
	turfratechain->SetBranchAddress("phiTrigMask",&phiTrigMask);
	turfratechain->SetBranchAddress("realTime",&realTime_turfrate);
	//turfratechain->SetBranchAddress("iev",&iturfevent);
	turfratechain->BuildIndex("realTime");
	//tindex = (TTreeIndex*)turfratechai,n->GetTree()->GetTreeIndex() ;
	//tindex = new TTreeIndex(turfratechain, "realTime", "iev");
	//    tindex=new TTreeIndex((TTree*)turfratechain,"realTime","iev");
	
	turfratechain->GetEvent(0);
	realTime_tr_min=realTime_turfrate; // realTime of first event in the file
	turfratechain->GetEvent(turfratechain->GetEntries()-1);
	realTime_tr_max=realTime_turfrate; // realTime of last event in file
	
    }
    
    
    
    
    for (int i=0;i<10000;i++) {
	latitude_bn_anitalite[i]=0;
	longitude_bn_anitalite[i]=0;
	altitude_bn_anitalite[i]=0;
    }
    NPOINTS=0;
    REDUCEBALLOONPOSITIONS=100;
    
    
    
    if (WHICHPATH!=6) {
	// if it's not the anita path, you won't be able to read in thresholds and masks
	for (int i=0;i<9;i++) {
	    for (int j=0;j<2;j++) {
		surfTrigBandMask[i][j]=0;
	    }
	    for (int j=0;j<32;j++) {
		powerthresh[i][j]=3.E-12;
		meanp[i][j]=1.E-12;
	    }
	    
	}
	
	
    }
}
double Balloon::GetBalloonSpin(double heading) { // get the azimuth of the balloon
    
    
    double phi_spin;
    if (WHICHPATH==2 || WHICHPATH==6 || WHICHPATH==7)
	phi_spin=-1.*heading*RADDEG;
    else { 
	if (RANDOMIZE_BN_ORIENTATION==1) 
	    phi_spin=gRandom->Rndm()*2*PI;
	else
	    phi_spin=0.;
    }
    
    
    return phi_spin;
}


int Balloon::Getibnposition() {
    int ibnposition_tmp;
    if (WHICHPATH==1)
	ibnposition_tmp = (int)(r_bn.Lon() / 2);
    else if (WHICHPATH==2 || WHICHPATH==6 || WHICHPATH==7)
	ibnposition_tmp=(int)((double)igps/(double)REDUCEBALLOONPOSITIONS);
    else
	ibnposition_tmp=0;
    
    return ibnposition_tmp;
    
}
// this is called for each neutrino
void Balloon::PickBalloonPosition(IceModel *antarctica1,int inu) { // Settings *settings1,,Anita *anita1, r_bn_shadow=position of spot under the balloon on earth's surface
    
    //flightdatatree->SetBranchAddress("surfTrigBandMask",&surfTrigBandMask);
    //unsigned short test[9][2];
    
    //double latitude,longitude,altitude;
    
    
    //  cout << "I'm in pickballoonposition. whichpath is " << WHICHPATH << "\n";
    //Pick balloon position
    if (WHICHPATH==2 || WHICHPATH==6 || WHICHPATH==7) { // anita-lite or anita-I or -II path
        
	if (WHICHPATH==2) {
	    igps=NPOINTS_MIN+(igps_previous+1-NPOINTS_MIN)%(NPOINTS_MAX-NPOINTS_MIN); //Note: ignore last 140 points, where balloon is falling - Stephen
	    flatitude=(float)latitude_bn_anitalite[igps];
	    flongitude=(float)longitude_bn_anitalite[igps];
	    faltitude=(float)altitude_bn_anitalite[igps];
	    fheading=(float)heading_bn_anitalite[igps];
	    
	    
	}
	else if (WHICHPATH==6 || WHICHPATH==7) {  // For Anita 1 and Anita 2:  
	    
	    igps=(igps_previous+1)%flightdatachain->GetEntries(); // pick which event in the tree we want
	    
	    
	    flightdatachain->GetEvent(igps); // this grabs the balloon position data for this event
	    //cout<<"realTime is "<<realTime_flightdata<<"\n";
	    
	    while (faltitude<MINALTITUDE || fheading<0) { // if the altitude is too low, pick another event.
		igps++; // increment by 1
		igps=igps%flightdatachain->GetEntries(); // make sure it's not beyond the maximum entry number
		
		flightdatachain->GetEvent(igps);	  // get new event
	    }
	    if (WHICHPATH==7)  // this is for Anita 2
		// get phi masking	
		setphiTrigMask();// set phiTrigMask, and public variable of Balloon class
	    
	}
	igps_previous=igps;
	
	latitude=(double)flatitude;
	longitude=(double)flongitude;
	altitude=(double)faltitude;
	heading=(double)fheading;
	
	setr_bn(latitude,longitude); // sets theta_bn, phi_bn and r_bn.  r_bn is a unit vector pointing in the right direction
	
	
	
	if (WHICHPATH==2)
	    altitude_bn=altitude*12.*CMINCH/100.;
	else if (WHICHPATH==6 || WHICHPATH==7)
	    altitude_bn=altitude; // get the altitude of the balloon in the right units
	
	surface_under_balloon = antarctica1->Surface(r_bn); // get altitude of the surface under the balloon
	
	r_bn_shadow = surface_under_balloon * r_bn.Unit(); // this is a vector pointing to spot just under the balloon on the surface (its shadow at high noon)
	r_bn = (antarctica1->Geoid(r_bn)+altitude_bn) * r_bn.Unit();
	//r_bn = (antarctica->Surface(r_bn)+altitude_bn) * r_bn.Unit(); //this points to balloon position (not a unit vector)
    } //if (ANITA-lite path) or anita 1 or anita 2
    
    
    
    else if (WHICHPATH==1){ // pick random phi at 80 deg S
	phi_bn=gRandom->Rndm()*TWOPI;
	
	r_bn = Position(theta_bn,phi_bn);
	surface_under_balloon = antarctica1->Surface(r_bn);
	
	r_bn_shadow = surface_under_balloon * r_bn.Unit();
	//    r_bn = (antarctica->Geoid(r_bn)+altitude_bn) * r_bn.Unit();
	r_bn = (antarctica1->Surface(r_bn)+altitude_bn) * r_bn.Unit();
	igps=0;
    } //else if(random position at 80 deg S)      
    else if (WHICHPATH==0){
	
	igps=0; // set gps position to be zero - we're at a fixed balloon position
	
    }
    else if (WHICHPATH==5) { // for slac?
	igps=0;
	theta_bn=1.*RADDEG; // 1deg
	phi_bn=1.*RADDEG; // 1deg
	r_bn=Position(theta_bn,phi_bn); // sets r_bn
	if (BN_ALTITUDE!=0)
	    altitude_bn=BN_ALTITUDE*12.*CMINCH/100.; // set the altitude of the balloon to be what you pick.  This isn't in time for CreateHorizons though!  
	surface_under_balloon = antarctica1->Surface(r_bn);
	r_bn_shadow = surface_under_balloon * r_bn.Unit();
	//    r_bn = (antarctica->Geoid(r_bn)+altitude_bn) * r_bn.Unit();     
	
	r_bn = (antarctica1->Surface(r_bn)+altitude_bn) * r_bn.Unit();
    } // you pick it     
    else if (WHICHPATH==4){ // for making comparison with Peter's event
	//phi_bn=0.767945;
	// phi_bn=0.95993;
	//166.73 longitude
	//    phi_bn=-1.339;// this is McM phi
	// 156.73
	phi_bn=-1.1646;// this is McM phi
	//phi_bn=-0.523;
	// theta_bn has been set to 10 degrees in setdefaultballoonposition
	// we want theta_bn to be on a path along same longitude as
	// McM but up to 10 deg further south in latitude
	// 12.12 deg. is McM
	double theta_max=0.2115; // this is McM theta
	double theta_min=0.2115-0.110; // this approximately one horizon north of McM
	//This is 701.59 km along the surface of the earth
	int npositions=1000; // number of possible positions between those two bounds
	theta_bn=theta_min+(theta_max-theta_min)*(double)(inu%npositions)/(double)npositions;
        
	r_bn = Position(theta_bn,phi_bn);
	surface_under_balloon = antarctica1->Surface(r_bn);
	
	r_bn_shadow = surface_under_balloon * r_bn.Unit();
	//    r_bn = (antarctica->Geoid(r_bn)+altitude_bn) * r_bn.Unit();
	
	r_bn = (antarctica1->Surface(r_bn)+altitude_bn) * r_bn.Unit();
	igps=0;
    } // comparing with Peter's event
    
    
    
    ibnposition=Getibnposition();
    /*
    if (!settings1->UNBIASED_SELECTION)
	dtryingposition=antarctica1->GetBalloonPositionWeight(ibnposition);
    else
	dtryingposition=1.;
    
    phi_spin=GetBalloonSpin(heading); // get the azimuth of the balloon.
    
    // normalized balloon position
    n_bn = r_bn.Unit();
    
    // finding which direction is east under the balloon
    n_east = Vector(sin(phi_bn), -1*cos(phi_bn), 0.);
    
    if (settings1->SLAC) 
	AdjustSlacBalloonPosition(inu); // move payload around like we did at slac
    // find position of each antenna boresight
    if (settings1->BORESIGHTS) 
	GetBoresights(settings1,anita1);
    
    // now finding north
    n_north = n_bn.Cross(n_east);
    
    // these coordinates are for filling ntuples.
    horizcoord_bn=r_bn[0]/1000.; //m to km
    vertcoord_bn=r_bn[1]/1000.;
    
    calculate_antenna_positions(settings1, anita1);
    */
    
} // end PickBalloonPosition
void Balloon::AdjustSlacBalloonPosition(int inu) { // move payload around like we did at slac
    
    if (inu<=MAX_POSITIONS)  // use the event number for the position if we c an
	islacposition=inu;
    else
	islacposition=0; // otherwise just use the first positions
    // adjust the payload position
    r_bn+=slacpositions[islacposition].GetX()*n_east
    +slacpositions[islacposition].GetY()*n_bn;
    
}
void Balloon::CenterPayload(Settings *settings1,Anita *anita1,Vector n_exit2bn,Vector n_pol,int whichlayer,int centerphi) {
    
    double hitangle_e,hitangle_h,e_component,h_component;
    Vector ant_normal;
    
    // allow an option to rotate the payload so the signal is
    // always on the boresight of one phi sector
    GetHitAngles(settings1,anita1,0,n_exit2bn,n_pol,whichlayer,centerphi,
		 hitangle_e,hitangle_h,e_component,h_component,ant_normal);   
    phi_spin-=hitangle_e;
    
}

void Balloon::GetHitAngles(Settings *settings1,Anita *anita1,int inu,const Vector &n_exit2bn,const Vector &n_pol,int ilayer,int ifold,
			   double& hitangle_e,double& hitangle_h,double& e_component,double& h_component,Vector &ant_normal) {
    
    // this find the angles that the signal hits the antenna relative to the e and h polarizations
    Vector n_eplane(0,0,1);
    Vector n_hplane(0,-1,0);
    Vector n_normal(1,0,0); // face of antenna starts out relative to +x because phi is relative to +x
    
    // now rotate for antenna's orientation on the payload.
    
    if(settings1->WHICH==6 || settings1->WHICH==8) {
	n_eplane = n_eplane.RotateY(anita1->ANTENNA_DOWN[ilayer][ifold]);
	n_hplane = n_hplane.RotateY(anita1->ANTENNA_DOWN[ilayer][ifold]);
	n_normal = n_normal.RotateY(anita1->ANTENNA_DOWN[ilayer][ifold]);
    }
    else {
	n_eplane = n_eplane.RotateY(anita1->THETA_ZENITH[ilayer] - PI/2);
	n_hplane = n_hplane.RotateY(anita1->THETA_ZENITH[ilayer] - PI/2);
	n_normal = n_normal.RotateY(anita1->THETA_ZENITH[ilayer] - PI/2);
    }
    
    
    double phi;
    // rotate about z axis for phi
    if (settings1->CYLINDRICALSYMMETRY==1) {
	phi=(double)ifold/(double)anita1->NRX_PHI[ilayer]*2*PI + anita1->PHI_OFFSET[ilayer] + phi_spin;
    }
    else
	phi=anita1->PHI_EACHLAYER[ilayer][ifold] + anita1->PHI_OFFSET[ilayer]+phi_spin;
    
    
    
    n_eplane = n_eplane.RotateZ(phi);
    n_hplane = n_hplane.RotateZ(phi);
    n_normal = n_normal.RotateZ(phi);
    //when phi_spin=0, n_normal is in the x-z plane for antenna 0 (0-31)
    
    //   now rotate to balloon's position on the continent
    //n_eplane = n_eplane.ChangeCoord(r_bn);
    //n_hplane = n_hplane.ChangeCoord(r_bn);
    //n_normal = n_normal.ChangeCoord(r_bn); // These two versions of ChangeCoord gave different results. ChangeCoord(balloon_pos) is simpler and I didn't find any mistakes with it so I'll continue to use it for the ANITA Monte Carlo at OSU. We should make the two simulations use the same version of ChangeCoord after someone finds the bug in one of them. ---Brian
    n_eplane = n_eplane.ChangeCoord(n_north,-1.*n_east);
    n_hplane = n_hplane.ChangeCoord(n_north,-1.*n_east);
    n_normal = n_normal.ChangeCoord(n_north,-1.*n_east);
    // n_normal points northwards for antenna 0 (0-31)
    
    // find component along the e-plane for the purpose of finding hit angles
    e_component = -(n_exit2bn*n_eplane);
    // find component along the h-plane for the purpose of finding hit angles
    h_component = -(n_exit2bn*n_hplane);
    // find the component normal 
    double n_component = -(n_exit2bn * n_normal);
    
    
    
    hitangle_e=atan(h_component/n_component);
    
    if (n_component<0 && h_component<0)
	hitangle_e-=PI;
    if (n_component<0 && h_component>0)
	hitangle_e+=PI;  
    
    hitangle_h=atan(e_component/n_component);
    
    if (n_component<0 && e_component<0)
	hitangle_h-=PI;
    if (n_component<0 && e_component>0)
	hitangle_h+=PI;
    
    //Record antenna normal for plotting - Stephen
    ant_normal = n_normal;  
    
    e_component = n_pol*n_eplane;
    
    // find component along the h-plane  
    h_component = n_pol*n_hplane;
    
    
    if (settings1->REMOVEPOLARIZATION) {
	//Trying to remove effects of polarization at antenna. Stephen
	e_component = n_pol*n_pol;
	h_component = 0.001;
	n_component = 0.001;
    } //if
    
} //end void GetHitAngles


void Balloon::setphiTrigMask() {
    if (realTime_flightdata<realTime_tr_min || realTime_flightdata>realTime_tr_max) {
	phiTrigMask=0; // if the realTime for this balloon position is out of range then just set mask to 0
	
    }
    else { // if it's in range
	
	
	iturf=turfratechain->GetEntryNumberWithBestIndex(realTime_flightdata); // find entry in turfratechain that is closest to this realTime_flightdata
	
	
	if (iturf<0) // if it didn't find one
	    phiTrigMask=0; // set to zero
	
	else 
	    turfratechain->GetEvent(iturf);
	
	
    } // end if it's in range
    
}

void Balloon::setr_bn(double latitude,double longitude) {
    
    // latitude is between -90 and 0.  
    // theta_bn measured from the SP and is between 0 and PI/2.
    theta_bn=(90+latitude)*RADDEG;
    
    // this is the payload's longitude, not the azimuth of the balloon like it sounds.
    // longitude is between -180 to 180 with 0 at prime meridian
    // phi is from 0 to 360 with 0 at +90 longitude
    phi_bn=(-1*longitude+90.);
    
    if (phi_bn<0)
	phi_bn+=360.; 
    
    phi_bn*=RADDEG;
    
    r_bn = Position(theta_bn,phi_bn);  //r_bn is a unit vector pointing in the right direction
}

void Balloon::PickDownwardInteractionPoint(Interaction *interaction1,Anita *anita1,Settings *settings1,IceModel *antarctica1,int inu, 
					   
					   
					   Ray *ray1,
					   TH1F *prob_eachbin, TH1F *vol_eachbin, TH1F *prob_eachbin_weighted, TH1F *prob_eachbnposition,
					   TH1F *diff_ilon_bn_ibnposition, int &beyondhorizon, int& k, int& k1) {
    
    double distance=1.E7;
    double phi=0,theta=0;
    double lon=0;
    double latfromSP=0;
    
    
    if (settings1->UNBIASED_SELECTION==1) {
	
      if (antarctica1->PickUnbiased(inu,interaction1,antarctica1,k,k1)) { // pick neutrino direction and interaction point 
	    interaction1->dtryingdirection=1.;
	    interaction1->iceinteraction=1;
	}
      else{
	    interaction1->iceinteraction=0;
      }
    }//unbiased selection
    else {
	interaction1->iceinteraction=1;
	if (WHICHPATH==3) { //Force interaction point if we want to make a banana plot 
	    interaction1->posnu = interaction1->nu_banana;
	} //if (making banana plot)
	else if (WHICHPATH==4) {// Force interaction point for comparison with Peter.
	    
	    // want interaction location at Taylor Dome
	    // According to lab book it's 77 deg, 52.818' S=77.8803
	    // 158 deg, 27.555' east=158.45925
	    latfromSP=12.1197; // this is degrees latitude from the south pole 
	    //lon=180.+166.73;
	    lon=180.+158.45925;
	    //lon=180.+120.20;
	    
	    phi=EarthModel::LongtoPhi_0is180thMeridian(lon);
	    
	    theta = latfromSP*RADDEG;
	    
	    
	    
	    double elevation=antarctica1->SurfaceAboveGeoid(lon,latfromSP)-500.;
	    
	    interaction1->posnu = Vector((elevation+antarctica1->Geoid(latfromSP))*sin(theta)*cos(phi),(elevation+antarctica1->Geoid(latfromSP))*sin(theta)*sin(phi),(elevation+antarctica1->Geoid(latfromSP))*cos(theta));
	    
	    
	    
	} //if (WHICHPATH==4)
	
	else if (settings1->SLAC) {
	    
	    Vector zaxis(0.,0.,1.); // start with vector pointing in the +z direction
	    
	    interaction1->posnu=zaxis.RotateY(r_bn.Theta()-settings1->SLAC_HORIZDIST/EarthModel::R_EARTH); // rotate to theta of balloon, less the distance from the interaction to the balloon
	    
	    interaction1->posnu=interaction1->posnu.RotateZ(r_bn.Phi()); // rotate to phi of the balloon
	    
	    interaction1->posnu=(antarctica1->Surface(interaction1->posnu)-settings1->SLAC_DEPTH)*interaction1->posnu; // put the interaction position at depth settings1->SLAC_DEPTH in the ice 
	    
	}//slac
	else
	    interaction1->posnu = antarctica1->PickInteractionLocation(ibnposition,inu);
    }//else unbiased


    if (interaction1->iceinteraction) {
	//cout << "posnu is ";interaction1->posnu.Print();
    }
    // first guess at the rf exit point is just the point on the surface directly above the interaction
    ray1->rfexit[0] = antarctica1->Surface(interaction1->posnu) * interaction1->posnu.Unit();
    
    
    // unit vector pointing to antenna from exit point.
    ray1->n_exit2bn[0] = (r_bn - ray1->rfexit[0]).Unit();
    
    if (settings1->BORESIGHTS) {
	for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
	    for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
		ray1->rfexit_eachboresight[0][ilayer][ifold]=antarctica1->Surface(interaction1->posnu) * interaction1->posnu.Unit();// this first guess rfexit is the same for all antennas too
		ray1->n_exit2bn_eachboresight[0][ilayer][ifold]=(r_boresights[ilayer][ifold]- ray1->rfexit_eachboresight[0][ilayer][ifold]).Unit();
		//cout << "ilayer, ifold, n_exit2bn are " << ilayer << "\t" << ifold << " ";
	    }
	}
    }
    
    // first pass at direction of vector from interaction to exit point
    // just make the ray go radially outward away from center of earth.    
    // still for first guess
    // should incorporate this into PickInteractionPoint	 
    ray1->nrf_iceside[0] = interaction1->posnu.Unit();
    
    if (settings1->BORESIGHTS) {
	// this is the same for all of the antennas too
	for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) { // loop over layers on the payload
	    for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
		
		ray1->nrf_iceside_eachboresight[0][ilayer][ifold]=interaction1->posnu.Unit();
	    } // end loop over fold
	} // end loop over payload layers
    } // end if boresights  
    
    
    
    
    double r_down = 2*(antarctica1->Surface(interaction1->posnu)-antarctica1->IceThickness(interaction1->posnu))-interaction1->posnu.Mag();
    interaction1->posnu_down = r_down * interaction1->posnu.Unit();
    //position of the mirror point of interaction
    
    //interaction1->posnu is downward interaction1->posnu.
    distance=interaction1->posnu.Distance(r_bn);
    
    
    // depth of interaction
    // gets distance between interaction and exit point, this time it's same as depth
    // because our first guess at exit point is radially outward from interaction.
    // negative means below surface
    interaction1->altitude_int=-1*ray1->rfexit[0].Distance(interaction1->posnu);
    interaction1->altitude_int_mirror=-1*ray1->rfexit[0].Distance(interaction1->posnu_down);//get depth of mirror point
    
    
    interaction1->r_fromballoon[0]=r_bn.Distance(interaction1->posnu);
    
    //distance from the mirror point to the balloon because it is equal to the path that signals pass
    interaction1->r_fromballoon[1]=r_bn.Distance(interaction1->posnu_down);
    
    
    
    if (ray1->n_exit2bn[0].Angle(interaction1->posnu) > PI/2  && !(WHICHPATH==3 || WHICHPATH==4)) 
	beyondhorizon = 1;
    
    
    return;
}//PickDownwardInteractionPoint

void Balloon::calculate_antenna_positions(Settings *settings1, Anita *anita1){
    number_all_antennas = 0;
    Vector antenna_position;
    
    for (int ilayer = 0; ilayer < settings1->NLAYERS; ilayer++){
	for (int ifold = 0; ifold < anita1->NRX_PHI[ilayer]; ifold++){
	    double phi = 0;
	    if (settings1->WHICH==6 || settings1->WHICH==8){ //If payload is either 
		antenna_position = anita1->ANTENNA_POSITION_START[ilayer][ifold];	
	    }else{
		if (settings1->CYLINDRICALSYMMETRY==1){ // for timing code
		    // phi is 0 for antenna 0 (0-31) and antenna 16 (0-31)
		    // antenna 1 (1-32) and antenna 18 (1-32)
		    phi = (double) ifold / (double) anita1->NRX_PHI[ilayer] * 2 * PI + anita1->PHI_OFFSET[ilayer];
		}else{
		    phi = anita1->PHI_EACHLAYER[ilayer][ifold] + anita1->PHI_OFFSET[ilayer];
		}
		antenna_position = Vector(anita1->RRX[ilayer]*cos(phi) + anita1->LAYER_HPOSITION[ilayer]*cos(anita1->LAYER_PHIPOSITION[ilayer]), anita1->RRX[ilayer]*sin(phi)+anita1->LAYER_HPOSITION[ilayer]*sin(anita1->LAYER_PHIPOSITION[ilayer]), anita1->LAYER_VPOSITION[ilayer]);
		// antenna_position is pointing in the +x direction for
		// antenna 0 (0-31) and antenna 16 (0-31)
	    }
	    
	    antenna_position = antenna_position.RotateZ(phi_spin);
	    antenna_position = antenna_position.ChangeCoord(n_north, -1 * n_east);
	    
	    antenna_positions[number_all_antennas] = antenna_position;
	    number_all_antennas++;
	}
    }
    return;
}



void Balloon::GetSlacPositions(Anita *anita1) {
    
    sslacpositions[0]="phi 4-5";
    slacpositions[0]=Vector(0.3,179.7,223.8);
    
    sslacpositions[1]="phi 1, start location";
    slacpositions[1]=Vector(3.1,158.8,222.3);
    
    sslacpositions[2]="phi 13, start location";
    slacpositions[2]=Vector(3.5,158.4,226.2);
    
    sslacpositions[3]="phi 9, start location";
    slacpositions[3]=Vector(8.9,158.2,225.6);
    
    sslacpositions[4]="phi 5, start location";
    slacpositions[4]=Vector(6.0,157.2,219.6);
    
    sslacpositions[5]="phi 1, X=0, Y=2m";
    slacpositions[5]=Vector(1.3,233.8,221.0);
    
    sslacpositions[6]="phi 13 X=0, Y=0";
    slacpositions[6]=Vector(5.4,154.4,226.2);
    
    sslacpositions[7]="phi 13, X=3m, Y=0";
    slacpositions[7]=Vector(107.3,150.5,226.7);
    
    sslacpositions[8]="phi 5, X=3m, Y=0";
    slacpositions[8]=Vector(110.4,144.6,220.1);
    
    sslacpositions[9]="phi 1, left ant. at X=0";
    slacpositions[9]=Vector(104.7,147.8,223.2);
    
    sslacpositions[10]="phi 1, Y=2m, left ant at X=0";
    slacpositions[10]=Vector(104.8,232.4,220.8);
    
    sslacpositions[11]="phi 13, Y=2m, left ant at X=0";
    slacpositions[11]=Vector(104.1,233.1,225.5);
    
    sslacpositions[12]="phi 1, Y=0, right ant at X=0";
    slacpositions[12]=Vector(-102.5,155.8,222.5);
    
    sslacpositions[13]="phi 13, right ant at X=0";
    slacpositions[13]=Vector(-101.6,157.0,230.4);
    
    sslacpositions[14]="phi 13, X=0, Y=-2m";
    slacpositions[14]=Vector(1.7,81.3,230.6);
    
    sslacpositions[15]="phi 13, X=3m, Y=-2m";
    slacpositions[15]=Vector(103.1,82.4,231.2);
    
    sslacpositions[16]="phi 3, X=0, Y=0";
    slacpositions[16]=Vector(-0.5,156.1,225.8);
    
    sslacpositions[17]="phi 7, X=0, Y=0";
    slacpositions[17]=Vector(3.3,157.3,226.7);
    
    sslacpositions[18]="phi 11, X=0, Y=0";
    slacpositions[18]=Vector(2.2,159.2,230.0);
    
    sslacpositions[19]="phi 15, X=0, Y=0";
    slacpositions[19]=Vector(-1.8,158.8,228.5);
    
    sslacpositions[20]="phi 13, X=0, Y=2m";
    slacpositions[20]=Vector(-2.5,237.0,229.7);
    
    sslacpositions[21]="phi 13, X=-3m, Y=2m";
    slacpositions[21]=Vector(-109.8,237.4,229.2);
    
    sslacpositions[22]="phi 13, X=-3m, Y=-2m";
    slacpositions[22]=Vector(-101.8,88.7,231.4);
    
    sslacpositions[22]="phi 14, X=0m, Y=0m";
    slacpositions[22]=Vector(-0.8,160.4,228.2);
    
    sslacpositions[23]="phi 2, X=0m, Y=0m";
    slacpositions[23]=Vector(2.4,158.9,225.7);
    
    sslacpositions[24]="phi 6, X=0m, Y=0m";
    slacpositions[24]=Vector(5.5,158.2,223.6);
    
    sslacpositions[25]="phi 10, X=0m, Y=0m";
    slacpositions[25]=Vector(5.9,160.3,228.2);
    
    sslacpositions[26]="phi 12, X=0m, Y=0m";
    slacpositions[26]=Vector(3.1,160.5,230.4);
    
    sslacpositions[27]="phi 8, X=0m, Y=0m";
    slacpositions[27]=Vector(6.7,159.7,226.7);
    
    sslacpositions[28]="phi 4, X=0m, Y=0m";
    slacpositions[28]=Vector(4.5,157.8,224.1);
    
    sslacpositions[29]="phi 16, X=0m, Y=0m";
    slacpositions[29]=Vector(1.4,159.4,227.3);
    
    sslacpositions[30]="phi 13, X=0m, Y=-6m";
    slacpositions[30]=Vector(3.5,21.1,231.1);
    
    // if you add any positions, be sure the change the number in the next line!!
    for (int i=0;i<=30;i++) {
	
	slacpositions[i]=Vector(slacpositions[i].GetX()*CMINCH/100.,
				slacpositions[i].GetY()*CMINCH/100.,
				slacpositions[i].GetZ()*CMINCH/100.);// inches to meters
	slacpositions[i].SetY(slacpositions[i].GetY()-BN_ALTITUDE*12.*CMINCH/100.-0.5584-anita1->LAYER_VPOSITION[2]+0.4); // subtract balloon altitude, add depth of beam relative to surface.  Use 0.4 m for the height of an antenna.
	
    }
    
}

void GetBoresights(Settings *settings1,Anita *anita1,Position r_bn,double phi_spin,Position r_boresights[Anita::NLAYERS_MAX][Anita::NPHI_MAX],Vector n_north,Vector n_east) {
    
    // this fills r_boresights with the position of the antenna boresights.
    for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
	for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
            
	    Vector n_boresight(1,0,0); // this will be a unit vector that points from the center axis of the payload to the boresight on each level of the payload (it will stay in the x-y plane)
            
	    double phi;
	    // rotate about z axis for phi
	    if (settings1->CYLINDRICALSYMMETRY==1) {
		phi=(double)ifold/(double)anita1->NRX_PHI[ilayer]*2*PI + anita1->PHI_OFFSET[ilayer]+phi_spin;
	    }
	    else
		phi=anita1->PHI_EACHLAYER[ilayer][ifold] + anita1->PHI_OFFSET[ilayer]+phi_spin;
	    
	    
	    n_boresight = n_boresight.RotateZ(phi);
	    
	    // start with a vector that points in the +z direction but with the same length as r_bn
	    r_boresights[ilayer][ifold].SetX(0.); // this will be a unit vector that points from the center axis of the payload to the boresight on each level of the payload (it will stay in the x-y plane)
	    r_boresights[ilayer][ifold].SetY(0.);
	    r_boresights[ilayer][ifold].SetZ(r_bn.Mag());
            
	    Vector zaxis(0.,0.,1.);
	    Vector xaxis(1.,0.,0.);
	    Vector yaxis(0.,1.,0.);
	    
	    // Add the positions of the antennas relative to the center of the payload
	    r_boresights[ilayer][ifold]+=anita1->RRX[ilayer]*n_boresight
	    + anita1->LAYER_VPOSITION[ilayer]*zaxis
	    + anita1->LAYER_HPOSITION[ilayer]*cos(anita1->LAYER_PHIPOSITION[ilayer])*xaxis
	    + anita1->LAYER_HPOSITION[ilayer]*sin(anita1->LAYER_PHIPOSITION[ilayer])*yaxis;
            
	    
	    //   now rotate to balloon's position on the continent
	    r_boresights[ilayer][ifold] = r_boresights[ilayer][ifold].ChangeCoord(n_north,-1.*n_east);
	    
	    
	}
    }
}



void Balloon::GetBoresights(Settings *settings1,Anita *anita1) {
    
    // this fills r_boresights with the position of the antenna boresights.
    for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
	for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
            
	    Vector n_boresight(1,0,0); // this will be a unit vector that points from the center axis of the payload to the boresight on each level of the payload (it will stay in the x-y plane)
            
	    double phi;
	    // rotate about z axis for phi
	    if (settings1->CYLINDRICALSYMMETRY==1) {
		phi=(double)ifold/(double)anita1->NRX_PHI[ilayer]*2*PI + anita1->PHI_OFFSET[ilayer]+phi_spin;
	    }
	    else
		phi=anita1->PHI_EACHLAYER[ilayer][ifold] + anita1->PHI_OFFSET[ilayer]+phi_spin;
	    
	    
	    n_boresight = n_boresight.RotateZ(phi);
	    
	    // start with a vector that points in the +z direction but with the same length as r_bn
	    r_boresights[ilayer][ifold].SetX(0.); // this will be a unit vector that points from the center axis of the payload to the boresight on each level of the payload (it will stay in the x-y plane)
	    r_boresights[ilayer][ifold].SetY(0.);
	    r_boresights[ilayer][ifold].SetZ(r_bn.Mag());
            
	    Vector zaxis(0.,0.,1.);
	    Vector xaxis(1.,0.,0.);
	    Vector yaxis(0.,1.,0.);
	    
	    // Add the positions of the antennas relative to the center of the payload
	    r_boresights[ilayer][ifold]+=anita1->RRX[ilayer]*n_boresight
	    + anita1->LAYER_VPOSITION[ilayer]*zaxis
	    + anita1->LAYER_HPOSITION[ilayer]*cos(anita1->LAYER_PHIPOSITION[ilayer])*xaxis
	    + anita1->LAYER_HPOSITION[ilayer]*sin(anita1->LAYER_PHIPOSITION[ilayer])*yaxis;
            
	    
	    //   now rotate to balloon's position on the continent
	    r_boresights[ilayer][ifold] = r_boresights[ilayer][ifold].ChangeCoord(n_north,-1.*n_east);
	    
	    
	}
    }
}
