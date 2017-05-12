#include "Constants.h"
#include "TRandom3.h"
#include "Settings.h"
#include "earthmodel.hh"
#include "icemodel.hh"
#include <cmath>
#include "Tools.h"
#include "vector.hh"
#include "position.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include "signal.hh"
#include "Primaries.h"
#include "secondaries.hh"


using std::cout;
using std::endl;
using std::ios;
using std::fstream;


const double EarthModel::COASTLINE(30.);
const double EarthModel::MAXTHETA(180.);
const int EarthModel::ILAT_COASTLINE((int)((COASTLINE/MAXTHETA)*(double)NLAT+0.00001)); // corresponding latitude bin to "coastline"
const double EarthModel::GEOID_MAX(6.378137E6); // parameters of geoid model
const double EarthModel::GEOID_MIN(6.356752E6); // from Geodetic Reference System 1980, Bulletin Geodesique, Vol 54:395,1980. // The previous reference gave issue number 3 instead of page number 395


EarthModel::EarthModel(int model,int WEIGHTABSORPTION_SETTING) {
    
    radii[0]=1.2e13;
    radii[1]=(EarthModel::R_EARTH-4.0E4)*(EarthModel::R_EARTH-4.0E4);
    radii[2]=(EarthModel::R_EARTH*EarthModel::R_EARTH); // average radii of boundaries between earth layers
    
    
    //  cout << "In EarthModel, model is " << model << "\n";
    weightabsorption= WEIGHTABSORPTION_SETTING;
    
    CONSTANTICETHICKNESS = (int) (model / 1000);
    model -= CONSTANTICETHICKNESS * 1000;
    
    CONSTANTCRUST = (int) (model / 100);
    model -= CONSTANTCRUST * 100;
    
    FIXEDELEVATION = (int) (model / 10);
    model -= FIXEDELEVATION * 10;
    
    EARTH_MODEL = model;
    //cout<<"CONSTICETHK = "<<CONSTANTICETHICKNESS<<", CNSTCRST = "<<CONSTANTCRUST<<", FIXDELV = "<<FIXEDELEVATION<<", EARTH_MODEL = "<<EARTH_MODEL<<endl;
    for (int i=0;i<NLON;i++) {
	
	Tools::Zero(elevationarray[i],NLAT);
	Tools::Zero(waterthkarray[i],NLAT);
	Tools::Zero(icethkarray[i],NLAT);
	Tools::Zero(softsedthkarray[i],NLAT);
	Tools::Zero(hardsedthkarray[i],NLAT);
	Tools::Zero(uppercrustthkarray[i],NLAT);
	Tools::Zero(middlecrustthkarray[i],NLAT);
	Tools::Zero(lowercrustthkarray[i],NLAT);
	Tools::Zero(crustthkarray[i],NLAT);
	
	
	Tools::Zero(surfacer[i],NLAT);
	Tools::Zero(icer[i],NLAT);
	Tools::Zero(waterr[i],NLAT);
	Tools::Zero(softsedr[i],NLAT);
	Tools::Zero(hardsedr[i],NLAT);
	Tools::Zero(uppercrustr[i],NLAT);
	Tools::Zero(middlecrustr[i],NLAT);
	Tools::Zero(lowercrustr[i],NLAT);
	
	Tools::Zero(waterdensityarray[i],NLAT);
	Tools::Zero(icedensityarray[i],NLAT);
	Tools::Zero(softseddensityarray[i],NLAT);
	Tools::Zero(hardseddensityarray[i],NLAT);
	Tools::Zero(uppercrustdensityarray[i],NLAT);
	Tools::Zero(middlecrustdensityarray[i],NLAT);
	Tools::Zero(lowercrustdensityarray[i],NLAT);
        
    } //Zero Earth model arrays
    
    // see monte carlo note #17
    for (int i=0;i<NLAT;i++) {
	geoid[i]=GEOID_MIN*GEOID_MAX/sqrt(pow(GEOID_MIN,2.)-(pow(GEOID_MIN,2.)-pow(GEOID_MAX,2.))*pow(cos(dGetTheta(i)),2.));
    } //for
    
    
    // Crust 2.0 is binned in 2deg x 2deg bins, area of bin depends on latitude.
    // calculating surface area of bins
    phistep=2*PI/(double)NLON;
    thetastep=(MAXTHETA*RADDEG)/NLAT;
    for (int i=0;i<NLAT;i++) {
	area[i]=phistep*(cos(dGetTheta(i))-cos(dGetTheta(i+1)))*pow(geoid[i],2.);
    } //for
    
    
    if (EARTH_MODEL == 0)
	ReadCrust(crust20_in);
    else {
	cout<<"Error!  Unknown Earth model requested!  Defaulting to Crust 2.0 model.\n";
	ReadCrust(crust20_in);
    } //else
    
} //EarthModel constructor (int mode)

EarthModel::~EarthModel() {} //EarthModel destructor - no dynamic variables, nothing to delete


double EarthModel::LongtoPhi_0isPrimeMeridian(double longitude) {
    
    double phi;
    // convert longitude (-180 to 180) to phi (0 to 2pi) wrt +x
    // in radians
    phi=(90-longitude); 
    if (phi<0.)
	phi+=360.;
    
    phi=phi*RADDEG;
    
    return phi;
}
double EarthModel::LongtoPhi_0is180thMeridian(double longitude) {
    
    double phi;
    // convert longitude (0 to 360) to phi (0 to 2pi) wrt +x
    
    phi=(270.-longitude); 
    if (phi<0.)
	phi+=360.;
    
    phi=phi*RADDEG;
    // in radians
    
    return phi;
}

double EarthModel::Geoid(double latitude) {
    // latitude here is 0 at the south pole and 180 at the north pole
    
    return (GEOID_MIN*GEOID_MAX/
	    sqrt(GEOID_MIN*GEOID_MIN-(GEOID_MIN*GEOID_MIN-GEOID_MAX*GEOID_MAX)*
		 cos(latitude*RADDEG)*cos(latitude*RADDEG)));
} //Geoid(lat)

double EarthModel::Geoid(const Position &pos) {
    return Geoid(pos.Lat());
} //Geoid(Position)

double EarthModel::IceThickness(double lon,double lat) {
    return icethkarray[(int)(lon/2)][(int)(lat/2)]*1000.;
} //IceThickness(lon,lat)

double EarthModel::IceThickness(const Position& pos) {
    return IceThickness(pos.Lon(),pos.Lat());
} //IceThickness(Position)
int EarthModel::InFirn(const Position& pos) {
    if (pos.Mag()-Surface(pos)<FIRNDEPTH)
	return 0;
    return 1;
} //InFirn(Position)
double EarthModel::SurfaceDeepIce(const Position& pos) { // surface of the deep ice (where you reach the firn)
    return  surfacer[(int)(pos.Lon()/2)][(int)(pos.Lat()/2)] + geoid[(int)(pos.Lat()/2)] + FIRNDEPTH;
} //Surface(lon,lat)

double EarthModel::Surface(double lon,double lat) {
    return surfacer[(int)(lon/2)][(int)(lat/2)] + geoid[(int)(lat/2)];
} //Surface(lon,lat)

double EarthModel::Surface(const Position& pos) {
    return surfacer[(int)(pos.Lon()/2)][(int)(pos.Lat()/2)] + geoid[(int)(pos.Lat()/2)];
} //Surface(Position)

double EarthModel::RockSurface(double lon,double lat) {
    return (Surface(lon,lat) - IceThickness(lon,lat) - WaterDepth(lon,lat));
} //RockSurface(lon,lat)

double EarthModel::RockSurface(const Position& pos) {
    return RockSurface(pos.Lon(),pos.Lat());
} //RockSurface(lon,lat)

double EarthModel::SurfaceAboveGeoid(double lon,double lat) {
    return surfacer[(int)(lon/2)][(int)(lat/2)];
} //SurfaceAboveGeoid(lon,lat)

double EarthModel::SurfaceAboveGeoid(const Position& pos) {
    return surfacer[(int)(pos.Lon()/2)][(int)(pos.Lat()/2)];
} //SurfaceAboveGeoid(Position)

double EarthModel::WaterDepth(double lon,double lat) {
    return waterthkarray[(int)(lon/2)][(int)(lat/2)]*1000;
} //WaterDepth(lon,lat)

double EarthModel::WaterDepth(const Position& pos) {
    return WaterDepth(pos.Lon(),pos.Lat());
} //WaterDepth(Position)

double EarthModel::GetLat(double theta) {
    return theta*DEGRAD;
} //GetLat

double EarthModel::GetLon(double phi) {
    // input is phi in radians wrt +x
    double phi_deg = phi*DEGRAD; 
    if (phi_deg > 270)   
	phi_deg = phi_deg - 360.; // this puts it from -90 to 270
    
    return (360.*3./4. - phi_deg); // returns 0 to 360 degrees (going from -180 to 180 deg longitude like Crust 2.0 does)
} //GetLon

double EarthModel::GetDensity(double altitude, const Position earth_in, const Position posnu,
			      int& crust_entered, // 1 or 0
			      int& mantle_entered, // 1 or 0
			      int& core_entered){
  
        Position where = earth_in;
	
	double lon = where.Lon();
	double lat = where.Lat();
	//cout <<"Lon and Lat are "<<lon<<","<<lat<<"\n";

	int ilon = (int)(lon/2);
	int ilat = (int)(lat/2);

	double ddensity =0; //initilize ddensity

	double surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point

	double local_icethickness = this->IceThickness(lon,lat);
	double local_waterdepth = WaterDepth(lon,lat);

	

	if(altitude>surface_elevation+0.1){ // if it is above the surface, it's messed up
	  
	    }
	if(altitude>surface_elevation+0.1){
	  ddensity=1.25;
	  //cout <<"density is air! \n";
	}
	if (altitude<=surface_elevation+0.1 && altitude>(surface_elevation-local_icethickness)) // the 0.1 is just to take care of precision issues.   It could have been 0.01 or 0.001.
		ddensity=icedensityarray[ilon][ilat]*1000;
	  
	else if (altitude<=(surface_elevation-local_icethickness) && altitude>(surface_elevation-local_icethickness-local_waterdepth))
		ddensity=waterdensityarray[ilon][ilat]*1000;
	else if (altitude<=(surface_elevation-local_icethickness-local_waterdepth) && altitude>softsedr[ilon][ilat]) {
		ddensity=softseddensityarray[ilon][ilat]*1000;
	       	crust_entered=1; //Switch that lets us know we've penetrated into the crust
		 } //end if
	else if (altitude<=softsedr[ilon][ilat] && altitude>hardsedr[ilon][ilat])
		ddensity=hardseddensityarray[ilon][ilat]*1000;
	else if (altitude<=hardsedr[ilon][ilat] && altitude>uppercrustr[ilon][ilat])
		ddensity=uppercrustdensityarray[ilon][ilat]*1000;
	else if (altitude<=uppercrustr[ilon][ilat] && altitude>middlecrustr[ilon][ilat])
		ddensity=middlecrustdensityarray[ilon][ilat]*1000;
	else if (altitude<=middlecrustr[ilon][ilat] && altitude>lowercrustr[ilon][ilat])
		ddensity=lowercrustdensityarray[ilon][ilat]*1000;
	else if (altitude<=lowercrustr[ilon][ilat])
		ddensity=densities[1];
	
	    return ddensity;

}//Get Density




Vector EarthModel::GetSurfaceNormal(const Position &r_out)
{
    Vector n_surf = r_out.Unit();
    if (FLATSURFACE)
	return n_surf;
    
    double theta=r_out.Theta();
    
    int ilon,ilat;
    GetILonILat(r_out,ilon,ilat);
    
    int ilon_previous=ilon-1;
    if (ilon_previous<0)
	ilon_previous=NLON-1;
    
    int ilon_next=ilon+1;
    if (ilon_next==NLON)
	ilon_next=0;
    
    double r=(geoid[ilat]+surfacer[ilon][ilat])*sin(theta);
    
    double slope_phi=(surfacer[ilon_next][ilat]-surfacer[ilon_previous][ilat])/(r*2*phistep);
    
    int ilat_previous=ilat-1;
    if (ilat_previous<0)
	ilat_previous=0;
    
    int ilat_next=ilat+1;
    if (ilat_next==NLAT)
	ilat_next=NLAT-1;
    
    double slope_costheta=(surfacer[ilon][ilat_next]-surfacer[ilon][ilat_previous])/((geoid[ilat]+surfacer[ilon][ilat])*2*thetastep);
    
    // first rotate n_surf according to tilt in costheta and position on continent - rotate around the y axis.
    double angle=atan(slope_costheta);
    
    n_surf = n_surf.RotateY(angle);
    
    // now rotate n_surf according to tilt in phi - rotate around the z axis.
    angle=atan(slope_phi);
    
    n_surf = n_surf.RotateZ(angle);
    
    return n_surf;
    
} //method GetSurfaceNormal

double EarthModel::SmearPhi(int ilon)  {
    
    
    double phi=((double)(360.*3./4.-((double)ilon+gRandom->Rndm())*360/180))*RADDEG;
    if (phi<0 && phi>-1*PI/2)
	phi+=2*PI;
    
    
    return phi;
} //SmearPhi

double EarthModel::SmearTheta(int ilat)  {
    
    // remember that we should smear it evenly in cos(theta).
    // first get the cos(theta)'s at the boundaries.
    
    double theta1=dGetTheta(ilat)-PI/(double)NLAT/2.;
    double theta2=dGetTheta(ilat+1)-PI/(double)NLAT/2.;
    
    double costheta1=cos(theta1);
    double costheta2=cos(theta2);
    
    
    
    double costheta=gRandom->Rndm()*(costheta2-costheta1)+costheta1;
    
    double theta=acos(costheta);
    
    return theta;
} //SmearTheta

void EarthModel::ReadCrust(string test) {
    
    // reads in altitudes of 7 layers of crust, ice and water
    // puts data in arrays
    
    fstream infile(test.c_str(),ios::in);
    
    string thisline; // for reading in file
    string slon; //longitude as a string
    string slat; // latitude as a string
    string selev; // elevation (km relative to geoid)
    string sdepth; // depth (km)
    string sdensity; // density (g/cm^3)
    double dlon,dlat; // longitude, latitude as double
    int endindex; // index along thisline for parsing
    int beginindex; // same
    
    int indexlon=0; // 180 bins in longitude
    int indexlat=0; // 90 bins in latitude
    
    string layertype; // water, ice, etc.
    
    while(!infile.eof()) {
	getline(infile,thisline,'\n'); 
	
	int loc=thisline.find("type, latitude, longitude,"); 
	
	if (loc!=(int)(string::npos)) {      
	    
	    beginindex=thisline.find_first_not_of(" ",57);
	    
	    endindex=thisline.find_first_of(" ",61);
	    
	    slat=thisline.substr(beginindex,endindex-beginindex);
	    dlat=(double)atof(slat.c_str());
	    
	    beginindex=thisline.find_first_not_of(" ",68);
	    endindex=thisline.find_first_of(" ",72);
	    
	    slon=thisline.substr(beginindex,endindex-beginindex);
	    dlon=(double)atof(slon.c_str());
	    
	    
	    indexlon=(int)((dlon+180)/2);
	    indexlat=(int)((90+dlat)/2);
	    
	    beginindex=thisline.find_first_not_of(" ",78);
	    endindex=thisline.find_first_of(" ",83);
	    
	    selev=thisline.substr(beginindex,endindex-beginindex);
	    elevationarray[indexlon][indexlat]=(double)atof(selev.c_str());
	    
	} //if
	
	for (int i=0;i<4;i++) {
	    getline(infile,thisline,'\n');
	} //for
	
	for (int i=0;i<7;i++) {
	    getline(infile,thisline,'\n');
	    
	    endindex=thisline.length()-1;
	    beginindex=thisline.find_last_of("0123456789",1000);
	    layertype=thisline.substr(beginindex+3,endindex-beginindex);
	    
	    
	    beginindex=thisline.find_first_not_of(" ",0);
	    endindex=thisline.find_first_of(" ",beginindex);
	    
	    sdepth=thisline.substr(beginindex,endindex-beginindex-1);
	    
	    
	    // fills arrays of thicknesses of each layer
	    if (layertype.substr(0,5)=="water") 
		waterthkarray[indexlon][indexlat]=(double)atof(sdepth.c_str()); 
	    if (layertype.substr(0,3)=="ice") 
		icethkarray[indexlon][indexlat]=(double)atof(sdepth.c_str());
	    if (layertype.substr(0,8)=="soft sed") 
		softsedthkarray[indexlon][indexlat]=(double)atof(sdepth.c_str());
	    if (layertype.substr(0,8)=="hard sed") 
		hardsedthkarray[indexlon][indexlat]=(double)atof(sdepth.c_str());
	    if (layertype.substr(0,11)=="upper crust") 
		uppercrustthkarray[indexlon][indexlat]=(double)atof(sdepth.c_str());
	    if (layertype.substr(0,12)=="middle crust") 
		middlecrustthkarray[indexlon][indexlat]=(double)atof(sdepth.c_str());
	    if (layertype.substr(0,11)=="lower crust") 
		lowercrustthkarray[indexlon][indexlat]=(double)atof(sdepth.c_str());
	    
	    //      cout << "indexlon, indexlat, icethkarray " << indexlon << " " << indexlat << " " << icethkarray[indexlon][indexlat] << "\n";
	    
	    // region where Ross Ice Shelf was not accounted for in Crust 2.0
	    // add it in by hand
	    if (indexlat==5 && (indexlon<=5 || indexlon>=176)) // Ross Ice Shelf
		icethkarray[indexlon][indexlat]=0.5;
	    
	    beginindex=thisline.find_first_not_of(" ",endindex);
	    endindex=thisline.find_first_of(" ",beginindex);
	    
	    
	    beginindex=thisline.find_first_not_of(" ",endindex);
	    endindex=thisline.find_first_of(" ",beginindex);
	    
	    beginindex=thisline.find_first_not_of(" ",endindex);
	    endindex=thisline.find_first_of(" ",beginindex);
	    
	    
	    sdensity=thisline.substr(beginindex,endindex-beginindex);
	    
	    double ddensity=(double)atof(sdensity.c_str());
	    
	    
	    // fills arrays of densities of each layer
	    if (layertype.substr(0,5)=="water") 
		waterdensityarray[indexlon][indexlat]=ddensity; 
	    if (layertype.substr(0,3)=="ice") 
		icedensityarray[indexlon][indexlat]=ddensity;
	    if (layertype.substr(0,8)=="soft sed") 
		softseddensityarray[indexlon][indexlat]=ddensity;
	    if (layertype.substr(0,8)=="hard sed") 
		hardseddensityarray[indexlon][indexlat]=ddensity;
	    if (layertype.substr(0,11)=="upper crust") 
		uppercrustdensityarray[indexlon][indexlat]=ddensity;
	    if (layertype.substr(0,12)=="middle crust")
		middlecrustdensityarray[indexlon][indexlat]=ddensity;
	    if (layertype.substr(0,11)=="lower crust") 
		lowercrustdensityarray[indexlon][indexlat]=ddensity;
	} //for (reading all lines for one location given in Crust 2.0 input file)
	
	if (CONSTANTCRUST) {
	    softsedthkarray[indexlon][indexlat]=40.;
	    hardsedthkarray[indexlon][indexlat]=0;
	    uppercrustthkarray[indexlon][indexlat]=0;
	    middlecrustthkarray[indexlon][indexlat]=0;
	    lowercrustthkarray[indexlon][indexlat]=0;
	    crustthkarray[indexlon][indexlat]=0;
	    softseddensityarray[indexlon][indexlat]=2.9;
	} //if (set crust thickness to constant everywhere)
	if (CONSTANTICETHICKNESS) {
	    icethkarray[indexlon][indexlat]=3.;
	    waterthkarray[indexlon][indexlat]=0.;
	} //if (set ice thickness to constant everywhere)
	
	// adds up total thickness of crust
	crustthkarray[indexlon][indexlat]=softsedthkarray[indexlon][indexlat]+
	hardsedthkarray[indexlon][indexlat]+
	uppercrustthkarray[indexlon][indexlat]+
	middlecrustthkarray[indexlon][indexlat]+
	lowercrustthkarray[indexlon][indexlat];
	
	if (indexlon==179 && indexlat==0)
	    break;
    }  // done reading file
    
    for (int i=0;i<NLON;i++) {
	for (int j=0;j<NLAT;j++) {
	    
	    if (FIXEDELEVATION) 	
		elevationarray[i][j]=icethkarray[i][j]*1000;
	    
	    if (waterthkarray[i][j] != 0) 
		surfacer[i][j]=elevationarray[i][j]+waterthkarray[i][j]*1000+icethkarray[i][j]*1000;	  
	    else
		surfacer[i][j]=elevationarray[i][j];
	    
	    if (fabs(surfacer[i][j])<1.E-10)
		surfacer[i][j] = 0;	
	    
	    // reminder- waterr is elevation at *bottom* of water layer, etc. 
	    // in units of m
	    waterr[i][j]=surfacer[i][j]-(icethkarray[i][j]+waterthkarray[i][j])*1000;
	    if ((double)fabs(waterr[i][j])<1.E-10)
		waterr[i][j]=0;
	    icer[i][j]=waterr[i][j]+
	    waterthkarray[i][j]*1000;
	    softsedr[i][j]=waterr[i][j]-
	    softsedthkarray[i][j]*1000;
	    hardsedr[i][j]=waterr[i][j]-
	    (softsedthkarray[i][j]+
	     hardsedthkarray[i][j])*1000;
	    uppercrustr[i][j]=waterr[i][j]-
	    (softsedthkarray[i][j]+
	     hardsedthkarray[i][j]+
	     uppercrustthkarray[i][j])*1000;
	    middlecrustr[i][j]=waterr[i][j]-
	    (softsedthkarray[i][j]+
	     hardsedthkarray[i][j]+
	     uppercrustthkarray[i][j]+
	     middlecrustthkarray[i][j])*1000;
	    lowercrustr[i][j]=waterr[i][j]-
	    (softsedthkarray[i][j]+
	     hardsedthkarray[i][j]+
	     uppercrustthkarray[i][j]+
	     middlecrustthkarray[i][j]+
	     lowercrustthkarray[i][j])*1000;
	} //for (latitude bins)
    } //for (longitude bins)
    
    // calculate ice volume for comparison with expectation
    volume=0;
    ice_area=0.;
    double sumarea=0; // sum of surface area of ice
    average_iceth = 0;
    max_icevol_perbin=0.;
    max_icethk_perbin=0.;
    for (int j=0;j<ILAT_COASTLINE;j++) {
	for (int i=0;i<NLON;i++) {
	    volume+=icethkarray[i][j]*1000.*area[j];
	    if (icethkarray[i][j]>0)
		ice_area+=area[j];
	    if (icethkarray[i][j]*area[j]>max_icevol_perbin)
		max_icevol_perbin=icethkarray[i][j]*area[j];
	    if (icethkarray[i][j]>max_icethk_perbin)
		max_icethk_perbin=icethkarray[i][j];
	    
	    /*
	     // j=6 corresponds to 80deg S
	     if (j==6) {
	     // fill output file, just for plotting
	     outfile << surfacer[i][j] << "\t" << waterr[i][j] << "\t" << icer[i][j] << "\t" << icethkarray[i][j] << " " << waterr[i][j]-icer[i][j] << " " << (waterr[i][j]-icer[i][j])*area[j] << " " << area[j] << " " << volume << "\n";
	     }//ifx
	     */
	    
	    // find average ice thickness
	    average_iceth+=(surfacer[i][j]-icer[i][j])*area[j];
	    sumarea+=area[j];
	} //for
    } //for
    average_iceth=average_iceth/sumarea; 
    
    // find the place where the crust is the deepest.
    // for finding where to start stepping in Getchord
    MIN_ALTITUDE_CRUST=1.E6;
    //MAX_VOL=-1.E6;
    for (int i=0;i<NLON;i++) {
	for (int j=0;j<NLAT;j++) {
	    if (elevationarray[i][j]-(crustthkarray[i][j])*1000<MIN_ALTITUDE_CRUST) {
		if (waterthkarray[i][j]==0)
		    MIN_ALTITUDE_CRUST=elevationarray[i][j]-(crustthkarray[i][j]+icethkarray[i][j])*1000;
		else
		    MIN_ALTITUDE_CRUST=elevationarray[i][j]-crustthkarray[i][j]*1000;
	    }//if
	    //if (icethkarray[i][j]*1000.*area[j]>MAX_VOL) 
	    //MAX_VOL=icethkarray[i][j]*1000.*area[j];      
	}//for
    }//for
    
    //record depth of crust-mantle interface
    radii[1]=(Geoid(0.)+MIN_ALTITUDE_CRUST)*(Geoid(0.)+MIN_ALTITUDE_CRUST);
    
}//ReadCrust


Vector EarthModel::PickPosnuForaLonLat(double lon,double lat,double theta,double phi) {
    
    
    double surface_above_geoid = this->SurfaceAboveGeoid(lon,lat);
    double local_ice_thickness = this->IceThickness(lon,lat);
    
    double rnd3=gRandom->Rndm();
    
    
    
    double elevation = surface_above_geoid - rnd3*local_ice_thickness; // elevation of interaction
    
    
    
    //cout << "Inside PickInteractionLocation, lon, lat, local_ice_thickness are " << lon << " " << lat << " " << local_ice_thickness << "\n";
    
    if (elevation > surface_above_geoid || elevation < (surface_above_geoid - local_ice_thickness))
	cout<<"elevation > surface_above_geoid || elevation < (surface_above_geoid - local_ice_thickness)\n";
    
    Vector posnu((elevation+this->Geoid(lat))*sin(theta)*cos(phi),(elevation+this->Geoid(lat))*sin(theta)*sin(phi),(elevation+this->Geoid(lat))*cos(theta));
    
    if (((this->Geoid(lat) + surface_above_geoid)) - posnu.Mag() < 0) {
	cout<<"/nYikes!  (Geoid(lat) + Surface(lon,lat)) - sqrt(dSquare(posnu) = "<<((this->Geoid(lat) + surface_above_geoid)) - posnu.Mag()<<"/nlon, lat: "<<lon<<" , "<<lat<< " " << endl<<endl;
	cout << "local_ice_thickness is " << local_ice_thickness << "\n";
    }
    
    return posnu;
}

double EarthModel::dGetTheta(int ilat) {
    return (((double)ilat+0.5)/(double)NLAT*MAXTHETA)*RADDEG;
} //dGetTheta(int)

double EarthModel::dGetPhi(int ilon) {
    // this takes as an input the crust 2.0 index 0=-180 deg longitude to 179=+180 deg longitude
    // its output is phi in radians
    // from ~ -pi/2 to 3*pi/2 
    return (double)(-1*((double)ilon+0.5)+(double)NLON)*2*PI/(double)NLON-PI/2;
} //dGetPhi(int)

void EarthModel::GetILonILat(const Position &p,int& ilon,int& ilat) {
    // Phi function outputs from 0 to 2*pi wrt +x
    double phi_deg=p.Phi()*DEGRAD;
    
    if (phi_deg>270)
	phi_deg=phi_deg-360;
    // now it's from -90 to 270
    
    ilon=(int)((360.*3./4.-phi_deg)*180./360.); // ilon is from 0 (at -180 longitude) to 180 (at 180 longitude)
    
    ilat=(int)((p.Theta()*DEGRAD)/2.);
    
} //method GetILonILat
void EarthModel::EarthCurvature(double *array,double depth_temp) {
    
    Position parray;
    parray.SetXYZ(array[0],array[1],array[2]);
    
    // adjust array coordinates so that it fits to a curved earth surface at a specific depth 
    double length=Surface(parray)-depth_temp; // length=distance from center of earth
    
    double rxposx=array[0]; // x coordinate of antenna position
    double rxposy=array[1]; // y coordinate of antenna position
    double rxdr=sqrt(rxposx*rxposx+rxposy*rxposy); // distance in horizontal plane from the center of the detector to the antenna
    if (Tools::dSquare(array)==0) cout << "Attempt of divide by zero in Earth curvature!!\n";
    double rxdtheta=asin(rxdr/sqrt(Tools::dSquare(array)));
    double rxdphi=atan2(rxposy,rxposx);
    
    array[0]=length*sin(rxdtheta)*cos(rxdphi);// have the array sit on a sphere of radius "length"
    array[1]=length*sin(rxdtheta)*sin(rxdphi);
    array[2]=length*cos(rxdtheta);
    
}
Position EarthModel::WhereDoesItEnter(const Position &posnu,const Vector &nnu) {
    // now get neutrino entry point...
    double p = posnu.Mag(); // radius of interaction
    double costheta = (nnu*posnu) / p; // theta of neutrino at interaction position
    double sintheta = sqrt(1-costheta*costheta);
    
    double lon = posnu.Lon();
    double lat = posnu.Lat();
    
    double a=0; // length of chord
    
    double R = Surface(lon,lat);
    double delta = R - p; // depth of the interaction
    // if interaction occurs below surface, as it should
    
    if (delta>-0.001) {
	a=p*costheta+sqrt(R*R*costheta*costheta+2*delta*R*sintheta*sintheta); // chord length
	if (a<0) {
	    cout << "Negative chord length: " << a << "\n";
	} //end if
    } //end if (interaction below surface)  
    else if (delta<=-0.001) {
	
	cout << "lon, lat from WhereDoesItEnter is " << lon << " " << lat << "\n";
	cout << "geoid, surface, p, surface-p are " << Geoid(lat) << " " << Surface(lon,lat) << " " << p << " , "<<(Surface(lon,lat)-p)<<"\n";
	
    } //else if: error: interaction takes place above the surface


    // first approx
    // find where nnu intersects spherical earth surface
    Position r_in = posnu - a*nnu;

    int iter = 0;
    // now do correction 3 times
    //for (iter=0; iter<3; iter++) {
    //    delta = r_in.Mag() - antarctica1->Surface( r_in );
    //    r_in = r_in + (delta * nnu);
    //}
    
    // how far away r_in is from surface, along line to center of earth
    delta = r_in.Mag() - Surface( r_in );
    
    while ( fabs(delta) >= 0.1 ) {  //if it's greater than 10 cm
      r_in = r_in + (delta * nnu); // move distance delta along nnu  for the next iteration
      delta = r_in.Mag() - Surface( r_in ); // find new delta
        iter++;
        if ( iter > 10 ) { // stop after 10 iterations and just revert to simple answer
            //cout<<"\n r_in iteration more than 10 times!!! delta : "<<delta<<". now set r_in as before."<<endl;
            r_in = Surface( r_in ) * r_in.Unit();   // the way before
            delta = r_in.Mag() - Surface( r_in );
        }
    }

    
    //lon = r_in.Lon();
    //lat = r_in.Lat();
    
    //r_in = antarctica1->Surface(lon,lat) * r_in.Unit();
    
    
    return r_in;
} //method WhereDoesItEnter
