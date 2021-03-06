
#include "TRandom3.h"
#include "Constants.h"
#include "vector.hh"
#include "position.hh"
#include "TF1.h"
#include <fstream>
#include <iostream>
#include "Settings.h"
#include "earthmodel.hh"
#include "icemodel.hh"
#include "Primaries.h"
#include "Settings.h"
#include "counting.hh"

#include <cmath>



#include "TH2D.h"
#include "TCanvas.h"

Primaries::Primaries(){//constructor

  // This is for parametrizations in Connolly et al. 2011  
  //in the form of [i][j] where i is neutrino type(nu_nubar) and j is current type, "nc" vs "cc".
  //[nu_nubar][currentint]
  //[0=nu, 1=nubar][0=neutral current, 1=charged current]
  //[0][0]->[nu][neutral current]
  //[0][1]->[nu][charged current]
  //[1][0]->[nubar][neutral current]
  //[1][1]->[nubar][charged current]
  
  //[nu][neutral current]
  c0[0][0]=-1.826;
  c1[0][0]=-17.31;
  c2[0][0]=-6.448; 
  c3[0][0]=1.431;
  c4[0][0]=-18.61;
  
  //[nu][charged current]
  c0[0][1]=-1.826;
  c1[0][1]=-17.31;
  c2[0][1]=-6.406; 
  c3[0][1]=1.431;
  c4[0][1]=-17.91;
  
  //[nubar][neutral current]	
  c0[1][0]=-1.033;
  c1[1][0]=-15.95;
  c2[1][0]= -7.296; 
  c3[1][0]=1.569;
  c4[1][0]=-18.30;
  
  //[nubar][charged current]
  c0[1][1]=-1.033;
  c1[1][1]=-15.95;
  c2[1][1]=-7.247; 
  c3[1][1]=1.569;
  c4[1][1]=-17.72;
  
  char ch[50];
  string stmp;
  string sbase="fsigma";
  for(int i=0; i<=1;i++){ // nu, nubar
    for(int j=0; j<=1; j++){ // nc, cc
      sprintf(ch,"%d%d",i,j);
      stmp=ch;	
      m_fsigma[i][j]=new TF1((sbase+stmp).c_str(),"pow(10, [1]+[2]*log(x-[0])+[3]*pow(log(x-[0]),2)+[4]/log(x-[0]))", 4., 21.);//check bounds. they're in log10 GeV.
      //x=log10(pnu/GeV).
      m_fsigma[i][j]->SetParameters(c0[i][j], c1[i][j], c2[i][j], c3[i][j], c4[i][j]);
      //"fsigma00"->[nu][neutral current]
      //"fsigma01"->[nu][charged current]
      //"fsigma10"->[nubar][neutral current]
      //"fsigma11"->[nubar][charged current]
    }		
  }
  m_csigma=new TCanvas("m_csigma","m_csigma title",1000, 700);
  m_hsigma=new TH2D("hsigma","title hsigma", 600, 7., 12., 600, -40., -30.);
  
  m_hsigma->SetTitle("log10 (pnu) vs.log10 Cross Section Sigma");
  m_hsigma->GetXaxis()->SetTitle("Log10(Ev/ GeV)");
  m_hsigma->GetYaxis()->SetTitle("log10(Cross Section/ m^2)");
  
  m_hsigma->Draw("scat");
  m_hsigma->SetMarkerStyle(7);
  m_hsigma->SetMarkerSize(3);
  
  // again y distributions from Connolly et al. 2011
  m_myY=new Y();
  
  //Low y/////////////////////////////From Table V. Connolly Calc 2011.
  //A_low[4];//same for any [i]nu_nubar and [j]currentint.
  A_low[0]=0.0;
  A_low[1]=0.0941;
  A_low[2]=4.72;
  A_low[3]=0.456;
  //high y////////////////////////////
  //[0=nu, 1=nubar][0=neutral current, 1=charged current]
  //[nu_bar][currentint];
  A0_high[0][0]=-0.005;
  A1_high[0][0]=0.23;
  A2_high[0][0]=3.0;
  A3_high[0][0]=1.7;
  
  A0_high[0][1]=-0.008;
  A1_high[0][1]=0.26;
  A2_high[0][1]=3.0;
  A3_high[0][1]=1.7;
  
  A0_high[1][0]=-0.005;
  A1_high[1][0]=0.23;
  A2_high[1][0]=3.0;
  A3_high[1][0]=1.7;
  
  A0_high[1][1]=-0.0026;
  A1_high[1][1]=0.085;
  A2_high[1][1]=4.1;
  A3_high[1][1]=1.7;
  
  b0=2.55;
  b1=-0.0949; //C2_low=b0+b1*epsilon;
  
  ymin_low=2.E-5;
  ymax_low=1.E-3;
  ymin_high=ymax_low;
  ymax_high=1.;
  
  run_old_code=0;//for GetSigma() & Gety() runs of the old code if 1, else runs current code.

  // 0=Reno
  // 1=Connolly et al. 2011
  mine[0]=1.2E15;
  mine[1]=1.E4;// minimum energy for cross section parametrizations
  maxe[0]=1.E21;
  maxe[1]=1.E21; // use the same upper limit for reno as for connolly et al.


}

Primaries::~Primaries(){//default deconstructor
  m_hsigma->Draw("same");
  m_csigma->Print("sigmaCrossSection.pdf");
  delete m_hsigma;
  delete m_myY;
  for(int i=0; i<=1;i++){ // nu, nubar
    for(int j=0; j<=1; j++){ // nc, cc
      delete m_fsigma[i][j];
    }
  }



}//deconstructor

int Primaries::GetSigma(double pnu,double& sigma,double &len_int_kgm2,Settings *settings1,int nu_nubar,int currentint){
  // calculate cross section
  if (pnu<mine[settings1->SIGMAPARAM] || pnu>maxe[settings1->SIGMAPARAM]) {
    cout <<  "Need a parameterization for this energy region.\n";
    return 0;
  } //if
  else {
   
    nu_nubar=1;//default.
    //nu=0, nubar=1
    if(nu_nubar!=0 && nu_nubar!=1){   
      cout<<"nu_nubar is not defined correctly!\n";
      return 0;
    }
    if (currentint!=0 && currentint!=1){//default "cc"
      cout<<"Current is not cc or nc!\n";
      return 0;
    }
    
    if(settings1->SIGMAPARAM==0){ // Reno
      // fit to cross sections calculated by M.H. Reno using the same method as Gandhi et al, but with the CTEQ6-DIS parton distribution functions instead of the CTEQ4-DIS distribution functions
      sigma=(2.501E-39)*pow(pnu/1.E9,0.3076)*settings1->SIGMA_FACTOR; // 10^18 eV - 10^21 eV(use this one for ANITA)
      //sigma=(1.2873E-39)*pow(pnu/1.E9,0.33646)*SIGMA_FACTOR; // 10^17 eV - 10^20 eV (use this one for SalSA)
    }//old code
    else if (settings1->SIGMAPARAM==1) {//Connolly et al.
      double pnuGeV=pnu/1.E9;//Convert eV to GeV.
      double epsilon=log10(pnuGeV);
      sigma=settings1->SIGMA_FACTOR*(m_fsigma[nu_nubar][currentint]->Eval(epsilon))/1.E4;//convert cm to meters. multiply by (1m^2/10^4 cm^2).
      
      
      if(m_hsigma->GetEntries()<2000){
	m_hsigma->Fill(epsilon, log10(sigma));
      }
    }//else current code
  }//if
  // interaction length in kg/m^2
  
  len_int_kgm2=M_NUCL/sigma; // kg/m^2
  return 1;
} //GetSigma


// The interaction
double Primaries::Gety(Settings *settings1,double pnu,int nu_nubar,int currentint) {
  // THIS IS A ROUGH PARAMETRIZATION OF PLOT 6 FROM 
  //  Ghandhi,Reno,Quigg,Sarcevic  hep-ph/9512364
  //  (the curves are not in their later article.)
  //  There is also a slow energy dependence.
   //cout << "I'm here.\n";
 
  if(settings1->YPARAM==0){
  	double rnd;
  	double x = 0;
  	const double R1=0.36787944;  // 1/e
	const double R2=0.63212056;  // 1-r1
  
  	// generate according to Ghandi fig. 6 
  	// adjust exponent until looks like the curve
  	//  and has right mean.
  	//  (Note this is not the fcn, but the inverse of the integral...)
  


  	rnd = gRandom->Rndm(1); // (0,1)
  	//  cout << "R1, R2, rnd are " << R1 << " " << R2 << " " << rnd << "\n";

	x=pow(-log(R1+rnd*R2),2.5); 


  	return x;   
  }//old Gety

  else { //use prescription in Connolly et al.2011
  	nu_nubar=0;
	double pnuGeV=pnu/1.E9;
	double elast_y=m_myY->pickY(nu_nubar,currentint,pnuGeV);
	return elast_y;   
  }//current Gety

} //Gety

double Primaries::Getyweight(double pnu, double y, int nu_nubar, int currentint){
	//from Connolly Calc 2011, Equations 9, 10, 11, 16, and 17.
	double dy=0.;//default
	//Ev, cc or nc, nu or nubar.
	
	double C0_highbar, C0_lowbar,C0_high, C0_low;//these C0's are normalization factors.
	double dNdy=0.;//default
	double U, W, B, T;//are added in to help with readability of equations.
	double C1_low, C2_low, C1_high, C1_highbar;
	double weighty;
	double epsilon=log10(pnu/1.E9);

	C2_low=b0+b1*epsilon;//Eq(17)
	C1_low=A_low[0]+A_low[1]*(-exp(-(epsilon-A_low[2])/A_low[3]));//Eq(16)
	nu_nubar=0;
	C1_high=A0_high[nu_nubar][currentint]+A1_high[nu_nubar][currentint]*(-exp(-(epsilon-A2_high[nu_nubar][currentint])/A3_high[nu_nubar][currentint]));

	nu_nubar=1;
	C1_highbar=A0_high[nu_nubar][currentint]+A1_high[nu_nubar][currentint]*(-exp(-(epsilon-A2_high[nu_nubar][currentint])/A3_high[nu_nubar][currentint]));//Eq(16)
	if(nu_nubar==0) {
		U=1-1/C2_low;
		W=fabs( (ymax_high-C1_high)/(ymin_high-C1_high));
		B=(pow(ymax_low-C1_low, 1/C2_low)/(ymax_low-C1_high));
		T=B*((pow(ymax_low-C1_low, U)-pow(ymin_low-C1_low, U) )/U)+log(W);
		C0_high=1/T;	
		C0_low=C0_high*(pow(ymax_low-C1_low, 1/C2_low))/(ymax_low-C1_high);
		if(y<ymax_low){//Eq(9)
				dy=0.00002;
				dNdy=C0_low/pow(y-C1_low, 1/C2_low);//Eq(10)
		}
		else if(y>=ymax_low && y<1.){//Eq(9)
				dy=0.001;
				dNdy=C0_highbar/(y-C1_highbar);//Eq(10)
		}
		else{
			dNdy=0.;
      			cout<<"y value is outside of the domain of y.\n";
		}
	}
	else if(nu_nubar==1){
		U=1-1/C2_low;
		W=fabs( (ymax_high-C1_highbar)/(ymin_high-C1_highbar));
		B=(pow(ymax_low-C1_low, 1/C2_low)/(ymax_low-C1_highbar));
		T=B*((pow(ymax_low-C1_low, U)-pow(ymin_low-C1_low, U) )/U)+log(W);
		C0_highbar=1/T;	
		C0_lowbar=C0_highbar*(pow(ymax_low-C1_low, 1/C2_low))/(ymax_low-C1_highbar);
		if(y<ymax_low){
			dy=0.00002;
			dNdy=C0_lowbar/pow(y-C1_low, 1/C2_low);
		}
		else if(y>=ymax_low && y<1.){
			dy=0.001;
			dNdy=C0_highbar/(y-C1_highbar);
		}
		else{
			dNdy=0;
      cout<<"y value is outside of the domain of y.\n";
		}
	}
	else{
	    cout<<"Nu_nubar is not defined!\n";
	}		
	weighty=dNdy;
	return weighty;
}//Getyweight

string Primaries::GetCurrent() {
  // choose CC or NC
  //  get from ratios in Ghandi etal paper
  // updated for the CTEQ6-DIS parton distribution functions
  string current;
  double rnd=gRandom->Rndm();
  if (rnd<=0.6865254) // 10^18 eV - 10^21 eV (use this one for ANITA)
//if (rnd<=0.6893498) // 10^17 eV - 10^20 eV (use this one for SalSA)
    current="cc";
  else
    current="nc"; 

  //hack!
  //current = "cc";
  return current;
} //GetCurrent

string Primaries::GetNuFlavor() {
    // pick a neutrino type, flavor ratio 1:1:1
  
  string nuflavor;

  double rnd=gRandom->Rndm();

  if (rnd<=(1./3.)) {  
    nuflavor="nue";
  } //if
  else if(rnd<=(2./3.)) { 
    nuflavor="numu";
  } //else if
  else if(rnd<=(1.)) { 
    nuflavor="nutau";
  } //else if
  else
    cout << "unable to pick nu flavor\n";
  
  //hack!
  // nuflavor="nutau";
  
  return nuflavor;
} //GetNuFlavor

Interaction::Interaction(string inttype,Primaries *primary1,Settings *settings1,int whichray,Counting *count1) : banana_flavor("numu"), banana_current("nc"),  nu_banana(Position(theta_nu_banana,phi_nu_banana)) {

  noway=0;
  wheredoesitleave_err=0;
  neverseesice=0;

    wheredoesitenterice_err=0;
    toohigh=0;
    toolow=0;

  iceinteraction=0;
  dtryingdirection=0.;
  dnutries=0.;
   //set to zero so the weight given to each event is correct
  /*weight_nu=0.;//Weight for neutrino that survives to posnu
  weight_nu_prob=0.;//Weight for neutrino that survives to posnu and interacts in the ice
  weight_tau=0.;//Weight for a tau neutrino to interact, create a tau, and the tau to survive to the ice. NOT USED
  weight_tau_prob=0.;//Weight for tau neutrino to interact, create a tau, tau survives and decays in the ice.
  */
    if (inttype=="banana") {
      nu_banana = (surface_over_banana_nu+altitude_nu_banana) * nu_banana;

      //Set neutrino direction
      nnu_banana = Vector(nu_banana_theta_angle + PI,nu_banana_phi_angle);
      nnu_banana = nnu_banana.ChangeCoord(nu_banana);
           
      current = banana_current;

      nuflavor = banana_flavor;


    }
    else {
    setNuFlavor(primary1,settings1,whichray,count1);
    setCurrent(primary1);
    //    setnu_nubar(primary1);//same function for inttype "banna" or otherwise.
    }
}



     
void Interaction::PickAnyDirection() {
  double rndlist[2];
  gRandom->RndmArray(2,rndlist);
  
  costheta_nutraject=2*rndlist[0]-1;

 
  // pick a neutrino azimuthal angle
  phi_nutraject=2*PI*rndlist[1];
  
  // check that these give the right result
  double thetanu=acos(costheta_nutraject);
  
  double sinthetanu=sin(thetanu);
  
  // find direction vector of neutrino
  // **** are cosine and sine flipped?
  nnu.SetX(sinthetanu*cos(phi_nutraject));
  nnu.SetY(sinthetanu*sin(phi_nutraject));
  nnu.SetZ(costheta_nutraject);


}


void  Interaction::setNuFlavor(Primaries *primary1,Settings *settings1,int whichray,Counting *counting1) {
     // pick the neutrino flavor,  type of tau decay when relevant,
      //  lpm energy.
      nuflavor=primary1->GetNuFlavor();
      
      if (settings1->MINRAY==whichray) {// only increment neutrino flavor on first ray so you don't count twice
      
	if (nuflavor=="nue")
	  counting1->nnu_e++;      
	if (nuflavor=="numu")
	  counting1->nnu_mu++;      
	if (nuflavor=="nutau")
	  counting1->nnu_tau++;      
      }
      
      if (settings1->FORSECKEL==1) // For making array of signal vs. freq, viewangle, just use muon neutrinos
	nuflavor="nue";
      if (nuflavor=="nue")  //For outputting to file
	nuflavorint=1;
      else if (nuflavor=="numu")
	nuflavorint=2;
      else if (nuflavor=="nutau")
	nuflavorint=3;
      else 
	cout<<"nuflavor is "<<nuflavor<<"\n";

}


void  Interaction::setCurrent(Primaries *primary1) {
     // pick whether it is neutral current
      // or charged current
      current=primary1->GetCurrent();


      if (current=="cc")   //For outputting to file
	currentint=1;
      else if(current=="nc")
	currentint=2;   
  }//setCurrent

///////////////////// Y ////////////////////
Y::Y() { // Constructor

  ffrac=new TF1("ffrac","[0]*sin([1]*(x-[2]))",7.,12.); // This is the fraction of the distribution in the low y region given by Equation 18. 
  
  ffrac->FixParameter(0,0.128); // These parameters are the same for all interaction types
  ffrac->FixParameter(1,-0.197);
  ffrac->FixParameter(2,21.8);

  string sbase="C1_high";
  char which[50];
  for (int i=0;i<2;i++) {
    for (int j=0;j<2;j++) {
      sprintf(which,"%d%d",i,j);
      string sname=sbase+which;
      fC1_high[i][j]=new TF1(sname.c_str(),"[0]+[1]*(-exp(-(x-[2])/[3]))",7.,12.); // parameterization of parameter C1 in the high y region according to Equation 16
    }
  }

// parameter A_0 in Table V for the high y region
  fC1_high[0][0]->FixParameter(0,-0.0026);//nubar, CC
  fC1_high[1][0]->FixParameter(0,-0.008); //nu,    CC
  fC1_high[0][1]->FixParameter(0,-0.005); //nubar, NC
  fC1_high[1][1]->FixParameter(0,-0.005); //nu,    NC

  // parameter A_1 in Table V for the high y region
  fC1_high[0][0]->FixParameter(1,0.085); // nubar, CC
  fC1_high[1][0]->FixParameter(1,0.26); // nu, CC
  fC1_high[0][1]->FixParameter(1,0.23); // nubar, NC
  fC1_high[1][1]->FixParameter(1,0.23); // nu, NC


  // parameter A_2 in Table V for the high y region
  fC1_high[0][0]->FixParameter(2,4.1); // nubar, CC
  fC1_high[1][0]->FixParameter(2,3.0); // nu, CC   
  fC1_high[0][1]->FixParameter(2,3.0); // nubar, NC
  fC1_high[1][1]->FixParameter(2,3.0); // nu, NC

  // parameter A_3 in Table V for the high y region.  This parameter is the same for all four interaction types
  for (int i=0;i<2;i++) { // nu, nubar
    for (int j=0;j<2;j++) { // CC, NC
      fC1_high[i][j]->FixParameter(3,1.7);
    }
  }

  fC1_low=new TF1("C1_low","[0]+[1]*(-exp(-(x-[2])/[3]))",7.,12.); // parameterization of parameter C1 in the low y region according to Equation 16.
  // This parameterization is the same for all interaction types.
  
  fC1_low->FixParameter(0,0.);
  fC1_low->FixParameter(1,0.0941);
  fC1_low->FixParameter(2,4.72);
  fC1_low->FixParameter(3,0.456);

  fC2=new TF1("C2","[0]+[1]*x",7.,12.); // parameterization of parameter C2 in the low y region according to Equation 17.
  // This parameterization is the same for all interaction types.
  fC2->FixParameter(0,2.55);
  fC2->FixParameter(1,-9.49E-2);

  // For picking inelasticity in low y region according to Equation 14.
   fy0_low=new TF3("fy0_low","x+(z*([1]-x)^(-1./y+1)+(1-z)*([0]-x)^(-1./y+1))^(y/(y-1))"); // x=C_1, y=C_2, z=R
   fy0_low->SetParameter(0,0.00002);  // y_min
   fy0_low->SetParameter(1,0.001); // y_max
   
   // For picking inelasticity in high y region according to Equation 15.
   fy0_high=new TF2("fy0_high","([1]-x)^y/([0]-x)^(y-1.)+x"); // x=C_1, y=R
   fy0_high->SetParameter(0,0.001); // y_min
   fy0_high->SetParameter(1,1.); // y_max

}//Y Constructor

double Y::pickY(int NU,int CURRENT,double e) {
	//e is in GeV.
  // Select a value of y that follows the appropriate distribution according to the prescription outlined in 
  // A. Connolly, R. Thorne and D. Waters, arXiv:1102.0691 [hep-ph]. 

  // pick a y region
 
  //double R1=Rand3Y.Rndm(); // choose our first random number
  double R1=gRandom->Rndm();

  int iyregion=0; // 0 for high y region 
  if (R1<ffrac->Eval(e)) // Is it going to be in low y region?
    iyregion=1; // 1 for low y region
 
  double C1_this;
  if (iyregion==0) // high y region
    C1_this=fC1_high[NU][CURRENT]->Eval(e); // C1 for this event
  else // low y region
    C1_this=fC1_low->Eval(e); // C1 for this event

  double C2_this=fC2->Eval(e); // C2 for this event
  
  // pick another random number
	double R2=gRandom->Rndm();//  double R2=Rand3Y.Rndm();
  double y0=0.;
  if (iyregion==0)  // high y region
    y0=fy0_high->Eval(C1_this,R2); // pick y0 according to Equation 15
  else if (iyregion==1)  // low y region
    y0=fy0_low->Eval(C1_this,C2_this,R2); // pick y0 according to Equation 14

  return y0;
}//pickY
