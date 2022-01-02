#include "xscns.h"
#include <math.h>
#include <iostream>

// All energies in MeV


const double A2forcm2=8.43103e-45; // ((G^2)/(2  Pi)) *hbarcinmeters^-4*(100)^2
const double hbarcincm=197.327e-13; // MeV-cm
const double Anuelforcm2 = 1.7233e-44; // ((2 G^2 me )/(Pi)) *hbarcinmeters^-4*(100)^2

/////////////////////////
// Differential cross section in cm^2 MeV^-1  for combination with GV, GA
// Input energy in MeV

double diffxscnvec(double knu, double mass, double erec) {

  double diffxscn = A2forcm2*mass*(2-mass*erec/(knu*knu)-2*erec/knu+erec*erec/(knu*knu));
        //      std::cout << erec<<" "<<q<<" "<<F<<" "<<std::endl;
  return diffxscn;

}

// Need cm^2

double diffangdist(double knu,  double cth) {

  // Only higher order terms kinematically
  double diffxscn = A2forcm2*knu*knu*(1+cth);
        //      std::cout << erec<<" "<<q<<" "<<F<<" "<<std::endl;
  return diffxscn;

}



double diffxscnaxial(double knu, double mass, double erec) {
  
        //      std::cout << erec<<" "<<q<<" "<<F<<" "<<std::endl;
  return A2forcm2*mass*(2+mass*erec/(knu*knu)-2*erec/knu+erec*erec/(knu*knu));

}

double diffxscninterf(double knu, double mass, double erec) {
  
        //      std::cout << erec<<" "<<q<<" "<<F<<" "<<std::endl;
  return A2forcm2*mass*(4*erec/knu-2*erec*erec/(knu*knu));

}


double diffxscnmag(double knu, double erec) {

  double me=0.51099895;
  double alpha = 0.0072973525664;
  

  return M_PI*pow(alpha,2)*pow(hbarcincm,2)/pow(me,2)*(1/erec-1/knu + erec/(4*knu*knu));

}


///////////////

// Neutrino-electron scattering cross section

double diffnuelectronxscn(int flav,  double knu, double erec) {

  // flav = 1: electron
  // flav = 2: muon
  // flav = 3: tau
  // flav is negative for antineutrinos



  double gV, gA;
  const double ssthW = 0.231; 
  
  switch (flav) {
  case 1:
    gV = 2*ssthW+0.5;
    gA = 0.5;
    break;
  case 2: 
    gV = 2*ssthW-0.5;
    gA = -0.5;
    break;
  case 3: 
    gV = 2*ssthW-0.5;
    gA = -0.5;
    break;
  case -1: 
    gV = 2*ssthW+0.5;
    gA = -0.5;
    break;
  case -2: 
    gV = 2*ssthW-0.5;
    gA = 0.5;
    break;
  case -3: 
    gV = 2*ssthW-0.5;
    gA = 0.5;
    break;
    std::cout<< "Wrong flavor "<<std::endl;
  }

  double me=0.51099895;

  // std::cout << "no pol "<<Anuelforcm2/4.*(pow(gA+gV,2))<<" pow(gA+gV,2) "<<pow(gA+gV,2)<<" "<<std::endl;
  return Anuelforcm2/4.*(pow(gA+gV,2)+ pow(gV-gA,2)*(1-erec/knu)*(1-erec/knu)+(gV*gV+gA*gA)*me*erec/(knu*knu));

  
}
			  


///////////////

// Neutrino-electron scattering cross section for nonzero mass, from Barranco et al.

double diffnuelectronxscn2(int flav,  double knu, double erec, double spar, double mnu, int maj) {

  // flav = 1: electron
  // flav = 2: muon
  // flav = 3: tau
  // flav is negative for antineutrinos

  // spar is parallel component of spin, nominally -1 for zero mass

  // maj = 1: Majorana
  // maj = 0: Dirac

  double gV, gA;
  const double ssthW = 0.231; 
  
  switch (flav) {
  case 1:
    gV = 2*ssthW+0.5;
    gA = 0.5;
    break;
  case 2: 
    gV = 2*ssthW-0.5;
    gA = -0.5;
    break;
  case 3: 
    gV = 2*ssthW-0.5;
    gA = -0.5;
    break;
  case -1: 
    gV = 2*ssthW+0.5;
    gA = -0.5;
    break;
  case -2: 
    gV = 2*ssthW-0.5;
    gA = 0.5;
    break;
  case -3: 
    gV = 2*ssthW-0.5;
    gA = 0.5;
    break;
    std::cout<< "Wrong flavor "<<std::endl;
  }

  double me=0.51099895;
  double xscn=0.;
  double pnu = pow(knu*knu-mnu*mnu, 0.5);
  if (maj == 0) {
    //   std::cout <<" pnu "<<pnu<< "pol "<<Anuelforcm2/8.*(pow(gA+gV,2)*knu*(knu-spar*pnu)/(pnu*pnu))<<" factor "<<
    //knu*(knu-spar*pnu)/(pnu*pnu)<<" pow(gA+gV,2) "<<pow(gA+gV,2)<<std::endl;
    xscn = Anuelforcm2/8.*(pow(gA+gV,2)*knu*(knu-spar*pnu)/(pnu*pnu)
			   + (pow(gV-gA,2)*pow(knu-erec,2)/(pnu*pnu)+(gA*gA-gV*gV)*me*erec/(pnu*pnu))
			   *(1-spar*knu/pnu)
			   + (spar/pnu*pow(gA-gV,2)*(knu-erec)*(1+erec/me)+(gA*gA-gV*gV)*(1-spar*erec/pnu))
			   * pow(mnu/pnu,2));
  } else {
    xscn = Anuelforcm2/(4.*pnu*pnu)*(2*(2*gA*gA-gV*gV)*mnu*mnu
			     - 4*knu/pnu*gA*gV*spar*erec*(knu+mnu*mnu/me)*(1-erec/(2*knu))
			     + (gA*gA-gV*gV)*me*erec
			     + (gA*gA+gV*gV)*(2*knu*knu*(1-erec/knu)+erec*erec*(1+mnu*mnu/(me*erec))));
      
  }

  return xscn;
}
			  


///////////////////////

// SM parameters

void sm_vector_couplings(int pdgyear, double* gv) {

  // Default is 2015 PDG, from Erler and Su paper.  This is for mu flavor.  Erler and Su has charge radius correction for mu flavor, so remove it for default not correction.

  // This is subtracting the charge radius correction for the 
  // proton coupling; these are values in the table for mu flavor,
  // and include -phinulWmu*2 (factor of 2 for V=L+R)= 
  // correction only for proton

  //  double chgradcorr = 2.*0.00571903;

    // For no charge radius correction
  //double gVp = 0.0298-chgradcorr;

  // This is for Giunti couplings
  double gVp = 0.0227;
  double gVn= -0.5117;

  if (pdgyear == 0) {
    // Custom couplings
    std::cout << "Custom couplings requested"<<std:: endl;
    exit(0);
  }
  
  if (pdgyear < 2004){
    gVp = 0.0152;
    gVn = -0.5122;
  }

  if (pdgyear >= 2004 && pdgyear < 2011){
    gVp = 0.0304;
    gVn = -0.5122;
  }

  if (pdgyear >= 2011 && pdgyear < 2012){
    gVp = 0.0306;
    gVn = -0.5120;
  }

  if (pdgyear >= 2012 && pdgyear < 2014){
    gVp = 0.0307;
    gVn = -0.5120;
  }

  if (pdgyear >= 2014 && pdgyear < 2020){

    // Erler and Su with Giunti charge correction removed for numu
    gVp = 0.01836;
    gVn= -0.5117;
  }

//    if (pdgyear < 2015) {
//      double Szhat2 = 0.23875;
//      //     double Szhat2 = 0.23120;
//      double rhoNeuNucNC = 1.0086;
//      double khatNeuNuc = 0.9978;
//      double lambdaUL = -0.0031;
//      double lambdaDL = -0.0025;
//      double lambdaDR = 7.5e-5;
//      double lambdaUR = lambdaDR/2.0;


//     gVp = rhoNeuNucNC*(0.5 - 2.0 * khatNeuNuc * Szhat2) + 2.0*lambdaUL + 2.0*lambdaUR + lambdaDL + lambdaDR;

// //    //      gVn = -0.5*rhoNeuNucNC + 2.*lambdaUL + lambdaUR + 2.0*lambdaDL + 2.0*lambdaDR;

//        gVn = -0.5*rhoNeuNucNC +  lambdaUL + lambdaUR + 2.0*lambdaDL + 2.0*lambdaDR;

//    }

  gv[0] = gVp;
  gv[1] = gVn;
  

}


double chgradcorr(int flavor, int type) {

  double gvpcorr = 0.;
  int sign = fabs(flavor)/flavor;
  flavor = fabs(flavor);

  if (type == 1) {
  // Correction to PDG from Erler
    if (flavor == 1) {
      gvpcorr = 0.0196964;
    } else if (flavor == 2) {
      gvpcorr = 0.0114381;
    }
    else if (flavor == 3) {
      gvpcorr = 0.00706633;
    }

  } else if (type==2) {

   // Updated Giunti but no sign change for antineutrinos, and slightly tweaked numbers according to the reference
    
    if (flavor == 1) {
      gvpcorr = 0.0212196;
    } else if (flavor == 2) {
      gvpcorr = 0.0122716;
    }
    else if (flavor == 3) {
      gvpcorr = 0.00766972;
    }

  }
  else if (type == 3) {

    // Giunti corrections

    if (flavor == 1) {
      gvpcorr = 0.02191;
    } else if (flavor == 2) {
      gvpcorr = 0.01267;
    }
    else if (flavor == 3) {
      gvpcorr = 0.00792;
    }

    // Opposite sign for antineutrinos-- note this is wrong (see Errata 2020 for Giunti) but included as legacy
    gvpcorr *= sign;

  } 

  return gvpcorr;
 
}

// For use in Tomalak charge radius correction
double Pifunc(double Q, double mf, double mu) {

  // Q: momentum transfer
  // Mf: lefton mass 
  // Mu: renormalization scale in GeV

  double Q2 = Q*Q;
  double mf2 = mf*mf;
  double A = sqrt(1+ 4*mf2/Q2);

  double Pival=0;
  if (Q2!=0) {

    Pival = (1/3.)*(2*log(mu)-log(mf2))
    +5./9. - 4*mf2/(3*Q2) + (1/3.)*(1-2*mf2/Q2)*A*(log(A-1)-log(A+1));

  } else {
    // Limit for small Q
    Pival = (1./3.)*(2*log(mu)-log(mf2))-Q2/(15.*mf2);

  }


  return Pival;

}

double chgradcorr_tomalak(double Q, int flav){

  // return Q-dependent flavor correction to gVP
  // Based on Tomalak et al., arXiv:2011.05960v2


    //A = alpha/(2*M_PI)/(sqrt(2)*GF)

  // Q-independent part
  // qcdcorr= alpha*deltaQCD/(2.*M_PI)
  
  const double qcdcorr = -0.00498612;
  const double A = 72.3779;
  double deltaf = 0;
  const double cLee = 2.41198e-05;
  const double cLemu = -8.8704e-06;
  const double cLmumu = cLee;
  const double cLtautau = cLee;
  const double cLetau = cLemu;
  const double cLmue = cLemu;
  const double cLmutau = cLemu;
  const double cLtaue = cLemu;
  const double cLtaumu = cLemu;
  const double cRl = 7.62469e-06;

  const double mu = 2000.;

  const double me = 0.510998;
  const double mmu = 105.6583755;
  const double mtau = 1776.86;  
  
  double chgradcorr = qcdcorr;
  if (abs(flav) ==1 ) {
    deltaf=A*((cLee+cRl)*Pifunc(Q,me,mu)+ (cLemu+cRl)*Pifunc(Q,mmu,mu)+(cLetau+cRl)*Pifunc(0,mtau,mu));
  } 

  if (abs(flav) ==2 ) {
    deltaf=A*((cLmue+cRl)*Pifunc(Q,me,mu)+ (cLmumu+cRl)*Pifunc(Q,mmu,mu)+(cLmutau+cRl)*Pifunc(0,mtau,mu));
  }

  if (abs(flav) ==3 ) {
    deltaf=A* ((cLetau+cRl)*Pifunc(Q,me,mu)+ (cLtaumu+cRl)*Pifunc(0,mmu,mu)+(cLtautau+cRl)*Pifunc(Q,mtau,mu));
  } 

  chgradcorr += deltaf;
  return chgradcorr;

}




void sm_axial_couplings(int pdgyear, int flav, double* ga) {

  // Default is 2015 PDG, from Erler paper

  // flav = 1: electron
  // flav = 2: muon
  // Negative sign means antineutrino
  // Actually this doesn't care the flavor if no NSIs, just cares about sign
  
  // Ignoring nutau in initial state

  //  double sstw = 0.231;

  if (pdgyear == 0) {
    // Custom couplings
    std::cout << "Custom couplings requested"<<std:: endl;
    exit(0);
  }
  
  double gAp =  0.4995*flav/fabs(flav);
  double gAn = -0.5121*flav/fabs(flav);

  // Actually I don't know where the 1.27 comes from-- Phil's paper

  if (pdgyear < 2004) {
    double dS = -0.15;
    
    gAp = (1.27-dS)/2.*flav/fabs(flav);
    gAn = -(1.27-dS)/2.*flav/fabs(flav);

  }

  if (pdgyear >= 2004 && pdgyear < 2011){
    gAp = 0.4955*flav/fabs(flav);
    gAn = -0.5125*flav/fabs(flav);
  }

  if (pdgyear >= 2011 && pdgyear < 2012){
    gAp = 0.4942*flav/fabs(flav);
    gAn = -0.5123*flav/fabs(flav);
  }

  if (pdgyear >= 2012 && pdgyear < 2014){
    gAp = 0.4953*flav/fabs(flav);
    gAn = -0.5124*flav/fabs(flav);
  }


  ga[0] = gAp;
  ga[1] = gAn;
  
}



// Not used
double GV_SM(int pdgyear, int Z, int N) {

  // Default is 2015 PDG, from Erler paper

  double gv[2];
  sm_vector_couplings(pdgyear, gv);
  double gVp = gv[0];
  double gVn = gv[1];

  return (Z*gVp + N*gVn);

}


// Not used
double GA_SM(int pdgyear, int flav, int Z, int N, int Zdiff, int Ndiff){

  double ga[2];
  sm_axial_couplings(pdgyear, flav, ga);
  double gAp = ga[0];
  double gAn = ga[1];
  
  return (gAp*Zdiff+gAn*Ndiff);

}

// These are from poly fit to Seghal paper charge ratios.  Not quite right at very low Q, take const value at 0.1.. It's wrt electron

double mufactor(double Q) {

  double logQ2 = log10(Q*Q);

  double mufact;
  if (logQ2<-1.) {
    logQ2=-1.;
  }

  if (logQ2>5.) {
    mufact = 1.;
  } else  {
    mufact = 1.0344 -0.00355405*logQ2 - 0.00146375*logQ2*logQ2 + 8.59952e-6*pow(logQ2,3) + 4.64804e-05*pow(logQ2,4)-2.95199e-06*pow(logQ2,5);
  }
    
  return mufact;
}

double taufactor(double Q) {

  double logQ2 = log10(Q*Q);

  double taufact;
  if (logQ2<-1.) {
    logQ2=-1.;
  }

  if (logQ2>6.) {
    taufact = 1.;
  } else  {
    taufact = 1.05292 -0.0034373*logQ2 - 0.00133684*logQ2*logQ2 + 3.41384e-05*pow(logQ2,3) + 4.28903e-05*pow(logQ2,4)-4.74854e-06*pow(logQ2,5);
  }
    
  return taufact;
}


//////////////////////

// NSI couplings... first two routines not used


// This is the NSI factor for the T<<Enu approximation
double GV_nsi_nonuniv(int flav, int Z, int N, double eeeuV, double eeedV, double eetauuV, double eetaudV, double eemuuV, double eemudV, double emumuuV, double emumudV, double emutauuV, double emutaudV){

  // flav = 1: electron
  // flav = 2: muon

  // Ignoring nutau in initial state

  //  cout << "Qw "<<Qw<<endl;
  
  double nonuniv=0;


  if (flav == 1 ) {
    nonuniv = Z*(2*eeeuV+ eeedV)+N*(eeeuV+2*eeedV);
  }
  else if (flav== 2) {
    nonuniv = Z*(2*emumuuV+ emumudV)+N*(emumuuV+2*emumudV);
  }
    

  return nonuniv;

}


// Already squared.
double GV_nsi_fc2(int flav, int Z, int N, double eeeuV, double eeedV, double eetauuV, double eetaudV, double eemuuV, double eemudV, double emumuuV, double emumudV, double emutauuV, double emutaudV){

  // flav = 1: electron
  // flav = 2: muon

  // Ignoring nutau in initial state

  
  double fc=0;

  //double nsi_factor = 0;
  if (flav == 1 ) {
    fc = Z*(2*eemuuV+eemudV)+N*(eemuuV+2*eemudV)
      +Z*(2*eetauuV+eetaudV)+N*(eetauuV+2*eetaudV);
  }
  else if (flav== 2) {
    fc = pow(Z*(2*eemuuV+eemudV)+N*(eemuuV+2*eemudV),2)
      +pow(Z*(2*emutauuV+emutaudV)+N*(emutauuV+2*emutaudV),2);
  }
    

  return fc;

}

void nsi_vector_couplings(double* nsi_gv, double epsilonu, double epsilond) {

  // Output for proton and neutrons separately
 
  double nsip;
  double nsin;

  nsip = 2*epsilonu + epsilond;
  nsin = epsilonu + 2*epsilond;
  

  nsi_gv[0] = nsip;
  nsi_gv[1] = nsin;
  

}

// Ignoring axial NSIs
