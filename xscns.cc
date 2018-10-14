#include "xscns.h"
#include <math.h>

// All energies in MeV

double A2forcm2=8.43103e-45; // ((G^2)/(2  Pi)) *hbarcinmeters^-4*(100)^2

/////////////////////////
// Differential cross section in cm^2 MeV^-1  for combination with GV, GA
// Input energy in MeV

double diffxscnvec(double knu, double mass, double erec) {

  double diffxscn = A2forcm2*mass*(2-mass*erec/(knu*knu)-2*erec/knu+erec*erec/(knu*knu));
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

///////////////////////

// SM parameters

void sm_vector_couplings(int pdgyear, double* gv) {

  // Default is 2015 PDG, from Erler and Su paper.  This is for mu flavor

  // This is subtracting the charge radius correction for the 
  // proton coupling; these are values in the table for mu flavor,
  // and include -phinulWmu*2 (factor of 2 for V=L+R)= 
  // correction only for proton

  //  double chgradcorr = 2.*0.00571903;

    // For no charge radius correction
  //double gVp = 0.0298-chgradcorr;

  double gVp = 0.01836;
  double gVn= -0.5117;

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
  

//   if (pdgyear < 2015) {
//     //   double Szhat2 = 0.23875;
//     double Szhat2 = 0.23120;
//     double rhoNeuNucNC = 1.0086;
//     double khatNeuNuc = 0.9978;
//     double lambdaUL = -0.0031;
//     double lambdaDL = -0.0025;
//     double lambdaDR = 7.5e-5;
//     double lambdaUR = lambdaDR/2.0;


//    gVp = rhoNeuNucNC*(0.5 - 2.0 * khatNeuNuc * Szhat2) + 2.0*lambdaUL + 2.0*lambdaUR + lambdaDL + lambdaDR;

//    //      gVn = -0.5*rhoNeuNucNC + 2.*lambdaUL + lambdaUR + 2.0*lambdaDL + 2.0*lambdaDR;

//       gVn = -0.5*rhoNeuNucNC +  lambdaUL + lambdaUR + 2.0*lambdaDL + 2.0*lambdaDR;

//   }

  gv[0] = gVp;
  gv[1] = gVn;
  

}


double chgradcorr(int flavor) {

  // Correction to PDG from Erler
  double gvpcorr = 0.;
  flavor = fabs(flavor);
  if (flavor == 1) {
    gvpcorr = 0.0196964;
  } else if (flavor == 2) {
    gvpcorr = 0.0114381;
  }
  else if (flavor == 3) {
    gvpcorr = 0.00706633;
  }

  return gvpcorr;
 
}

void sm_axial_couplings(int pdgyear, int flav, double* ga) {

  // Default is 2015 PDG, from Erler paper

  // flav = 1: electron
  // flav = 2: muon
  // Negative sign means antineutrino
  // Actually this doesn't care the flavor if no NSIs, just cares about sign
  
  // Ignoring nutau in initial state

  //  double sstw = 0.231;

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


double GV_SM(int pdgyear, int Z, int N) {

  // Default is 2015 PDG, from Erler paper

  double gv[2];
  sm_vector_couplings(pdgyear, gv);
  double gVp = gv[0];
  double gVn = gv[1];

  return (Z*gVp + N*gVn);

}


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

