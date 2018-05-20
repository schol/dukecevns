#include "xscns.h"
#include <math.h>

// All energies in MeV

double A2forcm2=8.43103e-45; // ((G^2)/(2  Pi)) *hbarcinmeters^-4*(100)^2

/////////////////////////
// Differential cross section in cm^2 MeV^-1  in combination with GV, GA
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

double GV_SM(int Z, int N) {

  // 2015 PDG, from Erler paper
  //  double sstw = 0.231;
  double gVp = 0.0298;

  double gVn= -0.5117;

  return (Z*gVp + N*gVn);

}


double GA_SM(int flav, int Z, int N, int Zdiff, int Ndiff){

  // flav = 1: electron
  // flav = 2: muon
  // Negative sign means antineutrino
  // Actually this doesn't care the flavor if no NSIs, just cares about sign
  
  // Ignoring nutau in initial state

  //  double sstw = 0.231;

  // Actually I don't know where the 1.27 comes from-- Phil's paper
  double gAp =  0.4995*flav/fabs(flav);
  double gAn = -0.5121*flav/fabs(flav);


  return (gAp*Zdiff+gAn*Ndiff);

}

