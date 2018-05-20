#include <iostream>
#include <fstream>

#include "FormFactor.h"
#include "NuFlux.h"
#include "DetectorResponse.h"
#include "TFile.h"
#include "TGraph.h"
#include "TObjArray.h"
#include "TString.h"
#include "TMath.h"
#include "TH2D.h"

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>
#include <math.h>
#include <map>

#include "xscns.h"

int main(int argc, char * argv[] )
{

#include "isomaps.h"
  // Info for relevant mixtures
#include "mixtures.h"

  
  if (argc<3) {
    std::cout << "Usage:  ./diff_rates [target material] [form factor] [Rnfac]"<<std::endl;
    exit(0);
  }


  // Set up the form factor

  const char * ffname = argv[2];

  double rnfac=1.;
  if (argc >=4) {
     rnfac = (double)atof(argv[3]);
  }


  // Array of pointers to form factors
  FormFactor** ff;
  ff = new FormFactor*[max_components];

    // Set up the flux  (here, fluence)

 //  NumericalFlux* livermore = new NumericalFlux();
//   livermore->SetFluxFilename("livermore.dat");
//   livermore->ReadFluxFile();
//   double kmax = livermore->maxEnu();
//   std::cout << "Max neutrino energy: "<<kmax<<std::endl;

  PinchedThermal* pinched = new PinchedThermal();

  double alpha[6]={3.,2.5,2.5,3.,2.5,2.5};
  double avgen[6]={10.,15.,15.,14.,15.,15.};
  double lumi[6]={1.6e52,1.6e52,1.6e52,1.6e52,1.6e52,1.6e52};
  pinched->SetAlpha(alpha);
  pinched->SetLuminosity(lumi);
  pinched->SetAvgEn(avgen);
  double kmax = pinched->maxEnu();

    // Set up the material

  std::string material = argv[1];

  std::ofstream outfile;
  std::string outfilename;
  outfilename = "supernova_diff_rates-"+material+"-"+std::string(ffname)+".out";
  outfile.open(outfilename);


  std::ofstream phoutfile;
  std::string phoutfilename;
  phoutfilename = "supernova_diff_rates-"+material+"-"+std::string(ffname)+"-photons.out";
  phoutfile.open(phoutfilename);


  // Set up a detector quenching factor

   DetectorResponse* arqf = new DetectorResponse();
   arqf->SetQFPolyFilename("arpoly.txt");
   arqf->ReadQFPolyFile();
  
   // Array for quenched total rates

   const int maxiq = 10000;
   double Eee[maxiq];
   double dNdEee[maxiq];


  //  std::cout << "Material "<<material << std::endl;

  //  std::string material = "Ar";

  double M;
  double Delta;
  int Nn,Z,A;
  int Zdiff, Ndiff;

  std::string matname = material;

 // These are defined in mixtures include
  std::vector<double> fraction = molar_fraction[material];
  std::vector<std::string> isotope_component= isotopes[material];

  std::cout << "Material "<<matname<<std::endl;


  int is=0;
  std::vector<std::string>::iterator v = isotope_component.begin();
  std::string isotope;

  // First get the total mass.  Get also the maximum recoil values

  double erecmaxvals[max_components];
  double Mtot = 0;
  v = isotope_component.begin();

  double minM = 1.e10;
  while( v != isotope_component.end()) {

    isotope = *v;
        std::cout << "isotope"<< isotope << std::endl;

    Z = Zs[std::string(isotope)];
    Nn = Ns[std::string(isotope)];
    Delta = Deltas[std::string(isotope)];
    M = (Z+Nn)*amu - Z*me + Delta;
    if (M<minM) {minM=M;}
    Mtot += M*fraction[is];
    erecmaxvals[is] = 2*kmax*kmax/(M+2*kmax);

 


    // Set up the form factor for this isotope

    //       std::cout << "Mass "<<M<<std::endl;


    if (strcmp(ffname, "helm")==0) {
      Helm* helmff= new Helm();
      ff[is] = helmff;
      helmff->Setsval(0.9);
    }
    else if (strcmp(ffname, "klein")==0) {
      Klein* kleinff = new Klein();
      ff[is] = kleinff;
      kleinff->Setakval(0.7);
    } 
    else if  (strcmp(ffname, "horowitz")==0){
      Horowitz* horowitzff = new Horowitz();
      ff[is] = horowitzff;
      std::string isoname = std::string(isotope);
      std::transform(isoname.begin(), isoname.end(),isoname.begin(), ::toupper);
      std::string horowitz_filename = isoname+".FF";
      //      std::cout << horowitz_filename << std::endl;
      horowitzff->SetFFfilename(horowitz_filename.c_str());
      horowitzff->ReadFFfile();
      horowitzff->SetRnfac(rnfac);
    }

    A = Nn + Z;
    ff[is]->SetA(A);
 

    v++; is++;
  }



  // Use the mass of the lightest component
  double erecmaxall = 2*kmax*kmax/(minM+2*kmax);
  
  double erecstart = 0.;
  double erecend = erecmaxall;
  double erecstep = 0.0001;

   // Now compute the differential recoil spectra

     
    double Erec;
    double knu;

    std::cout << "erecmaxall "<<erecmaxall<<std::endl;

   double knustep = 0.0001;

   // The totals
   double toterecoil = 0.;
   double totevents = 0.;
   

   int iq=0;
   // Loop over recoil energy
   for (Erec=erecstart+erecstep;Erec<=erecend; Erec+=erecstep) {
     

	double diffrate_e_vec = 0;
	double diffrate_ebar_vec = 0;
	double diffrate_mu_vec = 0;
	double diffrate_mubar_vec = 0;
	double diffrate_tau_vec = 0;
	double diffrate_taubar_vec = 0;

	double diffrate_e_axial = 0;
	double diffrate_ebar_axial = 0;
	double diffrate_mu_axial = 0;
	double diffrate_mubar_axial = 0;
	double diffrate_tau_axial = 0;
	double diffrate_taubar_axial = 0;

	double diffrate_e_interf = 0;
	double diffrate_ebar_interf = 0;
	double diffrate_mu_interf = 0;
	double diffrate_mubar_interf = 0;
	double diffrate_tau_interf = 0;
	double diffrate_taubar_interf = 0;


	v = isotope_component.begin();
	// Now loop over components
	is=0;
	while( v != isotope_component.end()) {
	  
	  isotope = *v;
	  //	  std::cout << "isotope"<< isotope << std::endl;
	  
	  Z = Zs[std::string(isotope)];
	  Nn = Ns[std::string(isotope)];
	  Delta = Deltas[std::string(isotope)];
	  M = (Z+Nn)*amu - Z*me + Delta;
    


	  Zdiff = Zdiffs[std::string(isotope)];
	  Ndiff = Ndiffs[std::string(isotope)];
	    
	  mass_fraction[is] = M/Mtot*fraction[is];
	    
	  A = Nn + Z;
	    //	  std::cout << " Z "<<Z<<" N "<<Nn<<" A "<<A<<" M "<<M << " "<<mass_fraction[is]<<std::endl;

     // Loop over neutrino energy contributions
	  
	  // Minimum neutrino energy contributing to a given recoil energy

	  double knumin = 0.5*(Erec+sqrt(Erec*Erec+2*M*Erec));
	  Double_t hbarc = 197.327; // MeV-fm, convert for Q in MeV for ff
	  double Q = sqrt(2*M*Erec); // MeV
	  double qq = Q/hbarc;
	  //    double ff2 = helmff->FFval(qq);

	  double ff2 = ff[is]->FFval(qq);

	  // SM Couplings

	  double GV_sm = GV_SM(Z,Nn);
	  double GA_sm = GA_SM(1,Z,Nn,Zdiff,Ndiff);
	  double GA_sm_bar = GA_SM(-1,Z,Nn,Zdiff,Ndiff);

	// Normalize for one ton of material
	// Weight by mass fraction
	
	  double Nt = 1.e6/(M/amu)*6.022e23;

	// A2: G^2/(2Pi) * hbarcinmeters^-4 
	  double norm = Nt;

	  double drate_e_vec=0;
	  double drate_ebar_vec=0;
	  double drate_mu_vec=0;
	  double drate_mubar_vec=0;
	  double drate_tau_vec=0;
	  double drate_taubar_vec=0;


	  double drate_e_axial=0;
	  double drate_ebar_axial=0;
	  double drate_mu_axial=0;
	  double drate_mubar_axial=0;
	  double drate_tau_axial=0;
	  double drate_taubar_axial=0;


	  double drate_e_interf=0;
	  double drate_ebar_interf=0;
	  double drate_mu_interf=0;
	  double drate_mubar_interf=0;
	  double drate_tau_interf=0;
	  double drate_taubar_interf=0;

	  for (knu=knumin;knu<=kmax;knu+=knustep) {

	  	  //Energy in GeV
	    //	    double energy = knu*0.001;

	    drate_e_vec += diffxscnvec(knu,M,Erec)*pinched->fluxval(knu,1,knustep);
	    drate_ebar_vec += diffxscnvec(knu,M,Erec)*pinched->fluxval(knu,-1,knustep);
	    drate_mu_vec += diffxscnvec(knu,M,Erec)*pinched->fluxval(knu,2,knustep);
	    drate_mubar_vec += diffxscnvec(knu,M,Erec)*pinched->fluxval(knu,-2,knustep);
	    drate_tau_vec +=  diffxscnvec(knu,M,Erec)*pinched->fluxval(knu,3,knustep);
	    drate_taubar_vec +=  diffxscnvec(knu,M,Erec)*pinched->fluxval(knu,-3,knustep);


	    drate_e_axial += diffxscnaxial(knu,M,Erec)*pinched->fluxval(knu,1,knustep);
	    drate_ebar_axial += diffxscnaxial(knu,M,Erec)*pinched->fluxval(knu,-1,knustep);
	    drate_mu_axial += diffxscnaxial(knu,M,Erec)*pinched->fluxval(knu,2,knustep);
	    drate_mubar_axial += diffxscnaxial(knu,M,Erec)*pinched->fluxval(knu,-2,knustep);
	    drate_tau_axial +=  diffxscnaxial(knu,M,Erec)*pinched->fluxval(knu,3,knustep);
	    drate_taubar_axial +=  diffxscnaxial(knu,M,Erec)*pinched->fluxval(knu,-3,knustep);


	    drate_e_interf += diffxscninterf(knu,M,Erec)*pinched->fluxval(knu,1,knustep);
	    drate_ebar_interf += diffxscninterf(knu,M,Erec)*pinched->fluxval(knu,-1,knustep);
	    drate_mu_interf += diffxscninterf(knu,M,Erec)*pinched->fluxval(knu,2,knustep);
	    drate_mubar_interf += diffxscninterf(knu,M,Erec)*pinched->fluxval(knu,-2,knustep);
	    drate_tau_interf +=  diffxscninterf(knu,M,Erec)*pinched->fluxval(knu,3,knustep);
	    drate_taubar_interf += diffxscninterf(knu,M,Erec)*pinched->fluxval(knu,-3,knustep);

	      
	    //	    std::cout << Erec << " "<<knu<<" "<<std::endl;

	  } // End of loop over neutrino energy contributions
	    

	  // Now multiply by target-dependent factors and add up this recoil energy bin

	  diffrate_e_vec += norm*pow(GV_sm,2)*ff2*mass_fraction[is]*drate_e_vec;
	  diffrate_ebar_vec += norm*pow(GV_sm,2)*ff2*mass_fraction[is]*drate_ebar_vec;
	  diffrate_mu_vec += norm*pow(GV_sm,2)*ff2*mass_fraction[is]*drate_mu_vec;
	  diffrate_mubar_vec += norm*pow(GV_sm,2)*ff2*mass_fraction[is]*drate_mubar_vec;
	  diffrate_tau_vec +=  norm*pow(GV_sm,2)*ff2*mass_fraction[is]*drate_tau_vec;
	  diffrate_taubar_vec += norm*pow(GV_sm,2)*ff2*mass_fraction[is]*drate_taubar_vec;


	  diffrate_e_axial += norm*pow(GA_sm,2)*ff2*mass_fraction[is]*drate_e_axial;
	  diffrate_ebar_axial += norm*pow(GA_sm_bar,2)*ff2*mass_fraction[is]*drate_ebar_axial;
	  diffrate_mu_axial += norm*pow(GA_sm,2)*ff2*mass_fraction[is]*drate_mu_axial;
	  diffrate_mubar_axial += norm*pow(GA_sm_bar,2)*ff2*mass_fraction[is]*drate_mubar_axial;
	  diffrate_tau_axial +=  norm*pow(GA_sm,2)*ff2*mass_fraction[is]*drate_tau_axial;
	  diffrate_taubar_axial +=  norm*pow(GA_sm_bar,2)*ff2*mass_fraction[is]*drate_taubar_axial;


	  diffrate_e_interf += norm*GV_sm*GA_sm*ff2*mass_fraction[is]*drate_e_interf;
	  diffrate_ebar_interf += norm*GV_sm*GA_sm_bar*ff2*mass_fraction[is]*drate_ebar_interf;
	  diffrate_mu_interf += norm*GV_sm*GA_sm*ff2*mass_fraction[is]*drate_mu_interf;
	  diffrate_mubar_interf += norm*GV_sm*GA_sm_bar*ff2*mass_fraction[is]*drate_mubar_interf;
	  diffrate_tau_interf +=  norm*GV_sm*GA_sm*ff2*mass_fraction[is]*drate_tau_interf;
	  diffrate_taubar_interf +=  norm*GV_sm*GA_sm_bar*ff2*mass_fraction[is]*drate_taubar_interf;

	    //	    std::cout << is<<" "<<"Erec "<<Erec<<" mass frac "<<mass_fraction[is]<<" "<<ff2<<" "<<diffrate_e<<std::endl;

	  v++;is++;

	  
	} // End of loop over material components

	// This is events per MeV

	double events=0;
	events = diffrate_e_vec + diffrate_ebar_vec + diffrate_mu_vec+ diffrate_mubar_vec+ diffrate_tau_vec + diffrate_taubar_vec;
	events += diffrate_e_axial + diffrate_ebar_axial + diffrate_mu_axial+ diffrate_mubar_axial+ diffrate_tau_axial + diffrate_taubar_axial;
        events += diffrate_e_interf + diffrate_ebar_interf + diffrate_mu_interf+ diffrate_mubar_interf+ diffrate_tau_interf + diffrate_taubar_interf;

	std::cout << Erec<<" "<<events<<" "<<diffrate_e_vec<<" "<<diffrate_ebar_vec<<" "<<diffrate_mu_vec<<" "<<diffrate_mubar_vec<<" "<<diffrate_tau_vec<<" "<<diffrate_taubar_vec<<" "<<diffrate_e_axial<<" "<<diffrate_ebar_axial<<" "<<diffrate_mu_axial<<" "<<diffrate_mubar_axial<<" "<<diffrate_tau_axial<<" "<<diffrate_taubar_axial<<" "<<diffrate_e_interf<<" "<<diffrate_ebar_interf<<" "<<diffrate_mu_interf<<" "<<diffrate_mubar_interf<<" "<<diffrate_tau_interf<<" "<<diffrate_taubar_interf <<std::endl;
	outfile << Erec<<" "<<events<<" "<<diffrate_e_vec<<" "<<diffrate_ebar_vec<<" "<<diffrate_mu_vec<<" "<<diffrate_mubar_vec<<" "<<diffrate_tau_vec<<" "<<diffrate_taubar_vec<<" "<<diffrate_e_axial<<" "<<diffrate_ebar_axial<<" "<<diffrate_mu_axial<<" "<<diffrate_mubar_axial<<" "<<diffrate_tau_axial<<" "<<diffrate_taubar_axial<<" "<<diffrate_e_interf<<" "<<diffrate_ebar_interf<<" "<<diffrate_mu_interf<<" "<<diffrate_mubar_interf<<" "<<diffrate_tau_interf<<" "<<diffrate_taubar_interf <<std::endl;
   
	// Now get the total for the distribution taking into account the quenching factor
	// Note this will be per MeVee, but uneven bins

	Eee[iq] = arqf->qfpoly(Erec)*Erec;
	// events is dNderec
	double qfderiv = abs(arqf->qfpolyderiv(Erec));

	if (qfderiv>0) {
	  dNdEee[iq] = events/qfderiv;
	} else {
	  dNdEee[iq] = 0.;
	}

	phoutfile  << Erec<< " "<<events<<" "<<Eee[iq]<<" "<<dNdEee[iq]<<std::endl;
	iq++;

	totevents+=events*erecstep;

	toterecoil += events*Erec*erecstep;

  } // End of loop over Erec

   double time_interval = 10.; // Assume flux is over 10 seconds
   std::cout << "Total events:  "<< totevents*time_interval<< std::endl;
   std::cout << "Total recoil energy deposited:  "<< toterecoil*time_interval<< std::endl;

   outfile.close();
   phoutfile.close();

   // Interpolate TGraph of quenched differential spectrum to get evenly spaced Eee bins
   int nquenched = iq;
   TGraph* quenchedspec = new TGraph(nquenched,Eee,dNdEee);

   std::ofstream phoutfile2;
   std::string phoutfilename2;
   phoutfilename2 = "supernova_diff_rates-"+material+"-"+std::string(ffname)+"-photons2.out";
   phoutfile2.open(phoutfilename2);

   double totquencheden=0;
   double enee;
   double eneestep=0.0001;
   for (enee=0;enee<=Eee[nquenched-1];enee+=eneestep) {
     std::cout << enee << " "<<quenchedspec->Eval(enee)<<std::endl;
     phoutfile2 << enee << " "<<quenchedspec->Eval(enee)<<std::endl;
     totquencheden += enee*eneestep*quenchedspec->Eval(enee);
     
   }

   std::cout<< "Total quenched energy deposited in MeV "<<totquencheden*time_interval<<" photons: "<<totquencheden*24000.*time_interval<<std::endl;
   phoutfile2.close();
   

   return 0;

}


