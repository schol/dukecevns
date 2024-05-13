#include <iostream>
#include "TFile.h"
#include "TMath.h"
#include "TH2D.h"
#include "TString.h"
#include <string>
#include <fstream>


////////

bool fileExists(const TString& filename) {
    std::ifstream file(filename.Data());
    return file.good();
}


void get_flavor_weight(Int_t convolved, Double_t start_time, Double_t end_time, Double_t* teffic_params, Double_t* numu_weight, Double_t* numubar_weight, Double_t* nue_weight)
{

// Return the relative flux weighting by flavor for a given time 
// interval

// For use with assumption of perfect stopped pion spectrum
// Ignore nuebar

//  Double_t numu_weight, numubar_weight, nue_weight;

  TString filename;
  if (convolved == 1) {
    filename = "sns_out_BERT_convolved.root";
    if (fileExists(filename)) {
      std::cout<<" Reading file "<<filename<<std::endl;
    } else {
      std::cout<<"File "<<filename<<" does not exist"<<std::endl;
      return;
    }
  } else if (convolved == 2) {
    filename = "sns_out_BERT.root";
    if (fileExists(filename)) {
      std::cout<<" Reading file "<<filename<<std::endl;
    } else {
      std::cout<<"File "<<filename<<" does not exist"<<std::endl;
      return;
    }
  } else {
    std::cout << "Not using flux file "<<std::endl;
    return;
  }

  TFile f(filename);
  TH2D* nue;
  TH2D* nuebar;
  TH2D* numu;
  TH2D* numubar;


  if (convolved == 1) {

    nue = (TH2D*)f.Get("convolved_energy_time_of_nu_e");
    nuebar= (TH2D*)f.Get("convolved_energy_time_of_anti_nu_e");
    numu = (TH2D*)f.Get("convolved_energy_time_of_nu_mu");
    numubar = (TH2D*)f.Get("convolved_energy_time_of_anti_nu_mu");

  } else {


    nue = (TH2D*)f.Get("initial_energy_time_of_nu_e");
    nuebar= (TH2D*)f.Get("initial_energy_time_of_anti_nu_e");
    numu = (TH2D*)f.Get("initial_energy_time_of_nu_mu");
    numubar = (TH2D*)f.Get("initial_energy_time_of_anti_nu_mu");

  }

  Double_t nue_nbinx = nue->GetNbinsX();
  Double_t nue_nbiny = nue->GetNbinsY();

  Int_t ibinx, ibiny, bin;
  Double_t bin_content, bin_center;

  const Int_t maxenvals = 200;
  const Int_t minbinx = 1000;
  const Int_t maxbinx = 20000;
  const Int_t xstep = 1;
  const Int_t minbiny = 0;
  const Int_t maxbiny = 100;

  Double_t xval,yval;

  Double_t tot_nue = 0;  
  Double_t tot_nuebar = 0;  
  Double_t tot_numu = 0;  
  Double_t tot_numubar = 0;  
  Double_t tot_nue_intw = 0;  
  Double_t tot_nuebar_intw = 0;  
  Double_t tot_numu_intw = 0;  
  Double_t tot_numubar_intw = 0;  

  Int_t num_tbins = 0;

  //  nue->Print();

  // Loop over time
  //  for (ibinx=minbinx;ibinx<maxbinx;ibinx+=xstep) {
  for (ibinx=1;ibinx<nue_nbinx;ibinx+=xstep) {
    xval=nue->GetXaxis()->GetBinCenter(ibinx);

    // Loop over energy and sum
    for (ibiny=1;ibiny<nue_nbiny;ibiny++) {

      // Nue

      bin = nue->GetBin(ibinx,ibiny);
      yval=nue->GetYaxis()->GetBinCenter(ibiny);


      // Assume all histos have same binning
      bin_content = nue->GetBinContent(ibinx,ibiny);
      tot_nue += bin_content;

      // Time-dependent efficiency

      Double_t time_effic = 1.;

      Double_t offset = teffic_params[0];
      Double_t a=teffic_params[1];
      Double_t b=teffic_params[2];

      Double_t t = (xval-offset)/1000.; //microsecs
      if (t>a) {
	time_effic = TMath::Exp(-b*(t-a));
      }

       // Check if within time window

      
      if (xval>=start_time && xval<=end_time) { 
	tot_nue_intw +=bin_content*time_effic;
      }


      bin_content = nuebar->GetBinContent(ibinx,ibiny);
      tot_nuebar += bin_content;
      if (xval>=start_time && xval<=end_time) { 
	tot_nuebar_intw +=bin_content*time_effic;
      }


      bin_content = numu->GetBinContent(ibinx,ibiny);
      tot_numu += bin_content;
      if (xval>=start_time && xval<=end_time) { 
	tot_numu_intw +=bin_content*time_effic;
      }

      bin_content = numubar->GetBinContent(ibinx,ibiny);
      tot_numubar += bin_content;
      if (xval>=start_time && xval<=end_time) { 
	tot_numubar_intw +=bin_content*time_effic;
      }

    } // End of energy loop

  } // End of time loop

  // Ignoring nuebar

  
  *nue_weight = tot_nue_intw/tot_nue;
  *numu_weight = tot_numu_intw/tot_numu;
  *numubar_weight = tot_numubar_intw/tot_numubar;

  //  std::cout << "nue " << nue_weight<<" numu_weight "<< numu_weight<<" numubar_weight "<<numubar_weight<<std::endl;


}


