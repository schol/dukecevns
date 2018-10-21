#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>
#include "DetectorResponse.h"

// Numerical QF file related methods
// Energies in MeV

void DetectorResponse::ReadQFFile()
{

  double erec;
  double qf;
  std::ifstream qffile;
  std::string filename = qffilename;
  qffile.open(filename.c_str());
  if (!qffile) {
    std::cout << "File "<<filename<<" does not exist!" <<std::endl;
    exit(-1);
  } else {
    while(! qffile.eof() ) 
      {
        qffile >> erec >> qf;

        if (! qffile.eof()) {
          _qfmap[erec] = qf;
        }
      }
  }

  qffile.close();
}


////

double DetectorResponse::qfnum(double erec) 
{
  double qf = 1;


    //http://www.bnikolic.co.uk/blog/cpp-map-interp.html
    // Interpolate from the map.  Must have been initalized for output to make sense

  typedef std::map<double, double>::const_iterator i_t;

  //  std::map<double, double> _qfmap;
  
  i_t i=_qfmap.upper_bound(erec);

  if(i==_qfmap.end())
    {
      return (--i)->second;
    }
  if (i==_qfmap.begin())
    {
      return i->second;
    }
  i_t l=i; --l;
  
  const double delta=(erec- l->first)/(i->first - l->first);
  qf= delta*i->second +(1-delta)*l->second;

  if (isnan(qf)) {qf=0.;}

  return qf;

}

void DetectorResponse::SetQFFilename(const char * fname) {
  strcpy(qffilename, fname);
}

const char * DetectorResponse::GetQFFilename() {
  return qffilename;
}

double DetectorResponse::maxErec() 
{

  // Return the maximum energy in MeV.   For the numerical file 

  typedef std::map<double, double>::const_reverse_iterator i_t;
  i_t it = _qfmap.rbegin();
  double maxErec = it->first;

  return maxErec;

}

/////// 

// Polynomial QF-related methods

void DetectorResponse::ReadQFPolyFile() {

  double coeff;
  std::ifstream qfpolyfile;
  std::string filename = qfpolyfilename;
  qfpolyfile.open(filename.c_str());
  if (!qfpolyfile) {
    std::cout << "File "<<filename<<" does not exist!" <<std::endl;
    exit(-1);
  } else {
    qfpolyfile >> qfpolyrange[0]>>qfpolyrange[1];

    while(! qfpolyfile.eof() ) 
      {
	qfpolyfile>>coeff;
        if (! qfpolyfile.eof()) {
	  qfpolycoeff.push_back(coeff);
        }
      }
  }

  qfpolyfile.close();

}


double DetectorResponse::qfpoly(double erec) {

  // erec in MeV
  double qf=0;
  if (erec>=qfpolyrange[0]  && erec<=qfpolyrange[1] ) {
    for (int i=0; i<qfpolycoeff.size();i++) {
      qf += qfpolycoeff[i]*pow(erec,i);
    }
  }

  return qf;
} 

double DetectorResponse::qfpolyderiv(double erec) {

  // Return the value of the derivative of the polynomial times Erec (assume Eee = qf(Erec)*Erec
  // This useful for binning quenched distributions, dN/dEee = dN/dEr*dEr/dEee

  // erec in MeV
  double qfderiv = 0;
  if (erec>=qfpolyrange[0]  && erec<=qfpolyrange[1] ) {
    for (int i=0; i<qfpolycoeff.size();i++) {
      qfderiv += (i+1)*qfpolycoeff[i]*pow(erec,i);
    }
  }

  return qfderiv;
} 


void DetectorResponse::SetQFPolyFilename(const char * fname) {
  strcpy(qfpolyfilename, fname);
}

const char * DetectorResponse::GetQFPolyFilename() {
  return qfpolyfilename;
}


void DetectorResponse::SetQFPolyRange(double* range) {

  qfpolyrange[0] = range[0];
  qfpolyrange[1] = range[1];

}

double* DetectorResponse::GetQFPolyRange() { return qfpolyrange;}

// For Gaussian smearing

// Polynomial GS-related methods

void DetectorResponse::ReadGSPolyFile() {

  double coeff;
  std::ifstream gspolyfile;
  std::string filename = gspolyfilename;
  gspolyfile.open(filename.c_str());
  if (!gspolyfile) {
    std::cout << "File "<<filename<<" does not exist!" <<std::endl;
    exit(-1);
  } else {
    gspolyfile >> gspolyrange[0]>>gspolyrange[1];

    while(! gspolyfile.eof() ) 
      {
	gspolyfile>>coeff;
        if (! gspolyfile.eof()) {
	  gspolycoeff.push_back(coeff);
        }
      }
  }

  gspolyfile.close();

}


double DetectorResponse::gspoly(double en) {

  // returns the sigma as a function of energy (usually Eee)

  // erec in MeV
  double gs=0;
  if (en>=gspolyrange[0]  && en<=gspolyrange[1] ) {
    for (int i=0; i<gspolycoeff.size();i++) {
      gs += gspolycoeff[i]*pow(en,i);
    }
  }

  return gs;
} 



void DetectorResponse::SetGSPolyFilename(const char * fname) {
  strcpy(gspolyfilename, fname);
}

const char * DetectorResponse::GetGSPolyFilename() {
  return gspolyfilename;
}


void DetectorResponse::SetGSPolyRange(double* range) {

  gspolyrange[0] = range[0];
  gspolyrange[1] = range[1];

}

double* DetectorResponse::GetGSPolyRange() { return gspolyrange;}



// Generic methods

DetectorResponse::DetectorResponse(){}

DetectorResponse::DetectorResponse(const char * type)
{
   strcpy(detectortype,type);
}


void DetectorResponse::Setdetectortype(const char * type) {
  strcpy(detectortype, type);
}

const char * DetectorResponse::Getdetectortype() {
  return detectortype;
}

// Numerical efficiency file related methods

void DetectorResponse::ReadEfficFile()
{

  double erec;
  double effic;
  std::ifstream efficfile;
  std::string filename = efficfilename;
  efficfile.open(filename.c_str());
  if (!efficfile) {
    std::cout << "File "<<filename<<" does not exist!" <<std::endl;
    exit(-1);
  } else {
    while(! efficfile.eof() ) 
      {
        efficfile >> erec >> effic;

        if (! efficfile.eof()) {
          _efficmap[erec] = effic;
        }
      }
  }

  efficfile.close();
}


////

double DetectorResponse::efficnum(double erec) 
{
  double effic = 1;


    //http://www.bnikolic.co.uk/blog/cpp-map-interp.html
    // Interpolate from the map.  Must have been initalized for output to make sense

  typedef std::map<double, double>::const_iterator i_t;

  //  std::map<double, double> _efficmap;
  
  i_t i=_efficmap.upper_bound(erec);

  if(i==_efficmap.end())
    {
      return (--i)->second;
    }
  if (i==_efficmap.begin())
    {
      return i->second;
    }
  i_t l=i; --l;
  
  const double delta=(erec- l->first)/(i->first - l->first);
  effic= delta*i->second +(1-delta)*l->second;

  if (isnan(effic)) {effic=0.;}

  return effic;

}

void DetectorResponse::SetEfficFilename(const char * fname) {
  strcpy(efficfilename, fname);
}

const char * DetectorResponse::GetEfficFilename() {
  return efficfilename;
}

double DetectorResponse::maxEfficErec() 
{

  // Return the maximum energy in MeV.   For the numerical file 

  typedef std::map<double, double>::const_reverse_iterator i_t;
  i_t it = _efficmap.rbegin();
  double maxErec = it->first;

  return maxErec;

}


void DetectorResponse::SetStepThresh(double thresh) {

  step_thresh = thresh;
 
}

double DetectorResponse::GetStepThresh() { return step_thresh;}


/////// 
