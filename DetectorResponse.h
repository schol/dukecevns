#ifndef _DetectorResponse_
#define _DetectorResponse_

#include <map>
#include <vector>
#include <fstream>
#include <math.h>
#include <string>

// Energies in MeV
class DetectorResponse
{

 protected:

  char detectortype[80];

  // For QF in numerical format
  std::map<double,double> _qfmap;
  char qffilename[80];

  // For QF in polynomial format
  char qfpolyfilename[80];
  std::vector<double> polycoeff;
  double polyrange[2]; // Range of validity for polynomial

  // For efficiency in numerical format
  std::map<double,double> _efficmap;
  char efficfilename[80];



 public: 
  DetectorResponse();
  DetectorResponse(const char *);
  ~DetectorResponse(){};

  void Setdetectortype(const char *);
  const char * Getdetectortype();

  // For QF in numerical format
  void SetQFFilename(const char * qffilename);
  const char * GetQFFilename();
  void ReadQFFile();
  double qfnum(double);
  double maxErec();


  // For QF in polynomial format
  void SetPolyRange(double*);
  double* GetPolyRange();

  void SetQFPolyFilename(const char * qfpolyfilename);
  const char * GetQFPolyFilename();
  void ReadQFPolyFile();
  double qfpoly(double);
  double qfpolyderiv(double);

  // For efficiency as a function of Erec, file in numerical format

  void SetEfficFilename(const char * qffilename);
  const char * GetEfficFilename();
  void ReadEfficFile();
  double efficnum(double);
  double maxEfficErec();


};



#endif
