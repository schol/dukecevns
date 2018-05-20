#ifndef _FormFactor_
#define _FormFactor_

#include <map>
#include <vector>
#include <fstream>
#include <math.h>
#include <string>

using namespace std;

class FormFactor
{

 protected:

  int A; 
  double Rnfac;

  char fftype[80];
 

 public: 
  FormFactor();
  FormFactor(const char *);
  ~FormFactor(){};

  virtual double FFval(double) = 0;

  void SetA(int);
  int GetA();

  // Variation of Rn (as fraction of nominal)
  void SetRnfac(double);
  double GetRnfac();

  void Setfftype(const char *);
  const char * Getfftype();

};

class Helm: public FormFactor {

 protected:
  double sval;

 public:
  Helm() : FormFactor("helm") {}
  double FFval(double);
  void Setsval(double);
  double Getsval(); 

};

class Klein: public FormFactor {

 protected:
  double akval;

 public:
  Klein() : FormFactor("klein") {}
  double FFval(double);
  void Setakval(double);
  double Getakval(); 

};


class Horowitz: public FormFactor {

 protected:
  std::map<double,double> _ffmap;
  char filename[80];

 public:
  Horowitz() : FormFactor("horowitz") {}
  double FFval(double);
  void SetFFfilename(const char * filename);
  const char * GetFFfilename();
  void ReadFFfile();


};

#endif
