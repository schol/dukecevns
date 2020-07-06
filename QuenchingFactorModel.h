#include <cstdlib>
#include <ctime>
#include <cstring>
#include <sstream>
#include <iostream>
#include <math.h>
#include <map>
#include <vector>
class QuenchingFactorModel{
  //abstract base class for all QF models
protected:
  double qfrange[2]; // Range of validity
  char qffilename[80];
public:
  virtual ~QuenchingFactorModel(){};
  virtual double qf(double erec) = 0;
  virtual double qfderiv(double erec) = 0;//function to help convert Enr binning to Eee binning
  virtual void SetQFFilename(const char * fname){
    strcpy(qffilename, fname);
  }
  virtual void ReadQFFile() = 0;
  virtual void debuginfo(){
    std::cout<<"Lindhard QF at 10keV is:"<<qf(0.01)<<" derivative is: "<<qfderiv(0.01)<<std::endl;
  }
};

class QFNumModel : public QuenchingFactorModel{
protected:
  std::map<double,double> _qfmap;

public:
  // Numerical QF file related methods
  // Energies in MeV

  void ReadQFFile(){
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

  double qf(double erec){
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

  double qfderiv(double erec){
    double qfnumderiv = 0.;

    // Compute the derivative of the numerical map at value erec... it's derivative of erec*qf

    typedef std::map<double, double>::const_iterator i_t;

    i_t i=_qfmap.upper_bound(erec);

    double delta;
    double er1,er2;
    double qf1,qf2;
    double rise;
    double run;

    i_t l;
    if(i==_qfmap.end()) {
      // Actually same as normal case
        i_t np = i; --np;
        er1 = np->first;
        er2 = i->first;
        qf1 = np->second;
        qf2 = i->second;

        //      rise = qf2 * er2- qf1 * er1;
        //run = er2-er1;
    }
    else if (i==_qfmap.begin()){
        i_t nl = i; nl++;
        er1 = i->first;
        er2 = nl->first;
        qf1 = i->second;
        qf2 = nl->second;


    } else {
      l=i; --l;
      er1 = l->first;
      er2 = i->first;

      qf1 = l->second;
      qf2 = i->second;


    }

    //double qferec = qf1+(qf2-qf1)/(er2-er1)*(erec-er1);
    // This gives problems if erec=er1
    //  rise = qferec * erec - qf1 * er1;
    rise = qf2 * er2 - qf1 * er1;
    run = er2-er1;

    //  std::cout << er1<<" "<<er2<<" "<<erec<<" "<<qf1<<" "<<qf2<<" "<<qferec<<std::endl;
    delta= rise/run;

    qfnumderiv= delta;

    if (isnan(qfnumderiv)) {qfnumderiv=0.;}

    return qfnumderiv;
  }
};

class QFPolyModel : public QuenchingFactorModel{
protected:
  std::vector<double> qfpolycoeff;

public:
  void ReadQFFile(){
    double coeff;
    std::ifstream qfpolyfile;
    std::string filename = qffilename;
    qfpolyfile.open(filename.c_str());
    if (!qfpolyfile) {
      std::cout << "File "<<filename<<" does not exist!" <<std::endl;
      exit(-1);
    } else {
      qfpolyfile >> qfrange[0]>>qfrange[1];

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

  double qf(double erec){
    // erec in MeV
    double qf=0;
    for (size_t i=0; i<qfpolycoeff.size();i++) {
      if (erec>=qfrange[0]  && erec<=qfrange[1] ) {
        qf += qfpolycoeff[i]*pow(erec,i);
      } else if (erec<qfrange[0]) {
        qf += qfpolycoeff[i]*pow(qfrange[0],i);
      } else {
        qf += qfpolycoeff[i]*pow(qfrange[1],i);

      }
    }

    return qf;
  }

  double qfderiv(double erec){
    // Return the value of the derivative of the polynomial times Erec (assume Eee = qf(Erec)*Erec
    // This useful for binning quenched distributions, dN/dEee = dN/dEr*dEr/dEee

    // erec in MeV
    double qfderiv = 0;
    for (size_t i=0; i<qfpolycoeff.size();i++) {

      if (erec>=qfrange[0]  && erec<=qfrange[1] ) {

        qfderiv += (i+1)*qfpolycoeff[i]*pow(erec,i);
      } else if (erec<qfrange[0]) {
        qfderiv += (i+1)*qfpolycoeff[i]*pow(qfrange[0],i);
      } else {
        qfderiv += (i+1)*qfpolycoeff[i]*pow(qfrange[1],i);

      }

    }

    return qfderiv;
  }
};

class QFLindModel : public QuenchingFactorModel{
protected:
  double k;
  double pre_eps;//this parameter times Erec get epsilon


public:
  void ReadQFFile(){
    double Z,A,f;//atomic number Z,mass number A and f in lindhard k=fZA
    std::ifstream qffile;
    std::string filename = qffilename;
    qffile.open(filename.c_str());
    if (!qffile) {
      std::cout << "File "<<filename<<" does not exist!" <<std::endl;
      exit(-1);
    }
    else {
      qffile>>Z>>A>>f; //order in qf files
    }
    k=f*pow(Z,2.0/3)*pow(A,-0.5);
    pre_eps=11.5*1000*pow(Z,-7.0/3);//Erec unit is in MeV
    qffile.close();
  }

  double qf(double erec){
    double eps=pre_eps*erec;
    double ge=3*pow(eps,0.15)+0.7*pow(eps,0.6)+eps;
    return k*ge/(1+k*ge);
  }

  double qfderiv(double erec){
    //dEee/dErec=qf+Erec*dqf/dErec=qf+Erec*(1+kge)^-2*dkge/deps*deps/dErec
    double eps=pre_eps*erec;
    double ge=3*pow(eps,0.15)+0.7*pow(eps,0.6)+eps;
    double ge_deriv=0.45*pow(eps,-0.85)+0.42*pow(eps,-0.4)+1;
    return k*ge/(1+k*ge)+erec*pow(1+k*ge,-2)*k*ge_deriv*pre_eps;
  }

};
