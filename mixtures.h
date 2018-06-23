// Mixtures.  Assume add up to 100

const int max_components = 15;
double amu = 931.414;
double me=0.511;

  std::map<std::string, std::vector<std::string>> isotopes;
  std::map<std::string, std::vector<double>> molar_fraction;

double mass_fraction[max_components];

  isotopes["C"].push_back("C12");
  molar_fraction["C"].push_back(0.9893);
  isotopes["C"].push_back("C13");
  molar_fraction["C"].push_back(0.0107);


  isotopes["C12"].push_back("C12");
  molar_fraction["C"].push_back(1.0);
  isotopes["C13"].push_back("C13");
  molar_fraction["C13"].push_back(1.0);

  isotopes["Ne"].push_back("Ne20");
  molar_fraction["Ne"].push_back(0.9048);
  isotopes["Ne"].push_back("Ne21");
  molar_fraction["Ne"].push_back(0.0027);
  isotopes["Ne"].push_back("Ne22");
  molar_fraction["Ne"].push_back(0.0925);


  isotopes["Ne20"].push_back("Ne20");
  molar_fraction["Ne20"].push_back(1.0);


  isotopes["CsI"].push_back("Cs133");
  isotopes["CsI"].push_back("I127");

  molar_fraction["CsI"].push_back(0.5);
  molar_fraction["CsI"].push_back(0.5);

  isotopes["Cs133"].push_back("Cs133");
  molar_fraction["Cs133"].push_back(1.0);

  isotopes["I127"].push_back("I127");
  molar_fraction["I127"].push_back(1.0);


////

  isotopes["NaI"].push_back("Na23");
  isotopes["NaI"].push_back("I127");

  molar_fraction["NaI"].push_back(0.5);
  molar_fraction["NaI"].push_back(0.5);

  isotopes["Na23"].push_back("Na23");
  molar_fraction["Na23"].push_back(1.0);


///

  isotopes["Ar"].push_back("Ar36");
  molar_fraction["Ar"].push_back(0.003336);
  isotopes["Ar"].push_back("Ar38");
  molar_fraction["Ar"].push_back(0.000629);
  isotopes["Ar"].push_back("Ar40");
  molar_fraction["Ar"].push_back(0.996035);

// For single-isotope target
  isotopes["Ar36"].push_back("Ar36");
  molar_fraction["Ar36"].push_back(1.0);
  isotopes["Ar38"].push_back("Ar38");
  molar_fraction["Ar38"].push_back(1.0);
  isotopes["Ar40"].push_back("Ar40");
  molar_fraction["Ar40"].push_back(1.0);


//

  isotopes["Ge"].push_back("Ge70");
  molar_fraction["Ge"].push_back(0.2057);
  isotopes["Ge"].push_back("Ge72");
  molar_fraction["Ge"].push_back(0.2745);
  isotopes["Ge"].push_back("Ge73");
  molar_fraction["Ge"].push_back(0.0775);
  isotopes["Ge"].push_back("Ge74");
  molar_fraction["Ge"].push_back(0.3650);
  isotopes["Ge"].push_back("Ge76");
  molar_fraction["Ge"].push_back(0.0773);

  isotopes["Ge72"].push_back("Ge72");
  molar_fraction["Ge72"].push_back(1.);

  isotopes["Ge74"].push_back("Ge74");
  molar_fraction["Ge74"].push_back(1.);

  isotopes["Ge76"].push_back("Ge76");
  molar_fraction["Ge76"].push_back(1.);
//

  isotopes["Xe"].push_back("Xe124");
  molar_fraction["Xe"].push_back(0.000952);
  isotopes["Xe"].push_back("Xe126");
  molar_fraction["Xe"].push_back(0.000890);
  isotopes["Xe"].push_back("Xe128");
  molar_fraction["Xe"].push_back(0.019102);
  isotopes["Xe"].push_back("Xe129");
  molar_fraction["Xe"].push_back(0.264006);
  isotopes["Xe"].push_back("Xe130");
  molar_fraction["Xe"].push_back(0.04071);
  isotopes["Xe"].push_back("Xe131");
  molar_fraction["Xe"].push_back(0.21232);
  isotopes["Xe"].push_back("Xe132");
  molar_fraction["Xe"].push_back(0.269086);
  isotopes["Xe"].push_back("Xe134");
  molar_fraction["Xe"].push_back(0.104357);
  isotopes["Xe"].push_back("Xe136");
  molar_fraction["Xe"].push_back(0.088573);


// For single-isotope target
  isotopes["Xe132"].push_back("Xe132");
  molar_fraction["Xe132"].push_back(1.0);

   



