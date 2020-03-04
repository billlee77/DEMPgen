/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Asymmetry Class.
  Loads asymmetry parameterization into TF2 object.
  Persistent between events to minimize disk I/O.
  Handles fitting of Monte-Carlo asymmetry data,
  extrapolation between fits.
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef Asymmetry_H
#define Asymmetry_H

#include <vector>
#include "TF1.h"

using namespace std;

class Asymmetry
{
private:

  char * root_filepath;     // Name of .root file with Raws and/or
                            // parameters
  char * AsyNameStr;        // Name of asymmetry
                            // (must match branch name)

  int nQsq;                 // Number of Qsq values
  vector<double> Qsq_Vec;   // Qsq values
  vector<TF1*> AsyFunction; // Asymmetry Functions

  char * FuncForm;          // General form of function to fit data

  int nPars;                // Number of pars in fitting function


  double Extrap(double x0, double x1, double x2,
                double y1, double y2);
                            // Extrapolation function

public:
  Asymmetry(char * in_AsyName, char * in_Func,
            vector<double> in_Qsq = vector<double>(),
            bool refit = true);

  int SetPars(vector<double> in_Qsq);            // Set Pars from file

  int Parameterize();       // Set pars from fitting

  int Parameterize(vector<double> in_Qsq);
                            // Use only specified Qsq values

  double GetAsyAmp(double Qsq, double tp);
                            // Return asymmetry amplitude

  // void SetFunc(const string& fitfunc, int n);
                            // Sets the fitting function from string

  void PrintPars();

};
#endif
