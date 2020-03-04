#ifndef SigmaL_H
#define SigmaL_H

#include <vector>

#include "TTree.h"

class SigmaL
{
  TTree* Raw;
  TTree* tPar;
  TTree* Q2Par;

  vector< TF2* > Functions;

  vector< double > WVals;

public:
  void Parameterize();

#endif
