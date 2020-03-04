/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Takes particle and event data, and adds data addresses to to
  a TTree to be saved.
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef TreeBuilder_H
#define TreeBuilder_H

#include "Particle.hxx"
#include "DEMPEvent.hxx"

#include "TTree.h"
#include "TFile.h"

#include <vector>

class TreeBuilder
{
  int nParticles;
  std::vector<Particle*> ParticleList;
  std::vector< std::vector< double* >* > PropList;

  int nVars;

  const char* tree_name;
  const char* file_name;

  TTree * Tree_Out;
  TFile * File_Out;

  void Retrieve();
public:

  TreeBuilder(const char * name);

  void AddParticle(Particle * p);
  void AddDouble(double * x, const char* name);
  void AddEvent(DEMPEvent * event);
  void Fill();
  void Save();
};


#endif
