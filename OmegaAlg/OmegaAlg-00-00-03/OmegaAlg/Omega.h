#ifndef Physics_Analysis_Omega_H
#define Physics_Analysis_Omega_H

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
//#include "VertexFit/ReadBeamParFromDb.h"

class Omega : public Algorithm
{
public:
  Omega(const std::string &name, ISvcLocator *pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private:
  double m_vr0cut;
  double m_vz0cut;
  double m_energyThreshold;
  double m_gammaPhiCut;
  double m_gammaThetaCut;
  double m_gammaAngleCut;
  int m_test4C;
  int m_test5C;
  int m_checkDedx;
  int m_checkTof;
  double m_energy;

  Ntuple::Tuple *m_tuplet; // Topology
  NTuple::Item<int> runID;
  NTuple::Item<int> eventID;
  NTuple::Item<int> m_idxmc;
  NTuple::Array<int> m_pdgid;
  NTuple::Array<int> m_motheridx;

  NTuple::Tuple *m_tuple4; // Omega 4C
  NTuple::Item<double> m_chisq_fit;
  NTuple::Item<double> m_omega_fit;
  NTuple::Item<double> m_pi01_fit;
  NTuple::Item<double> m_pi02_fit;
  NTuple::Item<double> m_pi03_fit;

  NTuple::Tuple *m_tuple5; // Omega 5C
  NTuple::Item<double> m_chisq;
  NTuple::Item<double> m_mpi01;
  NTuple::Item<double> m_mpi02;
  NTuple::Item<double> m_mpi03;
  NTuple::Item<double> m_momega;
};

#endif
