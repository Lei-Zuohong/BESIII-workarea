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

  NTuple::Tuple *m_tuplet; // Topology
  NTuple::Item<int> runID;
  NTuple::Item<int> eventID;
  NTuple::Item<int> m_idxmc;
  NTuple::Array<int> m_pdgid;
  NTuple::Array<int> m_motheridx;

  NTuple::Tuple *m_tuple4; // Omega 4C
  NTuple::Item<double> m_chisq_4c;
  NTuple::Item<double> m_omega_4c;
  NTuple::Item<double> m_pi01_4c;
  NTuple::Item<double> m_pi02_4c;
  NTuple::Item<double> m_pi03_4c;
  NTuple::Item<int> runID_4c;
  NTuple::Item<int> eventID_4c;
  NTuple::Item<int> m_idxmc_4c;
  NTuple::Array<int> m_pdgid_4c;
  NTuple::Array<int> m_motheridx_4c;

  NTuple::Tuple *m_tuple5; // Omega 5C
  NTuple::Item<double> m_chisq_5c;
  NTuple::Item<double> m_mpi01_5c;
  NTuple::Item<double> m_mpi02_5c;
  NTuple::Item<double> m_mpi03_5c;
  NTuple::Item<double> m_momega_5c;
  NTuple::Item<int> runID_5c;
  NTuple::Item<int> eventID_5c;
  NTuple::Item<int> m_idxmc_5c;
  NTuple::Array<int> m_pdgid_5c;
  NTuple::Array<int> m_motheridx_5c;
};

#endif
