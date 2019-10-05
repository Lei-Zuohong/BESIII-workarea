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
  int m_test4C;
  int m_test5C;
  int m_checkDedx;
  int m_checkTof;
  double m_energy;

  NTuple::Tuple *m_tuple4; // Omega 4C
  NTuple::Item<double> m_chisq_4c;
  NTuple::Item<double> m_chisq_3pi;
  NTuple::Item<double> m_omega_4c;
  NTuple::Item<double> m_pi01_4c;
  NTuple::Item<double> m_pi02_4c;
  NTuple::Item<double> m_pi03_4c;
  NTuple::Item<double> m_f0_4c;
  NTuple::Item<double> m_b0_4c;

  NTuple::Item<int> runID;
  NTuple::Item<int> eventID;
  NTuple::Item<int> m_idxmc;
  NTuple::Array<int> m_pdgid;
  NTuple::Array<int> m_motheridx;

  NTuple::Tuple *m_tuple5; // Omega 5C
  NTuple::Item<double> m_chisq_5c;
  NTuple::Item<double> m_omega_5c;
  NTuple::Item<double> m_pi01_5c;
  NTuple::Item<double> m_pi02_5c;
  NTuple::Item<double> m_pi03_5c;
  NTuple::Item<double> m_f0_5c;
  NTuple::Item<double> m_b0_5c;
};

#endif
