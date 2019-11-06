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
  double m_energy;

  NTuple::Tuple *m_tuple4; // Omega 4C
  NTuple::Item<double> m_chisq_4c;
  NTuple::Item<double> m_chisq_3pi;
  NTuple::Item<double> m_pi01;
  NTuple::Item<double> m_pi02;
  NTuple::Item<double> m_pi03;
  NTuple::Item<double> m_omega;
  NTuple::Item<double> m_omegapi02;
  NTuple::Item<double> m_omegapi03;
  NTuple::Item<double> m_pi02pi03;
  NTuple::Item<double> m_pi01pi02;
  NTuple::Item<double> m_pi01pi03;
  

  NTuple::Item<int> runID;
  NTuple::Item<int> eventID;
  NTuple::Item<int> m_idxmc;
  NTuple::Array<int> m_pdgid;
  NTuple::Array<int> m_motheridx;
};

#endif
