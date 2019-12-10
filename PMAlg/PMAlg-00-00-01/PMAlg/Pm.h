#ifndef Physics_Analysis_Pm_H
#define Physics_Analysis_pm_H

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
//#include "VertexFit/ReadBeamParFromDb.h"

class Pm : public Algorithm
{
public:
  Pm(const std::string &name, ISvcLocator *pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private:
  int m_test4C;
  double m_energy;

  NTuple::Tuple *m_tuple4; // Pm 4C
  NTuple::Item<double> m_chisq_4c;
  NTuple::Item<double> m_pi0;
  NTuple::Item<double> m_omega;
  NTuple::Item<double> m_omegapip;
  NTuple::Item<double> m_omegapim;
  NTuple::Item<double> m_pipm;
};

#endif
