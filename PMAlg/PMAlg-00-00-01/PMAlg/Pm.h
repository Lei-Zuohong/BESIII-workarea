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
  NTuple::Item<double> m_chisq_3pi;
  NTuple::Item<double> m_pi01;
  NTuple::Item<double> m_pi02;
  NTuple::Item<double> m_pi03;
  NTuple::Item<double> m_Pm;
  NTuple::Item<double> m_Pmpi02;
  NTuple::Item<double> m_Pmpi03;
  NTuple::Item<double> m_pi02pi03;
  NTuple::Item<double> m_pi01pi02;
  NTuple::Item<double> m_pi01pi03;
};

#endif
