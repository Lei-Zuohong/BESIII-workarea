#ifndef Physics_Analysis_Pm_H
#define Physics_Analysis_Pm_H

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
  double m_energy;
  // tree for "truth"
  NTuple::Tuple *m_tuple1;
  // isr
  NTuple::Item<double> t1_misr, t1_aisr, t1_pisr, t1_eisr, t1_pxisr, t1_pyisr, t1_pzisr;
  // pip2 pim2
  NTuple::Item<double> t1_mpip2, t1_apip2, t1_ppip2, t1_epip2, t1_pxpip2, t1_pypip2, t1_pzpip2;
  NTuple::Item<double> t1_mpim2, t1_apim2, t1_ppim2, t1_epim2, t1_pxpim2, t1_pypim2, t1_pzpim2;
  // pi0 pip1 pim1 w wpip2 wpim2 pip2pip3
  NTuple::Item<double> t1_mpi0, t1_api0, t1_ppi0, t1_epi0, t1_pxpi0, t1_pypi0, t1_pzpi0;
  NTuple::Item<double> t1_mpip1, t1_apip1, t1_ppip1, t1_epip1, t1_pxpip1, t1_pypip1, t1_pzpip1;
  NTuple::Item<double> t1_mpim1, t1_apim1, t1_ppim1, t1_epim1, t1_pxpim1, t1_pypim1, t1_pzpim1;
  NTuple::Item<double> t1_momega, t1_aomega, t1_pomega, t1_eomega, t1_pxomega, t1_pyomega, t1_pzomega;
  NTuple::Item<double> t1_momegapip2, t1_aomegapip2, t1_pomegapip2;
  NTuple::Item<double> t1_momegapim2, t1_aomegapim2, t1_pomegapim2;
  NTuple::Item<double> t1_mpip2pim2, t1_apip2pim2, t1_ppip2pim2;
  // tree for "fit4c"
  NTuple::Tuple *m_tuple4;
  // chisq
  NTuple::Item<double> t4_chisq_4c;
  // pip2 pim2
  NTuple::Item<double> t4_mpip2, t4_apip2, t4_ppip2, t4_epip2, t4_pxpip2, t4_pypip2, t4_pzpip2;
  NTuple::Item<double> t4_mpim2, t4_apim2, t4_ppim2, t4_epim2, t4_pxpim2, t4_pypim2, t4_pzpim2;
  // pi0 pip1 pim1 w wpip2 wpim2 pip2pip3
  NTuple::Item<double> t4_mpi0, t4_api0, t4_ppi0, t4_epi0, t4_pxpi0, t4_pypi0, t4_pzpi0;
  NTuple::Item<double> t4_mpip1, t4_apip1, t4_ppip1, t4_epip1, t4_pxpip1, t4_pypip1, t4_pzpip1;
  NTuple::Item<double> t4_mpim1, t4_apim1, t4_ppim1, t4_epim1, t4_pxpim1, t4_pypim1, t4_pzpim1;
  NTuple::Item<double> t4_momega, t4_aomega, t4_pomega, t4_eomega, t4_pxomega, t4_pyomega, t4_pzomega;
  NTuple::Item<double> t4_momegapip2, t4_aomegapip2, t4_pomegapip2;
  NTuple::Item<double> t4_momegapim2, t4_aomegapim2, t4_pomegapim2;
  NTuple::Item<double> t4_mpip2pim2, t4_apip2pim2, t4_ppip2pim2;
  // topo信息
  NTuple::Item<int> runID;
  NTuple::Item<int> eventID;
  NTuple::Item<int> m_idxmc;
  NTuple::Array<int> m_pdgid;
  NTuple::Array<int> m_motheridx;
};

#endif
