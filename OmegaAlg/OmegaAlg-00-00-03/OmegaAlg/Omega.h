#ifndef Physics_Analysis_Omega_H
#define Physics_Analysis_Omega_H

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"

class Omega : public Algorithm
{
public:
  Omega(const std::string &name, ISvcLocator *pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private:
  double m_energy;
  // tree for "truth"
  NTuple::Tuple *m_tuple1;
  // isr pi+ pi-
  NTuple::Item<double> t1_misr, t1_aisr, t1_pisr, t1_eisr, t1_pxisr, t1_pyisr, t1_pzisr;
  NTuple::Item<double> t1_mpip, t1_apip, t1_ppip, t1_epip, t1_pxpip, t1_pypip, t1_pzpip;
  NTuple::Item<double> t1_mpim, t1_apim, t1_ppim, t1_epim, t1_pxpim, t1_pypim, t1_pzpim;
  // pi01 pi02 pi03 w wpi02 wpi03 pi02pi03
  NTuple::Item<double> t1_mpi01, t1_api01, t1_ppi01, t1_epi01, t1_pxpi01, t1_pypi01, t1_pzpi01;
  NTuple::Item<double> t1_mpi02, t1_api02, t1_ppi02, t1_epi02, t1_pxpi02, t1_pypi02, t1_pzpi02;
  NTuple::Item<double> t1_mpi03, t1_api03, t1_ppi03, t1_epi03, t1_pxpi03, t1_pypi03, t1_pzpi03;
  NTuple::Item<double> t1_momega, t1_aomega, t1_pomega, t1_eomega, t1_pxomega, t1_pyomega, t1_pzomega;
  NTuple::Item<double> t1_momegapi02, t1_aomegapi02, t1_pomegapi02;
  NTuple::Item<double> t1_momegapi03, t1_aomegapi03, t1_pomegapi03;
  NTuple::Item<double> t1_mpi02pi03, t1_api02pi03, t1_ppi02pi03;

  // tree for "fit4c"
  NTuple::Tuple *m_tuple4;
  NTuple::Item<double> t4_chisq_4c;
  NTuple::Item<double> t4_chisq_3pi;
  // pi+ pi-
  NTuple::Item<double> t4_mpip, t4_apip, t4_ppip, t4_epip, t4_pxpip, t4_pypip, t4_pzpip;
  NTuple::Item<double> t4_mpim, t4_apim, t4_ppim, t4_epim, t4_pxpim, t4_pypim, t4_pzpim;
  // pi01 pi02 pi03 w wpi02 wpi03 pi02pi03
  NTuple::Item<double> t4_mpi01, t4_api01, t4_ppi01, t4_epi01, t4_pxpi01, t4_pypi01, t4_pzpi01;
  NTuple::Item<double> t4_mpi02, t4_api02, t4_ppi02, t4_epi02, t4_pxpi02, t4_pypi02, t4_pzpi02;
  NTuple::Item<double> t4_mpi03, t4_api03, t4_ppi03, t4_epi03, t4_pxpi03, t4_pypi03, t4_pzpi03;
  NTuple::Item<double> t4_momega, t4_aomega, t4_pomega, t4_eomega, t4_pxomega, t4_pyomega, t4_pzomega;
  NTuple::Item<double> t4_momegapi02, t4_aomegapi02, t4_pomegapi02;
  NTuple::Item<double> t4_momegapi03, t4_aomegapi03, t4_pomegapi03;
  NTuple::Item<double> t4_mpi02pi03, t4_api02pi03, t4_ppi02pi03;
  // topo信息
  NTuple::Item<int> runID;
  NTuple::Item<int> eventID;
  NTuple::Item<int> m_idxmc;
  NTuple::Array<int> m_pdgid;
  NTuple::Array<int> m_motheridx;
};

#endif
