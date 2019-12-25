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
  double m_energy;

  // tree for "truth"
  NTuple::Tuple *m_tuple1;
  // gamma信息
  NTuple::Item<double> t1_misr;
  NTuple::Item<double> t1_aisr;
  NTuple::Item<double> t1_pisr;
  NTuple::Item<double> t1_mpip;
  NTuple::Item<double> t1_apip;
  NTuple::Item<double> t1_ppip;
  NTuple::Item<double> t1_mpim;
  NTuple::Item<double> t1_apim;
  NTuple::Item<double> t1_ppim;
  NTuple::Item<double> t1_mgamma1;
  NTuple::Item<double> t1_agamma1;
  NTuple::Item<double> t1_pgamma1;
  NTuple::Item<double> t1_mgamma2;
  NTuple::Item<double> t1_agamma2;
  NTuple::Item<double> t1_pgamma2;
  NTuple::Item<double> t1_mgamma3;
  NTuple::Item<double> t1_agamma3;
  NTuple::Item<double> t1_pgamma3;
  NTuple::Item<double> t1_mgamma4;
  NTuple::Item<double> t1_agamma4;
  NTuple::Item<double> t1_pgamma4;
  NTuple::Item<double> t1_mgamma5;
  NTuple::Item<double> t1_agamma5;
  NTuple::Item<double> t1_pgamma5;
  NTuple::Item<double> t1_mgamma6;
  NTuple::Item<double> t1_agamma6;
  NTuple::Item<double> t1_pgamma6;
  // 粒子信息
  NTuple::Item<double> t1_mpi01;
  NTuple::Item<double> t1_api01;
  NTuple::Item<double> t1_ppi01;
  NTuple::Item<double> t1_mpi02;
  NTuple::Item<double> t1_api02;
  NTuple::Item<double> t1_ppi02;
  NTuple::Item<double> t1_mpi03;
  NTuple::Item<double> t1_api03;
  NTuple::Item<double> t1_ppi03;
  NTuple::Item<double> t1_momega;
  NTuple::Item<double> t1_aomega;
  NTuple::Item<double> t1_pomega;
  NTuple::Item<double> t1_momegapi02;
  NTuple::Item<double> t1_aomegapi02;
  NTuple::Item<double> t1_pomegapi02;
  NTuple::Item<double> t1_momegapi03;
  NTuple::Item<double> t1_aomegapi03;
  NTuple::Item<double> t1_pomegapi03;
  NTuple::Item<double> t1_mpi02pi03;
  NTuple::Item<double> t1_api02pi03;
  NTuple::Item<double> t1_ppi02pi03;

  // tree for "fit4c"
  NTuple::Tuple *m_tuple4;
  NTuple::Item<double> t4_chisq_4c;
  NTuple::Item<double> t4_chisq_3pi;
  // gamma信息
  NTuple::Item<double> t4_mpip;
  NTuple::Item<double> t4_apip;
  NTuple::Item<double> t4_ppip;
  NTuple::Item<double> t4_mpim;
  NTuple::Item<double> t4_apim;
  NTuple::Item<double> t4_ppim;
  NTuple::Item<double> t4_mgamma1;
  NTuple::Item<double> t4_agamma1;
  NTuple::Item<double> t4_pgamma1;
  NTuple::Item<double> t4_mgamma2;
  NTuple::Item<double> t4_agamma2;
  NTuple::Item<double> t4_pgamma2;
  NTuple::Item<double> t4_mgamma3;
  NTuple::Item<double> t4_agamma3;
  NTuple::Item<double> t4_pgamma3;
  NTuple::Item<double> t4_mgamma4;
  NTuple::Item<double> t4_agamma4;
  NTuple::Item<double> t4_pgamma4;
  NTuple::Item<double> t4_mgamma5;
  NTuple::Item<double> t4_agamma5;
  NTuple::Item<double> t4_pgamma5;
  NTuple::Item<double> t4_mgamma6;
  NTuple::Item<double> t4_agamma6;
  NTuple::Item<double> t4_pgamma6;
  // 粒子信息
  NTuple::Item<double> t4_mpi01;
  NTuple::Item<double> t4_api01;
  NTuple::Item<double> t4_ppi01;
  NTuple::Item<double> t4_mpi02;
  NTuple::Item<double> t4_api02;
  NTuple::Item<double> t4_ppi02;
  NTuple::Item<double> t4_mpi03;
  NTuple::Item<double> t4_api03;
  NTuple::Item<double> t4_ppi03;
  NTuple::Item<double> t4_momega;
  NTuple::Item<double> t4_aomega;
  NTuple::Item<double> t4_pomega;
  NTuple::Item<double> t4_momegapi02;
  NTuple::Item<double> t4_aomegapi02;
  NTuple::Item<double> t4_pomegapi02;
  NTuple::Item<double> t4_momegapi03;
  NTuple::Item<double> t4_aomegapi03;
  NTuple::Item<double> t4_pomegapi03;
  NTuple::Item<double> t4_mpi02pi03;
  NTuple::Item<double> t4_api02pi03;
  NTuple::Item<double> t4_ppi02pi03;
  // topo信息
  NTuple::Item<int> runID;
  NTuple::Item<int> eventID;
  NTuple::Item<int> m_idxmc;
  NTuple::Array<int> m_pdgid;
  NTuple::Array<int> m_motheridx;
};

#endif
