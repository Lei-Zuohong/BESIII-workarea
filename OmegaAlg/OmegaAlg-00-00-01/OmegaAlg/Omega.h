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
    // ************************************************************************************************************************
    double m_energy;
    int Do_compare_6_57;
    // ************************************************************************************************************************
    // topo信息
    NTuple::Item<int> runID;
    NTuple::Item<int> eventID;
    NTuple::Item<int> m_idxmc;
    NTuple::Array<int> m_pdgid;
    NTuple::Array<int> m_motheridx;
    // ************************************************************************************************************************
    // TREE - "truth"
    NTuple::Tuple *m_tuple1;
    // isr
    NTuple::Item<double> truth_misr, truth_aisr, truth_pisr, truth_eisr, truth_pxisr, truth_pyisr, truth_pzisr;
    // pip pim pi01~3 omega
    NTuple::Item<double> truth_mpip, truth_apip, truth_ppip, truth_epip, truth_pxpip, truth_pypip, truth_pzpip;
    NTuple::Item<double> truth_mpim, truth_apim, truth_ppim, truth_epim, truth_pxpim, truth_pypim, truth_pzpim;
    NTuple::Item<double> truth_mpi01, truth_api01, truth_ppi01, truth_epi01, truth_pxpi01, truth_pypi01, truth_pzpi01;
    NTuple::Item<double> truth_mpi02, truth_api02, truth_ppi02, truth_epi02, truth_pxpi02, truth_pypi02, truth_pzpi02;
    NTuple::Item<double> truth_mpi03, truth_api03, truth_ppi03, truth_epi03, truth_pxpi03, truth_pypi03, truth_pzpi03;
    NTuple::Item<double> truth_momega, truth_aomega, truth_pomega, truth_eomega, truth_pxomega, truth_pyomega, truth_pzomega;
    // omegapi02 omegapi03 pi02pi03 3pi0
    NTuple::Item<double> truth_momegapi02, truth_aomegapi02, truth_pomegapi02;
    NTuple::Item<double> truth_momegapi03, truth_aomegapi03, truth_pomegapi03;
    NTuple::Item<double> truth_mpi02pi03, truth_api02pi03, truth_ppi02pi03;
    NTuple::Item<double> truth_m3pi0, truth_a3pi0, truth_p3pi0;
    // ************************************************************************************************************************
    // TREE - "charge"
    NTuple::Tuple *m_tuple2;
    NTuple::Item<double> charge_ngood;
    NTuple::Item<double> charge_ncharge;
    // ************************************************************************************************************************
    // TREE - "vertex"
    NTuple::Tuple *m_tuple3;
    NTuple::Item<double> vertex_chisq;
    // ************************************************************************************************************************
    // TREE - "fit4c"
    NTuple::Tuple *m_tuple4;
    // chisq
    NTuple::Item<double> fit4c_chisq_6g;
    // pip pim pi01~3 omega
    NTuple::Item<double> fit4c_mpip, fit4c_apip, fit4c_ppip, fit4c_epip, fit4c_pxpip, fit4c_pypip, fit4c_pzpip;
    NTuple::Item<double> fit4c_mpim, fit4c_apim, fit4c_ppim, fit4c_epim, fit4c_pxpim, fit4c_pypim, fit4c_pzpim;
    NTuple::Item<double> fit4c_mpi01, fit4c_api01, fit4c_ppi01, fit4c_epi01, fit4c_pxpi01, fit4c_pypi01, fit4c_pzpi01;
    NTuple::Item<double> fit4c_mpi02, fit4c_api02, fit4c_ppi02, fit4c_epi02, fit4c_pxpi02, fit4c_pypi02, fit4c_pzpi02;
    NTuple::Item<double> fit4c_mpi03, fit4c_api03, fit4c_ppi03, fit4c_epi03, fit4c_pxpi03, fit4c_pypi03, fit4c_pzpi03;
    NTuple::Item<double> fit4c_momega, fit4c_aomega, fit4c_pomega, fit4c_eomega, fit4c_pxomega, fit4c_pyomega, fit4c_pzomega;
    // omegapi02 omegapi03 pi02pi03 3pi0
    NTuple::Item<double> fit4c_momegapi02, fit4c_aomegapi02, fit4c_pomegapi02;
    NTuple::Item<double> fit4c_momegapi03, fit4c_aomegapi03, fit4c_pomegapi03;
    NTuple::Item<double> fit4c_mpi02pi03, fit4c_api02pi03, fit4c_ppi02pi03;
    NTuple::Item<double> fit4c_m3pi0, fit4c_a3pi0, fit4c_p3pi0;
    // combination
    NTuple::Item<double> fit4c_acombination;
    NTuple::Item<double> fit4c_nphoton;
    NTuple::Item<double> fit4c_eloss;
    // ************************************************************************************************************************
};
#endif
