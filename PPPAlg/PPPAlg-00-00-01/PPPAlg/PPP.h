#ifndef Physics_Analysis_Omega_H
#define Physics_Analysis_Omega_H

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"

class PPP : public Algorithm
{
public:
    PPP(const std::string &name, ISvcLocator *pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

private:
    // ****************************************
    double m_energy;
    // ****************************************
    // TREE - "truth"
    NTuple::Tuple *m_tuple1;
    // isr
    NTuple::Item<double> truth_misr, truth_aisr, truth_pisr, truth_eisr, truth_pxisr, truth_pyisr, truth_pzisr;
    // pip pim
    NTuple::Item<double> truth_mpip, truth_apip, truth_ppip, truth_epip, truth_pxpip, truth_pypip, truth_pzpip;
    NTuple::Item<double> truth_mpim, truth_apim, truth_ppim, truth_epim, truth_pxpim, truth_pypim, truth_pzpim;
    // gamma
    NTuple::Item<double> truth_mgamma1, truth_agamma1, truth_pgamma1, truth_egamma1, truth_pxgamma1, truth_pygamma1, truth_pzgamma1;
    NTuple::Item<double> truth_mgamma2, truth_agamma2, truth_pgamma2, truth_egamma2, truth_pxgamma2, truth_pygamma2, truth_pzgamma2;
    NTuple::Item<double> truth_mgamma3, truth_agamma3, truth_pgamma3, truth_egamma3, truth_pxgamma3, truth_pygamma3, truth_pzgamma3;
    NTuple::Item<double> truth_mgamma4, truth_agamma4, truth_pgamma4, truth_egamma4, truth_pxgamma4, truth_pygamma4, truth_pzgamma4;
    NTuple::Item<double> truth_mgamma5, truth_agamma5, truth_pgamma5, truth_egamma5, truth_pxgamma5, truth_pygamma5, truth_pzgamma5;
    NTuple::Item<double> truth_mgamma6, truth_agamma6, truth_pgamma6, truth_egamma6, truth_pxgamma6, truth_pygamma6, truth_pzgamma6;
    // pi01 pi02 pi03 w wpi02 wpi03 pi02pi03
    NTuple::Item<double> truth_mpi01, truth_api01, truth_ppi01, truth_epi01, truth_pxpi01, truth_pypi01, truth_pzpi01;
    NTuple::Item<double> truth_mpi02, truth_api02, truth_ppi02, truth_epi02, truth_pxpi02, truth_pypi02, truth_pzpi02;
    NTuple::Item<double> truth_mpi03, truth_api03, truth_ppi03, truth_epi03, truth_pxpi03, truth_pypi03, truth_pzpi03;
    NTuple::Item<double> truth_momega, truth_aomega, truth_pomega, truth_eomega, truth_pxomega, truth_pyomega, truth_pzomega;
    NTuple::Item<double> truth_momegapi02, truth_aomegapi02, truth_pomegapi02;
    NTuple::Item<double> truth_momegapi03, truth_aomegapi03, truth_pomegapi03;
    NTuple::Item<double> truth_mpi02pi03, truth_api02pi03, truth_ppi02pi03;
    NTuple::Item<double> truth_m3pi0;
    // hgamma
    NTuple::Item<double> truth_hgamma1, truth_hgamma2, truth_hgamma3, truth_hgamma4, truth_hgamma5, truth_hgamma6;
    // ****************************************
    // TREE - "charge"
    NTuple::Tuple *m_tuple2;
    NTuple::Item<double> charge_ngood;
    NTuple::Item<double> charge_ncharge;
    // ****************************************
    // TREE - "vertex"
    NTuple::Tuple *m_tuple3;
    NTuple::Item<double> vertex_chisq;
    // ****************************************
    // TREE - "fit4c"
    NTuple::Tuple *m_tuple4;
    // chisq
    NTuple::Item<double> fit4c_chisq_4c;
    NTuple::Item<double> fit4c_chisq_3pi;
    // pip pim
    NTuple::Item<double> fit4c_mpip, fit4c_apip, fit4c_ppip, fit4c_epip, fit4c_pxpip, fit4c_pypip, fit4c_pzpip;
    NTuple::Item<double> fit4c_mpim, fit4c_apim, fit4c_ppim, fit4c_epim, fit4c_pxpim, fit4c_pypim, fit4c_pzpim;
    // gamma
    NTuple::Item<double> fit4c_mgamma1, fit4c_agamma1, fit4c_pgamma1, fit4c_egamma1, fit4c_pxgamma1, fit4c_pygamma1, fit4c_pzgamma1;
    NTuple::Item<double> fit4c_mgamma2, fit4c_agamma2, fit4c_pgamma2, fit4c_egamma2, fit4c_pxgamma2, fit4c_pygamma2, fit4c_pzgamma2;
    NTuple::Item<double> fit4c_mgamma3, fit4c_agamma3, fit4c_pgamma3, fit4c_egamma3, fit4c_pxgamma3, fit4c_pygamma3, fit4c_pzgamma3;
    NTuple::Item<double> fit4c_mgamma4, fit4c_agamma4, fit4c_pgamma4, fit4c_egamma4, fit4c_pxgamma4, fit4c_pygamma4, fit4c_pzgamma4;
    NTuple::Item<double> fit4c_mgamma5, fit4c_agamma5, fit4c_pgamma5, fit4c_egamma5, fit4c_pxgamma5, fit4c_pygamma5, fit4c_pzgamma5;
    NTuple::Item<double> fit4c_mgamma6, fit4c_agamma6, fit4c_pgamma6, fit4c_egamma6, fit4c_pxgamma6, fit4c_pygamma6, fit4c_pzgamma6;
    // pi01 pi02 pi03 w wpi02 wpi03 pi02pi03
    NTuple::Item<double> fit4c_mpi01, fit4c_api01, fit4c_ppi01, fit4c_epi01, fit4c_pxpi01, fit4c_pypi01, fit4c_pzpi01;
    NTuple::Item<double> fit4c_mpi02, fit4c_api02, fit4c_ppi02, fit4c_epi02, fit4c_pxpi02, fit4c_pypi02, fit4c_pzpi02;
    NTuple::Item<double> fit4c_mpi03, fit4c_api03, fit4c_ppi03, fit4c_epi03, fit4c_pxpi03, fit4c_pypi03, fit4c_pzpi03;
    NTuple::Item<double> fit4c_momega, fit4c_aomega, fit4c_pomega, fit4c_eomega, fit4c_pxomega, fit4c_pyomega, fit4c_pzomega;
    NTuple::Item<double> fit4c_momegapi02, fit4c_aomegapi02, fit4c_pomegapi02;
    NTuple::Item<double> fit4c_momegapi03, fit4c_aomegapi03, fit4c_pomegapi03;
    NTuple::Item<double> fit4c_mpi02pi03, fit4c_api02pi03, fit4c_ppi02pi03;
    NTuple::Item<double> fit4c_m3pi0;
    // hgamma
    NTuple::Item<double> fit4c_hgamma1, fit4c_hgamma2, fit4c_hgamma3, fit4c_hgamma4, fit4c_hgamma5, fit4c_hgamma6;
    // combination
    NTuple::Item<double> fit4c_acombination;
    NTuple::Item<double> fit4c_nphoton;

    
    // topo信息
    NTuple::Item<int> runID;
    NTuple::Item<int> eventID;
    NTuple::Item<int> m_idxmc;
    NTuple::Array<int> m_pdgid;
    NTuple::Array<int> m_motheridx;
};

#endif
