#ifndef Physics_Analysis_Omega_H
#define Physics_Analysis_Omega_H

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"

class Pppmpz : public Algorithm
{
public:
    Pppmpz(const std::string &name, ISvcLocator *pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

private:
    // ****************************************
    double m_energy;
    // ****************************************
    /*
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
    */
    // ****************************************
    // TREE - "fit4c"
    NTuple::Tuple *m_tuple4;
    // chisq
    NTuple::Item<double> fit4c_chisq_4c;
    // pip pim
    NTuple::Item<double> fit4c_mpip, fit4c_apip, fit4c_ppip;
    NTuple::Item<double> fit4c_mpim, fit4c_apim, fit4c_ppim;
    NTuple::Item<double> fit4c_mpiz, fit4c_apiz, fit4c_ppiz;
    NTuple::Item<double> fit4c_mpipz, fit4c_apipz, fit4c_ppipz;
    NTuple::Item<double> fit4c_mpimz, fit4c_apimz, fit4c_ppimz;
    NTuple::Item<double> fit4c_mpipm, fit4c_apipm, fit4c_ppipm;

    // topo信息
    NTuple::Item<int> runID;
    NTuple::Item<int> eventID;
    NTuple::Item<int> m_idxmc;
    NTuple::Array<int> m_pdgid;
    NTuple::Array<int> m_motheridx;
};

#endif
