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
    // Parameter from joboption
    double job_energy;
    int job_flag3;

    int job_do_truth;
    int job_do_pid;
    int job_do_4c;
    int job_do_4c_0;
    int job_do_4c_1;
    int job_do_4c_3;
    int job_do_4c_4;
    // ****************************************
    // TREE - "truth"
    NTuple::Tuple *m_tuple_truth;
    NTuple::Item<double> truth_mgamma1, truth_agamma1, truth_pgamma1;
    NTuple::Item<double> truth_mgamma2, truth_agamma2, truth_pgamma2;
    NTuple::Item<double> truth_mpip, truth_apip, truth_ppip;
    NTuple::Item<double> truth_mpim, truth_apim, truth_ppim;
    NTuple::Item<double> truth_mpiz, truth_apiz, truth_ppiz;
    NTuple::Item<double> truth_mpipz, truth_apipz, truth_ppipz;
    NTuple::Item<double> truth_mpimz, truth_apimz, truth_ppimz;
    NTuple::Item<double> truth_mpipm, truth_apipm, truth_ppipm;
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
    // Topology信息
    NTuple::Item<int> runNo;
    NTuple::Item<int> event;
    NTuple::Item<int> m_idxmc;
    NTuple::Array<int> m_pdgid;
    NTuple::Array<int> m_motheridx;
    NTuple::Item<int> flag1;
    NTuple::Item<int> flag2;
    NTuple::Item<int> flag3;
    // For Cut
    NTuple::Item<double> fit4c_chisq;
    NTuple::Item<double> fit4c_chisq_0g;
    NTuple::Item<double> fit4c_chisq_1g;
    NTuple::Item<double> fit4c_chisq_3g;
    NTuple::Item<double> fit4c_chisq_4g;
    // For track
    NTuple::Item<double> fit4c_egamma1;
    NTuple::Item<double> fit4c_egamma2;
    NTuple::Item<double> fit4c_hgamma1;
    NTuple::Item<double> fit4c_hgamma2;
    NTuple::Item<double> fit4c_agamma;
    NTuple::Item<double> fit4c_mpip, fit4c_apip, fit4c_ppip;
    NTuple::Item<double> fit4c_mpim, fit4c_apim, fit4c_ppim;
    NTuple::Item<double> fit4c_mpiz, fit4c_apiz, fit4c_ppiz;
    NTuple::Item<double> fit4c_mpipz, fit4c_apipz, fit4c_ppipz;
    NTuple::Item<double> fit4c_mpimz, fit4c_apimz, fit4c_ppimz;
    NTuple::Item<double> fit4c_mpipm, fit4c_apipm, fit4c_ppipm;
    // ****************************************
    // TREE - "fit5c"
    NTuple::Tuple *m_tuple5;
    NTuple::Item<double> fit5c_chisq;
    NTuple::Item<double> fit5c_mpip, fit5c_apip, fit5c_ppip;
    NTuple::Item<double> fit5c_mpim, fit5c_apim, fit5c_ppim;
    NTuple::Item<double> fit5c_mpiz, fit5c_apiz, fit5c_ppiz;
    NTuple::Item<double> fit5c_mpipz, fit5c_apipz, fit5c_ppipz;
    NTuple::Item<double> fit5c_mpimz, fit5c_apimz, fit5c_ppimz;
    NTuple::Item<double> fit5c_mpipm, fit5c_apipm, fit5c_ppipm;
};

#endif
