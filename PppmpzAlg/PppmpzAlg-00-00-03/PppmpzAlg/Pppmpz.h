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
    int Do_Truth;
    int Do_PID;
    int Do_Vertexfit;
    int Do_4c;
    int Do_Compare_2_1;
    int Do_Compare_2_3;
    int Do_1c;
    int Do_5c;
    // ****************************************
    // TREE - "truth"
    NTuple::Tuple *m_tuple1;
    NTuple::Item<double> truth_mpip, truth_apip, truth_ppip;
    NTuple::Item<double> truth_mpim, truth_apim, truth_ppim;
    NTuple::Item<double> truth_mpiz, truth_apiz, truth_ppiz;
    NTuple::Item<double> truth_mpipz, truth_apipz, truth_ppipz;
    NTuple::Item<double> truth_mpimz, truth_apimz, truth_ppimz;
    NTuple::Item<double> truth_mpipm, truth_apipm, truth_ppipm;
    // ****************************************
    // TREE - "charge"
    // NTuple::Tuple *m_tuple2;
    // NTuple::Item<double> charge_ngood;
    // NTuple::Item<double> charge_ncharge;
    // ****************************************
    // TREE - "vertex"
    // NTuple::Tuple *m_tuple3;
    // NTuple::Item<double> vertex_chisq;
    // ****************************************
    // ****************************************
    // topo信息
    NTuple::Item<int> runID;
    NTuple::Item<int> eventID;
    NTuple::Item<int> m_idxmc;
    NTuple::Array<int> m_pdgid;
    NTuple::Array<int> m_motheridx;
    // ****************************************
    // TREE - "fit4c"
    NTuple::Tuple *m_tuple4;
    NTuple::Item<double> fit4c_chisq;
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
    // ****************************************
    // TREE - "fit1c"
    NTuple::Tuple *m_tuple6;
    NTuple::Item<double> fit1c_chisq;
    NTuple::Item<double> fit1c_mpip, fit1c_apip, fit1c_ppip;
    NTuple::Item<double> fit1c_mpim, fit1c_apim, fit1c_ppim;
    NTuple::Item<double> fit1c_mpiz, fit1c_apiz, fit1c_ppiz;
    NTuple::Item<double> fit1c_mpipz, fit1c_apipz, fit1c_ppipz;
    NTuple::Item<double> fit1c_mpimz, fit1c_apimz, fit1c_ppimz;
    NTuple::Item<double> fit1c_mpipm, fit1c_apipm, fit1c_ppipm;
};

#endif
