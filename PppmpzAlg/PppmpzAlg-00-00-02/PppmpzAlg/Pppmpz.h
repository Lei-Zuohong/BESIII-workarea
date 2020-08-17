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
    int job_do_4c;
    int job_do_4c_0;
    int job_do_4c_1;
    int job_do_4c_3;
    int job_do_4c_4;
    int job_do_4c_e;
    int job_do_4c_k;
    int job_do_4c_m;

    int job_selection;
    // Parameter for topology
    NTuple::Item<int> runNo;
    NTuple::Item<int> event;
    NTuple::Item<int> m_idxmc;
    NTuple::Array<int> m_pdgid;
    NTuple::Array<int> m_motheridx;
    NTuple::Item<int> flag1;
    NTuple::Item<int> flag2;
    NTuple::Item<int> flag3;
    // TREE - "truth"
    NTuple::Tuple *m_tuple_truth;
    NTuple::Item<double> truth_pip_m, truth_pip_p, truth_pip_a, truth_pip_pe, truth_pip_px, truth_pip_py, truth_pip_pz;
    NTuple::Item<double> truth_pim_m, truth_pim_p, truth_pim_a, truth_pim_pe, truth_pim_px, truth_pim_py, truth_pim_pz;
    NTuple::Item<double> truth_gamma1_m, truth_gamma1_p, truth_gamma1_a, truth_gamma1_pe, truth_gamma1_px, truth_gamma1_py, truth_gamma1_pz;
    NTuple::Item<double> truth_gamma2_m, truth_gamma2_p, truth_gamma2_a, truth_gamma2_pe, truth_gamma2_px, truth_gamma2_py, truth_gamma2_pz;
    NTuple::Item<double> truth_piz_m, truth_piz_p, truth_piz_a, truth_piz_pe, truth_piz_px, truth_piz_py, truth_piz_pz;
    NTuple::Item<double> truth_pipm_m, truth_pipm_p, truth_pipm_a, truth_pipm_pe, truth_pipm_px, truth_pipm_py, truth_pipm_pz;
    NTuple::Item<double> truth_pipz_m, truth_pipz_p, truth_pipz_a, truth_pipz_pe, truth_pipz_px, truth_pipz_py, truth_pipz_pz;
    NTuple::Item<double> truth_pimz_m, truth_pimz_p, truth_pimz_a, truth_pimz_pe, truth_pimz_px, truth_pimz_py, truth_pimz_pz;
    // TREE - "charge"
    NTuple::Tuple *m_tuple_charge;
    NTuple::Item<double> charge_ngood;
    NTuple::Item<double> charge_ncharge;
    // ****************************************
    // TREE - "fit4c"
    NTuple::Tuple *m_tuple_fit4c;
    NTuple::Tuple *m_tuple_selection;
    // For Cut
    NTuple::Item<double> fit4c_chisq;
    NTuple::Item<double> fit4c_chisq_0g;
    NTuple::Item<double> fit4c_chisq_1g;
    NTuple::Item<double> fit4c_chisq_3g;
    NTuple::Item<double> fit4c_chisq_4g;
    NTuple::Item<double> fit4c_chisq_e;
    NTuple::Item<double> fit4c_chisq_k;
    NTuple::Item<double> fit4c_chisq_m;
    // For 00 track
    NTuple::Item<double> fit4c_gamma1_heli;
    NTuple::Item<double> fit4c_gamma2_heli;
    NTuple::Item<double> fit4c_a_pippim;
    // For +- track
    NTuple::Item<double> fit4c_pip_ep;
    NTuple::Item<double> fit4c_pim_ep;
    // For track
    NTuple::Item<double> fit4c_pip_m, fit4c_pip_p, fit4c_pip_a, fit4c_pip_pe, fit4c_pip_px, fit4c_pip_py, fit4c_pip_pz;
    NTuple::Item<double> fit4c_pim_m, fit4c_pim_p, fit4c_pim_a, fit4c_pim_pe, fit4c_pim_px, fit4c_pim_py, fit4c_pim_pz;
    NTuple::Item<double> fit4c_gamma1_m, fit4c_gamma1_p, fit4c_gamma1_a, fit4c_gamma1_pe, fit4c_gamma1_px, fit4c_gamma1_py, fit4c_gamma1_pz;
    NTuple::Item<double> fit4c_gamma2_m, fit4c_gamma2_p, fit4c_gamma2_a, fit4c_gamma2_pe, fit4c_gamma2_px, fit4c_gamma2_py, fit4c_gamma2_pz;
    NTuple::Item<double> fit4c_piz_m, fit4c_piz_p, fit4c_piz_a, fit4c_piz_pe, fit4c_piz_px, fit4c_piz_py, fit4c_piz_pz;
    NTuple::Item<double> fit4c_pipm_m, fit4c_pipm_p, fit4c_pipm_a, fit4c_pipm_pe, fit4c_pipm_px, fit4c_pipm_py, fit4c_pipm_pz;
    NTuple::Item<double> fit4c_pipz_m, fit4c_pipz_p, fit4c_pipz_a, fit4c_pipz_pe, fit4c_pipz_px, fit4c_pipz_py, fit4c_pipz_pz;
    NTuple::Item<double> fit4c_pimz_m, fit4c_pimz_p, fit4c_pimz_a, fit4c_pimz_pe, fit4c_pimz_px, fit4c_pimz_py, fit4c_pimz_pz;
};

#endif
