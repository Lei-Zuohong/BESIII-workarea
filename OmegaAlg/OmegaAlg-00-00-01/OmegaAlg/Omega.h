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
    // 容器信息
    double job_energy;
    int job_truth;
    int job_do_567;
    // topo信息
    NTuple::Item<int> runID;
    NTuple::Item<int> eventID;
    NTuple::Item<int> m_idxmc;
    NTuple::Array<int> m_pdgid;
    NTuple::Array<int> m_motheridx;
    // TREE - "truth"
    NTuple::Tuple *m_tuple_truth;
    NTuple::Item<double> truth_isr_m, truth_isr_p, truth_isr_a, truth_isr_pe, truth_isr_px, truth_isr_py, truth_isr_pz;
    NTuple::Item<double> truth_pip_m, truth_pip_p, truth_pip_a, truth_pip_pe, truth_pip_px, truth_pip_py, truth_pip_pz;
    NTuple::Item<double> truth_pim_m, truth_pim_p, truth_pim_a, truth_pim_pe, truth_pim_px, truth_pim_py, truth_pim_pz;
    NTuple::Item<double> truth_pi01_m, truth_pi01_p, truth_pi01_a, truth_pi01_pe, truth_pi01_px, truth_pi01_py, truth_pi01_pz;
    NTuple::Item<double> truth_pi02_m, truth_pi02_p, truth_pi02_a, truth_pi02_pe, truth_pi02_px, truth_pi02_py, truth_pi02_pz;
    NTuple::Item<double> truth_pi03_m, truth_pi03_p, truth_pi03_a, truth_pi03_pe, truth_pi03_px, truth_pi03_py, truth_pi03_pz;
    NTuple::Item<double> truth_omega_m, truth_omega_p, truth_omega_a, truth_omega_pe, truth_omega_px, truth_omega_py, truth_omega_pz;
    NTuple::Item<double> truth_omegapi02_m, truth_omegapi02_p, truth_omegapi02_a, truth_omegapi02_pe, truth_omegapi02_px, truth_omegapi02_py, truth_omegapi02_pz;
    NTuple::Item<double> truth_omegapi03_m, truth_omegapi03_p, truth_omegapi03_a, truth_omegapi03_pe, truth_omegapi03_px, truth_omegapi03_py, truth_omegapi03_pz;
    NTuple::Item<double> truth_pi02pi03_m, truth_pi02pi03_p, truth_pi02pi03_a, truth_pi02pi03_pe, truth_pi02pi03_px, truth_pi02pi03_py, truth_pi02pi03_pz;
    // TREE - "charge"
    NTuple::Tuple *m_tuple_charge;
    NTuple::Item<double> charge_ngood;
    NTuple::Item<double> charge_ncharge;
    // TREE - "vertex"
    NTuple::Tuple *m_tuple_vertex;
    NTuple::Item<double> vertex_chisq;
    // ************************************************************************************************************************
    // TREE - "fit4c"
    NTuple::Tuple *m_tuple_fit4c;
    NTuple::Item<double> fit4c_chisq;
    NTuple::Item<double> fit4c_pip_m, fit4c_pip_p, fit4c_pip_a, fit4c_pip_pe, fit4c_pip_px, fit4c_pip_py, fit4c_pip_pz;
    NTuple::Item<double> fit4c_pim_m, fit4c_pim_p, fit4c_pim_a, fit4c_pim_pe, fit4c_pim_px, fit4c_pim_py, fit4c_pim_pz;
    NTuple::Item<double> fit4c_pi01_m, fit4c_pi01_p, fit4c_pi01_a, fit4c_pi01_pe, fit4c_pi01_px, fit4c_pi01_py, fit4c_pi01_pz;
    NTuple::Item<double> fit4c_pi02_m, fit4c_pi02_p, fit4c_pi02_a, fit4c_pi02_pe, fit4c_pi02_px, fit4c_pi02_py, fit4c_pi02_pz;
    NTuple::Item<double> fit4c_pi03_m, fit4c_pi03_p, fit4c_pi03_a, fit4c_pi03_pe, fit4c_pi03_px, fit4c_pi03_py, fit4c_pi03_pz;
    NTuple::Item<double> fit4c_omega_m, fit4c_omega_p, fit4c_omega_a, fit4c_omega_pe, fit4c_omega_px, fit4c_omega_py, fit4c_omega_pz;
    NTuple::Item<double> fit4c_omegapi02_m, fit4c_omegapi02_p, fit4c_omegapi02_a, fit4c_omegapi02_pe, fit4c_omegapi02_px, fit4c_omegapi02_py, fit4c_omegapi02_pz;
    NTuple::Item<double> fit4c_omegapi03_m, fit4c_omegapi03_p, fit4c_omegapi03_a, fit4c_omegapi03_pe, fit4c_omegapi03_px, fit4c_omegapi03_py, fit4c_omegapi03_pz;
    NTuple::Item<double> fit4c_pi02pi03_m, fit4c_pi02pi03_p, fit4c_pi02pi03_a, fit4c_pi02pi03_pe, fit4c_pi02pi03_px, fit4c_pi02pi03_py, fit4c_pi02pi03_pz;

    NTuple::Item<double> fit4c_pif1_m, fit4c_pif2_m, fit4c_pif3_m, fit4c_pif4_m;
    // ************************************************************************************************************************
};
#endif
