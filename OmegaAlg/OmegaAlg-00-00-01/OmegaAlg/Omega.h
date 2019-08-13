#ifndef Physics_Analysis_Omega_H
#define Physics_Analysis_Omega_H 

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
//#include "VertexFit/ReadBeamParFromDb.h"

class Omega : public Algorithm {
public:
  Omega(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();  

private:

  double m_vr0cut;
  double m_vz0cut;
  double m_energyThreshold;
  double m_gammaPhiCut;
  double m_gammaThetaCut;
  double m_gammaAngleCut;
  int m_test4C;
  int m_test5C;
  int m_checkDedx;
  int m_checkTof;
  double m_energy;

  // define Ntuples here

  NTuple::Tuple*  m_tuple1;      // charged track vertex
  NTuple::Item<double>  m_vx0;
  NTuple::Item<double>  m_vy0;
  NTuple::Item<double>  m_vz0;
  NTuple::Item<double>  m_vr0;
  NTuple::Item<double>  m_rvxy0;
  NTuple::Item<double>  m_rvz0;
  NTuple::Item<double>  m_rvphi0;

  NTuple::Tuple*  m_tuple2;      // fake photon
  NTuple::Item<double>  m_dthe;
  NTuple::Item<double>  m_dphi;
  NTuple::Item<double>  m_dang;
  NTuple::Item<double>  m_eraw;

  NTuple::Tuple*  m_tuple3;     // Omega: raw mgg, etot
  NTuple::Item<double>  m_m2gg;
  NTuple::Item<double>  m_etot;

  NTuple::Tuple*  m_tuple4;     // Omega 4C
  NTuple::Item<double>  m_chi1;
  NTuple::Item<double>  m_mpi0;

  NTuple::Tuple*  m_tuple5;     // Omega 5C
  NTuple::Item<double>  m_chi2;
  //NTuple::Item<double>  m_mrh0;
  //NTuple::Item<double>  m_mrhp;
  //NTuple::Item<double>  m_mrhm;
  NTuple::Item<double>  m_mpi01;
  NTuple::Item<double>  m_mpi02;
  NTuple::Item<double>  m_mpi03;
  NTuple::Item<double>  m_momega;

  NTuple::Tuple*  m_tuple6;    // photons
  NTuple::Item<double>  m_fcos;
  NTuple::Item<double>  m_elow;

  NTuple::Tuple* m_tuple7;    // dE/dx
  NTuple::Item<double> m_ptrk;
  NTuple::Item<double> m_chie;
  NTuple::Item<double> m_chimu;
  NTuple::Item<double> m_chipi;
  NTuple::Item<double> m_chik;
  NTuple::Item<double> m_chip;
  NTuple::Item<double> m_probPH;
  NTuple::Item<double> m_normPH;
  NTuple::Item<double> m_ghit;
  NTuple::Item<double> m_thit;

  NTuple::Tuple* m_tuple8;   // endcap tof
  NTuple::Item<double> m_ptot_etof;
  NTuple::Item<double> m_cntr_etof;
  NTuple::Item<double> m_te_etof;
  NTuple::Item<double> m_tmu_etof;
  NTuple::Item<double> m_tpi_etof;
  NTuple::Item<double> m_tk_etof;
  NTuple::Item<double> m_tp_etof;
  NTuple::Item<double> m_ph_etof;
  NTuple::Item<double> m_rhit_etof;
  NTuple::Item<double> m_qual_etof;

  NTuple::Tuple* m_tuple9;  // barrel inner tof
  NTuple::Item<double> m_ptot_btof1;
  NTuple::Item<double> m_cntr_btof1;
  NTuple::Item<double> m_te_btof1;
  NTuple::Item<double> m_tmu_btof1;
  NTuple::Item<double> m_tpi_btof1;
  NTuple::Item<double> m_tk_btof1;
  NTuple::Item<double> m_tp_btof1;
  NTuple::Item<double> m_ph_btof1;
  NTuple::Item<double> m_zhit_btof1;
  NTuple::Item<double> m_qual_btof1;

  NTuple::Tuple* m_tuple10;  // barrel outer tof
  NTuple::Item<double> m_ptot_btof2;
  NTuple::Item<double> m_cntr_btof2;
  NTuple::Item<double> m_te_btof2;
  NTuple::Item<double> m_tmu_btof2;
  NTuple::Item<double> m_tpi_btof2;
  NTuple::Item<double> m_tk_btof2;
  NTuple::Item<double> m_tp_btof2;
  NTuple::Item<double> m_ph_btof2;
  NTuple::Item<double> m_zhit_btof2;
  NTuple::Item<double> m_qual_btof2;
  
  NTuple::Tuple* m_tuple11;  // Particle ID info.
  NTuple::Item<double> m_ptrk_pid;
  NTuple::Item<double> m_cost_pid;
  NTuple::Item<double> m_dedx_pid;
  NTuple::Item<double> m_tof1_pid;
  NTuple::Item<double> m_tof2_pid;
  NTuple::Item<double> m_prob_pid;

};

#endif 
