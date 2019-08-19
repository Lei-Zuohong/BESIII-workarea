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

  NTuple::Tuple*  m_tuple4;          // Omega 4C
  NTuple::Item<double>  m_chisq_4;
  NTuple::Item<double>  m_omega_4;
  NTuple::Item<double>  m_pi01_4;
  NTuple::Item<double>  m_pi02_4;
  NTuple::Item<double>  m_pi03_4;
  NTuple::Tuple*  m_tuple5;          // Omega 5C
  NTuple::Item<double>  m_chi2;
  NTuple::Item<double>  m_mpi01;
  NTuple::Item<double>  m_mpi02;
  NTuple::Item<double>  m_mpi03;
  NTuple::Item<double>  m_momega;

  NTuple::Tuple*  m_tuple6;          // photons
  NTuple::Item<double>  m_fcos;
  NTuple::Item<double>  m_elow;


};

#endif 