//*********************************************************************************************************
//***                                                代码引用                                            ***
//*********************************************************************************************************
// include McTruth
#include "McTruth/McParticle.h"
// include GaudiKernel
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
// include other
#include "EventModel/EventModel.h"
#include "EventModel/EventHeader.h"
#include "EventModel/Event.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"
#include "TMath.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"
using CLHEP::Hep2Vector;
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include <vector>
#include "PmAlg/Pm.h"
//*********************************************************************************************************
//***                                                变量设置                                            ***
//*********************************************************************************************************
const double mpi = 0.13957;
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
std::vector<int> SeriesRun;
std::vector<int> SeriesNum;
int Ncut0, Ncut1, Ncut2, Ncut3, Ncut4, Ncut5, number, checkexit, checki, firstrun;
//*********************************************************************************************************
//***                                                声明容器                                            ***
//*********************************************************************************************************
Pm::Pm(const std::string &name, ISvcLocator *pSvcLocator) : Algorithm(name, pSvcLocator)
{
	declareProperty("Energy", m_energy = 2.125);
}
//*********************************************************************************************************
//***                                               initialize                                          ***
//*********************************************************************************************************
StatusCode Pm::initialize()
{
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in initialize()" << endmsg;
	StatusCode status;
	// tree for "truth"
	if (1 == 1)
	{
		NTuplePtr nt1(ntupleSvc(), "FILE1/truth");
		if (nt1)
		{
			m_tuple1 = nt1;
		}
		else
		{
			m_tuple1 = ntupleSvc()->book("FILE1/truth", CLID_ColumnWiseTuple, "ks N-Tuple example");
			if (m_tuple1)
			{
				// isr
				status = m_tuple1->addItem("misr", t1_misr);
				status = m_tuple1->addItem("aisr", t1_aisr);
				status = m_tuple1->addItem("pisr", t1_pisr);
				status = m_tuple1->addItem("eisr", t1_eisr);
				status = m_tuple1->addItem("pxisr", t1_pxisr);
				status = m_tuple1->addItem("pyisr", t1_pyisr);
				status = m_tuple1->addItem("pzisr", t1_pzisr);
				// pip2 pim2
				status = m_tuple1->addItem("mpip2", t1_mpip2);
				status = m_tuple1->addItem("apip2", t1_apip2);
				status = m_tuple1->addItem("ppip2", t1_ppip2);
				status = m_tuple1->addItem("epip2", t1_epip2);
				status = m_tuple1->addItem("pxpip2", t1_pxpip2);
				status = m_tuple1->addItem("pypip2", t1_pypip2);
				status = m_tuple1->addItem("pzpip2", t1_pzpip2);
				status = m_tuple1->addItem("mpim2", t1_mpim2);
				status = m_tuple1->addItem("apim2", t1_apim2);
				status = m_tuple1->addItem("ppim2", t1_ppim2);
				status = m_tuple1->addItem("epim2", t1_epim2);
				status = m_tuple1->addItem("pxpim2", t1_pxpim2);
				status = m_tuple1->addItem("pypim2", t1_pypim2);
				status = m_tuple1->addItem("pzpim2", t1_pzpim2);
				// pi0 pip1 pim1 w wpip2 wpim2 pip2pip3
				status = m_tuple1->addItem("mpi0", t1_mpi0);
				status = m_tuple1->addItem("api0", t1_api0);
				status = m_tuple1->addItem("ppi0", t1_ppi0);
				status = m_tuple1->addItem("epi0", t1_epi0);
				status = m_tuple1->addItem("pxpi0", t1_pxpi0);
				status = m_tuple1->addItem("pypi0", t1_pypi0);
				status = m_tuple1->addItem("pzpi0", t1_pzpi0);
				status = m_tuple1->addItem("mpip1", t1_mpip1);
				status = m_tuple1->addItem("apip1", t1_apip1);
				status = m_tuple1->addItem("ppip1", t1_ppip1);
				status = m_tuple1->addItem("epip1", t1_epip1);
				status = m_tuple1->addItem("pxpip1", t1_pxpip1);
				status = m_tuple1->addItem("pypip1", t1_pypip1);
				status = m_tuple1->addItem("pzpip1", t1_pzpip1);
				status = m_tuple1->addItem("mpim1", t1_mpim1);
				status = m_tuple1->addItem("apim1", t1_apim1);
				status = m_tuple1->addItem("ppim1", t1_ppim1);
				status = m_tuple1->addItem("epim1", t1_epim1);
				status = m_tuple1->addItem("pxpim1", t1_pxpim1);
				status = m_tuple1->addItem("pypim1", t1_pypim1);
				status = m_tuple1->addItem("pzpim1", t1_pzpim1);
				status = m_tuple1->addItem("momega", t1_momega);
				status = m_tuple1->addItem("aomega", t1_aomega);
				status = m_tuple1->addItem("pomega", t1_pomega);
				status = m_tuple1->addItem("eomega", t1_eomega);
				status = m_tuple1->addItem("pxomega", t1_pxomega);
				status = m_tuple1->addItem("pyomega", t1_pyomega);
				status = m_tuple1->addItem("pzomega", t1_pzomega);
				status = m_tuple1->addItem("momegapip2", t1_momegapip2);
				status = m_tuple1->addItem("aomegapip2", t1_aomegapip2);
				status = m_tuple1->addItem("pomegapip2", t1_pomegapip2);
				status = m_tuple1->addItem("momegapim2", t1_momegapim2);
				status = m_tuple1->addItem("aomegapim2", t1_aomegapim2);
				status = m_tuple1->addItem("pomegapim2", t1_pomegapim2);
				status = m_tuple1->addItem("mpip2pim2", t1_mpip2pim2);
				status = m_tuple1->addItem("apip2pim2", t1_apip2pim2);
				status = m_tuple1->addItem("ppip2pim2", t1_ppip2pim2);
			}
			else
			{
				log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
	// tree for "fit4c"
	if (1 == 1)
	{
		NTuplePtr nt4(ntupleSvc(), "FILE1/fit4c");
		if (nt4)
		{
			m_tuple4 = nt4;
		}
		else
		{
			m_tuple4 = ntupleSvc()->book("FILE1/fit4c", CLID_ColumnWiseTuple, "ks N-Tuple example");
			if (m_tuple4)
			{
				status = m_tuple4->addItem("chisq", t4_chisq_4c);
				// topo 信息
				status = m_tuple4->addItem("runID", runID);
				status = m_tuple4->addItem("eventID", eventID);
				status = m_tuple4->addItem("indexmc", m_idxmc, 0, 100);
				status = m_tuple4->addIndexedItem("pdgid", m_idxmc, m_pdgid);
				status = m_tuple4->addIndexedItem("motheridx", m_idxmc, m_motheridx);
				// pip2 pim2
				status = m_tuple4->addItem("mpip2", t4_mpip2);
				status = m_tuple4->addItem("apip2", t4_apip2);
				status = m_tuple4->addItem("ppip2", t4_ppip2);
				status = m_tuple4->addItem("epip2", t4_epip2);
				status = m_tuple4->addItem("pxpip2", t4_pxpip2);
				status = m_tuple4->addItem("pypip2", t4_pypip2);
				status = m_tuple4->addItem("pzpip2", t4_pzpip2);
				status = m_tuple4->addItem("mpim2", t4_mpim2);
				status = m_tuple4->addItem("apim2", t4_apim2);
				status = m_tuple4->addItem("ppim2", t4_ppim2);
				status = m_tuple4->addItem("epim2", t4_epim2);
				status = m_tuple4->addItem("pxpim2", t4_pxpim2);
				status = m_tuple4->addItem("pypim2", t4_pypim2);
				status = m_tuple4->addItem("pzpim2", t4_pzpim2);
				// pi01 pi02 pi03 w wpi02 wpi03 pi02pi03
				status = m_tuple4->addItem("mpi0", t4_mpi0);
				status = m_tuple4->addItem("api0", t4_api0);
				status = m_tuple4->addItem("ppi0", t4_ppi0);
				status = m_tuple4->addItem("epi0", t4_epi0);
				status = m_tuple4->addItem("pxpi0", t4_pxpi0);
				status = m_tuple4->addItem("pypi0", t4_pypi0);
				status = m_tuple4->addItem("pzpi0", t4_pzpi0);
				status = m_tuple4->addItem("mpip1", t4_mpip1);
				status = m_tuple4->addItem("apip1", t4_apip1);
				status = m_tuple4->addItem("ppip1", t4_ppip1);
				status = m_tuple4->addItem("epip1", t4_epip1);
				status = m_tuple4->addItem("pxpip1", t4_pxpip1);
				status = m_tuple4->addItem("pypip1", t4_pypip1);
				status = m_tuple4->addItem("pzpip1", t4_pzpip1);
				status = m_tuple4->addItem("mpim1", t4_mpim1);
				status = m_tuple4->addItem("apim1", t4_apim1);
				status = m_tuple4->addItem("ppim1", t4_ppim1);
				status = m_tuple4->addItem("epim1", t4_epim1);
				status = m_tuple4->addItem("pxpim1", t4_pxpim1);
				status = m_tuple4->addItem("pypim1", t4_pypim1);
				status = m_tuple4->addItem("pzpim1", t4_pzpim1);
				status = m_tuple4->addItem("momega", t4_momega);
				status = m_tuple4->addItem("aomega", t4_aomega);
				status = m_tuple4->addItem("pomega", t4_pomega);
				status = m_tuple4->addItem("eomega", t4_eomega);
				status = m_tuple4->addItem("pxomega", t4_pxomega);
				status = m_tuple4->addItem("pyomega", t4_pyomega);
				status = m_tuple4->addItem("pzomega", t4_pzomega);
				status = m_tuple4->addItem("momegapip2", t4_momegapip2);
				status = m_tuple4->addItem("aomegapip2", t4_aomegapip2);
				status = m_tuple4->addItem("pomegapip2", t4_pomegapip2);
				status = m_tuple4->addItem("momegapim2", t4_momegapim2);
				status = m_tuple4->addItem("aomegapim2", t4_aomegapim2);
				status = m_tuple4->addItem("pomegapim2", t4_pomegapim2);
				status = m_tuple4->addItem("mpip2pim2", t4_mpip2pim2);
				status = m_tuple4->addItem("apip2pim2", t4_apip2pim2);
				status = m_tuple4->addItem("ppip2pim2", t4_ppip2pim2);
			}
			else
			{
				log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple4) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
	log << MSG::INFO << "successfully return from initialize()" << endmsg;
	return StatusCode::SUCCESS;
}
//*********************************************************************************************************
//***                                                execute                                            ***
//*********************************************************************************************************
StatusCode Pm::execute()																//
{																						//
	MsgStream log(msgSvc(), name());													//
	log << MSG::INFO << "in execute()" << endreq;										//
	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");		//
	int runNo = eventHeader->runNumber();												// 读取runNo：runnumber
	int event = eventHeader->eventNumber();												// 读取event：eventnumber
	runID = runNo;																		// 变量：topo
	eventID = event;																	//
	log << MSG::DEBUG << "run, evtnum = "												//
		<< runNo << " , "																//
		<< event << endreq;																//
	Ncut0++;																			//
	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent); //
	log << MSG::DEBUG << "ncharg, nneu, tottks = "										//
		<< evtRecEvent->totalCharged() << " , "											//
		<< evtRecEvent->totalNeutral() << " , "											//
		<< evtRecEvent->totalTracks() << endreq;										//
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
	//*********************************************************************************
	// Selection 0: Topology
	//*********************************************************************************
	if (eventHeader->runNumber() < 0)															 //
	{																							 //
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol"); //
		int m_numParticle = 1;																	 //
		if (!mcParticleCol)																		 //
		{																						 //
			cout << "Could not retrieve McParticelCol" << endl;									 //
			return StatusCode::FAILURE;															 //
		}																						 //
		else																					 //
		{																						 //
			int nprmary = 0;																	 //
			Event::McParticleCol::iterator iter_mc1 = mcParticleCol->begin();					 //
			for (; iter_mc1 != mcParticleCol->end(); iter_mc1++)								 //
			{																					 //
				if (!(*iter_mc1)->decayFromGenerator())											 //
					continue;																	 //
				if ((*iter_mc1)->primaryParticle())												 //
				{																				 //
					nprmary++;																	 //
				}																				 //
			}																					 //
			Event::McParticleCol::iterator iter_mc2 = mcParticleCol->begin();					 //
			if (nprmary == 1)																	 //
			{																					 //
				m_numParticle = 0;																 //
				for (; iter_mc2 != mcParticleCol->end(); iter_mc2++)							 //
				{																				 //
					if (!(*iter_mc2)->decayFromGenerator())										 //
						continue;																 //
					if ((*iter_mc2)->primaryParticle())											 //
					{																			 //
						m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();				 // 变量：topo
						m_motheridx[m_numParticle] = 0;											 // 变量：topo
					}																			 //
					else																		 //
					{																			 //
						m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();				 // 变量：topo
						m_motheridx[m_numParticle] = ((*iter_mc2)->mother()).trackIndex();		 // 变量：topo
					}																			 //
					m_numParticle += 1;															 //
				}																				 //
				m_idxmc = m_numParticle;														 // 变量：topo
			}																					 //
			if (nprmary > 1)																	 //
			{																					 //
				m_numParticle = 1;																 //
				for (; iter_mc2 != mcParticleCol->end(); iter_mc2++)							 //
				{																				 //
					if (!(*iter_mc2)->decayFromGenerator())										 //
						continue;																 //
					if ((*iter_mc2)->primaryParticle())											 //
					{																			 //
						m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();				 // 变量：topo
						m_motheridx[m_numParticle] = 0;											 // 变量：topo
					}																			 //
					else																		 //
					{																			 //
																								 //
						m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();				 // 变量：topo
						m_motheridx[m_numParticle] = ((*iter_mc2)->mother()).trackIndex() + 1;   // 变量：topo
					}																			 //
					m_numParticle += 1;															 //
					m_pdgid[0] = 11111;															 // 变量：topo
					m_motheridx[0] = 0;															 // 变量：topo
				}																				 //
				m_idxmc = m_numParticle;														 // 变量：topo
			}																					 //
		}																						 //
	}
	//*********************************************************************************
	// Selection 0: Truth
	//*********************************************************************************
	int truth_check = 0;
	if (eventHeader->runNumber() < 0)
	{
		// 设定变量
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
		HepLorentzVector isr_track;
		HepLorentzVector pip1_track;
		HepLorentzVector pip2_track;
		HepLorentzVector pim1_track;
		HepLorentzVector pim2_track;
		HepLorentzVector pi0_track;
		int nisr = 0;
		int npip1 = 0;
		int npip2 = 0;
		int npim1 = 0;
		int npim2 = 0;
		int npi0 = 0;
		HepLorentzVector omega_track;
		int nomega = 0;
		HepLorentzVector omegapip2_track;
		HepLorentzVector omegapim2_track;
		HepLorentzVector pip2pim2_track;
		Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
		// 统计 omega, pi+, pi-, pi0
		for (; iter_mc != mcParticleCol->end(); iter_mc++)
		{
			if (!(*iter_mc)->decayFromGenerator())
				continue;
			HepLorentzVector mctrue_track = (*iter_mc)->initialFourMomentum();
			int mctrue_index = (*iter_mc)->trackIndex();
			// 统计omega
			if ((*iter_mc)->particleProperty() == 223)
			{
				omega_track = mctrue_track;
				nomega += 1;
			}
			// 统计pi+
			if ((*iter_mc)->particleProperty() == 211)
			{
				if (((*iter_mc)->mother()).particleProperty() == 223)
				{
					pip1_track = mctrue_track;
					npip1 += 1;
				}
				else
				{
					pip2_track = mctrue_track;
					npip2 += 1;
				}
			}
			// 统计pi-
			if ((*iter_mc)->particleProperty() == -211)
			{
				if (((*iter_mc)->mother()).particleProperty() == 223)
				{
					pim1_track = mctrue_track;
					npim1 += 1;
				}
				else
				{
					pim2_track = mctrue_track;
					npim2 += 1;
				}
			}
			// 统计pi0
			if ((*iter_mc)->particleProperty() == 111)
			{
				if (((*iter_mc)->mother()).particleProperty() == 223)
				{
					pi0_track = mctrue_track;
					npi0 += 1;
				}
			}
			// 统计isr
			if ((*iter_mc)->particleProperty() == 22)
			{
				if (((*iter_mc)->mother()).particleProperty() != 111)
				{
					isr_track = mctrue_track;
					nisr += 1;
				}
			}
		}
		if (npip1 == 1 && npip2 == 1 && npim1 == 1 && npim2 == 1 && npi0 == 1 && nomega == 1)
		{
			truth_check = 1;
		}
		//填入信息
		if (truth_check == 1)
		{
			omega_track = pip1_track + pim1_track + pi0_track;
			omegapip2_track = omega_track + pip2_track;
			omegapim2_track = omega_track + pim2_track;
			pip2pim2_track = pip2_track + pim2_track;
			// 输出gamma信息
			if (1 == 1)
			{
				//isr
				t1_misr = isr_track.m();
				t1_aisr = cos(isr_track.theta());
				t1_pisr = isr_track.rho();
				t1_eisr = isr_track.e();
				t1_pxisr = isr_track.px();
				t1_pyisr = isr_track.py();
				t1_pzisr = isr_track.pz();
				// pip2 pim2
				t1_mpip2 = pip2_track.m();
				t1_apip2 = cos(pip2_track.theta());
				t1_ppip2 = pip2_track.rho();
				t1_epip2 = pip2_track.e();
				t1_pxpip2 = pip2_track.px();
				t1_pypip2 = pip2_track.py();
				t1_pzpip2 = pip2_track.pz();
				t1_mpim2 = pim2_track.m();
				t1_apim2 = cos(pim2_track.theta());
				t1_ppim2 = pim2_track.rho();
				t1_epim2 = pim2_track.e();
				t1_pxpim2 = pim2_track.px();
				t1_pypim2 = pim2_track.py();
				t1_pzpim2 = pim2_track.pz();
				// pi0 pip1 pim1 w wpip2 wpim2 pip2pip3
				t1_mpi0 = pi0_track.m();
				t1_api0 = cos(pi0_track.theta());
				t1_ppi0 = pi0_track.rho();
				t1_epi0 = pi0_track.e();
				t1_pxpi0 = pi0_track.px();
				t1_pypi0 = pi0_track.py();
				t1_pzpi0 = pi0_track.pz();
				t1_mpip1 = pip1_track.m();
				t1_apip1 = cos(pip1_track.theta());
				t1_ppip1 = pip1_track.rho();
				t1_epip1 = pip1_track.e();
				t1_pxpip1 = pip1_track.px();
				t1_pypip1 = pip1_track.py();
				t1_pzpip1 = pip1_track.pz();
				t1_mpim1 = pim1_track.m();
				t1_apim1 = cos(pim1_track.theta());
				t1_ppim1 = pim1_track.rho();
				t1_epim1 = pim1_track.e();
				t1_pxpim1 = pim1_track.px();
				t1_pypim1 = pim1_track.py();
				t1_pzpim1 = pim1_track.pz();
				t1_momega = omega_track.m();
				t1_aomega = cos(omega_track.theta());
				t1_pomega = omega_track.rho();
				t1_eomega = omega_track.e();
				t1_pxomega = omega_track.px();
				t1_pyomega = omega_track.py();
				t1_pzomega = omega_track.pz();
				t1_momegapip2 = omegapip2_track.m();
				t1_aomegapip2 = cos(omegapip2_track.theta());
				t1_pomegapip2 = omegapip2_track.rho();
				t1_momegapim2 = omegapim2_track.m();
				t1_aomegapim2 = cos(omegapim2_track.theta());
				t1_pomegapim2 = omegapim2_track.rho();
				t1_mpip2pim2 = pip2pim2_track.m();
				t1_apip2pim2 = cos(pip2pim2_track.theta());
				t1_ppip2pim2 = pip2pim2_track.rho();
			}
			m_tuple1->write();
			Ncut1 += 1;
		}
	}
	//*********************************************************************************
	// Selection 0: Check runnumber
	//*********************************************************************************
	if (eventHeader->runNumber() > 0)
	{
		int m_status = 0;
		if (runNo == 34326 || runNo == 34334 || runNo == 34478 || runNo == 34818 || runNo == 34982 ||
			runNo == 35101 || runNo == 40459 || runNo == 40460 || runNo == 40461 || runNo == 40462 ||
			runNo == 41408 || runNo == 41416 || runNo == 41902)
		{
			m_status = -1;
		}
		if (runNo == 34018 || runNo == 34020 || runNo == 34021 ||
			runNo == 34023 || runNo == 34027 || runNo == 34040 || runNo == 34042 || runNo == 34046 ||
			runNo == 34048 || runNo == 34058 || runNo == 34061 || runNo == 34065 || runNo == 34099 ||
			runNo == 34100 || runNo == 34110 || runNo == 34117 || runNo == 34123 || runNo == 34124 ||
			runNo == 34145 || runNo == 34152 || runNo == 34153 || runNo == 34168 || runNo == 34169 ||
			runNo == 34170 || runNo == 34171 || runNo == 34172 || runNo == 34173 || runNo == 34174 ||
			runNo == 34177 || runNo == 34178 || runNo == 34211 || runNo == 34212 || runNo == 34217 ||
			runNo == 34227 || runNo == 34228 || runNo == 34237 || runNo == 34238 || runNo == 34261 ||
			runNo == 34262 || runNo == 34282 || runNo == 34283 || runNo == 34301 || runNo == 34308 ||
			runNo == 34309 || runNo == 34336 || runNo == 34337 || runNo == 34338 || runNo == 34352 ||
			runNo == 34353 || runNo == 34377 || runNo == 34378 || runNo == 34383 || runNo == 34392 ||
			runNo == 34401 || runNo == 34404 || runNo == 34406 || runNo == 34429 || runNo == 34432 ||
			runNo == 34434 || runNo == 34440 || runNo == 34441 || runNo == 34464 || runNo == 34465 ||
			runNo == 34467 || runNo == 34490 || runNo == 34492 || runNo == 34516 || runNo == 34517 ||
			runNo == 34518 || runNo == 34519 || runNo == 34520 || runNo == 34521 || runNo == 34522 ||
			runNo == 34523 || runNo == 34525 || runNo == 34527 || runNo == 34528 || runNo == 34529 ||
			runNo == 34537 || runNo == 34544 || runNo == 34546 || runNo == 34547 || runNo == 34548 ||
			runNo == 34549 || runNo == 34576 || runNo == 34577 || runNo == 34578 || runNo == 34600 ||
			runNo == 34602 || runNo == 34625 || runNo == 34626 || runNo == 34628 || runNo == 34650 ||
			runNo == 34653 || runNo == 34676 || runNo == 34678 || runNo == 34700 || runNo == 34711 ||
			runNo == 34720 || runNo == 34722 || runNo == 34744 || runNo == 34745 || runNo == 34746 ||
			runNo == 34748 || runNo == 34770 || runNo == 34771 || runNo == 34772 || runNo == 34795 ||
			runNo == 34798 || runNo == 34815 || runNo == 34817 || runNo == 34848 || runNo == 34849 ||
			runNo == 34853 || runNo == 34854 || runNo == 34870 || runNo == 34872 || runNo == 34875 ||
			runNo == 34880 || runNo == 34882 || runNo == 34900 || runNo == 34904 || runNo == 34907 ||
			runNo == 34916 || runNo == 34917 || runNo == 34936 || runNo == 34951 || runNo == 34954 ||
			runNo == 34975 || runNo == 34977 || runNo == 34979 || runNo == 34985 || runNo == 34986 ||
			runNo == 34988 || runNo == 34989 || runNo == 34990 || runNo == 35006 || runNo == 35014 ||
			runNo == 35021 || runNo == 35023 || runNo == 35043 || runNo == 35044 || runNo == 35045 ||
			runNo == 35046 || runNo == 35052 || runNo == 35066 || runNo == 35071 || runNo == 35074 ||
			runNo == 35080 || runNo == 35091 || runNo == 35092 || runNo == 35096 || runNo == 35103 ||
			runNo == 35117 || runNo == 35118)
		{
			m_status = -1;
		}
		if (m_status == -1)
		{
			return StatusCode::SUCCESS;
		}
	}
	//*********************************************************************************
	// Selection 0: Statistic Runnumber
	//*********************************************************************************
	if (firstrun == 0)						   //
	{										   //
		number = 0;							   //
		checkexit = 0;						   //
		checki = 0;							   //
		SeriesRun.clear();					   //
		SeriesNum.clear();					   //
	}										   //
	firstrun = 1;							   //
	checkexit = 0;							   //
	for (int i = 0; i < SeriesRun.size(); i++) //
	{										   //
		if (runNo == SeriesRun[i])			   //
		{									   //
			checkexit = 1;					   //
			checki = i;						   //
		}									   //
	}										   //
	if (checkexit == 1)						   //
	{										   //
		SeriesNum[checki] += 1;				   //
	}										   //
	else									   //
	{										   //
		SeriesRun.push_back(runNo);			   //
		SeriesNum.push_back(1);				   //
	}
	//*********************************************************************************
	// Selection 1: Good Charged Track Selection
	//*********************************************************************************
	Vint iGood, ipip, ipim;															   //
	iGood.clear();																	   // 变量：iGood[]（参数为good-track序号，内容为track编号）
	ipip.clear();																	   // 变量：ipip[]（参数为good-track+序号，内容为track编号）
	ipim.clear();																	   // 变量：ipim[]（参数为good-track-序号，内容为track编号）
	Vp4 ppip, ppim;																	   //
	ppip.clear();																	   // 变量：ppip[]（good-track+的四动量）
	ppim.clear();																	   // 变量：ppim[]（good-track-的四动量）
	int nCharge;																	   //
	nCharge = 0;																	   //
	Hep3Vector xorigin(0, 0, 0);													   //
	IVertexDbSvc *vtxsvc;															   //
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);							   //
	if (vtxsvc->isVertexValid())													   //
	{																				   //
		double *dbv = vtxsvc->PrimaryVertex();										   //
		double *vv = vtxsvc->SigmaPrimaryVertex();									   //
		xorigin.setX(dbv[0]);														   //
		xorigin.setY(dbv[1]);														   //
		xorigin.setZ(dbv[2]);														   //
	}																				   //
	for (int i = 0; i < evtRecEvent->totalCharged(); i++)							   //
	{																				   //
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;						   //
		if (!(*itTrk)->isMdcTrackValid())											   //
			continue;																   //
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();									   //
		double pch = mdcTrk->p();													   //
		double x0 = mdcTrk->x();													   //
		double y0 = mdcTrk->y();													   //
		double z0 = mdcTrk->z();													   //
		double phi0 = mdcTrk->helix(1);												   //表示螺旋线参数
																					   //0: d0 -> 螺旋线到对撞顶点的最小距离
																					   //1: phi0 -> 最小距离的xy平面相角
																					   //2: kappa
																					   //3: d
																					   //4: tan(lamda)
		double xv = xorigin.x();													   //
		double yv = xorigin.y();													   //
		double Rxy = (x0 - xv) * cos(phi0) + (y0 - yv) * sin(phi0);					   //
		HepVector a = mdcTrk->helix();												   //
		HepSymMatrix Ea = mdcTrk->err();											   //
		HepPoint3D point0(0., 0., 0.);												   //
		HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);							   //
		VFHelix helixip(point0, a, Ea);												   //
		helixip.pivot(IP);															   //
		HepVector vecipa = helixip.a();												   //
		double Rvxy0 = fabs(vecipa[0]);												   //
		double Rvz0 = vecipa[3];													   //
		double Rvphi0 = vecipa[1];													   //
		if (fabs(cos(mdcTrk->theta())) > 0.93)										   // 选择：cos(theta)
			continue;																   //
		if (fabs(Rvz0) >= 10)														   // 选择：Rvz0
			continue;																   //
		if (fabs(Rvxy0) >= 1)														   // 选择：Rvxy0
			continue;																   //
		iGood.push_back(i);															   //
		nCharge += mdcTrk->charge();												   //
	}																				   //
	int nGood = iGood.size();														   // 变量：nGood（带电track数目）
	log << MSG::DEBUG << "ngood, totcharge = " << nGood << " , " << nCharge << endreq; // 变量：nCharge（带电track总电量）
	if ((nGood != 4) || (nCharge != 0))												   // 选择：nGood
	{																				   // 选择：nCharge
		return StatusCode::SUCCESS;													   //
	}																				   //
	//*********************************************************************************
	// Selection 2: Good Photon Selection
	//*********************************************************************************
	Vint iGam;																							  //
	iGam.clear();																						  //
	for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks(); i++)						  //
	{																									  //
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;											  //
		if (!(*itTrk)->isEmcShowerValid())																  //
			continue;																					  //
		RecEmcShower *emcTrk = (*itTrk)->emcShower();													  //
		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());										  //
		double dthe = 200.;																				  //
		double dphi = 200.;																				  //
		double dang = 200.;																				  //
		for (int j = 0; j < evtRecEvent->totalCharged(); j++)											  //
		{																								  //
			EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;										  //
			if (!(*jtTrk)->isExtTrackValid())															  //
				continue;																				  //
			RecExtTrack *extTrk = (*jtTrk)->extTrack();													  //
			if (extTrk->emcVolumeNumber() == -1)														  //
				continue;																				  //
			Hep3Vector extpos = extTrk->emcPosition();													  //
			double angd = extpos.angle(emcpos);															  //
			double thed = extpos.theta() - emcpos.theta();												  //
			double phid = extpos.deltaPhi(emcpos);														  //
			thed = fmod(thed + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) - CLHEP::pi;				  //
			phid = fmod(phid + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) - CLHEP::pi;				  //
			if (angd < dang)																			  //
			{																							  //
				dang = angd;																			  //
				dthe = thed;																			  //
				dphi = phid;																			  //
			}																							  //
		}																								  //
		if (dang >= 200)																				  //
			continue;																					  //
		double eraw = emcTrk->energy();																	  //
		dthe = dthe * 180 / (CLHEP::pi);																  //
		dphi = dphi * 180 / (CLHEP::pi);																  //
		dang = dang * 180 / (CLHEP::pi);																  //
		double ctht = cos(emcTrk->theta());																  //
		if ((emcTrk->module() == 1) && (fabs(ctht) > 0.8))												  // 选择：E-endcap
			continue;																					  //
		if ((emcTrk->module() == 0 || emcTrk->module() == 2) && (fabs(ctht) > 0.92 || fabs(ctht) < 0.86)) // 选择：E-endcap
			continue;																					  //
		if (emcTrk->time() < 0 || emcTrk->time() > 14)													  // 选择：TDC
			continue;																					  //
		if ((emcTrk->module() == 1) && eraw < 0.025)													  // 选择：E-barrel
			continue;																					  //
		if ((emcTrk->module() == 0 || emcTrk->module() == 2) && eraw < 0.050)							  // 选择：E-barrel
			continue;																					  //
		if (fabs(dang) < 10)																			  // 选择：
			continue;																					  //
		iGam.push_back(i);																				  // 变量：iGam[]（参数为good-track序号，内容为track编号）
	}																									  //
	int nGam = iGam.size();																				  // 变量：nGam（中性track数量）
	log << MSG::DEBUG << "num Good Photon " << nGam << " , " << evtRecEvent->totalNeutral() << endreq;	//
	if (nGam < 2 || nGam > 55)																			  // 选择：nGam
	{																									  //
		return StatusCode::SUCCESS;																		  //
	}
	//*********************************************************************************
	// Calculation 2: 4-momentum to each photon
	//*********************************************************************************
	Vp4 pGam;														 // 变量：pGam[]（参数为good-track序号，内容为动量）
	pGam.clear();													 // 变量：iGam[]（参数为good-track序号，内容为track编号）
	for (int i = 0; i < nGam; i++)									 // 变量：nGam（中性track数量）
	{																 //
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGam[i]; //
		RecEmcShower *emcTrk = (*itTrk)->emcShower();				 //
		double eraw = emcTrk->energy();								 //
		double phi = emcTrk->phi();									 //
		double the = emcTrk->theta();								 //
		HepLorentzVector ptrk;										 //
		ptrk.setPx(eraw * sin(the) * cos(phi));						 //
		ptrk.setPy(eraw * sin(the) * sin(phi));						 //
		ptrk.setPz(eraw * cos(the));								 //
		ptrk.setE(eraw);											 //
		pGam.push_back(ptrk);										 //
	}																 //
	//*********************************************************************************
	// Calculation 2: 4-momentum to each charged track
	//*********************************************************************************
	ParticleID *pid = ParticleID::instance();											   //
	for (int i = 0; i < nGood; i++)														   //
	{																					   //
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];					   //
		pid->init();																	   // 对于Likelihood方法
		pid->setMethod(pid->methodProbability());										   // pid->setMethod(pid->methodLikelihood());
		pid->setChiMinCut(4);															   //
		pid->setRecTrack(*itTrk);														   //
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE()); // 选择pid的粒子
		pid->identify(pid->onlyPion() | pid->onlyKaon());								   // pid->identify(pid->onlyPion());
		pid->calculate();																   // pid->identify(pid->onlyKaon());
		if (!(pid->IsPidInfoValid()))													   //
			continue;																	   // 对于Likelihood方法(0=electron 1=muon 2=pion 3=kaon 4=proton)
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();										   // if(pid->pdf(2) < pid->pdf(3))
		if (pid->probPion() < pid->probKaon())											   // continue;
			continue;																	   //
		RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();							   // 对于ParticleID, 用RecMdcKalTrack代替RecMdcTrack
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);								   // PID可以设定为electron, muon, pion, kaon and proton;The default setting is pion
		if (mdcKalTrk->charge() > 0)													   //
		{																				   //
			ipip.push_back(iGood[i]);													   //
			HepLorentzVector ptrk;														   //
			ptrk.setPx(mdcKalTrk->px());												   //
			ptrk.setPy(mdcKalTrk->py());												   //
			ptrk.setPz(mdcKalTrk->pz());												   //
			double p3 = ptrk.mag();														   //
			ptrk.setE(sqrt(p3 * p3 + mpi * mpi));										   //
			ppip.push_back(ptrk);														   // 变量：ppip[]（值为pi+动量）
		}																				   // 变量：ppim[]（值为pi-动量）
		else																			   //
		{																				   // 变量：ipip[]（值为pi+编号）
			ipim.push_back(iGood[i]);													   // 变量：ipim[]（值为pi-编号）
			HepLorentzVector ptrk;														   //
			ptrk.setPx(mdcKalTrk->px());												   //
			ptrk.setPy(mdcKalTrk->py());												   //
			ptrk.setPz(mdcKalTrk->pz());												   //
			double p3 = ptrk.mag();														   //
			ptrk.setE(sqrt(p3 * p3 + mpi * mpi));										   //
			ppim.push_back(ptrk);														   //
		}																				   //
	}																					   //
	int npip = ipip.size();																   // 变量：npip（pi+数目）
	int npim = ipim.size();																   // 变量：npim（pi-数目）
	if (npip * npim != 4)																   //
		return SUCCESS;																	   //
	Ncut3++;																			   //
	//*********************************************************************************
	// Selection 3: Vertex fit Selection, check ppi0, pTot
	//*********************************************************************************
	RecMdcKalTrack *pipTrk1 = (*(evtRecTrkCol->begin() + ipip[0]))->mdcKalTrack(); // 读取track
	RecMdcKalTrack *pipTrk2 = (*(evtRecTrkCol->begin() + ipip[1]))->mdcKalTrack(); //
	RecMdcKalTrack *pimTrk1 = (*(evtRecTrkCol->begin() + ipim[0]))->mdcKalTrack(); //
	RecMdcKalTrack *pimTrk2 = (*(evtRecTrkCol->begin() + ipim[1]))->mdcKalTrack(); //
	WTrackParameter wvpipTrk1;													   //
	WTrackParameter wvpipTrk2;													   //
	WTrackParameter wvpimTrk1;													   //
	WTrackParameter wvpimTrk2;													   //
	wvpipTrk1 = WTrackParameter(mpi, pipTrk1->getZHelix(), pipTrk1->getZError());  //
	wvpipTrk2 = WTrackParameter(mpi, pipTrk2->getZHelix(), pipTrk2->getZError());  //
	wvpimTrk1 = WTrackParameter(mpi, pimTrk1->getZHelix(), pimTrk1->getZError());  //
	wvpimTrk2 = WTrackParameter(mpi, pimTrk2->getZHelix(), pimTrk2->getZError());  //
	HepPoint3D vx(0., 0., 0.);													   // 设定顶点
	HepSymMatrix Evx(3, 0);														   //
	double bx = 1E+6;															   //
	double by = 1E+6;															   //
	double bz = 1E+6;															   //
	Evx[0][0] = bx * bx;														   //
	Evx[1][1] = by * by;														   //
	Evx[2][2] = bz * bz;														   //
	VertexParameter vxpar;														   //
	vxpar.setVx(vx);															   //
	vxpar.setEvx(Evx);															   //
	VertexFit *vtxfit = VertexFit::instance();									   // 进行顶点拟合
	vtxfit->init();																   //
	vtxfit->AddTrack(0, wvpipTrk1);												   //
	vtxfit->AddTrack(1, wvpipTrk2);												   //
	vtxfit->AddTrack(2, wvpimTrk1);												   //
	vtxfit->AddTrack(3, wvpimTrk2);												   //
	vtxfit->AddVertex(0, vxpar, 0, 1);											   //
	if (!vtxfit->Fit(0))														   //
		return SUCCESS;															   //
	vtxfit->Swim(0);															   //
	//*********************************************************************************
	// Selection 7: 4~5C Selection
	//*********************************************************************************
	int selecto[4][4] = {{0, 2, 1, 3},
						 {0, 3, 1, 2},
						 {1, 2, 0, 3},
						 {1, 3, 0, 2}};
	//*********************************************************************************
	// Selection 7-1: 4C Selection
	//*********************************************************************************
	WTrackParameter wpip1 = vtxfit->wtrk(0);																						 //
	WTrackParameter wpip2 = vtxfit->wtrk(1);																						 //
	WTrackParameter wpim1 = vtxfit->wtrk(2);																						 //
	WTrackParameter wpim2 = vtxfit->wtrk(3);																						 //
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();																		 //
	HepLorentzVector ecms(0.034 * m_energy / 3.097, 0, 0, m_energy);																 //
	HepLorentzVector ptrackp1, n1ptrackp1;																							 //
	HepLorentzVector ptrackp2, n1ptrackp2;																							 //
	HepLorentzVector ptrackm1, n1ptrackm1;																							 //
	HepLorentzVector ptrackm2, n1ptrackm2;																							 //
	HepLorentzVector ptrackn1;																										 //
	HepLorentzVector ptrackn2;																										 //
	double chisq_4c_2g = 9999;																										 //
	double chisq_w = 9999;																											 //
	for (int i1 = 0; i1 < nGam; i1++)																								 // 得到2Gamma-chisq
	{																																 //
		RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();													 //
		for (int i2 = i1; i2 < nGam; i2++)																							 //
		{																															 //
			if (i2 == i1)																											 //
			{																														 //
				continue;																											 //
			}																														 //
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();												 //
			kmfit->init();																											 //
			kmfit->AddTrack(0, wpip1);																								 //
			kmfit->AddTrack(1, wpip2);																								 //
			kmfit->AddTrack(2, wpim1);																								 //
			kmfit->AddTrack(3, wpim2);																								 //
			kmfit->AddTrack(4, 0.0, g1Trk);																							 //
			kmfit->AddTrack(5, 0.0, g2Trk);																							 //
			kmfit->AddFourMomentum(0, ecms);																						 //
			bool oksq = kmfit->Fit();																								 //
			if (oksq)																												 //
			{																														 //
				double chi2 = kmfit->chisq();																						 //
				if (chi2 <= chisq_4c_2g)																							 //
				{																													 //
					chisq_4c_2g = chi2;																								 //
					ptrackp1 = kmfit->pfit(0);																						 //
					ptrackp2 = kmfit->pfit(1);																						 //
					ptrackm1 = kmfit->pfit(2);																						 //
					ptrackm2 = kmfit->pfit(3);																						 //
					ptrackn1 = kmfit->pfit(4);																						 //
					ptrackn2 = kmfit->pfit(5);																						 //
				}																													 //
			}																														 //
		}																															 //
	}																																 //
	int g2 = 0;																														 //
	if (chisq_4c_2g < 200)																											 // 2Gamma-chisq < 200
	{																																 //
		g2 = 1;																														 //
	}																																 //
	if (g2 == 1)																													 //
	{																																 //
		if (1 == 1)																													 // 重建w
		{																															 //
			HepLorentzVector ptrackg1[4] = {ptrackp1,																				 //
											ptrackp2,																				 //
											ptrackm1,																				 //
											ptrackm2};																				 //
			for (int i = 0; i < 4; i++)																								 //
			{																														 //
				double chisq_momega = pow((ptrackg1[selecto[i][0]] + ptrackg1[selecto[i][1]] + ptrackn1 + ptrackn2).m() - 0.782, 2); //
				if (chisq_momega < chisq_w)																							 //
				{																													 //
					chisq_w = chisq_momega;																							 //
					n1ptrackp1 = ptrackg1[selecto[i][0]];																			 //
					n1ptrackm1 = ptrackg1[selecto[i][1]];																			 //
					n1ptrackp2 = ptrackg1[selecto[i][2]];																			 //
					n1ptrackm2 = ptrackg1[selecto[i][3]];																			 //
				}																													 //
			}																														 //
		}																															 //
		HepLorentzVector out_pi0 = ptrackn1 + ptrackn2;																				 //
		HepLorentzVector out_pip1 = n1ptrackp1;																						 //
		HepLorentzVector out_pip2 = n1ptrackp2;																						 //
		HepLorentzVector out_pim1 = n1ptrackm1;																						 //
		HepLorentzVector out_pim2 = n1ptrackm2;																						 //
		HepLorentzVector out_omega = n1ptrackp1 + n1ptrackm1 + out_pi0;																 //
		HepLorentzVector out_omegapip2 = out_omega + n1ptrackp2;																	 //
		HepLorentzVector out_omegapim2 = out_omega + n1ptrackm2;																	 //
		HepLorentzVector out_pip2pim2 = n1ptrackp2 + n1ptrackm2;																	 //
		if (1 == 1)																													 //
		{
			t4_chisq_4c = chisq_4c_2g;
			// pip2 pim2
			t4_mpip2 = out_pip2.m();
			t4_apip2 = cos(out_pip2.theta());
			t4_ppip2 = out_pip2.rho();
			t4_epip2 = out_pip2.e();
			t4_pxpip2 = out_pip2.px();
			t4_pypip2 = out_pip2.py();
			t4_pzpip2 = out_pip2.pz();
			t4_mpim2 = out_pim2.m();
			t4_apim2 = cos(out_pim2.theta());
			t4_ppim2 = out_pim2.rho();
			t4_epim2 = out_pim2.e();
			t4_pxpim2 = out_pim2.px();
			t4_pypim2 = out_pim2.py();
			t4_pzpim2 = out_pim2.pz();
			// pi0 pip1 pim1 w wpip2 wpim2 pip2pip3
			t4_mpi0 = out_pi0.m();
			t4_api0 = cos(out_pi0.theta());
			t4_ppi0 = out_pi0.rho();
			t4_epi0 = out_pi0.e();
			t4_pxpi0 = out_pi0.px();
			t4_pypi0 = out_pi0.py();
			t4_pzpi0 = out_pi0.pz();
			t4_mpip1 = out_pip1.m();
			t4_apip1 = cos(out_pip1.theta());
			t4_ppip1 = out_pip1.rho();
			t4_epip1 = out_pip1.e();
			t4_pxpip1 = out_pip1.px();
			t4_pypip1 = out_pip1.py();
			t4_pzpip1 = out_pip1.pz();
			t4_mpim1 = out_pim1.m();
			t4_apim1 = cos(out_pim1.theta());
			t4_ppim1 = out_pim1.rho();
			t4_epim1 = out_pim1.e();
			t4_pxpim1 = out_pim1.px();
			t4_pypim1 = out_pim1.py();
			t4_pzpim1 = out_pim1.pz();
			t4_momega = out_omega.m();
			t4_aomega = cos(out_omega.theta());
			t4_pomega = out_omega.rho();
			t4_eomega = out_omega.e();
			t4_pxomega = out_omega.px();
			t4_pyomega = out_omega.py();
			t4_pzomega = out_omega.pz();
			t4_momegapip2 = out_omegapip2.m();
			t4_aomegapip2 = cos(out_omegapip2.theta());
			t4_pomegapip2 = out_omegapip2.rho();
			t4_momegapim2 = out_omegapim2.m();
			t4_aomegapim2 = cos(out_omegapim2.theta());
			t4_pomegapim2 = out_omegapim2.rho();
			t4_mpip2pim2 = out_pip2pim2.m();
			t4_apip2pim2 = cos(out_pip2pim2.theta());
			t4_ppip2pim2 = out_pip2pim2.rho();
			m_tuple4->write();
			Ncut5++;
		}
	}
	//*********************************************************************************
	// Selection 8: End mission
	//*********************************************************************************
	return StatusCode::SUCCESS;
}
//*********************************************************************************************************
//***                                               finalize                                            ***
//*********************************************************************************************************
StatusCode Pm::finalize()
{
	cout << "Energy:               " << m_energy << endl;
	cout << "total number:         " << Ncut0 << endl;
	cout << "Pass charge:          " << Ncut1 << endl;
	cout << "Pass photon:          " << Ncut2 << endl;
	cout << "Pass Pid:             " << Ncut3 << endl;
	cout << "Pass Vertex:          " << Ncut4 << endl;
	cout << "Pass 4C:              " << Ncut5 << endl;
	for (int i = 0; i < SeriesRun.size(); i++)
	{
		cout << "oooooooo" << SeriesRun[i] << "oooooooo" << SeriesNum[i] << "oooooooo" << endl;
	}
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}
