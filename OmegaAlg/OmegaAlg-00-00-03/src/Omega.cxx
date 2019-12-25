//*********************************************************************************************************
//***                                                代码引用                                            ***
//*********************************************************************************************************
#include "McTruth/McParticle.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "VertexFit/IVertexDbSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"
#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"
#include "EventModel/EventHeader.h"
#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep2Vector;
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
//#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"
#include <vector>
#include "OmegaAlg/Omega.h"
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
Omega::Omega(const std::string &name, ISvcLocator *pSvcLocator) : Algorithm(name, pSvcLocator)
{
	declareProperty("Energy", m_energy = 2.125);
}
//*********************************************************************************************************
//***                                               initialize                                          ***
//*********************************************************************************************************
StatusCode Omega::initialize()
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
				//gamma信息
				if (1 == 1)
				{
					status = m_tuple1->addItem("misr", t1_misr);
					status = m_tuple1->addItem("aisr", t1_aisr);
					status = m_tuple1->addItem("pisr", t1_pisr);
					status = m_tuple1->addItem("mpip", t1_mpip);
					status = m_tuple1->addItem("apip", t1_apip);
					status = m_tuple1->addItem("ppip", t1_ppip);
					status = m_tuple1->addItem("mpim", t1_mpim);
					status = m_tuple1->addItem("apim", t1_apim);
					status = m_tuple1->addItem("ppim", t1_ppim);
				}
				// 粒子信息
				if (1 == 1)
				{
					status = m_tuple1->addItem("mpi01", t1_mpi01);
					status = m_tuple1->addItem("api01", t1_api01);
					status = m_tuple1->addItem("ppi01", t1_ppi01);
					status = m_tuple1->addItem("mpi02", t1_mpi02);
					status = m_tuple1->addItem("api02", t1_api02);
					status = m_tuple1->addItem("ppi02", t1_ppi02);
					status = m_tuple1->addItem("mpi03", t1_mpi03);
					status = m_tuple1->addItem("api03", t1_api03);
					status = m_tuple1->addItem("ppi03", t1_ppi03);
					status = m_tuple1->addItem("momega", t1_momega);
					status = m_tuple1->addItem("aomega", t1_aomega);
					status = m_tuple1->addItem("pomega", t1_pomega);
					status = m_tuple1->addItem("momegapi02", t1_momegapi02);
					status = m_tuple1->addItem("aomegapi02", t1_aomegapi02);
					status = m_tuple1->addItem("pomegapi02", t1_pomegapi02);
					status = m_tuple1->addItem("momegapi03", t1_momegapi03);
					status = m_tuple1->addItem("aomegapi03", t1_aomegapi03);
					status = m_tuple1->addItem("pomegapi03", t1_pomegapi03);
					status = m_tuple1->addItem("mpi02pi03", t1_mpi02pi03);
					status = m_tuple1->addItem("api02pi03", t1_api02pi03);
					status = m_tuple1->addItem("ppi02pi03", t1_ppi02pi03);
				}
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
				status = m_tuple4->addItem("chisq_pi", t4_chisq_3pi);
				// topo 信息
				status = m_tuple4->addItem("runID", runID);
				status = m_tuple4->addItem("eventID", eventID);
				status = m_tuple4->addItem("indexmc", m_idxmc, 0, 100);
				status = m_tuple4->addIndexedItem("pdgid", m_idxmc, m_pdgid);
				status = m_tuple4->addIndexedItem("motheridx", m_idxmc, m_motheridx);
				//gamma信息
				if (1 == 1)
				{
					status = m_tuple4->addItem("mpip", t4_mpip);
					status = m_tuple4->addItem("apip", t4_apip);
					status = m_tuple4->addItem("ppip", t4_ppip);
					status = m_tuple4->addItem("mpim", t4_mpim);
					status = m_tuple4->addItem("apim", t4_apim);
					status = m_tuple4->addItem("ppim", t4_ppim);

					status = m_tuple4->addItem("epip", t4_epip);
					status = m_tuple4->addItem("pxpip", t4_pxpip);
					status = m_tuple4->addItem("pypip", t4_pypip);
					status = m_tuple4->addItem("pzpip", t4_pzpip);
					status = m_tuple4->addItem("epim", t4_epim);
					status = m_tuple4->addItem("pxpim", t4_pxpim);
					status = m_tuple4->addItem("pypim", t4_pypim);
					status = m_tuple4->addItem("pzpim", t4_pzpim);

				}
				// 粒子信息
				if (1 == 1)
				{
					status = m_tuple4->addItem("mpi01", t4_mpi01);
					status = m_tuple4->addItem("api01", t4_api01);
					status = m_tuple4->addItem("ppi01", t4_ppi01);
					status = m_tuple4->addItem("mpi02", t4_mpi02);
					status = m_tuple4->addItem("api02", t4_api02);
					status = m_tuple4->addItem("ppi02", t4_ppi02);
					status = m_tuple4->addItem("mpi03", t4_mpi03);
					status = m_tuple4->addItem("api03", t4_api03);
					status = m_tuple4->addItem("ppi03", t4_ppi03);
					status = m_tuple4->addItem("momega", t4_momega);
					status = m_tuple4->addItem("aomega", t4_aomega);
					status = m_tuple4->addItem("pomega", t4_pomega);

					status = m_tuple4->addItem("epi01", t4_epi01);
					status = m_tuple4->addItem("pxpi01", t4_pxpi01);
					status = m_tuple4->addItem("pypi01", t4_pypi01);
					status = m_tuple4->addItem("pzpi01", t4_pzpi01);
					status = m_tuple4->addItem("epi02", t4_epi02);
					status = m_tuple4->addItem("pxpi02", t4_pxpi02);
					status = m_tuple4->addItem("pypi02", t4_pypi02);
					status = m_tuple4->addItem("pzpi02", t4_pzpi02);
					status = m_tuple4->addItem("epi03", t4_epi03);
					status = m_tuple4->addItem("pxpi03", t4_pxpi03);
					status = m_tuple4->addItem("pypi03", t4_pypi03);
					status = m_tuple4->addItem("pzpi03", t4_pzpi03);
					status = m_tuple4->addItem("eomega", t4_eomega);
					status = m_tuple4->addItem("pxomega", t4_pxomega);
					status = m_tuple4->addItem("pyomega", t4_pyomega);
					status = m_tuple4->addItem("pzomega", t4_pzomega);

					status = m_tuple4->addItem("momegapi02", t4_momegapi02);
					status = m_tuple4->addItem("aomegapi02", t4_aomegapi02);
					status = m_tuple4->addItem("pomegapi02", t4_pomegapi02);
					status = m_tuple4->addItem("momegapi03", t4_momegapi03);
					status = m_tuple4->addItem("aomegapi03", t4_aomegapi03);
					status = m_tuple4->addItem("pomegapi03", t4_pomegapi03);
					status = m_tuple4->addItem("mpi02pi03", t4_mpi02pi03);
					status = m_tuple4->addItem("api02pi03", t4_api02pi03);
					status = m_tuple4->addItem("ppi02pi03", t4_ppi02pi03);
				}
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
StatusCode Omega::execute()																//
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
		HepLorentzVector pip_track;
		HepLorentzVector pim_track;
		HepLorentzVector gamma1_track;
		HepLorentzVector gamma2_track;
		HepLorentzVector gamma3_track;
		HepLorentzVector gamma4_track;
		HepLorentzVector gamma5_track;
		HepLorentzVector gamma6_track;
		int nisr = 0;
		int npip = 0;
		int npim = 0;
		int ngamma1 = 0;
		int ngamma2 = 0;
		int ngamma3 = 0;
		int ngamma4 = 0;
		int ngamma5 = 0;
		int ngamma6 = 0;
		HepLorentzVector pi01_track;
		HepLorentzVector pi02_track, medium_pi02_track;
		HepLorentzVector pi03_track, medium_pi03_track;
		HepLorentzVector omega_track;
		HepLorentzVector omegapi02_track;
		HepLorentzVector omegapi03_track;
		HepLorentzVector pi02pi03_track;
		int npi01 = 0;
		int npi02 = 0;
		int npi03 = 0;
		int pi01_index;
		int pi02_index, medium_pi02_index;
		int pi03_index, medium_pi03_index;
		int nomega = 0;
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
			if ((*iter_mc)->particleProperty() == 211 && ((*iter_mc)->mother()).particleProperty() == 223)
			{
				pip_track = mctrue_track;
				npip += 1;
			}
			// 统计pi-
			if ((*iter_mc)->particleProperty() == -211 && ((*iter_mc)->mother()).particleProperty() == 223)
			{
				pim_track = mctrue_track;
				npim += 1;
			}
			// 统计pi0
			if ((*iter_mc)->particleProperty() == 111)
			{
				if (((*iter_mc)->mother()).particleProperty() == 223)
				{
					pi01_track = mctrue_track;
					pi01_index = mctrue_index;
					npi01 += 1;
				}
				else if (npi02 == 0)
				{
					pi02_track = mctrue_track;
					pi02_index = mctrue_index;
					npi02 += 1;
				}
				else
				{
					pi03_track = mctrue_track;
					pi03_index = mctrue_index;
					npi03 += 1;
				}
			}
		}
		// 定位 pi02 pi03
		if (1 == 1)
		{
			omegapi02_track = omega_track + pi02_track;
			omegapi03_track = omega_track + pi03_track;
			if (omegapi02_track.m() > omegapi03_track.m())
			{
				medium_pi02_track = pi03_track;
				medium_pi02_index = pi03_index;
				medium_pi03_track = pi02_track;
				medium_pi03_index = pi02_index;
				pi02_track = medium_pi02_track;
				pi02_index = medium_pi02_index;
				pi03_track = medium_pi03_track;
				pi03_index = medium_pi03_index;
			}
			omegapi02_track = omega_track + pi02_track;
			omegapi03_track = omega_track + pi03_track;
			pi02pi03_track = pi02_track + pi03_track;
		}
		//统计gamma
		iter_mc = mcParticleCol->begin();
		for (; iter_mc != mcParticleCol->end(); iter_mc++)
		{
			if (!(*iter_mc)->decayFromGenerator())
				continue;
			HepLorentzVector mctrue_track = (*iter_mc)->initialFourMomentum();
			int mcture_index = (*iter_mc)->trackIndex();
			if ((*iter_mc)->particleProperty() == 22)
			{
				//统计isr gamma
				if (((*iter_mc)->mother()).particleProperty() != 111)
				{
					isr_track = mctrue_track;
					nisr += 1;
				}
				// 统计pi01 gamma
				if (((*iter_mc)->mother()).trackIndex() == pi01_index)
				{
					if (ngamma1 == 0)
					{
						gamma1_track = mctrue_track;
						ngamma1 += 1;
					}
					else
					{
						gamma2_track = mctrue_track;
						ngamma2 += 1;
					}
				}
				// 统计pi02 gamma
				if (((*iter_mc)->mother()).trackIndex() == pi02_index)
				{
					if (ngamma3 == 0)
					{
						gamma3_track = mctrue_track;
						ngamma3 += 1;
					}
					else
					{
						gamma4_track = mctrue_track;
						ngamma4 += 1;
					}
				}
				// 统计pi03 gamma
				if (((*iter_mc)->mother()).trackIndex() == pi03_index)
				{
					if (ngamma5 == 0)
					{
						gamma5_track = mctrue_track;
						ngamma5 += 1;
					}
					else
					{
						gamma6_track = mctrue_track;
						ngamma6 += 1;
					}
				}
			}
		}
		if (npip == 1 && npim == 1 && npi01 == 1 && npi02 == 1 && npi03 == 1 && nomega == 1 && ngamma1 == 1 && ngamma2 == 1 && ngamma3 == 1 && ngamma4 == 1 && ngamma5 == 1 && ngamma6 == 1)
		{
			truth_check = 1;
		}
		//填入信息
		if (truth_check == 1)
		{
			// 输出gamma信息
			if (1 == 1)
			{
				//isr
				t1_misr = isr_track.m();
				t1_aisr = cos(isr_track.theta());
				t1_pisr = isr_track.rho();
				//pip pim
				t1_mpip = pip_track.m();
				t1_apip = cos(pip_track.theta());
				t1_ppip = pip_track.rho();
				t1_mpim = pim_track.m();
				t1_apim = cos(pim_track.theta());
				t1_ppim = pim_track.rho();
			}
			// 输出粒子信息
			if (1 == 1)
			{
				t1_mpi01 = pi01_track.m();
				t1_api01 = cos(pi01_track.theta());
				t1_ppi01 = pi01_track.rho();
				t1_mpi02 = pi02_track.m();
				t1_api02 = cos(pi02_track.theta());
				t1_ppi02 = pi02_track.rho();
				t1_mpi03 = pi03_track.m();
				t1_api03 = cos(pi03_track.theta());
				t1_ppi03 = pi03_track.rho();
				t1_momega = omega_track.m();
				t1_aomega = cos(omega_track.theta());
				t1_pomega = omega_track.rho();
				t1_momegapi03 = omegapi03_track.m();
				t1_aomegapi03 = cos(omegapi03_track.theta());
				t1_pomegapi03 = omegapi03_track.rho();
				t1_momegapi02 = omegapi02_track.m();
				t1_aomegapi02 = cos(omegapi02_track.theta());
				t1_pomegapi02 = omegapi02_track.rho();
				t1_mpi02pi03 = pi02pi03_track.m();
				t1_api02pi03 = cos(pi02pi03_track.theta());
				t1_ppi02pi03 = pi02pi03_track.rho();
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
	if ((nGood != 2) || (nCharge != 0))												   // 选择：nGood
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
	if (nGam < 6 || nGam > 55)																			  // 选择：nGam
	{																									  //
		return StatusCode::SUCCESS;																		  //
	}																									  //
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
	if (npip * npim != 1)																   //
		return SUCCESS;																	   //
	Ncut3++;																			   //
	//*********************************************************************************
	// Selection 3: Vertex fit Selection, check ppi0, pTot
	//*********************************************************************************
	RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin() + ipip[0]))->mdcKalTrack(); // Default is pion, for other particles:
	RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin() + ipim[0]))->mdcKalTrack(); // wvppTrk = WTrackParameter(mp, pipTrk->getZHelixP(), pipTrk->getZErrorP()); proton
	WTrackParameter wvpipTrk, wvpimTrk;											  // wvmupTrk = WTrackParameter(mmu, pipTrk->getZHelixMu(), pipTrk->getZErrorMu()); muon
	wvpipTrk = WTrackParameter(mpi, pipTrk->getZHelix(), pipTrk->getZError());	// wvepTrk = WTrackParameter(me, pipTrk->getZHelixE(), pipTrk->getZErrorE()); electron
	wvpimTrk = WTrackParameter(mpi, pimTrk->getZHelix(), pimTrk->getZError());	// wvkpTrk = WTrackParameter(mk, pipTrk->getZHelixK(), pipTrk->getZErrorK()); kaon
	HepPoint3D vx(0., 0., 0.);													  //
	HepSymMatrix Evx(3, 0);														  //
	double bx = 1E+6;															  //
	double by = 1E+6;															  //
	double bz = 1E+6;															  //
	Evx[0][0] = bx * bx;														  //
	Evx[1][1] = by * by;														  //
	Evx[2][2] = bz * bz;														  //
	VertexParameter vxpar;														  //
	vxpar.setVx(vx);															  //
	vxpar.setEvx(Evx);															  //
	VertexFit *vtxfit = VertexFit::instance();									  //
	vtxfit->init();																  //
	vtxfit->AddTrack(0, wvpipTrk);												  // 设定track0
	vtxfit->AddTrack(1, wvpimTrk);												  // 设定track1
	vtxfit->AddVertex(0, vxpar, 0, 1);											  // 设定顶点0
	if (!vtxfit->Fit(0))														  //
		return SUCCESS;															  //
	vtxfit->Swim(0);															  //
	//*********************************************************************************
	// Selection 7: 4 Selection
	//*********************************************************************************
	int combine[15][6] = {{1, 2, 3, 4, 5, 6},  //
						  {1, 2, 3, 5, 4, 6},  //
						  {1, 2, 3, 6, 4, 5},  //
						  {1, 3, 2, 4, 5, 6},  //
						  {1, 3, 2, 5, 4, 6},  //
						  {1, 3, 2, 6, 4, 5},  //
						  {1, 4, 3, 2, 5, 6},  //
						  {1, 4, 3, 5, 2, 6},  //
						  {1, 4, 3, 6, 2, 5},  //
						  {1, 5, 3, 4, 2, 6},  //
						  {1, 5, 3, 2, 4, 6},  //
						  {1, 5, 3, 6, 4, 2},  //
						  {1, 6, 3, 4, 5, 2},  //
						  {1, 6, 3, 5, 4, 2},  //
						  {1, 6, 3, 2, 4, 5}}; //
	int selecto[3][3] = {{1, 2, 3},
						 {2, 3, 1},
						 {3, 1, 2}};
	int selectb[2][2] = {{1, 2},
						 {2, 1}};
	//*********************************************************************************
	// Selection 7-1: 4C Selection
	//*********************************************************************************
	WTrackParameter wpip = vtxfit->wtrk(0);																			  //
	WTrackParameter wpim = vtxfit->wtrk(1);																			  //
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();														  //
	HepLorentzVector ecms(0.034 * m_energy / 3.097, 0, 0, m_energy);												  //
	HepLorentzVector ptrackp;																						  //
	HepLorentzVector ptrackm;																						  //
	HepLorentzVector ptrack1, n1ptrack1, n2ptrack1, n3ptrack1;														  //
	HepLorentzVector ptrack2, n1ptrack2, n2ptrack2, n3ptrack2;														  //
	HepLorentzVector ptrack3, n1ptrack3, n2ptrack3, n3ptrack3;														  //
	HepLorentzVector ptrack4, n1ptrack4, n2ptrack4, n3ptrack4;														  //
	HepLorentzVector ptrack5, n1ptrack5, n2ptrack5, n3ptrack5;														  //
	HepLorentzVector ptrack6, n1ptrack6, n2ptrack6, n3ptrack6;														  //
	double chisq_4c_5g = 9999;																						  //
	double chisq_4c_6g = 9999;																						  //
	double chisq_4c_7g = 9999;																						  //
	double chisq_4c_pi = 9999;																						  //
	double chisq_4c_om = 9999;																						  //
	double chisq_4c_b0 = 9999;																						  //
	for (int i1 = 0; i1 < nGam; i1++)																				  // 得到6Gamma-chisq
	{																												  //
		RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();									  //
		for (int i2 = i1; i2 < nGam; i2++)																			  //
		{																											  //
			if (i2 == i1)																							  //
			{																										  //
				continue;																							  //
			}																										  //
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();								  //
			for (int i3 = i2; i3 < nGam; i3++)																		  //
			{																										  //
				if (i3 == i1 || i3 == i2)																			  //
				{																									  //
					continue;																						  //
				}																									  //
				RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[i3]))->emcShower();							  //
				for (int i4 = i3; i4 < nGam; i4++)																	  //
				{																									  //
					if (i4 == i1 || i4 == i2 || i4 == i3)															  //
					{																								  //
						continue;																					  //
					}																								  //
					RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[i4]))->emcShower();						  //
					for (int i5 = i4; i5 < nGam; i5++)																  //
					{																								  //
						if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4)											  //
						{																							  //
							continue;																				  //
						}																							  //
						RecEmcShower *g5Trk = (*(evtRecTrkCol->begin() + iGam[i5]))->emcShower();					  //
						for (int i6 = i5; i6 < nGam; i6++)															  //
						{																							  //
							if (i6 == i1 || i6 == i2 || i6 == i3 || i6 == i4 || i6 == i5)							  //
							{																						  //
								continue;																			  //
							}																						  //
							RecEmcShower *g6Trk = (*(evtRecTrkCol->begin() + iGam[i6]))->emcShower();				  //
							kmfit->init();																			  //
							kmfit->AddTrack(0, wpip);																  //
							kmfit->AddTrack(1, wpim);																  //
							kmfit->AddTrack(2, 0.0, g1Trk);															  //
							kmfit->AddTrack(3, 0.0, g2Trk);															  //
							kmfit->AddTrack(4, 0.0, g3Trk);															  //
							kmfit->AddTrack(5, 0.0, g4Trk);															  //
							kmfit->AddTrack(6, 0.0, g5Trk);															  //
							kmfit->AddTrack(7, 0.0, g6Trk);															  //
							kmfit->AddFourMomentum(0, ecms);														  //
							bool oksq = kmfit->Fit();																  //
							if (oksq)																				  //
							{																						  //
								double chi2 = kmfit->chisq();														  //
								if (chi2 <= chisq_4c_6g)															  // 选择：最小chi-4c
								{																					  //
									chisq_4c_6g = chi2;																  //
									ptrackp = kmfit->pfit(0);														  //
									ptrackm = kmfit->pfit(1);														  //
									ptrack1 = kmfit->pfit(2);														  //
									ptrack2 = kmfit->pfit(3);														  //
									ptrack3 = kmfit->pfit(4);														  //
									ptrack4 = kmfit->pfit(5);														  //
									ptrack5 = kmfit->pfit(6);														  //
									ptrack6 = kmfit->pfit(7);														  //
								}																					  //
							}																						  //
						}																							  //
					}																								  //
				}																									  //
			}																										  //
		}																											  //
	}																												  //
	for (int i1 = 0; i1 < nGam; i1++)																				  // 得到5Gamma-chisq
	{																												  //
		RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();									  //
		for (int i2 = i1; i2 < nGam; i2++)																			  //
		{																											  //
			if (i2 == i1)																							  //
			{																										  //
				continue;																							  //
			}																										  //
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();								  //
			for (int i3 = i2; i3 < nGam; i3++)																		  //
			{																										  //
				if (i3 == i1 || i3 == i2)																			  //
				{																									  //
					continue;																						  //
				}																									  //
				RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[i3]))->emcShower();							  //
				for (int i4 = i3; i4 < nGam; i4++)																	  //
				{																									  //
					if (i4 == i1 || i4 == i2 || i4 == i3)															  //
					{																								  //
						continue;																					  //
					}																								  //
					RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[i4]))->emcShower();						  //
					for (int i5 = i4; i5 < nGam; i5++)																  //
					{																								  //
						if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4)											  //
						{																							  //
							continue;																				  //
						}																							  //
						RecEmcShower *g5Trk = (*(evtRecTrkCol->begin() + iGam[i5]))->emcShower();					  //
						kmfit->init();																				  //
						kmfit->AddTrack(0, wpip);																	  //
						kmfit->AddTrack(1, wpim);																	  //
						kmfit->AddTrack(2, 0.0, g1Trk);																  //
						kmfit->AddTrack(3, 0.0, g2Trk);																  //
						kmfit->AddTrack(4, 0.0, g3Trk);																  //
						kmfit->AddTrack(5, 0.0, g4Trk);																  //
						kmfit->AddTrack(6, 0.0, g5Trk);																  //
						kmfit->AddFourMomentum(0, ecms);															  //
						bool oksq = kmfit->Fit();																	  //
						if (oksq)																					  //
						{																							  //
							double chi2 = kmfit->chisq();															  //
							if (chi2 < chisq_4c_5g)																	  //
							{																						  //
								chisq_4c_5g = chi2;																	  //
							}																						  //
						}																							  //
					}																								  //
				}																									  //
			}																										  //
		}																											  //
	}																												  //
	if (nGam > 6)																									  // 得到7Gamma-chisq
	{																												  //
		for (int i1 = 0; i1 < nGam; i1++)																			  //
		{																											  //
			RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();								  //
			for (int i2 = i1; i2 < nGam; i2++)																		  //
			{																										  //
				if (i2 == i1)																						  //
				{																									  //
					continue;																						  //
				}																									  //
				RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();							  //
				for (int i3 = i2; i3 < nGam; i3++)																	  //
				{																									  //
					if (i3 == i1 || i3 == i2)																		  //
					{																								  //
						continue;																					  //
					}																								  //
					RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[i3]))->emcShower();						  //
					for (int i4 = i3; i4 < nGam; i4++)																  //
					{																								  //
						if (i4 == i1 || i4 == i2 || i4 == i3)														  //
						{																							  //
							continue;																				  //
						}																							  //
						RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[i4]))->emcShower();					  //
						for (int i5 = i4; i5 < nGam; i5++)															  //
						{																							  //
							if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4)										  //
							{																						  //
								continue;																			  //
							}																						  //
							RecEmcShower *g5Trk = (*(evtRecTrkCol->begin() + iGam[i5]))->emcShower();				  //
							for (int i6 = i5; i6 < nGam; i6++)														  //
							{																						  //
								if (i6 == i1 || i6 == i2 || i6 == i3 || i6 == i4 || i6 == i5)						  //
								{																					  //
									continue;																		  //
								}																					  //
								RecEmcShower *g6Trk = (*(evtRecTrkCol->begin() + iGam[i6]))->emcShower();			  //
								for (int i7 = i6; i7 < nGam; i7++)													  //
								{																					  //
									if (i7 == i1 || i7 == i2 || i7 == i3 || i7 == i4 || i7 == i5 || i7 == i6)		  //
									{																				  //
										continue;																	  //
									}																				  //
									RecEmcShower *g7Trk = (*(evtRecTrkCol->begin() + iGam[i7]))->emcShower();		  //
									kmfit->init();																	  //
									kmfit->AddTrack(0, wpip);														  //
									kmfit->AddTrack(1, wpim);														  //
									kmfit->AddTrack(2, 0.0, g1Trk);													  //
									kmfit->AddTrack(3, 0.0, g2Trk);													  //
									kmfit->AddTrack(4, 0.0, g3Trk);													  //
									kmfit->AddTrack(5, 0.0, g4Trk);													  //
									kmfit->AddTrack(6, 0.0, g5Trk);													  //
									kmfit->AddTrack(7, 0.0, g6Trk);													  //
									kmfit->AddTrack(8, 0.0, g7Trk);													  //
									kmfit->AddFourMomentum(0, ecms);												  //
									bool oksq = kmfit->Fit();														  //
									if (oksq)																		  //
									{																				  //
										double chi2 = kmfit->chisq();												  //
										if (chi2 < chisq_4c_7g)														  //
										{																			  //
											chisq_4c_7g = chi2;														  //
										}																			  //
									}																				  //
								}																					  //
							}																						  //
						}																							  //
					}																								  //
				}																									  //
			}																										  //
		}																											  //
	}																												  //
	int g6 = 0;																										  //
	if (chisq_4c_6g < 200)																							  // 选择：chisq-6Gamma
	{																												  //
		g6 = 1;																										  //
	}																												  //
	if (chisq_4c_5g < chisq_4c_6g)																					  // 选择：chisq-5Gamma
	{																												  //
		g6 = 0;																										  //
	}																												  //
	if (nGam > 6)																									  // 选择：chisq-7Gamma
	{																												  //
		if (chisq_4c_7g < chisq_4c_6g)																				  //
		{																											  //
			g6 = 0;																									  //
		}																											  //
	}																												  //
	if (g6 == 1)																									  //
	{																												  //
		if (1 == 1)																									  //
		{																											  //
			HepLorentzVector ptrackg1[6] = {ptrack1,																  //
											ptrack2,																  //
											ptrack3,																  //
											ptrack4,																  //
											ptrack5,																  //
											ptrack6};																  //
			for (int i = 0; i < 15; i++)																			  //
			{																										  //
				double chisq_mpi01 = pow((ptrackg1[combine[i][0] - 1] + ptrackg1[combine[i][1] - 1]).m() - 0.135, 2); //
				double chisq_mpi02 = pow((ptrackg1[combine[i][2] - 1] + ptrackg1[combine[i][3] - 1]).m() - 0.135, 2); //
				double chisq_mpi03 = pow((ptrackg1[combine[i][4] - 1] + ptrackg1[combine[i][5] - 1]).m() - 0.135, 2); //
				double chisq_mpi0 = (chisq_mpi01 + chisq_mpi02 + chisq_mpi03) / 3;									  //
				if (chisq_mpi0 < chisq_4c_pi)																		  //
				{																									  //
					chisq_4c_pi = chisq_mpi0;																		  //
					n1ptrack1 = ptrackg1[combine[i][0] - 1];														  //
					n1ptrack2 = ptrackg1[combine[i][1] - 1];														  //
					n1ptrack3 = ptrackg1[combine[i][2] - 1];														  //
					n1ptrack4 = ptrackg1[combine[i][3] - 1];														  //
					n1ptrack5 = ptrackg1[combine[i][4] - 1];														  //
					n1ptrack6 = ptrackg1[combine[i][5] - 1];														  //
				}																									  //
			}																										  //
		}																											  //
		if (1 == 1)																									  //
		{																											  //
			HepLorentzVector ptrackg2[6] = {n1ptrack1,																  //
											n1ptrack2,																  //
											n1ptrack3,																  //
											n1ptrack4,																  //
											n1ptrack5,																  //
											n1ptrack6};																  //
			for (int i = 0; i < 3; i++)
			{
				double chisq_momega = pow((ptrackp + ptrackm + ptrackg2[2 * selecto[i][0] - 2] + ptrackg2[2 * selecto[i][0] - 1]).m() - 0.782, 2);
				if (chisq_momega < chisq_4c_om)
				{
					chisq_4c_om = chisq_momega;
					n2ptrack1 = ptrackg2[2 * selecto[i][0] - 2];
					n2ptrack2 = ptrackg2[2 * selecto[i][0] - 1];
					n2ptrack3 = ptrackg2[2 * selecto[i][1] - 2];
					n2ptrack4 = ptrackg2[2 * selecto[i][1] - 1];
					n2ptrack5 = ptrackg2[2 * selecto[i][2] - 2];
					n2ptrack6 = ptrackg2[2 * selecto[i][2] - 1];
				}
			}
		}
		if (1 == 1)
		{
			HepLorentzVector ptrackg3[6] = {n2ptrack1,
											n2ptrack2,
											n2ptrack3,
											n2ptrack4,
											n2ptrack5,
											n2ptrack6};
			for (int i = 0; i < 2; i++)
			{
				double chisq_mb = pow((ptrackp + ptrackm + ptrackg3[0] + ptrackg3[1] + ptrackg3[2 * selectb[i][0]] + ptrackg3[2 * selectb[i][0] + 1]).m(), 2);
				if (chisq_mb < chisq_4c_b0)
				{
					chisq_4c_b0 = chisq_mb;
					n3ptrack1 = ptrackg3[0];
					n3ptrack2 = ptrackg3[1];
					n3ptrack3 = ptrackg3[2 * selectb[i][0] + 0];
					n3ptrack4 = ptrackg3[2 * selectb[i][0] + 1];
					n3ptrack5 = ptrackg3[2 * selectb[i][1] + 0];
					n3ptrack6 = ptrackg3[2 * selectb[i][1] + 1];
				}
			}
		}
		// gamma 信息
		HepLorentzVector out_pip = ptrackp;
		HepLorentzVector out_pim = ptrackm;
		HepLorentzVector out_gamma1 = n3ptrack1;
		HepLorentzVector out_gamma2 = n3ptrack2;
		HepLorentzVector out_gamma3 = n3ptrack3;
		HepLorentzVector out_gamma4 = n3ptrack4;
		HepLorentzVector out_gamma5 = n3ptrack5;
		HepLorentzVector out_gamma6 = n3ptrack6;
		// 粒子信息
		HepLorentzVector out_pi01 = out_gamma1 + out_gamma2;
		HepLorentzVector out_pi02 = out_gamma3 + out_gamma4;
		HepLorentzVector out_pi03 = out_gamma5 + out_gamma6;
		HepLorentzVector out_omega = out_pip + out_pim + out_pi01;
		HepLorentzVector out_omegapi02 = out_omega + out_pi02;
		HepLorentzVector out_omegapi03 = out_omega + out_pi03;
		HepLorentzVector out_pi02pi03 = out_pi02 + out_pi03;
		// 填入信息
		if (1 == 1)
		{
			t4_chisq_4c = chisq_4c_6g;
			t4_chisq_3pi = chisq_4c_pi;
			// 输出gamma信息
			if (1 == 1)
			{
				//pip pim
				t4_mpip = out_pip.m();
				t4_apip = cos(out_pip.theta());
				t4_ppip = out_pip.rho();
				t4_mpim = out_pim.m();
				t4_apim = cos(out_pim.theta());
				t4_ppim = out_pim.rho();
				//4-momentum of pip pim
				t4_epip = out_pip.e();
				t4_pxpip = out_pip.px();
				t4_pypip = out_pip.py();
				t4_pzpip = out_pip.pz();
				t4_epim = out_pim.e();
				t4_pxpim = out_pim.px();
				t4_pypim = out_pim.py();
				t4_pzpim = out_pim.pz();
			}
			// 输出粒子信息
			if (1 == 1)
			{
				//
				t4_mpi01 = out_pi01.m();
				t4_api01 = cos(out_pi01.theta());
				t4_ppi01 = out_pi01.rho();
				t4_mpi02 = out_pi02.m();
				t4_api02 = cos(out_pi02.theta());
				t4_ppi02 = out_pi02.rho();
				t4_mpi03 = out_pi03.m();
				t4_api03 = cos(out_pi03.theta());
				t4_ppi03 = out_pi03.rho();
				t4_momega = out_omega.m();
				t4_aomega = cos(out_omega.theta());
				t4_pomega = out_omega.rho();
				//
				t4_epi01 = out_pi01.e();
				t4_pxpi01 = out_pi01.px();
				t4_pypi01 = out_pi01.py();
				t4_pzpi01 = out_pi01.pz();
				t4_epi02 = out_pi02.e();
				t4_pxpi02 = out_pi02.px();
				t4_pypi02 = out_pi02.py();
				t4_pzpi02 = out_pi02.pz();
				t4_epi03 = out_pi01.e();
				t4_pxpi03 = out_pi03.px();
				t4_pypi03 = out_pi03.py();
				t4_pzpi03 = out_pi03.pz();
				t4_eomega = out_omega.e();
				t4_pxomega = out_omega.px();
				t4_pyomega = out_omega.py();
				t4_pzomega = out_omega.pz();
				//
				t4_momegapi03 = out_omegapi03.m();
				t4_aomegapi03 = cos(out_omegapi03.theta());
				t4_pomegapi03 = out_omegapi03.rho();
				t4_momegapi02 = out_omegapi02.m();
				t4_aomegapi02 = cos(out_omegapi02.theta());
				t4_pomegapi02 = out_omegapi02.rho();
				t4_mpi02pi03 = out_pi02pi03.m();
				t4_api02pi03 = cos(out_pi02pi03.theta());
				t4_ppi02pi03 = out_pi02pi03.rho();
			}
			m_tuple4->write();
			Ncut4++;
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
StatusCode Omega::finalize()
{
	cout << "能量为                 " << m_energy << endl;
	cout << "total number:         " << Ncut0 << endl;
	cout << "Pass truth:           " << Ncut1 << endl;
	cout << "Pass Pid:             " << Ncut3 << endl;
	cout << "Pass 4C:              " << Ncut4 << endl;
	for (int i = 0; i < SeriesRun.size(); i++)
	{
		cout << "oooooooo" << SeriesRun[i] << "oooooooo" << SeriesNum[i] << "oooooooo" << endl;
	}
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}
