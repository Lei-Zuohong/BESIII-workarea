#pragma region 调用头文件
// include common
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
#pragma endregion
#pragma region 调用类型对象
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
using CLHEP::Hep2Vector;
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
#include <vector>
#include "OmegaAlg/Omega.h"
#include "OmegaAlg/Myfunc.h"
#pragma endregion
#pragma region 定义全局变量
// 定义常数
const double mpi = 0.13957;
// 定义调用类型
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
// 定义统计run的参数
std::vector<int> SeriesRun;
std::vector<int> SeriesNum;
int number, checkexit, checki, firstrun;
// 定义计数参数
int Ncut0, Ncut1, Ncut2, Ncut3, Ncut4, Ncut5;
#pragma endregion
#pragma region 引入变量容器
Omega::Omega(const std::string &name, ISvcLocator *pSvcLocator) : Algorithm(name, pSvcLocator)
{
	declareProperty("Energy", m_energy = 0);
}
#pragma endregion
#pragma region 初始化
StatusCode Omega::initialize()
{
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in initialize()" << endmsg;
	StatusCode status;
#pragma region 新建树：truth
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
			status = m_tuple1->addItem("misr", truth_misr);
			status = m_tuple1->addItem("aisr", truth_aisr);
			status = m_tuple1->addItem("pisr", truth_pisr);
			status = m_tuple1->addItem("eisr", truth_eisr);
			status = m_tuple1->addItem("pxisr", truth_pxisr);
			status = m_tuple1->addItem("pyisr", truth_pyisr);
			status = m_tuple1->addItem("pzisr", truth_pzisr);
			// pip
			status = m_tuple1->addItem("mpip", truth_mpip);
			status = m_tuple1->addItem("apip", truth_apip);
			status = m_tuple1->addItem("ppip", truth_ppip);
			status = m_tuple1->addItem("epip", truth_epip);
			status = m_tuple1->addItem("pxpip", truth_pxpip);
			status = m_tuple1->addItem("pypip", truth_pypip);
			status = m_tuple1->addItem("pzpip", truth_pzpip);
			// pim
			status = m_tuple1->addItem("mpim", truth_mpim);
			status = m_tuple1->addItem("apim", truth_apim);
			status = m_tuple1->addItem("ppim", truth_ppim);
			status = m_tuple1->addItem("epim", truth_epim);
			status = m_tuple1->addItem("pxpim", truth_pxpim);
			status = m_tuple1->addItem("pypim", truth_pypim);
			status = m_tuple1->addItem("pzpim", truth_pzpim);
			// gamma1
			status = m_tuple1->addItem("mgamma1", truth_mgamma1);
			status = m_tuple1->addItem("agamma1", truth_agamma1);
			status = m_tuple1->addItem("pgamma1", truth_pgamma1);
			status = m_tuple1->addItem("egamma1", truth_egamma1);
			status = m_tuple1->addItem("pxgamma1", truth_pxgamma1);
			status = m_tuple1->addItem("pygamma1", truth_pygamma1);
			status = m_tuple1->addItem("pzgamma1", truth_pzgamma1);
			// gamma2
			status = m_tuple1->addItem("mgamma2", truth_mgamma2);
			status = m_tuple1->addItem("agamma2", truth_agamma2);
			status = m_tuple1->addItem("pgamma2", truth_pgamma2);
			status = m_tuple1->addItem("egamma2", truth_egamma2);
			status = m_tuple1->addItem("pxgamma2", truth_pxgamma2);
			status = m_tuple1->addItem("pygamma2", truth_pygamma2);
			status = m_tuple1->addItem("pzgamma2", truth_pzgamma2);
			// gamma3
			status = m_tuple1->addItem("mgamma3", truth_mgamma3);
			status = m_tuple1->addItem("agamma3", truth_agamma3);
			status = m_tuple1->addItem("pgamma3", truth_pgamma3);
			status = m_tuple1->addItem("egamma3", truth_egamma3);
			status = m_tuple1->addItem("pxgamma3", truth_pxgamma3);
			status = m_tuple1->addItem("pygamma3", truth_pygamma3);
			status = m_tuple1->addItem("pzgamma3", truth_pzgamma3);
			// gamma4
			status = m_tuple1->addItem("mgamma4", truth_mgamma4);
			status = m_tuple1->addItem("agamma4", truth_agamma4);
			status = m_tuple1->addItem("pgamma4", truth_pgamma4);
			status = m_tuple1->addItem("egamma4", truth_egamma4);
			status = m_tuple1->addItem("pxgamma4", truth_pxgamma4);
			status = m_tuple1->addItem("pygamma4", truth_pygamma4);
			status = m_tuple1->addItem("pzgamma4", truth_pzgamma4);
			// gamma5
			status = m_tuple1->addItem("mgamma5", truth_mgamma5);
			status = m_tuple1->addItem("agamma5", truth_agamma5);
			status = m_tuple1->addItem("pgamma5", truth_pgamma5);
			status = m_tuple1->addItem("egamma5", truth_egamma5);
			status = m_tuple1->addItem("pxgamma5", truth_pxgamma5);
			status = m_tuple1->addItem("pygamma5", truth_pygamma5);
			status = m_tuple1->addItem("pzgamma5", truth_pzgamma5);
			// gamma6
			status = m_tuple1->addItem("mgamma6", truth_mgamma6);
			status = m_tuple1->addItem("agamma6", truth_agamma6);
			status = m_tuple1->addItem("pgamma6", truth_pgamma6);
			status = m_tuple1->addItem("egamma6", truth_egamma6);
			status = m_tuple1->addItem("pxgamma6", truth_pxgamma6);
			status = m_tuple1->addItem("pygamma6", truth_pygamma6);
			status = m_tuple1->addItem("pzgamma6", truth_pzgamma6);
			// pi01
			status = m_tuple1->addItem("mpi01", truth_mpi01);
			status = m_tuple1->addItem("api01", truth_api01);
			status = m_tuple1->addItem("ppi01", truth_ppi01);
			status = m_tuple1->addItem("epi01", truth_epi01);
			status = m_tuple1->addItem("pxpi01", truth_pxpi01);
			status = m_tuple1->addItem("pypi01", truth_pypi01);
			status = m_tuple1->addItem("pzpi01", truth_pzpi01);
			// pi02
			status = m_tuple1->addItem("mpi02", truth_mpi02);
			status = m_tuple1->addItem("api02", truth_api02);
			status = m_tuple1->addItem("ppi02", truth_ppi02);
			status = m_tuple1->addItem("epi02", truth_epi02);
			status = m_tuple1->addItem("pxpi02", truth_pxpi02);
			status = m_tuple1->addItem("pypi02", truth_pypi02);
			status = m_tuple1->addItem("pzpi02", truth_pzpi02);
			// pi03
			status = m_tuple1->addItem("mpi03", truth_mpi03);
			status = m_tuple1->addItem("api03", truth_api03);
			status = m_tuple1->addItem("ppi03", truth_ppi03);
			status = m_tuple1->addItem("epi03", truth_epi03);
			status = m_tuple1->addItem("pxpi03", truth_pxpi03);
			status = m_tuple1->addItem("pypi03", truth_pypi03);
			status = m_tuple1->addItem("pzpi03", truth_pzpi03);
			// omega
			status = m_tuple1->addItem("momega", truth_momega);
			status = m_tuple1->addItem("aomega", truth_aomega);
			status = m_tuple1->addItem("pomega", truth_pomega);
			status = m_tuple1->addItem("eomega", truth_eomega);
			status = m_tuple1->addItem("pxomega", truth_pxomega);
			status = m_tuple1->addItem("pyomega", truth_pyomega);
			status = m_tuple1->addItem("pzomega", truth_pzomega);
			// omegapi02
			status = m_tuple1->addItem("momegapi02", truth_momegapi02);
			status = m_tuple1->addItem("aomegapi02", truth_aomegapi02);
			status = m_tuple1->addItem("pomegapi02", truth_pomegapi02);
			// omegapi03
			status = m_tuple1->addItem("momegapi03", truth_momegapi03);
			status = m_tuple1->addItem("aomegapi03", truth_aomegapi03);
			status = m_tuple1->addItem("pomegapi03", truth_pomegapi03);
			// pi02pi03
			status = m_tuple1->addItem("mpi02pi03", truth_mpi02pi03);
			status = m_tuple1->addItem("api02pi03", truth_api02pi03);
			status = m_tuple1->addItem("ppi02pi03", truth_ppi02pi03);
			// boost
			status = m_tuple1->addItem("hgamma1", truth_hgamma1);
			status = m_tuple1->addItem("hgamma2", truth_hgamma2);
			status = m_tuple1->addItem("hgamma3", truth_hgamma3);
			status = m_tuple1->addItem("hgamma4", truth_hgamma4);
			status = m_tuple1->addItem("hgamma5", truth_hgamma5);
			status = m_tuple1->addItem("hgamma6", truth_hgamma6);
		}
		else
		{
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1) << endmsg;
			return StatusCode::FAILURE;
		}
	}
#pragma endregion
#pragma region 新建树：fit4c
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
			status = m_tuple4->addItem("chisq", fit4c_chisq_4c);
			// topo 信息
			status = m_tuple4->addItem("runID", runID);
			status = m_tuple4->addItem("eventID", eventID);
			status = m_tuple4->addItem("indexmc", m_idxmc, 0, 100);
			status = m_tuple4->addIndexedItem("pdgid", m_idxmc, m_pdgid);
			status = m_tuple4->addIndexedItem("motheridx", m_idxmc, m_motheridx);
			// pip
			status = m_tuple4->addItem("mpip", fit4c_mpip);
			status = m_tuple4->addItem("apip", fit4c_apip);
			status = m_tuple4->addItem("ppip", fit4c_ppip);
			status = m_tuple4->addItem("epip", fit4c_epip);
			status = m_tuple4->addItem("pxpip", fit4c_pxpip);
			status = m_tuple4->addItem("pypip", fit4c_pypip);
			status = m_tuple4->addItem("pzpip", fit4c_pzpip);
			// pim
			status = m_tuple4->addItem("mpim", fit4c_mpim);
			status = m_tuple4->addItem("apim", fit4c_apim);
			status = m_tuple4->addItem("ppim", fit4c_ppim);
			status = m_tuple4->addItem("epim", fit4c_epim);
			status = m_tuple4->addItem("pxpim", fit4c_pxpim);
			status = m_tuple4->addItem("pypim", fit4c_pypim);
			status = m_tuple4->addItem("pzpim", fit4c_pzpim);
			// gamma1
			status = m_tuple4->addItem("mgamma1", fit4c_mgamma1);
			status = m_tuple4->addItem("agamma1", fit4c_agamma1);
			status = m_tuple4->addItem("pgamma1", fit4c_pgamma1);
			status = m_tuple4->addItem("egamma1", fit4c_egamma1);
			status = m_tuple4->addItem("pxgamma1", fit4c_pxgamma1);
			status = m_tuple4->addItem("pygamma1", fit4c_pygamma1);
			status = m_tuple4->addItem("pzgamma1", fit4c_pzgamma1);
			// gamma2
			status = m_tuple4->addItem("mgamma2", fit4c_mgamma2);
			status = m_tuple4->addItem("agamma2", fit4c_agamma2);
			status = m_tuple4->addItem("pgamma2", fit4c_pgamma2);
			status = m_tuple4->addItem("egamma2", fit4c_egamma2);
			status = m_tuple4->addItem("pxgamma2", fit4c_pxgamma2);
			status = m_tuple4->addItem("pygamma2", fit4c_pygamma2);
			status = m_tuple4->addItem("pzgamma2", fit4c_pzgamma2);
			// gamma3
			status = m_tuple4->addItem("mgamma3", fit4c_mgamma3);
			status = m_tuple4->addItem("agamma3", fit4c_agamma3);
			status = m_tuple4->addItem("pgamma3", fit4c_pgamma3);
			status = m_tuple4->addItem("egamma3", fit4c_egamma3);
			status = m_tuple4->addItem("pxgamma3", fit4c_pxgamma3);
			status = m_tuple4->addItem("pygamma3", fit4c_pygamma3);
			status = m_tuple4->addItem("pzgamma3", fit4c_pzgamma3);
			// gamma4
			status = m_tuple4->addItem("mgamma4", fit4c_mgamma4);
			status = m_tuple4->addItem("agamma4", fit4c_agamma4);
			status = m_tuple4->addItem("pgamma4", fit4c_pgamma4);
			status = m_tuple4->addItem("egamma4", fit4c_egamma4);
			status = m_tuple4->addItem("pxgamma4", fit4c_pxgamma4);
			status = m_tuple4->addItem("pygamma4", fit4c_pygamma4);
			status = m_tuple4->addItem("pzgamma4", fit4c_pzgamma4);
			// gamma5
			status = m_tuple4->addItem("mgamma5", fit4c_mgamma5);
			status = m_tuple4->addItem("agamma5", fit4c_agamma5);
			status = m_tuple4->addItem("pgamma5", fit4c_pgamma5);
			status = m_tuple4->addItem("egamma5", fit4c_egamma5);
			status = m_tuple4->addItem("pxgamma5", fit4c_pxgamma5);
			status = m_tuple4->addItem("pygamma5", fit4c_pygamma5);
			status = m_tuple4->addItem("pzgamma5", fit4c_pzgamma5);
			// gamma6
			status = m_tuple4->addItem("mgamma6", fit4c_mgamma6);
			status = m_tuple4->addItem("agamma6", fit4c_agamma6);
			status = m_tuple4->addItem("pgamma6", fit4c_pgamma6);
			status = m_tuple4->addItem("egamma6", fit4c_egamma6);
			status = m_tuple4->addItem("pxgamma6", fit4c_pxgamma6);
			status = m_tuple4->addItem("pygamma6", fit4c_pygamma6);
			status = m_tuple4->addItem("pzgamma6", fit4c_pzgamma6);
			// pi01
			status = m_tuple4->addItem("mpi01", fit4c_mpi01);
			status = m_tuple4->addItem("api01", fit4c_api01);
			status = m_tuple4->addItem("ppi01", fit4c_ppi01);
			status = m_tuple4->addItem("epi01", fit4c_epi01);
			status = m_tuple4->addItem("pxpi01", fit4c_pxpi01);
			status = m_tuple4->addItem("pypi01", fit4c_pypi01);
			status = m_tuple4->addItem("pzpi01", fit4c_pzpi01);
			// pi02
			status = m_tuple4->addItem("mpi02", fit4c_mpi02);
			status = m_tuple4->addItem("api02", fit4c_api02);
			status = m_tuple4->addItem("ppi02", fit4c_ppi02);
			status = m_tuple4->addItem("epi02", fit4c_epi02);
			status = m_tuple4->addItem("pxpi02", fit4c_pxpi02);
			status = m_tuple4->addItem("pypi02", fit4c_pypi02);
			status = m_tuple4->addItem("pzpi02", fit4c_pzpi02);
			// pi03
			status = m_tuple4->addItem("mpi03", fit4c_mpi03);
			status = m_tuple4->addItem("api03", fit4c_api03);
			status = m_tuple4->addItem("ppi03", fit4c_ppi03);
			status = m_tuple4->addItem("epi03", fit4c_epi03);
			status = m_tuple4->addItem("pxpi03", fit4c_pxpi03);
			status = m_tuple4->addItem("pypi03", fit4c_pypi03);
			status = m_tuple4->addItem("pzpi03", fit4c_pzpi03);
			// omega
			status = m_tuple4->addItem("momega", fit4c_momega);
			status = m_tuple4->addItem("aomega", fit4c_aomega);
			status = m_tuple4->addItem("pomega", fit4c_pomega);
			status = m_tuple4->addItem("eomega", fit4c_eomega);
			status = m_tuple4->addItem("pxomega", fit4c_pxomega);
			status = m_tuple4->addItem("pyomega", fit4c_pyomega);
			status = m_tuple4->addItem("pzomega", fit4c_pzomega);
			// omegapi02
			status = m_tuple4->addItem("momegapi02", fit4c_momegapi02);
			status = m_tuple4->addItem("aomegapi02", fit4c_aomegapi02);
			status = m_tuple4->addItem("pomegapi02", fit4c_pomegapi02);
			// omegapi03
			status = m_tuple4->addItem("momegapi03", fit4c_momegapi03);
			status = m_tuple4->addItem("aomegapi03", fit4c_aomegapi03);
			status = m_tuple4->addItem("pomegapi03", fit4c_pomegapi03);
			// pi02pi03
			status = m_tuple4->addItem("mpi02pi03", fit4c_mpi02pi03);
			status = m_tuple4->addItem("api02pi03", fit4c_api02pi03);
			status = m_tuple4->addItem("ppi02pi03", fit4c_ppi02pi03);
			// boost
			status = m_tuple4->addItem("hgamma1", fit4c_hgamma1);
			status = m_tuple4->addItem("hgamma2", fit4c_hgamma2);
			status = m_tuple4->addItem("hgamma3", fit4c_hgamma3);
			status = m_tuple4->addItem("hgamma4", fit4c_hgamma4);
			status = m_tuple4->addItem("hgamma5", fit4c_hgamma5);
			status = m_tuple4->addItem("hgamma6", fit4c_hgamma6);
		}
		else
		{
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple4) << endmsg;
			return StatusCode::FAILURE;
		}
	}
#pragma endregion
	log << MSG::INFO << "successfully return from initialize()" << endmsg;
	return StatusCode::SUCCESS;
}
#pragma endregion
#pragma region 循环执行事例
StatusCode Omega::execute() //
{
#pragma region section_初始化
	MsgStream log(msgSvc(), name());														   //
	log << MSG::INFO << "in execute()" << endreq;											   //
	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");			   //
	int runNo = eventHeader->runNumber();													   // 读取runNo：runnumber
	int event = eventHeader->eventNumber();													   // 读取event：eventnumber
	runID = runNo;																			   // 变量：topo
	eventID = event;																		   //
	Ncut0++;																				   // Ncut0
	log << MSG::DEBUG << "run, evtnum = "													   //
		<< runNo << " , "																	   //
		<< event << endreq;																	   //
	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);		   //
	log << MSG::DEBUG << "ncharg, nneu, tottks = "											   //
		<< evtRecEvent->totalCharged() << " , "												   //
		<< evtRecEvent->totalNeutral() << " , "												   //
		<< evtRecEvent->totalTracks() << endreq;											   //
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol); //
#pragma endregion
#pragma region section_topoloty
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
						m_motheridx[m_numParticle] = ((*iter_mc2)->mother()).trackIndex() + 1;	 // 变量：topo
					}																			 //
					m_numParticle += 1;															 //
					m_pdgid[0] = 11111;															 // 变量：topo
					m_motheridx[0] = 0;															 // 变量：topo
				}																				 //
				m_idxmc = m_numParticle;														 // 变量：topo
			}																					 //
		}																						 //
	}																							 //
#pragma endregion
#pragma region section_truth
	int truth_check = 0;
	if (eventHeader->runNumber() < 0)
	{
		// 设定变量
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
		// 设定直接轨迹变量 - track
		HepLorentzVector track_isr, track_pip, track_pim;
		HepLorentzVector track_gamma1, track_gamma2, track_gamma3, track_gamma4, track_gamma5, track_gamma6;
		HepLorentzVector track_pi01, track_pi02, track_pi03, track_omega;
		// 设定直接轨迹变量 - number
		int n_isr = 0, n_pip = 0, n_pim = 0;
		int n_gamma1 = 0, n_gamma2 = 0, n_gamma3 = 0, n_gamma4 = 0, n_gamma5 = 0, n_gamma6 = 0;
		int n_pi01 = 0, n_pi02 = 0, n_pi03 = 0, n_omega = 0;
		// 设定直接轨迹变量 - index
		int index_isr, index_pip, index_pim;
		int index_gamma1, index_gamma2, index_gamma3, index_gamma4, index_gamma5, index_gamma6;
		int index_pi01, index_pi02, index_pi03, index_omega;
		// 设定间接轨迹变量
		HepLorentzVector track_omegapi02;
		HepLorentzVector track_omegapi03;
		HepLorentzVector track_pi02pi03;
		// 设定暂存交换变量
		HepLorentzVector track_medium;
		int index_medium;
		// 开始筛选
		Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
		// 统计 omega, pi+, pi-, pi0
		for (; iter_mc != mcParticleCol->end(); iter_mc++)
		{
			if (!(*iter_mc)->decayFromGenerator())
				continue;
			HepLorentzVector mctrue_track = (*iter_mc)->initialFourMomentum();
			int mctrue_index = (*iter_mc)->trackIndex();
			// 统计omega,标记233
			if ((*iter_mc)->particleProperty() == 223)
			{
				track_omega = mctrue_track;
				index_omega = mctrue_index;
				n_omega += 1;
			}
			// 统计pi+,标记211,且来自223
			if ((*iter_mc)->particleProperty() == 211 && ((*iter_mc)->mother()).particleProperty() == 223)
			{
				track_pip = mctrue_track;
				index_pip = mctrue_index;
				n_pip += 1;
			}
			// 统计pi-,标记-211,且来自223
			if ((*iter_mc)->particleProperty() == -211 && ((*iter_mc)->mother()).particleProperty() == 223)
			{
				track_pim = mctrue_track;
				index_pim = mctrue_index;
				n_pim += 1;
			}
			// 统计pi0,标记111,来自223定为pi01,否则定为pi02,pi03
			if ((*iter_mc)->particleProperty() == 111)
			{
				if (((*iter_mc)->mother()).particleProperty() == 223)
				{
					track_pi01 = mctrue_track;
					index_pi01 = mctrue_index;
					n_pi01 += 1;
				}
				else if (n_pi02 == 0)
				{
					track_pi02 = mctrue_track;
					index_pi02 = mctrue_index;
					n_pi02 += 1;
				}
				else
				{
					track_pi03 = mctrue_track;
					index_pi03 = mctrue_index;
					n_pi03 += 1;
				}
			}
			// 统计gamma,标记22
			if ((*iter_mc)->particleProperty() == 22 && ((*iter_mc)->mother()).particleProperty() != 111)
			{
				track_isr = mctrue_track;
				index_isr = mctrue_index;
				n_isr += 1;
			}
		}
		// 定位 pi02 pi03
		if (1 == 1)
		{
			track_omegapi02 = track_omega + track_pi02;
			track_omegapi03 = track_omega + track_pi03;
			if (track_omegapi02.m() > track_omegapi03.m())
			{
				track_medium = track_pi03;
				track_pi03 = track_pi02;
				track_pi02 = track_medium;
				index_medium = index_pi03;
				index_pi03 = index_pi02;
				index_pi02 = index_medium;
			}
			track_omegapi02 = track_omega + track_pi02;
			track_omegapi03 = track_omega + track_pi03;
			track_pi02pi03 = track_pi02 + track_pi03;
		}
		// 再次开始筛选
		iter_mc = mcParticleCol->begin();
		for (; iter_mc != mcParticleCol->end(); iter_mc++)
		{
			if (!(*iter_mc)->decayFromGenerator())
				continue;
			HepLorentzVector mctrue_track = (*iter_mc)->initialFourMomentum();
			int mctrue_index = (*iter_mc)->trackIndex();
			if ((*iter_mc)->particleProperty() == 22 && ((*iter_mc)->mother()).particleProperty() == 111)
			{
				if (((*iter_mc)->mother()).trackIndex() == index_pi01)
				{
					if (n_gamma1 == 0)
					{
						track_gamma1 = mctrue_track;
						index_gamma1 = mctrue_index;
						n_gamma1++;
					}
					else
					{
						track_gamma2 = mctrue_track;
						index_gamma2 = mctrue_index;
						n_gamma2++;
					}
				}
				if (((*iter_mc)->mother()).trackIndex() == index_pi02)
				{
					if (n_gamma3 == 0)
					{
						track_gamma3 = mctrue_track;
						index_gamma3 = mctrue_index;
						n_gamma3++;
					}
					else
					{
						track_gamma4 = mctrue_track;
						index_gamma4 = mctrue_index;
						n_gamma4++;
					}
				}
				if (((*iter_mc)->mother()).trackIndex() == index_pi03)
				{
					if (n_gamma5 == 0)
					{
						track_gamma5 = mctrue_track;
						index_gamma5 = mctrue_index;
						n_gamma5++;
					}
					else
					{
						track_gamma6 = mctrue_track;
						index_gamma6 = mctrue_index;
						n_gamma6++;
					}
				}
			}
		}
		// 判断
		if (n_pip == 1 && n_pim == 1 && n_pi01 == 1 && n_pi02 == 1 && n_pi03 == 1 && n_omega == 1)
		{
			if (n_gamma1 == 1 && n_gamma2 == 1 && n_gamma3 == 1 && n_gamma4 == 1 && n_gamma5 == 1 && n_gamma6 == 1)
			{
				truth_check = 1;
			}
		}
		HepLorentzVector track_bgamma1 = track_gamma1.boost(-track_pi01.boostVector());
		HepLorentzVector track_bgamma2 = track_gamma2.boost(-track_pi01.boostVector());
		HepLorentzVector track_bgamma3 = track_gamma3.boost(-track_pi02.boostVector());
		HepLorentzVector track_bgamma4 = track_gamma4.boost(-track_pi02.boostVector());
		HepLorentzVector track_bgamma5 = track_gamma5.boost(-track_pi03.boostVector());
		HepLorentzVector track_bgamma6 = track_gamma6.boost(-track_pi03.boostVector());
		double px1 = track_bgamma1.px();
		double py1 = track_bgamma1.py();
		double pz1 = track_bgamma1.pz();
		double px2 = track_bgamma2.px();
		double py2 = track_bgamma2.py();
		double pz2 = track_bgamma2.pz();
		double px3 = track_bgamma3.px();
		double py3 = track_bgamma3.py();
		double pz3 = track_bgamma3.pz();
		double px4 = track_bgamma4.px();
		double py4 = track_bgamma4.py();
		double pz4 = track_bgamma4.pz();
		double px5 = track_bgamma5.px();
		double py5 = track_bgamma5.py();
		double pz5 = track_bgamma5.pz();
		double px6 = track_bgamma6.px();
		double py6 = track_bgamma6.py();
		double pz6 = track_bgamma6.pz();
		double px01 = track_pi01.px();
		double py01 = track_pi01.py();
		double pz01 = track_pi01.pz();
		double px02 = track_pi02.px();
		double py02 = track_pi02.py();
		double pz02 = track_pi02.pz();
		double px03 = track_pi03.px();
		double py03 = track_pi03.py();
		double pz03 = track_pi03.pz();
		// 填入信息
		if (truth_check == 1)
		{
			// 输出gamma信息
			if (1 == 1)
			{
				//isr
				truth_misr = track_isr.m();
				truth_aisr = track_isr.cosTheta();
				truth_pisr = track_isr.rho();
				truth_eisr = track_isr.e();
				truth_pxisr = track_isr.px();
				truth_pyisr = track_isr.py();
				truth_pzisr = track_isr.pz();
				// pip
				truth_mpip = track_pip.m();
				truth_apip = track_pip.cosTheta();
				truth_ppip = track_pip.rho();
				truth_epip = track_pip.e();
				truth_pxpip = track_pip.px();
				truth_pypip = track_pip.py();
				truth_pzpip = track_pip.pz();
				// pim
				truth_mpim = track_pim.m();
				truth_apim = track_pim.cosTheta();
				truth_ppim = track_pim.rho();
				truth_epim = track_pim.e();
				truth_pxpim = track_pim.px();
				truth_pypim = track_pim.py();
				truth_pzpim = track_pim.pz();
				// gamma1
				truth_mgamma1 = track_gamma1.m();
				truth_agamma1 = track_gamma1.cosTheta();
				truth_pgamma1 = track_gamma1.rho();
				truth_egamma1 = track_gamma1.e();
				truth_pxgamma1 = track_gamma1.px();
				truth_pygamma1 = track_gamma1.py();
				truth_pzgamma1 = track_gamma1.pz();
				// gamma2
				truth_mgamma2 = track_gamma2.m();
				truth_agamma2 = track_gamma2.cosTheta();
				truth_pgamma2 = track_gamma2.rho();
				truth_egamma2 = track_gamma2.e();
				truth_pxgamma2 = track_gamma2.px();
				truth_pygamma2 = track_gamma2.py();
				truth_pzgamma2 = track_gamma2.pz();
				// gamma3
				truth_mgamma3 = track_gamma3.m();
				truth_agamma3 = track_gamma3.cosTheta();
				truth_pgamma3 = track_gamma3.rho();
				truth_egamma3 = track_gamma3.e();
				truth_pxgamma3 = track_gamma3.px();
				truth_pygamma3 = track_gamma3.py();
				truth_pzgamma3 = track_gamma3.pz();
				// gamma4
				truth_mgamma4 = track_gamma4.m();
				truth_agamma4 = track_gamma4.cosTheta();
				truth_pgamma4 = track_gamma4.rho();
				truth_egamma4 = track_gamma4.e();
				truth_pxgamma4 = track_gamma4.px();
				truth_pygamma4 = track_gamma4.py();
				truth_pzgamma4 = track_gamma4.pz();
				// gamma5
				truth_mgamma5 = track_gamma5.m();
				truth_agamma5 = track_gamma5.cosTheta();
				truth_pgamma5 = track_gamma5.rho();
				truth_egamma5 = track_gamma5.e();
				truth_pxgamma5 = track_gamma5.px();
				truth_pygamma5 = track_gamma5.py();
				truth_pzgamma5 = track_gamma5.pz();
				// gamma6
				truth_mgamma6 = track_gamma6.m();
				truth_agamma6 = track_gamma6.cosTheta();
				truth_pgamma6 = track_gamma6.rho();
				truth_egamma6 = track_gamma6.e();
				truth_pxgamma6 = track_gamma6.px();
				truth_pygamma6 = track_gamma6.py();
				truth_pzgamma6 = track_gamma6.pz();
				// pi01
				truth_mpi01 = track_pi01.m();
				truth_api01 = track_pi01.cosTheta();
				truth_ppi01 = track_pi01.rho();
				truth_epi01 = track_pi01.e();
				truth_pxpi01 = track_pi01.px();
				truth_pypi01 = track_pi01.py();
				truth_pzpi01 = track_pi01.pz();
				// pi02
				truth_mpi02 = track_pi02.m();
				truth_api02 = track_pi02.cosTheta();
				truth_ppi02 = track_pi02.rho();
				truth_epi02 = track_pi02.e();
				truth_pxpi02 = track_pi02.px();
				truth_pypi02 = track_pi02.py();
				truth_pzpi02 = track_pi02.pz();
				// pi03
				truth_mpi03 = track_pi03.m();
				truth_api03 = track_pi03.cosTheta();
				truth_ppi03 = track_pi03.rho();
				truth_epi03 = track_pi03.e();
				truth_pxpi03 = track_pi03.px();
				truth_pypi03 = track_pi03.py();
				truth_pzpi03 = track_pi03.pz();
				// omega
				truth_momega = track_omega.m();
				truth_aomega = track_omega.cosTheta();
				truth_pomega = track_omega.rho();
				truth_eomega = track_omega.e();
				truth_pxomega = track_omega.px();
				truth_pyomega = track_omega.py();
				truth_pzomega = track_omega.pz();
				// omegapi02
				truth_momegapi02 = track_omegapi02.m();
				truth_aomegapi02 = track_omegapi02.cosTheta();
				truth_pomegapi02 = track_omegapi02.rho();
				// omegapi03
				truth_momegapi03 = track_omegapi03.m();
				truth_aomegapi03 = track_omegapi03.cosTheta();
				truth_pomegapi03 = track_omegapi03.rho();
				// pi02pi03
				truth_mpi02pi03 = track_pi02pi03.m();
				truth_api02pi03 = track_pi02pi03.cosTheta();
				truth_ppi02pi03 = track_pi02pi03.rho();
				// boost gamma
				truth_hgamma1 = my_helicityangle(track_gamma1, track_pi01);
				truth_hgamma1 = my_helicityangle(track_gamma2, track_pi01);
				truth_hgamma1 = my_helicityangle(track_gamma3, track_pi02);
				truth_hgamma1 = my_helicityangle(track_gamma4, track_pi02);
				truth_hgamma1 = my_helicityangle(track_gamma5, track_pi03);
				truth_hgamma1 = my_helicityangle(track_gamma6, track_pi03);
				m_tuple1->write();
				Ncut1 += 1;
			}
		}
	}
#pragma endregion
#pragma region section_runnunber筛选
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
#pragma endregion
#pragma region section_runnumber统计
	if (firstrun == 0)						   // 第一个事例进行初始化
	{										   //
		number = 0;							   //
		checkexit = 0;						   //
		checki = 0;							   //
		SeriesRun.clear();					   //
		SeriesNum.clear();					   //
	}										   //
	firstrun = 1;							   //
	checkexit = 0;							   // run号在SeriesRun中，则为1，在第几位放入checki
	for (int i = 0; i < SeriesRun.size(); i++) // run号不在SeriesRun中，则为0
	{										   //
		if (runNo == SeriesRun[i])			   //
		{									   //
			checkexit = 1;					   //
			checki = i;						   //
		}									   //
	}										   //
	if (checkexit == 1)						   // 若为1，则统计加一
	{										   //
		SeriesNum[checki] += 1;				   //
	}										   //
	else									   //
	{										   //
		SeriesRun.push_back(runNo);			   // 若为0，则加一统计变量
		SeriesNum.push_back(1);				   //
	}										   //
#pragma endregion
#pragma region section_charged track
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
		double phi0 = mdcTrk->helix(1);												   // 表示螺旋线参数
		double xv = xorigin.x();													   // 0: d0 -> 螺旋线到对撞顶点的最小距离
		double yv = xorigin.y();													   // 1: phi0 -> 最小距离的xy平面相角
		double Rxy = (x0 - xv) * cos(phi0) + (y0 - yv) * sin(phi0);					   // 2: kappa
		HepVector a = mdcTrk->helix();												   // 3: d
		HepSymMatrix Ea = mdcTrk->err();											   // 4: tan(lamda)
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
#pragma endregion
#pragma region section_neutral track
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
		if ((emcTrk->module() == 1) && eraw < 0.025)													  // 选择：E-barrel
			continue;																					  //
		if ((emcTrk->module() == 0 || emcTrk->module() == 2) && eraw < 0.050)							  // 选择：E-barrel
			continue;																					  //
		if (emcTrk->time() < 0 || emcTrk->time() > 14)													  // 选择：TDC
			continue;																					  //
		//if (fabs(dang) < 10)																			  // 选择：dang
		//	continue;																					  //
		iGam.push_back(i);			// 变量：iGam[]（参数为good-track序号，内容为track编号）
	}								//
	int nGam = iGam.size();			// 变量：nGam（中性track数量）
	if (nGam < 6 || nGam > 55)		// 选择：nGam
	{								//
		return StatusCode::SUCCESS; //
	}								//
#pragma endregion
#pragma region section_neutral track momentum
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
#pragma endregion
#pragma region section_charged track momentum
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
	Ncut2++;																			   // Ncut2
#pragma endregion
#pragma region section_vertex fit
	RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin() + ipip[0]))->mdcKalTrack(); // Default is pion, for other particles:
	RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin() + ipim[0]))->mdcKalTrack(); // wvppTrk = WTrackParameter(mp, pipTrk->getZHelixP(), pipTrk->getZErrorP()); proton
	WTrackParameter wvpipTrk, wvpimTrk;											  // wvmupTrk = WTrackParameter(mmu, pipTrk->getZHelixMu(), pipTrk->getZErrorMu()); muon
	wvpipTrk = WTrackParameter(mpi, pipTrk->getZHelix(), pipTrk->getZError());	  // wvepTrk = WTrackParameter(me, pipTrk->getZHelixE(), pipTrk->getZErrorE()); electron
	wvpimTrk = WTrackParameter(mpi, pimTrk->getZHelix(), pimTrk->getZError());	  // wvkpTrk = WTrackParameter(mk, pipTrk->getZHelixK(), pipTrk->getZErrorK()); kaon
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
#pragma endregion
#pragma region section_four constrain
#pragma region fourc_循环数组定义
	int selectb[2][2] = {{1, 2},  // selectb
						 {2, 1}}; // 12顺序遍历
#pragma endregion
#pragma region fourc_初始参数定义
	WTrackParameter wpip = vtxfit->wtrk(0);							 //
	WTrackParameter wpim = vtxfit->wtrk(1);							 //
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();		 //
	HepLorentzVector ecms(0.034 * m_energy / 3.097, 0, 0, m_energy); //
	HepLorentzVector ptrackp;										 // ptrack? 代表最初拟合得到的track
	HepLorentzVector ptrackm;										 // n?ptrack? 代表经过?次拟合得到的track
	HepLorentzVector ptrack1;										 //
	HepLorentzVector ptrack2;										 //
	HepLorentzVector ptrack3;										 //
	HepLorentzVector ptrack4;										 //
	HepLorentzVector ptrack5;										 //
	HepLorentzVector ptrack6;										 //
	double chisq_4c_5g = 9999;										 //
	double chisq_4c_6g = 9999;										 //
	double chisq_4c_7g = 9999;										 //
#pragma endregion
#pragma region fourc_6gamma拟合
	for (int i1 = 0; i1 < nGam; i1++)																  // 得到6Gamma-chisq
	{																								  //
		RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();					  //
		for (int i2 = i1; i2 < nGam; i2++)															  //
		{																							  //
			if (i2 == i1)																			  //
			{																						  //
				continue;																			  //
			}																						  //
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();				  //
			for (int i3 = i2; i3 < nGam; i3++)														  //
			{																						  //
				if (i3 == i1 || i3 == i2)															  //
				{																					  //
					continue;																		  //
				}																					  //
				RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[i3]))->emcShower();			  //
				for (int i4 = i3; i4 < nGam; i4++)													  //
				{																					  //
					if (i4 == i1 || i4 == i2 || i4 == i3)											  //
					{																				  //
						continue;																	  //
					}																				  //
					RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[i4]))->emcShower();		  //
					for (int i5 = i4; i5 < nGam; i5++)												  //
					{																				  //
						if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4)							  //
						{																			  //
							continue;																  //
						}																			  //
						RecEmcShower *g5Trk = (*(evtRecTrkCol->begin() + iGam[i5]))->emcShower();	  //
						for (int i6 = i5; i6 < nGam; i6++)											  //
						{																			  //
							if (i6 == i1 || i6 == i2 || i6 == i3 || i6 == i4 || i6 == i5)			  //
							{																		  //
								continue;															  //
							}																		  //
							RecEmcShower *g6Trk = (*(evtRecTrkCol->begin() + iGam[i6]))->emcShower(); //
							kmfit->init();															  //
							kmfit->AddTrack(0, wpip);												  //
							kmfit->AddTrack(1, wpim);												  //
							kmfit->AddTrack(2, 0.0, g1Trk);											  //
							kmfit->AddTrack(3, 0.0, g2Trk);											  //
							kmfit->AddTrack(4, 0.0, g3Trk);											  //
							kmfit->AddTrack(5, 0.0, g4Trk);											  //
							kmfit->AddTrack(6, 0.0, g5Trk);											  //
							kmfit->AddTrack(7, 0.0, g6Trk);											  //
							kmfit->AddFourMomentum(0, ecms);										  //
							bool oksq = kmfit->Fit();												  //
							if (oksq)																  //
							{																		  //
								double chi2 = kmfit->chisq();										  //
								if (chi2 <= chisq_4c_6g)											  // 选择：最小chi-4c
								{																	  //
									chisq_4c_6g = chi2;												  //
									ptrackp = kmfit->pfit(0);										  //
									ptrackm = kmfit->pfit(1);										  //
									ptrack1 = kmfit->pfit(2);										  //
									ptrack2 = kmfit->pfit(3);										  //
									ptrack3 = kmfit->pfit(4);										  //
									ptrack4 = kmfit->pfit(5);										  //
									ptrack5 = kmfit->pfit(6);										  //
									ptrack6 = kmfit->pfit(7);										  //
								}																	  //
							}																		  //
						}																			  //
					}																				  //
				}																					  //
			}																						  //
		}																							  //
	}																								  //
#pragma endregion
#pragma region fourc_5gamma拟合
	for (int i1 = 0; i1 < nGam; i1++)															  // 得到5Gamma-chisq
	{																							  //
		RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();				  //
		for (int i2 = i1; i2 < nGam; i2++)														  //
		{																						  //
			if (i2 == i1)																		  //
			{																					  //
				continue;																		  //
			}																					  //
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();			  //
			for (int i3 = i2; i3 < nGam; i3++)													  //
			{																					  //
				if (i3 == i1 || i3 == i2)														  //
				{																				  //
					continue;																	  //
				}																				  //
				RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[i3]))->emcShower();		  //
				for (int i4 = i3; i4 < nGam; i4++)												  //
				{																				  //
					if (i4 == i1 || i4 == i2 || i4 == i3)										  //
					{																			  //
						continue;																  //
					}																			  //
					RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[i4]))->emcShower();	  //
					for (int i5 = i4; i5 < nGam; i5++)											  //
					{																			  //
						if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4)						  //
						{																		  //
							continue;															  //
						}																		  //
						RecEmcShower *g5Trk = (*(evtRecTrkCol->begin() + iGam[i5]))->emcShower(); //
						kmfit->init();															  //
						kmfit->AddTrack(0, wpip);												  //
						kmfit->AddTrack(1, wpim);												  //
						kmfit->AddTrack(2, 0.0, g1Trk);											  //
						kmfit->AddTrack(3, 0.0, g2Trk);											  //
						kmfit->AddTrack(4, 0.0, g3Trk);											  //
						kmfit->AddTrack(5, 0.0, g4Trk);											  //
						kmfit->AddTrack(6, 0.0, g5Trk);											  //
						kmfit->AddFourMomentum(0, ecms);										  //
						bool oksq = kmfit->Fit();												  //
						if (oksq)																  //
						{																		  //
							double chi2 = kmfit->chisq();										  //
							if (chi2 < chisq_4c_5g)												  //
							{																	  //
								chisq_4c_5g = chi2;												  //
							}																	  //
						}																		  //
					}																			  //
				}																				  //
			}																					  //
		}																						  //
	}																							  //
#pragma endregion
#pragma region fourc_7gamma拟合
	if (nGam > 6)																							  // 得到7Gamma-chisq
	{																										  //
		for (int i1 = 0; i1 < nGam; i1++)																	  //
		{																									  //
			RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();						  //
			for (int i2 = i1; i2 < nGam; i2++)																  //
			{																								  //
				if (i2 == i1)																				  //
				{																							  //
					continue;																				  //
				}																							  //
				RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();					  //
				for (int i3 = i2; i3 < nGam; i3++)															  //
				{																							  //
					if (i3 == i1 || i3 == i2)																  //
					{																						  //
						continue;																			  //
					}																						  //
					RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[i3]))->emcShower();				  //
					for (int i4 = i3; i4 < nGam; i4++)														  //
					{																						  //
						if (i4 == i1 || i4 == i2 || i4 == i3)												  //
						{																					  //
							continue;																		  //
						}																					  //
						RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[i4]))->emcShower();			  //
						for (int i5 = i4; i5 < nGam; i5++)													  //
						{																					  //
							if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4)								  //
							{																				  //
								continue;																	  //
							}																				  //
							RecEmcShower *g5Trk = (*(evtRecTrkCol->begin() + iGam[i5]))->emcShower();		  //
							for (int i6 = i5; i6 < nGam; i6++)												  //
							{																				  //
								if (i6 == i1 || i6 == i2 || i6 == i3 || i6 == i4 || i6 == i5)				  //
								{																			  //
									continue;																  //
								}																			  //
								RecEmcShower *g6Trk = (*(evtRecTrkCol->begin() + iGam[i6]))->emcShower();	  //
								for (int i7 = i6; i7 < nGam; i7++)											  //
								{																			  //
									if (i7 == i1 || i7 == i2 || i7 == i3 || i7 == i4 || i7 == i5 || i7 == i6) //
									{																		  //
										continue;															  //
									}																		  //
									RecEmcShower *g7Trk = (*(evtRecTrkCol->begin() + iGam[i7]))->emcShower(); //
									kmfit->init();															  //
									kmfit->AddTrack(0, wpip);												  //
									kmfit->AddTrack(1, wpim);												  //
									kmfit->AddTrack(2, 0.0, g1Trk);											  //
									kmfit->AddTrack(3, 0.0, g2Trk);											  //
									kmfit->AddTrack(4, 0.0, g3Trk);											  //
									kmfit->AddTrack(5, 0.0, g4Trk);											  //
									kmfit->AddTrack(6, 0.0, g5Trk);											  //
									kmfit->AddTrack(7, 0.0, g6Trk);											  //
									kmfit->AddTrack(8, 0.0, g7Trk);											  //
									kmfit->AddFourMomentum(0, ecms);										  //
									bool oksq = kmfit->Fit();												  //
									if (oksq)																  //
									{																		  //
										double chi2 = kmfit->chisq();										  //
										if (chi2 < chisq_4c_7g)												  //
										{																	  //
											chisq_4c_7g = chi2;												  //
										}																	  //
									}																		  //
								}																			  //
							}																				  //
						}																					  //
					}																						  //
				}																							  //
			}																								  //
		}																									  //
	}																										  //
#pragma endregion
#pragma region fourc_判断6gamma是否成功
	int g6 = 0;							   //
	if (chisq_4c_6g < 200)				   // 选择：chisq-6Gamma
	{									   //
		g6 = 1;							   //
		Ncut3++;						   // Ncut3
	}									   //
	if (1 == 1)							   // 是否比较6gamma与57gamma
	{									   //
		if (chisq_4c_5g < chisq_4c_6g)	   // 选择：chisq-5Gamma
		{								   //
			g6 = 0;						   //
		}								   //
		if (nGam > 6)					   // 选择：chisq-7Gamma
		{								   //
			if (chisq_4c_7g < chisq_4c_6g) //
			{							   //
				g6 = 0;					   //
			}							   //
		}								   //
	}									   //
#pragma endregion
#pragma region fourc_粒子重建
	HepLorentzVector *tracklist;
	tracklist = my_newheplorentzvector(8);
	tracklist[0] = ptrackp;
	tracklist[1] = ptrackm;
	tracklist[2] = ptrack1;
	tracklist[3] = ptrack2;
	tracklist[4] = ptrack3;
	tracklist[5] = ptrack4;
	tracklist[6] = ptrack5;
	tracklist[7] = ptrack6;
	if (g6 == 1)
	{
		if (1 == 1)
		{
			tracklist = my_recon_3pi(tracklist, 0.135, 0.135, 0.135);
		}
		if (1 == 1)
		{
			tracklist = my_recon_omega(tracklist, 0.782);
		}
		if (1 == 1)
		{
			tracklist = my_recon_lower(tracklist);
		}
#pragma endregion
#pragma region fourc_信息输出
		// gamma 信息
		HepLorentzVector out_pip = tracklist[0];
		HepLorentzVector out_pim = tracklist[1];
		HepLorentzVector out_gamma1 = tracklist[2];
		HepLorentzVector out_gamma2 = tracklist[3];
		HepLorentzVector out_gamma3 = tracklist[4];
		HepLorentzVector out_gamma4 = tracklist[5];
		HepLorentzVector out_gamma5 = tracklist[6];
		HepLorentzVector out_gamma6 = tracklist[7];
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
			fit4c_chisq_4c = chisq_4c_6g;
			// pip
			fit4c_mpip = out_pip.m();
			fit4c_apip = out_pip.cosTheta();
			fit4c_ppip = out_pip.rho();
			fit4c_epip = out_pip.e();
			fit4c_pxpip = out_pip.px();
			fit4c_pypip = out_pip.py();
			fit4c_pzpip = out_pip.pz();
			// pim
			fit4c_mpim = out_pim.m();
			fit4c_apim = out_pim.cosTheta();
			fit4c_ppim = out_pim.rho();
			fit4c_epim = out_pim.e();
			fit4c_pxpim = out_pim.px();
			fit4c_pypim = out_pim.py();
			fit4c_pzpim = out_pim.pz();
			// gamma1
			fit4c_mgamma1 = out_gamma1.m();
			fit4c_agamma1 = out_gamma1.cosTheta();
			fit4c_pgamma1 = out_gamma1.rho();
			fit4c_egamma1 = out_gamma1.e();
			fit4c_pxgamma1 = out_gamma1.px();
			fit4c_pygamma1 = out_gamma1.py();
			fit4c_pzgamma1 = out_gamma1.pz();
			// gamma2
			fit4c_mgamma2 = out_gamma2.m();
			fit4c_agamma2 = out_gamma2.cosTheta();
			fit4c_pgamma2 = out_gamma2.rho();
			fit4c_egamma2 = out_gamma2.e();
			fit4c_pxgamma2 = out_gamma2.px();
			fit4c_pygamma2 = out_gamma2.py();
			fit4c_pzgamma2 = out_gamma2.pz();
			// gamma3
			fit4c_mgamma3 = out_gamma3.m();
			fit4c_agamma3 = out_gamma3.cosTheta();
			fit4c_pgamma3 = out_gamma3.rho();
			fit4c_egamma3 = out_gamma3.e();
			fit4c_pxgamma3 = out_gamma3.px();
			fit4c_pygamma3 = out_gamma3.py();
			fit4c_pzgamma3 = out_gamma3.pz();
			// gamma4
			fit4c_mgamma4 = out_gamma4.m();
			fit4c_agamma4 = out_gamma4.cosTheta();
			fit4c_pgamma4 = out_gamma4.rho();
			fit4c_egamma4 = out_gamma4.e();
			fit4c_pxgamma4 = out_gamma4.px();
			fit4c_pygamma4 = out_gamma4.py();
			fit4c_pzgamma4 = out_gamma4.pz();
			// gamma5
			fit4c_mgamma5 = out_gamma5.m();
			fit4c_agamma5 = out_gamma5.cosTheta();
			fit4c_pgamma5 = out_gamma5.rho();
			fit4c_egamma5 = out_gamma5.e();
			fit4c_pxgamma5 = out_gamma5.px();
			fit4c_pygamma5 = out_gamma5.py();
			fit4c_pzgamma5 = out_gamma5.pz();
			// gamma6
			fit4c_mgamma6 = out_gamma6.m();
			fit4c_agamma6 = out_gamma6.cosTheta();
			fit4c_pgamma6 = out_gamma6.rho();
			fit4c_egamma6 = out_gamma6.e();
			fit4c_pxgamma6 = out_gamma6.px();
			fit4c_pygamma6 = out_gamma6.py();
			fit4c_pzgamma6 = out_gamma6.pz();
			// pi01
			fit4c_mpi01 = out_pi01.m();
			fit4c_api01 = out_pi01.cosTheta();
			fit4c_ppi01 = out_pi01.rho();
			fit4c_epi01 = out_pi01.e();
			fit4c_pxpi01 = out_pi01.px();
			fit4c_pypi01 = out_pi01.py();
			fit4c_pzpi01 = out_pi01.pz();
			// pi02
			fit4c_mpi02 = out_pi02.m();
			fit4c_api02 = out_pi02.cosTheta();
			fit4c_ppi02 = out_pi02.rho();
			fit4c_epi02 = out_pi02.e();
			fit4c_pxpi02 = out_pi02.px();
			fit4c_pypi02 = out_pi02.py();
			fit4c_pzpi02 = out_pi02.pz();
			// pi03
			fit4c_mpi03 = out_pi03.m();
			fit4c_api03 = out_pi03.cosTheta();
			fit4c_ppi03 = out_pi03.rho();
			fit4c_epi03 = out_pi03.e();
			fit4c_pxpi03 = out_pi03.px();
			fit4c_pypi03 = out_pi03.py();
			fit4c_pzpi03 = out_pi03.pz();
			// omega
			fit4c_momega = out_omega.m();
			fit4c_aomega = out_omega.cosTheta();
			fit4c_pomega = out_omega.rho();
			fit4c_eomega = out_omega.e();
			fit4c_pxomega = out_omega.px();
			fit4c_pyomega = out_omega.py();
			fit4c_pzomega = out_omega.pz();
			// omegapi02
			fit4c_momegapi02 = out_omegapi02.m();
			fit4c_aomegapi02 = out_omegapi02.cosTheta();
			fit4c_pomegapi02 = out_omegapi02.rho();
			// omegapi03
			fit4c_momegapi03 = out_omegapi03.m();
			fit4c_aomegapi03 = out_omegapi03.cosTheta();
			fit4c_pomegapi03 = out_omegapi03.rho();
			// pi02pi03
			fit4c_mpi02pi03 = out_pi02pi03.m();
			fit4c_api02pi03 = out_pi02pi03.cosTheta();
			fit4c_ppi02pi03 = out_pi02pi03.rho();
			// boost gamma
			fit4c_hgamma1 = my_helicityangle(out_gamma1, out_pi01);
			fit4c_hgamma2 = my_helicityangle(out_gamma2, out_pi01);
			fit4c_hgamma3 = my_helicityangle(out_gamma3, out_pi02);
			fit4c_hgamma4 = my_helicityangle(out_gamma4, out_pi02);
			fit4c_hgamma5 = my_helicityangle(out_gamma5, out_pi03);
			fit4c_hgamma6 = my_helicityangle(out_gamma6, out_pi03);
			m_tuple4->write();
			Ncut4++;
		}
#pragma endregion
	}
	my_constant a_constant;
	cout << a_constant.pi << endl;
#pragma endregion
	return StatusCode::SUCCESS;
}
#pragma endregion
#pragma region 结束输出
//*********************************************************************************************************
//***                                               finalize                                            ***
//*********************************************************************************************************
StatusCode Omega::finalize()
{
	cout << "energy point:         " << m_energy << endl;
	cout << "total number:         " << Ncut0 << endl;
	cout << "Pass truth:           " << Ncut1 << endl;
	cout << "Pass Pid:             " << Ncut2 << endl;
	cout << "Pass 4c-6gamma:       " << Ncut3 << endl;
	cout << "Pass 4C:              " << Ncut4 << endl;
	for (int i = 0; i < SeriesRun.size(); i++)
	{
		cout << "oooooooo" << SeriesRun[i] << "oooooooo" << SeriesNum[i] << "oooooooo" << endl;
	}
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}
#pragma endregion
