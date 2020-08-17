#pragma region 1.1.include
#include "McTruth/McParticle.h"
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
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
using CLHEP::Hep2Vector;
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
#include <vector>
#include "PppmpzAlg/Pppmpz.h"
#include "PppmpzAlg/Myfunc.h"
#pragma endregion
#pragma region 1.2.parameter_init
// 定义类型
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
// 定义常数
my_constant use_constant;
const double mpiz = use_constant.mpiz;
const double mpipm = use_constant.mpipm;
int Ncut0, Ncut1, Ncut2, Ncut3, Ncut4, Ncut5;
Vint SeriesRun;
Vint SeriesNum;
Vint SeriesNum1;
Vint SeriesNum2;
Vint SeriesNum3;
Vint SeriesNum4;
Vint SeriesNum5;
int firstrun = 0;
Pppmpz::Pppmpz(const std::string &name, ISvcLocator *pSvcLocator) : Algorithm(name, pSvcLocator)
{
	// Determine parameter
	declareProperty("Energy", job_energy = 0);
	declareProperty("flag3", job_flag3 = 0);
	// Determine process
	declareProperty("do_truth", job_do_truth = 1);
	declareProperty("do_pid", job_do_pid = 0);
	declareProperty("do_4c", job_do_4c = 1);
	declareProperty("do_4c_0", job_do_4c_0 = 1);
	declareProperty("do_4c_1", job_do_4c_1 = 1);
	declareProperty("do_4c_3", job_do_4c_3 = 1);
	declareProperty("do_4c_4", job_do_4c_4 = 1);
}
#pragma endregion
#pragma region 3.tuple
StatusCode Pppmpz::initialize()
{
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in initialize()" << endmsg;
	StatusCode status;
#pragma region buildtruth
	if (job_do_truth)
	{
		NTuplePtr nt1(ntupleSvc(), "FILE1/truth");
		if (nt1)
		{
			m_tuple_truth = nt1;
		}
		else
		{
			m_tuple_truth = ntupleSvc()->book("FILE1/truth", CLID_ColumnWiseTuple, "ks N-Tuple example");
			if (m_tuple_truth)
			{
				status = m_tuple_truth->addItem("mpip", truth_mpip);
				status = m_tuple_truth->addItem("apip", truth_apip);
				status = m_tuple_truth->addItem("ppip", truth_ppip);
				status = m_tuple_truth->addItem("mpim", truth_mpim);
				status = m_tuple_truth->addItem("apim", truth_apim);
				status = m_tuple_truth->addItem("ppim", truth_ppim);
				status = m_tuple_truth->addItem("mpiz", truth_mpiz);
				status = m_tuple_truth->addItem("apiz", truth_apiz);
				status = m_tuple_truth->addItem("ppiz", truth_ppiz);
				status = m_tuple_truth->addItem("mpipz", truth_mpipz);
				status = m_tuple_truth->addItem("apipz", truth_apipz);
				status = m_tuple_truth->addItem("ppipz", truth_ppipz);
				status = m_tuple_truth->addItem("mpimz", truth_mpimz);
				status = m_tuple_truth->addItem("apimz", truth_apimz);
				status = m_tuple_truth->addItem("ppimz", truth_ppimz);
				status = m_tuple_truth->addItem("mpipm", truth_mpipm);
				status = m_tuple_truth->addItem("apipm", truth_apipm);
				status = m_tuple_truth->addItem("ppipm", truth_ppipm);
			}
			else
			{
				log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple_truth) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
#pragma endregion
#pragma region 新建树：charge
	NTuplePtr nt2(ntupleSvc(), "FIL		1/charge");
	if (nt2)
	{
		m_tuple2 = nt2;
	}
	else
	{
		m_tuple2 = ntupleSvc()->book("FILE1/charge", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if (m_tuple2)
		{
			status = m_tuple2->addItem("ngood", charge_ngood);
			status = m_tuple2->addItem("ncharge", charge_ncharge);
		}
		else
		{
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple2) << endmsg;
			return StatusCode::FAILURE;
		}
	}
#pragma endregion
#pragma region 新建树：fit4c
	if (job_do_4c)
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
				// Topology信息
				status = m_tuple4->addItem("runID", runNo);
				status = m_tuple4->addItem("eventID", event);
				status = m_tuple4->addItem("indexmc", m_idxmc, 0, 100);
				status = m_tuple4->addIndexedItem("pdgid", m_idxmc, m_pdgid);
				status = m_tuple4->addIndexedItem("motheridx", m_idxmc, m_motheridx);
				status = m_tuple4->addItem("flag1", flag1);
				status = m_tuple4->addItem("flag2", flag2);
				status = m_tuple4->addItem("flag3", flag3);
				// For Cut
				status = m_tuple4->addItem("chisq", fit4c_chisq);
				status = m_tuple4->addItem("chisq_0g", fit4c_chisq_0g);
				status = m_tuple4->addItem("chisq_1g", fit4c_chisq_1g);
				status = m_tuple4->addItem("chisq_3g", fit4c_chisq_3g);
				status = m_tuple4->addItem("chisq_4g", fit4c_chisq_4g);
				// For track
				status = m_tuple4->addItem("egamma1", fit4c_egamma1);
				status = m_tuple4->addItem("egamma2", fit4c_egamma2);
				status = m_tuple4->addItem("hgamma1", fit4c_hgamma1);
				status = m_tuple4->addItem("hgamma2", fit4c_hgamma2);
				status = m_tuple4->addItem("agamma", fit4c_agamma);
				status = m_tuple4->addItem("mpip", fit4c_mpip);
				status = m_tuple4->addItem("apip", fit4c_apip);
				status = m_tuple4->addItem("ppip", fit4c_ppip);
				status = m_tuple4->addItem("mpim", fit4c_mpim);
				status = m_tuple4->addItem("apim", fit4c_apim);
				status = m_tuple4->addItem("ppim", fit4c_ppim);
				status = m_tuple4->addItem("mpiz", fit4c_mpiz);
				status = m_tuple4->addItem("apiz", fit4c_apiz);
				status = m_tuple4->addItem("ppiz", fit4c_ppiz);
				status = m_tuple4->addItem("mpipz", fit4c_mpipz);
				status = m_tuple4->addItem("apipz", fit4c_apipz);
				status = m_tuple4->addItem("ppipz", fit4c_ppipz);
				status = m_tuple4->addItem("mpimz", fit4c_mpimz);
				status = m_tuple4->addItem("apimz", fit4c_apimz);
				status = m_tuple4->addItem("ppimz", fit4c_ppimz);
				status = m_tuple4->addItem("mpipm", fit4c_mpipm);
				status = m_tuple4->addItem("apipm", fit4c_apipm);
				status = m_tuple4->addItem("ppipm", fit4c_ppipm);
			}
			else
			{
				log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple4) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
#pragma endregion
	log << MSG::INFO << "successfully return from initialize()" << endmsg;
	return StatusCode::SUCCESS;
}
#pragma endregion
#pragma region 循环执行事例
StatusCode Pppmpz::execute() //
{
#pragma region section_初始化
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;
	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
	runNo = eventHeader->runNumber();
	event = eventHeader->eventNumber();
	if (runNo < 0)
	{
		flag1 = eventHeader->flag1();
		flag2 = eventHeader->flag2();
	}
	else
	{
		flag1 = 0;
		flag2 = 0;
	}
	flag3 = job_flag3;
	Ncut0++;
	log << MSG::DEBUG << "run, evtnum = "
		<< runNo << " , "
		<< event << endreq;
	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	log << MSG::DEBUG << "ncharg, nneu, tottks = "
		<< evtRecEvent->totalCharged() << " , "
		<< evtRecEvent->totalNeutral() << " , "
		<< evtRecEvent->totalTracks() << endreq;
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
#pragma endregion
#pragma region section_runnunber_筛选
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
#pragma region section_checkrun_初始化
	my_seriesinit(SeriesRun,
				  SeriesNum,
				  SeriesNum1,
				  SeriesNum2,
				  SeriesNum3,
				  SeriesNum4,
				  SeriesNum5,
				  runNo,
				  firstrun);
#pragma region section_truth
	int truth_check = 0;
	if (eventHeader->runNumber() < 0)
	{														 
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
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
	if (eventHeader->runNumber() < 0 && job_do_truth)
	{
		// 设定变量
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
		// 设定直接轨迹变量 - track
		HepLorentzVector truth_isr;
		HepLorentzVector truth_pip;
		HepLorentzVector truth_pim;
		HepLorentzVector truth_gamma1;
		HepLorentzVector truth_gamma2;
		HepLorentzVector truth_piz;
		// 设定直接轨迹变量 - number
		int n_isr = 0, index_isr;
		int n_pip = 0, index_pip;
		int n_pim = 0, index_pim;
		int n_gamma1 = 0, index_gamma1;
		int n_gamma2 = 0, index_gamma2;
		int n_piz = 0, index_piz;
		// 设定间接轨迹变量
		HepLorentzVector truth_pipm;
		HepLorentzVector truth_pipz;
		HepLorentzVector truth_pimz;
		// 设定暂存交换变量
		HepLorentzVector truth_medium;
		int index_medium;
		// 开始筛选
		Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
		// 统计粒子
		for (; iter_mc != mcParticleCol->end(); iter_mc++)
		{
			if (!(*iter_mc)->decayFromGenerator())
				continue;
			HepLorentzVector mctrue_track = (*iter_mc)->initialFourMomentum();
			int mctrue_index = (*iter_mc)->trackIndex();
			// 统计pi+,标记211
			if ((*iter_mc)->particleProperty() == 211)
			{
				truth_pip = mctrue_track;
				n_pip += 1;
			}
			// 统计pi-,标记-211
			if ((*iter_mc)->particleProperty() == -211)
			{
				truth_pim = mctrue_track;
				n_pim += 1;
			}
			// 统计pi0,标记111
			if ((*iter_mc)->particleProperty() == 111)
			{
				truth_piz = mctrue_track;
				n_piz += 1;
			}
			// 统计gamma,标记22
			if ((*iter_mc)->particleProperty() == 22 && ((*iter_mc)->mother()).particleProperty() != 111)
			{
				truth_isr = mctrue_track;
				n_isr += 1;
			}
			// 统计gamma,标记22
			if ((*iter_mc)->particleProperty() == 22 && ((*iter_mc)->mother()).particleProperty() == 111)
			{
				if (n_gamma1 == 0)
				{
					truth_gamma1 = mctrue_track;
					n_gamma1++;
				}
				else
				{
					truth_gamma2 = mctrue_track;
					n_gamma2++;
				}
			}
		}
		// 填入信息
		if (n_pip == 1 && n_pim == 1 && n_piz == 1 && n_gamma1 == 1 && n_gamma2 == 1)
		{
			// 输出gamma信息
			if (1 == 1)
			{
				truth_mpip = truth_pip.m();
				truth_apip = truth_pip.cosTheta();
				truth_ppip = truth_pip.rho();
				truth_mpim = truth_pim.m();
				truth_apim = truth_pim.cosTheta();
				truth_ppim = truth_pim.rho();
				truth_mpiz = truth_piz.m();
				truth_apiz = truth_piz.cosTheta();
				truth_ppiz = truth_piz.rho();
				truth_mpipm = truth_pipm.m();
				truth_apipm = truth_pipm.cosTheta();
				truth_ppipm = truth_pipm.rho();
				truth_mpipz = truth_pipz.m();
				truth_apipz = truth_pipz.cosTheta();
				truth_ppipz = truth_pipz.rho();
				truth_mpimz = truth_pimz.m();
				truth_apimz = truth_pimz.cosTheta();
				truth_ppimz = truth_pimz.rho();
				m_tuple_truth->write();
				Ncut1 += 1;
				my_seriescount(SeriesRun, SeriesNum1, runNo);
			}
		}
		else
		{
			cout << "Truth error with:" << n_pip << n_pim << n_piz << n_gamma1 << n_gamma2 << endl;
		}
	}
#pragma endregion
#pragma region section_charged track
	Vint iGood, ipip, ipim;															   // ****************************************
	iGood.clear();																	   // 变量：iGood[]（参数为good-track序号，内容为track编号）
	ipip.clear();																	   // 变量：ipip[]（参数为good-track+序号，内容为tra			ck编号）
	ipim.clear();																	   // 变量：ipim[]（参数为good-track-序号，内容为track编号）
	Vp4 ppip, ppim;																	   //
	ppip.clear();																	   // 变量：ppip[]（good-track+的四动量）
	ppim.clear();																	   // 变量：ppim[]（good-track-的四动量）
	int nCharge;																	   // ****************************************
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
		double z0 = mdcTrk->z();													   // ****************************************
		double phi0 = mdcTrk->helix(1);												   // 表示螺旋线参数
		double xv = xorigin.x();													   // 0: d0 -> 螺旋线到对撞顶点的最小距离
		double yv = xorigin.y();													   // 1: phi0 -> 最小距离的xy平面相角
		double Rxy = (x0 - xv) * cos(phi0) + (y0 - yv) * sin(phi0);					   // 2: kappa
		HepVector a = mdcTrk->helix();												   // 3: d
		HepSymMatrix Ea = mdcTrk->err();											   // 4: tan(lamda)
		HepPoint3D point0(0., 0., 0.);												   // ****************************************
		HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);							   //
		VFHelix helixip(point0, a, Ea);												   //
		helixip.pivot(IP);															   //
		HepVector vecipa = helixip.a();												   //
		double Rvxy0 = fabs(vecipa[0]);												   //
		double Rvz0 = vecipa[3];													   //
		double Rvphi0 = vecipa[1];													   // ****************************************
		if (fabs(cos(mdcTrk->theta())) > 0.93)										   // 选择：cos(theta)
			continue;																   //
		if (fabs(Rvz0) >= 10)														   // 选择：Rvz0
			continue;																   //
		if (fabs(Rvxy0) >= 1)														   // 选择：Rvxy0
			continue;																   // ****************************************
		if ((*itTrk)->isEmcShowerValid())											   //
		{																			   //
			RecEmcShower *emcTrk = (*itTrk)->emcShower();							   //
			if ((emcTrk->energy() / mdcTrk->p()) > 0.85)							   //
				continue;															   //
		}																			   //
		iGood.push_back(i);															   //
		nCharge += mdcTrk->charge();												   //
	}																				   //
	int nGood = iGood.size();														   //
	log << MSG::DEBUG << "ngood, totcharge = " << nGood << " , " << nCharge << endreq; // ****************************************
	if (1 == 1)																		   // 输出2：charge
	{																				   // ****************************************
		charge_ngood = nGood;														   //
		charge_ncharge = nCharge;													   //
		m_tuple2->write();															   //
	}																				   //
	if (1 == 1)																		   // ****************************************
	{																				   // 变量：nGood（带电track数目）
		if ((nGood != 2) || (nCharge != 0))											   // 变量：nCharge（带电track总电量）
		{																			   // 选择1：nGood = 2, nCharge = 0
			return StatusCode::SUCCESS;												   // 选择2：nGood >= 2
		}																			   // ****************************************
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
		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());												  //
		double dthe = 200.;																			  //
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
	if (nGam < 2 || nGam > 10)		// 选择：nGam
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
		ptrk.setPy(eraw * sin(the) * sin(phi));							 //
		ptrk.setPz(eraw * cos(the));								 //
		ptrk.setE(eraw);											 //
		pGam.push_back(ptrk);										 //
	}																 //
#pragma endregion
#pragma region section_charged track momentum
	ParticleID *pid = ParticleID::instance();												 //
	int pid_np = 0;																			 //
	int pid_nm = 0;																			 //
	for (int i = 0; i < nGood; i++)															 //
	{																						 //
		int pid_check = 1;																	 //
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];						 // ****************************************
		pid->init();																		 // 对于Likelihood方法
		pid->setMethod(pid->methodProbability());											 // pid->setMethod(pid->methodLikelihood());
		pid->setChiMinCut(4);																 //
		pid->setRecTrack(*itTrk);															 //
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE()); // 选择pid的粒子
		pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyProton());				 // pid->identify(pid->onlyPion());
		pid->calculate();																	 // pid->identify(pid->onlyKaon());
		if (!(pid->IsPidInfoValid()))														 //
			pid_check = 0;																	 // 对于Likelihood方法(0=electron 1=muon 2=pion 3=kaon 4=proton)
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();											 // if(pid->pdf(2) < pid->pdf(3))
		if (pid->probPion() < pid->probKaon())												 // continue;
			pid_check = 0;																	 //
		if (pid->probPion() < pid->probProton())											 //
			pid_check = 0;																	 //
		RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();								 // 对于ParticleID, 用RecMdcKalTrack代替RecMdcTrack
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);									 // PID可以设定为electron, muon, pion, kaon and proton;The default setting is pion
		HepLorentzVector ptrk;																 // ****************************************
		ptrk.setPx(mdcKalTrk->px());														 // 变量：ppip[]（值为pi+动量）
		ptrk.setPy(mdcKalTrk->py());														 // 变量：ppim[]（值为pi-动量）
		ptrk.setPz(mdcKalTrk->pz());														 // 变量：ipip[]（值为pi+编号）
		double p3 = ptrk.mag();																 // 变量：ipim[]（值为pi-编号）
		ptrk.setE(sqrt(p3 * p3 + mpipm * mpipm));											 // ****************************************
		if (mdcKalTrk->charge() > 0)														 //
		{																					 //
			ipip.push_back(iGood[i]);														 //
			ppip.push_back(ptrk);															 //
			if (pid_check)																	 //
			{																				 //
				pid_np += 1;																 //
			}																				 //
		}																					 //
		else																				 //
		{																					 //
			ipim.push_back(iGood[i]);														 //
			ppim.push_back(ptrk);															 //
			if (pid_check)																	 //
			{																				 //
				pid_nm += 1;																 //
			}																				 //
		}																					 //
	}																						 //
	if (pid_np * pid_nm != 1)																 //
	{																						 //
		if (job_do_pid)																			 //
		{																					 //
			return SUCCESS;																	 //
		}																					 //
	}																						 //
	Ncut2++;																				 //
	my_seriescount(SeriesRun, SeriesNum2, runNo);											 //
#pragma endregion
#pragma region section_vertex fit
	RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin() + ipip[0]))->mdcKalTrack(); // ****************************************
	RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin() + ipim[0]))->mdcKalTrack(); // Default is pion, for other particles:
	WTrackParameter wvpipTrk, wvpimTrk;											  // wvppTrk = WTrackParameter(mp, pipTrk->getZHelixP(), pipTrk->getZErrorP()); proton
	wvpipTrk = WTrackParameter(mpipm, pipTrk->getZHelix(), pipTrk->getZError());  // wvepTrk = WTrackParameter(me, pipTrk->getZHelixE(), pipTrk->getZErrorE()); electron
	wvpimTrk = WTrackParameter(mpipm, pimTrk->getZHelix(), pimTrk->getZError());  // wvkpTrk = WTrackParameter(mk, pipTrk->getZHelixK(), pipTrk->getZErrorK()); kaon
																				  // wvmupTrk = WTrackParameter(mmu, pipTrk->getZHelixMu(), pipTrk-			>getZErrorMu());	 	m	uon
																				  //****************************************
#pragma endregion
#pragma region section_4C
#pragma region fourc_初始参数定义
	WTrackParameter wpip = wvpipTrk;
	WTrackParameter wpim = wvpimTrk;
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
	HepLorentzVector ecms(0.010978 * job_energy, 0, 0, job_energy);
	HepLorentzVector ptrackp;
	HepLorentzVector ptrackm;
	HepLorentzVector ptrack1;
	HepLorentzVector ptrack2;
	double chisq_4c_0g = 9999;
	double chisq_4c_1g = 9999;
	double chisq_4c_2g = 9999;
	double chisq_4c_3g = 9999;
	double chisq_4c_4g = 9999;
	int index_gamma1 = 0;
	int index_gamma2 = 0;
	int index_gamma3 = 0;
	int index_gamma4 = 0;
#pragma endregion
#pragma region fourc_2gamma拟合
	for (int i1 = 0; i1 < nGam; i1++)
	{
		RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();
		for (int i2 = i1; i2 < nGam; i2++)
		{
			if (i2 == i1)
				continue;
			RecEmcShower *g2Trk =(*(evtRecTrkCol->begin() +	iGam[i2]))->emcShower();
			kmfit->init();
			kmfit->AddTrack(0, wpip);
			kmfit->AddTrack(1, wpim);
			kmfit->AddTrack(2, 0.0, g1Trk);
			kmfit->AddTrack(3, 0.0, g2Trk);
			kmfit->AddFourMomentum(0, ecms);
			bool oksq = kmfit->Fit();
			if (oksq)
			{
				double chi2 = kmfit->chisq();
				if (chi2 <= chisq_4c_2g)
				{
					chisq_4c_2g = chi2;
					ptrackp = kmfit->pfit(0);
					ptrackm = kmfit->pfit(1);
					ptrack1 = kmfit->pfit(2);
					ptrack2 = kmfit->pfit(3);
					index_gamma1 = i1;
					index_gamma2 = i2;
				}
			}
		}
	}
#pragma endregion
#pragma re	gio	n fourc_0gamma	拟合
		if(job_do_4c_0)
	{
		kmfit->init();
		kmfit->AddTrack(0, wpip);
		kmfit->AddTrack(1, wpip);
		kmfit->AddFourMomentum(0, ecms);
		bool oksq = kmfit->Fit();
		if (oksq)
		{
			double chi2 = kmfit->chisq();
			if (chi2 < chisq_4c_0g)
			{
				chisq_4c_0g = chi2;
			}
		}
	}
#pragma endregion
#pragma region fourc_1gamma拟合
	if (job_do_4c_1)
	{
		for (int i1 = 0; i1 < nGam; i1++)
		{
			RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();
			kmfit->init();
			kmfit->AddTrack(0, wpip);
			kmfit->AddTrack(1, wpim);
			kmfit->AddTrack(2, 0.0, g1Trk);
			kmfit->AddFourMomentum(0, ecms);
			bool oksq = kmfit->Fit();
			if (oksq)
			{
				double chi2 = kmfit->chisq();
				if (chi2 <= chisq_4c_1g)
				{
					chisq_4c_1g = chi2;
				}
			}
		}
	}
#pragma endregion
#pragma region fourc_3gamma拟合
	if (job_do_4c_3 && nGam > 2)
	{
		for (int i1 = 0; i1 < nGam; i1++)
		{
			RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();
			for (int i2 = i1; i2 < nGam; i2++)
			{
				if (i2 == i1)
					continue;
				RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();
				for (int i3 = i2; i3 < nGam; i3++)
				{
					if (i3 == i1 || i3 == i2)
						continue;
					RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[i3]))->emcShower();
					kmfit->init();
					kmfit->AddTrack(0, wpip);
					kmfit->AddTrack(1, wpim);
					kmfit->AddTrack(2, 0.0, g1Trk);
					kmfit->AddTrack(3, 0.0, g2Trk);
					kmfit->AddTrack(4, 0.0, g3Trk);
					kmfit->AddFourMomentum(0, ecms);
					bool oksq = kmfit->Fit();
					if (oksq)
					{
						double chi2 = kmfit->chisq();
						if (chi2 < 	chisq_4c_3g)
							{
								chisq_4c_3g = chi2;
						}
					}
				}
			}
		}
	}
#pragma endregion
#pragma region fourc_judge
	int g2 = 0;
	if (chisq_4c_2g < 200)
	{
		g2 = 1;
	}
#pragma endregion
#pragma region fourc_4gamma特殊拟合
	if (job_do_4c_4 && nGam > 3 && g2 != 0)
	{
		double diff_gamma34 = 9999;
		for (int i1 = 0; i1 < nGam; i1++)
		{
			if (i1 == index_gamma1 || i1 == index_gamma2)
				continue;
			for (int i2 = 0; i2 < nGam; i2++)
			{
				if (i2 == i1 || i2 == index_gamma1 || i2 == index_gamma2)
					continue;
				double diff = fabs((pGam[i1] + pGam[i2]).m() - mpiz);
				if (diff < diff_gamma34)
				{
					diff_gamma34 = diff;
					index_gamma3 = i1;
					index_gamma4 = i2;
				}
			}
		}
		RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[index_gamma1]))->emcShower();
		RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[index_gamma2]))->emcShower();
		RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[index_gamma3]))->emcShower();
		RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[index_gamma4]))->emcShower();
		kmfit->init();
		kmfit->AddTrack(0, wpip);
		kmfit->AddTrack(1, wpim);
		kmfit->AddTrack(2, 0.0, g1Trk);
		kmfit->AddTrack(3, 0.0, g2Trk);
		kmfit->AddTrack(4, 0.0, g3Trk);
		kmfit->AddTrack(5, 0.0, g4Trk);
		kmfit->AddFourMomentum(0, ecms);
		bool oksq = kmfit->Fit();
		if (oksq)
		{
			chisq_4c_4g = kmfit->chisq();
		}
	}
#pragma endregion
#pragma region fourc_信息输出
	if (g2 != 0)
	{
		Ncut3++;
		my_seriescount(SeriesRun, SeriesNum3, runNo);

		HepLorentzVector fit4c_pip = ptrackp;
		HepLorentzVector fit4c_pim = ptrackm;
		HepLorentzVector fit4c_gamma1 = ptrack1;
		HepLorentzVector fit4c_gamma2 = ptrack2;
		HepLorentzVector fit4c_piz = fit4c_gamma1 + fit4c_gamma2;
		HepLorentzVector fit4c_pipm = fit4c_pip + fit4c_pim;
		HepLorentzVector fit4c_pipz = fit4c_pip + fit4c_piz;
		HepLorentzVector fit4c_pimz = fit4c_pim + fit4c_piz;
		fit4c_chisq = chisq_4c_2g;
		fit4c_chisq_0g = chisq_4c_0g;
		fit4c_chisq_1g = chisq_4c_1g;
		fit4c_chisq_3g = chisq_4c_3g;
		fit4c_chisq_4g = chisq_4c_4g;

		fit4c_egamma1 = fit4c_gamma1.e();
		fit4c_egamma2 = fit4c_gamma2.e();
		fit4c_hgamma1 = my_helicityangle(fit4c_gamma1, fit4c_piz);
		fit4c_hgamma2 = my_helicityangle(fit4c_gamma2, fit4c_piz);
		fit4c_agamma =	my_angle(fit4c_pip, fit4c_pim);

		fit4c_mpip = fit4c_pip.m();
		fit4c_apip = fit4c_pip.cosTheta();
		fit4c_ppip = fit4c_pip.rho();
		fit4c_mpim = fit4c_pim.m();
		fit4c_apim = fit4c_pim.cosTheta();
		fit4c_ppim = fit4c_pim.rho();
		fit4c_mpiz = fit4c_piz.m();
		fit4c_apiz = fit4c_piz.cosTheta();
		fit4c_ppiz = fit4c_piz.rho();
		fit4c_mpipz = fit4c_pipz.m();
		fit4c_apipz = fit4c_pipz.cosTheta();
		fit4c_ppipz = fit4c_pipz.rho();
		fit4c_mpimz = fit4c_pimz.m();
		fit4c_apimz = fit4c_pimz.cosTheta();
		fit4c_ppimz = fit4c_pimz.rho();
		fit4c_mpipm = fit4c_pipm.m();
		fit4c_apipm = fit4c_pipm.cosTheta();
		fit4c_ppipm = fit4c_pipm.rho();

		m_tuple4->write();
		Ncut4++;
		my_seriescount(SeriesRun, SeriesNum4, runNo);
	}
#pragma endregion
#pragma endregion
	return StatusCode::SUCCESS;
}
#pragma endregion
#pragma region 结束输出
//*******************************************************************************************************
//***                                               finalize                                            ***
//*********************************************************************************************************
StatusCode Pppmpz::finalize()
{
	cout << "energy point:         " << job_energy << endl;
	cout << "total number:         " << Ncut0 << endl;
	cout << "Pass truth:           " << Ncut1 << endl;
	cout << "Pass Pid:             " << Ncut2 << endl;
	cout << "Pass 4C:              " << Ncut3 << endl;
	cout << "Pass 5C:              " << Ncut4 << endl;
	// Start my output
	cout << "****************************************************" << endl;
	cout << "**********Exporting Run numbers and events**********" << endl;
	for (int i = 0; i < SeriesRun.size(); i++)
	{
		cout << "Run:" << SeriesRun[i] << "for" << SeriesNum[i] << "times" << endl;
	}
	cout << "****************************************************" << endl;
	for (int i = 0; i < SeriesRun.size(); i++)
	{
		cout << "Run:" << SeriesRun[i] << "get signal1" << SeriesNum1[i] << "times" << endl;
	}
	cout << "****************************************************" << endl;
	for (int i = 0; i < SeriesRun.size(); i++)
	{
		cout << "Run:" << SeriesRun[i] << "get signal2" << SeriesNum2[i] << "times" << endl;
	}
	cout << "****************************************************" << endl;
	for (int i = 0; i < SeriesRun.size(); i++)
	{
		cout << "Run:" << SeriesRun[i] << "get signal3" << SeriesNum3[i] << "times" << endl;
	}
	cout << "****************************************************" << endl;
	for (int i = 0; i < SeriesRun.size(); i++)
	{
		cout << "Run:" << SeriesRun[i] << "get signal4" << SeriesNum4[i] << "times" << endl;
	}
	cout << "****************************************************" << endl;
	for (int i = 0; i < SeriesRun.size(); i++)
	{
		cout << "Run:" << SeriesRun[i] << "get signal5" << SeriesNum5[i] << "times" << endl;
	}
	cout << "****************************************************" << endl;
	cout << "Finish script" << endl;
	// End my output
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}
#pragma endregion
