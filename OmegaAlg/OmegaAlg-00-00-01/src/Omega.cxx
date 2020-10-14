#pragma region 准备：调用头文件
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
#pragma endregion
#pragma region 1.2.调用类型
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
#pragma region 1.3.定义全局变量
typedef std::vector<int> Vint;				  // 定义类型
typedef std::vector<double> Vdouble;		  //
typedef std::vector<HepLorentzVector> Vp4;	  //
my_constant use_constant;					  // 定义常数
const double mpiz = use_constant.mpiz;		  //
const double mpipm = use_constant.mpipm;	  //
int Ncut0, Ncut1, Ncut2, Ncut3, Ncut4, Ncut5; // 定义统计总数的参数
Vint SeriesRun;								  //
Vint SeriesNum;								  //
Vint SeriesNum1;							  //
Vint SeriesNum2;							  //
Vint SeriesNum3;							  //
Vint SeriesNum4;							  //
Vint SeriesNum5;							  //
int firstrun = 0;							  //
#pragma endregion
#pragma region 1.4.定义输入变量
Omega::Omega(const std::string &name, ISvcLocator *pSvcLocator) : Algorithm(name, pSvcLocator)
{
	declareProperty("Energy", job_energy = 0);
	declareProperty("truth", job_truth = 1);
	declareProperty("do_567", job_do_567 = 1);
}
#pragma endregion
#pragma region 1.5.定义输出
StatusCode Omega::initialize()
{
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in initialize()" << endmsg;
	StatusCode status;
#pragma region 新建树：truth
	if (job_truth == 1)
	{
		NTuplePtr nt_truth(ntupleSvc(), "FILE1/truth");
		if (nt_truth)
		{
			m_tuple_truth = nt_truth;
		}
		else
		{
			m_tuple_truth = ntupleSvc()->book("FILE1/truth", CLID_ColumnWiseTuple, "ks N-Tuple example");
			if (m_tuple_truth)
			{
				status = m_tuple_truth->addItem("isr_m", truth_isr_m);
				status = m_tuple_truth->addItem("isr_p", truth_isr_p);
				status = m_tuple_truth->addItem("isr_a", truth_isr_a);
				status = m_tuple_truth->addItem("isr_pe", truth_isr_pe);
				status = m_tuple_truth->addItem("isr_px", truth_isr_px);
				status = m_tuple_truth->addItem("isr_py", truth_isr_py);
				status = m_tuple_truth->addItem("isr_pz", truth_isr_pz);
				status = m_tuple_truth->addItem("pip_m", truth_pip_m);
				status = m_tuple_truth->addItem("pip_p", truth_pip_p);
				status = m_tuple_truth->addItem("pip_a", truth_pip_a);
				status = m_tuple_truth->addItem("pip_pe", truth_pip_pe);
				status = m_tuple_truth->addItem("pip_px", truth_pip_px);
				status = m_tuple_truth->addItem("pip_py", truth_pip_py);
				status = m_tuple_truth->addItem("pip_pz", truth_pip_pz);
				status = m_tuple_truth->addItem("pim_m", truth_pim_m);
				status = m_tuple_truth->addItem("pim_p", truth_pim_p);
				status = m_tuple_truth->addItem("pim_a", truth_pim_a);
				status = m_tuple_truth->addItem("pim_pe", truth_pim_pe);
				status = m_tuple_truth->addItem("pim_px", truth_pim_px);
				status = m_tuple_truth->addItem("pim_py", truth_pim_py);
				status = m_tuple_truth->addItem("pim_pz", truth_pim_pz);
				status = m_tuple_truth->addItem("pi01_m", truth_pi01_m);
				status = m_tuple_truth->addItem("pi01_p", truth_pi01_p);
				status = m_tuple_truth->addItem("pi01_a", truth_pi01_a);
				status = m_tuple_truth->addItem("pi01_pe", truth_pi01_pe);
				status = m_tuple_truth->addItem("pi01_px", truth_pi01_px);
				status = m_tuple_truth->addItem("pi01_py", truth_pi01_py);
				status = m_tuple_truth->addItem("pi01_pz", truth_pi01_pz);
				status = m_tuple_truth->addItem("pi02_m", truth_pi02_m);
				status = m_tuple_truth->addItem("pi02_p", truth_pi02_p);
				status = m_tuple_truth->addItem("pi02_a", truth_pi02_a);
				status = m_tuple_truth->addItem("pi02_pe", truth_pi02_pe);
				status = m_tuple_truth->addItem("pi02_px", truth_pi02_px);
				status = m_tuple_truth->addItem("pi02_py", truth_pi02_py);
				status = m_tuple_truth->addItem("pi02_pz", truth_pi02_pz);
				status = m_tuple_truth->addItem("pi03_m", truth_pi03_m);
				status = m_tuple_truth->addItem("pi03_p", truth_pi03_p);
				status = m_tuple_truth->addItem("pi03_a", truth_pi03_a);
				status = m_tuple_truth->addItem("pi03_pe", truth_pi03_pe);
				status = m_tuple_truth->addItem("pi03_px", truth_pi03_px);
				status = m_tuple_truth->addItem("pi03_py", truth_pi03_py);
				status = m_tuple_truth->addItem("pi03_pz", truth_pi03_pz);
				status = m_tuple_truth->addItem("omega_m", truth_omega_m);
				status = m_tuple_truth->addItem("omega_p", truth_omega_p);
				status = m_tuple_truth->addItem("omega_a", truth_omega_a);
				status = m_tuple_truth->addItem("omega_pe", truth_omega_pe);
				status = m_tuple_truth->addItem("omega_px", truth_omega_px);
				status = m_tuple_truth->addItem("omega_py", truth_omega_py);
				status = m_tuple_truth->addItem("omega_pz", truth_omega_pz);
				status = m_tuple_truth->addItem("omegapi02_m", truth_omegapi02_m);
				status = m_tuple_truth->addItem("omegapi02_p", truth_omegapi02_p);
				status = m_tuple_truth->addItem("omegapi02_a", truth_omegapi02_a);
				status = m_tuple_truth->addItem("omegapi02_pe", truth_omegapi02_pe);
				status = m_tuple_truth->addItem("omegapi02_px", truth_omegapi02_px);
				status = m_tuple_truth->addItem("omegapi02_py", truth_omegapi02_py);
				status = m_tuple_truth->addItem("omegapi02_pz", truth_omegapi02_pz);
				status = m_tuple_truth->addItem("omegapi03_m", truth_omegapi03_m);
				status = m_tuple_truth->addItem("omegapi03_p", truth_omegapi03_p);
				status = m_tuple_truth->addItem("omegapi03_a", truth_omegapi03_a);
				status = m_tuple_truth->addItem("omegapi03_pe", truth_omegapi03_pe);
				status = m_tuple_truth->addItem("omegapi03_px", truth_omegapi03_px);
				status = m_tuple_truth->addItem("omegapi03_py", truth_omegapi03_py);
				status = m_tuple_truth->addItem("omegapi03_pz", truth_omegapi03_pz);
				status = m_tuple_truth->addItem("pi02pi03_m", truth_pi02pi03_m);
				status = m_tuple_truth->addItem("pi02pi03_p", truth_pi02pi03_p);
				status = m_tuple_truth->addItem("pi02pi03_a", truth_pi02pi03_a);
				status = m_tuple_truth->addItem("pi02pi03_pe", truth_pi02pi03_pe);
				status = m_tuple_truth->addItem("pi02pi03_px", truth_pi02pi03_px);
				status = m_tuple_truth->addItem("pi02pi03_py", truth_pi02pi03_py);
				status = m_tuple_truth->addItem("pi02pi03_pz", truth_pi02pi03_pz);
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
	NTuplePtr nt_charge(ntupleSvc(), "FILE1/charge");
	if (nt_charge)
	{
		m_tuple_charge = nt_charge;
	}
	else
	{
		m_tuple_charge = ntupleSvc()->book("FILE1/charge", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if (m_tuple_charge)
		{
			status = m_tuple_charge->addItem("ngood", charge_ngood);
			status = m_tuple_charge->addItem("ncharge", charge_ncharge);
		}
		else
		{
			log << MSG::ERROR << "Cannot book N-tuple:" << long(m_tuple_charge) << endmsg;
			return StatusCode::FAILURE;
		}
	}
#pragma endregion
#pragma region 新建树：vertex
	NTuplePtr nt_vertex(ntupleSvc(), "FILE1/vertex");
	if (nt_vertex)
	{
		m_tuple_vertex = nt_vertex;
	}
	else
	{
		m_tuple_vertex = ntupleSvc()->book("FILE1/vertex", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if (m_tuple_vertex)
		{
			status = m_tuple_vertex->addItem("chisq", vertex_chisq);
		}
		else
		{
			log << MSG::ERROR << "Cannot book N-tuple:" << long(m_tuple_vertex) << endmsg;
			return StatusCode::FAILURE;
		}
	}
#pragma endregion
#pragma region 新建树：fit4c
	NTuplePtr nt_fit4c(ntupleSvc(), "FILE1/fit4c");
	if (nt_fit4c)
	{
		m_tuple_fit4c = nt_fit4c;
	}
	else
	{
		m_tuple_fit4c = ntupleSvc()->book("FILE1/fit4c", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if (m_tuple_fit4c)
		{
			// topo 信息
			status = m_tuple_fit4c->addItem("runID", runID);
			status = m_tuple_fit4c->addItem("eventID", eventID);
			status = m_tuple_fit4c->addItem("indexmc", m_idxmc, 0, 100);
			status = m_tuple_fit4c->addIndexedItem("pdgid", m_idxmc, m_pdgid);
			status = m_tuple_fit4c->addIndexedItem("motheridx", m_idxmc, m_motheridx);
			//
			status = m_tuple_fit4c->addItem("chisq", fit4c_chisq);
			status = m_tuple_fit4c->addItem("pip_m", fit4c_pip_m);
			status = m_tuple_fit4c->addItem("pip_p", fit4c_pip_p);
			status = m_tuple_fit4c->addItem("pip_a", fit4c_pip_a);
			status = m_tuple_fit4c->addItem("pip_pe", fit4c_pip_pe);
			status = m_tuple_fit4c->addItem("pip_px", fit4c_pip_px);
			status = m_tuple_fit4c->addItem("pip_py", fit4c_pip_py);
			status = m_tuple_fit4c->addItem("pip_pz", fit4c_pip_pz);
			status = m_tuple_fit4c->addItem("pim_m", fit4c_pim_m);
			status = m_tuple_fit4c->addItem("pim_p", fit4c_pim_p);
			status = m_tuple_fit4c->addItem("pim_a", fit4c_pim_a);
			status = m_tuple_fit4c->addItem("pim_pe", fit4c_pim_pe);
			status = m_tuple_fit4c->addItem("pim_px", fit4c_pim_px);
			status = m_tuple_fit4c->addItem("pim_py", fit4c_pim_py);
			status = m_tuple_fit4c->addItem("pim_pz", fit4c_pim_pz);
			status = m_tuple_fit4c->addItem("pi01_m", fit4c_pi01_m);
			status = m_tuple_fit4c->addItem("pi01_p", fit4c_pi01_p);
			status = m_tuple_fit4c->addItem("pi01_a", fit4c_pi01_a);
			status = m_tuple_fit4c->addItem("pi01_pe", fit4c_pi01_pe);
			status = m_tuple_fit4c->addItem("pi01_px", fit4c_pi01_px);
			status = m_tuple_fit4c->addItem("pi01_py", fit4c_pi01_py);
			status = m_tuple_fit4c->addItem("pi01_pz", fit4c_pi01_pz);
			status = m_tuple_fit4c->addItem("pi02_m", fit4c_pi02_m);
			status = m_tuple_fit4c->addItem("pi02_p", fit4c_pi02_p);
			status = m_tuple_fit4c->addItem("pi02_a", fit4c_pi02_a);
			status = m_tuple_fit4c->addItem("pi02_pe", fit4c_pi02_pe);
			status = m_tuple_fit4c->addItem("pi02_px", fit4c_pi02_px);
			status = m_tuple_fit4c->addItem("pi02_py", fit4c_pi02_py);
			status = m_tuple_fit4c->addItem("pi02_pz", fit4c_pi02_pz);
			status = m_tuple_fit4c->addItem("pi03_m", fit4c_pi03_m);
			status = m_tuple_fit4c->addItem("pi03_p", fit4c_pi03_p);
			status = m_tuple_fit4c->addItem("pi03_a", fit4c_pi03_a);
			status = m_tuple_fit4c->addItem("pi03_pe", fit4c_pi03_pe);
			status = m_tuple_fit4c->addItem("pi03_px", fit4c_pi03_px);
			status = m_tuple_fit4c->addItem("pi03_py", fit4c_pi03_py);
			status = m_tuple_fit4c->addItem("pi03_pz", fit4c_pi03_pz);
			status = m_tuple_fit4c->addItem("omega_m", fit4c_omega_m);
			status = m_tuple_fit4c->addItem("omega_p", fit4c_omega_p);
			status = m_tuple_fit4c->addItem("omega_a", fit4c_omega_a);
			status = m_tuple_fit4c->addItem("omega_pe", fit4c_omega_pe);
			status = m_tuple_fit4c->addItem("omega_px", fit4c_omega_px);
			status = m_tuple_fit4c->addItem("omega_py", fit4c_omega_py);
			status = m_tuple_fit4c->addItem("omega_pz", fit4c_omega_pz);
			status = m_tuple_fit4c->addItem("omegapi02_m", fit4c_omegapi02_m);
			status = m_tuple_fit4c->addItem("omegapi02_p", fit4c_omegapi02_p);
			status = m_tuple_fit4c->addItem("omegapi02_a", fit4c_omegapi02_a);
			status = m_tuple_fit4c->addItem("omegapi02_pe", fit4c_omegapi02_pe);
			status = m_tuple_fit4c->addItem("omegapi02_px", fit4c_omegapi02_px);
			status = m_tuple_fit4c->addItem("omegapi02_py", fit4c_omegapi02_py);
			status = m_tuple_fit4c->addItem("omegapi02_pz", fit4c_omegapi02_pz);
			status = m_tuple_fit4c->addItem("omegapi03_m", fit4c_omegapi03_m);
			status = m_tuple_fit4c->addItem("omegapi03_p", fit4c_omegapi03_p);
			status = m_tuple_fit4c->addItem("omegapi03_a", fit4c_omegapi03_a);
			status = m_tuple_fit4c->addItem("omegapi03_pe", fit4c_omegapi03_pe);
			status = m_tuple_fit4c->addItem("omegapi03_px", fit4c_omegapi03_px);
			status = m_tuple_fit4c->addItem("omegapi03_py", fit4c_omegapi03_py);
			status = m_tuple_fit4c->addItem("omegapi03_pz", fit4c_omegapi03_pz);
			status = m_tuple_fit4c->addItem("pi02pi03_m", fit4c_pi02pi03_m);
			status = m_tuple_fit4c->addItem("pi02pi03_p", fit4c_pi02pi03_p);
			status = m_tuple_fit4c->addItem("pi02pi03_a", fit4c_pi02pi03_a);
			status = m_tuple_fit4c->addItem("pi02pi03_pe", fit4c_pi02pi03_pe);
			status = m_tuple_fit4c->addItem("pi02pi03_px", fit4c_pi02pi03_px);
			status = m_tuple_fit4c->addItem("pi02pi03_py", fit4c_pi02pi03_py);
			status = m_tuple_fit4c->addItem("pi02pi03_pz", fit4c_pi02pi03_pz);

			status = m_tuple_fit4c->addItem("pif1_m", fit4c_pif1_m);
			status = m_tuple_fit4c->addItem("pif2_m", fit4c_pif2_m);
			status = m_tuple_fit4c->addItem("pif3_m", fit4c_pif3_m);
			status = m_tuple_fit4c->addItem("pif4_m", fit4c_pif4_m);
		}
		else
		{
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple_fit4c) << endmsg;
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
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;
	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
	int runNo = eventHeader->runNumber();
	int event = eventHeader->eventNumber();
	runID = runNo;
	eventID = event;
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
	if (runNo > 0)
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
#pragma endregion
#pragma region section_topoloty
	if (eventHeader->runNumber() < 0)
	{
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
		int m_numParticle = 1;
		if (!mcParticleCol)
		{
			cout << "Could not retrieve McParticelCol" << endl;
			return StatusCode::FAILURE;
		}
		else
		{
			int nprmary = 0;
			Event::McParticleCol::iterator iter_mc1 = mcParticleCol->begin();
			for (; iter_mc1 != mcParticleCol->end(); iter_mc1++)
			{
				if (!(*iter_mc1)->decayFromGenerator())
					continue;
				if ((*iter_mc1)->primaryParticle())
				{
					nprmary++;
				}
			}
			Event::McParticleCol::iterator iter_mc2 = mcParticleCol->begin();
			if (nprmary == 1)
			{
				m_numParticle = 0;
				for (; iter_mc2 != mcParticleCol->end(); iter_mc2++)
				{
					if (!(*iter_mc2)->decayFromGenerator())
						continue;
					if ((*iter_mc2)->primaryParticle())
					{
						m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();
						m_motheridx[m_numParticle] = 0;
					}
					else
					{
						m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();
						m_motheridx[m_numParticle] = ((*iter_mc2)->mother()).trackIndex();
					}
					m_numParticle += 1;
				}
				m_idxmc = m_numParticle;
			}
			if (nprmary > 1)
			{
				m_numParticle = 1;
				for (; iter_mc2 != mcParticleCol->end(); iter_mc2++)
				{
					if (!(*iter_mc2)->decayFromGenerator())
						continue;
					if ((*iter_mc2)->primaryParticle())
					{
						m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();
						m_motheridx[m_numParticle] = 0;
					}
					else
					{

						m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();
						m_motheridx[m_numParticle] = ((*iter_mc2)->mother()).trackIndex() + 1;
					}
					m_numParticle += 1;
					m_pdgid[0] = 11111;
					m_motheridx[0] = 0;
				}
				m_idxmc = m_numParticle;
			}
		}
	}
#pragma endregion
#pragma region section_truth
	int truth_check = 0;
	if (eventHeader->runNumber() < 0)
	{
		// 设定变量
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
		// 设定直接轨迹变量 - track
		HepLorentzVector truth_isr, truth_pip, truth_pim;
		HepLorentzVector truth_gamma1, truth_gamma2, truth_gamma3, truth_gamma4, truth_gamma5, truth_gamma6;
		HepLorentzVector truth_pi01, truth_pi02, truth_pi03, truth_omega;
		HepLorentzVector truth_omegapi02;
		HepLorentzVector truth_omegapi03;
		HepLorentzVector truth_pi02pi03;
		// 设定直接轨迹变量 - number
		int n_isr = 0, n_pip = 0, n_pim = 0;
		int n_gamma1 = 0, n_gamma2 = 0, n_gamma3 = 0, n_gamma4 = 0, n_gamma5 = 0, n_gamma6 = 0;
		int n_pi01 = 0, n_pi02 = 0, n_pi03 = 0, n_omega = 0;
		// 设定直接轨迹变量 - index
		int index_isr, index_pip, index_pim;
		int index_gamma1, index_gamma2, index_gamma3, index_gamma4, index_gamma5, index_gamma6;
		int index_pi01, index_pi02, index_pi03, index_omega;
		// 设定暂存交换变量
		HepLorentzVector truth_medium;
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
				truth_omega = mctrue_track;
				index_omega = mctrue_index;
				n_omega += 1;
			}
			// 统计pi+,标记211,且来自223
			if ((*iter_mc)->particleProperty() == 211 && ((*iter_mc)->mother()).particleProperty() == 223)
			{
				truth_pip = mctrue_track;
				index_pip = mctrue_index;
				n_pip += 1;
			}
			// 统计pi-,标记-211,且来自223
			if ((*iter_mc)->particleProperty() == -211 && ((*iter_mc)->mother()).particleProperty() == 223)
			{
				truth_pim = mctrue_track;
				index_pim = mctrue_index;
				n_pim += 1;
			}
			// 统计pi0,标记111,来自223定为pi01,否则定为pi02,pi03
			if ((*iter_mc)->particleProperty() == 111)
			{
				if (((*iter_mc)->mother()).particleProperty() == 223)
				{
					truth_pi01 = mctrue_track;
					index_pi01 = mctrue_index;
					n_pi01 += 1;
				}
				else if (n_pi02 == 0)
				{
					truth_pi02 = mctrue_track;
					index_pi02 = mctrue_index;
					n_pi02 += 1;
				}
				else
				{
					truth_pi03 = mctrue_track;
					index_pi03 = mctrue_index;
					n_pi03 += 1;
				}
			}
			// 统计gamma,标记22
			if ((*iter_mc)->particleProperty() == 22 && ((*iter_mc)->mother()).particleProperty() != 111)
			{
				truth_isr = mctrue_track;
				index_isr = mctrue_index;
				n_isr += 1;
			}
		}
		// 定位 pi02 pi03
		if (1 == 1)
		{
			truth_omegapi02 = truth_omega + truth_pi02;
			truth_omegapi03 = truth_omega + truth_pi03;
			if (truth_omegapi02.m() > truth_omegapi03.m())
			{
				truth_medium = truth_pi03;
				truth_pi03 = truth_pi02;
				truth_pi02 = truth_medium;
				index_medium = index_pi03;
				index_pi03 = index_pi02;
				index_pi02 = index_medium;
			}
			truth_omegapi02 = truth_omega + truth_pi02;
			truth_omegapi03 = truth_omega + truth_pi03;
			truth_pi02pi03 = truth_pi02 + truth_pi03;
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
						truth_gamma1 = mctrue_track;
						index_gamma1 = mctrue_index;
						n_gamma1++;
					}
					else
					{
						truth_gamma2 = mctrue_track;
						index_gamma2 = mctrue_index;
						n_gamma2++;
					}
				}
				if (((*iter_mc)->mother()).trackIndex() == index_pi02)
				{
					if (n_gamma3 == 0)
					{
						truth_gamma3 = mctrue_track;
						index_gamma3 = mctrue_index;
						n_gamma3++;
					}
					else
					{
						truth_gamma4 = mctrue_track;
						index_gamma4 = mctrue_index;
						n_gamma4++;
					}
				}
				if (((*iter_mc)->mother()).trackIndex() == index_pi03)
				{
					if (n_gamma5 == 0)
					{
						truth_gamma5 = mctrue_track;
						index_gamma5 = mctrue_index;
						n_gamma5++;
					}
					else
					{
						truth_gamma6 = mctrue_track;
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
		// 填入信息
		if (truth_check == 1)
		{
			// 输出gamma信息
			if (1 == 1)
			{
				my_tracktovalue(truth_isr, truth_isr_m, truth_isr_p, truth_isr_a, truth_isr_pe, truth_isr_px, truth_isr_py, truth_isr_pz);
				my_tracktovalue(truth_pip, truth_pip_m, truth_pip_p, truth_pip_a, truth_pip_pe, truth_pip_px, truth_pip_py, truth_pip_pz);
				my_tracktovalue(truth_pim, truth_pim_m, truth_pim_p, truth_pim_a, truth_pim_pe, truth_pim_px, truth_pim_py, truth_pim_pz);
				my_tracktovalue(truth_pi01, truth_pi01_m, truth_pi01_p, truth_pi01_a, truth_pi01_pe, truth_pi01_px, truth_pi01_py, truth_pi01_pz);
				my_tracktovalue(truth_pi02, truth_pi02_m, truth_pi02_p, truth_pi02_a, truth_pi02_pe, truth_pi02_px, truth_pi02_py, truth_pi02_pz);
				my_tracktovalue(truth_pi03, truth_pi03_m, truth_pi03_p, truth_pi03_a, truth_pi03_pe, truth_pi03_px, truth_pi03_py, truth_pi03_pz);
				my_tracktovalue(truth_omega, truth_omega_m, truth_omega_p, truth_omega_a, truth_omega_pe, truth_omega_px, truth_omega_py, truth_omega_pz);
				my_tracktovalue(truth_omegapi02, truth_omegapi02_m, truth_omegapi02_p, truth_omegapi02_a, truth_omegapi02_pe, truth_omegapi02_px, truth_omegapi02_py, truth_omegapi02_pz);
				my_tracktovalue(truth_omegapi03, truth_omegapi03_m, truth_omegapi03_p, truth_omegapi03_a, truth_omegapi03_pe, truth_omegapi03_px, truth_omegapi03_py, truth_omegapi03_pz);
				my_tracktovalue(truth_pi02pi03, truth_pi02pi03_m, truth_pi02pi03_p, truth_pi02pi03_a, truth_pi02pi03_pe, truth_pi02pi03_px, truth_pi02pi03_py, truth_pi02pi03_pz);
				m_tuple_truth->write();
				Ncut1 += 1;
				my_seriescount(SeriesRun, SeriesNum1, runNo);
			}
		}
	}
#pragma endregion
#pragma region section_charged track
	Vint iGood, ipip, ipim;															   // ****************************************
	iGood.clear();																	   // 变量：iGood[]（参数为good-track序号，内容为track编号）
	ipip.clear();																	   // 变量：ipip[]（参数为good-track+序号，内容为track编号）
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
		iGood.push_back(i);															   //
		nCharge += mdcTrk->charge();												   //
	}																				   //
	int nGood = iGood.size();														   //
	log << MSG::DEBUG << "ngood, totcharge = " << nGood << " , " << nCharge << endreq; // ****************************************
	if (1 == 1)																		   // 输出2：charge
	{																				   // ****************************************
		charge_ngood = nGood;														   //
		charge_ncharge = nCharge;													   //
		m_tuple_charge->write();													   //
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
	if (nGam < 6 || nGam > 10)		// 选择：nGam
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
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];					   // ****************************************
		pid->init();																	   // 对于Likelihood方法
		pid->setMethod(pid->methodProbability());										   // pid->setMethod(pid->methodLikelihood());
		pid->setChiMinCut(4);															   //
		pid->setRecTrack(*itTrk);														   //
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE()); // 选择pid的粒子
		pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyProton());			   // pid->identify(pid->onlyPion());
		pid->calculate();																   // pid->identify(pid->onlyKaon());
		if (!(pid->IsPidInfoValid()))													   //
			continue;																	   // 对于Likelihood方法(0=electron 1=muon 2=pion 3=kaon 4=proton)
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();										   // if(pid->pdf(2) < pid->pdf(3))
		if (pid->probPion() < pid->probKaon())											   // continue;
			continue;																	   //
		if (pid->probPion() < pid->probProton())										   //
			continue;																	   //
		RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();							   // 对于ParticleID, 用RecMdcKalTrack代替RecMdcTrack
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);								   // PID可以设定为electron, muon, pion, kaon and proton;The default setting is pion
		if (mdcKalTrk->charge() > 0)													   // ****************************************
		{																				   //
			ipip.push_back(iGood[i]);													   // ****************************************
			HepLorentzVector ptrk;														   // 变量：ppip[]（值为pi+动量）
			ptrk.setPx(mdcKalTrk->px());												   // 变量：ppim[]（值为pi-动量）
			ptrk.setPy(mdcKalTrk->py());												   // 变量：ipip[]（值为pi+编号）
			ptrk.setPz(mdcKalTrk->pz());												   // 变量：ipim[]（值为pi-编号）
			double p3 = ptrk.mag();														   // ****************************************
			ptrk.setE(sqrt(p3 * p3 + mpipm * mpipm));									   //
			ppip.push_back(ptrk);														   //
		}																				   //
		else																			   //
		{																				   //
			ipim.push_back(iGood[i]);													   //
			HepLorentzVector ptrk;														   //
			ptrk.setPx(mdcKalTrk->px());												   //
			ptrk.setPy(mdcKalTrk->py());												   //
			ptrk.setPz(mdcKalTrk->pz());												   //
			double p3 = ptrk.mag();														   //
			ptrk.setE(sqrt(p3 * p3 + mpipm * mpipm));									   //
			ppim.push_back(ptrk);														   //
		}																				   //
	}																					   // ****************************************
	int npip = ipip.size();																   // 变量：npip（pi+数目）
	int npim = ipim.size();																   // 变量：npim（pi-数目）
	if (npip * npim != 1)																   // ****************************************
		return SUCCESS;																	   //
	Ncut2++;																			   // Ncut2
	my_seriescount(SeriesRun, SeriesNum2, runNo);										   // Series2
#pragma endregion
#pragma region section_vertex fit
	RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin() + ipip[0]))->mdcKalTrack(); //
	RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin() + ipim[0]))->mdcKalTrack(); // Default is pion, for other particles:
	WTrackParameter wvpipTrk, wvpimTrk;											  // wvppTrk = WTrackParameter(mp, pipTrk->getZHelixP(), pipTrk->getZErrorP()); proton
	wvpipTrk = WTrackParameter(mpipm, pipTrk->getZHelix(), pipTrk->getZError());  // wvepTrk = WTrackParameter(me, pipTrk->getZHelixE(), pipTrk->getZErrorE()); electron
	wvpimTrk = WTrackParameter(mpipm, pimTrk->getZHelix(), pimTrk->getZError());  // wvkpTrk = WTrackParameter(mk, pipTrk->getZHelixK(), pipTrk->getZErrorK()); kaon
	HepPoint3D vx(0., 0., 0.);													  // wvmupTrk = WTrackParameter(mmu, pipTrk->getZHelixMu(), pipTrk->getZErrorMu()); muon
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
	vtxfit->AddTrack(0, wvpipTrk);												  //
	vtxfit->AddTrack(1, wvpimTrk);												  //
	vtxfit->AddVertex(0, vxpar, 0, 1);											  //
	if (!vtxfit->Fit(0))														  //
		return SUCCESS;															  //
	double chisq_vertex = vtxfit->chisq(0);										  //
	if (chisq_vertex > 100)														  //
		return SUCCESS;															  //
	if (1 == 1)																	  //
	{																			  //
		vertex_chisq = chisq_vertex;											  //
		m_tuple_vertex->write();												  //
	}																			  //
	vtxfit->Swim(0);															  //
#pragma endregion
#pragma region fourc_初始参数定义
	WTrackParameter wpip = vtxfit->wtrk(0);
	WTrackParameter wpim = vtxfit->wtrk(1);
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
	HepLorentzVector ecms(0.010978 * job_energy, 0, 0, job_energy);
	HepLorentzVector ptrackp;
	HepLorentzVector ptrackm;
	HepLorentzVector ptrack1;
	HepLorentzVector ptrack2;
	HepLorentzVector ptrack3;
	HepLorentzVector ptrack4;
	HepLorentzVector ptrack5;
	HepLorentzVector ptrack6;
	double chisq_4c_5g = 9999;
	double chisq_4c_6g = 9999;
	double chisq_4c_7g = 9999;
#pragma endregion
#pragma region fourc_6gamma拟合
	for (int i1 = 0; i1 < nGam; i1++)
	{
		RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();
		for (int i2 = i1; i2 < nGam; i2++)
		{
			if (i2 == i1)
			{
				continue;
			}
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();
			for (int i3 = i2; i3 < nGam; i3++)
			{
				if (i3 == i1 || i3 == i2)
				{
					continue;
				}
				RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[i3]))->emcShower();
				for (int i4 = i3; i4 < nGam; i4++)
				{
					if (i4 == i1 || i4 == i2 || i4 == i3)
					{
						continue;
					}
					RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[i4]))->emcShower();
					for (int i5 = i4; i5 < nGam; i5++)
					{
						if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4)
						{
							continue;
						}
						RecEmcShower *g5Trk = (*(evtRecTrkCol->begin() + iGam[i5]))->emcShower();
						for (int i6 = i5; i6 < nGam; i6++)
						{
							if (i6 == i1 || i6 == i2 || i6 == i3 || i6 == i4 || i6 == i5)
							{
								continue;
							}
							RecEmcShower *g6Trk = (*(evtRecTrkCol->begin() + iGam[i6]))->emcShower();
							kmfit->init();
							kmfit->AddTrack(0, wpip);
							kmfit->AddTrack(1, wpim);
							kmfit->AddTrack(2, 0.0, g1Trk);
							kmfit->AddTrack(3, 0.0, g2Trk);
							kmfit->AddTrack(4, 0.0, g3Trk);
							kmfit->AddTrack(5, 0.0, g4Trk);
							kmfit->AddTrack(6, 0.0, g5Trk);
							kmfit->AddTrack(7, 0.0, g6Trk);
							kmfit->AddFourMomentum(0, ecms);
							bool oksq = kmfit->Fit();
							if (oksq)
							{
								double chi2 = kmfit->chisq();
								if (chi2 <= chisq_4c_6g)
								{
									chisq_4c_6g = chi2;
									ptrackp = kmfit->pfit(0);
									ptrackm = kmfit->pfit(1);
									ptrack1 = kmfit->pfit(2);
									ptrack2 = kmfit->pfit(3);
									ptrack3 = kmfit->pfit(4);
									ptrack4 = kmfit->pfit(5);
									ptrack5 = kmfit->pfit(6);
									ptrack6 = kmfit->pfit(7);
								}
							}
						}
					}
				}
			}
		}
	}
#pragma endregion
#pragma region fourc_5gamma拟合
	for (int i1 = 0; i1 < nGam; i1++)
	{
		RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();
		for (int i2 = i1; i2 < nGam; i2++)
		{
			if (i2 == i1)
			{
				continue;
			}
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();
			for (int i3 = i2; i3 < nGam; i3++)
			{
				if (i3 == i1 || i3 == i2)
				{
					continue;
				}
				RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[i3]))->emcShower();
				for (int i4 = i3; i4 < nGam; i4++)
				{
					if (i4 == i1 || i4 == i2 || i4 == i3)
					{
						continue;
					}
					RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[i4]))->emcShower();
					for (int i5 = i4; i5 < nGam; i5++)
					{
						if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4)
						{
							continue;
						}
						RecEmcShower *g5Trk = (*(evtRecTrkCol->begin() + iGam[i5]))->emcShower();
						kmfit->init();
						kmfit->AddTrack(0, wpip);
						kmfit->AddTrack(1, wpim);
						kmfit->AddTrack(2, 0.0, g1Trk);
						kmfit->AddTrack(3, 0.0, g2Trk);
						kmfit->AddTrack(4, 0.0, g3Trk);
						kmfit->AddTrack(5, 0.0, g4Trk);
						kmfit->AddTrack(6, 0.0, g5Trk);
						kmfit->AddFourMomentum(0, ecms);
						bool oksq = kmfit->Fit();
						if (oksq)
						{
							double chi2 = kmfit->chisq();
							if (chi2 < chisq_4c_5g)
							{
								chisq_4c_5g = chi2;
							}
						}
					}
				}
			}
		}
	}
#pragma endregion
#pragma region fourc_7gamma拟合
	if (nGam > 6)
	{
		for (int i1 = 0; i1 < nGam; i1++)
		{
			RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();
			for (int i2 = i1; i2 < nGam; i2++)
			{
				if (i2 == i1)
				{
					continue;
				}
				RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();
				for (int i3 = i2; i3 < nGam; i3++)
				{
					if (i3 == i1 || i3 == i2)
					{
						continue;
					}
					RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[i3]))->emcShower();
					for (int i4 = i3; i4 < nGam; i4++)
					{
						if (i4 == i1 || i4 == i2 || i4 == i3)
						{
							continue;
						}
						RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[i4]))->emcShower();
						for (int i5 = i4; i5 < nGam; i5++)
						{
							if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4)
							{
								continue;
							}
							RecEmcShower *g5Trk = (*(evtRecTrkCol->begin() + iGam[i5]))->emcShower();
							for (int i6 = i5; i6 < nGam; i6++)
							{
								if (i6 == i1 || i6 == i2 || i6 == i3 || i6 == i4 || i6 == i5)
								{
									continue;
								}
								RecEmcShower *g6Trk = (*(evtRecTrkCol->begin() + iGam[i6]))->emcShower();
								for (int i7 = i6; i7 < nGam; i7++)
								{
									if (i7 == i1 || i7 == i2 || i7 == i3 || i7 == i4 || i7 == i5 || i7 == i6)
									{
										continue;
									}
									RecEmcShower *g7Trk = (*(evtRecTrkCol->begin() + iGam[i7]))->emcShower();
									kmfit->init();
									kmfit->AddTrack(0, wpip);
									kmfit->AddTrack(1, wpim);
									kmfit->AddTrack(2, 0.0, g1Trk);
									kmfit->AddTrack(3, 0.0, g2Trk);
									kmfit->AddTrack(4, 0.0, g3Trk);
									kmfit->AddTrack(5, 0.0, g4Trk);
									kmfit->AddTrack(6, 0.0, g5Trk);
									kmfit->AddTrack(7, 0.0, g6Trk);
									kmfit->AddTrack(8, 0.0, g7Trk);
									kmfit->AddFourMomentum(0, ecms);
									bool oksq = kmfit->Fit();
									if (oksq)
									{
										double chi2 = kmfit->chisq();
										if (chi2 < chisq_4c_7g)
										{
											chisq_4c_7g = chi2;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
#pragma endregion
#pragma region fourc_判断6gamma是否成功
	int check_6g = 0;
	if (chisq_4c_6g < 200)
	{
		check_6g = 1;
		Ncut3++;
		my_seriescount(SeriesRun, SeriesNum3, runNo);
	}
	if (check_6g && job_do_567)
	{
		if (chisq_4c_5g < chisq_4c_6g || chisq_4c_7g < chisq_4c_6g)
		{
			check_6g = 0;
		}
	}
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
	if (check_6g == 1)
	{
		if (1 == 1)
		{
			tracklist = my_omega_3pi(tracklist, mpiz, mpiz, mpiz);
		}
		if (1 == 1)
		{
			tracklist = my_omega_omega(tracklist, 0.782);
		}
		if (1 == 0)
		{
			tracklist = my_omega_dalitz(tracklist);
		}
		if (1 == 1)
		{
			tracklist = my_omega_lower(tracklist);
		}
#pragma endregion
#pragma region fourc_信息输出
		// gamma 信息
		HepLorentzVector fit4c_pip = tracklist[0];
		HepLorentzVector fit4c_pim = tracklist[1];
		HepLorentzVector fit4c_gamma1 = tracklist[2];
		HepLorentzVector fit4c_gamma2 = tracklist[3];
		HepLorentzVector fit4c_gamma3 = tracklist[4];
		HepLorentzVector fit4c_gamma4 = tracklist[5];
		HepLorentzVector fit4c_gamma5 = tracklist[6];
		HepLorentzVector fit4c_gamma6 = tracklist[7];
		// 粒子信息
		HepLorentzVector fit4c_pi01 = fit4c_gamma1 + fit4c_gamma2;
		HepLorentzVector fit4c_pi02 = fit4c_gamma3 + fit4c_gamma4;
		HepLorentzVector fit4c_pi03 = fit4c_gamma5 + fit4c_gamma6;
		HepLorentzVector fit4c_omega = fit4c_pip + fit4c_pim + fit4c_pi01;
		HepLorentzVector fit4c_omegapi02 = fit4c_omega + fit4c_pi02;
		HepLorentzVector fit4c_omegapi03 = fit4c_omega + fit4c_pi03;
		HepLorentzVector fit4c_pi02pi03 = fit4c_pi02 + fit4c_pi03;
		// 填入信息
		fit4c_chisq = chisq_4c_6g;
		my_tracktovalue(fit4c_pip, fit4c_pip_m, fit4c_pip_p, fit4c_pip_a, fit4c_pip_pe, fit4c_pip_px, fit4c_pip_py, fit4c_pip_pz);
		my_tracktovalue(fit4c_pim, fit4c_pim_m, fit4c_pim_p, fit4c_pim_a, fit4c_pim_pe, fit4c_pim_px, fit4c_pim_py, fit4c_pim_pz);
		my_tracktovalue(fit4c_pi01, fit4c_pi01_m, fit4c_pi01_p, fit4c_pi01_a, fit4c_pi01_pe, fit4c_pi01_px, fit4c_pi01_py, fit4c_pi01_pz);
		my_tracktovalue(fit4c_pi02, fit4c_pi02_m, fit4c_pi02_p, fit4c_pi02_a, fit4c_pi02_pe, fit4c_pi02_px, fit4c_pi02_py, fit4c_pi02_pz);
		my_tracktovalue(fit4c_pi03, fit4c_pi03_m, fit4c_pi03_p, fit4c_pi03_a, fit4c_pi03_pe, fit4c_pi03_px, fit4c_pi03_py, fit4c_pi03_pz);
		my_tracktovalue(fit4c_omega, fit4c_omega_m, fit4c_omega_p, fit4c_omega_a, fit4c_omega_pe, fit4c_omega_px, fit4c_omega_py, fit4c_omega_pz);
		my_tracktovalue(fit4c_omegapi02, fit4c_omegapi02_m, fit4c_omegapi02_p, fit4c_omegapi02_a, fit4c_omegapi02_pe, fit4c_omegapi02_px, fit4c_omegapi02_py, fit4c_omegapi02_pz);
		my_tracktovalue(fit4c_omegapi03, fit4c_omegapi03_m, fit4c_omegapi03_p, fit4c_omegapi03_a, fit4c_omegapi03_pe, fit4c_omegapi03_px, fit4c_omegapi03_py, fit4c_omegapi03_pz);
		my_tracktovalue(fit4c_pi02pi03, fit4c_pi02pi03_m, fit4c_pi02pi03_p, fit4c_pi02pi03_a, fit4c_pi02pi03_pe, fit4c_pi02pi03_px, fit4c_pi02pi03_py, fit4c_pi02pi03_pz);
		//
		fit4c_pif1_m = (fit4c_gamma1 + fit4c_gamma3).m();
		fit4c_pif2_m = (fit4c_gamma1 + fit4c_gamma4).m();
		fit4c_pif3_m = (fit4c_gamma1 + fit4c_gamma5).m();
		fit4c_pif4_m = (fit4c_gamma1 + fit4c_gamma6).m();
		//
		m_tuple_fit4c->write();
		Ncut4++;
		my_seriescount(SeriesRun, SeriesNum4, runNo);
		Ncut5++;
		my_seriescount(SeriesRun, SeriesNum5, runNo);
#pragma endregion
	}
	return StatusCode::SUCCESS;
}
#pragma endregion
#pragma region 结束输出
//*********************************************************************************************************
//***                                               finalize                                            ***
//*********************************************************************************************************
StatusCode Omega::finalize()
{
	cout << "energy point:         " << job_energy << endl;
	cout << "total number:         " << Ncut0 << endl;
	cout << "Pass truth:           " << Ncut1 << endl;
	cout << "Pass Pid:             " << Ncut2 << endl;
	cout << "Pass 4c-6gamma:       " << Ncut3 << endl;
	cout << "Pass 4C:              " << Ncut4 << endl;
	cout << "Pass 4C:              " << Ncut5 << endl;
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
