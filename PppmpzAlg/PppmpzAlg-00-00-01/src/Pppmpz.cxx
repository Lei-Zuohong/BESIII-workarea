#pragma region 1.1.引用函数库
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
#include "PppmpzAlg/headc/bes.h"
#include "PppmpzAlg/headc/common.h"
#pragma endregion
#pragma region 1.2.初始化变量
// 定义类型
typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;
// 定义常数
bes::CONST use_const;
bes::RUN use_run;
Pppmpz::Pppmpz(const std::string &name, ISvcLocator *pSvcLocator) : Algorithm(name, pSvcLocator)
{
	declareProperty("Energy", job_energy = 0);
	declareProperty("flag3", job_flag3 = 0);

	declareProperty("do_truth", job_do_truth = 1);
	declareProperty("do_4c", job_do_4c = 1);
	declareProperty("do_4c_0", job_do_4c_0 = 0);
	declareProperty("do_4c_1", job_do_4c_1 = 0);
	declareProperty("do_4c_3", job_do_4c_3 = 0);
	declareProperty("do_4c_4", job_do_4c_4 = 0);
}
#pragma endregion
#pragma region 1.3.初始化输出
StatusCode Pppmpz::initialize()
{
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in initialize()" << endmsg;
	StatusCode status;
#pragma region 新建树：truth
	if (job_do_truth)
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
				status = m_tuple_truth->addItem("gamma1_m", truth_gamma1_m);
				status = m_tuple_truth->addItem("gamma1_p", truth_gamma1_p);
				status = m_tuple_truth->addItem("gamma1_a", truth_gamma1_a);
				status = m_tuple_truth->addItem("gamma1_pe", truth_gamma1_pe);
				status = m_tuple_truth->addItem("gamma1_px", truth_gamma1_px);
				status = m_tuple_truth->addItem("gamma1_py", truth_gamma1_py);
				status = m_tuple_truth->addItem("gamma1_pz", truth_gamma1_pz);
				status = m_tuple_truth->addItem("gamma2_m", truth_gamma2_m);
				status = m_tuple_truth->addItem("gamma2_p", truth_gamma2_p);
				status = m_tuple_truth->addItem("gamma2_a", truth_gamma2_a);
				status = m_tuple_truth->addItem("gamma2_pe", truth_gamma2_pe);
				status = m_tuple_truth->addItem("gamma2_px", truth_gamma2_px);
				status = m_tuple_truth->addItem("gamma2_py", truth_gamma2_py);
				status = m_tuple_truth->addItem("gamma2_pz", truth_gamma2_pz);
				status = m_tuple_truth->addItem("piz_m", truth_piz_m);
				status = m_tuple_truth->addItem("piz_p", truth_piz_p);
				status = m_tuple_truth->addItem("piz_a", truth_piz_a);
				status = m_tuple_truth->addItem("piz_pe", truth_piz_pe);
				status = m_tuple_truth->addItem("piz_px", truth_piz_px);
				status = m_tuple_truth->addItem("piz_py", truth_piz_py);
				status = m_tuple_truth->addItem("piz_pz", truth_piz_pz);
				status = m_tuple_truth->addItem("pipm_m", truth_pipm_m);
				status = m_tuple_truth->addItem("pipm_p", truth_pipm_p);
				status = m_tuple_truth->addItem("pipm_a", truth_pipm_a);
				status = m_tuple_truth->addItem("pipm_pe", truth_pipm_pe);
				status = m_tuple_truth->addItem("pipm_px", truth_pipm_px);
				status = m_tuple_truth->addItem("pipm_py", truth_pipm_py);
				status = m_tuple_truth->addItem("pipm_pz", truth_pipm_pz);
				status = m_tuple_truth->addItem("pipz_m", truth_pipz_m);
				status = m_tuple_truth->addItem("pipz_p", truth_pipz_p);
				status = m_tuple_truth->addItem("pipz_a", truth_pipz_a);
				status = m_tuple_truth->addItem("pipz_pe", truth_pipz_pe);
				status = m_tuple_truth->addItem("pipz_px", truth_pipz_px);
				status = m_tuple_truth->addItem("pipz_py", truth_pipz_py);
				status = m_tuple_truth->addItem("pipz_pz", truth_pipz_pz);
				status = m_tuple_truth->addItem("pimz_m", truth_pimz_m);
				status = m_tuple_truth->addItem("pimz_p", truth_pimz_p);
				status = m_tuple_truth->addItem("pimz_a", truth_pimz_a);
				status = m_tuple_truth->addItem("pimz_pe", truth_pimz_pe);
				status = m_tuple_truth->addItem("pimz_px", truth_pimz_px);
				status = m_tuple_truth->addItem("pimz_py", truth_pimz_py);
				status = m_tuple_truth->addItem("pimz_pz", truth_pimz_pz);
			}
			else
			{
				log << MSG::ERROR << "Cannot book N-tuple:" << long(m_tuple_truth) << endmsg;
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
			status = m_tuple_charge->addItem("nCharge", charge_nCharge);
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
	if (job_do_4c)
	{
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
				// Topology信息
				status = m_tuple_fit4c->addItem("runID", runNo);
				status = m_tuple_fit4c->addItem("eventID", eventNo);
				status = m_tuple_fit4c->addItem("indexmc", m_idxmc, 0, 100);
				status = m_tuple_fit4c->addIndexedItem("pdgid", m_idxmc, m_pdgid);
				status = m_tuple_fit4c->addIndexedItem("motheridx", m_idxmc, m_motheridx);
				status = m_tuple_fit4c->addItem("flag1", flag1);
				status = m_tuple_fit4c->addItem("flag2", flag2);
				status = m_tuple_fit4c->addItem("flag3", flag3);
				// Charged track
				status = m_tuple_fit4c->addItem("pip_ep", fit4c_pip_ep);
				status = m_tuple_fit4c->addItem("pim_ep", fit4c_pim_ep);
				status = m_tuple_fit4c->addItem("pip_pid_pi", fit4c_pip_pid_pi);
				status = m_tuple_fit4c->addItem("pim_pid_pi", fit4c_pim_pid_pi);
				status = m_tuple_fit4c->addItem("pip_pid_mu", fit4c_pip_pid_mu);
				status = m_tuple_fit4c->addItem("pim_pid_mu", fit4c_pim_pid_mu);
				status = m_tuple_fit4c->addItem("pip_pid_e", fit4c_pip_pid_e);
				status = m_tuple_fit4c->addItem("pim_pid_e", fit4c_pim_pid_e);
				// Neutral track
				status = m_tuple_fit4c->addItem("ngamma", fit4c_ngamma);
				// Vertex fit
				status = m_tuple_fit4c->addItem("vertex", fit4c_vertex);
				// 4C fit
				status = m_tuple_fit4c->addItem("chisq", fit4c_chisq);
				status = m_tuple_fit4c->addItem("chisq_0g", fit4c_chisq_0g);
				status = m_tuple_fit4c->addItem("chisq_1g", fit4c_chisq_1g);
				status = m_tuple_fit4c->addItem("chisq_3g", fit4c_chisq_3g);
				status = m_tuple_fit4c->addItem("chisq_4g", fit4c_chisq_4g);
				// Momentum transform
				status = m_tuple_fit4c->addItem("gamma1_heli", fit4c_gamma1_heli);
				status = m_tuple_fit4c->addItem("gamma2_heli", fit4c_gamma2_heli);
				status = m_tuple_fit4c->addItem("a_pippim", fit4c_a_pippim);
				status = m_tuple_fit4c->addItem("b_pippim", fit4c_b_pippim);
				status = m_tuple_fit4c->addItem("dalitz_pm", fit4c_dalitz_pm);
				status = m_tuple_fit4c->addItem("dalitz_pz", fit4c_dalitz_pz);
				status = m_tuple_fit4c->addItem("dalitz_mz", fit4c_dalitz_mz);
				// Momentum
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
				status = m_tuple_fit4c->addItem("gamma1_m", fit4c_gamma1_m);
				status = m_tuple_fit4c->addItem("gamma1_p", fit4c_gamma1_p);
				status = m_tuple_fit4c->addItem("gamma1_a", fit4c_gamma1_a);
				status = m_tuple_fit4c->addItem("gamma1_pe", fit4c_gamma1_pe);
				status = m_tuple_fit4c->addItem("gamma1_px", fit4c_gamma1_px);
				status = m_tuple_fit4c->addItem("gamma1_py", fit4c_gamma1_py);
				status = m_tuple_fit4c->addItem("gamma1_pz", fit4c_gamma1_pz);
				status = m_tuple_fit4c->addItem("gamma2_m", fit4c_gamma2_m);
				status = m_tuple_fit4c->addItem("gamma2_p", fit4c_gamma2_p);
				status = m_tuple_fit4c->addItem("gamma2_a", fit4c_gamma2_a);
				status = m_tuple_fit4c->addItem("gamma2_pe", fit4c_gamma2_pe);
				status = m_tuple_fit4c->addItem("gamma2_px", fit4c_gamma2_px);
				status = m_tuple_fit4c->addItem("gamma2_py", fit4c_gamma2_py);
				status = m_tuple_fit4c->addItem("gamma2_pz", fit4c_gamma2_pz);
				status = m_tuple_fit4c->addItem("piz_m", fit4c_piz_m);
				status = m_tuple_fit4c->addItem("piz_p", fit4c_piz_p);
				status = m_tuple_fit4c->addItem("piz_a", fit4c_piz_a);
				status = m_tuple_fit4c->addItem("piz_pe", fit4c_piz_pe);
				status = m_tuple_fit4c->addItem("piz_px", fit4c_piz_px);
				status = m_tuple_fit4c->addItem("piz_py", fit4c_piz_py);
				status = m_tuple_fit4c->addItem("piz_pz", fit4c_piz_pz);
				status = m_tuple_fit4c->addItem("pipm_m", fit4c_pipm_m);
				status = m_tuple_fit4c->addItem("pipm_p", fit4c_pipm_p);
				status = m_tuple_fit4c->addItem("pipm_a", fit4c_pipm_a);
				status = m_tuple_fit4c->addItem("pipm_pe", fit4c_pipm_pe);
				status = m_tuple_fit4c->addItem("pipm_px", fit4c_pipm_px);
				status = m_tuple_fit4c->addItem("pipm_py", fit4c_pipm_py);
				status = m_tuple_fit4c->addItem("pipm_pz", fit4c_pipm_pz);
				status = m_tuple_fit4c->addItem("pipz_m", fit4c_pipz_m);
				status = m_tuple_fit4c->addItem("pipz_p", fit4c_pipz_p);
				status = m_tuple_fit4c->addItem("pipz_a", fit4c_pipz_a);
				status = m_tuple_fit4c->addItem("pipz_pe", fit4c_pipz_pe);
				status = m_tuple_fit4c->addItem("pipz_px", fit4c_pipz_px);
				status = m_tuple_fit4c->addItem("pipz_py", fit4c_pipz_py);
				status = m_tuple_fit4c->addItem("pipz_pz", fit4c_pipz_pz);
				status = m_tuple_fit4c->addItem("pimz_m", fit4c_pimz_m);
				status = m_tuple_fit4c->addItem("pimz_p", fit4c_pimz_p);
				status = m_tuple_fit4c->addItem("pimz_a", fit4c_pimz_a);
				status = m_tuple_fit4c->addItem("pimz_pe", fit4c_pimz_pe);
				status = m_tuple_fit4c->addItem("pimz_px", fit4c_pimz_px);
				status = m_tuple_fit4c->addItem("pimz_py", fit4c_pimz_py);
				status = m_tuple_fit4c->addItem("pimz_pz", fit4c_pimz_pz);
			}
			else
			{
				log << MSG::ERROR << "Cannot book N-tuple:" << long(m_tuple_fit4c) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
#pragma endregion
	log << MSG::INFO << "successfully return from initialize()" << endmsg;
	return StatusCode::SUCCESS;
}
#pragma endregion
#pragma region 2.0.循环执行事例
StatusCode Pppmpz::execute() //
{
#pragma region section_初始化
	// 初始化运行
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;
	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
	runNo = eventHeader->runNumber();
	eventNo = eventHeader->eventNumber();
	use_run.INIT(eventHeader->runNumber());
	// 处理flag信息
	if (eventHeader->runNumber() < 0)
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
	log << MSG::DEBUG << "run, evtnum = "
		<< eventHeader->runNumber() << " , "
		<< eventHeader->eventNumber() << endreq;
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
		if (eventHeader->runNumber() == 34326 || eventHeader->runNumber() == 34334 || eventHeader->runNumber() == 34478 || eventHeader->runNumber() == 34818 || eventHeader->runNumber() == 34982 ||
			eventHeader->runNumber() == 35101 || eventHeader->runNumber() == 40459 || eventHeader->runNumber() == 40460 || eventHeader->runNumber() == 40461 || eventHeader->runNumber() == 40462 ||
			eventHeader->runNumber() == 41408 || eventHeader->runNumber() == 41416 || eventHeader->runNumber() == 41902)
		{
			m_status = -1;
		}
		if (m_status == -1)
		{
			return StatusCode::SUCCESS;
		}
	}
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
	if (eventHeader->runNumber() < 0 && job_do_truth)
	{
		// 设定变量
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
		// 设定直接轨迹变量 - track
		HepLorentzVector truth_isr;
		HepLorentzVector truth_pip;
		HepLorentzVector truth_pim;
		HepLorentzVector truth_piz;
		HepLorentzVector truth_gamma1;
		HepLorentzVector truth_gamma2;
		// 设定直接轨迹变量 - number
		int n_isr = 0, index_isr;
		int n_pip = 0, index_pip;
		int n_pim = 0, index_pim;
		int n_piz = 0, index_piz;
		int n_gamma1 = 0, index_gamma1;
		int n_gamma2 = 0, index_gamma2;
		// 设定间接轨迹变量
		HepLorentzVector truth_pipm;
		HepLorentzVector truth_pipz;
		HepLorentzVector truth_pimz;
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
			truth_pipm = truth_pip + truth_pim;
			truth_pipz = truth_pip + truth_piz;
			truth_pimz = truth_pim + truth_piz;
			bes::tracktovalue(truth_pip, truth_pip_m, truth_pip_p, truth_pip_a, truth_pip_pe, truth_pip_px, truth_pip_py, truth_pip_pz);
			bes::tracktovalue(truth_pim, truth_pim_m, truth_pim_p, truth_pim_a, truth_pim_pe, truth_pim_px, truth_pim_py, truth_pim_pz);
			bes::tracktovalue(truth_gamma1, truth_gamma1_m, truth_gamma1_p, truth_gamma1_a, truth_gamma1_pe, truth_gamma1_px, truth_gamma1_py, truth_gamma1_pz);
			bes::tracktovalue(truth_gamma2, truth_gamma2_m, truth_gamma2_p, truth_gamma2_a, truth_gamma2_pe, truth_gamma2_px, truth_gamma2_py, truth_gamma2_pz);
			bes::tracktovalue(truth_piz, truth_piz_m, truth_piz_p, truth_piz_a, truth_piz_pe, truth_piz_px, truth_piz_py, truth_piz_pz);
			bes::tracktovalue(truth_pipm, truth_pipm_m, truth_pipm_p, truth_pipm_a, truth_pipm_pe, truth_pipm_px, truth_pipm_py, truth_pipm_pz);
			bes::tracktovalue(truth_pipz, truth_pipz_m, truth_pipz_p, truth_pipz_a, truth_pipz_pe, truth_pipz_px, truth_pipz_py, truth_pipz_pz);
			bes::tracktovalue(truth_pimz, truth_pimz_m, truth_pimz_p, truth_pimz_a, truth_pimz_pe, truth_pimz_px, truth_pimz_py, truth_pimz_pz);
			m_tuple_truth->write();
			use_run.COUNT(eventHeader->runNumber(), 1);
		}
	}
#pragma endregion
#pragma region section_charged track
	Vint iCharge, ipip, ipim;										//
	iCharge.clear();												// 变量：iCharge[]（参数为good-track序号，内容为track编号）
	ipip.clear();													// 变量：ipip[]（参数为good-track+序号，内容为track编号）
	ipim.clear();													// 变量：ipim[]（参数为good-track-序号，内容为track编号）
	Vdouble eppip, eppim;											//
	eppip.clear();													// 变量：eppip[]（参数为good-track+序号，内容为e/p ratio）
	eppim.clear();													// 变量：eppim[]（参数为good-track-序号，内容为e/p ratio）
	int nCharge = 0;												//
	int npip = 0;													//
	int npim = 0;													//
	Hep3Vector xorigin(0, 0, 0);									//
	IVertexDbSvc *vtxsvc;											//
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);			//
	if (vtxsvc->isVertexValid())									//
	{																//
		double *dbv = vtxsvc->PrimaryVertex();						//
		double *vv = vtxsvc->SigmaPrimaryVertex();					//
		xorigin.setX(dbv[0]);										//
		xorigin.setY(dbv[1]);										//
		xorigin.setZ(dbv[2]);										//
	}																//
	for (int i = 0; i < evtRecEvent->totalCharged(); i++)			//
	{																//
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;		//
		if (!(*itTrk)->isMdcTrackValid())							//
			continue;												//
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();					//
		double pch = mdcTrk->p();									//
		double x0 = mdcTrk->x();									//
		double y0 = mdcTrk->y();									//
		double z0 = mdcTrk->z();									//
		double phi0 = mdcTrk->helix(1);								// 表示螺旋线参数
		double xv = xorigin.x();									// 0: d0 -> 螺旋线到对撞顶点的最小距离
		double yv = xorigin.y();									// 1: phi0 -> 最小距离的xy平面相角
		double Rxy = (x0 - xv) * cos(phi0) + (y0 - yv) * sin(phi0); // 2: kappa
		HepVector a = mdcTrk->helix();								// 3: d
		HepSymMatrix Ea = mdcTrk->err();							// 4: tan(lamda)
		HepPoint3D point0(0., 0., 0.);								//
		HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);			//
		VFHelix helixip(point0, a, Ea);								//
		helixip.pivot(IP);											//
		HepVector vecipa = helixip.a();								//
		double Rvxy0 = fabs(vecipa[0]);								//
		double Rvz0 = vecipa[3];									//
		double Rvphi0 = vecipa[1];									//
		if (fabs(cos(mdcTrk->theta())) > 0.93)						// 选择：cos(theta)
			continue;												//
		if (fabs(Rvz0) >= 10)										// 选择：Rvz0
			continue;												//
		if (fabs(Rvxy0) >= 1)										// 选择：Rvxy0
			continue;												//
		double epratio = 9999;										//
		if ((*itTrk)->isEmcShowerValid())							//
		{															//
			RecEmcShower *emcTrk = (*itTrk)->emcShower();			//
			epratio = emcTrk->energy() / mdcTrk->p();				//
		}															//
		if (mdcTrk->charge() > 0)									//
		{															//
			iCharge.push_back(i);									//
			ipip.push_back(i);										//
			eppip.push_back(epratio);								//
		}															//
		else if (mdcTrk->charge() < 0)								//
		{															//
			iCharge.push_back(i);									//
			ipim.push_back(i);										//
			eppim.push_back(epratio);								//
		}															//
	}																//
	nCharge = iCharge.size();										//
	npip = ipip.size();												//
	npim = ipim.size();												//
	if (1 == 1)														//
	{																//
		charge_nCharge = nCharge;									//
		m_tuple_charge->write();									//
	}																//
	if ((npip != 1) || (npim != 1))									//
	{																//
		return StatusCode::SUCCESS;									//
	}																//
#pragma endregion
#pragma region section_charged track momentum
	Vp4 ppip, ppim;																									  //
	ppip.clear();																									  // 变量：ppip[]（good-track+的四动量）
	ppim.clear();																									  // 变量：ppim[]（good-track-的四动量）
	Vdouble pip_pid_pi, pim_pid_pi;																					  //
	Vdouble pip_pid_mu, pim_pid_mu;																					  //
	Vdouble pip_pid_e, pim_pid_e;																					  //
	pip_pid_pi.clear();																								  // 变量：pip_pid_?（track+的PID概率）
	pim_pid_pi.clear();																								  // 变量：pim_pid_?（track-的PID概率）
	pip_pid_mu.clear();																								  //
	pim_pid_mu.clear();																								  //
	pip_pid_e.clear();																								  //
	pim_pid_e.clear();																								  //
	ParticleID *pid = ParticleID::instance();																		  //
	for (int i = 0; i < nCharge; i++)																				  //
	{																												  //
		int pid_check = 1;																							  //
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iCharge[i];												  //
		pid->init();																								  //
		pid->setMethod(pid->methodProbability());																	  // 对于Likelihood方法：pid->setMethod(pid->methodLikelihood());
		pid->setChiMinCut(4);																						  //
		pid->setRecTrack(*itTrk);																					  //
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE());							  //
		pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyProton() | pid->onlyElectron() | pid->onlyMuon()); //
		pid->calculate();																							  //
		if (!(pid->IsPidInfoValid()))																				  //
			pid_check = 0;																							  // 对于Likelihood方法(0=electron 1=muon 2=pion 3=kaon 4=proton)
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();																	  // if(pid->pdf(2) < pid->pdf(3))
		if (pid->probPion() < pid->probKaon())																		  //
			pid_check = 0;																							  //
		if (pid->probPion() < pid->probProton())																	  //
			pid_check = 0;																							  //
		if (pid->probPion() < pid->probElectron())																	  //
			pid_check = 0;																							  //
		RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();														  // 对于ParticleID, 用RecMdcKalTrack代替RecMdcTrack
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);															  // PID可以设定为electron, muon, pion, kaon and proton;The default setting is pion
		HepLorentzVector ptrk;																						  //
		ptrk.setPx(mdcKalTrk->px());																				  //
		ptrk.setPy(mdcKalTrk->py());																				  //
		ptrk.setPz(mdcKalTrk->pz());																				  //
		double p3 = ptrk.mag();																						  //
		ptrk.setE(sqrt(p3 * p3 + use_const.mpipm * use_const.mpipm));												  //
		if (mdcKalTrk->charge() > 0)																				  //
		{																											  //
			ppip.push_back(ptrk);																					  //
			pip_pid_pi.push_back(pid->probPion());																	  //
			pip_pid_mu.push_back(pid->probMuon());																	  //
			pip_pid_e.push_back(pid->probElectron());																  //
		}																											  //
		else																										  //
		{																											  //
			ppim.push_back(ptrk);																					  //
			pim_pid_pi.push_back(pid->probPion());																	  //
			pim_pid_mu.push_back(pid->probMuon());																	  //
			pim_pid_e.push_back(pid->probElectron());																  //
		}																											  //
	}																												  //
	use_run.COUNT(eventHeader->runNumber(), 2);																		  //
#pragma endregion
#pragma region section_neutral track
	Vint iGam;																							  //
	iGam.clear();																						  //
	Vdouble aGam;																						  //
	aGam.clear();																						  //
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
		if (fabs(dang) < 10)																			  //
			continue;																					  //
		iGam.push_back(i);																				  // 变量：iGam[]（参数为good-track序号，内容为track编号）
		aGam.push_back(fabs(dang));																		  // 变量：aGam[]（参数为good-track序号，内容为shower和charged最小夹角）
	}																									  //
	int nGam = iGam.size();																				  // 变量：nGam（中性track数量）
	if (nGam < 2 || nGam > 20)																			  //
	{																									  //
		return StatusCode::SUCCESS;																		  //
	}																									  //
#pragma endregion
#pragma region section_neutral track momentum
	Vp4 pGam;														 // 变量：pGam[]（参数为good-track序号，内容为动量）
	pGam.clear();													 //
	for (int i = 0; i < nGam; i++)									 //
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
	use_run.COUNT(eventHeader->runNumber(), 3);						 //
#pragma endregion
#pragma region section_vertex fit
	RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin() + ipip[0]))->mdcKalTrack();	   // ****************************************
	RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin() + ipim[0]))->mdcKalTrack();	   // Default is pion, for other particles:
	WTrackParameter wpip, wpim;														   // wvppTrk = WTrackParameter(mp, pipTrk->getZHelixP(), pipTrk->getZErrorP()); proton
	wpip = WTrackParameter(use_const.mpipm, pipTrk->getZHelix(), pipTrk->getZError()); // wvepTrk = WTrackParameter(me, pipTrk->getZHelixE(), pipTrk->getZErrorE()); electron
	wpim = WTrackParameter(use_const.mpipm, pimTrk->getZHelix(), pimTrk->getZError()); // wvkpTrk = WTrackParameter(mk, pipTrk->getZHelixK(), pipTrk->getZErrorK()); kaon
	HepPoint3D vx(0., 0., 0.);														   // wvmupTrk = WTrackParameter(mmu, pipTrk->getZHelixMu(), pipTrk->getZErrorMu());	muon
	HepSymMatrix Evx(3, 0);															   //****************************************
	double bx = 1E+6;																   //
	double by = 1E+6;																   //
	double bz = 1E+6;																   //
	Evx[0][0] = bx * bx;															   //
	Evx[1][1] = by * by;															   //
	Evx[2][2] = bz * bz;															   //
	VertexParameter vxpar;															   //
	vxpar.setVx(vx);																   //
	vxpar.setEvx(Evx);																   //
	VertexFit *vtxfit = VertexFit::instance();										   //
	vtxfit->init();																	   //
	vtxfit->AddTrack(0, wpip);														   //
	vtxfit->AddTrack(1, wpim);														   //
	vtxfit->AddVertex(0, vxpar, 0, 1);												   //
	double chisq_vertex = 9999;														   //
	if (vtxfit->Fit(0))																   //
		chisq_vertex = vtxfit->chisq(0);											   //
	if (1 == 1)																		   //
	{																				   //
		vertex_chisq = chisq_vertex;												   //
		m_tuple_vertex->write();													   //
	}																				   //
	fit4c_vertex = chisq_vertex;													   //
	vtxfit->Swim(0);																   //
	use_run.COUNT(eventHeader->runNumber(), 4);										   //
#pragma endregion
#pragma region fourc_初始参数定义
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
#pragma endregion
#pragma region fourc_2gamma拟合
	for (int i1 = 0; i1 < nGam; i1++)
	{
		RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();
		for (int i2 = i1; i2 < nGam; i2++)
		{
			if (i2 == i1)
				continue;
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();
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
#pragma region fourc_0gamma拟合
	if (job_do_4c_0 && g2 != 0)
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
	if (job_do_4c_1 && g2 != 0)
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
	if (job_do_4c_3 && nGam > 2 && g2 != 0)
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
						if (chi2 < chisq_4c_3g)
						{
							chisq_4c_3g = chi2;
						}
					}
				}
			}
		}
	}
#pragma endregion
#pragma region fourc_4gamma拟合
	if (job_do_4c_3 && nGam > 2 && g2 != 0)
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
					for (int i4 = i3; i4 < nGam; i4++)
					{
						if (i4 == i1 || i4 == i2 || i4 == i3)
							continue;
						RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[i4]))->emcShower();
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
				}
			}
		}
	}
#pragma endregion
#pragma region fourc_信息输出
	if (g2 != 0)
	{
		// ?.带电track
		fit4c_pip_ep = eppip[0];
		fit4c_pim_ep = eppim[0];
		fit4c_pip_pid_pi = pip_pid_pi[0];
		fit4c_pim_pid_pi = pim_pid_pi[0];
		fit4c_pip_pid_mu = pip_pid_mu[0] / pip_pid_pi[0];
		fit4c_pim_pid_mu = pim_pid_mu[0] / pim_pid_pi[0];
		fit4c_pip_pid_e = pip_pid_e[0] / pip_pid_pi[0];
		fit4c_pim_pid_e = pim_pid_e[0] / pim_pid_pi[0];
		// ?.中性track
		fit4c_ngamma = nGam;
		// ?.事例选择
		fit4c_chisq = chisq_4c_2g;
		fit4c_chisq_0g = chisq_4c_0g;
		fit4c_chisq_1g = chisq_4c_1g;
		fit4c_chisq_3g = chisq_4c_3g;
		fit4c_chisq_4g = chisq_4c_4g;
		// ?.重建track
		HepLorentzVector fit4c_pip = ptrackp;
		HepLorentzVector fit4c_pim = ptrackm;
		HepLorentzVector fit4c_gamma1 = ptrack1;
		HepLorentzVector fit4c_gamma2 = ptrack2;
		HepLorentzVector fit4c_piz = fit4c_gamma1 + fit4c_gamma2;
		HepLorentzVector fit4c_pipm = fit4c_pip + fit4c_pim;
		HepLorentzVector fit4c_pipz = fit4c_pip + fit4c_piz;
		HepLorentzVector fit4c_pimz = fit4c_pim + fit4c_piz;
		// ?.track计算
		fit4c_gamma1_heli = bes::helicityangle(fit4c_gamma1, fit4c_piz);
		fit4c_gamma2_heli = bes::helicityangle(fit4c_gamma2, fit4c_piz);
		fit4c_a_pippim = bes::angle_boost(fit4c_pip, fit4c_pim) * 180. / 3.1415926;
		fit4c_b_pippim = bes::angle(fit4c_pip, fit4c_pim) * 180. / 3.1415926;
		// ?.track取值
		bes::tracktovalue(fit4c_pip, fit4c_pip_m, fit4c_pip_p, fit4c_pip_a, fit4c_pip_pe, fit4c_pip_px, fit4c_pip_py, fit4c_pip_pz);
		bes::tracktovalue(fit4c_pim, fit4c_pim_m, fit4c_pim_p, fit4c_pim_a, fit4c_pim_pe, fit4c_pim_px, fit4c_pim_py, fit4c_pim_pz);
		bes::tracktovalue(fit4c_gamma1, fit4c_gamma1_m, fit4c_gamma1_p, fit4c_gamma1_a, fit4c_gamma1_pe, fit4c_gamma1_px, fit4c_gamma1_py, fit4c_gamma1_pz);
		bes::tracktovalue(fit4c_gamma2, fit4c_gamma2_m, fit4c_gamma2_p, fit4c_gamma2_a, fit4c_gamma2_pe, fit4c_gamma2_px, fit4c_gamma2_py, fit4c_gamma2_pz);
		bes::tracktovalue(fit4c_piz, fit4c_piz_m, fit4c_piz_p, fit4c_piz_a, fit4c_piz_pe, fit4c_piz_px, fit4c_piz_py, fit4c_piz_pz);
		bes::tracktovalue(fit4c_pipm, fit4c_pipm_m, fit4c_pipm_p, fit4c_pipm_a, fit4c_pipm_pe, fit4c_pipm_px, fit4c_pipm_py, fit4c_pipm_pz);
		bes::tracktovalue(fit4c_pipz, fit4c_pipz_m, fit4c_pipz_p, fit4c_pipz_a, fit4c_pipz_pe, fit4c_pipz_px, fit4c_pipz_py, fit4c_pipz_pz);
		bes::tracktovalue(fit4c_pimz, fit4c_pimz_m, fit4c_pimz_p, fit4c_pimz_a, fit4c_pimz_pe, fit4c_pimz_px, fit4c_pimz_py, fit4c_pimz_pz);
		// 变量处理：动量归一化
		fit4c_pip_p = fit4c_pip_p * 2 / job_energy;
		fit4c_pim_p = fit4c_pim_p * 2 / job_energy;
		fit4c_gamma1_p = fit4c_gamma1_p * 2 / job_energy;
		fit4c_gamma2_p = fit4c_gamma2_p * 2 / job_energy;
		fit4c_piz_p = fit4c_piz_p * 2 / job_energy;
		fit4c_pipm_p = fit4c_pipm_p * 2 / job_energy;
		fit4c_pipz_p = fit4c_pipz_p * 2 / job_energy;
		fit4c_pimz_p = fit4c_pimz_p * 2 / job_energy;
		// 变量处理：dalitz变量
		fit4c_dalitz_pm = fit4c_pipm_m * fit4c_pipm_m;
		fit4c_dalitz_pz = fit4c_pipz_m * fit4c_pipz_m;
		fit4c_dalitz_mz = fit4c_pimz_m * fit4c_pimz_m;

		if (fit4c_pim_ep < 1.2 && fit4c_pip_ep < 1.2)
		{
			m_tuple_fit4c->write();
			use_run.COUNT(eventHeader->runNumber(), 5);
		}
	}
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
	cout << "Energy point:         " << job_energy << endl;
	cout << "Total number:         " << use_run.countnum[0] << endl;
	cout << "Pass truth:           " << use_run.countnum[1] << endl;
	cout << "Pass charged track:   " << use_run.countnum[2] << endl;
	cout << "Pass pid:             " << use_run.countnum[3] << endl;
	cout << "Pass vertex fit:      " << use_run.countnum[4] << endl;
	cout << "Pass 4C               " << use_run.countnum[5] << endl;
	// Start my output
	cout << "****************************************************" << endl;
	cout << "**********Exporting Run numbers and events**********" << endl;
	for (int i = 0; i < use_run.seriesrun.size(); i++)
	{
		cout << "Run: " << use_run.seriesrun[i] << " for " << use_run.seriesnum[0][i] << " times" << endl;
	}
	cout << "****************************************************" << endl;
	for (int i = 0; i < use_run.seriesrun.size(); i++)
	{
		cout << "Run: " << use_run.seriesrun[i] << " signal1 " << use_run.seriesnum[1][i] << " times" << endl;
	}
	cout << "****************************************************" << endl;
	for (int i = 0; i < use_run.seriesrun.size(); i++)
	{
		cout << "Run: " << use_run.seriesrun[i] << " signal2 " << use_run.seriesnum[2][i] << " times" << endl;
	}
	cout << "****************************************************" << endl;
	for (int i = 0; i < use_run.seriesrun.size(); i++)
	{
		cout << "Run: " << use_run.seriesrun[i] << " signal3 " << use_run.seriesnum[3][i] << " times" << endl;
	}
	cout << "****************************************************" << endl;
	for (int i = 0; i < use_run.seriesrun.size(); i++)
	{
		cout << "Run: " << use_run.seriesrun[i] << " signal4 " << use_run.seriesnum[4][i] << " times" << endl;
	}
	cout << "****************************************************" << endl;
	for (int i = 0; i < use_run.seriesrun.size(); i++)
	{
		cout << "Run: " << use_run.seriesrun[i] << " signal5 " << use_run.seriesnum[5][i] << " times" << endl;
	}
	cout << "****************************************************" << endl;
	cout << "Finish script" << endl;
	// End my output
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}
#pragma endregion
