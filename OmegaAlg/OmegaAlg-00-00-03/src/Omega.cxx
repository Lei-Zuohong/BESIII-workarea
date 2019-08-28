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
const double xmass[5] = {0.000511, 0.105658, 0.139570, 0.493677, 0.938272};
const double velc = 299.792458;
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
int Ncut0, Ncut1, Ncut2, Ncut3, Ncut4, Ncut5, Ncut6;
//*********************************************************************************************************
//***                                                声明容器                                            ***
//*********************************************************************************************************
Omega::Omega(const std::string &name, ISvcLocator *pSvcLocator) : Algorithm(name, pSvcLocator)
{
	declareProperty("Test4C", m_test4C = 1);
	declareProperty("Test5C", m_test5C = 1);
	declareProperty("CheckDedx", m_checkDedx = 1);
	declareProperty("CheckTof", m_checkTof = 1);
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
	// initialize-data-in-topo
	if (1 == 1)
	{
		NTuplePtr ntt(ntupleSvc(), "FILE1/topo");
		if (ntt)
		{
			m_tuplet = ntt;
		}
		else
		{
			m_tuplet = ntupleSvc()->book("FILE1/topo", CLID_ColumnWiseTuple, "ks N-Tuple example");
			if (m_tuplet)
			{
				status = m_tuplet->addItem("runID", runID);
				status = m_tuplet->addItem("eventID", eventID);
				status = m_tuplet->addItem("indexmc", m_idxmc, 0, 100);
				status = m_tuplet->addIndexedItem("pdgid", m_idxmc, m_pdgid);
				status = m_tuplet->addIndexedItem("motheridx", m_idxmc, m_motheridx);
			}
			else
			{
				log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuplet) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
	// initialize-data-in-fit4c
	if (m_test4C == 1)
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
				status = m_tuple4->addItem("chisq", m_chisq_4c);
				status = m_tuple4->addItem("mb0", m_b0_4c);
				status = m_tuple4->addItem("momega", m_omega_4c);
				status = m_tuple4->addItem("mpi01", m_pi01_4c);
				status = m_tuple4->addItem("mpi02", m_pi02_4c);
				status = m_tuple4->addItem("mpi03", m_pi03_4c);

				status = m_tuple4->addItem("runID", runID_4c);
				status = m_tuple4->addItem("eventID", eventID_4c);
				status = m_tuple4->addItem("indexmc", m_idxmc_4c, 0, 100);
				status = m_tuple4->addIndexedItem("pdgid", m_idxmc_4c, m_pdgid_4c);
				status = m_tuple4->addIndexedItem("motheridx", m_idxmc_4c, m_motheridx_4c);
			}
			else
			{
				log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple4) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
	//initialize-data-in-fit5c
	if (m_test5C == 1)
	{
		NTuplePtr nt5(ntupleSvc(), "FILE1/fit5c");
		if (nt5)
		{
			m_tuple5 = nt5;
		}
		else
		{
			m_tuple5 = ntupleSvc()->book("FILE1/fit5c", CLID_ColumnWiseTuple, "ks N-Tuple example");
			if (m_tuple5)
			{
				status = m_tuple5->addItem("chisq", m_chisq_5c);
				status = m_tuple5->addItem("mb0", m_b0_5c);
				status = m_tuple5->addItem("momega", m_momega_5c);
				status = m_tuple5->addItem("mpi01", m_mpi01_5c);
				status = m_tuple5->addItem("mpi02", m_mpi02_5c);
				status = m_tuple5->addItem("mpi03", m_mpi03_5c);

				status = m_tuple5->addItem("runID", runID_5c);
				status = m_tuple5->addItem("eventID", eventID_5c);
				status = m_tuple5->addItem("indexmc", m_idxmc_5c, 0, 100);
				status = m_tuple5->addIndexedItem("pdgid", m_idxmc_5c, m_pdgid_5c);
				status = m_tuple5->addIndexedItem("motheridx", m_idxmc_5c, m_motheridx_5c);
			}
			else
			{
				log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple5) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
	//end-line
	log << MSG::INFO << "successfully return from initialize()" << endmsg;
	return StatusCode::SUCCESS;
}
//*********************************************************************************************************
//***                                                execute                                            ***
//*********************************************************************************************************
StatusCode Omega::execute()																	   //
{																							   //
	std::cout << "execute()" << std::endl;													   //
	MsgStream log(msgSvc(), name());														   //
	log << MSG::INFO << "in execute()" << endreq;											   //
	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");			   //
	int runNo = eventHeader->runNumber();													   // 读取runNo：runnumber
	int event = eventHeader->eventNumber();													   // 读取event：eventnumber
	runID = runNo;																			   //
	eventID = event;																		   //
	runID_4c = runNo;																		   //
	eventID_4c = event;																		   //
	runID_5c = runNo;																		   //
	eventID_5c = event;																		   //
	log << MSG::DEBUG << "run, evtnum = "													   //
		<< runNo << " , "																	   //
		<< event << endreq;																	   //
	Ncut0++;																				   //
	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);		   //
	log << MSG::DEBUG << "ncharg, nneu, tottks = "											   //
		<< evtRecEvent->totalCharged() << " , "												   //
		<< evtRecEvent->totalNeutral() << " , "												   //
		<< evtRecEvent->totalTracks() << endreq;											   //
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol); //
	//*********************************************************************************
	// Selection 0: Topology
	//*********************************************************************************
	if (eventHeader->runNumber() < 0)															  //
	{																							  //
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");  //
		int m_numParticle = 1;																	  //
		if (!mcParticleCol)																		  //
		{																						  //
			cout << "Could not retrieve McParticelCol" << endl;									  //
			return StatusCode::FAILURE;															  //
		}																						  //
		else																					  //
		{																						  //
			int nprmary = 0;																	  //
			Event::McParticleCol::iterator iter_mc1 = mcParticleCol->begin();					  //
			for (; iter_mc1 != mcParticleCol->end(); iter_mc1++)								  //
			{																					  //
				if (!(*iter_mc1)->decayFromGenerator())											  //
					continue;																	  //
				if ((*iter_mc1)->primaryParticle())												  //
				{																				  //
					nprmary++;																	  //
				}																				  //
			}																					  //
			Event::McParticleCol::iterator iter_mc2 = mcParticleCol->begin();					  //
			if (nprmary == 1)																	  //
			{																					  //
				m_numParticle = 0;																  //
				for (; iter_mc2 != mcParticleCol->end(); iter_mc2++)							  //
				{																				  //
					if (!(*iter_mc2)->decayFromGenerator())										  //
						continue;																  //
					if ((*iter_mc2)->primaryParticle())											  //
					{																			  //
						m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();				  //
						m_pdgid_4c[m_numParticle] = (*iter_mc2)->particleProperty();			  //
						m_pdgid_5c[m_numParticle] = (*iter_mc2)->particleProperty();			  //
						m_motheridx[m_numParticle] = 0;											  //
						m_motheridx_4c[m_numParticle] = 0;										  //
						m_motheridx_5c[m_numParticle] = 0;										  //
					}																			  //
					else																		  //
					{																			  //
						m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();				  //
						m_pdgid_4c[m_numParticle] = (*iter_mc2)->particleProperty();			  //
						m_pdgid_5c[m_numParticle] = (*iter_mc2)->particleProperty();			  //
						m_motheridx[m_numParticle] = ((*iter_mc2)->mother()).trackIndex();		  //
						m_motheridx_4c[m_numParticle] = ((*iter_mc2)->mother()).trackIndex();	 //
						m_motheridx_5c[m_numParticle] = ((*iter_mc2)->mother()).trackIndex();	 //
					}																			  //
					m_numParticle += 1;															  //
				}																				  //
				m_idxmc = m_numParticle;														  //
				m_idxmc_4c = m_numParticle;														  //
				m_idxmc_5c = m_numParticle;														  //
			}																					  //
			if (nprmary > 1)																	  //
			{																					  //
				m_numParticle = 1;																  //
				for (; iter_mc2 != mcParticleCol->end(); iter_mc2++)							  //
				{																				  //
					if (!(*iter_mc2)->decayFromGenerator())										  //
						continue;																  //
					if ((*iter_mc2)->primaryParticle())											  //
					{																			  //
						m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();				  //
						m_pdgid_4c[m_numParticle] = (*iter_mc2)->particleProperty();			  //
						m_pdgid_5c[m_numParticle] = (*iter_mc2)->particleProperty();			  //
						m_motheridx[m_numParticle] = 0;											  //
						m_motheridx_4c[m_numParticle] = 0;										  //
						m_motheridx_5c[m_numParticle] = 0;										  //
					}																			  //
					else																		  //
					{																			  //
																								  //
						m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();				  //
						m_pdgid_4c[m_numParticle] = (*iter_mc2)->particleProperty();			  //
						m_pdgid_5c[m_numParticle] = (*iter_mc2)->particleProperty();			  //
						m_motheridx[m_numParticle] = ((*iter_mc2)->mother()).trackIndex() + 1;	//
						m_motheridx_4c[m_numParticle] = ((*iter_mc2)->mother()).trackIndex() + 1; //
						m_motheridx_5c[m_numParticle] = ((*iter_mc2)->mother()).trackIndex() + 1; //
					}																			  //
					m_numParticle += 1;															  //
					m_pdgid[0] = 11111;															  //
					m_pdgid_4c[0] = 11111;														  //
					m_pdgid_5c[0] = 11111;														  //
					m_motheridx[0] = 0;															  //
					m_motheridx_4c[0] = 0;														  //
					m_motheridx_5c[0] = 0;														  //
				}																				  //
				m_idxmc = m_numParticle;														  //
				m_idxmc_4c = m_numParticle;														  //
				m_idxmc_5c = m_numParticle;														  //
			}																					  //
		}																						  //
	}																							  //
	//*********************************************************************************
	// Selection 1: Good Charged Track Selection
	//*********************************************************************************
	Vint iGood, ipip, ipim;															   //
	iGood.clear();																	   //存good的charge track的编号
	ipip.clear();																	   //存good的Pi+的编号
	ipim.clear();																	   //存good的Pi-的编号
	Vp4 ppip, ppim;																	   //
	ppip.clear();																	   //存Pi+的四动量
	ppim.clear();																	   //存Pi-的四动量
	int nCharge;																	   //
	nCharge = 0;																	   //存带电总量
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
	}																				   //存平均对撞顶点
	for (int i = 0; i < evtRecEvent->totalCharged(); i++)							   //循环所有带电track
	{																				   //
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;						   //
		if (!(*itTrk)->isMdcTrackValid())											   //
			continue;																   //
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();									   //
		double pch = mdcTrk->p();													   //track开始的动量
		double x0 = mdcTrk->x();													   //track开始的x坐标
		double y0 = mdcTrk->y();													   //track开始的y坐标
		double z0 = mdcTrk->z();													   //track开始的z坐标
		double phi0 = mdcTrk->helix(1);												   //表示螺旋线参数
																					   //0: d0 -> 螺旋线到对撞顶点的最小距离
																					   //1: phi0 -> 最小距离的xy平面相角
																					   //2: kappa
																					   //3: d
																					   //4: tan(lamda)
		double xv = xorigin.x();													   //对撞顶点的x坐标
		double yv = xorigin.y();													   //对撞顶点的y坐标
		double Rxy = (x0 - xv) * cos(phi0) + (y0 - yv) * sin(phi0);					   //计算顶点到track开始位置的xy平面距离
		HepVector a = mdcTrk->helix();												   //helix存入a
		HepSymMatrix Ea = mdcTrk->err();											   //helix-err存入Ea
		HepPoint3D point0(0., 0., 0.);												   //探测器原点
		HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);							   //对撞原点
		VFHelix helixip(point0, a, Ea);												   //组合数据{(0,0,0),helix(),err()}
		helixip.pivot(IP);															   //调整数据{(0,0,0),helix()-IP,err()-IP}
		HepVector vecipa = helixip.a();												   //提取新helix
		double Rvxy0 = fabs(vecipa[0]);												   //
		double Rvz0 = vecipa[3];													   //
		double Rvphi0 = vecipa[1];													   //
		if (fabs(cos(mdcTrk->theta())) > 0.93)										   //
			continue;																   //
		if (fabs(Rvz0) >= 10)														   //
			continue;																   //
		if (fabs(Rvxy0) >= 1)														   //
			continue;																   //
		iGood.push_back(i);															   //
		nCharge += mdcTrk->charge();												   //
	}																				   //
	int nGood = iGood.size();														   //ngood记录good trach数量，ncharge记录带电量
	log << MSG::DEBUG << "ngood, totcharge = " << nGood << " , " << nCharge << endreq; //
	if ((nGood != 2) || (nCharge != 0))												   //参数设定：带电总数，带电总量
	{																				   //
		return StatusCode::SUCCESS;													   //
	}																				   //
	Ncut1++;																		   //
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
		if ((emcTrk->module() == 1) && (fabs(ctht) > 0.8))												  //
			continue;																					  //
		if ((emcTrk->module() == 0 || emcTrk->module() == 2) && (fabs(ctht) > 0.92 || fabs(ctht) < 0.86)) //
			continue;																					  //
		if (emcTrk->time() < 0 || emcTrk->time() > 14)													  //
			continue;																					  //
		if (eraw < 0.025)																				  //
			continue;																					  //
		if (fabs(dang) < 10)																			  //
			continue;																					  //
		iGam.push_back(i);																				  //
	}																									  //
	int nGam = iGam.size();																				  //
	log << MSG::DEBUG << "num Good Photon " << nGam << " , " << evtRecEvent->totalNeutral() << endreq;	//
	if (nGam < 6)																						  //
	{																									  //
		return StatusCode::SUCCESS;																		  //
	}																									  //
	Ncut2++;																							  //
	//*********************************************************************************
	// Calculation 1: 4-momentum to each photon
	//*********************************************************************************
	Vp4 pGam;														 //变量：pGam[i],i对应光子新编号,值为光子动量,类型HepLorentzVector
	pGam.clear();													 //变量：iGam[i],i对应光子新编号,值为光子原始编号,类型Vint
	for (int i = 0; i < nGam; i++)									 //变量：nGam,值为iGam长度,类型int
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
		pid->init();																	   //对于Likelihood方法
		pid->setMethod(pid->methodProbability());										   //pid->setMethod(pid->methodLikelihood());
		pid->setChiMinCut(4);															   //
		pid->setRecTrack(*itTrk);														   //
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE()); //选择pid的粒子
		pid->identify(pid->onlyPion() | pid->onlyKaon());								   //pid->identify(pid->onlyPion());
		pid->calculate();																   //pid->identify(pid->onlyKaon());
		if (!(pid->IsPidInfoValid()))													   //
			continue;																	   //对于Likelihood方法(0=electron 1=muon 2=pion 3=kaon 4=proton)
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();										   //if(pid->pdf(2)<pid->pdf(3))
		if (pid->probPion() < 0.001)													   //continue;
			continue;																	   //
		RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();							   //对于ParticleID, 用RecMdcKalTrack代替RecMdcTrack
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);								   //PID可以设定为electron, muon, pion, kaon and proton;The default setting is pion
		if (mdcKalTrk->charge() > 0)													   //
		{																				   //
			ipip.push_back(iGood[i]);													   //
			HepLorentzVector ptrk;														   //
			ptrk.setPx(mdcKalTrk->px());												   //
			ptrk.setPy(mdcKalTrk->py());												   //
			ptrk.setPz(mdcKalTrk->pz());												   //
			double p3 = ptrk.mag();														   //
			ptrk.setE(sqrt(p3 * p3 + mpi * mpi));										   //
			ppip.push_back(ptrk);														   //变量：ppip,值为pi+动量,类型为HepLorentzVector
		}																				   //变量：ppim,值为pi-动量,类型为HepLorentzVector
		else																			   //
		{																				   //
			ipim.push_back(iGood[i]);													   //
			HepLorentzVector ptrk;														   //
			ptrk.setPx(mdcKalTrk->px());												   //
			ptrk.setPy(mdcKalTrk->py());												   //
			ptrk.setPz(mdcKalTrk->pz());												   //
			double p3 = ptrk.mag();														   //
			ptrk.setE(sqrt(p3 * p3 + mpi * mpi));										   //
			ppim.push_back(ptrk);														   //
		}																				   //
	}																					   //
	int npip = ipip.size();																   //
	int npim = ipim.size();																   //
	if (npip * npim != 1)																   //
		return SUCCESS;																	   //
	Ncut3++;																			   //
	cout << "pass pid" << endl;															   //
																						   //*********************************************************************************
																						   // Selection 3: Vertex fit Selection, check ppi0, pTot
																						   //*********************************************************************************
	RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin() + ipip[0]))->mdcKalTrack();		   //Default is pion, for other particles:
	RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin() + ipim[0]))->mdcKalTrack();		   //wvppTrk = WTrackParameter(mp, pipTrk->getZHelixP(), pipTrk->getZErrorP()); proton
	WTrackParameter wvpipTrk, wvpimTrk;													   //wvmupTrk = WTrackParameter(mmu, pipTrk->getZHelixMu(), pipTrk->getZErrorMu()); muon
	wvpipTrk = WTrackParameter(mpi, pipTrk->getZHelix(), pipTrk->getZError());			   //wvepTrk = WTrackParameter(me, pipTrk->getZHelixE(), pipTrk->getZErrorE()); electron
	wvpimTrk = WTrackParameter(mpi, pimTrk->getZHelix(), pimTrk->getZError());			   //wvkpTrk = WTrackParameter(mk, pipTrk->getZHelixK(), pipTrk->getZErrorK()); kaon
	HepPoint3D vx(0., 0., 0.);															   //
	HepSymMatrix Evx(3, 0);																   //
	double bx = 1E+6;																	   //
	double by = 1E+6;																	   //
	double bz = 1E+6;																	   //
	Evx[0][0] = bx * bx;																   //
	Evx[1][1] = by * by;																   //
	Evx[2][2] = bz * bz;																   //
	VertexParameter vxpar;																   //
	vxpar.setVx(vx);																	   //
	vxpar.setEvx(Evx);																	   //
	VertexFit *vtxfit = VertexFit::instance();											   //
	vtxfit->init();																		   //
	vtxfit->AddTrack(0, wvpipTrk);														   //设定track0
	vtxfit->AddTrack(1, wvpimTrk);														   //设定track1
	vtxfit->AddVertex(0, vxpar, 0, 1);													   //设定顶点0
	if (!vtxfit->Fit(0))																   //
		return SUCCESS;																	   //
	vtxfit->Swim(0);																	   //
	Ncut6++;																			   //
	cout << "pass VFit" << endl;														   //
	//*********************************************************************************
	// Selection 7: 4~5C Selection
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
	int select[3][3] = {{0, 1, 2},			   //
						{1, 0, 2},			   //
						{2, 0, 1}};			   //
	//*********************************************************************************
	// Selection 7-1: 4C Selection
	//*********************************************************************************
	cout << "before 4c" << endl;										 //准备4C声明
	WTrackParameter wpip = vtxfit->wtrk(0);								 //
	WTrackParameter wpim = vtxfit->wtrk(1);								 //
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();			 //
	if (m_test4C == 1)													 //
	{																	 //
		HepLorentzVector ecms(0.034 * m_energy / 3.097, 0, 0, m_energy); //
		double chisq_4c = 9999;											 //
		HepLorentzVector ptrack0;										 //记录ptrack_fit0~7
		HepLorentzVector ptrack1;										 //
		HepLorentzVector ptrack2;										 //
		HepLorentzVector ptrack3;										 //
		HepLorentzVector ptrack4;										 //
		HepLorentzVector ptrack5;										 //
		HepLorentzVector ptrack6;										 //
		HepLorentzVector ptrack7;										 //
		int igamma0 = -1;												 //记录igamma0~5
		int igamma1 = -1;												 //
		int igamma2 = -1;												 //
		int igamma3 = -1;												 //
		int igamma4 = -1;												 //
		int igamma5 = -1;												 //
		for (int i1 = 0; i1 < nGam - 5; i1++)
		{
			RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();
			for (int i2 = i1 + 1; i2 < nGam - 4; i2++)
			{
				RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();
				for (int i3 = i2 + 1; i3 < nGam - 3; i3++)
				{
					RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[i3]))->emcShower();
					for (int i4 = i3 + 1; i4 < nGam - 2; i4++)
					{
						RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[i4]))->emcShower();
						for (int i5 = i4 + 1; i5 < nGam - 1; i5++)
						{
							RecEmcShower *g5Trk = (*(evtRecTrkCol->begin() + iGam[i5]))->emcShower();
							for (int i6 = i5 + 1; i6 < nGam - 0; i6++)
							{
								RecEmcShower *g6Trk = (*(evtRecTrkCol->begin() + iGam[i6]))->emcShower();
								kmfit->init();					  //
								kmfit->AddTrack(0, wpip);		  //
								kmfit->AddTrack(1, wpim);		  //
								kmfit->AddTrack(2, 0.0, g1Trk);   //
								kmfit->AddTrack(3, 0.0, g2Trk);   //
								kmfit->AddTrack(4, 0.0, g3Trk);   //
								kmfit->AddTrack(5, 0.0, g4Trk);   //
								kmfit->AddTrack(6, 0.0, g5Trk);   //
								kmfit->AddTrack(7, 0.0, g6Trk);   //
								kmfit->AddFourMomentum(0, ecms);  //
								bool oksq = kmfit->Fit();		  //
								if (oksq)						  //
								{								  //
									double chi2 = kmfit->chisq(); //
									if (chi2 < chisq_4c)		  //
									{							  //
										chisq_4c = chi2;		  //
										ptrack0 = kmfit->pfit(0); //
										ptrack1 = kmfit->pfit(1); //
										ptrack2 = kmfit->pfit(2); //
										ptrack3 = kmfit->pfit(3); //
										ptrack4 = kmfit->pfit(4); //
										ptrack5 = kmfit->pfit(5); //
										ptrack6 = kmfit->pfit(6); //
										ptrack7 = kmfit->pfit(7); //
										igamma0 = i1;			  //
										igamma1 = i2;			  //
										igamma2 = i3;			  //
										igamma3 = i4;			  //
										igamma4 = i5;			  //
										igamma5 = i6;			  //
									}							  //
								}								  //
							}									  //
						}										  //
					}											  //
				}												  //
			}													  //
		}														  //
		if (chisq_4c < 200)										  //
		{														  //
			double chisq = 9999;								  //
			double mb0;											  //
			double momega;										  //
			double mpi01;										  //
			double mpi02;										  //
			double mpi03;										  //
			HepLorentzVector ptrack[7] = {ptrack0,				  //
										  ptrack2,				  //
										  ptrack3,				  //
										  ptrack4,				  //
										  ptrack5,				  //
										  ptrack6,				  //
										  ptrack7};				  //
			double ib0;											  //
			double iomega;										  //
			double ipi01;										  //
			double ipi02;										  //
			double ipi03;										  //
			double chisqb;										  //
			double chisqo;										  //
			double chisq1;										  //
			double chisq2;										  //
			double chisq3;										  //
			for (int i = 0; i < 15; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					ib0 = (ppip[0] +
						   ppim[0] +
						   ptrack[combine[i][2 * select[j][0]]] +
						   ptrack[combine[i][2 * select[j][0] + 1]] +
						   ptrack[combine[i][2 * select[j][1]]] +
						   ptrack[combine[i][2 * select[j][1] + 1]])
							  .m();
					iomega = (ppip[0] +
							  ppim[0] +
							  ptrack[combine[i][2 * select[j][0]]] +
							  ptrack[combine[i][2 * select[j][0] + 1]])
								 .m();
					ipi01 = (ptrack[combine[i][2 * select[j][0]]] +
							 ptrack[combine[i][2 * select[j][0] + 1]])
								.m();
					ipi02 = (ptrack[combine[i][2 * select[j][1]]] +
							 ptrack[combine[i][2 * select[j][1] + 1]])
								.m();
					ipi03 = (ptrack[combine[i][2 * select[j][2]]] +
							 ptrack[combine[i][2 * select[j][2] + 1]])
								.m();
					chisqb = pow((ib0 - 1.235), 2);
					chisqo = pow((iomega - 0.782), 2);
					chisq1 = pow((ipi01 - 0.135), 2);
					chisq2 = pow((ipi02 - 0.135), 2);
					chisq3 = pow((ipi03 - 0.135), 2);
					double chi2 = chisqb / 10 + chisqo + (chisq1 + chisq2 + chisq3) / 3;
					if (chi2 < chisq)
					{
						chisq = chi2;
						mb0 = ib0;
						momega = iomega;
						mpi01 = ipi01;
						mpi02 = ipi02;
						mpi03 = ipi03;
					}
				}
			}
			if (1 == 1)
			{
				m_chisq_4c = chisq;
				m_b0_4c = mb0;
				m_omega_4c = momega;
				m_pi01_4c = mpi01;
				m_pi02_4c = mpi02;
				m_pi03_4c = mpi03;
				m_tuple4->write();
				Ncut4++;
			}
		}
	}
	//*********************************************************************************
	// Selection 7-2: 5C Selection
	//*********************************************************************************
	cout << "before 5c" << endl; //准备5C声明
	if (m_test5C == 1)
	{
		HepLorentzVector ecms(0.034 * m_energy / 3.097, 0, 0, m_energy);								  //
		double chisq = 9999;																			  //
		int ig1 = -1;																					  //
		int ig2 = -1;																					  //
		int ig3 = -1;																					  //
		int ig4 = -1;																					  //
		int ig5 = -1;																					  //
		int ig6 = -1;																					  //
		int io[6] = {-1, -1, -1, -1, -1, -1};															  //
		for (int i1 = 0; i1 < nGam - 5; i1++)															  //
		{																								  //
			RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();					  //
			for (int i2 = i1 + 1; i2 < nGam - 4; i2++)													  //
			{																							  //
				RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();				  //
				for (int i3 = i2 + 1; i3 < nGam - 3; i3++)												  //
				{																						  //
					RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[i3]))->emcShower();			  //
					for (int i4 = i3 + 1; i4 < nGam - 2; i4++)											  //
					{																					  //
						RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[i4]))->emcShower();		  //
						for (int i5 = i4 + 1; i5 < nGam - 1; i5++)										  //
						{																				  //
							RecEmcShower *g5Trk = (*(evtRecTrkCol->begin() + iGam[i5]))->emcShower();	 //
							for (int i6 = i5 + 1; i6 < nGam - 0; i6++)									  //
							{																			  //
								RecEmcShower *g6Trk = (*(evtRecTrkCol->begin() + iGam[i6]))->emcShower(); //
								for (int i = 0; i < 15; i++)											  //
								{																		  //
									kmfit->init();														  //
									kmfit->AddTrack(0, wpip);											  //
									kmfit->AddTrack(1, wpim);											  //
									kmfit->AddTrack(2, 0.0, g1Trk);										  //
									kmfit->AddTrack(3, 0.0, g2Trk);										  //
									kmfit->AddTrack(4, 0.0, g3Trk);										  //
									kmfit->AddTrack(5, 0.0, g4Trk);										  //
									kmfit->AddTrack(6, 0.0, g5Trk);										  //
									kmfit->AddTrack(7, 0.0, g6Trk);										  //
									kmfit->AddResonance(0, 0.135, combine[i][0] + 1, combine[i][1] + 1);  //
									kmfit->AddResonance(1, 0.135, combine[i][2] + 1, combine[i][3] + 1);  //
									kmfit->AddResonance(2, 0.135, combine[i][4] + 1, combine[i][5] + 1);  //
									kmfit->AddFourMomentum(3, ecms);									  //
									bool oksq = kmfit->Fit();											  //
									if (oksq)															  //
									{																	  //
										double chi2 = kmfit->chisq();									  //
										if (chi2 < chisq)												  //
										{																  //
											chisq = chi2;												  //
											int ig[6] = {i1, i2, i3, i4, i5, i6};						  //
											io[0] = ig[combine[i][0] - 1];								  //
											io[1] = ig[combine[i][1] - 1];								  //
											io[2] = ig[combine[i][2] - 1];								  //
											io[3] = ig[combine[i][3] - 1];								  //
											io[4] = ig[combine[i][4] - 1];								  //
											io[5] = ig[combine[i][5] - 1];								  //
										}																  //
									}																	  //
								}																		  //
							}																			  //
						}																				  //
					}																					  //
				}																						  //
			}																							  //
		}																								  //
		log << MSG::INFO << " chisq = " << chisq << endreq;												  //
		if (chisq < 200)																				  //
		{																								  //
			double chisq_o = 9999;																		  //
			double mb0 = -1;																			  //
			double momega = -1;																			  //
			double mpi01 = -1;																			  //
			double mpi02 = -1;																			  //
			double mpi03 = -1;																			  //
			RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[io[0]]))->emcShower();				  //
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[io[1]]))->emcShower();				  //
			RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[io[2]]))->emcShower();				  //
			RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[io[3]]))->emcShower();				  //
			RecEmcShower *g5Trk = (*(evtRecTrkCol->begin() + iGam[io[4]]))->emcShower();				  //
			RecEmcShower *g6Trk = (*(evtRecTrkCol->begin() + iGam[io[5]]))->emcShower();				  //
			kmfit->init();																				  //
			kmfit->AddTrack(0, wpip);																	  //
			kmfit->AddTrack(1, wpim);																	  //
			kmfit->AddTrack(2, 0.0, g1Trk);																  //
			kmfit->AddTrack(3, 0.0, g2Trk);																  //
			kmfit->AddTrack(4, 0.0, g3Trk);																  //
			kmfit->AddTrack(5, 0.0, g4Trk);																  //
			kmfit->AddTrack(6, 0.0, g5Trk);																  //
			kmfit->AddTrack(7, 0.0, g6Trk);																  //
			kmfit->AddResonance(0, 0.135, 2, 3);														  //
			kmfit->AddResonance(1, 0.135, 4, 5);														  //
			kmfit->AddResonance(2, 0.135, 6, 7);														  //
			kmfit->AddFourMomentum(3, ecms);															  //
			bool oksq = kmfit->Fit();																	  //
			if (oksq)																					  //
			{																							  //
				for (int j = 0; j < 3; j++)																  //
				{																						  //
					double chi2_1 = pow((kmfit->pfit(0) +
										 kmfit->pfit(1) +
										 kmfit->pfit(2 * select[j][0] + 2) +
										 kmfit->pfit(2 * select[j][0] + 3))
												.m() -
											1.235,
										2);
					double chi2_2 = pow((kmfit->pfit(0) +
										 kmfit->pfit(1) +
										 kmfit->pfit(2 * select[j][0] + 2) +
										 kmfit->pfit(2 * select[j][0] + 3) +
										 kmfit->pfit(2 * select[j][1] + 2) +
										 kmfit->pfit(2 * select[j][1] + 3))
												.m() -
											0.782,
										2);
					double chi2 = chi2_1 + chi2_2 / 100;
					if (chi2 < chisq_o)
					{
						chisq_o = chi2;
						mb0 = (kmfit->pfit(0) +
							   kmfit->pfit(1) +
							   kmfit->pfit(2 * select[j][0] + 2) +
							   kmfit->pfit(2 * select[j][0] + 3))
								  .m();
						momega = (kmfit->pfit(0) +
								  kmfit->pfit(1) +
								  kmfit->pfit(2 * select[j][0] + 2) +
								  kmfit->pfit(2 * select[j][0] + 3))
									 .m();
						mpi01 = (kmfit->pfit(2 * select[j][0] + 2) +
								 kmfit->pfit(2 * select[j][0] + 3))
									.m();
						mpi02 = (kmfit->pfit(2 * select[j][1] + 2) +
								 kmfit->pfit(2 * select[j][1] + 3))
									.m();
						mpi03 = (kmfit->pfit(2 * select[j][2] + 2) +
								 kmfit->pfit(2 * select[j][2] + 3))
									.m();
					}
				}
				if (1 == 1)
				{
					m_chisq_5c = kmfit->chisq();
					m_b0_5c = mb0;
					m_momega_5c = momega;
					m_mpi01_5c = mpi01;
					m_mpi02_5c = mpi02;
					m_mpi03_5c = mpi03;
					m_tuple5->write();
					Ncut5++;
				}
			}
		}
	}
	return StatusCode::SUCCESS;
}
//*********************************************************************************************************
//***                                               finalize                                            ***
//*********************************************************************************************************
StatusCode Omega::finalize()
{
	cout << "total number:         " << Ncut0 << endl;
	cout << "nGood==2, nCharge==0: " << Ncut1 << endl;
	cout << "nGam>=6:              " << Ncut2 << endl;
	cout << "Pass Pid:             " << Ncut3 << endl;
	cout << "Pass VFit:            " << Ncut6 << endl;
	cout << "Pass 4C:              " << Ncut4 << endl;
	cout << "Pass 5C:              " << Ncut5 << endl;
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}
