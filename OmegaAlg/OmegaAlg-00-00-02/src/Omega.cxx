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
int Ncut0, Ncut3, Ncut4, Ncut5;
//*********************************************************************************************************
//***                                                声明容器                                            ***
//*********************************************************************************************************
Omega::Omega(const std::string &name, ISvcLocator *pSvcLocator) : Algorithm(name, pSvcLocator)
{
	declareProperty("Test4C", m_test4C = 1);
	declareProperty("Test5C", m_test5C = 0);
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
	// initialize-data-in-fit4c
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
				status = m_tuple4->addItem("chisq", m_chisq_4c);
				status = m_tuple4->addItem("chisq_pi", m_chisq_3pi);
				status = m_tuple4->addItem("mpi01", m_pi01);
				status = m_tuple4->addItem("mpi02", m_pi02);
				status = m_tuple4->addItem("mpi03", m_pi03);
				status = m_tuple4->addItem("momega", m_omega);
				status = m_tuple4->addItem("momegapi02", m_omegapi02);
				status = m_tuple4->addItem("momegapi03", m_omegapi03);
				status = m_tuple4->addItem("mpi02pi03", m_pi02pi03);

				status = m_tuple4->addItem("runID", runID);
				status = m_tuple4->addItem("eventID", eventID);
				status = m_tuple4->addItem("indexmc", m_idxmc, 0, 100);
				status = m_tuple4->addIndexedItem("pdgid", m_idxmc, m_pdgid);
				status = m_tuple4->addIndexedItem("motheridx", m_idxmc, m_motheridx);
			}
			else
			{
				log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple4) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
	//
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
	runID = runNo;																			   // 变量：topo
	eventID = event;																		   //
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
	}																							 //
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
		if ((emcTrk->module() == 1) && (fabs(ctht) > 0.8))												  //
			continue;																					  //
		if ((emcTrk->module() == 0 || emcTrk->module() == 2) && (fabs(ctht) > 0.92 || fabs(ctht) < 0.86)) // 选择：E-endcap
			continue;																					  //
		if (emcTrk->time() < 0 || emcTrk->time() > 14)													  // 选择：TDC
			continue;																					  //
		if (eraw < 0.025)																				  // 选择：E-barrel
			continue;																					  //
		if (fabs(dang) < 10)																			  // 选择：
			continue;																					  //
		iGam.push_back(i);																				  // 变量：iGam[]（参数为good-track序号，内容为track编号）
	}																									  //
	int nGam = iGam.size();																				  // 变量：nGam（中性track数量）
	log << MSG::DEBUG << "num Good Photon " << nGam << " , " << evtRecEvent->totalNeutral() << endreq;	//
	if (nGam < 6)																						  // 选择：nGam
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
	cout << "pass pid" << endl;
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
	cout << "pass VFit" << endl;												  //
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
	if (m_test4C == 1)																								  //
	{																												  //
		HepLorentzVector ecms(0.034 * m_energy / 3.097, 0, 0, m_energy);											  //
		HepLorentzVector ptrackp;																					  //
		HepLorentzVector ptrackm;																					  //
		HepLorentzVector ptrack1, nptrack1;																			  //
		HepLorentzVector ptrack2, nptrack2;																			  //
		HepLorentzVector ptrack3, nptrack3;																			  //
		HepLorentzVector ptrack4, nptrack4;																			  //
		HepLorentzVector ptrack5, nptrack5;																			  //
		HepLorentzVector ptrack6, nptrack6;																			  //
		double chisq_4c_5g = 9999;																					  //
		double chisq_4c_6g = 9999;																					  //
		double chisq_4c_7g = 9999;																					  //
		double chisq_4c_pi = 9999;																					  //
		double chisq_4c_om = 9999;																					  //
		double chisq_4c_b0 = 9999;																					  //
		for (int i1 = 0; i1 < nGam; i1++)																			  // 得到6Gamma-chisq
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
								kmfit->init();																		  //
								kmfit->AddTrack(0, wpip);															  //
								kmfit->AddTrack(1, wpim);															  //
								kmfit->AddTrack(2, 0.0, g1Trk);														  //
								kmfit->AddTrack(3, 0.0, g2Trk);														  //
								kmfit->AddTrack(4, 0.0, g3Trk);														  //
								kmfit->AddTrack(5, 0.0, g4Trk);														  //
								kmfit->AddTrack(6, 0.0, g5Trk);														  //
								kmfit->AddTrack(7, 0.0, g6Trk);														  //
								kmfit->AddFourMomentum(0, ecms);													  //
								bool oksq = kmfit->Fit();															  //
								if (oksq)																			  //
								{																					  //
									double chi2 = kmfit->chisq();													  //
									if (chi2 <= chisq_4c_6g)														  // 选择：最小chi-4c
									{																				  //
										chisq_4c_6g = chi2;															  //
										ptrackp = kmfit->pfit(0);													  //
										ptrackm = kmfit->pfit(1);													  //
										ptrack1 = kmfit->pfit(2);													  //
										ptrack2 = kmfit->pfit(3);													  //
										ptrack3 = kmfit->pfit(4);													  //
										ptrack4 = kmfit->pfit(5);													  //
										ptrack5 = kmfit->pfit(6);													  //
										ptrack6 = kmfit->pfit(7);													  //
									}																				  //
								}																					  //
							}																						  //
						}																							  //
					}																								  //
				}																									  //
			}																										  //
		}																											  //
		for (int i1 = 0; i1 < nGam; i1++)																			  // 得到5Gamma-chisq
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
							kmfit->init();																			  //
							kmfit->AddTrack(0, wpip);																  //
							kmfit->AddTrack(1, wpim);																  //
							kmfit->AddTrack(2, 0.0, g1Trk);															  //
							kmfit->AddTrack(3, 0.0, g2Trk);															  //
							kmfit->AddTrack(4, 0.0, g3Trk);															  //
							kmfit->AddTrack(5, 0.0, g4Trk);															  //
							kmfit->AddTrack(6, 0.0, g5Trk);															  //
							kmfit->AddFourMomentum(0, ecms);														  //
							bool oksq = kmfit->Fit();																  //
							if (oksq)																				  //
							{																						  //
								double chi2 = kmfit->chisq();														  //
								if (chi2 < chisq_4c_5g)																  //
								{																					  //
									chisq_4c_5g = chi2;																  //
								}																					  //
							}																						  //
						}																							  //
					}																								  //
				}																									  //
			}																										  //
		}																											  //
		if (nGam > 6)																								  // 得到7Gamma-chisq
		{																											  //
			for (int i1 = 0; i1 < nGam; i1++)																		  //
			{																										  //
				RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();							  //
				for (int i2 = i1; i2 < nGam; i2++)																	  //
				{																									  //
					if (i2 == i1)																					  //
					{																								  //
						continue;																					  //
					}																								  //
					RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();						  //
					for (int i3 = i2; i3 < nGam; i3++)																  //
					{																								  //
						if (i3 == i1 || i3 == i2)																	  //
						{																							  //
							continue;																				  //
						}																							  //
						RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[i3]))->emcShower();					  //
						for (int i4 = i3; i4 < nGam; i4++)															  //
						{																							  //
							if (i4 == i1 || i4 == i2 || i4 == i3)													  //
							{																						  //
								continue;																			  //
							}																						  //
							RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[i4]))->emcShower();				  //
							for (int i5 = i4; i5 < nGam; i5++)														  //
							{																						  //
								if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4)									  //
								{																					  //
									continue;																		  //
								}																					  //
								RecEmcShower *g5Trk = (*(evtRecTrkCol->begin() + iGam[i5]))->emcShower();			  //
								for (int i6 = i5; i6 < nGam; i6++)													  //
								{																					  //
									if (i6 == i1 || i6 == i2 || i6 == i3 || i6 == i4 || i6 == i5)					  //
									{																				  //
										continue;																	  //
									}																				  //
									RecEmcShower *g6Trk = (*(evtRecTrkCol->begin() + iGam[i6]))->emcShower();		  //
									for (int i7 = i6; i7 < nGam; i7++)												  //
									{																				  //
										if (i7 == i1 || i7 == i2 || i7 == i3 || i7 == i4 || i7 == i5 || i7 == i6)	 //
										{																			  //
											continue;																  //
										}																			  //
										RecEmcShower *g7Trk = (*(evtRecTrkCol->begin() + iGam[i7]))->emcShower();	 //
										kmfit->init();																  //
										kmfit->AddTrack(0, wpip);													  //
										kmfit->AddTrack(1, wpim);													  //
										kmfit->AddTrack(2, 0.0, g1Trk);												  //
										kmfit->AddTrack(3, 0.0, g2Trk);												  //
										kmfit->AddTrack(4, 0.0, g3Trk);												  //
										kmfit->AddTrack(5, 0.0, g4Trk);												  //
										kmfit->AddTrack(6, 0.0, g5Trk);												  //
										kmfit->AddTrack(7, 0.0, g6Trk);												  //
										kmfit->AddTrack(8, 0.0, g7Trk);												  //
										kmfit->AddFourMomentum(0, ecms);											  //
										bool oksq = kmfit->Fit();													  //
										if (oksq)																	  //
										{																			  //
											double chi2 = kmfit->chisq();											  //
											if (chi2 < chisq_4c_7g)													  //
											{																		  //
												chisq_4c_7g = chi2;													  //
											}																		  //
										}																			  //
									}																				  //
								}																					  //
							}																						  //
						}																							  //
					}																								  //
				}																									  //
			}																										  //
		}																											  //
		int g6 = 0;																									  //
		if (chisq_4c_6g < 200)																						  // 选择：chisq-6Gamma
		{																											  //
			g6 = 1;																									  //
		}																											  //
		if (chisq_4c_5g < chisq_4c_6g)																				  // 选择：chisq-5Gamma
		{																											  //
			g6 = 0;																									  //
		}																											  //
		if (nGam > 6)																								  // 选择：chisq-7Gamma
		{																											  //
			if (chisq_4c_7g < chisq_4c_6g)																			  //
			{																										  //
				g6 = 0;																								  //
			}																										  //
		}																											  //
		if (g6 == 1)																								  //
		{																											  //
			if (1 == 1)																								  //
			{																										  //
				HepLorentzVector ptrack[6] = {ptrack1,																  //
											  ptrack2,																  //
											  ptrack3,																  //
											  ptrack4,																  //
											  ptrack5,																  //
											  ptrack6};																  //
				for (int i = 0; i < 15; i++)																		  //
				{																									  //
					double chisq_mpi01 = pow((ptrack[combine[i][0] - 1] + ptrack[combine[i][0] - 1]).m() - 0.135, 2); //
					double chisq_mpi02 = pow((ptrack[combine[i][2] - 1] + ptrack[combine[i][3] - 1]).m() - 0.135, 2); //
					double chisq_mpi03 = pow((ptrack[combine[i][4] - 1] + ptrack[combine[i][5] - 1]).m() - 0.135, 2); //
					double chisq_mpi0 = (chisq_mpi01 + chisq_mpi02 + chisq_mpi03) / 3;								  //
					if (chisq_mpi0 < chisq_4c_pi)																	  //
					{																								  //
						chisq_4c_pi = chisq_mpi0;																	  //
						nptrack1 = ptrack[combine[i][0] - 1];														  //
						nptrack2 = ptrack[combine[i][1] - 1];														  //
						nptrack3 = ptrack[combine[i][2] - 1];														  //
						nptrack4 = ptrack[combine[i][3] - 1];														  //
						nptrack5 = ptrack[combine[i][4] - 1];														  //
						nptrack6 = ptrack[combine[i][5] - 1];														  //
					}																								  //
				}																									  //
			}																										  //
			if (1 == 1)
			{
				HepLorentzVector ptrack[6] = {nptrack1,  //
											  nptrack2,  //
											  nptrack3,  //
											  nptrack4,  //
											  nptrack5,  //
											  nptrack6}; //
				for (int i = 0; i < 3; i++)
				{
					double chisq_momega = pow((ptrackp + ptrackm + ptrack[2 * selecto[i][0] - 2] + ptrack[2 * selecto[i][0] - 1]).m() - 0.782, 2);
					if (chisq_momega < chisq_4c_om)
					{
						chisq_4c_om = chisq_momega;
						nptrack1 = ptrack[2 * selecto[i][0] - 2];
						nptrack2 = ptrack[2 * selecto[i][0] - 1];
						nptrack3 = ptrack[2 * selecto[i][1] - 2];
						nptrack4 = ptrack[2 * selecto[i][1] - 1];
						nptrack5 = ptrack[2 * selecto[i][2] - 2];
						nptrack6 = ptrack[2 * selecto[i][2] - 1];
					}
				}
			}
			if (1 == 0)
			{
				HepLorentzVector ptrack[6] = {nptrack1,  //
											  nptrack2,  //
											  nptrack3,  //
											  nptrack4,  //
											  nptrack5,  //
											  nptrack6}; //
				for (int i = 0; i < 2; i++)
				{
					double chisq_mb = pow((ptrack[2 * selectb[i][0]] + ptrack[2 * selectb[i][0] + 1]).m(), 2);
					if (chisq_mb < chisq_4c_b0)
					{
						chisq_4c_b0 = chisq_mb;
						nptrack1 = ptrack[0];
						nptrack2 = ptrack[1];
						nptrack3 = ptrack[2 * selectb[i][0] + 0];
						nptrack4 = ptrack[2 * selectb[i][0] + 1];
						nptrack5 = ptrack[2 * selectb[i][1] + 0];
						nptrack6 = ptrack[2 * selectb[i][1] + 1];
					}
				}
			}
			HepLorentzVector out_pi01 = nptrack1 + nptrack2;
			HepLorentzVector out_pi02 = nptrack3 + nptrack4;
			HepLorentzVector out_pi03 = nptrack5 + nptrack6;
			HepLorentzVector out_omega = ptrackp + ptrackm + out_pi01;
			HepLorentzVector out_omegapi02 = out_omega + out_pi02;
			HepLorentzVector out_omegapi03 = out_omega + out_pi03;
			HepLorentzVector out_pi02pi03 = out_pi02 + out_pi03;
			double out_mpi01 = out_pi01.m();
			double out_mpi02 = out_pi02.m();
			double out_mpi03 = out_pi03.m();
			double out_momega = out_omega.m();
			double out_momegapi02 = out_omegapi02.m();
			double out_momegapi03 = out_omegapi03.m();
			double out_mpi02pi03 = out_pi02pi03.m();
			if (1 == 1)
			{
				m_chisq_4c = chisq_4c_6g;
				m_chisq_3pi = chisq_4c_pi;
				m_pi01 = out_mpi01;
				m_pi02 = out_mpi02;
				m_pi03 = out_mpi03;
				m_omega = out_momega;
				m_omegapi02 = out_momegapi02;
				m_omegapi03 = out_momegapi03;
				m_pi02pi03 = out_mpi02pi03;
				m_tuple4->write();
				Ncut4++;
			}
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
	cout << "Pass Pid:             " << Ncut3 << endl;
	cout << "Pass 4C:              " << Ncut4 << endl;
	cout << "Pass 5C:              " << Ncut5 << endl;
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}
