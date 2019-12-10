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
	declareProperty("Test4C", m_test4C = 1);
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
				status = m_tuple4->addItem("mpi0", m_pi0);
				status = m_tuple4->addItem("momega", m_omega);
				status = m_tuple4->addItem("momegapip", m_omegapip);
				status = m_tuple4->addItem("momegapim", m_omegapim);
				status = m_tuple4->addItem("mpipm", m_pipm);
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
	Ncut1++;
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
		if (emcTrk->time() < 0 || emcTrk->time() > 15)													  // 选择：TDC
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
	if (nGam < 2 || nGam > 55)																			  // 选择：nGam
	{																									  //
		return StatusCode::SUCCESS;																		  //
	}																									  //
	Ncut2++;
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
	RecMdcKalTrack *pipTrk1 = (*(evtRecTrkCol->begin() + ipip[1]))->mdcKalTrack(); //
	RecMdcKalTrack *pimTrk2 = (*(evtRecTrkCol->begin() + ipim[0]))->mdcKalTrack(); //
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
	Ncut4++;
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
	WTrackParameter wpip1 = vtxfit->wtrk(0);																							 //
	WTrackParameter wpip2 = vtxfit->wtrk(1);																							 //
	WTrackParameter wpim1 = vtxfit->wtrk(2);																							 //
	WTrackParameter wpim2 = vtxfit->wtrk(3);																							 //
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();																			 //
	if (m_test4C == 1)																													 //
	{																																	 //
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
						ptrackp2 = kmfit->pfit(0);																						 //
						ptrackm1 = kmfit->pfit(0);																						 //
						ptrackm2 = kmfit->pfit(0);																						 //
						ptrackn1 = kmfit->pfit(0);																						 //
						ptrackn2 = kmfit->pfit(0);																						 //
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
			HepLorentzVector out_omega = n1ptrackp1 + n1ptrackm1 + out_pi0;																 //
			HepLorentzVector out_omegapip = out_omega + n1ptrackp2;																		 //
			HepLorentzVector out_omegapim = out_omega + n1ptrackm2;																		 //
			HepLorentzVector out_pipm = n1ptrackp2 + n1ptrackm2;																		 //
			double out_mpi0 = out_pi0.m();																								 //
			double out_momega = out_omega.m();																							 //
			double out_momegapip = out_omegapip.m();																					 //
			double out_momegapim = out_omegapim.m();																					 //
			double out_mpipm = out_pipm.m();																							 //
			if (1 == 1)																													 //
			{																															 //
				m_chisq_4c = chisq_4c_2g;																								 //
				m_pi0 = out_mpi0;																										 //
				m_omega = out_momega;																									 //
				m_omegapip = out_momegapip;																								 //
				m_omegapim = out_momegapim;																								 //
				m_pipm = out_mpipm;																										 //
				m_tuple4->write();																										 //
				Ncut5++;																												 //
			}																															 //
		}																																 //
	}																																	 //
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
