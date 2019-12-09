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
int Ncut0, Ncut1, Ncut3, Ncut4, Ncut5, number, checkexit, checki, firstrun;
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
				status = m_tuple4->addItem("chisq_pi", m_chisq_3pi);
				status = m_tuple4->addItem("mpi01", m_pi01);
				status = m_tuple4->addItem("mpi02", m_pi02);
				status = m_tuple4->addItem("mpi03", m_pi03);
				status = m_tuple4->addItem("mPm", m_Pm);
				status = m_tuple4->addItem("mPmpi02", m_Pmpi02);
				status = m_tuple4->addItem("mPmpi03", m_Pmpi03);
				status = m_tuple4->addItem("mpi02pi03", m_pi02pi03);
				status = m_tuple4->addItem("mpi01pi02", m_pi01pi02);
				status = m_tuple4->addItem("mpi01pi03", m_pi01pi03);
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
	RecMdcKalTrack *pipTrk1 = (*(evtRecTrkCol->begin() + ipip[0]))->mdcKalTrack(); //
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
	HepPoint3D vx(0., 0., 0.);													   //
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
	VertexFit *vtxfit = VertexFit::instance();									   //
	vtxfit->init();																   //
	vtxfit->AddTrack(0, wvpipTrk);												   // 设定track0
	vtxfit->AddTrack(1, wvpimTrk);												   // 设定track1
	vtxfit->AddVertex(0, vxpar, 0, 1);											   // 设定顶点0
	if (!vtxfit->Fit(0))														   //
		return SUCCESS;															   //
	vtxfit->Swim(0);															   //
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
	WTrackParameter wpip = vtxfit->wtrk(0);												  //
	WTrackParameter wpim = vtxfit->wtrk(1);												  //
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();							  //
	if (m_test4C == 1)																	  //
	{																					  //
		HepLorentzVector ecms(0.034 * m_energy / 3.097, 0, 0, m_energy);				  //
		HepLorentzVector ptrackp;														  //
		HepLorentzVector ptrackm;														  //
		HepLorentzVector ptrack1, n1ptrack1, n2ptrack1, n3ptrack1;						  //
		HepLorentzVector ptrack2, n1ptrack2, n2ptrack2, n3ptrack2;						  //
		double chisq_4c_5g = 9999;														  //
		double chisq_4c_6g = 9999;														  //
		double chisq_4c_7g = 9999;														  //
		double chisq_4c_pi = 9999;														  //
		double chisq_4c_om = 9999;														  //
		for (int i1 = 0; i1 < nGam; i1++)												  // 得到6Gamma-chisq
		{																				  //
			RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();	 //
			for (int i2 = i1; i2 < nGam; i2++)											  //
			{																			  //
				if (i2 == i1)															  //
				{																		  //
					continue;															  //
				}																		  //
				RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower(); //
																						  //
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
	//
	//
	//
	//
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
				double chisq_mPm = pow((ptrackp + ptrackm + ptrackg2[2 * selecto[i][0] - 2] + ptrackg2[2 * selecto[i][0] - 1]).m() - 0.782, 2);
				if (chisq_mPm < chisq_4c_om)
				{
					chisq_4c_om = chisq_mPm;
					n2ptrack1 = ptrackg2[2 * selecto[i][0] - 2];
					n2ptrack2 = ptrackg2[2 * selecto[i][0] - 1];
					n2ptrack3 = ptrackg2[2 * selecto[i][1] - 2];
					n2ptrack4 = ptrackg2[2 * selecto[i][1] - 1];
					n2ptrack5 = ptrackg2[2 * selecto[i][2] - 2];
					n2ptrack6 = ptrackg2[2 * selecto[i][2] - 1];
				}
			}
		}
		if (1 == 0)
		{
			HepLorentzVector ptrackg3[6] = {n2ptrack1,  //
											n2ptrack2,  //
											n2ptrack3,  //
											n2ptrack4,  //
											n2ptrack5,  //
											n2ptrack6}; //
			for (int i = 0; i < 2; i++)
			{
				double chisq_mb = pow((ptrackg3[2 * selectb[i][0]] + ptrackg3[2 * selectb[i][0] + 1]).m(), 2);
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
		HepLorentzVector out_pi01 = n2ptrack1 + n2ptrack2;
		HepLorentzVector out_pi02 = n2ptrack3 + n2ptrack4;
		HepLorentzVector out_pi03 = n2ptrack5 + n2ptrack6;
		HepLorentzVector out_Pm = ptrackp + ptrackm + out_pi01;
		HepLorentzVector out_Pmpi02 = out_Pm + out_pi02;
		HepLorentzVector out_Pmpi03 = out_Pm + out_pi03;
		HepLorentzVector out_pi02pi03 = out_pi02 + out_pi03;
		HepLorentzVector out_pi01pi02 = out_pi01 + out_pi02;
		HepLorentzVector out_pi01pi03 = out_pi01 + out_pi03;
		double out_mpi01 = out_pi01.m();
		double out_mpi02 = out_pi02.m();
		double out_mpi03 = out_pi03.m();
		double out_mPm = out_Pm.m();
		double out_mPmpi02 = out_Pmpi02.m();
		double out_mPmpi03 = out_Pmpi03.m();
		double out_mpi02pi03 = out_pi02pi03.m();
		double out_mpi01pi02 = out_pi01pi02.m();
		double out_mpi01pi03 = out_pi01pi03.m();
		if (1 == 1)
		{
			m_chisq_4c = chisq_4c_6g;
			m_chisq_3pi = chisq_4c_pi;
			m_pi01 = out_mpi01;
			m_pi02 = out_mpi02;
			m_pi03 = out_mpi03;
			m_Pm = out_mPm;
			m_Pmpi02 = out_mPmpi02;
			m_Pmpi03 = out_mPmpi03;
			m_pi02pi03 = out_mpi02pi03;
			m_pi01pi02 = out_mpi02pi03;
			m_pi01pi03 = out_mpi02pi03;
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
StatusCode Pm::finalize()
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
