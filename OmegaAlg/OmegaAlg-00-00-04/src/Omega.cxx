//*********************************************************************************************************
//***                                                代码引用                                            ***
//*********************************************************************************************************
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
	declareProperty("Vr0cut", m_vr0cut = 1.0);
	declareProperty("Vz0cut", m_vz0cut = 10.0);
	declareProperty("EnergyThreshold", m_energyThreshold = 0.025);
	declareProperty("GammaPhiCut", m_gammaPhiCut = 20.0);
	declareProperty("GammaThetaCut", m_gammaThetaCut = 20.0);
	declareProperty("GammaAngleCut", m_gammaAngleCut = 10.0);
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
	//initialize-start-line
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in initialize()" << endmsg;
	StatusCode status;
	//initialize-data-in-ori4c
	if (m_test4C == 1)
	{
		NTuplePtr nt3(ntupleSvc(), "FILE1/ori4c");
		if (nt3)
		{
			m_tuple3 = nt3;
		}
		else
		{
			m_tuple3 = ntupleSvc()->book("FILE1/ori4c", CLID_ColumnWiseTuple, "ks N-Tuple example");
			if (m_tuple3)
			{
				status = m_tuple3->addItem("chisq", m_chisq_ori);
				status = m_tuple3->addItem("momega", m_omega_ori);
				status = m_tuple3->addItem("mpi01", m_pi01_ori);
				status = m_tuple3->addItem("mpi02", m_pi02_ori);
				status = m_tuple3->addItem("mpi03", m_pi03_ori);
			}
			else
			{
				log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple3) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
	//initialize-data-in-fit4c
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
				status = m_tuple4->addItem("chisq", m_chisq_fit);
				status = m_tuple4->addItem("momega", m_omega_fit);
				status = m_tuple4->addItem("mpi01", m_pi01_fit);
				status = m_tuple4->addItem("mpi02", m_pi02_fit);
				status = m_tuple4->addItem("mpi03", m_pi03_fit);
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
				status = m_tuple5->addItem("chi2", m_chi2);
				status = m_tuple5->addItem("mpi01", m_mpi01);
				status = m_tuple5->addItem("mpi02", m_mpi02);
				status = m_tuple5->addItem("mpi03", m_mpi03);
				status = m_tuple5->addItem("momega", m_momega);
			}
			else
			{
				log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple5) << endmsg;
				return StatusCode::FAILURE;
			}
		}
		NTuplePtr nt6(ntupleSvc(), "FILE1/geff");
		if (nt6)
		{
			m_tuple6 = nt6;
		}
		else
		{
			m_tuple6 = ntupleSvc()->book("FILE1/geff", CLID_ColumnWiseTuple, "ks N-Tuple example");
			if (m_tuple6)
			{
				status = m_tuple6->addItem("fcos", m_fcos);
				status = m_tuple6->addItem("elow", m_elow);
			}
			else
			{
				log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple6) << endmsg;
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
StatusCode Omega::execute()
{
	std::cout << "execute()" << std::endl;
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;
	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
	int runNo = eventHeader->runNumber();
	int event = eventHeader->eventNumber();
	log << MSG::DEBUG << "run, evtnum = "
		<< runNo << " , "
		<< event << endreq;
	cout << "event " << event << endl;
	Ncut0++;
	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	log << MSG::DEBUG << "ncharg, nneu, tottks = "
		<< evtRecEvent->totalCharged() << " , "
		<< evtRecEvent->totalNeutral() << " , "
		<< evtRecEvent->totalTracks() << endreq;

	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
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
		if (fabs(Rvz0) >= m_vz0cut)													   //
			continue;																   //
		if (fabs(Rvxy0) >= m_vr0cut)												   //
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
	Vint iGam;																						   //
	iGam.clear();																					   //
	for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks(); i++)					   //
	{																								   //
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;										   //
		if (!(*itTrk)->isEmcShowerValid())															   //
			continue;																				   //
		RecEmcShower *emcTrk = (*itTrk)->emcShower();												   //
		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());									   //
		double dthe = 200.;																			   //
		double dphi = 200.;																			   //
		double dang = 200.;																			   //
		for (int j = 0; j < evtRecEvent->totalCharged(); j++)										   //
		{																							   //
			EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;									   //
			if (!(*jtTrk)->isExtTrackValid())														   //
				continue;																			   //
			RecExtTrack *extTrk = (*jtTrk)->extTrack();												   //
			if (extTrk->emcVolumeNumber() == -1)													   //
				continue;																			   //
			Hep3Vector extpos = extTrk->emcPosition();												   //
			double angd = extpos.angle(emcpos);														   //
			double thed = extpos.theta() - emcpos.theta();											   //
			double phid = extpos.deltaPhi(emcpos);													   //
			thed = fmod(thed + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) - CLHEP::pi;			   //
			phid = fmod(phid + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) - CLHEP::pi;			   //
			if (angd < dang)																		   //
			{																						   //
				dang = angd;																		   //
				dthe = thed;																		   //
				dphi = phid;																		   //
			}																						   //
		}																							   //
		if (dang >= 200)																			   //
			continue;																				   //
		double eraw = emcTrk->energy();																   //
		dthe = dthe * 180 / (CLHEP::pi);															   //
		dphi = dphi * 180 / (CLHEP::pi);															   //
		dang = dang * 180 / (CLHEP::pi);															   //
		if (eraw < m_energyThreshold)																   //筛选参数使用：m_energyThreshold
			continue;																				   //
		if (fabs(dang) < m_gammaAngleCut)															   //筛选参数使用：m_gammaAngleCut
			continue;																				   //
		iGam.push_back(i);																			   //
	}																								   //
	int nGam = iGam.size();																			   //
	log << MSG::DEBUG << "num Good Photon " << nGam << " , " << evtRecEvent->totalNeutral() << endreq; //
	if (nGam < 6)																					   //
	{																								   //
		return StatusCode::SUCCESS;																	   //
	}																								   //
	Ncut2++;																						   //
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
	RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin() + ipip[0]))->mdcKalTrack(); //Default is pion, for other particles:
	RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin() + ipim[0]))->mdcKalTrack(); //wvppTrk = WTrackParameter(mp, pipTrk->getZHelixP(), pipTrk->getZErrorP()); proton
	WTrackParameter wvpipTrk, wvpimTrk;											  //wvmupTrk = WTrackParameter(mmu, pipTrk->getZHelixMu(), pipTrk->getZErrorMu()); muon
	wvpipTrk = WTrackParameter(mpi, pipTrk->getZHelix(), pipTrk->getZError());	//wvepTrk = WTrackParameter(me, pipTrk->getZHelixE(), pipTrk->getZErrorE()); electron
	wvpimTrk = WTrackParameter(mpi, pimTrk->getZHelix(), pimTrk->getZError());	//wvkpTrk = WTrackParameter(mk, pipTrk->getZHelixK(), pipTrk->getZErrorK()); kaon
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
	vtxfit->AddTrack(0, wvpipTrk);												  //设定track0
	vtxfit->AddTrack(1, wvpimTrk);												  //设定track1
	vtxfit->AddVertex(0, vxpar, 0, 1);											  //设定顶点0
	if (!vtxfit->Fit(0))														  //
		return SUCCESS;															  //
	vtxfit->Swim(0);															  //更新信息
	Ncut6++;																	  //
	cout << "pass VFit" << endl;												  //
	//*********************************************************************************
	// Selection 7: 0C Selection ?
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
	//*********************************************************************************
	// Selection 7: 4C Selection
	//*********************************************************************************
	WTrackParameter wpip = vtxfit->wtrk(0);								 //
	WTrackParameter wpim = vtxfit->wtrk(1);								 //
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();			 //
	if (m_test4C == 1)													 //
	{																	 //
		HepLorentzVector ecms(0.034 * m_energy / 3.097, 0, 0, m_energy); //
		double chisq_4c = 9999;											 //
		HepLorentzVector ptrack0_fit;									 //记录ptrack_fit0~7
		HepLorentzVector ptrack1_fit;									 //
		HepLorentzVector ptrack2_fit;									 //
		HepLorentzVector ptrack3_fit;									 //
		HepLorentzVector ptrack4_fit;									 //
		HepLorentzVector ptrack5_fit;									 //
		HepLorentzVector ptrack6_fit;									 //
		HepLorentzVector ptrack7_fit;									 //
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
								kmfit->init();						  //
								kmfit->AddTrack(0, wpip);			  //
								kmfit->AddTrack(1, wpim);			  //
								kmfit->AddTrack(2, 0.0, g1Trk);		  //
								kmfit->AddTrack(3, 0.0, g2Trk);		  //
								kmfit->AddTrack(4, 0.0, g3Trk);		  //
								kmfit->AddTrack(5, 0.0, g4Trk);		  //
								kmfit->AddTrack(6, 0.0, g5Trk);		  //
								kmfit->AddTrack(7, 0.0, g6Trk);		  //
								kmfit->AddFourMomentum(0, ecms);	  //
								bool oksq = kmfit->Fit();			  //
								if (oksq)							  //
								{									  //
									double chi2 = kmfit->chisq();	 //
									if (chi2 < chisq_4c)			  //
									{								  //
										chisq_4c = chi2;			  //
										ptrack0_fit = kmfit->pfit(0); //
										ptrack1_fit = kmfit->pfit(1); //
										ptrack2_fit = kmfit->pfit(2); //
										ptrack3_fit = kmfit->pfit(3); //
										ptrack4_fit = kmfit->pfit(4); //
										ptrack5_fit = kmfit->pfit(5); //
										ptrack6_fit = kmfit->pfit(6); //
										ptrack7_fit = kmfit->pfit(7); //
										igamma0 = i1;				  //
										igamma1 = i2;				  //
										igamma2 = i3;				  //
										igamma3 = i4;				  //
										igamma4 = i5;				  //
										igamma5 = i6;				  //
									}
								}
							}
						}
					}
				}
			}
		}
		if (chisq_fit < 200)
		{
			double chisq_fit = 9999;						//
			double chisq_ori = 9999;						//
			double momega_fit;								//
			double momega_ori;								//
			double mpi01_fit;								//
			double mpi01_ori;								//
			double mpi02_fit;								//
			double mpi02_ori;								//
			double mpi03_fit;								//
			double mpi03_ori;								//
			HepLorentzVector ptrack0_ori = ppip[0];			//
			HepLorentzVector ptrack1_ori = ppim[0];			//
			HepLorentzVector ptrack2_ori = pGam[igamma0];   //
			HepLorentzVector ptrack3_ori = pGam[igamma1];   //
			HepLorentzVector ptrack4_ori = pGam[igamma2];   //
			HepLorentzVector ptrack5_ori = pGam[igamma3];   //
			HepLorentzVector ptrack6_ori = pGam[igamma4];   //
			HepLorentzVector ptrack7_ori = pGam[igamma5];   //
			HepLorentzVector ptrack_fit[7] = {ptrack0_fit,  //
											  ptrack2_fit,  //
											  ptrack3_fit,  //
											  ptrack4_fit,  //
											  ptrack5_fit,  //
											  ptrack6_fit,  //
											  ptrack7_fit}; //
			HepLorentzVector ptrack_ori[7] = {ptrack0_ori,  //
											  ptrack2_ori,  //
											  ptrack3_ori,  //
											  ptrack4_ori,  //
											  ptrack5_ori,  //
											  ptrack6_ori,  //
											  ptrack7_ori}; //
			double iomega_fit;								//
			double ipi01_fit;								//
			double ipi02_fit;								//
			double ipi03_fit;								//
			double chisqo_fit;								//
			double chisq1_fit;								//
			double chisq2_fit;								//
			double chisq3_fit;								//
			double iomega_ori;								//
			double ipi01_ori;								//
			double ipi02_ori;								//
			double ipi03_ori;								//
			double chisqo_ori;								//
			double chisq1_ori;								//
			double chisq2_ori;								//
			double chisq3_ori;								//
			for (int i = 0; i < 15; i++)
			{
				iomega_fit = (ptrack0_fit + ptrack1_fit + ptrack_fit[combine[i][0]] + ptrack_fit[combine[i][1]]).m();
				ipi01_fit = (ptrack_fit[combine[i][0]] + ptrack_fit[combine[i][1]]).m();
				ipi02_fit = (ptrack_fit[combine[i][2]] + ptrack_fit[combine[i][3]]).m();
				ipi03_fit = (ptrack_fit[combine[i][4]] + ptrack_fit[combine[i][5]]).m();
				iomega_ori = (ptrack0_ori + ptrack1_ori + ptrack_ori[combine[i][0]] + ptrack_ori[combine[i][1]]).m();
				ipi01_ori = (ptrack_ori[combine[i][0]] + ptrack_ori[combine[i][1]]).m();
				ipi02_ori = (ptrack_ori[combine[i][2]] + ptrack_ori[combine[i][3]]).m();
				ipi03_ori = (ptrack_ori[combine[i][4]] + ptrack_ori[combine[i][5]]).m();
				chisqo_fit = pow((iomega_fit - 0.782), 2);
				chisq1_fit = pow((ipi01_fit - 0.135), 2);
				chisq2_fit = pow((ipi02_fit - 0.135), 2);
				chisq3_fit = pow((ipi03_fit - 0.135), 2);
				chisqo_ori = pow((iomega_ori - 0.782), 2);
				chisq1_ori = pow((ipi01_ori - 0.135), 2);
				chisq2_ori = pow((ipi02_ori - 0.135), 2);
				chisq3_ori = pow((ipi03_ori - 0.135), 2);
				double chi2_fit = chisqo_fit;
				double chi2_ori = chisqo_ori;
				if (chi2_fit < chisq_fit)
				{
					chisq_fit = chi2_fit;
					momega_fit = iomega_fit;
					mpi01_fit = ipi01_fit;
					mpi02_fit = ipi01_fit;
					mpi03_fit = ipi01_fit;
				}
				if (chi2_ori < chisq_ori)
				{
					chisq_ori = chi2_ori;
					momega_ori = iomega_ori;
					mpi01_ori = ipi01_ori;
					mpi02_ori = ipi01_ori;
					mpi03_ori = ipi01_ori;
				}

			}
			if (1 == 1)
			{
				m_chisq_fit = chisq_fit;
				m_omega_fit = momega_fit;
				m_pi01_fit = mpi01_fit;
				m_pi02_fit = mpi02_fit;
				m_pi03_fit = mpi03_fit;
				m_tuple4->write();
				Ncut4++;
			}
			if (1 == 1)
			{
				m_chisq_ori = chisq_ori;
				m_omega_ori = momega_ori;
				m_pi01_ori = mpi01_ori;
				m_pi02_ori = mpi02_ori;
				m_pi03_ori = mpi03_ori;
				m_tuple3->write();
			}
		}
	}
	//*********************************************************************************
	// Selection 8: 5C Selection
	//     find the best combination over all possible pi+ pi- gamma gamma pair
	//*********************************************************************************
	cout << "before 5c" << endl; //准备5C声明
	if (m_test5C == 1)
	{
		HepLorentzVector ecms(0.034 * m_energy / 3.097, 0, 0, m_energy);
		double chisq = 9999.;
		int ig1 = -1;
		int ig2 = -1;
		int ig3 = -1;
		int ig4 = -1;
		int ig5 = -1;
		int ig6 = -1;
		for (int i1 = 0; i1 < nGam - 5; i1++) //循环所有6Gamma对
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
								for (int j = 0; j < 15; j++)
								{
									kmfit->init();
									kmfit->AddTrack(0, wpip);
									kmfit->AddTrack(1, wpim);
									kmfit->AddTrack(2, 0.0, g1Trk);
									kmfit->AddTrack(3, 0.0, g2Trk);
									kmfit->AddTrack(4, 0.0, g3Trk);
									kmfit->AddTrack(5, 0.0, g4Trk);
									kmfit->AddTrack(6, 0.0, g5Trk);
									kmfit->AddTrack(7, 0.0, g6Trk);
									kmfit->AddResonance(0, 0.135, combine[j][0] + 1, combine[j][1] + 1);
									kmfit->AddResonance(1, 0.135, combine[j][2] + 1, combine[j][3] + 1);
									kmfit->AddResonance(2, 0.135, combine[j][4] + 1, combine[j][5] + 1);
									kmfit->AddResonance(3, 0.782, 0, 1, combine[j][0] + 1, combine[j][1] + 1);
									kmfit->AddFourMomentum(4, ecms);
									if (!kmfit->Fit(0))
										continue;
									if (!kmfit->Fit(1))
										continue;
									if (!kmfit->Fit(2))
										continue;
									if (!kmfit->Fit(3))
										continue;
									if (!kmfit->Fit(4))
										continue;
									bool oksq = kmfit->Fit();
									if (oksq)
									{
										double chi2 = kmfit->chisq();
										if (chi2 < chisq)
										{
											int outcheck[6] = {i1, i2, i3, i4, i5, i6};
											chisq = chi2;
											ig1 = iGam[outcheck[combine[j][0] - 1]];
											ig2 = iGam[outcheck[combine[j][1] - 1]];
											ig3 = iGam[outcheck[combine[j][2] - 1]];
											ig4 = iGam[outcheck[combine[j][3] - 1]];
											ig5 = iGam[outcheck[combine[j][4] - 1]];
											ig6 = iGam[outcheck[combine[j][5] - 1]];
											cout << "通过5c拟合" << endl;
											cout << "通过5c拟合的chisq = " << chisq << endl;
										}
									}
								}
							}
						}
					}
				}
			}
		}
		log << MSG::INFO << " chisq = " << chisq << endreq;

		if (chisq < 200)
		{
			cout << "准备5c写入" << endl;
			RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + ig1))->emcShower();
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + ig2))->emcShower();
			RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + ig3))->emcShower();
			RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + ig4))->emcShower();
			RecEmcShower *g5Trk = (*(evtRecTrkCol->begin() + ig5))->emcShower();
			RecEmcShower *g6Trk = (*(evtRecTrkCol->begin() + ig6))->emcShower();
			kmfit->init();
			kmfit->AddTrack(0, wpip);
			kmfit->AddTrack(1, wpim);
			kmfit->AddTrack(2, 0.0, g1Trk);
			kmfit->AddTrack(3, 0.0, g2Trk);
			kmfit->AddTrack(4, 0.0, g3Trk);
			kmfit->AddTrack(5, 0.0, g4Trk);
			kmfit->AddTrack(6, 0.0, g5Trk);
			kmfit->AddTrack(7, 0.0, g6Trk);
			kmfit->AddResonance(0, 0.135, 2, 3);
			kmfit->AddResonance(1, 0.135, 4, 5);
			kmfit->AddResonance(2, 0.135, 6, 7);
			kmfit->AddResonance(3, 0.782, 0, 1, 2, 3);
			kmfit->AddFourMomentum(4, ecms);
			bool oksq = kmfit->Fit();
			if (oksq)
			{
				cout << "重新拟合成功" << endl;
				HepLorentzVector ppi01 = kmfit->pfit(2) + kmfit->pfit(3);
				HepLorentzVector ppi02 = kmfit->pfit(4) + kmfit->pfit(5);
				HepLorentzVector ppi03 = kmfit->pfit(6) + kmfit->pfit(7);
				HepLorentzVector pomega = kmfit->pfit(0) + kmfit->pfit(1) + kmfit->pfit(2) + kmfit->pfit(3);
				//HepLorentzVector prho0 = kmfit->pfit(0) + kmfit->pfit(1);
				//HepLorentzVector prhop = ppi0 + kmfit->pfit(0);
				//HepLorentzVector prhom = ppi0 + kmfit->pfit(1);
				//double eg1 = (kmfit->pfit(2)).e();
				//double eg2 = (kmfit->pfit(3)).e();
				//double fcos = abs(eg1 - eg2) / ppi0.rho();
				if (1 == 1)
				{
					m_chi2 = kmfit->chisq();
					//m_mrh0 = prho0.m();
					//m_mrhp = prhop.m();
					//m_mrhm = prhom.m();
					m_mpi01 = ppi01.m();
					m_mpi02 = ppi02.m();
					m_mpi03 = ppi03.m();
					m_momega = pomega.m();
					m_tuple5->write();
					Ncut5++;
					cout << "c5写入成功" << endl;
				}
				//
				//  Measure the photon detection efficiences via
				//          J/psi -> rho0 pi0
				//
				//if (fabs(prho0.m() - 0.770) < 0.150)
				if (1 == 0)
				{
					//if (fabs(fcos) < 0.99)
					if (1 == 0)
					{
						//m_fcos = (eg1 - eg2) / ppi0.rho();
						//m_elow = (eg1 < eg2) ? eg1 : eg2;
						m_fcos = 0;
						m_elow = 0;
						m_tuple6->write();
					}
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
