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
#include "OmegaAlg/Omega.h"
//#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"
#include <vector>

//*********************************************************************************************************
//***                                                变量设置                                            ***
//*********************************************************************************************************
//const double twopi = 6.2831853;
//const double pi = 3.1415927;
const double mpi = 0.13957;
const double xmass[5] = {0.000511, 0.105658, 0.139570, 0.493677, 0.938272};
//const double velc = 29.9792458;  tof_path unit in cm.
const double velc = 299.792458; // tof path unit in mm
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
				status = m_tuple4->addItem("chi2", m_chi1);
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
				//status = m_tuple5->addItem("mrh0", m_mrh0);
				//status = m_tuple5->addItem("mrhp", m_mrhp);
				//status = m_tuple5->addItem("mrhm", m_mrhm);
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
	//  log << MSG::INFO << "get event tag OK" << endreq;
	log << MSG::DEBUG << "ncharg, nneu, tottks = "
		<< evtRecEvent->totalCharged() << " , "
		<< evtRecEvent->totalNeutral() << " , "
		<< evtRecEvent->totalTracks() << endreq;

	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
	//*********************************************************************************
	// Selection 1: Good Charged Track Selection
	//*********************************************************************************
	Vint iGood, ipip, ipim;
	iGood.clear();        			     			     			     	        //存good的charge track的编号
	ipip.clear();        			     			     			     	        //存good的Pi+的编号
	ipim.clear();         			     			     			     	        //存good的Pi-的编号
	Vp4 ppip, ppim;
	ppip.clear();					     			     			     			//存Pi+的四动量
	ppim.clear();						     			     			     		//存Pi-的四动量
	int nCharge;
	nCharge = 0;			     			     			     			        //存带电总量
	Hep3Vector xorigin(0, 0, 0);
	IVertexDbSvc *vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if (vtxsvc->isVertexValid())
	{
		double *dbv = vtxsvc->PrimaryVertex();
		double *vv = vtxsvc->SigmaPrimaryVertex();
		//HepVector dbv = m_reader.PrimaryVertex(runNo);
		//HepVector vv = m_reader.SigmaPrimaryVertex(runNo);
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);
	}                                                                               //存平均对撞顶点
	for (int i = 0; i < evtRecEvent->totalCharged(); i++)                           //循环所有带电track
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
		if (!(*itTrk)->isMdcTrackValid())
			continue;
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
		double pch = mdcTrk->p();                                                   //track开始的动量
		double x0 = mdcTrk->x();                                                    //track开始的x坐标
		double y0 = mdcTrk->y();                                                    //track开始的y坐标
		double z0 = mdcTrk->z();                                                    //track开始的z坐标
		double phi0 = mdcTrk->helix(1);
		//表示螺旋线参数
		//0: d0 -> 螺旋线到对撞顶点的最小距离
		//1: phi0 -> 最小距离的xy平面相角
		//2: kappa
		//3: d
		//4: tan(lamda)
		double xv = xorigin.x();                                                    //对撞顶点的x坐标
		double yv = xorigin.y();                                                    //对撞顶点的y坐标
		double Rxy = (x0 - xv) * cos(phi0) + (y0 - yv) * sin(phi0);                 //计算顶点到track开始位置的xy平面距离
		HepVector a = mdcTrk->helix();					                            //helix存入a
		HepSymMatrix Ea = mdcTrk->err();				                            //helix-err存入Ea
		HepPoint3D point0(0., 0., 0.);					                            //探测器原点
		HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);                          //对撞原点
		VFHelix helixip(point0, a, Ea);					                            //组合数据{(0,0,0),helix(),err()}
		helixip.pivot(IP);								                            //调整数据{(0,0,0),helix()-IP,err()-IP}
		HepVector vecipa = helixip.a();					                            //提取新helix
		double Rvxy0 = fabs(vecipa[0]);					   
		double Rvz0 = vecipa[3];						   
		double Rvphi0 = vecipa[1];						   
		//if (1 == 1)										                        
		//{
		//	m_vx0 = x0;		                                                        //output-track开始位置x
		//	m_vy0 = y0;		                                                        //output-track开始位置y
		//	m_vz0 = z0;		                                                        //output-track开始位置z
		//	m_vr0 = Rxy;	                                                        //output-track开始位置-对撞原点位置的xy平面距离
		//	m_rvxy0 = Rvxy0;                                                        //output-相对位置的d0 -> 螺旋线到对撞顶点的最小距离
		//	m_rvz0 = Rvz0;	                                                        //output-相对位置的dz
		//	m_rvphi0 = Rvphi0;                                                      //相对位置的phi0
		//	m_tuple1->write();
		//}
		if(fabs(Rvz0) >= m_vz0cut) continue;
		if(fabs(Rvxy0) >= m_vr0cut) continue;
		//if (fabs(Rvz0) >= 10.0) continue;
		//if (fabs(Rvxy0) >= 1) continue;
		iGood.push_back(i);
		nCharge += mdcTrk->charge();
	}
	int nGood = iGood.size();                                                       //ngood记录good trach数量，ncharge记录带电量
	log << MSG::DEBUG << "ngood, totcharge = " << nGood << " , " << nCharge << endreq;
	if ((nGood != 2) || (nCharge != 0))                                             //参数设定：带电总数，带电总量
	{
		return StatusCode::SUCCESS;
	}
	Ncut1++;
	//*********************************************************************************
	// Selection 2: Good Photon Selection
	//*********************************************************************************
	Vint iGam;
	iGam.clear();
	for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks(); i++)  //循环所有中性track
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
		if (!(*itTrk)->isEmcShowerValid())
			continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
		// find the nearest charged track
		double dthe = 200.;
		double dphi = 200.;
		double dang = 200.;
		for (int j = 0; j < evtRecEvent->totalCharged(); j++)                       //循环所有带电track
		{
			EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
			if (!(*jtTrk)->isExtTrackValid())
				continue;
			RecExtTrack *extTrk = (*jtTrk)->extTrack();
			if (extTrk->emcVolumeNumber() == -1)
				continue;
			Hep3Vector extpos = extTrk->emcPosition();
			//double ctht = extpos.cosTheta(emcpos);
			double angd = extpos.angle(emcpos);
			double thed = extpos.theta() - emcpos.theta();
			double phid = extpos.deltaPhi(emcpos);
			thed = fmod(thed + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) - CLHEP::pi;
			phid = fmod(phid + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) - CLHEP::pi;
			if (angd < dang)
			{
				dang = angd;
				dthe = thed;
				dphi = phid;
			}
		}
		if (dang >= 200)
			continue;                                                               //排除中性流是有带电流辐射产生,动量方向夹角
		//********************
		//eraw = emcTrk->energy()
		//代表能量
		//dtheta,dphi,dangular
		//代表三维角度
		//********************
		double eraw = emcTrk->energy();
		dthe = dthe * 180 / (CLHEP::pi);
		dphi = dphi * 180 / (CLHEP::pi);
		dang = dang * 180 / (CLHEP::pi);                                            //单位转化
		//if (1 == 1)						                                        
		//{
		//	m_dthe = dthe;                                                          //dthe
		//	m_dphi = dphi;                                                          //dphi
		//	m_dang = dang;                                                          //dang
		//	m_eraw = eraw;                                                          //eraw
		//	m_tuple2->write();
		//}
		if (eraw < m_energyThreshold)                                               //参数设定：selection能量界限 > m_energyThreshold = 0.025
			continue; 
		//if((fabs(dthe) < m_gammaThetaCut) && (fabs(dphi)<m_gammaPhiCut) ) continue;
		if (fabs(dang) < m_gammaAngleCut)
			continue;                                                               //参数设定：selection角度 > m_gammaAngleCut = 10
		iGam.push_back(i);
	}
	int nGam = iGam.size(); //nGam记录good photon数量
	log << MSG::DEBUG << "num Good Photon " << nGam << " , " << evtRecEvent->totalNeutral() << endreq;
	if (nGam < 6)                                                                   //参数设定：nGam >= 6
	{
		return StatusCode::SUCCESS;
	}
	Ncut2++;
	//*********************************************************************************
	// Selection 3: dE/dx Selection
	//*********************************************************************************
	if (m_checkDedx == 1)
	{
		for (int i = 0; i < nGood; i++)
		{
			EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
			if (!(*itTrk)->isMdcTrackValid())
				continue;
			if (!(*itTrk)->isMdcDedxValid())
				continue;
			RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
			RecMdcDedx *dedxTrk = (*itTrk)->mdcDedx();
		}
	}
	//*********************************************************************************
	// Selection 4: Tof Selection
	//*********************************************************************************
	if (m_checkTof == 1)
	{
		for (int i = 0; i < nGood; i++)
		{
			EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
			if (!(*itTrk)->isMdcTrackValid())
				continue;
			if (!(*itTrk)->isTofTrackValid())
				continue;

			RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
			SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();

			double ptrk = mdcTrk->p();

			SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
			for (; iter_tof != tofTrkCol.end(); iter_tof++)
			{
				TofHitStatus *status = new TofHitStatus;
				status->setStatus((*iter_tof)->status());
				if (!(status->is_barrel()))
				{ //endcap
					if (!(status->is_counter()))
						continue; // ?
					if (status->layer() != 0)
						continue;					   //layer1
					double path = (*iter_tof)->path(); // ?
					double tof = (*iter_tof)->tof();
					double ph = (*iter_tof)->ph();
					double rhit = (*iter_tof)->zrhit();
					double qual = 0.0 + (*iter_tof)->quality();
					double cntr = 0.0 + (*iter_tof)->tofID();
					double texp[5];
					for (int j = 0; j < 5; j++)
					{
						double gb = ptrk / xmass[j];
						double beta = gb / sqrt(1 + gb * gb);
						texp[j] = 10 * path / beta / velc;
					}
				}
				else
				{
					if (!(status->is_counter()))
						continue;
					if (status->layer() == 1)
					{
						double path = (*iter_tof)->path();
						double tof = (*iter_tof)->tof();
						double ph = (*iter_tof)->ph();
						double rhit = (*iter_tof)->zrhit();
						double qual = 0.0 + (*iter_tof)->quality();
						double cntr = 0.0 + (*iter_tof)->tofID();
						double texp[5];
						for (int j = 0; j < 5; j++)
						{
							double gb = ptrk / xmass[j];
							double beta = gb / sqrt(1 + gb * gb);
							texp[j] = 10 * path / beta / velc;
						}
					}

					if (status->layer() == 2)
					{
						double path = (*iter_tof)->path();
						double tof = (*iter_tof)->tof();
						double ph = (*iter_tof)->ph();
						double rhit = (*iter_tof)->zrhit();
						double qual = 0.0 + (*iter_tof)->quality();
						double cntr = 0.0 + (*iter_tof)->tofID();
						double texp[5];
						for (int j = 0; j < 5; j++)
						{
							double gb = ptrk / xmass[j];
							double beta = gb / sqrt(1 + gb * gb);
							texp[j] = 10 * path / beta / velc;
						}
					}
				}
				delete status;
			}
		}
	}
	//*********************************************************************************
	// Calculation 1: 4-momentum to each photon
	//*********************************************************************************
	Vp4 pGam;                                                                       //存光子4动量
	pGam.clear();
	for (int i = 0; i < nGam; i++)                                                  //循环所有nGam
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGam[i];
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		double eraw = emcTrk->energy();
		double phi = emcTrk->phi();
		double the = emcTrk->theta();
		HepLorentzVector ptrk;
		ptrk.setPx(eraw * sin(the) * cos(phi));
		ptrk.setPy(eraw * sin(the) * sin(phi));
		ptrk.setPz(eraw * cos(the));
		ptrk.setE(eraw);
		//ptrk = ptrk.boost(-0.011,0,0);//boost to cms
		pGam.push_back(ptrk);
	}
	cout << "before pid" << endl; //准备pid声明
	//*********************************************************************************
	// Calculation 2: 4-momentum to each charged track
	//*********************************************************************************
	ParticleID *pid = ParticleID::instance();
	for (int i = 0; i < nGood; i++) //循环所有nGood
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
		//if(pid) delete pid;
		pid->init();
		pid->setMethod(pid->methodProbability());
		//pid->setMethod(pid->methodLikelihood());//for Likelihood Method
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE()); //use PID sub-system
		pid->identify(pid->onlyPion() | pid->onlyKaon());								   //seperater Pion/Kaon
		//  pid->identify(pid->onlyPion());
		//  pid->identify(pid->onlyKaon());
		pid->calculate();
		if (!(pid->IsPidInfoValid()))
			continue;
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
		//if(pid->probPion() < 0.001 || (pid->probPion() < pid->probKaon())) continue;
		if (pid->probPion() < 0.001)
			continue; //PID筛选条件
		//if(pid->pdf(2)<pid->pdf(3)) continue; //for Likelihood Method(0=electron 1=muon 2=pion 3=kaon 4=proton)
		RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack(); //对于ParticleID, 用RecMdcKalTrack代替RecMdcTrack
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);	 //PID can set to electron, muon, pion, kaon and proton;The default setting is pion
		if (mdcKalTrk->charge() > 0)						 //判断正负
		{
			ipip.push_back(iGood[i]);
			HepLorentzVector ptrk;
			ptrk.setPx(mdcKalTrk->px());
			ptrk.setPy(mdcKalTrk->py());
			ptrk.setPz(mdcKalTrk->pz());
			double p3 = ptrk.mag();
			ptrk.setE(sqrt(p3 * p3 + mpi * mpi));

			//ptrk = ptrk.boost(-0.011,0,0);//boost to cms

			ppip.push_back(ptrk);
		}
		else
		{
			ipim.push_back(iGood[i]);
			HepLorentzVector ptrk;
			ptrk.setPx(mdcKalTrk->px());
			ptrk.setPy(mdcKalTrk->py());
			ptrk.setPz(mdcKalTrk->pz());
			double p3 = ptrk.mag();
			ptrk.setE(sqrt(p3 * p3 + mpi * mpi));

			//ptrk = ptrk.boost(-0.011,0,0);//boost to cms

			ppim.push_back(ptrk);
		}
	}
	/*
	   for(int i = 0; i < nGood; i++) {//for Omega without PID
	   EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
	   RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack(); 
	   if(mdcTrk->charge() >0) {
	   ipip.push_back(iGood[i]);
	   HepLorentzVector ptrk;
	   ptrk.setPx(mdcTrk->px());
	   ptrk.setPy(mdcTrk->py());
	   ptrk.setPz(mdcTrk->pz());
	   double p3 = ptrk.mag();
	   ptrk.setE(sqrt(p3*p3+mpi*mpi));
	   ppip.push_back(ptrk);
	   } else {
	   ipim.push_back(iGood[i]);
	   HepLorentzVector ptrk;
	   ptrk.setPx(mdcTrk->px());
	   ptrk.setPy(mdcTrk->py());
	   ptrk.setPz(mdcTrk->pz());
	   double p3 = ptrk.mag();
	   ptrk.setE(sqrt(p3*p3+mpi*mpi));
	   ppim.push_back(ptrk);
	   }
	   }// without PID
	   */

	int npip = ipip.size();
	int npim = ipim.size();
	if (npip * npim != 1)
		return SUCCESS; //要求只有一个P+,P-
	Ncut3++;			//计数
	//*********************************************************************************
	// Selection 5: Gamma Pair Selection, check ppi0, pTot
	//*********************************************************************************
	HepLorentzVector pTot;
	for (int i = 0; i < nGam - 1; i++)
	{
		for (int j = i + 1; j < nGam; j++)
		{
			HepLorentzVector p2g = pGam[i] + pGam[j];
			pTot = ppip[0] + ppim[0];
			pTot += p2g;
		}
	}
	//*********************************************************************************
	// Selection 6: Vertex fit Selection, check ppi0, pTot
	//*********************************************************************************
	RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin() + ipip[0]))->mdcKalTrack();
	RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin() + ipim[0]))->mdcKalTrack();
	WTrackParameter wvpipTrk, wvpimTrk;
	wvpipTrk = WTrackParameter(mpi, pipTrk->getZHelix(), pipTrk->getZError());
	wvpimTrk = WTrackParameter(mpi, pimTrk->getZHelix(), pimTrk->getZError());
	/* 	Default is pion, for other particles:
	   	wvppTrk = WTrackParameter(mp, pipTrk->getZHelixP(), pipTrk->getZErrorP());//proton
	   	wvmupTrk = WTrackParameter(mmu, pipTrk->getZHelixMu(), pipTrk->getZErrorMu());//muon
	   	wvepTrk = WTrackParameter(me, pipTrk->getZHelixE(), pipTrk->getZErrorE());//electron
	   	wvkpTrk = WTrackParameter(mk, pipTrk->getZHelixK(), pipTrk->getZErrorK());//kaon     */
	HepPoint3D vx(0., 0., 0.);
	HepSymMatrix Evx(3, 0);
	double bx = 1E+6;
	double by = 1E+6;
	double bz = 1E+6;
	Evx[0][0] = bx * bx;
	Evx[1][1] = by * by;
	Evx[2][2] = bz * bz;
	VertexParameter vxpar;
	vxpar.setVx(vx);
	vxpar.setEvx(Evx);
	VertexFit *vtxfit = VertexFit::instance();
	vtxfit->init();
	vtxfit->AddTrack(0, wvpipTrk);	 //设定track0
	vtxfit->AddTrack(1, wvpimTrk);	 //设定track1
	vtxfit->AddVertex(0, vxpar, 0, 1); //设定顶点0
	if (!vtxfit->Fit(0))
		return SUCCESS;
	vtxfit->Swim(0); //更新信息
	//*********************************************************************************
	// Selection 7: 4C Selection
	//*********************************************************************************
	WTrackParameter wpip = vtxfit->wtrk(0);
	WTrackParameter wpim = vtxfit->wtrk(1);
	//KinematicFit * kmfit = KinematicFit::instance();
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
	cout << "before 4c" << endl;
	if (m_test4C == 1)
	{
		HepLorentzVector ecms(0.034 * m_energy / 3.097, 0, 0, m_energy);
		double chisq = 9999.; //这个是ka^2
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
					                if (chi2 < chisq)
					                {
						                chisq = chi2;
						                ig1 = iGam[i1];
						                ig2 = iGam[i2];
						                ig3 = iGam[i3];
						                ig4 = iGam[i4];
						                ig5 = iGam[i5];
						                ig6 = iGam[i6];
					                }
				                }
							}
						}
					}
				}
			}
		}

		if (chisq < 200)
		{
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
			kmfit->AddFourMomentum(0, ecms);
			bool oksq = kmfit->Fit();
			if (oksq)
			{
				m_chi1 = kmfit->chisq();
				m_tuple4->write();
				Ncut4++;
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
		//double ecms = 3.097;
		HepLorentzVector ecms(0.034 * m_energy / 3.097, 0, 0, m_energy);
		double chisq = 9999.;
		int ig1 = -1;
		int ig2 = -1;
		int ig3 = -1;
		int ig4 = -1;
		int ig5 = -1;
		int ig6 = -1;
		int combine[15][6] = {{1,2,3,4,5,6},
		{1,2,3,5,4,6},
		{1,2,3,6,4,5},
		{1,3,2,4,5,6},
		{1,3,2,5,4,6},
		{1,3,2,6,4,5},
		{1,4,3,2,5,6},
		{1,4,3,5,2,6},
		{1,4,3,6,2,5},
		{1,5,3,4,2,6},
		{1,5,3,2,4,6},
		{1,5,3,6,4,2},
		{1,6,3,4,5,2},
		{1,6,3,5,4,2},
		{1,6,3,2,4,5}};
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
				                	kmfit->AddResonance(0, 0.135, combine[j][0]+1, combine[j][1]+1);
				                	kmfit->AddResonance(1, 0.135, combine[j][2]+1, combine[j][3]+1);
				                	kmfit->AddResonance(2, 0.135, combine[j][4]+1, combine[j][5]+1);
				                	kmfit->AddResonance(3, 0.782, 0, 1, combine[j][0]+1, combine[j][1]+1);
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
											int outcheck[6] = {i1,i2,i3,i4,i5,i6};
						            	    chisq = chi2;
						            	    ig1 = iGam[outcheck[combine[j][0]]-1];
						            	    ig2 = iGam[outcheck[combine[j][1]]-1];
						            	    ig3 = iGam[outcheck[combine[j][2]]-1];
						            	    ig4 = iGam[outcheck[combine[j][3]]-1];
						            	    ig5 = iGam[outcheck[combine[j][4]]-1];
						            	    ig6 = iGam[outcheck[combine[j][5]]-1];
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
					cout << "写入成功" << endl;
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
						Ncut6++;
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
	cout << "nGam>=2:              " << Ncut2 << endl;
	cout << "Pass Pid:             " << Ncut3 << endl;
	cout << "Pass 4C:              " << Ncut4 << endl;
	cout << "Pass 5C:              " << Ncut5 << endl;
	cout << "J/psi->rho0 pi0:      " << Ncut6 << endl;
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}
