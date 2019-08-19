#容器写入规范
if (1)
{
    Omega::Omega(const std::string &name, ISvcLocator *pSvcLocator) : Algorithm(name, pSvcLocator)
    {
        declareProperty("Vr0cut", m_vr0cut = 1.0);
        declareProperty("Vz0cut", m_vz0cut = 10.0);
    }
}

#变量写入规范
if (1)
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

#数据写入规范
if (1)
{
    m_chi2 = kmfit->chisq();
    m_mpi01 = ppi01.m();
    m_mpi02 = ppi02.m();
    m_mpi03 = ppi03.m();
    m_momega = pomega.m();
    m_tuple5->write();
    Ncut5++;
    cout << "c5写入成功" << endl;
}




if (m_test111C == 1)
	{
		HepLorentzVector ecms(0.034 * m_energy / 3.097, 0, 0, m_energy); //
		double chisq_1 = 9999;											 //
		double chisq_2 = 9999;											 //
		double chisq_3 = 9999;											 //
		int ig1 = -1;													 //
		int ig2 = -1;													 //
		int ig3 = -1;													 //
		int ig4 = -1;													 //
		int ig5 = -1;													 //
		int ig6 = -1;													 //
		int process = 0;
		for (int i1 = 0; i1 < nGam - 1; i1++)
		{
			RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();
			for (int i2 = i1 + 1; i2 < nGam; i2++)
			{
				RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();
				kmfit->init();
				kmfit->AddTrack(0, 0.0, g1Trk);
				kmfit->AddTrack(1, 0.0, g2Trk);
				kmfit->AddResonance(0, 0.135, 0, 1);
				bool oksq = kmfit->Fit();
				if (oksq)
				{
					double chi2 = kmfit->chisq();
					if (chi2 < chisq_1)
					{
						chisq_1 = chi2;
						ig1 = i1;
						ig2 = i2;
						process = 1;
					}
				}
			}
		}
		if (process == 1)
		{
			for (int i1 = 0; i1 < nGam - 1; i1++)
			{
				if (i1 == ig1 || i1 == ig2)
				{
					continue;
				}
				RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();
				for (int i2 = i1 + 1; i2 < nGam; i2++)
				{
					if (i2 == ig1 || i2 == ig2)
					{
						continue;
					}
					RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();
					kmfit->init();
					kmfit->AddTrack(0, 0.0, g1Trk);
					kmfit->AddTrack(1, 0.0, g2Trk);
					kmfit->AddResonance(0, 0.135, 0, 1);
					bool oksq = kmfit->Fit();
					if (oksq)
					{
						double chi2 = kmfit->chisq();
						if (chi2 < chisq_2)
						{
							chisq_2 = chi2;
							ig3 = i1;
							ig4 = i2;
							process = 2;
						}
					}
				}
			}
		}
		if (process == 2)
		{
			for (int i1 = 0; i1 < nGam - 1; i1++)
			{
				if (i1 == ig1 || i1 == ig2 || i1 == ig3 || i1 == ig4)
				{
					continue;
				}
				RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i1]))->emcShower();
				for (int i2 = i1 + 1; i2 < nGam; i2++)
				{
					if (i2 == ig1 || i2 == ig2 || i2 == ig3 || i2 == ig4)
					{
						continue;
					}
					RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[i2]))->emcShower();
					kmfit->init();
					kmfit->AddTrack(0, 0.0, g1Trk);
					kmfit->AddTrack(1, 0.0, g2Trk);
					kmfit->AddResonance(0, 0.135, 0, 1);
					bool oksq = kmfit->Fit();
					if (oksq)
					{
						double chi2 = kmfit->chisq();
						if (chi2 < chisq_3)
						{
							chisq_3 = chi2;
							ig5 = i1;
							ig6 = i2;
							process = 3;
						}
					}
				}
			}
		}
		if (process == 3)
		{
			HepLorentzVector ptrack0 = pGam[ig1];
			HepLorentzVector ptrack1 = pGam[ig2];
			HepLorentzVector ptrack2 = pGam[ig3];
			HepLorentzVector ptrack3 = pGam[ig4];
			HepLorentzVector ptrack4 = pGam[ig5];
			HepLorentzVector ptrack5 = pGam[ig6];
			double chisq_o = 9999;
			double momega_111c = -1;
			double mpi01_111c = -1;
			double mpi02_111c = -1;
			double mpi03_111c = -1;
			double chisq_o1 = pow((ptrack0 + ptrack1 + ppip[0] + ppim[0]).m() - 0.782, 2);
			double chisq_o2 = pow((ptrack2 + ptrack3 + ppip[0] + ppim[0]).m() - 0.782, 2);
			double chisq_o3 = pow((ptrack4 + ptrack5 + ppip[0] + ppim[0]).m() - 0.782, 2);
			if (chisq_o1 < chisq_o)
			{
				chisq_o = chisq_o1;
				momega_111c = (ptrack0 + ptrack1 + ppip[0] + ppim[0]).m();
				mpi01_111c = (ptrack0 + ptrack1).m();
				mpi02_111c = (ptrack2 + ptrack3).m();
				mpi03_111c = (ptrack4 + ptrack5).m();
				process = 4;
			}
			if (chisq_o2 < chisq_o)
			{
				chisq_o = chisq_o2;
				momega_111c = (ptrack2 + ptrack3 + ppip[0] + ppim[0]).m();
				mpi01_111c = (ptrack2 + ptrack3).m();
				mpi02_111c = (ptrack0 + ptrack1).m();
				mpi03_111c = (ptrack4 + ptrack5).m();
				process = 4;
			}
			if (chisq_o3 < chisq_o)
			{
				chisq_o = chisq_o3;
				momega_111c = (ptrack4 + ptrack5 + ppip[0] + ppim[0]).m();
				mpi01_111c = (ptrack4 + ptrack5).m();
				mpi02_111c = (ptrack2 + ptrack3).m();
				mpi03_111c = (ptrack0 + ptrack1).m();
				process = 4;
			}
			if (process == 4)
			{
				m_omega_111c = momega_111c;
				m_pi01_111c = mpi01_111c;
				m_pi02_111c = mpi01_111c;
				m_pi03_111c = mpi01_111c;
				m_chisqo = chisq_o;
				m_chisq1 = chisq_1;
				m_chisq2 = chisq_2;
				m_chisq3 = chisq_3;
				m_tuple7->write();
				cout << "write tuple7" << endl;
			}
		}
	}