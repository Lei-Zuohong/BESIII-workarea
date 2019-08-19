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