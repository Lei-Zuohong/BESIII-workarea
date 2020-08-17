//using CLHEP::HepLorentzVector;
#pragma region 类定义
class my_constant
{
public:
    double pi;
    double mpiz;
    double mpipm;
    double me;
    double mkaon;
    double mmiu;
    my_constant()
    {
        pi = 3.1415926;
        mpiz = 0.13498;
        mpipm = 0.13957;
        me = 0.000511;
        mkaon = 0.493677;
        mmiu = 0.105658;
    }
};
#pragma endregion
#pragma region 变量操作
double *my_newlist1d(const int a,
                     double inital = 0)
{ //新建一个指定长度的double形数组

    double *group = new double[a];
    for (int i = 0; i < a; i++)
    {
        group[i] = inital;
    }
    return group;
}
double **my_newlist2d(const int a,
                      const int b,
                      double inital = 0)
{ //新建一个指定长度的double形二维数组
    double **group = new double *[a];
    for (int i = 0; i < a; i++)
    {
        group[i] = new double[b];
        for (int j = 0; j < b; j++)
        {
            group[i][j] = inital;
        }
    }
    return group;
}
HepLorentzVector *my_newheplorentzvector(const int a)
{ //新建一个指定长度的HepLorentzVector形数组
    HepLorentzVector *group = new HepLorentzVector[a];
    return group;
}
void my_tracktovalue(HepLorentzVector track,
                     NTuple::Item<double> &value_m,
                     NTuple::Item<double> &value_p,
                     NTuple::Item<double> &value_a,
                     NTuple::Item<double> &value_pe,
                     NTuple::Item<double> &value_px,
                     NTuple::Item<double> &value_py,
                     NTuple::Item<double> &value_pz)
{
    value_m = track.m();
    value_p = track.rho();
    value_a = track.cosTheta();
    value_pe = track.e();
    value_px = track.px();
    value_py = track.py();
    value_pz = track.pz();
}

#pragma endregion
#pragma region 函数计算
double my_helicityangle(HepLorentzVector track_gamma,
                        HepLorentzVector track_pi)
{ //计算两个track的helicity angle
    HepLorentzVector track_bgamma = track_gamma.boost(-track_pi.boostVector());
    double px1 = track_bgamma.px();
    double py1 = track_bgamma.py();
    double pz1 = track_bgamma.pz();
    double e1 = track_bgamma.e();
    double px2 = track_pi.px();
    double py2 = track_pi.py();
    double pz2 = track_pi.pz();
    double e2 = track_pi.e();
    double angle = (px1 * px2 + py1 * py2 + pz1 * pz2) / (std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1)) / (std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2));
    return angle;
}
double my_angle(HepLorentzVector track1,
                HepLorentzVector track2)
{
    double px1 = track1.px();
    double py1 = track1.py();
    double pz1 = track1.pz();
    double e1 = track1.e();
    double px2 = track2.px();
    double py2 = track2.py();
    double pz2 = track2.pz();
    double e2 = track2.e();
    double angle = (px1 * px2 + py1 * py2 + pz1 * pz2) / (std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1)) / (std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2));
    return angle;
}
#pragma endregion
#pragma region 事例号统计
void my_seriesinit(vector<int> &seriesrun,
                   vector<int> &seriesnum,
                   vector<int> &seriesnum1,
                   vector<int> &seriesnum2,
                   vector<int> &seriesnum3,
                   vector<int> &seriesnum4,
                   vector<int> &seriesnum5,
                   int runNo,
                   int &firstrun)
{
    if (firstrun == 0)
    {
        // first run, need clear data
        // 第一次运行，清除数组
        seriesrun.clear();
        seriesnum.clear();
        seriesnum1.clear();
        seriesnum2.clear();
        seriesnum3.clear();
        seriesnum4.clear();
        seriesnum5.clear();
        firstrun = 1;
    }
    int check_exist = 0;
    int check_i = 0;
    for (int i = 0; i < seriesrun.size(); i++)
    {
        // search in seriesrun
        // 在seriesrun中查找
        if (runNo == seriesrun[i])
        {
            check_exist = 1;
            check_i = i;
        }
    }
    if (check_exist == 1)
    {
        // exist, count
        // 存在，进行一次计数
        seriesnum[check_i] += 1;
    }
    else
    {
        // not exist, create
        // 不存在，加入容器
        seriesrun.push_back(runNo);
        seriesnum.push_back(1);
        seriesnum1.push_back(0);
        seriesnum2.push_back(0);
        seriesnum3.push_back(0);
        seriesnum4.push_back(0);
        seriesnum5.push_back(0);
    }
    return;
}
void my_seriescount(vector<int> &seriesrun,
                    vector<int> &seriescount,
                    int runNo)
{
    int check_i;
    for (int i = 0; i < seriesrun.size(); i++)
    {
        // check if runnumber exist
        if (runNo == seriesrun[i])
        {
            check_i = i;
        }
    }
    seriescount[check_i] += 1;
    return;
}
#pragma endregion