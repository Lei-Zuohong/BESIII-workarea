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
double my_angle_boost(HepLorentzVector track1,
                      HepLorentzVector track2)
{
    HepLorentzVector track1_new = track1.boost(-0.011, 0, 0);
    HepLorentzVector track2_new = track2.boost(-0.011, 0, 0);
    double px1 = track1_new.px();
    double py1 = track1_new.py();
    double pz1 = track1_new.pz();
    double e1 = track1_new.e();
    double px2 = track2_new.px();
    double py2 = track2_new.py();
    double pz2 = track2_new.pz();
    double e2 = track2_new.e();
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
#pragma region 事例选择
HepLorentzVector *my_omega_3pi(HepLorentzVector *tracklist,
                               double m1,
                               double m2,
                               double m3)
{ //将八个track中的后六个重建为三个指定质量的粒子
    int combine[15][6] = {{1, 2, 3, 4, 5, 6},
                          {1, 2, 3, 5, 4, 6},
                          {1, 2, 3, 6, 4, 5},
                          {1, 3, 2, 4, 5, 6},
                          {1, 3, 2, 5, 4, 6},
                          {1, 3, 2, 6, 4, 5},
                          {1, 4, 3, 2, 5, 6},
                          {1, 4, 3, 5, 2, 6},
                          {1, 4, 3, 6, 2, 5},
                          {1, 5, 3, 4, 2, 6},
                          {1, 5, 3, 2, 4, 6},
                          {1, 5, 3, 6, 4, 2},
                          {1, 6, 3, 4, 5, 2},
                          {1, 6, 3, 5, 4, 2},
                          {1, 6, 3, 2, 4, 5}};
    double chisq = 9999;
    HepLorentzVector *new_tracklist;
    new_tracklist = my_newheplorentzvector(8);
    for (int i = 0; i < 15; i++)
    {
        double chisq_mpi01 = pow((tracklist[combine[i][0] + 1] + tracklist[combine[i][1] + 1]).m() - m1, 2);
        double chisq_mpi02 = pow((tracklist[combine[i][2] + 1] + tracklist[combine[i][3] + 1]).m() - m2, 2);
        double chisq_mpi03 = pow((tracklist[combine[i][4] + 1] + tracklist[combine[i][5] + 1]).m() - m3, 2);
        double chisq_mpi0 = (chisq_mpi01 + chisq_mpi02 + chisq_mpi03) / 3;
        if (chisq_mpi0 < chisq)
        {
            chisq = chisq_mpi0;
            if (rand() % 10000 < 5000)
            {
                new_tracklist[0] = tracklist[0];
                new_tracklist[1] = tracklist[1];
                new_tracklist[2] = tracklist[combine[i][0] + 1];
                new_tracklist[3] = tracklist[combine[i][1] + 1];
                new_tracklist[4] = tracklist[combine[i][2] + 1];
                new_tracklist[5] = tracklist[combine[i][3] + 1];
                new_tracklist[6] = tracklist[combine[i][4] + 1];
                new_tracklist[7] = tracklist[combine[i][5] + 1];
            }
            else
            {
                new_tracklist[0] = tracklist[0];
                new_tracklist[1] = tracklist[1];
                new_tracklist[2] = tracklist[combine[i][1] + 1];
                new_tracklist[3] = tracklist[combine[i][0] + 1];
                new_tracklist[4] = tracklist[combine[i][3] + 1];
                new_tracklist[5] = tracklist[combine[i][2] + 1];
                new_tracklist[6] = tracklist[combine[i][5] + 1];
                new_tracklist[7] = tracklist[combine[i][4] + 1];
            }
        }
    }
    return new_tracklist;
}
HepLorentzVector *my_omega_omega(HepLorentzVector *tracklist,
                                 double mass)
{ //将八个track中的前四个重建为一个指定质量的粒子
    int combine[3][3] = {{1, 2, 3},
                         {2, 3, 1},
                         {3, 1, 2}};
    double chisq = 9999;
    HepLorentzVector *new_tracklist;
    new_tracklist = my_newheplorentzvector(8);
    for (int i = 0; i < 3; i++)
    {
        double chisq_m = pow((tracklist[0] + tracklist[1] + tracklist[2 * combine[i][0]] + tracklist[2 * combine[i][0] + 1]).m() - mass, 2);
        if (chisq_m < chisq)
        {
            chisq = chisq_m;
            new_tracklist[0] = tracklist[0];
            new_tracklist[1] = tracklist[1];
            new_tracklist[2] = tracklist[2 * combine[i][0]];
            new_tracklist[3] = tracklist[2 * combine[i][0] + 1];
            new_tracklist[4] = tracklist[2 * combine[i][1]];
            new_tracklist[5] = tracklist[2 * combine[i][1] + 1];
            new_tracklist[6] = tracklist[2 * combine[i][2]];
            new_tracklist[7] = tracklist[2 * combine[i][2] + 1];
        }
    }
    return new_tracklist;
}
HepLorentzVector *my_omega_lower(HepLorentzVector *tracklist)
{ //将八个track中的后四个按照判据排列
    int combine[2][2] = {{1, 2},
                         {2, 1}};
    double chisq = 9999;
    HepLorentzVector *new_tracklist;
    new_tracklist = my_newheplorentzvector(8);
    for (int i = 0; i < 2; i++)
    {
        HepLorentzVector track_omega = tracklist[0] + tracklist[1] + tracklist[2] + tracklist[3];
        double chisq_m = (track_omega + tracklist[2 * combine[i][0] + 2] + tracklist[2 * combine[i][0] + 3]).m();
        if (chisq_m < chisq)
        {
            chisq = chisq_m;
            new_tracklist[0] = tracklist[0];
            new_tracklist[1] = tracklist[1];
            new_tracklist[2] = tracklist[2];
            new_tracklist[3] = tracklist[3];
            new_tracklist[4] = tracklist[2 * combine[i][0] + 2];
            new_tracklist[5] = tracklist[2 * combine[i][0] + 3];
            new_tracklist[6] = tracklist[2 * combine[i][1] + 2];
            new_tracklist[7] = tracklist[2 * combine[i][1] + 3];
        }
    }
    return new_tracklist;
}
#pragma endregion