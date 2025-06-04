// simple_test.cpp - 简单的测试程序来调试main函数
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include "IMSim/EQSource_CB14PCA.h"

using namespace std;

// 从原main函数复制的lamda_M函数
double lamda_M(double M)
{
    double lamda;
    lamda = pow(10.0, -(M - 3.0));
    return lamda;
}

int main()
{
    ifstream EQSourcefile;
	ifstream Siteinfofile;

    EQSourcefile.open("EQSource.txt", ios::in);
    Siteinfofile.open("SiteFile.txt", ios::in);

    // 检查文件是否成功打开
    if (!EQSourcefile.is_open()) {
        cout << "错误：无法打开 EQSource.txt 文件！" << endl;
        return -1;
    }
    if (!Siteinfofile.is_open()) {
        cout << "错误：无法打开 SiteFile.txt 文件！" << endl;
        return -1;
    }

	stringstream str2any;
	string temp;

	//0/1是否输出中位值
	bool ifmedian; 
	EQSourcefile >> ifmedian;
	cout << "读取到的 ifmedian: " << ifmedian << endl;

	//震级
	double M;
	EQSourcefile >> M;
	cout << "读取到的 M: " << M << endl;

	//模拟次数
	int N_sim;
	EQSourcefile >> N_sim;
	cout << "读取到的 N_sim: " << N_sim << endl;

	//种子
	int seed; 
	EQSourcefile >> seed;
	cout << "读取到的 seed: " << seed << endl;
    default_random_engine p(seed);

	//设置震源
	double lon_0 = 0; EQSourcefile >> lon_0;
	double lat_0 = 0; EQSourcefile >> lat_0;
	EQSource_CB14PCA eqs(lon_0, lat_0);
	{
		eqs.set_randomengine(&p);
		//eqs.set_lamda_M(lamda_M, 5, 8);
		//断裂面参数
		double W; EQSourcefile >> W; eqs.set_W(W);
		double length; EQSourcefile >> length; eqs.set_length(length);
		double RuptureNormal_x, RuptureNormal_y, RuptureNormal_z;
		EQSourcefile >> RuptureNormal_x >> RuptureNormal_y >> RuptureNormal_z;
		eqs.set_RuptureNormal(RuptureNormal_x, RuptureNormal_y, RuptureNormal_z);
		double lambda; EQSourcefile >> lambda; eqs.set_lambda(lambda);
		int Fhw; EQSourcefile >> Fhw; eqs.set_Fhw(Fhw);
		double Zhyp; EQSourcefile >> Zhyp; eqs.set_Zhyp(Zhyp);
		int region; EQSourcefile >> region; eqs.set_region(region);
		int nPCs; EQSourcefile >> nPCs; eqs.set_nPCs(nPCs);
	}

	//输出场地数据
	int ID; double lon, lat, elevation_km, T0, Vs30, Z25;
	while (getline(Siteinfofile, temp))
	{
		str2any.clear();
		str2any.str(temp);
		str2any >> ID >> lon >> lat >> elevation_km >> T0 >> Vs30 >> Z25;
		eqs.register_site(ID, lon, lat, elevation_km, T0, Vs30, Z25);
	}

	EQSourcefile.close();
	Siteinfofile.close();

	// 模拟IM分布
	vector<double> M_;
	for (size_t i = 0; i < N_sim; i++)
	{
		M_.push_back(M);
	}
	eqs.SimulateIntensities(M_, ifmedian);

	// 输出
	eqs.io_IM("IM sim.txt");
	// eqs.io_IM_AllT("IM sim allT.txt");
	eqs.io_XY("XY coord.txt");
    
    
    return 0;
}
