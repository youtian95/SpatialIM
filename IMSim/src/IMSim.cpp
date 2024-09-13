// IMSim.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <sstream>
#include "IMSim/EQSource_CB14PCA.h"

double lamda_M(double M)
{
	double lamda;
	lamda = pow(10.0, -(M - 3.0));
	return lamda;
}

/// <summary>
/// 根据震源和场地输入模拟空间相关的烈度分布
/// </summary>
/// <param name="argc"></param>
/// <param name="argv">
/// 依次为震源信息和场地信息的文件名，例如 EQSource.txt，SiteFile.txt
/// 
/// EQSource.txt 每行依次为(空格分开)：
/// ifmedian - 0/1是否输出中位值
/// M - 震级
/// N_sim - 次数
/// seed - int, 随机数种子
/// lon_0，lat_0 - 震中经纬度，°
/// W - down-dip width of the fault rupture plane. 
///		if unknown, input: 999
/// length - 断裂面水平方向长度
/// RuptureNormal_x, RuptureNormal_y, RuptureNormal_z - 断裂面朝上的法线方向（向东为x,向北为y,向上为z）
/// lambda - rake angle (degree) - average angle of slip measured in the plane of rupture
/// Fhw - hanging wall effect
///		= 1 for including
///		= 0 for excluding
/// Zhyp - (km) Hypocentral depth of the earthquake measured from sea level
///		if unknown, input: 999
/// region
///		= 0 for global (incl. Taiwan)
///		= 1 for California
///		= 3 for China or Turkey
///		= 4 for Italy
/// nPCs - IM相关性PCA方法模拟考虑的主成分阶数，推荐大于等于5
/// 
/// SiteFile.txt 每行为一个场地的数据，每一行空格分开依次为:
/// ID - int
/// lon - 经度
/// lat - 纬度
/// elevation_km - 高程
/// period1 - 基本周期
/// Vs30_mpers - 剪切波速
/// Z25_km - Depth to the 2.5 km/s shear-wave velocity horizon (km) 
///		if in California or Japan and Z2.5 is unknow, then input: 999）
/// </param>
/// <returns></returns>
int main(int argc, char* argv[])
{
	ifstream EQSourcefile;
	ifstream Siteinfofile;

	if (argc == 1) {
		EQSourcefile.open("EQSource.txt", ios::in);
		Siteinfofile.open("SiteFile.txt", ios::in);
	}
	else if (argc == 3) {
		EQSourcefile.open(argv[1], ios::in);
		Siteinfofile.open(argv[2], ios::in);
	}
	else
	{
		cout << "输入参数错误！" << endl;
		return 1;
	}

	stringstream str2any;
	string temp;

	//0/1是否输出中位值
	bool ifmedian; 
	EQSourcefile >> ifmedian;

	//震级
	double M;
	EQSourcefile >> M;

	//模拟次数
	int N_sim;
	EQSourcefile >> N_sim;

	//种子
	int seed; 
	EQSourcefile >> seed;
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
}
