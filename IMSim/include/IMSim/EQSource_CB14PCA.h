#pragma once

#include <random>
#include <iostream>
#include <windows.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "Eigen/Dense"
#include "Eigen/Geometry"
#include "GMPE.h"
#include "SimulateResiduals.h"
#include "ResidualCE.h"

using namespace std;
using namespace Eigen;

class EQSource_CB14PCA
	//震源
	//设置震源参数->所有需要模拟烈度的场地注册信息->统一模拟烈度结果
	//所有震源参数必须准确输入, 不能输入 999
	//ps: 内部x,y以震源中心为原点, km
{
public:
	struct SiteInfo
	{
		int ID;
		double x, y, T, Rrup, Rjb, Rx, Ztor, Zbot, Vs30, Z25;
	};
public:
	EQSource_CB14PCA()
	{}
	EQSource_CB14PCA(double lon, double lat) :
		location(pair<double, double>(lon, lat))
	{
	}
	void set_randomengine(default_random_engine* p)
	{
		pRND = p;
		SR.set_random_engine(p);
	}
	void set_location(double lon, double lat)
	{
		location = make_pair(lon, lat);
	}
	void set_W(double w)
		//W - down - dip width of the fault rupture plane
		// if unknown, input: 999
	{
		W = w;
	}
	void set_length(double in)
	{
		length = in;
	}
	void set_RuptureNormal(double x, double y, double z)
		//设置断裂面朝上的法向向量
	{
		RuptureNormal << x, y, z;
		RuptureNormal.normalize();
		Vector3f z_unit(0, 0, 1);
		delta = acos(abs(RuptureNormal.dot(z_unit))) / pi * 180;
		assert(delta >= -0.1 && delta <= 90.1);
	}
	void set_lambda(double in)
		//rake angle(degree) - average angle of slip measured in
		//	the plance of rupture
	{
		lambda = in;
	}
	void set_Fhw(bool in)
		// hanging wall effect
		//	 = 1 for including
		//	 = 0 for excluding
	{
		Fhw = in;
	}
	void set_Zhyp(double in)
		// (km)Hypocentral depth of the earthquake measured from sea level
		//		if unknown, input: 999
	{
		Zhyp = in;
	}
	void set_region(int in)
		//= 0 for global(incl.Taiwan)
		//= 1 for California
		//= 2 for Japan
		//= 3 for China or Turkey
		//= 4 for Italy
	{
		region = in;
	}
	void set_nPCs(int in)
		//IM相关性PCA方法模拟考虑的主成分阶数，推荐大于等于5
	{
		nPCs = in;
	}
	void set_lamda_M(double (*f)(double M), double M_min_, double M_max_)
		//设置函数指针
	{
		assert(M_min_ < M_max_);
		lamda_M = f;
		M_min = M_min_;
		M_max = M_max_;
	}
	//注册场地和模拟烈度
	void register_site(int ID, double lon, double lat, double elevation_km,
		double period1, double Vs30_mpers, double Z25_km)
	{
		int ID_;
		double x_, y_, T_, Rrup_, Rjb_, Rx_, Ztor_, Zbot_, Vs30_, Z25_;
		{
			ID_ = ID;
			x_ = getX(lon);
			y_ = getY(lat);
			T_ = period1;
			Rrup_ = Rrup(lon, lat, elevation_km);
			Rjb_ = Rjb(lon, lat);
			Rx_ = Rx(lon, lat);
			Ztor_ = elevation_km + Zhyp - W * sin(delta / 180.0 * pi) / 2.0;
			Zbot_ = Ztor_ + W * sin(delta / 180.0 * pi);
			Vs30_ = Vs30_mpers;
			Z25_ = Z25_km;
		}
		SiteInfo Site = { ID_, x_, y_, T_, Rrup_, Rjb_, Rx_, Ztor_, Zbot_, Vs30_, Z25_ };
		SiteInfo_map.insert(make_pair(ID_, Site));
		firstSim = 1;
	}
	void clear_sites()
		//清除所有注册场地
	{
		SiteInfo_map.clear();
		firstSim = 1;
	}
	void SimulateIntensities(vector<double> M_, bool ifmedian = false, bool OutputIMAllPeriods = true)
	{
		//中间结果，所有模拟的周期
		VectorXf T_sim(19);
		T_sim << 0.010, 0.020, 0.030, 0.050, 0.075, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50,
			0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0;

		PGA_map.clear();
		Sa_map.clear();
		Sa_map_allPeriods.clear();
		SaMedian_map_allPeriods.clear();
		SaMedian_PlusBR_map_allPeriods.clear();
		
		for (auto& iter : SiteInfo_map)
		{
			PGA_map.insert(make_pair(iter.first, VectorXf::Zero(M_.size())));
			Sa_map.insert(make_pair(iter.first, VectorXf::Zero(M_.size())));
			vector<VectorXf> sa_allperiod_eachsite;
			vector<double> saMedian_allperiod_eachsite;
			for (size_t i = 0; i < T_sim.size(); i++)
			{
				sa_allperiod_eachsite.push_back(VectorXf::Zero(M_.size()));
				saMedian_allperiod_eachsite.push_back(0.0);
			}
			Sa_map_allPeriods.insert(make_pair(iter.first, sa_allperiod_eachsite));
			SaMedian_map_allPeriods.insert(make_pair(iter.first, saMedian_allperiod_eachsite));
			SaMedian_PlusBR_map_allPeriods.insert(make_pair(iter.first, sa_allperiod_eachsite));
		}

		vector<MatrixXf> IMNormResiduals; // n_sim <n_loc x n_T>
		MatrixXf IMNormInterResiduals; // n_sim x n_T
		if (!ifmedian)
		{
			VectorXf x(SiteInfo_map.size()), y(SiteInfo_map.size());
			int j_site = 0;
			for (map<int, SiteInfo>::iterator iter = SiteInfo_map.begin();
				iter != SiteInfo_map.end(); iter++)
			{
				x[j_site] = iter->second.x;
				y[j_site] = iter->second.y;
				j_site++;
			}
			if (firstSim)
				//第一次模拟需要注册场地
			{
				SR.RegisterSites(x, y);
				firstSim = 0;
			}
			//模拟 within-event residuals
			IMNormResiduals = SR.simulateResiduals(T_sim, M_.size());
			//模拟 between-event residuals
			IMNormInterResiduals = ResidualCE::Inter_Event_Residuals_Simulation(T_sim, M_.size(), pRND);
		}

		std::normal_distribution<> dist;

		GMPE gmpe;

		int j = 0;
		for (map<int, SiteInfo>::iterator iter = SiteInfo_map.begin();
			iter != SiteInfo_map.end(); iter++)
			// 每个场地的建筑
		{
			auto site = iter->second;

			//每个场地缓存一个结果，如果震级不变，可以直接用此结果
			VectorXf Median_PGA_buffer, sigma_PGA_buffer, tau_PGA_buffer,
				Median_Sa_buffer, sigma_Sa_buffer, tau_Sa_buffer, 
				Median_Sa_allT_buffer, sigma_Sa_allT_buffer, tau_Sa_allT_buffer;
			double M_buffer = 0;

			for (int i = 0; i < M_.size(); i++)
				// 每次地震
			{

				VectorXf Median_PGA, sigma_PGA, tau_PGA, period1;
				VectorXf Median_Sa, sigma_Sa, tau_Sa;
				VectorXf Median_Sa_allT, sigma_Sa_allT, tau_Sa_allT;


				if (M_[i] != M_buffer)
				{
					// 用0s处的Sa代表PGA
					gmpe.CB_2014_nga(&Median_PGA, &sigma_PGA, &tau_PGA, &period1,
						M_[i], VectorXf::Zero(1), site.Rrup, site.Rjb, site.Rx,
						W, site.Ztor, site.Zbot, delta, lambda, Fhw,
						site.Vs30, site.Z25, Zhyp, region);
					// Sa
					VectorXf T(1); T << site.T;
					gmpe.CB_2014_nga(&Median_Sa, &sigma_Sa, &tau_Sa, &period1,
						M_[i], T, site.Rrup, site.Rjb, site.Rx,
						W, site.Ztor, site.Zbot, delta, lambda, Fhw,
						site.Vs30, site.Z25, Zhyp, region);
					//Sa all T
					gmpe.CB_2014_nga(&Median_Sa_allT, &sigma_Sa_allT, &tau_Sa_allT, &period1,
						M_[i], T_sim, site.Rrup, site.Rjb, site.Rx,
						W, site.Ztor, site.Zbot, delta, lambda, Fhw,
						site.Vs30, site.Z25, Zhyp, region);

					//修改缓存
					M_buffer = M_[i];
					Median_PGA_buffer = Median_PGA;
					sigma_PGA_buffer = sigma_PGA;
					tau_PGA_buffer = tau_PGA;
					Median_Sa_buffer = Median_Sa;
					sigma_Sa_buffer = sigma_Sa;
					tau_Sa_buffer = tau_Sa;
					Median_Sa_allT_buffer = Median_Sa_allT;
					sigma_Sa_allT_buffer = sigma_Sa_allT;
					tau_Sa_allT_buffer = tau_Sa_allT;
				}
				else
				{
					Median_PGA = Median_PGA_buffer;
					sigma_PGA = sigma_PGA_buffer;
					tau_PGA = tau_PGA_buffer;
					Median_Sa = Median_Sa_buffer;
					sigma_Sa = sigma_Sa_buffer;
					tau_Sa = tau_Sa_buffer;
					Median_Sa_allT = Median_Sa_allT_buffer;
					sigma_Sa_allT = sigma_Sa_allT_buffer;
					tau_Sa_allT = tau_Sa_allT_buffer;
				}

				// 根据相关性模拟结果
				if (!ifmedian)
				{
					// within-event residuals
					float PGA_Res = GMPE::interp1(T_sim, IMNormResiduals[i].row(j), 0.001);
					// between-event residuals
					float PGA_Inter_Res = GMPE::interp1(T_sim, IMNormInterResiduals.row(i), 0.001);
					// 最后的模拟结果
					float PGA = exp(log(Median_PGA[0]) + PGA_Inter_Res * tau_PGA[0]
						+ sqrt(sigma_PGA[0] * sigma_PGA[0] - tau_PGA[0] * tau_PGA[0]) * PGA_Res);
					// within-event residuals
					float Sa_Res = GMPE::interp1(T_sim, IMNormResiduals[i].row(j), site.T);
					// between-event residuals
					float Sa_Inter_Res = GMPE::interp1(T_sim, IMNormInterResiduals.row(i), site.T);
					// 最后的模拟结果
					float Sa = exp(log(Median_Sa[0]) + Sa_Inter_Res * tau_Sa[0]
						+ sqrt(sigma_Sa[0] * sigma_Sa[0] - tau_Sa[0] * tau_Sa[0]) * Sa_Res);
					// 所有周期的结果
					for (size_t i_T_sim = 0; i_T_sim < T_sim.size(); i_T_sim++)
					{
						float Sa_T_sim = exp(log(Median_Sa_allT[i_T_sim]) 
							+ IMNormInterResiduals.row(i)[i_T_sim] * tau_Sa_allT[i_T_sim]
							+ sqrt(sigma_Sa_allT[i_T_sim] * sigma_Sa_allT[i_T_sim] 
							- tau_Sa_allT[i_T_sim] * tau_Sa_allT[i_T_sim])
							* IMNormResiduals[i].row(j)[i_T_sim]);
						Sa_map_allPeriods[site.ID][i_T_sim][i] = Sa_T_sim;
						SaMedian_map_allPeriods[site.ID][i_T_sim] = Median_Sa_allT[i_T_sim];
						SaMedian_PlusBR_map_allPeriods[site.ID][i_T_sim][i] =
							exp(log(Median_Sa_allT[i_T_sim])
								+ IMNormInterResiduals.row(i)[i_T_sim] * tau_Sa_allT[i_T_sim]);
					}
					//赋值
					PGA_map[site.ID][i] = PGA;
					Sa_map[site.ID][i] = Sa;
				}
				else
					// 输出中值
				{
					PGA_map[site.ID][i] = Median_PGA[0];
					Sa_map[site.ID][i] = Median_Sa[0];
				}
			}

			j++;
		}

		if (OutputIMAllPeriods && (!ifmedian))
		{
			for (size_t i_T_sim = 0; i_T_sim < T_sim.size(); i_T_sim++)
			{
				ofstream of_Sa;
				ofstream of_SaMedian;
				ofstream of_SaMedian_PlusBR;
				string period_temp = to_string(T_sim[i_T_sim]);
				of_Sa.open("IM sim with period "
					+ period_temp.erase(period_temp.find_last_not_of('0') + 1, std::string::npos)
					+ ".txt", ios::out | ios::trunc);
				of_SaMedian.open("IM median with period "
					+ period_temp.erase(period_temp.find_last_not_of('0') + 1, std::string::npos)
					+ ".txt", ios::out | ios::trunc);
				of_SaMedian_PlusBR.open("IM median plus Between-event residual with period "
					+ period_temp.erase(period_temp.find_last_not_of('0') + 1, std::string::npos)
					+ ".txt", ios::out | ios::trunc);

				auto iter = Sa_map_allPeriods.begin();
				while (iter != Sa_map_allPeriods.end()) {

					int ID = iter->first;
					of_Sa << ID << "  ";
					VectorXf SaVec = iter->second[i_T_sim];
					for (size_t i = 0; i < SaVec.size(); i++)
					{
						of_Sa << SaVec[i] << "  ";
					}
					of_Sa << endl;
					iter++;
				}

				auto iter_SaMedian = SaMedian_map_allPeriods.begin();
				while (iter_SaMedian != SaMedian_map_allPeriods.end()) {

					int ID = iter_SaMedian->first;
					of_SaMedian << ID << "  ";
					double SaMedian = iter_SaMedian->second[i_T_sim];
					of_SaMedian << SaMedian;
					of_SaMedian << endl;
					iter_SaMedian++;
				}

				auto iter_SaMedian_PlusBR = SaMedian_PlusBR_map_allPeriods.begin();
				while (iter_SaMedian_PlusBR != SaMedian_PlusBR_map_allPeriods.end()) {

					int ID = iter_SaMedian_PlusBR->first;
					of_SaMedian_PlusBR << ID << "  ";
					VectorXf SaVec = iter_SaMedian_PlusBR->second[i_T_sim];
					for (size_t i = 0; i < SaVec.size(); i++)
					{
						of_SaMedian_PlusBR << SaVec[i] << "  ";
					}
					of_SaMedian_PlusBR << endl;
					iter_SaMedian_PlusBR++;
				}

				of_Sa.close();
				of_SaMedian.close();
				of_SaMedian_PlusBR.close();
			}
		}
	}
	const VectorXf& get_PGA(int ID) const
	{
		return PGA_map.at(ID);
	}
	const VectorXf& get_Sa(int ID) const
	{
		return Sa_map.at(ID);
	}
	//按场地输出结果
	void io_IM(string filename) {
		// 每一行为一个场地的结果，依次为 [ID, 1 x N_sim]

		ofstream of;
		of.open(filename, ios::out | ios::trunc);
		auto iter = Sa_map.begin();
		while (iter != Sa_map.end()) {
			int ID = iter->first;
			of << ID << "  ";
			VectorXf SaVec = iter->second;
			for (size_t i = 0; i < SaVec.size(); i++)
			{
				of << SaVec[i] << "  ";
			}
			of << endl;
			iter++;
		}

		of.close();
	}
	void io_IM_AllT(string filename) {
		// 每一行为一个场地的结果，依次为 [ID,T,1 x N_sim], 对于同一个ID场地，T有几个周期就有几行

		//中间结果，所有模拟的周期
		VectorXf T_sim(19);
		T_sim << 0.010, 0.020, 0.030, 0.050, 0.075, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50,
			0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0;

		ofstream of;
		of.open(filename, ios::out | ios::trunc);
		auto iter = Sa_map_allPeriods.begin();
		// 每个场地
		while (iter != Sa_map_allPeriods.end()) {
			int ID = iter->first;
			// T=0的PGA
			of << ID << "  ";
			VectorXf PGA_1site = PGA_map[ID];
			for (size_t i = 0; i < PGA_1site.size(); i++)
			{
				of << PGA_1site[i] << "  ";
			}
			of << endl;
			// 每个周期的SA
			vector<VectorXf> SaVec_allT = iter->second;
			for (size_t i_T = 0; i_T < SaVec_allT.size(); i_T++)
			{
				of << ID << "  ";
				of << T_sim[i_T] << "  ";
				// 每次模拟
				VectorXf SaVec = SaVec_allT[i_T];
				for (size_t i = 0; i < SaVec.size(); i++)
				{
					of << SaVec[i] << "  ";
				}
				of << endl;
			}
			iter++;
		}

		of.close();
	}
	void io_XY(string filename) {
		// 每一行为一个场地的相对于震中的坐标，依次为 [ID, x ,y]

		ofstream of;
		of.open(filename, ios::out | ios::trunc);
		auto iter = SiteInfo_map.begin();
		while (iter != SiteInfo_map.end()) {
			int ID = iter->first;
			of << ID << "  ";
			SiteInfo Siteinfo = iter->second;
			of << Siteinfo.x << "  " << Siteinfo.y << "  ";
			of << endl;
			iter++;
		}

		of.close();
	}
	//按网格输出
	void io_PGA_Sa_Grid(vector<double> ML, bool ifmedian, double period,
		double elevation_km, double Vs30, double Z25_km,
		double bound_left, double bound_right, //经度
		double bound_down, double bound_up, //纬度
		int N_lat = 100, int N_lon = 100)
		//输出PGA Sa, 历史记录的场地信息和模拟结果会被消除
		//每一行是同一个纬度，经度和纬度均为从小到大排序
	{
		//注册场地
		SiteInfo_map.clear();
		int ID = 0;
		map<int, pair<double, double>> ID2XY;
		for (int n_lat = 0; n_lat <= N_lat; n_lat++)
		{
			for (int n_lon = 0; n_lon <= N_lon; n_lon++)
			{
				ID++;
				//当前点经度和纬度
				double lon = bound_left + n_lon * (bound_right - bound_left) / (double)N_lon;
				double lat = bound_down + n_lat * (bound_up - bound_down) / (double)N_lat;
				ID2XY.insert(make_pair(ID, make_pair(getX(lon), getY(lat))));
				register_site(ID, lon, lat, elevation_km, period, Vs30, Z25_km);
			}
		}

		//输出网格的x,y坐标
		cout << "输出网格x,y坐标... " << endl;
		{
			ofstream of;
			of.open("X_grid.txt", ios::out | ios::trunc);
			ID = 0;
			for (int n_lat = 0; n_lat <= N_lat; n_lat++)
			{
				for (int n_lon = 0; n_lon <= N_lon; n_lon++)
				{
					ID++;
					of << ID2XY.at(ID).first;
					if (n_lon < N_lon)
					{
						of << ",";
					}
					else
					{
						of << endl;
					}
				}
			}
			of.close();
		}
		{
			ofstream of;
			of.open("Y_grid.txt", ios::out | ios::trunc);
			ID = 0;
			for (int n_lat = 0; n_lat <= N_lat; n_lat++)
			{
				for (int n_lon = 0; n_lon <= N_lon; n_lon++)
				{
					ID++;
					of << ID2XY.at(ID).second;
					if (n_lon < N_lon)
					{
						of << ",";
					}
					else
					{
						of << endl;
					}
				}
			}
			of.close();
		}

		//输出网格的经纬度坐标
		cout << "输出网格经纬度坐标... " << endl;
		{
			ofstream of;
			of.open("longitude_grid.txt", ios::out | ios::trunc);
			ID = 0;
			for (int n_lat = 0; n_lat <= N_lat; n_lat++)
			{
				for (int n_lon = 0; n_lon <= N_lon; n_lon++)
				{
					ID++;
					//当前点经度和纬度
					double lon = bound_left + n_lon * (bound_right - bound_left) / (double)N_lon;
					double lat = bound_down + n_lat * (bound_up - bound_down) / (double)N_lat;
					of << fixed << setprecision(6) << lon;
					if (n_lon < N_lon)
					{
						of << ",";
					}
					else
					{
						of << endl;
					}
				}
			}
			of.close();
		}
		{
			ofstream of;
			of.open("latitude_grid.txt", ios::out | ios::trunc);
			ID = 0;
			for (int n_lat = 0; n_lat <= N_lat; n_lat++)
			{
				for (int n_lon = 0; n_lon <= N_lon; n_lon++)
				{
					ID++;
					//当前点经度和纬度
					double lon = bound_left + n_lon * (bound_right - bound_left) / N_lon;
					double lat = bound_down + n_lat * (bound_up - bound_down) / N_lat;
					of << fixed << setprecision(6) << lat;
					if (n_lon < N_lon)
					{
						of << ",";
					}
					else
					{
						of << endl;
					}
				}
			}
			of.close();
		}

		//模拟地震
		M_list = ML;
		//simulate_intensities(M_list, ifmedian);
		SimulateIntensities(M_list, ifmedian);

		//输出模拟的PGA
		cout << "输出PGA网格... " << endl;
		{
			vector<ofstream> outfile(M_list.size());
			for (int i = 0; i < outfile.size(); i++)
			{
				string filename = "PGA_grid_No_" + to_string(i) + ".txt";
				outfile[i].open(filename, ios::out | ios::trunc);
			}
			ID = 0;
			for (int n_lat = 0; n_lat <= N_lat; n_lat++)
			{
				for (int n_lon = 0; n_lon <= N_lon; n_lon++)
				{
					ID++;
					auto& IM_vec = get_PGA(ID);
					for (int i = 0; i < outfile.size(); i++)
					{
						outfile[i] << fixed << setprecision(4) << IM_vec[i];
						if (n_lon < N_lon)
						{
							outfile[i] << ",";
						}
						else
						{
							outfile[i] << endl;
						}
					}
				}
			}
			for (int i = 0; i < outfile.size(); i++)
			{
				outfile[i].close();
			}
		}

		//输出模拟的Sa
		cout << "输出Sa网格... " << endl;
		{
			vector<ofstream> outfile(M_list.size());
			for (int i = 0; i < outfile.size(); i++)
			{
				string filename = "Sa_grid_No_" + to_string(i) + ".txt";
				outfile[i].open(filename, ios::out | ios::trunc);
			}
			ID = 0;
			for (int n_lat = 0; n_lat <= N_lat; n_lat++)
			{
				for (int n_lon = 0; n_lon <= N_lon; n_lon++)
				{
					ID++;
					auto& IM_vec = get_Sa(ID);
					for (int i = 0; i < outfile.size(); i++)
					{
						outfile[i] << fixed << setprecision(4) << IM_vec[i];
						if (n_lon < N_lon)
						{
							outfile[i] << ",";
						}
						else
						{
							outfile[i] << endl;
						}
					}
				}
			}
			for (int i = 0; i < outfile.size(); i++)
			{
				outfile[i].close();
			}
		}

	}
	//按网格输出场地参数
	void io_Rrup(string filename, double elevation_km,
		double bound_left, double bound_right, //经度
		double bound_down, double bound_up, //纬度
		int N_lat = 100, int N_lon = 100) const
		//每一行是同一个纬度，经度和纬度均为从小到大排序
	{
		cout << "输出Rrup网格: " << filename << " ..." << endl;
		ofstream outfile;
		outfile.open(filename, ios::out | ios::trunc);

		for (int n_lat = 0; n_lat <= N_lat; n_lat++)
		{
			for (int n_lon = 0; n_lon < N_lon; n_lon++)
			{
				//当前点经度和纬度
				double lon = bound_left + n_lon * (bound_right - bound_left) / N_lon;
				double lat = bound_down + n_lat * (bound_up - bound_down) / N_lat;
				//输出
				outfile << Rrup(lon, lat, elevation_km) << ",";
			}
			double lat = bound_down + n_lat * (bound_up - bound_down) / N_lat;
			outfile << Rrup(bound_right, lat, elevation_km) << endl;
		}

		outfile.close();
	}
	void io_Rx(string filename,
		double bound_left, double bound_right, //经度
		double bound_down, double bound_up, //纬度
		int N_lat = 100, int N_lon = 100) const
		//每一行是同一个纬度，经度和纬度均为从小到大排序
	{
		cout << "输出Rx网格: " << filename << " ..." << endl;
		ofstream outfile;
		outfile.open(filename, ios::out | ios::trunc);

		for (int n_lat = 0; n_lat <= N_lat; n_lat++)
		{
			for (int n_lon = 0; n_lon < N_lon; n_lon++)
			{
				//当前点经度和纬度
				double lon = bound_left + n_lon * (bound_right - bound_left) / N_lon;
				double lat = bound_down + n_lat * (bound_up - bound_down) / N_lat;
				//输出

				outfile << Rx(lon, lat) << ",";
			}
			double lat = bound_down + n_lat * (bound_up - bound_down) / N_lat;
			outfile << Rx(bound_right, lat) << endl;
		}

		outfile.close();
	}
	void io_Rjb(string filename,
		double bound_left, double bound_right, //经度
		double bound_down, double bound_up, //纬度
		int N_lat = 100, int N_lon = 100) const
		//每一行是同一个纬度，经度和纬度均为从小到大排序
	{
		cout << "输出Rjb网格: " << filename << " ..." << endl;
		ofstream outfile;
		outfile.open(filename, ios::out | ios::trunc);

		for (int n_lat = 0; n_lat <= N_lat; n_lat++)
		{
			for (int n_lon = 0; n_lon < N_lon; n_lon++)
			{
				//当前点经度和纬度
				double lon = bound_left + n_lon * (bound_right - bound_left) / N_lon;
				double lat = bound_down + n_lat * (bound_up - bound_down) / N_lat;
				//输出
				outfile << Rjb(lon, lat) << ",";
			}
			double lat = bound_down + n_lat * (bound_up - bound_down) / N_lat;
			outfile << Rjb(bound_right, lat) << endl;
		}

		outfile.close();
	}
	//功能函数
	pair<double, double> get_location() const
	{
		return location;
	}
	double get_lamda_all() const
		//一年平均发生次数
	{
		if (!lamda_M)
		{
			cout << "Error: 没有输入震源概率信息." << endl;
			throw "";
		}
		return (lamda_M(M_min) - lamda_M(M_max));
	}
	double get_RandomMagnitude(default_random_engine* p) const
		//随机生成一个震级
	{
		if (!lamda_M)
		{
			cout << "Error: 没有输入震源概率信息." << endl;
			throw "";
		}

		//参数
		double error = 0.0005; int n_max = 10;

		double lamda_max = lamda_M(M_min);
		double lamda_min = lamda_M(M_max);
		uniform_real_distribution<double> dis(0, lamda_max);
		double lamda_random = dis(*p);
		if (lamda_random <= lamda_min)
		{
			return M_max;
		}
		else
		{
			return Finverse(lamda_random, lamda_M, error, n_max, M_min, M_max);
		}
	}
	vector<double> get_DiscreteMagnitude(int n_M) const
		//根据n_M划分震级范围, 返回震级向量
	{
		if (!lamda_M)
		{
			cout << "Error: 没有输入震源概率信息." << endl;
			throw "";
		}
		assert(n_M >= 2);
		vector<double> result(n_M);
		for (int i = 0; i < n_M; ++i)
		{
			result[i] = M_min + (M_max - M_min) / ((double)(n_M - 1)) * i;
		}
		return result;
	}
	double get_Integral_f_M_X_x(const map<double, double>& x_M) const
		//计算 Integral[f(M)x(M)]dM 积分，f(M)为M的概率密度函数, 
		//x(M)为外部输入函数 map<M, x> 必须从小到大
	{
		double integral = 0;
		if (x_M.size() < 2) throw "Error!";

		map<double, double>::const_iterator iter, iter_next, iter_last, iter_final;
		iter_final = x_M.end(); iter_final--;
		iter = x_M.begin();
		for (iter = x_M.begin(); iter != x_M.end(); iter++) {
			iter_next = iter; iter_last = iter;
			if (iter == x_M.begin())
			{
				iter_last = iter;
			}
			else
			{
				iter_last--;
			}
			if (iter == iter_final)
			{
				iter_next = iter;
			}
			else
			{
				iter_next++;
			}

			double M = iter->first;
			double x = iter->second;
			double f_M_dM;
			if (iter == x_M.begin())
			{
				f_M_dM = (get_lamdaFromM(M_min) - get_lamdaFromM((M + iter_next->first) / 2.0)) / get_lamdaFromM(M_min);
			}
			else if (iter == iter_final)
			{
				f_M_dM = get_lamdaFromM((M + iter_last->first) / 2.0) / get_lamdaFromM(M_min);
			}
			else
			{
				double M_1 = (M + iter_last->first) / 2.0;
				double lamda1 = get_lamdaFromM(M_1);
				double M_2 = (M + iter_next->first) / 2.0;
				double lamda2 = get_lamdaFromM(M_2);
				f_M_dM = (lamda1 - lamda2) / get_lamdaFromM(M_min);
			}
			integral += f_M_dM * x;
		}
		return integral;
	}
	double get_lamdaFromM(double M) const
		//查询lamda函数
	{
		if (!lamda_M)
		{
			cout << "Error: 没有输入震源概率信息." << endl;
			throw "";
		}
		return lamda_M(M);
	}
	double getX(double lon) const
		//根据经纬度获取某一点的X
	{
		double x_in = (lon - location.first) * pi / 180.0 * R * cos(location.second * pi / 180.0) / 1000.0;
		return x_in;
	}
	double getY(double lat) const
		//根据经纬度获取某一点的Y
	{
		double y_in = (lat - location.second) * pi / 180.0 * R / 1000.0;
		return y_in;
	}
	double getLon(double x) const
		//根据X获取某一点的经度
	{
		double lon = x * 1000.0 / cos(location.second * pi / 180.0) / R * 180.0 / pi + location.first;
		return lon;
	}
	double getLat(double y) const
		//根据y获取某一点的纬度
	{
		double lat = y * 1000.0 / R * 180.0 / pi + location.second;
		return lat;
	}

	static double Finverse(double y, double (*f)(double), double error, int n, double x_min, double x_max)
		//二分法求单调函数f(x)的反函数
		//y为计算的值,f为函数指针, error为误差, 最多计算n次
	{
		assert(error > 0);
		assert(x_min < x_max);
		if (f(x_min) == f(x_max)) return x_min; //相等时, 输出最小x
		if (f(x_min) <= f(x_max))
			//单增
		{
			assert(y >= f(x_min) && y <= f(x_max));
		}
		else
			//单减
		{
			assert(y >= f(x_max) && y <= f(x_min));
		}

		double x_left = x_min;
		double x_right = x_max;
		double x = (x_left + x_right) / 2.0;
		int i = 0;
		while (abs(f(x) - y) > error && i < n)
			//误差没达到或者计算次数没达到
		{
			if ((y - f(x_left)) * (f(x) - y) >= 0)
				//在x_left和x之间
			{
				x_right = x;
				x = (x_left + x_right) / 2.0;
			}
			else
				//在x_right和x之间
			{
				x_left = x;
				x = (x_left + x_right) / 2.0;
			}
			i++;
		}
		return x;
	}
	static int RecentM(double M, const vector<double>& M_list)
		//返回M_list中M对应的索引位置
		//之外返回-1
	{
		//中间
		for (int i = 1; i < M_list.size() - 1; i++)
		{
			if (((M_list[i - 1] + M_list[i]) / 2.0) <= M &&
				M <= ((M_list[i + 1] + M_list[i]) / 2.0))
			{
				return i;
			}
		}
		//i=0
		if (M_list[0] <= M &&
			M <= ((M_list[0] + M_list[1]) / 2.0))
		{
			return 0;
		}
		//i=end
		auto iter = M_list.end(); iter--; iter--;
		if (*M_list.rbegin() >= M &&
			M >= ((*iter + *M_list.rbegin()) / 2.0))
		{
			return (M_list.size() - 1);
		}
		return -1;
	}

protected:
	default_random_engine* pRND;
	//是否第一次模拟相关性
	bool firstSim = 1;
	//震源信息
	pair<double, double> location; //经纬度
	double R = 6371393; //地球半径m
	double pi = 3.1415926;
	double W = 0;
	double length = 0; //断裂面水平方向长度
	Vector3f RuptureNormal; 
	double delta = 90; //degree
	double lambda = 0; //degree
	bool Fhw = 1;
	double Zhyp = 0;
	int region = 3;
	int nPCs = 18;
	//场地信息
	map<int, SiteInfo> SiteInfo_map;
	//模拟地震
	vector<double> M_list;
	double (*lamda_M)(double M) = 0; //震级为M年均超越次数
	double M_min; double M_max; //考虑的震级上下限
	SimulateResiduals SR;
	//模拟结果
	map<int, VectorXf> PGA_map; //每次模拟的结果
	map<int, VectorXf> Sa_map; //每次模拟的结果
	map<int, vector<VectorXf>> Sa_map_allPeriods; //所有周期的Sa模拟结果
	map<int, vector<double>> SaMedian_map_allPeriods; //所有周期的Sa中值模拟结果
	map<int, vector<VectorXf>> SaMedian_PlusBR_map_allPeriods; //所有周期的 Sa中值+Between-event residul 模拟结果

	static bool FootpointInsideQuad(Vector3f p, vector<Vector3f> quad)
		//判断p的垂足是否在矩形之内,计算垂足在矩形平面的相对坐标, 矩形坐标的顺序必须连续
	{
		//矩形两条边
		Vector3f P12 = quad[1] - quad[0];
		Vector3f P23 = quad[2] - quad[1];
		//矩形中心
		Vector3f c = 0.25 * (quad[0] + quad[1] + quad[2] + quad[3]);
		//垂足
		Vector3f Normal = P12.cross(P23);
		Normal.normalize();
		if ((c - p).dot(Normal) < 0) Normal = -Normal;
		Vector3f fp = abs((c - p).dot(Normal)) * Normal + p;
		//矩形平面内的直线
		Vector3f l = fp - c;
		//相对坐标
		double x = abs(l.dot(P12.normalized()) / P12.norm());
		double y = abs(l.dot(P23.normalized()) / P23.norm());

		if (x > 0.5 || y > 0.5)
		{
			return false;
		}
		else
		{
			return true;
		}
	}
	static float Distance_PointAndLineSegment(Vector3f p, Vector3f p_1, Vector3f p_2)
		// return: 点到线段的最小距离
	{
		double dis;
		if ((p - p_1).dot(p_2 - p_1) < 0 || (p - p_2).dot(p_1 - p_2) < 0)
		{
			dis = min((p - p_1).norm(), (p - p_2).norm());
		}
		else
		{
			dis = (p_1 - p).cross(p_2 - p).norm() / (p_1 - p_2).norm();
		}
		return dis;
	}

	vector<Vector3f> Rupture4Points() const
		//断裂面的四个点, 上->下, 相邻排列
	{
		vector<Vector3f> result;

		Vector3f Z_unit(0, 0, 1);
		Vector3f Hor_unit = Z_unit.cross(RuptureNormal);
		Hor_unit.normalize();
		Vector3f Hor_half = length / 2.0 * Hor_unit;
		Vector3f ramp_up_half = RuptureNormal.cross(Hor_half);
		ramp_up_half.normalize();
		ramp_up_half = W / 2.0 * ramp_up_half;
		//四个点
		Vector3f P1 = ramp_up_half + Hor_half; 
		result.push_back(P1);
		Vector3f P2 = P1 - 2.0 * Hor_half;
		result.push_back(P2);
		Vector3f P3 = P2 - 2.0 * ramp_up_half;
		result.push_back(P3);
		Vector3f P4 = P3 + 2.0 * Hor_half;
		result.push_back(P4);

		return result;
	}
	double Rrup(double lon, double lat, double elevation_km) const
		//Rrup距离, 断裂面为矩形
		//到断裂面的最小距离
	{
		double x = getX(lon);
		double y = getY(lat);
		double z = elevation_km + Zhyp;
		Vector3f site(x,y,z);
		//断裂面矩形四个点坐标, 以断裂面中心为原点
		vector<Vector3f> Rupture4P = Rupture4Points();
		
		if (FootpointInsideQuad(site, Rupture4P))
		{
			//垂线距离
			double dist1 = (Rupture4P[0] - site).dot(-RuptureNormal);
			dist1 = abs(dist1);
			return dist1;
		}
		else
			//计算到四条边的最近距离
		{
			float dist1 = Distance_PointAndLineSegment(site, Rupture4P[0], Rupture4P[1]);
			float dist2 = Distance_PointAndLineSegment(site, Rupture4P[1], Rupture4P[2]);
			float dist3 = Distance_PointAndLineSegment(site, Rupture4P[2], Rupture4P[3]);
			float dist4 = Distance_PointAndLineSegment(site, Rupture4P[3], Rupture4P[0]);
			return min(min(dist1, dist2), min(dist3, dist4));
		}
	}
	double Rjb(double lon, double lat) const
		//到断裂面投影的最近距离
	{
		double x = getX(lon);
		double y = getY(lat);
		double z = 0;
		Vector3f site(x, y, z);
		//断裂面矩形四个点坐标, 以断裂面中心为原点
		vector<Vector3f> Rupture4P = Rupture4Points();
		for (int i = 0; i < Rupture4P.size(); i++)
		{
			Rupture4P[i][2] = 0;
		}

		if (FootpointInsideQuad(site, Rupture4P))
		{
			return 0;
		}
		else
		{
			float dist1 = Distance_PointAndLineSegment(site, Rupture4P[0], Rupture4P[1]);
			float dist2 = Distance_PointAndLineSegment(site, Rupture4P[1], Rupture4P[2]);
			float dist3 = Distance_PointAndLineSegment(site, Rupture4P[2], Rupture4P[3]);
			float dist4 = Distance_PointAndLineSegment(site, Rupture4P[3], Rupture4P[0]);
			return min(min(dist1, dist2), min(dist3, dist4));
		}
	}
	double Rx(double lon, double lat) const
		//到断裂面上边缘投影的垂线距离, Hanging wall side为正
	{
		double x = getX(lon);
		double y = getY(lat);
		double z = 0;
		Vector3f site(x, y, z);
		//断裂面矩形四个点坐标, 以断裂面中心为原点
		vector<Vector3f> Rupture4P = Rupture4Points();
		for (int i = 0; i < Rupture4P.size(); i++)
		{
			Rupture4P[i][2] = 0;
		}
		//垂线距离
		double result1 = (Rupture4P[0] - site).cross(Rupture4P[1] - site).norm() 
			/ (Rupture4P[0] - Rupture4P[1]).norm();
		assert(result1 >= 0);
		//判断Hanging wall side
		Vector3f n_xy = RuptureNormal; n_xy[2] = 0; n_xy.normalize();
		Vector3f center = 0.5 * (Rupture4P[0] + Rupture4P[1]);
		if ((site - center).dot(n_xy) > 0)
		{
			return result1;
		}
		else
		{
			return -result1;
		}
	}
};
