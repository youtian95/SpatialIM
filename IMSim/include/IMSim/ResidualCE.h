// 2021年3月24日16:24:52
// ***需要进一步测试 
// 
// 一次场景地震两个场地（距离为h）谱加速度的epsilon或者eta的相关系数rho
// rho_epsilon(h,T1,T2), rho_eta(T1,T2) （注意within-event与距离无关）, rho_total(h,T1,T2)

#pragma once
#define _USE_MATH_DEFINES 
#include <map>
#include <vector>
#include <cmath>
#include <random>
#include <assert.h>
#include "Eigen/Dense"
#include "Eigen/Geometry"
#include "GMPE.h"

using namespace std;
using namespace Eigen;

class ResidualCE
{
public:
	// 两个场地（距离为h）不同周期下的谱加速度的epsilon的相关系数
    static double rho_epsilon_loth_baker(double h, double T1, double T2)
        // Compute the spatial correlation of epsilons for the NGA ground motion models
        // 
        // The function is strictly empirical, fitted over the range the range 0.01s <= T1, T2 <= 10s
        // 
        // INPUT
        // 
        // T1, T2 = The two periods of interest.The periods may be equal,
        // and there is no restriction on which one is larger.
        // h = The separation distance between two sites(units of km)
        // 
        // OUTPUT
        // 
        // rho = The predicted correlation coefficient
        //
        // 参考文献：
        // Loth C, Baker JW. A spatial cross-correlation model of spectral accelerations 
        // at multiple periods. Earthquake Eng Struc. 2013;42:397-417. 
        // https://doi.org/10.1002/eqe.2212
    {
        // Verify the validity of input arguments
        if (min(T1, T2) < 0.01)
        {
            throw "The periods must be greater or equal to 0.01s";
        }
        if (max(T1, T2) > 10.0)
        { 
            throw "The periods must be less or equal to 10s";
        }
        if (h < 0.0)
        { 
            throw "The separation distance must be positive";
        }

        VectorXf Tlist(9);
        Tlist << 0.01, 0.1, 0.2, 0.5, 1, 2, 5, 7.5, 10.0001;

        ///////////////////////////////////////////////////////////////////////////////
        //
        // Model coefficients--note that these were updated in 2019 due to an error
        // in the original publication.The original publication's coefficients are
        // at the bottom of this function.
        //
        //Table II.Short range coregionalization matrix, B1
        Matrix<float, 9, 9> B1;
        B1 << 0.29, 0.25, 0.23, 0.23, 0.18, 0.1, 0.06, 0.06, 0.06,
            0.25, 0.30, 0.2, 0.16, 0.1, 0.04, 0.03, 0.04, 0.05,
            0.23, 0.20, 0.27, 0.18, 0.1, 0.03, 0, 0.01, 0.02,
            0.23, 0.16, 0.18, 0.31, 0.22, 0.14, 0.08, 0.07, 0.07,
            0.18, 0.10, 0.1, 0.22, 0.33, 0.24, 0.16, 0.13, 0.12,
            0.10, 0.04, 0.03, 0.14, 0.24, 0.33, 0.26, 0.21, 0.19,
            0.06, 0.03, 0, 0.08, 0.16, 0.26, 0.37, 0.3, 0.26,
            0.06, 0.04, 0.01, 0.07, 0.13, 0.21, 0.3, 0.28, 0.24,
            0.06, 0.05, 0.02, 0.07, 0.12, 0.19, 0.26, 0.24, 0.23;

        // Table III.Long range coregionalization matrix, B2
        Matrix<float, 9, 9> B2;
        B2 << 0.47, 0.4, 0.43, 0.35, 0.27, 0.15, 0.13, 0.09, 0.12,
            0.4, 0.42, 0.37, 0.25, 0.15, 0.03, 0.04, 0, 0.03,
            0.43, 0.37, 0.45, 0.36, 0.26, 0.15, 0.09, 0.05, 0.08,
            0.35, 0.25, 0.36, 0.42, 0.37, 0.29, 0.2, 0.16, 0.16,
            0.27, 0.15, 0.26, 0.37, 0.48, 0.41, 0.26, 0.21, 0.21,
            0.15, 0.03, 0.15, 0.29, 0.41, 0.55, 0.37, 0.33, 0.32,
            0.13, 0.04, 0.09, 0.2, 0.26, 0.37, 0.51, 0.49, 0.49,
            0.09, 0.00, 0.05, 0.16, 0.21, 0.33, 0.49, 0.62, 0.6,
            0.12, 0.03, 0.08, 0.16, 0.21, 0.32, 0.49, 0.6, 0.68;

        // Table IV.Nugget effect coregionalization matrix, B3
        Matrix<float, 9, 9> B3;
        B3 << 0.24, 0.22, 0.21, 0.09, -0.02, 0.01, 0.03, 0.02, 0.01,
            0.22, 0.28, 0.2, 0.04, -0.05, 0, 0.01, 0.01, -0.01,
            0.21, 0.20, 0.28, 0.05, -0.06, 0, 0.04, 0.03, 0.01,
            0.09, 0.04, 0.05, 0.26, 0.14, 0.05, 0.05, 0.05, 0.04,
            -0.02, -0.05, -0.06, 0.14, 0.20, 0.07, 0.05, 0.05, 0.05,
            0.01, 0.00, 0.00, 0.05, 0.07, 0.12, 0.08, 0.07, 0.06,
            0.03, 0.01, 0.04, 0.05, 0.05, 0.08, 0.12, 0.1, 0.08,
            0.02, 0.01, 0.03, 0.05, 0.05, 0.07, 0.1, 0.1, 0.09,
            0.01, -0.01, 0.01, 0.04, 0.05, 0.06, 0.08, 0.09, 0.09;

        // Find in which interval each input period is located
        int index1, index2;
        for (size_t i = 0; i < Tlist.size() - 1; i++)
        {
            if ((T1 - Tlist(i) >= 0.0) && (T1 - Tlist(i + 1) < 0.0))
                index1 = i;
            if ((T2 - Tlist(i) >= 0.0) && (T2 - Tlist(i + 1) < 0.0))
                index2 = i;
        }

        // Linearly interpolate the corresponding value of each coregionalization
        // matrix coefficient

        float B1coeff1 = B1(index1, index2) + (B1(index1 + 1, index2) - B1(index1, index2)) /
            (Tlist(index1 + 1) - Tlist(index1)) * (T1 - Tlist(index1));

        float B1coeff2 = B1(index1, index2 + 1) + (B1(index1 + 1, index2 + 1) - B1(index1, index2 + 1)) /
            (Tlist(index1 + 1) - Tlist(index1)) * (T1 - Tlist(index1));

        float B1coeff0 = B1coeff1 + (B1coeff2 - B1coeff1) /
            (Tlist(index2 + 1) - Tlist(index2)) * (T2 - Tlist(index2));


        float B2coeff1 = B2(index1, index2) + (B2(index1 + 1, index2) - B2(index1, index2)) /
            (Tlist(index1 + 1) - Tlist(index1)) * (T1 - Tlist(index1));

        float B2coeff2 = B2(index1, index2 + 1) + (B2(index1 + 1, index2 + 1) - B2(index1, index2 + 1)) /
            (Tlist(index1 + 1) - Tlist(index1)) * (T1 - Tlist(index1));

        float B2coeff0 = B2coeff1 + (B2coeff2 - B2coeff1) /
            (Tlist(index2 + 1) - Tlist(index2)) * (T2 - Tlist(index2));


        float B3coeff1 = B3(index1, index2) + (B3(index1 + 1, index2) - B3(index1, index2)) /
            (Tlist(index1 + 1) - Tlist(index1)) * (T1 - Tlist(index1));

        float B3coeff2 = B3(index1, index2 + 1) + (B3(index1 + 1, index2 + 1) - B3(index1, index2 + 1)) /
            (Tlist(index1 + 1) - Tlist(index1)) * (T1 - Tlist(index1));

        float B3coeff0 = B3coeff1 + (B3coeff2 - B3coeff1) /
            (Tlist(index2 + 1) - Tlist(index2)) * (T2 - Tlist(index2));


        // Compute the correlation coefficient(Equation 42)

        double rho = B1coeff0 * exp(-3.0 * h / 20.0) + B2coeff0 * exp(-3.0 * h / 70.0);

        if (h == 0.0)
        { 
            rho = B1coeff0 * exp(-3.0 * h / 20.0) + B2coeff0 * exp(-3.0 * h / 70.0) + B3coeff0;
        }

        return rho;

        /////////////////////////////////////////////////////////////////////////
        //
        // Below are the original tables from the 2013 manuscript, but they are
        // superceded now that we have found an error in the original calculations.
        // Uncomment theseand move them up in the function to reproduce the
        // original paper's results.
        //
        // %Table II.Short range coregionalization matrix, B1
        // B1 = [0.30 0.24 0.23 0.22 0.16 0.07 0.03 0 0;
        // 0.24 0.27 0.19 0.13 0.08 0 0 0 0;
        // 0.23 0.19 0.26 0.19 0.12 0.04 0 0 0;
        // 0.22 0.13 0.19 0.32 0.23 0.14 0.09 0.06 0.04;
        // 0.16 0.08 0.12 0.23 0.32 0.22 0.13 0.09 0.07;
        // 0.07 0 0.04 0.14 0.22 0.33 0.23 0.19 0.16;
        // 0.03 0 0 0.09 0.13 0.23 0.34 0.29 0.24;
        // 0 0 0 0.06 0.09 0.19 0.29 0.30 0.25;
        // 0 0 0 0.04 0.07 0.16 0.24 0.25 0.24];
        //
        // Table III.Long range coregionalization matrix, B2
        // B2 = [0.31 0.26 0.27 0.24 0.17 0.11 0.08 0.06 0.05;
        // 0.26 0.29 0.22 0.15 0.07 0 0 0 - 0.03;
        // 0.27 0.22 0.29 0.24 0.15 0.09 0.03 0.02 0;
        // 0.24 0.15 0.24 0.33 0.27 0.23 0.17 0.14 0.14;
        // 0.17 0.07 0.15 0.27 0.38 0.34 0.23 0.19 0.21;
        // 0.11 0 0.09 0.23 0.34 0.44 0.33 0.29 0.32;
        // 0.08 0 0.03 0.17 0.23 0.33 0.45 0.42 0.42;
        // 0.06 0 0.02 0.14 0.19 0.29 0.42 0.47 0.47;
        // 0.05 - 0.03 0 0.14 0.21 0.32 0.42 0.47 0.54];
        //
        // Table IV.Nugget effect coregionalization matrix, B3
        // B3 = [0.38 0.36 0.35 0.17 0.04 0.04 0 0.03 0.08;
        // 0.36 0.43 0.35 0.13 0 0.02 0 0.02 0.08;
        // 0.35 0.35 0.45 0.11 - 0.04 - 0.02 - 0.04 - 0.02 0.03;
        // 0.17 0.13 0.11 0.35 0.2 0.06 0.02 0.04 0.02;
        // 0.04 0 - 0.04 0.20 0.30 0.14 0.09 0.12 0.04;
        // 0.04 0.02 - 0.02 0.06 0.14 0.22 0.12 0.13 0.09;
        // 0 0 - 0.04 0.02 0.09 0.12 0.21 0.17 0.13;
        // 0.03 0.02 - 0.02 0.04 0.12 0.13 0.17 0.23 0.10;
        // 0.08 0.08 0.03 0.02 0.04 0.09 0.13 0.10 0.22];

	}

    // 单个场地不同周期下的谱加速度的total的相关系数（包含事件内和事件间的相关性）
    static double rho_total_baker_jayaram(double T1, double T2)
        // Compute the correlation of epsilons for the NGA ground motion models
        //
        // The function is strictly empirical, fitted over the range the range 0.01s <= T1, T2 <= 10s
        //
        // INPUT
        //
        // T1, T2 = The two periods of interest.The periods may be equal,
        // and there is no restriction on which one is larger.
        //
        // OUTPUT
        //
        // rho = The predicted correlation coefficient
        //
        // 参考文献：
        // Baker JW, Jayaram N. Correlation of spectral acceleration values 
        // from NGA ground motion models. Earthquake Spectra. 2008;24:299-317. 
        // https://doi.org/10.1193/1.2857544
    {
        double T_min = min(T1, T2);
        double T_max = max(T1, T2);
        double pi = 3.1415926;

        double C1, C2, C3, C4;
        
        C1 = (1.0 - std::cos(pi / 2.0 - std::log(T_max / max(T_min, 0.109)) * 0.366));

        if (T_max < 0.2)
        { 
            C2 = 1.0 - 0.105 * (1.0 - 1.0 / (1.0 + std::exp(100.0 * T_max - 5.0)))
                * (T_max - T_min) / (T_max - 0.0099);
        }

        if (T_max < 0.109)
        { 
            C3 = C2;
        }
        else
        { 
            C3 = C1;
        }

        C4 = C1 + 0.5 * (std::sqrt(C3) - C3) * (1.0 + std::cos(pi * (T_min) / (0.109)));

        double rho;

        if (T_max <= 0.109)
        {
            rho = C2;
        }
        else if (T_min > 0.109)
        { 
            rho = C1;
        }
        else if (T_max < 0.2)
        { 
            rho = min(C2, C4);
        }
        else
        { 
            rho = C4;
        }

        return rho;
    }

    // 不同周期下的谱加速度的etan的相关系数
    //
    // loth&baker(2013)的epsilon的相关系数 + baker&jayaram(2008)的total的相关系数 
    //  + Campbell&Bozorgnia(2014)的GMPE估计标准差 + Goda&Hong(2008)的方程
    // 来估计出etan的相关系数
    static double rho_eta_Combined_Method(double T1, double T2)
    {
        if (T1 <= 0.01) T1 = 0.01;
        if (T2 <= 0.01) T2 = 0.01;
        if (T1 == T2) return 1.0;

        VectorXf Sa;
        VectorXf sigma;
        VectorXf tau; 
        VectorXf period1;
        VectorXf T(1); 
        GMPE gmpe;

        T << T1;
        gmpe.CB_2014_nga(&Sa, &sigma, &tau, &period1,
            7, T, 20, 20, 20,
            10, 5, 15, 45, 0, 0,
            300, 999, 999, 0);
        double sigma_eta_T1 = tau(0);
        double sigma_epsilon_T1 = std::sqrt(sigma(0) * sigma(0) - tau(0) * tau(0));

        T(0) = T2;
        gmpe.CB_2014_nga(&Sa, &sigma, &tau, &period1,
            7, T, 20, 20, 20,
            10, 5, 15, 45, 0, 0,
            300, 999, 999, 0);
        double sigma_eta_T2 = tau(0);
        double sigma_epsilon_T2 = std::sqrt(sigma(0) * sigma(0) - tau(0) * tau(0));

        double rho_total_h_T1_T2 = rho_total_baker_jayaram(T1, T2); //h=0
        double rho_epsilon_h_T1_T2 = rho_epsilon_loth_baker(0, T1, T2);

        return rho_eta_T1_T2(T1, T2, rho_total_h_T1_T2, rho_epsilon_h_T1_T2,
            sigma_epsilon_T1, sigma_epsilon_T2,
            sigma_eta_T1, sigma_eta_T2);
    }
    static MatrixXf rho_eta_Combined_Method(VectorXf T)
        //返回相关系数矩阵
    {
        MatrixXf CE(T.size(), T.size());
        for (size_t i = 0; i < T.size(); i++)
        {
            for (size_t j = 0; j < T.size(); j++)
            {
                if (i == j)
                {
                    CE(i, j) = 1;
                }
                else if (i > j)
                    //左下
                {
                    CE(i, j) = CE(j, i);
                }
                else
                    //右上
                {
                    CE(i, j) = rho_eta_Combined_Method(T(i), T(j));
                }
            }
        }
        return CE;
    }

    //模拟生成 between-event residuals 
    static MatrixXf Inter_Event_Residuals_Simulation(VectorXf T, int N, default_random_engine* pRND)
        // 输入：
        // VectorXf T - 周期向量
        // int N - 模拟的次数
        // default_random_engine* pRND - 生成随机数的引擎指针
        //
        // 输出：
        // MatrixXf - (i_sim, i_T), 互相相关的标准正态分布变量
    {
        MatrixXf rho_eta = rho_eta_Combined_Method(T);

        MatrixXf L;
        //乔斯基分解
        Eigen::MatrixXf normTransform(T.size(), T.size());
        Eigen::LLT<Eigen::MatrixXf> cholSolver(rho_eta);
        // We can only use the cholesky decomposition if 
        // the covariance matrix is symmetric, pos-definite.
        // But a covariance matrix might be pos-semi-definite.
        // In that case, we'll go to an EigenSolver
        if (cholSolver.info() == Eigen::Success) {
            // Use cholesky solver
            normTransform = cholSolver.matrixL();
        }
        else {
            //cout << "警告：相关系数矩阵的cholesky分解失败，采用特征值分解，且将负特征值置为0！" << endl;
            // Use eigen solver
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigenSolver(rho_eta);
            VectorXf eigenvalues = eigenSolver.eigenvalues();
            for (size_t i = 0; i < eigenvalues.size(); i++)
            {
                if (eigenvalues(i) < 0) eigenvalues(i) = 0;
            }
            normTransform = eigenSolver.eigenvectors()
                * eigenvalues.cwiseSqrt().asDiagonal();
        }
        L = normTransform;

        // 生成相关的变量
        MatrixXf residuals(N,T.size());
        std::normal_distribution<> dist;
        for (int i_sim = 0; i_sim < N; i_sim++)
        {
            VectorXf randN(T.size());
            for (int j = 0; j < T.size(); j++)
            {
                randN(j) = dist(*pRND);
            }
            residuals.row(i_sim) = L * randN;
        }

        return residuals;
    }

    // rho_epsilon(h,T1,T2), rho_eta(T1,T2), rho_total(h,T1,T2)之间的关系
    // 
    // rho_total(h,T1,T2) * sigma_total(T1) * sigma_total(T2) 
    //      = rho_eta(T1,T2) * sigma_eta(T1) * sigma_eta(T2)
    //        + rho_epsilon(h,T1,T2) * sigma_epsilon(T1) * sigma_epsilon(T2)
    //
    // 参考文献：
    // Goda K, Hong HP. Spatial correlation of peak ground motions and 
    // response spectra. B Seismol Soc Am. 2008;98:354-65. 
    // https://doi.org/10.1785/0120070078
    static double rho_total_h_T1_T2(double h, double T1, double T2,
        double rho_epsilon_h_T1_T2, double rho_eta_T1_T2,
        double sigma_epsilon_T1, double sigma_epsilon_T2, 
        double sigma_eta_T1, double sigma_eta_T2)
    {
        double sigma_total_T1 = std::sqrt(sigma_eta_T1 * sigma_eta_T1 + sigma_epsilon_T1 * sigma_epsilon_T1);
        double sigma_total_T2 = std::sqrt(sigma_eta_T2 * sigma_eta_T2 + sigma_epsilon_T2 * sigma_epsilon_T2);

        double rho_total = (rho_epsilon_h_T1_T2 * sigma_epsilon_T1 * sigma_epsilon_T2
            + rho_eta_T1_T2 * sigma_eta_T1 * sigma_eta_T2)
            / (sigma_total_T1 * sigma_total_T2);

        return rho_total;
    }
    static double rho_eta_T1_T2(double T1, double T2,
        double rho_total_h_T1_T2, double rho_epsilon_h_T1_T2,
        double sigma_epsilon_T1, double sigma_epsilon_T2,
        double sigma_eta_T1, double sigma_eta_T2)
        // rho_total_h_T1_T2和rho_epsilon_h_T1_T2为相同h下的相关系数
    {
        double sigma_total_T1 = std::sqrt(sigma_eta_T1 * sigma_eta_T1 + sigma_epsilon_T1 * sigma_epsilon_T1);
        double sigma_total_T2 = std::sqrt(sigma_eta_T2 * sigma_eta_T2 + sigma_epsilon_T2 * sigma_epsilon_T2);

        double rho_eta_T1_T2 = (rho_total_h_T1_T2 * sigma_total_T1 * sigma_total_T2
            - rho_epsilon_h_T1_T2 * sigma_epsilon_T1 * sigma_epsilon_T2)
            / (sigma_eta_T1 * sigma_eta_T2);

        return rho_eta_T1_T2;
    }

};

