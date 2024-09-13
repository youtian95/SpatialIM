//Ground motion prediction equations

#pragma once
#include <cmath>
#include <map>
#include <vector>
#include <assert.h>
#include "Eigen/Dense"
#include "Eigen/Geometry"

using namespace std;
using namespace Eigen;

#ifndef M_PI
#define M_PI			3.14159265358979323846264338327950288419716939937510582   // pi
#endif

class GMPE
{
public:
	GMPE()
	{
		c0 << -4.365, -4.348, -4.024, -3.479, -3.293, -3.666, -4.866, -5.411, -5.962, -6.403, -7.566, -8.379, -9.841, -11.011, -12.469, -12.969, -13.306, -14.02, -14.558, -15.509, -15.975, -4.416, -2.895;
		c1 << 0.977, 0.976, 0.931, 0.887, 0.902, 0.993, 1.267, 1.366, 1.458, 1.528, 1.739, 1.872, 2.021, 2.180, 2.270, 2.271, 2.150, 2.132, 2.116, 2.223, 2.132, 0.984, 1.510;
		c2 << 0.533, 0.549, 0.628, 0.674, 0.726, 0.698, 0.510, 0.447, 0.274, 0.193, -0.020, -0.121, -0.042, -0.069, 0.047, 0.149, 0.368, 0.726, 1.027, 0.169, 0.367, 0.537, 0.270;
		c3 << -1.485, -1.488, -1.494, -1.388, -1.469, -1.572, -1.669, -1.750, -1.711, -1.770, -1.594, -1.577, -1.757, -1.707, -1.621, -1.512, -1.315, -1.506, -1.721, -0.756, -0.800, -1.499, -1.299;
		c4 << -0.499, -0.501, -0.517, -0.615, -0.596, -0.536, -0.490, -0.451, -0.404, -0.321, -0.426, -0.440, -0.443, -0.527, -0.630, -0.768, -0.890, -0.885, -0.878, -1.077, -1.282, -0.496, -0.453;
		c5 << -2.773, -2.772, -2.782, -2.791, -2.745, -2.633, -2.458, -2.421, -2.392, -2.376, -2.303, -2.296, -2.232, -2.158, -2.063, -2.104, -2.051, -1.986, -2.021, -2.179, -2.244, -2.773, -2.466;
		// c6 = 0.248,	0.247,	0.246,	0.240,	0.227,	0.210,	0.183,	0.182,	0.189,	0.195,	0.185,	0.186,	0.186,	0.169,	0.158,	0.158,	0.148,	0.135,	0.140,	0.178,	0.194,	0.248,	0.204; // these coefficients are in the PEER report
		c6 << 0.248, 0.247, 0.246, 0.240, 0.227, 0.210, 0.183, 0.182, 0.189, 0.195, 0.185, 0.186, 0.186, 0.169, 0.158, 0.158, 0.148, 0.135, 0.135, 0.165, 0.180, 0.248, 0.204; // these coefficients are in the Earthquake Spectra paper
		c7 << 6.753, 6.502, 6.291, 6.317, 6.861, 7.294, 8.031, 8.385, 7.534, 6.990, 7.012, 6.902, 5.522, 5.650, 5.795, 6.632, 6.759, 7.978, 8.538, 8.468, 6.564, 6.768, 5.837;
		c8 << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		c9 << -0.214, -0.208, -0.213, -0.244, -0.266, -0.229, -0.211, -0.163, -0.150, -0.131, -0.159, -0.153, -0.090, -0.105, -0.058, -0.028, 0, 0, 0, 0, 0, -0.212, -0.168;
		c10 << 0.720, 0.730, 0.759, 0.826, 0.815, 0.831, 0.749, 0.764, 0.716, 0.737, 0.738, 0.718, 0.795, 0.556, 0.480, 0.401, 0.206, 0.105, 0, 0, 0, 0.720, 0.305;
		c11 << 1.094, 1.149, 1.290, 1.449, 1.535, 1.615, 1.877, 2.069, 2.205, 2.306, 2.398, 2.355, 1.995, 1.447, 0.330, -0.514, -0.848, -0.793, -0.748, -0.664, -0.576, 1.090, 1.713;
		c12 << 2.191, 2.189, 2.164, 2.138, 2.446, 2.969, 3.544, 3.707, 3.343, 3.334, 3.544, 3.016, 2.616, 2.470, 2.108, 1.327, 0.601, 0.568, 0.356, 0.075, -0.027, 2.186, 2.602;
		c13 << 1.416, 1.453, 1.476, 1.549, 1.772, 1.916, 2.161, 2.465, 2.766, 3.011, 3.203, 3.333, 3.054, 2.562, 1.453, 0.657, 0.367, 0.306, 0.268, 0.374, 0.297, 1.420, 2.457;
		c14 << -0.0070, -0.0167, -0.0422, -0.0663, -0.0794, -0.0294, 0.0642, 0.0968, 0.1441, 0.1597, 0.1410, 0.1474, 0.1764, 0.2593, 0.2881, 0.3112, 0.3478, 0.3747, 0.3382, 0.3754, 0.3506, -0.0064, 0.1060;
		c15 << -0.207, -0.199, -0.202, -0.339, -0.404, -0.416, -0.407, -0.311, -0.172, -0.084, 0.085, 0.233, 0.411, 0.479, 0.566, 0.562, 0.534, 0.522, 0.477, 0.321, 0.174, -0.202, 0.332;
		c16 << 0.390, 0.387, 0.378, 0.295, 0.322, 0.384, 0.417, 0.404, 0.466, 0.528, 0.540, 0.638, 0.776, 0.771, 0.748, 0.763, 0.686, 0.691, 0.670, 0.757, 0.621, 0.393, 0.585;
		c17 << 0.0981, 0.1009, 0.1095, 0.1226, 0.1165, 0.0998, 0.0760, 0.0571, 0.0437, 0.0323, 0.0209, 0.0092, -0.0082, -0.0131, -0.0187, -0.0258, -0.0311, -0.0413, -0.0281, -0.0205, 0.0009, 0.0977, 0.0517;
		c18 << 0.0334, 0.0327, 0.0331, 0.0270, 0.0288, 0.0325, 0.0388, 0.0437, 0.0463, 0.0508, 0.0432, 0.0405, 0.0420, 0.0426, 0.0380, 0.0252, 0.0236, 0.0102, 0.0034, 0.0050, 0.0099, 0.0333, 0.0327;
		c19 << 0.00755, 0.00759, 0.00790, 0.00803, 0.00811, 0.00744, 0.00716, 0.00688, 0.00556, 0.00458, 0.00401, 0.00388, 0.00420, 0.00409, 0.00424, 0.00448, 0.00345, 0.00603, 0.00805, 0.00280, 0.00458, 0.00757, 0.00613;
		c20 << -0.0055, -0.0055, -0.0057, -0.0063, -0.0070, -0.0073, -0.0069, -0.0060, -0.0055, -0.0049, -0.0037, -0.0027, -0.0016, -0.0006, 0, 0, 0, 0, 0, 0, 0, -0.0055, -0.0017;
		
		Dc20_ << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		Dc20_JI << -0.0035, -0.0035, -0.0034, -0.0037, -0.0037, -0.0034, -0.003, -0.0031, -0.0033, -0.0035, -0.0034, -0.0034, -0.0032, -0.003, -0.0019, -0.0005, 0, 0, 0, 0, 0, -0.0035, -0.0006;
		Dc20_CH << 0.0036, 0.0036, 0.0037, 0.004, 0.0039, 0.0042, 0.0042, 0.0041, 0.0036, 0.0031, 0.0028, 0.0025, 0.0016, 0.0006, 0, 0, 0, 0, 0, 0, 0, 0.0036, 0.0017;
		a2 << 0.168, 0.166, 0.167, 0.173, 0.198, 0.174, 0.198, 0.204, 0.185, 0.164, 0.16, 0.184, 0.216, 0.596, 0.596, 0.596, 0.596, 0.596, 0.596, 0.596, 0.596, 0.167, 0.596;
		h1 << 0.242, 0.244, 0.246, 0.251, 0.26, 0.259, 0.254, 0.237, 0.206, 0.21, 0.226, 0.217, 0.154, 0.117, 0.117, 0.117, 0.117, 0.117, 0.117, 0.117, 0.117, 0.241, 0.117;
		h2 << 1.471, 1.467, 1.467, 1.449, 1.435, 1.449, 1.461, 1.484, 1.581, 1.586, 1.544, 1.554, 1.626, 1.616, 1.616, 1.616, 1.616, 1.616, 1.616, 1.616, 1.616, 1.474, 1.616;
		h3 << -0.714, -0.711, -0.713, -0.701, -0.695, -0.708, -0.715, -0.721, -0.787, -0.795, -0.77, -0.77, -0.78, -0.733, -0.733, -0.733, -0.733, -0.733, -0.733, -0.733, -0.733, -0.715, -0.733;
		h4 << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
		h5 << -0.336, -0.339, -0.338, -0.338, -0.347, -0.391, -0.449, -0.393, -0.339, -0.447, -0.525, -0.407, -0.371, -0.128, -0.128, -0.128, -0.128, -0.128, -0.128, -0.128, -0.128, -0.337, -0.128;
		h6 << -0.27, -0.263, -0.259, -0.263, -0.219, -0.201, -0.099, -0.198, -0.21, -0.121, -0.086, -0.281, -0.285, -0.756, -0.756, -0.756, -0.756, -0.756, -0.756, -0.756, -0.756, -0.27, -0.756;
		k1 << 865, 865, 908, 1054, 1086, 1032, 878, 748, 654, 587, 503, 457, 410, 400, 400, 400, 400, 400, 400, 400, 400, 865, 400;
		k2 << -1.186, -1.219, -1.273, -1.346, -1.471, -1.624, -1.931, -2.188, -2.381, -2.518, -2.657, -2.669, -2.401, -1.955, -1.025, -0.299, 0, 0, 0, 0, 0, -1.186, -1.955;
		k3 << 1.839, 1.84, 1.841, 1.843, 1.845, 1.847, 1.852, 1.856, 1.861, 1.865, 1.874, 1.883, 1.906, 1.929, 1.974, 2.019, 2.11, 2.2, 2.291, 2.517, 2.744, 1.839, 1.929;

		f1 << 0.734, 0.738, 0.747, 0.777, 0.782, 0.769, 0.769, 0.761, 0.744, 0.727, 0.69, 0.663, 0.606, 0.579, 0.541, 0.529, 0.527, 0.521, 0.502, 0.457, 0.441, 0.734, 0.655;
		f2 << 0.492, 0.496, 0.503, 0.52, 0.535, 0.543, 0.543, 0.552, 0.545, 0.568, 0.593, 0.611, 0.633, 0.628, 0.603, 0.588, 0.578, 0.559, 0.551, 0.546, 0.543, 0.492, 0.494;
		t1 << 0.404, 0.417, 0.446, 0.508, 0.504, 0.445, 0.382, 0.339, 0.34, 0.34, 0.356, 0.379, 0.43, 0.47, 0.497, 0.499, 0.5, 0.543, 0.534, 0.523, 0.466, 0.409, 0.317;
		t2 << 0.325, 0.326, 0.344, 0.377, 0.418, 0.426, 0.387, 0.338, 0.316, 0.3, 0.264, 0.263, 0.326, 0.353, 0.399, 0.4, 0.417, 0.393, 0.421, 0.438, 0.438, 0.322, 0.297;
		flnAF << 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3;
		rlnPGA_lnY << 1, 0.998, 0.986, 0.938, 0.887, 0.87, 0.876, 0.87, 0.85, 0.819, 0.743, 0.684, 0.562, 0.467, 0.364, 0.298, 0.234, 0.202, 0.184, 0.176, 0.154, 1, 0.684;
	}
	void CB_2014_nga(VectorXf* Sa, VectorXf* sigma, VectorXf* tau, VectorXf* period1,
		double M, VectorXf T, double Rrup, double Rjb, double Rx,
		double W, double Ztor, double Zbot, double delta, double lambda, bool Fhw,
		double Vs30, double Z25, double Zhyp, int region) const
		// Campbell and Bozorgnia 2014 ground motion prediciton model. Citation for the model :
		//		Campbell, K.W., and Bozorgnia, Y. (2014). "NGA-West2 Ground Motion Model 
		//		for the Average Horizontal Components of PGA, PGV, and 5 % Damped Linear
		//		Acceleration Response Spectra." Earthquake Spectra, 30(3), 1087-1115.
		//
		//	Provides ground - mtion prediction equations for computing mediansand
		//		standard deviations of average horizontal components of PGA, PGV and 5 %
		//		damped linear pseudo - absolute aceeleration response spectra
		// ------------------------------------------------------------------------------------
		// Input Variables
		// M			= Magnitude
		// T			= Period(sec);
		//					Use 1000 for output the array of Sa with period
		// Rrup			= Closest distance coseismic rupture(km)
		// Rjb			= Joyner - Boore distance(km)
		// Rx			= Closest distance to the surface projection of the
		//					coseismic fault rupture plane
		// W			= down - dip width of the fault rupture plane
		//					if unknown, input: 999
		// Ztor			= Depth to the top of coseismic rupture(km)
		//					if unknown, input : 999
		// Zbot			= Depth to the bottom of the seismogenic crust
		//					needed only when W is unknow;
		// delta		= average dip of the rupture place(degree)
		// lambda		= rake angle(degree) - average angle of slip measured in
		//					the plance of rupture
		// Fhw			= hanging wall effect
		//				= 1 for including
		//				= 0 for excluding
		// Vs30			= shear wave velocity averaged over top 30 m(m / s)
		// Z25			= Depth to the 2.5 km / s shear - wave velocity horizon(km)
		//					if in California or Japan and Z2.5 is unknow, then
		//					input: 999
		// Zhyp			= Hypocentral depth of the earthquake measured from sea level
		//					if unknown, input : 999
		// region		= 0 for global(incl.Taiwan)
		//				= 1 for California
		//				= 2 for Japan
		//				= 3 for China or Turkey
		//				= 4 for Italy
		//
		// Output Variables
		// Sa			= Median spectral acceleration prediction
		// sigma		= logarithmic standard deviation of spectral acceleration prediction
		// tau			= inter-events logarithmic standard deviation
		// period1		= period array, if input T = 1000. Outputs above will be arrays as well
		// ------------------------------------------------------------------------------------
	{
		// Period
		VectorXf period(23);
		period << 0.010, 0.020, 0.030, 0.050, 0.075, 0.10, 0.15, 0.20, 0.25, 0.30,
			0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 0, -1;

		// Set initial A1100 value to 999. 
		// A1100: median estimated value of PGA on rock with Vs30 = 1100m / s
		double A1100 = 999;

		// Style of faulting
		bool Frv = (lambda > 30 & lambda < 150);
		bool Fnm = (lambda > -150 & lambda < -30);

		// if Ztor is unknown..
		if (Ztor == 999)
		{
			if (Frv == 1)
				Ztor = pow(max(2.704 - 1.226 * max(M - 5.849, 0.0), 0.0), 2);
			else
				Ztor = pow(max(2.673 - 1.136 * max(M - 4.970, 0.0), 0.0), 2);
		}

		// if W is unknown...
		if (W == 999)
		{
			double Ztori;
			if (Frv)
				Ztori = pow(max(2.704 - 1.226 * max(M - 5.849, 0.0), 0.0), 2);
			else
				Ztori = pow(max(2.673 - 1.136 * max(M - 4.970, 0.0), 0.0), 2);

			W = min(sqrt(pow(10.0, (M - 4.07) / 0.98)), (Zbot - Ztori) / sin(M_PI / 180.0 * delta));

			Zhyp = 9;
		}

		// if Zhyp is unknown...
		if (Zhyp == 999 && W != 999)
		{
			double fdZM;
			if (M < 6.75)
				fdZM = -4.317 + 0.984 * M;
			else
				fdZM = 2.325;

			double fdZD;
			if (delta <= 40)
				fdZD = 0.0445 * (delta - 40);
			else
				fdZD = 0;

			double Ztori;
			if (Frv == 1)
				Ztori = pow(max(2.704 - 1.226 * max(M - 5.849, 0.0), 0.0), 2);
			else
				Ztori = pow(max(2.673 - 1.136 * max(M - 4.970, 0.0), 0.0), 2);

			double Zbor = Ztori + W * sin(M_PI / 180.0 * delta); // The depth to the bottom of the rupture plane
			double d_Z = exp(min(fdZM + fdZD, log(0.9 * (Zbor - Ztori))));
			Zhyp = d_Z + Ztori;

		}


		if (T.size() == 1 && T[0] == 1000)
			// Compute Sa and sigma with pre - defined period
		{
			*Sa = VectorXf::Zero(period.size() - 2);
			*sigma = VectorXf::Zero(period.size() - 2);
			*tau = VectorXf::Zero(period.size() - 2);
			*period1 = period.head(period.size() - 2);

			for (int ipT = 1; ipT <= period1->size(); ipT++)
			{
				double Sa_temp, sigma_temp, tau_temp, PGA_temp;
				CB_2014_nga_sub(&Sa_temp, &sigma_temp, &tau_temp,
					M, ipT, Rrup, Rjb, Rx, W, Ztor, Zbot, delta, Fhw, Vs30,
					Z25, Zhyp, lambda, Frv, Fnm, region, A1100);
				(*Sa)[ipT - 1] = Sa_temp;
				(*sigma)[ipT - 1] = sigma_temp;
				(*tau)[ipT - 1] = tau_temp;

				CB_2014_nga_sub(&PGA_temp, &sigma_temp, &tau_temp,
					M, 22, Rrup, Rjb, Rx, W, Ztor, Zbot, delta,
					Fhw, Vs30, Z25, Zhyp, lambda, Frv, Fnm, region, A1100);
				double PGA = PGA_temp;
				if ((*Sa)[ipT - 1] < PGA && (*period1)[ipT - 1] < 0.25)
					(*Sa)[ipT - 1] = PGA;
			}
		}
		else
			//Compute Sa and sigma with user-defined period 
		{
			*Sa = VectorXf::Zero(T.size());
			*sigma = VectorXf::Zero(T.size());
			*tau = VectorXf::Zero(T.size());
			*period1 = T;

			for (int i = 1; i <= T.size(); i++)
			{
				double Ti = T[i - 1];

				int ip_T = 0;
				bool period_has_Ti = 0;
				for (int i_period = 0; i_period < period.size(); i_period++)
				{
					if (abs(period[i_period] - Ti) < 0.0001)
					{
						ip_T = i_period + 1;
						period_has_Ti = 1;
						break;
					}
				}

				if (!period_has_Ti)
					// The user defined period requires interpolation
				{
					int ip_low = 1;
					double T_low = 0;
					{
						for (int i_period = 0; i_period < period.size(); i_period++)
						{
							if (period[i_period] < Ti)
							{
								if ((Ti - period[i_period]) < (Ti - T_low))
								{
									ip_low = i_period + 1;
									T_low = period[i_period];
								}
							}
						}
					}
					int ip_high = period.size();
					double T_high = period.maxCoeff();
					{
						for (int i_period = 0; i_period < period.size(); i_period++)
						{
							if (period[i_period] > Ti)
							{
								if ((period[i_period] - Ti) < (T_high - Ti))
								{
									ip_high = i_period + 1;
									T_high = period[i_period];
								}
							}
						}
					}

					double Sa_low, sigma_low, tau_low;
					CB_2014_nga_sub(&Sa_low, &sigma_low, &tau_low,
						M, ip_low, Rrup, Rjb, Rx, W, Ztor, Zbot,
						delta, Fhw, Vs30, Z25, Zhyp, lambda, Frv, Fnm, region, A1100);
					double Sa_high, sigma_high, tau_high;
					CB_2014_nga_sub(&Sa_high, &sigma_high, &tau_high,
						M, ip_high, Rrup, Rjb, Rx, W, Ztor, Zbot,
						delta, Fhw, Vs30, Z25, Zhyp, lambda, Frv, Fnm, region, A1100);
					double PGA, sigma_temp, tau_temp;
					CB_2014_nga_sub(&PGA, &sigma_temp, &tau_temp,
						M, 22, Rrup, Rjb, Rx, W, Ztor, Zbot,
						delta, Fhw, Vs30, Z25, Zhyp, lambda, Frv, Fnm, region, A1100);

					VectorXf x(2), Y_sa(2), Y_sigma(2), Y_tau(2);
					x << log(T_low), log(T_high);
					Y_sa << log(Sa_low), log(Sa_high);
					Y_sigma << sigma_low, sigma_high;
					Y_tau << tau_low, tau_high;
					(*Sa)[i - 1] = exp(interp1(x, Y_sa, log(Ti)));
					(*sigma)[i - 1] = interp1(x, Y_sigma, log(Ti));
					(*tau)[i - 1] = interp1(x, Y_tau, log(Ti));

					if ((*Sa)[i - 1] < PGA && (*period1)[i - 1] < 0.25)
						(*Sa)[i - 1] = PGA;
				}
				else
				{
					double Sa_temp, sigma_temp, tau_temp, PGA;
					CB_2014_nga_sub(&Sa_temp, &sigma_temp, &tau_temp,
						M, ip_T, Rrup, Rjb, Rx, W, Ztor, Zbot, delta,
						Fhw, Vs30, Z25, Zhyp, lambda, Frv, Fnm, region, A1100);
					(*Sa)[i - 1] = Sa_temp;
					(*sigma)[i - 1] = sigma_temp;
					(*tau)[i - 1] = tau_temp;

					CB_2014_nga_sub(&PGA, &sigma_temp, &tau_temp,
						M, 22, Rrup, Rjb, Rx, W, Ztor, Zbot,
						delta, Fhw, Vs30, Z25, Zhyp, lambda, Frv, Fnm, region, A1100);

					if ((*Sa)[i - 1] < PGA && (*period1)[i - 1] < 0.25)
						(*Sa)[i - 1] = PGA;
				}
			}
		}
	}

	static double interp1(VectorXf x, VectorXf y, double vx)
		//²åÖµ, xµ¥µ÷
	{
		assert(x.size() == y.size());
		if (vx < x[0])
		{
			return y[0] + (y[1] - y[0]) / (x[1] - x[0]) * (vx - x[0]);
		}
		else if (vx > x.tail(1)[0])
		{
			auto y_ = y.tail(2);
			auto x_ = x.tail(2);
			return y_[1] + (y_[1] - y_[0]) / (x_[1] - x_[0]) * (vx - x_[1]);
		}
		else
		{
			int i = 0;
			for (; i < x.size() - 1; i++)
			{
				if (x[i] <= vx && x[i + 1] >= vx) break;
			}
			return y[i] + (y[1 + i] - y[i]) / (x[1 + i] - x[i]) * (vx - x[i]);
		}
	}

private:
	VectorXf c0 = VectorXf(23);
	VectorXf c1 = VectorXf(23);
	VectorXf c2 = VectorXf(23);
	VectorXf c3 = VectorXf(23);
	VectorXf c4 = VectorXf(23);
	VectorXf c5 = VectorXf(23);
	VectorXf c6 = VectorXf(23);
	VectorXf c7 = VectorXf(23);
	VectorXf c8 = VectorXf(23);
	VectorXf c9 = VectorXf(23);
	VectorXf c10 = VectorXf(23);
	VectorXf c11 = VectorXf(23);
	VectorXf c12 = VectorXf(23);
	VectorXf c13 = VectorXf(23);
	VectorXf c14 = VectorXf(23);
	VectorXf c15 = VectorXf(23);
	VectorXf c16 = VectorXf(23);
	VectorXf c17 = VectorXf(23);
	VectorXf c18 = VectorXf(23);
	VectorXf c19 = VectorXf(23);
	VectorXf c20 = VectorXf(23);
	VectorXf Dc20_ = VectorXf(23);
	VectorXf Dc20_JI = VectorXf(23);
	VectorXf Dc20_CH = VectorXf(23);
	VectorXf a2 = VectorXf(23);
	VectorXf h1 = VectorXf(23);
	VectorXf h2 = VectorXf(23);
	VectorXf h3 = VectorXf(23);
	VectorXf h4 = VectorXf(23);
	VectorXf h5 = VectorXf(23);
	VectorXf h6 = VectorXf(23);
	VectorXf k1 = VectorXf(23);
	VectorXf k2 = VectorXf(23);
	VectorXf k3 = VectorXf(23);
	VectorXf f1= VectorXf(23);
	VectorXf f2= VectorXf(23);
	VectorXf t1= VectorXf(23);
	VectorXf t2= VectorXf(23);
	VectorXf flnAF = VectorXf(23);
	VectorXf rlnPGA_lnY = VectorXf(23);

	void CB_2014_nga_sub(double* Sa, double* sigma, double* tau,
		double M, int ip, double Rrup, double Rjb, double Rx,
		double W, double Ztor, double Zbot, double delta, bool Fhw,
		double Vs30, double Z25in, double Zhyp, double lambda,
		bool Frv, bool Fnm, int region, double A1100) const 
	{
		double c = 1.88;
		double n = 1.18;

		auto Dc20 = Dc20_;
		// Adjustment factor based on region
		if (region == 2)
		{
			Dc20 = Dc20_JI;
		}
		else if (region == 4)
		{
			Dc20 = Dc20_JI;
		}
		else if (region == 3)
		{
			Dc20 = Dc20_CH;
		}
		// if region is in Japan...
		bool Sj = (region == 2);

		double Z25, Z25A;
		// if Z2.5 is unknown...
		if (Z25in == 999)
		{
			if (region != 2) // if in California or other locations
			{
				Z25 = exp(7.089 - 1.144 * log(Vs30));
				Z25A = exp(7.089 - 1.144 * log(1100.0));
			}
			else if (region == 2) // if in Japan
			{
				Z25 = exp(5.359 - 1.102 * log(Vs30));
				Z25A = exp(5.359 - 1.102 * log(1100.0));
			}
		}
		else
		{
			// Assign Z2.5 from user input into Z25 and calc Z25A for Vs30 = 1100m / s
			if (region != 2) // if in California or other locations
			{
				Z25 = Z25in;
				Z25A = exp(7.089 - 1.144 * log(1100.0));
			}
			else if (region == 2) // if in Japan
			{
				Z25 = Z25in;
				Z25A = exp(5.359 - 1.102 * log(1100.0));
			}
		}

		double fmag;
		// Magnitude dependence
		if (M <= 4.5)
		{
			fmag = c0[ip - 1] + c1[ip - 1] * M;
		}
		else if (M <= 5.5)
		{
			fmag = c0[ip - 1] + c1[ip - 1] * M + c2[ip - 1] * (M - 4.5);
		}
		else if (M <= 6.5)
		{
			fmag = c0[ip - 1] + c1[ip - 1] * M + c2[ip - 1] * (M - 4.5) + c3[ip - 1] * (M - 5.5);
		}
		else
		{
			fmag = c0[ip - 1] + c1[ip - 1] * M + c2[ip - 1] * (M - 4.5) + c3[ip - 1] * (M - 5.5) + c4[ip - 1] * (M - 6.5);
		}

		// Geometric attenuation term
		double fdis = (c5[ip - 1] + c6[ip - 1] * M) * log(sqrt(Rrup * Rrup + c7[ip - 1] * c7[ip - 1]));

		// Style of faulting
		double F_fltm;
		if (M <= 4.5)
		{
			F_fltm = 0;
		}
		else if (M <= 5.5)
		{
			F_fltm = M - 4.5;
		}
		else
		{
			F_fltm = 1;
		}
		double fflt = ((c8[ip - 1] * Frv) + (c9[ip - 1] * Fnm)) * F_fltm;

		// Hanging-wall effects

		double R1 = W * cos(M_PI / 180.0 * delta);// W - downdip width
		double R2 = 62.0 * M - 350.0;

		double f1_Rx = h1[ip - 1] + h2[ip - 1] * (Rx / R1) + h3[ip - 1] * (Rx / R1) * (Rx / R1);
		double f2_Rx = h4[ip - 1] + h5[ip - 1] * ((Rx - R1) / (R2 - R1)) + h6[ip - 1] * pow((Rx - R1) / (R2 - R1), 2);

		double f_hngRx;
		if (Fhw == 0)
		{
			f_hngRx = 0;
		}
		else if (Rx < R1 && Fhw == 1)
		{
			f_hngRx = f1_Rx;
		}
		else if (Rx >= R1 && Fhw == 1)
		{
			f_hngRx = max(f2_Rx, 0.0);
		}

		double f_hngRup;
		if (Rrup == 0)
		{
			f_hngRup = 1;
		}
		else
		{
			f_hngRup = (Rrup - Rjb) / Rrup;
		}

		double f_hngM;
		if (M <= 5.5)
		{
			f_hngM = 0;
		}
		else if (M <= 6.5)
		{
			f_hngM = (M - 5.5) * (1.0 + a2[ip - 1] * (M - 6.5));
		}
		else
		{
			f_hngM = 1.0 + a2[ip - 1] * (M - 6.5);
		}

		double f_hngZ;
		if (Ztor <= 16.66)
		{
			f_hngZ = 1.0 - 0.06 * Ztor;
		}
		else
		{
			f_hngZ = 0;
		}

		double f_hngdelta = (90.0 - delta) / 45.0;

		double fhng = c10[ip - 1] * f_hngRx * f_hngRup * f_hngM * f_hngZ * f_hngdelta;

		// Site conditions

		double f_siteG;
		if (Vs30 <= k1[ip - 1])
		{
			if (A1100 == 999)
			{
				VectorXf Sa_temp, sigma_temp, tau_temp, period1_temp;
				VectorXf T_temp(1); T_temp << 0;
				CB_2014_nga(&Sa_temp, &sigma_temp, &tau_temp, &period1_temp,
					M, T_temp, Rrup, Rjb, Rx, W, Ztor, Zbot, delta, lambda, Fhw, 1100, Z25A, Zhyp, region);
				assert(Sa_temp.size() == 1);
				A1100 = Sa_temp[0];
			}
			f_siteG = c11[ip - 1] * log(Vs30 / k1[ip - 1]) + k2[ip - 1] * (log(A1100 + c * pow(Vs30 / k1[ip - 1], n)) - log(A1100 + c));
		}
		else if (Vs30 > k1(ip))
		{
			f_siteG = (c11[ip - 1] + k2[ip - 1] * n) * log(Vs30 / k1[ip - 1]);
		}

		double f_siteJ;
		if (Vs30 <= 200)
			f_siteJ = (c12[ip - 1] + k2[ip - 1] * n) * (log(Vs30 / k1[ip - 1]) - log(200 / k1[ip - 1])) * Sj;
		else
			f_siteJ = (c13[ip - 1] + k2[ip - 1] * n) * log(Vs30 / k1[ip - 1]) * Sj;

		double fsite = f_siteG + f_siteJ;

		// Basin Response Term - Sediment effects

		double fsed;
		if (Z25 <= 1)
			fsed = (c14[ip - 1] + c15[ip - 1] * Sj) * (Z25 - 1.0);
		else if (Z25 <= 3)
			fsed = 0;
		else if (Z25 > 3)
			fsed = c16[ip - 1] * k3[ip - 1] * exp(-0.75) * (1.0 - exp(-0.25 * (Z25 - 3.0)));

		// Hypocenteral Depth term
		double f_hypH;
		if (Zhyp <= 7)
			f_hypH = 0;
		else if (Zhyp <= 20)
			f_hypH = Zhyp - 7.0;
		else
			f_hypH = 13;

		double f_hypM;
		if (M <= 5.5)
			f_hypM = c17[ip - 1];
		else if (M <= 6.5)
			f_hypM = c17[ip - 1] + (c18[ip - 1] - c17[ip - 1]) * (M - 5.5);
		else
			f_hypM = c18[ip - 1];

		double fhyp = f_hypH * f_hypM;

		// Fault Dip term
		double f_dip;
		if (M <= 4.5)
			f_dip = c19[ip - 1] * delta;
		else if (M <= 5.5)
			f_dip = c19[ip - 1] * (5.5 - M) * delta;
		else
			f_dip = 0;

		// Anelastic Attenuation Term
		double f_atn;
		if (Rrup > 80)
			f_atn = (c20(ip) + Dc20(ip)) * (Rrup - 80);
		else
			f_atn = 0;

		// Median value
		*Sa = exp(fmag + fdis + fflt + fhng + fsite + fsed + fhyp + f_dip + f_atn);

		// Standard deviation computations
		double tau_lny, tau_lnPGA, phi_lny, phi_lnPGA;
		if (M <= 4.5)
		{
			tau_lny = t1[ip - 1];
			tau_lnPGA = t1[22 - 1]; // ip = PGA
			phi_lny = f1[ip - 1];
			phi_lnPGA = f1[22 - 1];
		}
		else if (M < 5.5)
		{
			tau_lny = t2[ip - 1] + (t1[ip - 1] - t2[ip - 1]) * (5.5 - M);
			tau_lnPGA = t2[22 - 1] + (t1[22 - 1] - t2[22 - 1]) * (5.5 - M);// ip = PGA
			phi_lny = f2[ip - 1] + (f1[ip - 1] - f2[ip - 1]) * (5.5 - M);
			phi_lnPGA = f2[22 - 1] + (f1[22 - 1] - f2[22 - 1]) * (5.5 - M);
		}
		else
		{
			tau_lny = t2[ip - 1];
			tau_lnPGA = t2[22 - 1];
			phi_lny = f2[ip - 1];
			phi_lnPGA = f2[22 - 1];
		}

		double tau_lnyB = tau_lny;
		double tau_lnPGAB = tau_lnPGA;
		double phi_lnyB = sqrt(phi_lny * phi_lny - flnAF[ip - 1] * flnAF[ip - 1]);
		double phi_lnPGAB = sqrt(phi_lnPGA * phi_lnPGA - flnAF[ip - 1] * flnAF[ip - 1]);

		double alpha;
		if (Vs30 < k1[ip])
			alpha = k2[ip - 1] * A1100 * (1.0 / (A1100 + c * pow(Vs30 / k1[ip - 1], n)) - 1.0 / (A1100 + c));
		else
			alpha = 0;


		*tau = sqrt(tau_lnyB * tau_lnyB + alpha * alpha * tau_lnPGAB * tau_lnPGAB + 2.0 * alpha * rlnPGA_lnY[ip - 1] * tau_lnyB * tau_lnPGAB);
		double phi = sqrt(phi_lnyB * phi_lnyB + flnAF[ip - 1] * flnAF[ip - 1] + alpha * alpha * phi_lnPGAB * phi_lnPGAB + 2.0 * alpha * rlnPGA_lnY[ip - 1] * phi_lnyB * phi_lnPGAB);

		*sigma = sqrt((*tau) * (*tau) + phi * phi);
	}
};

