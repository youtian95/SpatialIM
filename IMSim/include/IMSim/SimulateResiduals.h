//采用PCA方法模拟多个场地、不同周期下谱加速度Sa(T)的within-event residuals
//这些within-event residuals是相关的
//
//使用方法：注册场地->模拟
//
// 参考文献：
// Markhvida M, Ceferino L, Baker JW.Modeling spatially correlated spectral 
// accelerations at multiple periods using principal component analysis and 
// geostatistics.Earthquake Eng Struc. 2018; 47:1107 - 23. 
// https ://doi.org/10.1002/eqe.3007

#pragma once
#define _USE_MATH_DEFINES 
#include <map>
#include <vector>
#include <cmath>
#include <random>
#include <assert.h>
#include "Eigen/Dense"
#include "Eigen/Geometry"

using namespace std;
using namespace Eigen;

class SimulateResiduals
{
    struct ModelVario
    {
        float Cn = 0;
        float C1 = 0;
        float a1 = 0;
        float C2 = 0;
        float a2 = 0;
        string Type = "nug";
    };

public:
    int nPCs = 10;
	SimulateResiduals()
	{
        T_ << 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5;
        PCAcoefs_ << 0.270963956240108, -0.139418157111539, 0.0690420060759825, -0.106094866059058, -0.0922880747536406, -0.11348997612886, -0.188935371416588, 0.153956801827048, -0.160082932478587, -0.0485878661500656, 0.106169114128661, 0.0545367125260001, -0.0842347288819985, 0.00206507178166717, 0.233666515544382, -0.0444106080542835, -0.298766213268524, -0.527588859528322, -0.580349072958488,
            0.270185457409379, -0.141734438837191, 0.0770156687411159, -0.11639353409999, -0.103464378298138, -0.124082463392091, -0.199840301009291, 0.155452551301531, -0.157024101304166, -0.0511781532022337, 0.102685985245256, 0.053409178198703, -0.0785807732589693, 0.00538047460862058, 0.220317828901249, -0.039452593127257, -0.257172668625987, -0.15099488917774, 0.781868928002583,
            0.266716484131893, -0.150918021372557, 0.101241750461614, -0.144620230365225, -0.128327845057672, -0.150413273486784, -0.21751911461417, 0.15453312842206, -0.14455513348748, -0.0493784913365888, 0.0865380808702005, 0.0369034473261574, -0.0551975904000634, 0.00787249481774084, 0.149651509850332, -0.0232087969820605, -0.0284597879541604, 0.808901723567414, -0.226437334732122,
            0.251688240452975, -0.184642998841202, 0.178879968100826, -0.221328310597117, -0.17555752575313, -0.176668887066346, -0.1886513513666, 0.0424749058524636, -0.0455090101766787, -0.0291885726697761, -0.0315764520880226, -0.0605411454951765, 0.0935508235649957, 0.0224423933374744, -0.299181352126819, 0.0599168858817312, 0.754350402692246, -0.206472972801382, 0.023109344507298,
            0.236434540660266, -0.218922079202474, 0.23725418394185, -0.234559034122274, -0.133267087790244, -0.0431828093963935, 0.119447151365187, -0.272310554909179, 0.238192698302355, 0.100676333288767, -0.263315034080161, -0.121207177353811, 0.202769694434657, 0.00661936093850581, -0.493306767118467, 0.116246172684757, -0.475923822815555, 0.0367733765077595, -0.00643955732220411,
            0.232994643271854, -0.228087987254185, 0.230554572919478, -0.16044302411133, 0.0400564218872398, 0.181657484726543, 0.427112684239685, -0.324579279763868, 0.263780433255153, 0.142634796459082, -0.0813780347714383, 0.0465305509298986, -0.151801546733413, -0.0833183140197916, 0.534198178434453, -0.184596749983802, 0.210357916588565, -0.0028422556361022, 0.00331108666921881,
            0.238919759244457, -0.211905954003063, 0.132646222385762, 0.0820453502922968, 0.32794697293789, 0.393273105011282, 0.325836316409324, 0.162029624546621, -0.182164846060428, -0.138319895254022, 0.470111475270552, 0.177876087511963, -0.1112565000941, 0.0883177907230135, -0.291143253148395, 0.262494244929844, -0.00152291509197, 0.0154770015927024, 0.00150080927064067,
            0.247247200513419, -0.174053609784519, -0.00825743327819653, 0.277382297057791, 0.403271333843203, 0.220437620135438, -0.0837312940200531, 0.224796020087722, -0.171941472633753, -0.0292747130158189, -0.38152412132245, -0.237244495814365, 0.356271008578838, -0.0850494996785505, -0.0125434811410974, -0.442559354727306, 0.0151125187962444, 0.0110964548458275, 0.00194187520085148,
            0.253677096925723, -0.122375885147353, -0.148595585559558, 0.365271223153003, 0.253186443805678, -0.0612442510569459, -0.283389735157485, -0.0811437931924326, 0.212210668293107, 0.14336242671851, -0.275727618896954, -0.0411307164923357, -0.202014525375018, 0.022845903851548, 0.155257213551213, 0.632145331896681, 0.0455130187869974, 0.000690596761983894, 0.000469634928915885,
            0.254921191707501, -0.071319446353373, -0.237030888421264, 0.359073100317565, 0.0401080106567685, -0.248766601944991, -0.14185903779119, -0.286692239119504, 0.300971237941779, 0.0579993716202019, 0.328411991836789, 0.208361703130758, -0.194768507716927, 0.0324295008946133, -0.258822160901248, -0.477244327079575, 0.00191094743796522, 0.00662185258908696, 0.000188266912551749,
            0.252458254214951, 0.0125091293591534, -0.32712108086482, 0.226053196913636, -0.26129762020473, -0.216236975254179, 0.344080559560279, -0.121230620609225, -0.0602714322805403, -0.219189670580381, 0.211470671417798, -0.128634841134357, 0.576234521196323, -0.0550760430415841, 0.197333043823494, 0.20576645538163, 0.023662144998231, 0.00370370109158968, 0.000185197342559174,
            0.24594424065373, 0.0799604139803053, -0.358449872812475, 0.064099810496083, -0.341792253779399, 0.022496754565775, 0.388717982330275, 0.177122203733985, -0.255990758059885, -0.00644303562267801, -0.37539012296473, -0.0762061466572891, -0.502002035655987, 0.0183662355165662, -0.176351451500561, -0.0686345486566678, 0.0154233464320777, 0.00541822807684994, 0.00127868850390311,
            0.225758567264058, 0.191381035473567, -0.335176303224685, -0.216152633771253, -0.165178954634192, 0.423011571619087, -0.144255461671014, 0.187567292296388, 0.149360081833295, 0.530105300516035, 0.041796232067731, 0.326784764375549, 0.27460983615727, 0.0558724755215897, 0.00442737635788848, 0.0111352503084675, 0.0244045338837114, 0.000682889258092345, 0.00019866537927149,
            0.211097168683674, 0.259405649803992, -0.243643584807687, -0.325745718814885, 0.0763484285808871, 0.330279571319951, -0.220001696329049, -0.117395307490011, 0.271296590718256, -0.438774328341411, 0.148165305851337, -0.48454567943612, -0.143691433970383, -0.0389147466608019, 0.00893381613114769, -0.0205549231060488, -0.00611584426449596, 0.00380246484862132, -0.000389229000503269,
            0.188387436860863, 0.329799740851314, -0.0946692611641819, -0.273646378894757, 0.356651426359913, -0.153161129681719, -0.000682050969765767, -0.329897663448327, -0.267361749795915, -0.2793821559276, -0.26374050060882, 0.528987613257165, 0.0703281456138742, -0.0835258438506085, -0.0254953444060441, 0.0271139318663258, 0.0126194670232453, 0.00385357399236438, -0.000809502038431328,
            0.176395533106338, 0.357332294420881, 0.0554387508851632, -0.155161957798969, 0.354513035089274, -0.343041979311133, 0.161895880026785, -0.0275355960813715, -0.20785940912098, 0.506833313170992, 0.205450556183935, -0.413916086467618, -0.0407497250821277, 0.168472840930613, -0.00235465742383666, -0.00588317616785561, -0.00334222371106789, 0.00197868686362638, 0.00170392027445259,
            0.165469018345669, 0.360040619637271, 0.260392170105291, 0.0669803672081155, 0.0572701975680616, -0.220913199055637, 0.181062550106522, 0.519913777311151, 0.462086625105155, -0.104489884552583, -0.0219632037127792, 0.119011798891668, -0.00429988165275597, -0.417545575417374, -0.0401081045672218, 0.020060720387271, -0.00526025219733044, -0.00580167364426092, 0.00045823003408634,
            0.159580892256856, 0.347927159738324, 0.348346295207843, 0.239394140984691, -0.157211928082521, 0.0928779108728929, -0.00501602103568442, 0.0169759687957175, 0.109978222336614, -0.182603153808833, -0.121233489669211, 0.0711306930722663, 0.0620582733857731, 0.750450107378881, 0.0785420660886401, -0.0521102089376762, 0.0075193664688849, -0.00188268149870413, -0.00185190445021904,
            0.1488329207305, 0.33284769539218, 0.36509805168606, 0.3312276948947, -0.281334582456685, 0.283334016381439, -0.182848344502475, -0.325817997238157, -0.310946153889112, 0.128556113560638, 0.0837173663270596, -0.0703828760736262, -0.046192412575096, -0.441499963638099, -0.0405906261145499, 0.033716161820098, 0.00116291918248088, 0.0040160503301927, 0.000505570344725481;
        variance_scale_factor_ << 0.639841744964052, 0.846277138923352, 0.904533062923065, 0.933402823301903, 0.950154585181316, 0.960461560434709, 0.967972143777431, 0.973878205021816, 0.979297028245927, 0.983345801486615, 0.986639491844098, 0.989681397350252, 0.992367578022788, 0.994790440764041, 0.996929078974732, 0.998759738778621, 0.999764879971799, 0.999969807143844, 1;
        // Cn   C1  a1  C2  a2  Type
        modelVario_[0] = { 2.500000000000001, 4.520000000000002, 15, 6.780000000000003, 250, "iso nest" };
        modelVario_[1] = { 0.500000000000000, 1.400000000000000, 10, 2.600000000000001, 160, "iso nest" };
        modelVario_[2] = { 0.150000000000000, 0.420000000000000, 15, 0.630000000000000, 160, "iso nest" };
        modelVario_[3] = { 0.150000000000000, 0.225000000000000, 10, 0.225000000000000, 120, "iso nest" };
        modelVario_[4] = { 0.314321867136085,0,0,0,0, "nug" };
        modelVario_[5] = { 0.190749535511365,0,0,0,0, "nug" };
        modelVario_[6] = { 0.137846758971697,0,0,0,0, "nug" };
        modelVario_[7] = { 0.111283843493347,0,0,0,0, "nug" };
        modelVario_[8] = { 0.096499281204429,0,0,0,0, "nug" };
        modelVario_[9] = { 0.071736796680044,0,0,0,0, "nug" };
        modelVario_[10] = { 0.064816215163269,0,0,0,0, "nug" };
        modelVario_[11] = { 0.054076635653157,0,0,0,0, "nug" };
        modelVario_[12] = { 0.051188751201166,0,0,0,0, "nug" };
        modelVario_[13] = { 0.043316419835114,0,0,0,0, "nug" };
        modelVario_[14] = { 0.041398046055658,0,0,0,0, "nug" };
        modelVario_[15] = { 0.034663671578289,0,0,0,0, "nug" };
        modelVario_[16] = { 0.018796994070527,0,0,0,0, "nug" };
        modelVario_[17] = { 0.002856941187582,0,0,0,0, "nug" };
        modelVario_[18] = { 3.606453961200486e-04,0,0,0,0, "nug" };
    }
    void set_random_engine(default_random_engine* p)
    {
        pRND = p;
    }
    void RegisterSites(const VectorXf& x, const VectorXf& y)
        //注册场地, 乔斯基分解
    {
        VectorXf variance_scale_factor = variance_scale_factor_;
        vector<ModelVario> modelVario = modelVario_;

        int nLocs = x.size();

        // Create a matrix containing distances between locations
        MatrixXf distanceMatrix = getDistanceMatrix(x, y);

        // Scale variance if less than 19 principal components are used
        if (nPCs < 19)
        {
            for (int i = 0; i < modelVario.size(); i++)
            {
                if (modelVario[i].Type == "nug")
                {
                    modelVario[i].Cn = modelVario[i].Cn / variance_scale_factor[nPCs - 1];
                }
                else
                {
                    modelVario[i].Cn = modelVario[i].Cn / variance_scale_factor(nPCs);
                    modelVario[i].C1 = modelVario[i].C1 / variance_scale_factor(nPCs);
                    modelVario[i].C2 = modelVario[i].C2 / variance_scale_factor(nPCs);
                }
            }
        }

        L = vector<MatrixXf>(nPCs);

        // Create a covariance matricies for each of the principal components(PC's)
        vector<MatrixXf> covMatrix = vector<MatrixXf>(nPCs);
        for (int i = 0; i < nPCs; i++)
        {
            cout << i << endl;
            covMatrix[i] = getCovariance(distanceMatrix, modelVario[i]);

            //乔斯基分解

            Eigen::MatrixXf normTransform(nLocs, nLocs);

            Eigen::LLT<Eigen::MatrixXf> cholSolver(covMatrix[i]);

            // We can only use the cholesky decomposition if 
            // the covariance matrix is symmetric, pos-definite.
            // But a covariance matrix might be pos-semi-definite.
            // In that case, we'll go to an EigenSolver
            if (cholSolver.info() == Eigen::Success) {
                // Use cholesky solver
                normTransform = cholSolver.matrixL();
            }
            else {
                cout << "协方差矩阵的cholesky分解失败！" << endl;
                // Use eigen solver
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigenSolver(covMatrix[i]);
                normTransform = eigenSolver.eigenvectors()
                    * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
            }

            L[i] = normTransform;
        }
    }
    vector<MatrixXf> simulateResiduals(const VectorXf& T_sim, int nsims) const
        // Simulation of within - event ground motion residuals at multiple periods using PCA
        //
        // Further methodology documentation : M.Markhvida, L.Ceferino, and J.Baker.Modeling
        // spatially correlated spectral accelerations at multiple periods using
        // principal component analysisand geostatistics.Submitted to : Earthquake Engineering
        // & Structural Dynamics In Review(2017).
        //
        // INPUT
        // x, y     = x and y coordinates(in kilometers) for n locations of interest with size[n x 1]
        // T_sim    = a vector specifying periods to be simulated of size[1 x t]
        // nsims    = number of simulation to be performed
        // nPCs     = number of principal components to be used in simulation
        //              (recommended not less than 5)
        //
        // OUTPUT
        // sim_results  = [n x t] matrix of spatially cross - correlated within - event
        //                  normal residuals(mean 0 and variance 1) for n locations
        //                  and t periods
    {
        VectorXf T = T_;
        MatrixXf PCAcoefs = PCAcoefs_;

        int nLocs = L[0].rows();

        // Simulate each of the PC's
        vector<MatrixXf> sim_PCA = vector<MatrixXf>(nPCs);
        for (int i = 0; i < sim_PCA.size(); i++)
        {
            sim_PCA[i] = MatrixXf::Zero(nLocs, nsims);
        }
        VectorXf mu = VectorXf::Zero(nLocs);
        std::normal_distribution<> dist;
        for (int i_PC = 0; i_PC < nPCs; i_PC++)
        {
            MatrixXf samples(nsims, mu.size());
            for (int i_sim = 0; i_sim < nsims; i_sim++)
            {
                VectorXf randN(nLocs);
                for (int j = 0; j < nLocs; j++)
                {
                    randN(j) = dist(*pRND);
                }
                samples.row(i_sim) = L[i_PC] * randN + mu;
            }

            sim_PCA[i_PC] = samples.transpose();
        }

        // Transform simulated PC's to spectral acceleration residuals for desired periods

        VectorXi indexT = findmember(T_sim, T);
        vector<MatrixXf> sim_results(nsims);
        for (int i = 0; i < nsims; i++)
        {
            sim_results[i] = MatrixXf::Zero(nLocs, T_sim.size());
        }

        for (int i = 0; i < indexT.size(); i++)
        {
            for (int j = 0; j < nsims; j++)
            {
                MatrixXf temp_sim_PCA(nLocs, nPCs);
                for (int row = 0; row < nLocs; row++)
                {
                    for (int col = 0; col < nPCs; col++)
                    {
                        temp_sim_PCA(row, col) = sim_PCA[col](row, j);
                    }
                }
                if (indexT[i])
                {
                    MatrixXf temp = PCAcoefs.row(indexT[i] - 1).head(nPCs);
                    sim_results[j].col(i) = PCA_T_X(temp_sim_PCA,
                        temp,
                        VectorXf::Zero(1));
                }
                else
                {
                    MatrixXf extraPCAcoefs(1, nPCs);
                    for (int k = 0; k < nPCs; k++)
                    {
                        VectorXf PCAcoefs_T = PCAcoefs.col(k);
                        extraPCAcoefs(1, k) = interp1(T, PCAcoefs_T, T_sim[i]);
                    }
                    sim_results[j].col(i) = PCA_T_X(temp_sim_PCA, extraPCAcoefs, VectorXf::Zero(1));
                }
            }
        }
        return sim_results;
    }

private:
    default_random_engine* pRND;
    vector<MatrixXf> L; //每个PCA分量的 乔斯基分解 LL^T = COV
    VectorXf T_ = VectorXf(19);
    Matrix<float, 19, 19> PCAcoefs_;
    VectorXf variance_scale_factor_ = VectorXf(19);
    vector<ModelVario> modelVario_ = vector<ModelVario>(19);

    static MatrixXf getCovariance(const MatrixXf& DIST, ModelVario modelVario)
    {
        MatrixXf COV;
        if (modelVario.Type == "iso nest")
        {
            COV = getIsoNestedCOV(modelVario, DIST);
        }
        else if (modelVario.Type == "nug")
        {
            COV = getNugCOV(modelVario, DIST);
        }
        else
        {
            perror("NOT RIGHT MODEL");
        }
        return COV;
    }
    static MatrixXf getIsoNestedCOV(ModelVario varModel, const MatrixXf& DIST)
    {
        MatrixXf COV;
        float var = varModel.Cn + varModel.C1 + varModel.C2;
        COV = var - (varModel.Cn + varModel.C1 * (1.0 - (-3.0 * DIST / varModel.a1).array().exp()) + varModel.C2 * (1.0 - (-3.0 * DIST / varModel.a2).array().exp()));
        for (int row = 0; row < DIST.rows(); row++)
        {
            for (int col = 0; col < DIST.cols(); col++)
            {
                if (DIST(row, col) == 0)
                {
                    COV(row, col) = var;
                }
            }
        }
        return COV;
    }
    static MatrixXf getNugCOV(ModelVario varModel, const MatrixXf& DIST)
    {
        MatrixXf COV;
        COV = MatrixXf::Zero(DIST.rows(), DIST.cols());
        for (int row = 0; row < DIST.rows(); row++)
        {
            for (int col = 0; col < DIST.cols(); col++)
            {
                if (DIST(row, col) == 0)
                {
                    COV(row, col) = varModel.Cn;
                }
            }
        }
        return COV;
    }
    static MatrixXf getDistanceMatrix(const VectorXf& x, const VectorXf& y)
    {
        assert(x.size() == y.size());
        //MatrixXf DIST = MatrixXf::Zero(x.size(), x.size());

        MatrixXf X = x * (VectorXf::Ones(x.size()).transpose())
            - VectorXf::Ones(x.size()) * (x.transpose());
        MatrixXf Y = y * (VectorXf::Ones(y.size()).transpose())
            - VectorXf::Ones(y.size()) * (y.transpose());
        MatrixXf DIST = (X.array().square() + Y.array().square()).cwiseSqrt();

        /*for (int i = 0; i < x.size(); i++)
        {
            for (int j = 0; j < x.size(); j++)
            {
                DIST(i, j) = sqrt(pow(x[i] - x[j], 2.0) + pow(y[i] - y[j], 2.0));
                DIST(j, i) = DIST(i, j);
            }
        }*/
        return DIST;
    }
    static MatrixXf PCA_T_X(const MatrixXf& transformed, const MatrixXf& coef, const VectorXf& mu)
        // Function that takes from transformed space back to the original space
        // Input:
        // coef          Coefficient Matrix(txn)
        // mu            Mean of the original space(1xt)
        // transformed   Matrix in the transformed space(mxn)
        // Output :
        // original      Matrix in the original space(mxt)
    {
        MatrixXf original(transformed.rows(), mu.size());
        MatrixXf temp = transformed * (coef.transpose());
        for (int row = 0; row < original.rows(); row++)
        {
            original.row(row) = temp.row(row) + mu.transpose();
        }
        return original;
    }
    static VectorXi findmember(const VectorXf& x_sim, const VectorXf& x_all)
        //寻找x_sim在x_all中的位置, 从1开始
    {
        VectorXi result(x_sim.size());
        for (int i = 0; i < x_sim.size(); i++)
        {
            int index = 0;
            for (int j = 0; j < x_all.size(); j++)
            {
                if (x_all[j] == x_sim[i])
                {
                    index = j + 1;
                    break;
                }
            }
            result[i] = index;
        }
        return result;
    }
    static double interp1(VectorXf x, VectorXf y, double vx)
        //插值, x单调
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
};

