#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "IMSim/EQSource_CB14PCA.h"

// 声明lamda_M函数，它在IMSim.cpp中定义
double lamda_M(double M);

namespace py = pybind11;

PYBIND11_MODULE(spatialim, m) {
    m.doc() = "SpatialIM: 空间相关的烈度分布模拟库"; // 模块文档字符串
    
    // 导出EQSource_CB14PCA类
    py::class_<EQSource_CB14PCA>(m, "EQSource_CB14PCA")
        .def(py::init<double, double>(), py::arg("lon_0"), py::arg("lat_0"), "创建一个地震源，需要指定震中的经纬度(度)")
        
        // 设置断裂面参数
        .def("set_W", &EQSource_CB14PCA::set_W, py::arg("W"), "设置断裂面宽度W")
        .def("set_length", &EQSource_CB14PCA::set_length, py::arg("length"), "设置断裂面长度")
        .def("set_RuptureNormal", &EQSource_CB14PCA::set_RuptureNormal, py::arg("x"), py::arg("y"), py::arg("z"), "设置断裂面朝上的法线方向(向东为x,向北为y,向上为z)")
        .def("set_lambda", &EQSource_CB14PCA::set_lambda, py::arg("lambda_"), "设置rake角度(度)")
        .def("set_Fhw", &EQSource_CB14PCA::set_Fhw, py::arg("Fhw"), "设置hanging wall效应, 1: 包含, 0: 排除")
        .def("set_Zhyp", &EQSource_CB14PCA::set_Zhyp, py::arg("Zhyp"), "设置震源深度(km)")
        .def("set_region", &EQSource_CB14PCA::set_region, py::arg("region"), "设置区域: 0=全球(含台湾), 1=加州, 3=中国或土耳其, 4=意大利")
        .def("set_nPCs", &EQSource_CB14PCA::set_nPCs, py::arg("nPCs"), "设置IM相关性PCA方法模拟考虑的主成分阶数")
        
        // 随机数种子相关
        .def("set_seed", [](EQSource_CB14PCA &self, int seed) {
            static std::default_random_engine engine(seed);
            self.set_randomengine(&engine);
        }, py::arg("seed"), "设置随机数种子")
        
        // 注册场地
        .def("register_site", &EQSource_CB14PCA::register_site,
             py::arg("ID"), py::arg("lon"), py::arg("lat"), 
             py::arg("elevation_km"), py::arg("T0"), 
             py::arg("Vs30"), py::arg("Z25"),
             "注册一个场地点")
        
        // 模拟烈度 - 修复参数注释不匹配的问题
        .def("simulate_intensities", &EQSource_CB14PCA::SimulateIntensities,
             py::arg("magnitudes"), py::arg("ifmedian") = false, py::arg("OutputIMAllPeriods") = true,
             "模拟一系列震级下的烈度分布")
        
        // 输出结果
        .def("save_im", &EQSource_CB14PCA::io_IM, py::arg("filename"),
             "保存模拟的烈度到文件")
        .def("save_xy", &EQSource_CB14PCA::io_XY, py::arg("filename"),
             "保存场地坐标到文件")
          // 添加获取结果的方法(需要根据您的类实现添加)
        .def("get_results", [](EQSource_CB14PCA &self) {
            // 这里需要根据您的类结构实现获取结果的方法
            // 返回一个合适的Python结构，如字典或NumPy数组
            py::dict results;
            // TODO: 从self中提取结果数据并填充到results中
            return results;
        }, "获取模拟结果");
        
    // 导出辅助函数
    m.def("lamda_M", &lamda_M, py::arg("M"),
          "震级与地震发生率的转换函数");
    // 添加示例函数
    m.def("simulate_earthquake", [](
            double lon_0, double lat_0,    // 震中经纬度
            double M,                      // 震级
            int N_sim,                     // 模拟次数
            int seed,                      // 随机数种子
            double W,                      // 断裂面宽度
            double length,                 // 断裂面长度
            double normal_x, double normal_y, double normal_z, // 法线方向
            double lambda,                 // rake角度
            int Fhw,                       // hanging wall效应
            double Zhyp,                   // 震源深度
            int region,                    // 区域
            int nPCs,                      // 主成分数
            bool ifmedian,                 // 是否输出中位值
            py::list sites                 // 场地列表，每个场地为[ID, lon, lat, elevation_km, T0, Vs30, Z25]
        ) {
            // 创建地震源
            EQSource_CB14PCA eqs(lon_0, lat_0);
            
            // 设置随机数引擎
            std::default_random_engine p(seed);
            eqs.set_randomengine(&p);
            
            // 设置断裂面参数
            eqs.set_W(W);
            eqs.set_length(length);
            eqs.set_RuptureNormal(normal_x, normal_y, normal_z);
            eqs.set_lambda(lambda);
            eqs.set_Fhw(Fhw);
            eqs.set_Zhyp(Zhyp);
            eqs.set_region(region);
            eqs.set_nPCs(nPCs);
            
            // 注册场地
            for (auto site : sites) {
                py::list s = site.cast<py::list>();
                eqs.register_site(
                    s[0].cast<int>(),      // ID
                    s[1].cast<double>(),   // lon
                    s[2].cast<double>(),   // lat
                    s[3].cast<double>(),   // elevation_km
                    s[4].cast<double>(),   // T0
                    s[5].cast<double>(),   // Vs30
                    s[6].cast<double>()    // Z25
                );
            }
            
            // 创建震级列表
            std::vector<double> M_list(N_sim, M);
            
            // 模拟烈度
            eqs.SimulateIntensities(M_list, ifmedian);
            
            // 返回结果
            py::dict results;
            // TODO: 从eqs中提取结果数据并填充到results中
            // 这需要根据EQSource_CB14PCA类的实现添加获取结果的代码
            
            return results;
        },
        py::arg("lon_0"), py::arg("lat_0"),
        py::arg("M"),
        py::arg("N_sim"),
        py::arg("seed"),
        py::arg("W"),
        py::arg("length"),
        py::arg("normal_x"), py::arg("normal_y"), py::arg("normal_z"),
        py::arg("lambda"),
        py::arg("Fhw"),
        py::arg("Zhyp"),
        py::arg("region"),
        py::arg("nPCs"),
        py::arg("ifmedian") = false,
        py::arg("sites") = py::list(),
        "模拟地震烈度分布的便捷函数"
    );
}