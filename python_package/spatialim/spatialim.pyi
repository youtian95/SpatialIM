"""
SpatialIM类型存根文件
提供IntelliSense支持
"""

from typing import List, Dict, Any, Optional
import numpy as np

def lamda_M(M: float) -> float:
    """震级与地震发生率的转换函数"""
    ...

class EQSource_CB14PCA:
    """地震源类，用于模拟地震动强度空间分布"""
    
    def __init__(self, lon_0: float, lat_0: float) -> None:
        """
        创建一个地震源，需要指定震中的经纬度(度)
        
        Args:
            lon_0: 震中经度
            lat_0: 震中纬度
        """
        ...
    
    def set_W(self, W: float) -> None:
        """设置断裂面宽度W (km)"""
        ...
    
    def set_length(self, length: float) -> None:
        """设置断裂面长度 (km)"""
        ...
    
    def set_RuptureNormal(self, x: float, y: float, z: float) -> None:
        """
        设置断裂面朝上的法线方向
        
        Args:
            x: 向东为x
            y: 向北为y  
            z: 向上为z
        """
        ...
    
    def set_lambda(self, lambda_: float) -> None:
        """设置rake角度(度)"""
        ...
    
    def set_Fhw(self, Fhw: int) -> None:
        """
        设置hanging wall效应
        
        Args:
            Fhw: 1-包含, 0-排除
        """
        ...
    
    def set_Zhyp(self, Zhyp: float) -> None:
        """设置震源深度(km)"""
        ...
    
    def set_region(self, region: int) -> None:
        """
        设置区域
        
        Args:
            region: 0=全球(含台湾), 1=加州, 3=中国或土耳其, 4=意大利
        """
        ...
    
    def set_nPCs(self, nPCs: int) -> None:
        """设置IM相关性PCA方法模拟考虑的主成分阶数"""
        ...
    
    def set_seed(self, seed: int) -> None:
        """设置随机数种子"""
        ...
    
    def register_site(
        self, 
        ID: int, 
        lon: float, 
        lat: float,
        elevation_km: float, 
        T0: float, 
        Vs30: float, 
        Z25: float
    ) -> None:
        """
        注册一个场地点
        
        Args:
            ID: 场地ID
            lon: 经度
            lat: 纬度
            elevation_km: 海拔高程(km)
            T0: 周期(秒)
            Vs30: 30米平均剪切波速(m/s)
            Z25: 剪切波速2.5km/s的深度(km)
        """
        ...
    
    def simulate_intensities(
        self, 
        magnitudes: List[float], 
        ifmedian: bool = False,
        OutputIMAllPeriods: bool = True
    ) -> None:
        """
        模拟一系列震级下的烈度分布
        
        Args:
            magnitudes: 震级列表
            ifmedian: 是否输出中位值
            OutputIMAllPeriods: 是否输出所有周期的强度
        """
        ...
    
    def save_im(self, filename: str) -> None:
        """保存模拟的烈度到文件"""
        ...
    
    def save_xy(self, filename: str) -> None:
        """保存场地坐标到文件"""
        ...
    
    def get_results(self) -> Dict[str, Any]:
        """获取模拟结果"""
        ...
