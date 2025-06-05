"""
SpatialIM plotting module
Provides visualization functions for earthquake intensity spatial distribution
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any
import warnings
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.colors import Normalize
from scipy.interpolate import griddata


class SpatialIMPlotter:
    """SpatialIM result plotter"""
    
    def __init__(self):
        """Initialize the plotter"""
        # Set default font parameters
        plt.rcParams['axes.unicode_minus'] = False
    
    def _convert_fault_geometry_StrikeDip2RuptureNormal(self, strike: float, dip: float) -> Dict[str, np.ndarray]:
        """
        Convert fault strike and dip to geometric vectors
          Args:
            strike: Fault strike angle (degrees, 0=North, clockwise positive)
            dip: Fault dip angle (degrees, 90=vertical)
            
        Returns:
            Dictionary containing:
                - 'strike_vector': Unit vector along strike direction
                - 'dip_vector': Unit vector along dip direction  
                - 'normal_vector': Unit vector normal to fault plane
        """
        # Convert to radians
        strike_rad = np.radians(strike)
        dip_rad = np.radians(dip)
        
        # Strike vector (along fault length, horizontal)
        # Strike = 0° points North, increases clockwise
        strike_vector = np.array([
            np.cos(np.pi/2 - strike_rad),   # East component
            np.sin(np.pi/2 - strike_rad),   # North component  
            0.0                   # Vertical component
        ])
        
        # Dip vector (down-dip direction)
        # Points in the direction of maximum dip
        dip_azimuth = strike - 90.0  # Dip direction is perpendicular to strike
        dip_azimuth_rad = np.radians(dip_azimuth)
        
        dip_vector = np.array([
            np.cos(dip_azimuth_rad) * np.cos(dip_rad),  # East component
            np.sin(dip_azimuth_rad) * np.cos(dip_rad),  # North component
            -np.sin(dip_rad)                            # Vertical component (negative = downward)        
        ])
        
        # Normal vector to fault plane (rupture normal)
        # Pointing out of the fault surface
        normal_vector =  - np.cross(strike_vector, dip_vector)
        normal_vector = normal_vector / np.linalg.norm(normal_vector)  # Normalize
        
        return {
            'strike_vector': strike_vector,
            'dip_vector': dip_vector,
            'normal_vector': normal_vector
        }
    
    def _convert_fault_geometry_RuptureNormal2StrikeDip(self, rupture_normal: np.ndarray) -> Dict[str, float]:
        """
        Convert rupture normal vector to fault strike and dip angles
        
        Args:
            rupture_normal: Rupture normal vector as 3D numpy array [x, y, z]
            
        Returns:
            Dictionary containing:
                - 'strike': Fault strike angle (degrees, 0=North, clockwise positive)
                - 'dip': Fault dip angle (degrees, 90=vertical)
        """
        if rupture_normal.shape != (3,):
            raise ValueError("rupture_normal must be a 3D numpy array")
        
        # Normalize the rupture normal vector
        rupture_normal = rupture_normal / np.linalg.norm(rupture_normal)
        
        # Calculate dip from vertical component of normal
        dip_rad = np.arcsin(abs(rupture_normal[2]))
        dip = np.degrees(dip_rad)
        
        # Calculate strike from horizontal components
        horizontal_normal = rupture_normal[:2]
        if np.linalg.norm(horizontal_normal) > 1e-10:
            # Normalize horizontal component
            horizontal_normal = horizontal_normal / np.linalg.norm(horizontal_normal)
            horizontal_normal_rad = np.arctan2(horizontal_normal[1], horizontal_normal[0])  # atan2(y, x)
            # Strike is perpendicular to horizontal normal, and 90 degrees counter-clockwise from it
            strike_rad = horizontal_normal_rad + np.pi / 2.0
            # Adjust strike to be clockwise from North
            strike_rad = np.pi / 2.0 - strike_rad 
        else:
            # For vertical fault, strike can be arbitrary (set to 0)
            strike_rad = 0.0

        # Ensure positive angle and in [0, 360)
        strike = np.degrees(strike_rad)
        strike = strike % 360.0
        if strike < 0:
            strike += 360.0

        return {
            'strike': strike,
            'dip': dip
        }

    def _calculate_fault_projection(self, fault_width: float, fault_length: float, fault_strike: Optional[float] = None, fault_dip: Optional[float] = None, rupture_normal: Optional[np.ndarray] = None, epicenter_x: float = 0.0, epicenter_y: float = 0.0) -> List[Tuple[float, float]]:
        """
        Calculate fault surface projection rectangle corners
        
        Args:
            fault_width: Fault width (km)
            fault_length: Fault length (km) 
            fault_strike: Fault strike angle (degrees, 0=North, clockwise positive). If None, must provide rupture_normal
            fault_dip: Fault dip angle (degrees, 90=vertical). If None, must provide rupture_normal
            rupture_normal: Rupture normal vector as 3D numpy array [x, y, z]. If None, must provide fault_strike and fault_dip
            epicenter_x: Epicenter X coordinate (km)
            epicenter_y: Epicenter Y coordinate (km)
            
        Returns:
            List of (x, y) coordinates for rectangle corners
            
        Raises:
            ValueError: If neither (strike, dip) nor rupture_normal is provided
        """

        # Validate input parameters
        if rupture_normal is not None:
            # Convert rupture normal to strike and dip
            if isinstance(rupture_normal, np.ndarray) and rupture_normal.shape == (3,):
                fault_geometry = self._convert_fault_geometry_RuptureNormal2StrikeDip(rupture_normal)
                fault_strike = fault_geometry['strike']
                fault_dip = fault_geometry['dip']
                strike_rad = np.radians(fault_strike)
                dip_rad = np.radians(fault_dip)
            else:
                raise ValueError("rupture_normal must be a 3D numpy array")
        elif fault_strike is not None and fault_dip is not None:
            # Use provided strike and dip
            strike_rad = np.radians(fault_strike)
            dip_rad = np.radians(fault_dip)
        else:
            raise ValueError("Must provide either (fault_strike, fault_dip) or rupture_normal")
        
        # Calculate half dimensions
        half_length = fault_length / 2
        half_width = fault_width / 2
        
        # For surface projection, we need to consider the dip effect
        # The width on surface = fault_width * cos(dip)
        surface_width = fault_width * np.cos(dip_rad)
        half_surface_width = surface_width / 2
        
        # Calculate rectangle corners relative to epicenter
        # Strike direction (along fault length)
        strike_x = np.cos( - strike_rad + np.pi/2.0)
        strike_y = np.sin( - strike_rad + np.pi/2.0)
        
        # Perpendicular to strike (across fault width)
        perp_x = np.cos(strike_rad)
        perp_y = - np.sin(strike_rad)
        
        # Four corners of the fault rectangle
        corners = []
        for length_sign in [-1, 1]:
            for width_sign in [-1, 1]:
                x = (epicenter_x + 
                     length_sign * half_length * strike_x + 
                     width_sign * half_surface_width * perp_x)
                y = (epicenter_y + 
                     length_sign * half_length * strike_y + 
                     width_sign * half_surface_width * perp_y)
                corners.append((x, y))
        
        # Reorder corners to form a proper rectangle path
        # corners[0]: -length, -width
        # corners[1]: -length, +width  
        # corners[2]: +length, -width
        # corners[3]: +length, +width
        ordered_corners = [corners[0], corners[1], corners[3], corners[2], corners[0]]  # Close the rectangle
        
        return ordered_corners
    
    def plot_intensity_contour(
        self,
        im_file: str,
        coords_file: str,
        epicenter_x: float = 0.0,
        epicenter_y: float = 0.0,
        title: Optional[str] = None,
        save_path: Optional[str] = None,
        figsize: Tuple[int, int] = (12, 5),
        dpi: int = 300,
        grid_resolution: int = 100, # number of grid points for interpolation
        show_points: bool = False,
        show_epicenter: bool = True,
        colormap: str = 'viridis',        # 断层参数
        fault_width: Optional[float] = None,
        fault_length: Optional[float] = None,
        fault_strike: Optional[float] = None,  # 断层走向角度 (度)
        fault_dip: Optional[float] = None,    # 断层倾角 (度)
        rupture_normal: Optional[np.ndarray] = None,  # 破裂面法向量
        show_fault: bool = False
    ) -> Figure:
        """
        Plot earthquake intensity distribution contour
        
        Args:
            im_file: Intensity file path
            coords_file: Coordinates file path
            epicenter_x: Epicenter x
            epicenter_y: Epicenter y
            title: Plot title
            save_path: Save path
            figsize: Figure size
            dpi: Image resolution
            grid_resolution: Grid resolution
            show_points: Whether to show data points
            show_epicenter: Whether to show epicenter
            colormap: Color mapping
            fault_width: Fault width (km)
            fault_length: Fault length (km)
            fault_strike: Fault strike angle (degrees, 0=North, clockwise positive). If None, must provide rupture_normal
            fault_dip: Fault dip angle (degrees, 90=vertical). If None, must provide rupture_normal
            rupture_normal: Rupture normal vector as 3D numpy array [x, y, z]. If None, must provide fault_strike and fault_dip
            show_fault: Whether to show fault projection
            
        Returns:
            matplotlib.figure.Figure: Figure object
        """

        # 读取数据
        im_data = pd.read_csv(im_file, sep=r'\s+', header=None, names=['ID', 'Intensity'])
        coord_data = pd.read_csv(coords_file, sep=r'\s+', header=None, names=['ID', 'X', 'Y'])
        merged_data = pd.merge(im_data, coord_data, on='ID')
        
        x = merged_data['X'].values
        y = merged_data['Y'].values
        z = merged_data['Intensity'].values
        
        # 创建网格
        x_grid = np.linspace(x.min(), x.max(), grid_resolution)
        y_grid = np.linspace(y.min(), y.max(), grid_resolution)
        X, Y = np.meshgrid(x_grid, y_grid)
        
        # 插值
        Z = griddata((x, y), z, (X, Y), method='cubic')
        
        # 创建图形
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        
        # 散点图
        scatter = ax1.scatter(x, y, c=z, cmap=colormap, s=20, alpha=0.8, edgecolors='white', linewidth=0.5)
        ax1.set_xlabel('X (km)')
        ax1.set_ylabel('Y (km)')
        ax1.set_title(f'Intensity Measure - scatter')
        ax1.grid(True, alpha=0.3)
        ax1.set_aspect('equal')
        plt.colorbar(scatter, ax=ax1, label='Intensity')
        
        # 云图
        contour = ax2.contourf(X, Y, Z, levels=20, cmap=colormap, alpha=0.8)
        contour_lines = ax2.contour(X, Y, Z, levels=10, colors='black', alpha=0.6, linewidths=0.5)
        ax2.clabel(contour_lines, inline=True, fontsize=8, fmt='%.3f')
        
        if show_points:
            ax2.scatter(x, y, c='red', s=10, alpha=0.6, marker='x')
        
        ax2.set_xlabel('X (km)')
        ax2.set_ylabel('Y (km)')
        ax2.set_title(f'Intensity Measure - contour')
        ax2.grid(True, alpha=0.3)
        ax2.set_aspect('equal')
        plt.colorbar(contour, ax=ax2, label='Intensity')
        
        # 添加震中
        if show_epicenter:
            ax1.plot(epicenter_x, epicenter_y, 'r*', markersize=15, label='Epicenter', markeredgecolor='white', markeredgewidth=1)
            ax2.plot(epicenter_x, epicenter_y, 'r*', markersize=15, label='Epicenter', markeredgecolor='white', markeredgewidth=1)
            ax1.legend()
            ax2.legend()
        
        # 添加断层投影矩形
        if show_fault and fault_width is not None and fault_length is not None:
            fault_corners = self._calculate_fault_projection(
                fault_width=fault_width, 
                fault_length=fault_length,
                fault_strike=fault_strike,
                fault_dip=fault_dip,
                rupture_normal=rupture_normal,
                epicenter_x=epicenter_x, 
                epicenter_y=epicenter_y
            )
            
            # 提取x和y坐标
            fault_x = [corner[0] for corner in fault_corners]
            fault_y = [corner[1] for corner in fault_corners]
            
            # 分别绘制矩形的各条边
            # 矩形顶点顺序: [左上, 左下, 右下, 右上, 左上(闭合)]
            # 下边缘 (左下到右下): corners[1] to corners[2] - 使用虚线
            ax1.plot([fault_x[1], fault_x[2]], [fault_y[1], fault_y[2]], 'b--', linewidth=2, alpha=0.8)
            ax2.plot([fault_x[1], fault_x[2]], [fault_y[1], fault_y[2]], 'b--', linewidth=2, alpha=0.8)

            # 其他边缘 (右下到右上, 右上到左上，左上到左下): corners[2] to corners[3], corners[3] to corners[0], corners[0] to corners[1]
            ax1.plot([fault_x[2], fault_x[3]], [fault_y[2], fault_y[3]], 'b-', linewidth=2, alpha=0.8)
            ax2.plot([fault_x[2], fault_x[3]], [fault_y[2], fault_y[3]], 'b-', linewidth=2, alpha=0.8)
            ax1.plot([fault_x[3], fault_x[0]], [fault_y[3], fault_y[0]], 'b-', linewidth=2, alpha=0.8)
            ax2.plot([fault_x[3], fault_x[0]], [fault_y[3], fault_y[0]], 'b-', linewidth=2, alpha=0.8)
            ax1.plot([fault_x[0], fault_x[1]], [fault_y[0], fault_y[1]], 'b-', linewidth=2, alpha=0.8)
            ax2.plot([fault_x[0], fault_x[1]], [fault_y[0], fault_y[1]], 'b-', linewidth=2, alpha=0.8)
            
            # 更新图例
            if show_epicenter:
                ax1.legend()
                ax2.legend()
            else:
                ax1.legend(['Fault rupture projection'])
                ax2.legend(['Fault rupture projection'])
        
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        
        # 保存图像
        if save_path:
            plt.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
            print(f"Image saved to {save_path}")
        
        # 显示图像
        plt.show()
        
        return fig


