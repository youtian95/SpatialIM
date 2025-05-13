"""
SpatialIM: 空间相关的烈度分布模拟库
"""

# 从预编译的PYD模块导入所有内容
try:
    from .spatialim import *
except ImportError:
    # 如果导入失败，显示友好的错误信息
    import platform, sys
    python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
    system = platform.system()
    bit = platform.architecture()[0]
    
    error_msg = f"""
    无法导入SpatialIM模块。可能的原因：
    1. 您的Python版本({python_version})不是3.12
    2. 您的操作系统({system})不是Windows
    3. 您的系统架构({bit})不是64位
    
    此包仅支持Windows 64位平台上的Python 3.12。
    """
    raise ImportError(error_msg)
