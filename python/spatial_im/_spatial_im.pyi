from typing import Optional

def run_simulation(
    eq_source_path: str, 
    site_file_path: str, 
    gmpe_model: Optional[str] = None
) -> None:
    """
    运行模拟的主函数。

    Args:
        eq_source_path (str): 震源文件路径 (.json)
        site_file_path (str): 场地文件路径 (.csv)
        gmpe_model (Optional[str], optional): GMPE 模型名称，默认为 "CB14"。
    """
    ...
