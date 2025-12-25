import spatialim
import os
import sys

def main():
    # 获取当前脚本所在目录的绝对路径
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 定义测试数据路径
    eq_source_path = os.path.join(current_dir, "fixtures", "EQSource.json")
    site_file_path = os.path.join(current_dir, "fixtures", "SiteFile.csv")
    
    print(f"正在运行 spatial_im 演示...")
    print(f"震源文件: {eq_source_path}")
    print(f"场地文件: {site_file_path}")

    # 检查文件是否存在
    if not os.path.exists(eq_source_path):
        print(f"错误: 找不到震源文件: {eq_source_path}")
        return
    if not os.path.exists(site_file_path):
        print(f"错误: 找不到场地文件: {site_file_path}")
        return

    try:
        # 调用 Rust 编写的模块
        # 参数: 震源文件路径, 场地文件路径, GMPE模型(可选, 默认CB14)
        spatialim.run_simulation(
            eq_source_path, 
            site_file_path,
            None
        )
        print("模拟运行成功完成！")
        
    except RuntimeError as e:
        print(f"运行时发生错误: {e}")
    except Exception as e:
        print(f"发生未知错误: {e}")

if __name__ == "__main__":
    main()
