import os
import sys
from pathlib import Path

class PathManager:
    """管理文件路径的类"""
    
    def __init__(self, output_dir=None):
        """
        初始化路径管理器
        
        Args:
            output_dir (str, optional): 自定义输出目录。默认为scripts目录下的output。
        """
        # 获取脚本所在目录
        self.script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        
        # 设置基础目录
        if output_dir:
            self.output_dir = output_dir
        else:
            self.output_dir = os.path.join(self.script_dir, "output")
        
        # 设置日志目录
        self.logs_dir = os.path.join(self.script_dir, "logs")
        
        # 定义路径映射
        self.paths = {
            'script': self.script_dir,
            'output': self.output_dir,
            'logs': self.logs_dir
        }
    
    def ensure_dirs(self):
        """确保所有需要的目录都存在"""
        for dir_path in self.paths.values():
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)
                print(f"创建目录: {dir_path}")
    
    def get_path(self, key):
        """
        获取指定的路径
        
        Args:
            key (str): 路径键名 ('script', 'output', 'logs')
            
        Returns:
            str: 对应的路径
        """
        if key in self.paths:
            return self.paths[key]
        else:
            raise KeyError(f"未找到路径键: {key}") 