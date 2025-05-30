import pandas as pd
import sys
import os

def extract_significant_rows(file_path, p_value_column='p_value', threshold=0.05):
    """
    提取P值小于阈值的行
    
    参数:
    file_path: CSV文件路径
    p_value_column: 包含P值的列名，默认为'p_value'
    threshold: P值阈值，默认为0.05
    
    返回:
    significant_df: 仅包含显著性行的数据框，并保存到文件
    """
    try:
        # 读取CSV文件
        print(f"正在读取文件: {file_path}")
        df = pd.read_csv(file_path)
        
        # 确保P值列存在
        if p_value_column not in df.columns:
            print(f"错误: 列 '{p_value_column}' 在CSV文件中不存在")
            print("可用的列:", list(df.columns))
            return None
        
        # 将P值列转换为数值类型（处理可能的字符串表示）
        df[p_value_column] = pd.to_numeric(df[p_value_column], errors='coerce')
        
        # 提取P值小于阈值的行
        significant_df = df[df[p_value_column] < threshold]
        
        # 计算总行数和符合条件的行数
        total_rows = len(df)
        significant_rows = len(significant_df)
        
        # 输出结果到文件
        output_path = os.path.splitext(file_path)[0] + f'_p_less_than_{threshold}.csv'
        significant_df.to_csv(output_path, index=False)
        
        print(f"文件总行数: {total_rows}")
        print(f"P值<{threshold}的行数: {significant_rows} ({significant_rows/total_rows*100:.2f}%)")
        print(f"结果已保存到: {output_path}")
        
        return significant_df
    
    except Exception as e:
        print(f"处理文件时出错: {e}")
        return None

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("使用方法: python extract_significant_pvalues.py <csv文件路径> [P值列名] [阈值]")
        sys.exit(1)
    
    file_path = sys.argv[1]
    p_column = sys.argv[2] if len(sys.argv) > 2 else 'p_value'
    threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
    
    extract_significant_rows(file_path, p_column, threshold) 