import pandas as pd
import argparse
import os

# -r 25844 自己指定保留多少行
def clean_data(file_path, max_rows=25844, filter_zeros=True):
    """
    读取XLSX文件并只保留指定行数据，可选择性过滤零值
    
    Parameters:
    -----------
    file_path : str
        XLSX文件的路径
    max_rows : int
        要保留的最大行数
    filter_zeros : bool
        是否过滤掉全为零的行
        
    Returns:
    --------
    pd.DataFrame
        清理后的数据框
    """
    print(f"读取文件: {file_path}")
    # 读取XLSX文件
    df = pd.read_excel(file_path)
    
    total_rows = df.shape[0]
    print(f"原始数据: {total_rows} 行, {df.shape[1]} 列")
    
    # 只保留前max_rows行
    if max_rows and max_rows < total_rows:
        df_cleaned = df.iloc[:max_rows].copy()
        print(f"保留前 {max_rows} 行")
    else:
        df_cleaned = df.copy()
        print("保留所有行")
    
    # 提取计数列
    count_cols = [col for col in df_cleaned.columns if col.startswith('GF') or col.startswith('SCFA')]
    
    if count_cols:
        print(f"发现 {len(count_cols)} 个计数列: {', '.join(count_cols)}")
        
        # 如果需要过滤零值行
        if filter_zeros:
            print("开始过滤全零行...")
            before_filter = df_cleaned.shape[0]
            
            # 检查所有计数列是否都为零
            is_all_zeros = (df_cleaned[count_cols] == 0).all(axis=1)
            df_cleaned = df_cleaned[~is_all_zeros]
            
            print(f"过滤后: 从 {before_filter} 行减少到 {df_cleaned.shape[0]} 行")
            print(f"移除了 {before_filter - df_cleaned.shape[0]} 行全零数据")
    else:
        print("警告: 未找到计数列 (GF* 或 SCFA*)")
    
    # 检查数据基本信息
    print("\n数据统计:")
    for col in count_cols:
        zeros = (df_cleaned[col] == 0).sum()
        zeros_pct = zeros / len(df_cleaned) * 100
        print(f"列 {col}: {zeros} 个零值 ({zeros_pct:.2f}%)")
    
    return df_cleaned

def save_data(df, output_path, file_name=None):
    """
    保存清理后的数据
    
    Parameters:
    -----------
    df : pd.DataFrame
        数据框
    output_path : str
        输出路径 (文件或目录)
    file_name : str, optional
        如果output_path是目录，则用此文件名
    """
    # 确保输出目录存在
    if os.path.isdir(output_path) or not output_path.endswith(('.xlsx', '.csv')):
        # 输出路径是目录
        os.makedirs(output_path, exist_ok=True)
        if file_name is None:
            file_name = 'cleaned_data.xlsx'
        full_path = os.path.join(output_path, file_name)
    else:
        # 输出路径直接是文件
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        full_path = output_path
    
    # 根据文件扩展名保存
    if full_path.endswith('.csv'):
        df.to_csv(full_path, index=False)
        print(f"数据已保存为CSV: {full_path}")
    else:
        df.to_excel(full_path, index=False)
        print(f"数据已保存为Excel: {full_path}")
    
    return full_path

if __name__ == "__main__":
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='清理RNA-Seq数据XLSX文件')
    
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='输入Excel文件路径')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='输出文件或目录路径')
    parser.add_argument('-r', '--rows', type=int, default=25844,
                        help='要保留的最大行数（默认：25844）')
    parser.add_argument('-n', '--name', type=str, default=None,
                        help='输出文件名（当输出是目录时使用）')
    parser.add_argument('--keep-zeros', action='store_true',
                        help='保留全为零的行（默认会过滤掉）')
    
    args = parser.parse_args()
    
    # 清理数据
    cleaned_df = clean_data(args.input, args.rows, not args.keep_zeros)
    
    # 保存清理后的数据
    output_file = save_data(cleaned_df, args.output, args.name)
    
    print(f"\n数据清理完成!")
    print(f"原始数据: {args.input}")
    print(f"清理后数据: {output_file}")
    print(f"保留的行数: {len(cleaned_df)}")
