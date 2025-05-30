#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
从cle4.csv文件中提取差异表达基因（reject_null为TRUE的行）
"""

import pandas as pd
import os
import argparse

def extract_diff_genes(input_file, output_dir):
    """提取差异表达基因
    
    参数:
        input_file (str): 输入文件路径，如cle4.csv
        output_dir (str): 输出目录
    """
    print(f"正在从{input_file}中读取数据...")
    # 读取CSV文件
    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        print(f"读取文件时出错: {e}")
        return
    
    print(f"原始数据形状: {df.shape}")
    
    # 确保reject_null列存在
    if 'reject_null' not in df.columns:
        print("错误: 未找到'reject_null'列")
        print(f"可用的列: {df.columns.tolist()}")
        return
    
    # 提取reject_null为TRUE的行
    # 注意: 可能的情况包括'TRUE', 'True', True，都做判断
    filtered_df = df[df['reject_null'].astype(str).str.upper() == 'TRUE']
    
    # 如果过滤结果为空，尝试布尔值形式
    if len(filtered_df) == 0:
        filtered_df = df[df['reject_null'] == True]
    
    print(f"过滤后数据形状: {filtered_df.shape}，找到{len(filtered_df)}个差异基因")
    
    # 创建输出目录（如果不存在）
    os.makedirs(output_dir, exist_ok=True)
    
    # 生成输出文件名
    output_file = os.path.join(output_dir, 'diff_genes.csv')
    
    # 保存结果
    filtered_df.to_csv(output_file, index=False)
    print(f"差异基因数据已保存至: {output_file}")
    
    # 可选：根据基因类型或表达方向进一步分类
    if 'gene_type' in filtered_df.columns:
        # 按gene_type分组
        gene_types = filtered_df['gene_type'].unique()
        print(f"发现{len(gene_types)}个不同的基因类型: {gene_types}")
        
        for gene_type in gene_types:
            type_df = filtered_df[filtered_df['gene_type'] == gene_type]
            type_file = os.path.join(output_dir, f'diff_genes_{gene_type}.csv')
            type_df.to_csv(type_file, index=False)
            print(f"类型'{gene_type}'的基因({len(type_df)}个)已保存至: {type_file}")
    
    # 可选：如果有log2FC列，可以分别提取上调和下调基因
    if 'log2FC' in filtered_df.columns:
        # 上调基因(log2FC > 0)
        up_genes = filtered_df[filtered_df['log2FC'] > 0]
        up_file = os.path.join(output_dir, 'up_regulated.csv')
        up_genes.to_csv(up_file, index=False)
        print(f"上调基因({len(up_genes)}个)已保存至: {up_file}")
        
        # 下调基因(log2FC < 0)
        down_genes = filtered_df[filtered_df['log2FC'] < 0]
        down_file = os.path.join(output_dir, 'down_regulated.csv')
        down_genes.to_csv(down_file, index=False)
        print(f"下调基因({len(down_genes)}个)已保存至: {down_file}")

def main():
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='提取差异表达基因')
    parser.add_argument('-i', '--input', required=True, help='输入CSV文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出目录路径')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 提取差异基因
    extract_diff_genes(args.input, args.output)

if __name__ == "__main__":
    main() 