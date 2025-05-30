#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
处理enriched_gene_symbols.csv文件，进行基因统计分析
作者：数据分析团队
日期：2023

使用方法:
    python process_enriched_genes.py --input <输入文件路径> --output <输出目录路径>
    
    例如:
    python process_enriched_genes.py --input data/my_gene_data.csv --output results/my_analysis
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
from collections import Counter
import re
import argparse

# 默认文件路径
DEFAULT_INPUT_FILE = "trial/microbiome_analysis/data/sage_da/enriched_gene_symbols.csv"
DEFAULT_OUTPUT_DIR = "trial/microbiome_analysis/results"

def parse_arguments():
    """
    解析命令行参数
    
    返回:
        解析后的参数对象
    """
    parser = argparse.ArgumentParser(description='处理基因数据文件并进行统计分析')
    
    parser.add_argument('--input', '-i', 
                      type=str, 
                      default=DEFAULT_INPUT_FILE,
                      help='输入文件路径, 默认为: ' + DEFAULT_INPUT_FILE)
    
    parser.add_argument('--output', '-o', 
                      type=str, 
                      default=DEFAULT_OUTPUT_DIR,
                      help='输出目录路径, 默认为: ' + DEFAULT_OUTPUT_DIR)
    
    return parser.parse_args()

def load_enriched_genes(file_path):
    """
    加载并清理enriched_gene_symbols.csv文件
    
    参数:
        file_path: 输入文件路径
    返回:
        处理后的DataFrame
    """
    print(f"正在加载文件: {file_path}")
    df = pd.read_csv(file_path)
    
    # 显示基本信息
    print(f"数据集形状: {df.shape}")
    print("数据集列名:", df.columns.tolist())
    
    # 数据清理
    # 统计每列的缺失值
    missing_values = df.isnull().sum()
    print("\n缺失值统计:")
    for col, count in missing_values.items():
        print(f"{col}: {count} ({count/len(df)*100:.2f}%)")
    
    return df

def analyze_gene_distribution(df, output_dir):
    """
    分析基因分布情况
    
    参数:
        df: 包含基因信息的DataFrame
        output_dir: 输出目录路径
    """
    # 统计有基因ID和基因符号的行数
    has_gene_id = df['gene_id'].notna().sum()
    has_gene_symbol = df['mapping_gene_symb'].notna().sum()
    
    print(f"\n有基因ID的行数: {has_gene_id} ({has_gene_id/len(df)*100:.2f}%)")
    print(f"有基因符号的行数: {has_gene_symbol} ({has_gene_symbol/len(df)*100:.2f}%)")
    
    # 按p值分组统计
    p_value_bins = [0, 0.001, 0.01, 0.05, 1.0]
    p_value_labels = ['<0.001', '0.001-0.01', '0.01-0.05', '>0.05']
    df['p_value_group'] = pd.cut(df['p_value'], bins=p_value_bins, labels=p_value_labels, right=False)
    
    p_value_counts = df['p_value_group'].value_counts().sort_index()
    print("\nP值分布统计:")
    for group, count in p_value_counts.items():
        print(f"{group}: {count} ({count/len(df)*100:.2f}%)")
    
    # 绘制P值分布图
    plt.figure(figsize=(10, 6))
    plt.bar(p_value_counts.index, p_value_counts.values)
    plt.title('P值分布统计')
    plt.xlabel('P值范围')
    plt.ylabel('基因数量')
    plt.savefig(f"{output_dir}/p_value_distribution.png", dpi=300)
    
    # 只对有基因符号的行进行后续分析
    genes_with_symbols = df.dropna(subset=['mapping_gene_symb'])
    
    # 分析一个基因ID可能对应多个基因符号的情况
    gene_symbols_list = []
    for symbols in genes_with_symbols['mapping_gene_symb']:
        # 处理可能有多个基因符号的情况（用分号分隔）
        if isinstance(symbols, str):
            gene_symbols_list.extend([s.strip() for s in symbols.split(';')])
    
    # 统计最常见的基因符号
    gene_symbol_counter = Counter(gene_symbols_list)
    most_common_genes = gene_symbol_counter.most_common(20)
    
    print("\n最常见的20个基因符号:")
    for gene, count in most_common_genes:
        print(f"{gene}: {count}次")
    
    # 绘制最常见基因符号的柱状图
    if most_common_genes:
        plt.figure(figsize=(12, 8))
        genes, counts = zip(*most_common_genes)
        plt.barh(list(reversed(genes)), list(reversed(counts)))
        plt.title('最常见的基因符号')
        plt.xlabel('出现次数')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/most_common_genes.png", dpi=300)

def analyze_p_value_vs_gene_info(df, output_dir):
    """
    分析p值与基因信息的关系
    
    参数:
        df: 包含基因信息的DataFrame
        output_dir: 输出目录路径
    """
    # 创建一个新列，标记是否有基因符号
    df['has_gene_symbol'] = df['mapping_gene_symb'].notna()
    
    # 比较有无基因符号的p值分布
    plt.figure(figsize=(10, 6))
    sns.kdeplot(data=df, x='p_value', hue='has_gene_symbol', common_norm=False)
    plt.title('有无基因符号的p值分布')
    plt.xlabel('p值')
    plt.ylabel('密度')
    plt.savefig(f"{output_dir}/p_value_by_gene_symbol.png", dpi=300)
    
    # 计算有无基因符号的p值统计量
    with_symbol = df[df['has_gene_symbol']]['p_value']
    without_symbol = df[~df['has_gene_symbol']]['p_value']
    
    print("\n有基因符号的p值统计:")
    print(f"数量: {len(with_symbol)}")
    print(f"平均值: {with_symbol.mean():.6f}")
    print(f"中位数: {with_symbol.median():.6f}")
    print(f"标准差: {with_symbol.std():.6f}")
    
    print("\n无基因符号的p值统计:")
    print(f"数量: {len(without_symbol)}")
    print(f"平均值: {without_symbol.mean():.6f}")
    print(f"中位数: {without_symbol.median():.6f}")
    print(f"标准差: {without_symbol.std():.6f}")

def export_filtered_genes(df, output_dir):
    """
    导出筛选后的基因列表
    
    参数:
        df: 包含基因信息的DataFrame
        output_dir: 输出目录路径
    """
    # 筛选有基因符号且p值显著的基因
    significant_genes = df[(df['mapping_gene_symb'].notna()) & (df['p_value'] < 0.05)]
    
    # 按p值排序
    significant_genes = significant_genes.sort_values('p_value')
    
    # 导出到CSV文件
    significant_genes.to_csv(f"{output_dir}/significant_genes.csv", index=False)
    print(f"\n已导出{len(significant_genes)}个显著性基因到 {output_dir}/significant_genes.csv")
    
    # 统计前20个显著基因
    top_genes = significant_genes.head(20)
    print("\n前20个显著性基因:")
    for i, (_, row) in enumerate(top_genes.iterrows(), 1):
        print(f"{i}. {row['mapping_gene_symb']} (Gene ID: {row['gene_id']}, p-value: {row['p_value']:.6f})")

def main():
    """主函数"""
    # 解析命令行参数
    args = parse_arguments()
    
    input_file = args.input
    output_dir = args.output
    
    # 创建输出目录（如果不存在）
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 80)
    print(f"开始处理基因数据文件: {input_file}")
    print(f"结果将保存到: {output_dir}")
    print("=" * 80)
    
    # 加载数据
    df = load_enriched_genes(input_file)
    
    # 分析基因分布
    analyze_gene_distribution(df, output_dir)
    
    # 分析p值与基因信息的关系
    analyze_p_value_vs_gene_info(df, output_dir)
    
    # 导出筛选后的基因
    export_filtered_genes(df, output_dir)
    
    print("=" * 80)
    print(f"分析完成！结果已保存到 {output_dir} 目录")
    print("=" * 80)

if __name__ == "__main__":
    main() 