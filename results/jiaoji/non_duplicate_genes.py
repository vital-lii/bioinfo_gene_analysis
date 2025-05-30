#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
预处理cle4.csv文件，解决重复基因名的问题
"""

import pandas as pd
import numpy as np
import os
import argparse
from collections import defaultdict

def load_data(file_path):
    """加载CSV或Excel格式的基因表达数据"""
    # 根据文件扩展名选择读取方法
    file_ext = os.path.splitext(file_path)[1].lower()
    if file_ext == '.csv':
        print(f"检测到CSV文件，使用read_csv加载: {file_path}")
        df = pd.read_csv(file_path)
    else:
        print(f"检测到Excel文件，使用read_excel加载: {file_path}")
        df = pd.read_excel(file_path)
    
    print(f"从文件加载数据完成，形状: {df.shape}")
    print(f"列名: {df.columns.tolist()}")
    
    return df

def detect_duplicates(df, gene_column='gene_symb'):
    """检测DataFrame中的重复基因名"""
    
    # 确保指定的基因列存在
    if gene_column not in df.columns:
        print(f"错误: 未找到指定的基因列 '{gene_column}'")
        return None
    
    # 计数每个基因出现的次数
    gene_counts = df[gene_column].value_counts()
    
    # 找出出现多次的基因
    duplicate_genes = gene_counts[gene_counts > 1]
    
    if len(duplicate_genes) == 0:
        print("未发现重复基因名")
        return None
    
    print(f"发现{len(duplicate_genes)}个重复的基因名:")
    for gene, count in duplicate_genes.items():
        print(f"  - {gene}: 重复{count}次")
    
    return duplicate_genes

def handle_duplicates(df, method='merge', gene_column='gene_symb', value_columns=None):
    """
    处理DataFrame中的重复基因名
    
    参数:
    - df: 包含基因数据的DataFrame
    - method: 处理重复基因的方法，可选值：
        - 'merge': 合并重复基因行，对数值列取平均值
        - 'keep_first': 保留每个基因的第一次出现
        - 'suffix': 为重复基因添加后缀
        - 'best_pvalue': 保留P值最小的条目
    - gene_column: 基因名称所在的列名
    - value_columns: 需要特殊处理的值列名（如对数值列取平均等）
    
    返回:
    - 处理后的DataFrame
    """
    # 检测重复基因
    duplicate_genes = detect_duplicates(df, gene_column)
    
    if duplicate_genes is None or len(duplicate_genes) == 0:
        print("无需处理重复基因")
        return df
    
    # 如果没有指定值列，则使用所有数值列
    if value_columns is None:
        value_columns = df.select_dtypes(include=['number']).columns.tolist()
    
    # 为每个基因创建一个唯一ID，用于后续处理
    df['_temp_gene_id'] = range(len(df))
    
    # 根据选择的方法处理重复基因
    if method == 'merge':
        print("使用合并策略处理重复基因")
        # 创建一个新的DataFrame存储结果
        result_df = pd.DataFrame(columns=df.columns)
        unique_genes = df[gene_column].unique()
        
        for gene in unique_genes:
            gene_rows = df[df[gene_column] == gene]
            
            if len(gene_rows) == 1:
                # 如果只有一行，直接添加
                result_df = pd.concat([result_df, gene_rows])
            else:
                # 如果有多行，合并数据
                merged_row = gene_rows.iloc[0:1].copy()  # 创建基于第一行的副本
                
                # 对数值列取平均
                for col in value_columns:
                    if col in df.columns:
                        merged_row[col] = gene_rows[col].mean()
                
                # 对pvalue和padj列取最小值（更为显著的值）
                for col in ['pvalue', 'padj']:
                    if col in df.columns and not gene_rows[col].isna().all():
                        merged_row[col] = gene_rows[col].min()
                
                # 对log2FC使用与最显著P值对应的值
                if 'log2FC' in df.columns and 'padj' in df.columns:
                    # 检查padj列是否有非NaN值
                    if not gene_rows['padj'].isna().all():
                        most_sig_idx = gene_rows['padj'].idxmin()
                        merged_row['log2FC'] = df.loc[most_sig_idx, 'log2FC']
                    elif 'pvalue' in df.columns and not gene_rows['pvalue'].isna().all():
                        # 如果padj全为NaN但pvalue有值，则使用pvalue
                        most_sig_idx = gene_rows['pvalue'].idxmin()
                        merged_row['log2FC'] = df.loc[most_sig_idx, 'log2FC']
                    else:
                        # 两者都没有有效值，则使用log2FC的平均值
                        merged_row['log2FC'] = gene_rows['log2FC'].mean()
                
                # 合并字符串列内容（如果有的话）
                for col in df.select_dtypes(include=['object']).columns:
                    if col not in [gene_column, '_temp_gene_id']:  # 跳过基因名和临时ID列
                        unique_values = gene_rows[col].unique()
                        if len(unique_values) > 1:
                            merged_row[col] = '; '.join(str(x) for x in unique_values if pd.notna(x))
                
                result_df = pd.concat([result_df, merged_row])
        
        # 移除临时ID列
        result_df = result_df.drop('_temp_gene_id', axis=1)
        
        print(f"处理完成: 从{len(df)}行合并为{len(result_df)}行")
        return result_df
    
    elif method == 'keep_first':
        print("保留每个基因的第一次出现")
        # 依次处理每个重复基因，只保留第一次出现
        result_df = df.drop_duplicates(subset=[gene_column], keep='first')
        
        print(f"处理完成: 从{len(df)}行减少为{len(result_df)}行")
        return result_df.drop('_temp_gene_id', axis=1)
    
    elif method == 'suffix':
        print("为重复基因添加后缀")
        # 创建新的DataFrame和一个计数字典
        result_df = pd.DataFrame(columns=df.columns)
        gene_counter = defaultdict(int)
        
        # 重新组织数据，为重复基因添加后缀
        for _, row in df.iterrows():
            gene = row[gene_column]
            gene_counter[gene] += 1
            
            if gene_counter[gene] > 1:
                # 如果是重复基因，添加后缀
                new_row = row.copy()
                new_row[gene_column] = f"{gene}_{gene_counter[gene]}"
                result_df = pd.concat([result_df, pd.DataFrame([new_row])])
            else:
                # 如果是第一次出现，保持不变
                result_df = pd.concat([result_df, pd.DataFrame([row])])
        
        # 移除临时ID列
        result_df = result_df.drop('_temp_gene_id', axis=1)
        
        print(f"处理完成: 添加后缀后保留全部{len(df)}行数据")
        return result_df
    
    elif method == 'best_pvalue':
        print("保留每个基因中P值最小的条目")
        
        # 确保存在p值列
        p_col = None
        for col in ['padj', 'pvalue']:
            if col in df.columns:
                p_col = col
                break
        
        if p_col is None:
            print("未找到p值列，默认使用'keep_first'策略")
            return handle_duplicates(df, 'keep_first', gene_column, value_columns)
        
        # 对每个基因组，保留p值最小的行
        df_sorted = df.sort_values(by=[gene_column, p_col])
        result_df = df_sorted.drop_duplicates(subset=[gene_column], keep='first')
        
        # 移除临时ID列
        result_df = result_df.drop('_temp_gene_id', axis=1)
        
        print(f"处理完成: 从{len(df)}行减少为{len(result_df)}行")
        return result_df
    
    else:
        print(f"未知的处理方法: {method}，使用默认的'merge'策略")
        return handle_duplicates(df, 'merge', gene_column, value_columns)

def analyze_duplicates(df, gene_column='gene_symb'):
    """分析重复基因之间的差异程度"""
    
    # 检测重复基因
    duplicate_genes = detect_duplicates(df, gene_column)
    
    if duplicate_genes is None or len(duplicate_genes) == 0:
        print("无重复基因可分析")
        return
    
    # 数值列列表
    numeric_cols = df.select_dtypes(include=['number']).columns.tolist()
    
    # 分析重复基因之间的差异
    print("\n重复基因间的差异分析:")
    
    # 遍历每个重复基因
    for gene, count in duplicate_genes.items():
        gene_rows = df[df[gene_column] == gene]
        
        print(f"\n基因 {gene} (出现{count}次):")
        
        # 分析数值列
        for col in numeric_cols:
            values = gene_rows[col].values
            if len(values) > 0:
                min_val = np.min(values)
                max_val = np.max(values)
                mean_val = np.mean(values)
                std_val = np.std(values)
                
                # 计算变异系数(CV)，如果均值不为0
                cv = std_val / mean_val if mean_val != 0 else float('inf')
                
                # 只显示有显著差异的列
                if cv > 0.1:  # 变异系数大于10%
                    print(f"  - {col}: 范围[{min_val:.4f}, {max_val:.4f}], 均值={mean_val:.4f}, 标准差={std_val:.4f}, 变异系数={cv:.4f}")
        
        # 特别关注p值和log2FC
        for col in ['pvalue', 'padj', 'log2FC']:
            if col in df.columns:
                values = gene_rows[col].values
                if len(values) > 0:
                    min_val = np.min(values)
                    max_val = np.max(values)
                    print(f"  - {col}: 范围[{min_val:.6f}, {max_val:.6f}]")

def main():
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='预处理基因表达数据，解决重复基因名问题')
    
    # 添加参数
    parser.add_argument('-i', '--input', type=str, required=True, 
                        help='输入文件路径(CSV或Excel格式)')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='输出文件路径(将保存处理后的数据)')
    parser.add_argument('-m', '--method', type=str, default='merge',
                        choices=['merge', 'keep_first', 'suffix', 'best_pvalue'],
                        help='处理重复基因的方法: merge=合并取平均, keep_first=保留第一个, suffix=添加后缀, best_pvalue=保留最显著的')
    parser.add_argument('-g', '--gene_column', type=str, default='gene_symb',
                        help='基因名称所在的列名(默认: gene_symb)')
    parser.add_argument('-a', '--analyze', action='store_true',
                        help='分析重复基因之间的差异')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 加载数据
    df = load_data(args.input)
    
    # 分析重复基因差异
    if args.analyze:
        analyze_duplicates(df, args.gene_column)
    
    # 处理重复基因
    processed_df = handle_duplicates(df, args.method, args.gene_column)
    
    # 保存处理后的数据
    output_ext = os.path.splitext(args.output)[1].lower()
    if output_ext == '.csv':
        processed_df.to_csv(args.output, index=False)
    else:
        processed_df.to_excel(args.output, index=False)
    
    print(f"处理后的数据已保存至: {args.output}")
    print(f"处理前行数: {len(df)}, 处理后行数: {len(processed_df)}")

if __name__ == "__main__":
    main() 