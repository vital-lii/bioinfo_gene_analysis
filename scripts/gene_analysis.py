import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import argparse

# 设置绘图风格
plt.style.use('seaborn-v0_8-whitegrid')  
plt.rcParams['font.sans-serif'] = ['SimHei']  
plt.rcParams['axes.unicode_minus'] = False    

def load_data(file_path):
    """加载CSV或Excel格式的基因表达数据"""
    # 根据文件扩展名选择读取方法
    file_ext = os.path.splitext(file_path)[1].lower()
    if file_ext == '.csv':
        print(f"检测到CSV文件，使用read_csv加载")
        df = pd.read_csv(file_path)
    else:
        print(f"检测到Excel文件，使用read_excel加载")
        df = pd.read_excel(file_path)
    
    print(f"从文件加载数据完成，形状: {df.shape}")
    print(f"列名: {df.columns.tolist()}")
    
    return df

def analyze_differential_expression(df_counts, pvalue_threshold=0.05, fc_threshold=1.0):
    """进行差异表达分析，特别处理包含P值列的数据"""
    # 获取列名，检查数据类型
    columns = df_counts.columns.tolist()
    
    print(f"DEBUG: 可用列名: {columns}")
    
    # 检查是否包含常见的基因表达和P值列
    has_means = all(col in columns for col in ['GF_mean', 'SCFA_mean'])
    has_pvalues = 'pvalue' in columns or 'padj' in columns
    has_log2fc = 'log2FC' in columns
    
    if has_means and has_pvalues and has_log2fc:
        print("检测到完整的表达数据（包含均值, P值和log2FC列）")
        
        # 创建结果数据框
        df_stats = pd.DataFrame(index=df_counts.index)
        
        # 复制已有列
        for col in columns:
            df_stats[col] = df_counts[col]
        
        # 使用已有p值列
        if 'pvalue' in columns:
            print(f"使用原始文件中的pvalue列")
        else:
            print(f"未在原始文件中找到pvalue列，尝试添加默认值")
            df_stats['pvalue'] = 0.5
            
        # 使用已有校正p值列
        if 'padj' in columns:
            print(f"使用原始文件中的padj列")
        else:
            print(f"未在原始文件中找到padj列，使用pvalue列替代")
            df_stats['padj'] = df_stats['pvalue']
        
        # 确保log2FC列存在
        if 'log2FC' not in df_stats.columns:
            print("计算log2FC值")
            df_stats['log2FC'] = np.log2((df_stats['SCFA_mean'] + 1) / (df_stats['GF_mean'] + 1))
        
        # 确保数据类型正确
        for col in ['padj', 'pvalue', 'log2FC']:
            if col in df_stats.columns:
                print(f"{col}列类型: {df_stats[col].dtype}")
                if df_stats[col].dtype == 'object':
                    print(f"将{col}列转换为数值类型")
                    df_stats[col] = pd.to_numeric(df_stats[col], errors='coerce')
        
        # 设置显著性和上下调状态
        df_stats['significant'] = (df_stats['padj'] < pvalue_threshold)
        df_stats['regulation'] = 'Not significant'
        df_stats.loc[(df_stats['significant']) & (df_stats['log2FC'] > fc_threshold), 'regulation'] = 'Up'
        df_stats.loc[(df_stats['significant']) & (df_stats['log2FC'] < -fc_threshold), 'regulation'] = 'Down'
        
        if 'gene_type' in columns:
            print("检测到gene_type列，将用于后续分析")
            df_stats['is_specific'] = df_stats['gene_type'].str.contains('specific', na=False)
        
        return df_stats
    else:
        print("未检测到完整的表达数据，将进行简化分析")
        # 这里可以添加简化分析的逻辑，如果需要的话
        # 但对于当前数据，已有P值，不需要这部分
        return pd.DataFrame()

def create_visualizations(df_counts, df_stats, output_dir, top_genes=50, custom_genes=None):
    """创建多种可视化图表，安全处理混合数据类型"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 确保数据为数值类型
    df_counts_numeric = df_counts.select_dtypes(include=['number'])
    
    # 检查数据列以便后续分析
    print("可用数值列:", df_counts_numeric.columns.tolist())
    
    # 1. 火山图 
    if 'log2FC' in df_stats.columns and 'pvalue' in df_stats.columns:
        plt.figure(figsize=(10, 8))
        
        # 处理p值，避免log10(0)的问题
        min_pval = 1e-10
        plot_pval = df_stats['pvalue'].copy()
        plot_pval[plot_pval < min_pval] = min_pval  # 设置最小p值
        
        # 按调控状态分组绘图
        if 'regulation' in df_stats.columns:
            colors = {'Up': 'red', 'Down': 'blue', 'Not significant': 'gray'}
            sizes = {'Up': 60, 'Down': 60, 'Not significant': 30}
            alphas = {'Up': 0.8, 'Down': 0.8, 'Not significant': 0.4}
            
            # 首先绘制不显著的点，然后是下调和上调点（确保重要点在上层）
            for reg in ['Not significant', 'Down', 'Up']:
                subset = df_stats[df_stats['regulation'] == reg]
                if not subset.empty:
                    plt.scatter(subset['log2FC'], 
                              -np.log10(plot_pval[subset.index]),
                              c=colors[reg], 
                              s=sizes[reg], 
                              alpha=alphas[reg], 
                              label=reg)
        else:
            # 简单绘制所有点
            plt.scatter(df_stats['log2FC'], -np.log10(plot_pval), alpha=0.7)
        
        # 添加显著性和fold change阈值线
        plt.axhline(-np.log10(0.05), linestyle='--', color='gray', alpha=0.5, label='p = 0.05')
        plt.axvline(-1, linestyle='--', color='gray', alpha=0.5)
        plt.axvline(1, linestyle='--', color='gray', alpha=0.5)
        
        # 标注最显著的基因（分别从上调和下调中选择）
        try:
            # 选择最显著的上调和下调基因（使用pvalue排序）
            top_up = df_stats[(df_stats['regulation'] == 'Up')].nsmallest(8, 'pvalue')
            top_down = df_stats[(df_stats['regulation'] == 'Down')].nsmallest(8, 'pvalue')
            
            # 合并并标注
            for gene_idx in pd.concat([top_up, top_down]).index:
                x_val = df_stats.loc[gene_idx, 'log2FC']
                y_val = -np.log10(plot_pval[gene_idx])
                plt.annotate(str(gene_idx), 
                           (float(x_val), float(y_val)),
                           fontsize=8,
                           xytext=(5, 5),
                           textcoords='offset points',
                           bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
        except Exception as e:
            print(f"标注基因名时出错: {str(e)}，跳过基因标注")
        
        plt.xlabel('log2 Fold Change (SCFA/GF)')
        plt.ylabel('-log10 P-value')
        plt.ylim(0, 4)  # 设置y轴范围为0-4
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.title('Volcano Plot: SCFA vs GF')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'volcano_plot.png'), dpi=300, bbox_inches='tight')
        plt.close()
    
    # 2. 热图 
    if all(col in df_counts_numeric.columns for col in ['GF_mean', 'SCFA_mean']):
        # 选择最显著的基因，限制数量为指定的top_genes
        df_stats_sorted = df_stats.sort_values('padj')
        
        # 检测重复的基因名并去除重
        non_duplicate_genes = []
        seen_genes = set()
        
        
        for gene in df_stats_sorted.index:
            if gene not in seen_genes:
                non_duplicate_genes.append(gene)
                seen_genes.add(gene)
        
        top_sig_genes = non_duplicate_genes[:top_genes]
        print(f"已选择前{len(top_sig_genes)}个唯一基因（按P值排序）")
        
        # 准备热图数据列
        cols_for_heatmap = ['GF_mean', 'SCFA_mean']
        if 'log2FC' in df_counts_numeric.columns:
            cols_for_heatmap.append('log2FC')
        
        # 创建热图数据框
        try:
            # 创建一个空的DataFrame来存储热图数据
            heatmap_data = pd.DataFrame(index=[], columns=cols_for_heatmap)
            
            # 添加存在的基因及其数据
            existing_genes = []
            for gene in top_sig_genes:
                if gene in df_counts_numeric.index:
                    existing_genes.append(gene)
            
            if existing_genes:
                # 直接从原始数据提取所需基因的行和列，确保结果是2D
                heatmap_data = df_counts_numeric.loc[existing_genes, cols_for_heatmap].copy()
                
                # 确保数据是数值型
                for col in heatmap_data.columns:
                    heatmap_data[col] = pd.to_numeric(heatmap_data[col], errors='coerce')
                
                # 对非log2FC列进行对数转换
                print(f"热图数据形状: {heatmap_data.shape}")
                for col in heatmap_data.columns:
                    if col != 'log2FC' and (heatmap_data[col] > 0).any():
                        print(f"对{col}列进行log2(x+1)转换")
                        heatmap_data[col] = np.log2(heatmap_data[col] + 1)
                
                # 绘制热图 - 检查数据维度
                print(f"热图数据最终形状: {heatmap_data.shape}, 类型: {type(heatmap_data)}")
                if heatmap_data.shape[0] > 0 and heatmap_data.shape[1] > 0:
                    plt.figure(figsize=(10, 12))
                    sns.heatmap(heatmap_data, cmap='viridis')
                    plt.title(f'Top {len(existing_genes)} Genes by P-value')
                    plt.xticks(rotation=45, ha='right')
                    plt.tight_layout()
                    plt.savefig(os.path.join(output_dir, 'heatmap_top_genes.png'), dpi=300, bbox_inches='tight')
                    plt.close()
                else:
                    print("警告: 热图数据形状无效，无法创建热图")
            else:
                print("警告: 未找到有效的基因数据，无法创建热图")
        except Exception as e:
            print(f"创建热图时发生错误: {str(e)}")
            import traceback
            traceback.print_exc()
            
        # 直接为上下调基因创建热图
        try:
            # 获取上调和下调基因
            up_genes_all = df_stats[df_stats['regulation'] == 'Up'].index.tolist()
            down_genes_all = df_stats[df_stats['regulation'] == 'Down'].index.tolist()
            
            # 去重 - 使用字典来保持顺序，同时去除重复项
            up_genes = []
            up_genes_seen = set()
            for gene in up_genes_all:
                if gene not in up_genes_seen:
                    up_genes.append(gene)
                    up_genes_seen.add(gene)
            
            down_genes = []
            down_genes_seen = set()
            for gene in down_genes_all:
                if gene not in down_genes_seen:
                    down_genes.append(gene)
                    down_genes_seen.add(gene)
            
            # 限制每组至多25个基因，避免图表过大
            up_genes = up_genes[:25]
            down_genes = down_genes[:25]
            
            # 为所有调控基因创建唯一列表
            regulation_genes = list(set(up_genes + down_genes))
            print(f"找到{len(up_genes)}个唯一上调基因和{len(down_genes)}个唯一下调基因")
            
            # 检查是否有足够的基因
            if len(regulation_genes) > 0:
                # 创建调控基因热图数据
                valid_genes = [g for g in regulation_genes if g in df_counts_numeric.index]
                
                if valid_genes:
                    # 直接提取有效基因的数据
                    reg_data = df_counts_numeric.loc[valid_genes, cols_for_heatmap].copy()
                    
                    # 确保数据是数值型并进行对数转换
                    for col in reg_data.columns:
                        reg_data[col] = pd.to_numeric(reg_data[col], errors='coerce')
                        if col != 'log2FC' and (reg_data[col] > 0).any():
                            reg_data[col] = np.log2(reg_data[col] + 1)
                    
                    # 创建调控信息列表
                    regulation_info = []
                    for gene in reg_data.index:
                        if gene in up_genes:
                            regulation_info.append('Up')
                        elif gene in down_genes:
                            regulation_info.append('Down')
                        else:
                            regulation_info.append('Not significant')
                    
                    # 打印调试信息
                    print(f"调控基因热图数据形状: {reg_data.shape}")
                    print(f"调控信息长度: {len(regulation_info)}")
                    
                    # 绘制热图 - 使用普通heatmap代替clustermap以简化代码
                    plt.figure(figsize=(10, len(valid_genes) * 0.4 + 3))
                    
                    # 设置上下调基因的颜色
                    cmap = plt.cm.viridis
                    
                    # 绘制基本热图
                    ax = sns.heatmap(reg_data, cmap=cmap)
                    
                    # 添加颜色条表示上下调
                    # 在y轴左侧添加颜色条
                    for i, status in enumerate(regulation_info):
                        color = 'red' if status == 'Up' else 'blue' if status == 'Down' else 'gray'
                        ax.add_patch(plt.Rectangle((-0.5, i), 0.3, 1, color=color, clip_on=False, transform=ax.transData))
                    
                    plt.title('Differentially Expressed Genes')
                    plt.tight_layout()
                    plt.savefig(os.path.join(output_dir, 'heatmap_regulated_genes.png'), dpi=300, bbox_inches='tight')
                    plt.close()
                    
                    # 尝试创建单独的上调和下调基因热图 - 确保无重复
                    if len(up_genes) > 0:
                        valid_up = [g for g in up_genes if g in df_counts_numeric.index]
                        if valid_up:
                            # 确保没有重复项
                            valid_up = list(dict.fromkeys(valid_up))
                            print(f"绘制{len(valid_up)}个上调基因的热图")
                            
                            up_data = df_counts_numeric.loc[valid_up, cols_for_heatmap].copy()
                            for col in up_data.columns:
                                if col != 'log2FC' and (up_data[col] > 0).any():
                                    up_data[col] = np.log2(up_data[col] + 1)
                            plt.figure(figsize=(8, len(valid_up) * 0.4 + 3))
                            sns.heatmap(up_data, cmap='Reds')
                            plt.title('Up-regulated Genes')
                            plt.tight_layout()
                            plt.savefig(os.path.join(output_dir, 'heatmap_upregulated.png'), dpi=300, bbox_inches='tight')
                            plt.close()
                    
                    if len(down_genes) > 0:
                        valid_down = [g for g in down_genes if g in df_counts_numeric.index]
                        if valid_down:
                            # 确保没有重复项
                            valid_down = list(dict.fromkeys(valid_down))
                            print(f"绘制{len(valid_down)}个下调基因的热图")
                            
                            down_data = df_counts_numeric.loc[valid_down, cols_for_heatmap].copy()
                            for col in down_data.columns:
                                if col != 'log2FC' and (down_data[col] > 0).any():
                                    down_data[col] = np.log2(down_data[col] + 1)
                            plt.figure(figsize=(8, len(valid_down) * 0.4 + 3))
                            sns.heatmap(down_data, cmap='Blues')
                            plt.title('Down-regulated Genes')
                            plt.tight_layout()
                            plt.savefig(os.path.join(output_dir, 'heatmap_downregulated.png'), dpi=300, bbox_inches='tight')
                            plt.close()
                else:
                    print("警告: 没有找到有效的调控基因数据")
            else:
                print("没有找到显著的上调或下调基因")
        except Exception as e:
            print(f"创建调控基因热图时发生错误: {str(e)}")
            import traceback
            traceback.print_exc()
    
    # 3. 柱状图 - 显示上下调基因数量
    if 'regulation' in df_stats.columns:
        plt.figure(figsize=(10, 6))
        reg_counts = df_stats['regulation'].value_counts()
        reg_counts.plot(kind='bar', color=['gray', 'blue', 'red'])
        plt.ylabel('Number of Genes')
        plt.title('Differentially Expressed Genes')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'regulation_counts.png'), dpi=300, bbox_inches='tight')
        plt.close()
    
    # 4. 如果有special基因标识，创建特殊基因的图
    if 'gene_type' in df_stats.columns:
        plt.figure(figsize=(12, 8))
        type_counts = df_stats['gene_type'].value_counts()
        
        # 只显示频率最高的前10种类型
        if len(type_counts) > 10:
            type_counts = type_counts.head(10)
        
        type_counts.plot(kind='bar')
        plt.ylabel('Number of Genes')
        plt.title('Gene Type Distribution')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'gene_type_distribution.png'), dpi=300, bbox_inches='tight')
        plt.close()
    
    return "可视化完成，结果保存在 " + output_dir

def main(file_path, output_dir, pvalue_threshold=0.05, fc_threshold=1.0, top_genes=50, custom_genes=None):
    """主函数：执行完整的数据分析和可视化流程"""
    # 1. 加载数据
    df = load_data(file_path)
    
    # 2. 设置基因名为索引
    gene_col = 'gene_symb' if 'gene_symb' in df.columns else df.columns[0]
    df_counts = df.set_index(gene_col)
    
    # 3. 差异表达分析
    df_stats = analyze_differential_expression(df_counts, pvalue_threshold, fc_threshold)
    
    # 4. 创建可视化
    create_visualizations(df_counts, df_stats, output_dir, top_genes=top_genes, custom_genes=custom_genes)
    
    # 5. 保存分析结果
    result_file = os.path.join(output_dir, 'differential_expression_results.csv')
    df_stats.to_csv(result_file)
    print(f"分析结果已保存至 {result_file}")
    
    # 6. 返回差异表达结果的摘要
    sig_up = sum((df_stats['padj'] < pvalue_threshold) & (df_stats['log2FC'] > fc_threshold))
    sig_down = sum((df_stats['padj'] < pvalue_threshold) & (df_stats['log2FC'] < -fc_threshold))
    
    print(f"padj<{pvalue_threshold}的基因数量: {sum(df_stats['padj'] < pvalue_threshold)}")
    print(f"log2FC>{fc_threshold}的基因数量: {sum(df_stats['log2FC'] > fc_threshold)}")
    print(f"log2FC<-{fc_threshold}的基因数量: {sum(df_stats['log2FC'] < -fc_threshold)}")
    print(f"上调条件: (padj<{pvalue_threshold} & log2FC>{fc_threshold})的基因数量: {sig_up}")
    print(f"下调条件: (padj<{pvalue_threshold} & log2FC<-{fc_threshold})的基因数量: {sig_down}")
    
    if 'gene_type' in df_stats.columns:
        type_counts = df_stats['gene_type'].value_counts()
        print("\n基因类型统计:")
        for gene_type, count in type_counts.items():
            print(f" - {gene_type}: {count}个基因")
    
    print(f"\n分析完成！发现{sig_up}个上调基因和{sig_down}个下调基因（p < {pvalue_threshold}, |log2FC| > {fc_threshold}）")
    
    return df_stats

if __name__ == "__main__":
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='分析基因表达数据并生成可视化结果')
    
    # 添加参数
    parser.add_argument('-i', '--input', type=str, required=True, 
                        help='输入文件路径(CSV或Excel格式)')
    parser.add_argument('-o', '--output', type=str, default="results", 
                        help='输出结果目录')
    parser.add_argument('-p', '--pvalue', type=float, default=0.05,
                        help='差异表达的p值阈值（默认：0.05）')
    parser.add_argument('-f', '--fc', type=float, default=1.0,
                        help='差异表达的log2倍数变化阈值（默认：1.0）')
    parser.add_argument('-t', '--top', type=int, default=50,
                        help='热图中显示的顶部基因数量（默认：50）')
    parser.add_argument('-g', '--genes', type=str, 
                        help='指定关注的基因列表，用逗号分隔')

    args = parser.parse_args()
    
    custom_genes = None
    if args.genes:
        custom_genes = [g.strip() for g in args.genes.split(',')]
        print(f"将关注以下指定基因: {', '.join(custom_genes)}")
    
    main(args.input, args.output, args.pvalue, args.fc, args.top, custom_genes) 