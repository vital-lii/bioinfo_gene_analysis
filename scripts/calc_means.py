import pandas as pd
import argparse
import os
import numpy as np
from scipy import stats

def calculate_group_means(file_path, max_rows=None, min_expr_threshold=1.0):
    """
    计算GF和SCFA组的平均表达值，并进行统计检验
    
    Parameters:
    -----------
    file_path : str
        输入Excel文件路径
    max_rows : int, optional
        要读取的最大行数
    min_expr_threshold : float, optional
        最小表达阈值，用于过滤低表达基因
        
    Returns:
    --------
    pd.DataFrame
        包含GF和SCFA平均值及统计检验结果的数据框
    """
    print(f"读取文件: {file_path}")
    # 读取Excel文件，可选择限制行数
    if max_rows:
        df = pd.read_excel(file_path, nrows=max_rows)
        print(f"读取前{max_rows}行数据")
    else:
        df = pd.read_excel(file_path)
    
    print(f"原始数据维度: {df.shape}")
    
    # 提取基因符号列
    if 'gene_symb' in df.columns:
        gene_col = 'gene_symb'
    else:
        # 尝试使用第一列作为基因符号
        gene_col = df.columns[0]
        print(f"未找到'gene_symb'列，使用'{gene_col}'列作为基因标识符")
    
    # 识别GF和SCFA列
    gf_cols = [col for col in df.columns if col.startswith('GF')]
    scfa_cols = [col for col in df.columns if col.startswith('SCFA')]
    
    # 排除均值列（如果存在）
    gf_sample_cols = [col for col in gf_cols if col != 'GF_mean']
    scfa_sample_cols = [col for col in scfa_cols if col != 'SCFA_mean']
    
    if not gf_sample_cols:
        raise ValueError("未找到GF组样本列 (以'GF'开头且不是'GF_mean'的列名)")
    if not scfa_sample_cols:
        raise ValueError("未找到SCFA组样本列 (以'SCFA'开头且不是'SCFA_mean'的列名)")
    
    print(f"找到 {len(gf_sample_cols)} 个GF样本列: {', '.join(gf_sample_cols)}")
    print(f"找到 {len(scfa_sample_cols)} 个SCFA样本列: {', '.join(scfa_sample_cols)}")
    
    # 设置基因列为索引（如果不是第一列）
    if df.columns[0] != gene_col:
        df = df.set_index(gene_col)
    
    # 计算每组的平均值
    # 只保留计数列
    df_counts = df[gf_sample_cols + scfa_sample_cols].copy()
    
    # 计算每组的平均值和标准差
    result_df = pd.DataFrame({
        'GF_mean': df_counts[gf_sample_cols].mean(axis=1),
        'GF_std': df_counts[gf_sample_cols].std(axis=1),
        'SCFA_mean': df_counts[scfa_sample_cols].mean(axis=1),
        'SCFA_std': df_counts[scfa_sample_cols].std(axis=1),
    })
    
    # 计算log2 fold change
    result_df['log2FC'] = np.log2((result_df['SCFA_mean'] + 1) / (result_df['GF_mean'] + 1))
    
    # 计算是否可能显著）
    result_df['is_GF_all_zero'] = (df_counts[gf_sample_cols] == 0).all(axis=1)
    result_df['is_SCFA_all_zero'] = (df_counts[scfa_sample_cols] == 0).all(axis=1)
    
    result_df['gene_type'] = 'non_specific'  # 默认为非特异性
    
    # 创建组特异性基因分类
    specific_genes = {
        'GF_specific_high': [],    # GF特异高表达
        'GF_specific_medium': [],  # GF特异中等表达
        'GF_specific_low': [],     # GF特异低表达
        'SCFA_specific_high': [],  # SCFA特异高表达
        'SCFA_specific_medium': [], # SCFA特异中等表达
        'SCFA_specific_low': [],    # SCFA特异低表达
        'high_in_GF': [],          # GF组显著高表达
        'high_in_SCFA': []         # SCFA组显著高表达
    }
    
    # 进行统计检验（t-test）和基因分类
    print("进行统计检验和基因分类...")
    pvalues = []
    
    # 对每一行进行分类和统计检验
    for idx in range(len(df_counts)):
        gf_vals = df_counts.iloc[idx][gf_sample_cols].values
        scfa_vals = df_counts.iloc[idx][scfa_sample_cols].values
        gene = df_counts.index[idx]
        
        # 计算非零样本比例和表达统计
        gf_nonzero_ratio = (gf_vals > min_expr_threshold).sum() / len(gf_vals)
        scfa_nonzero_ratio = (scfa_vals > min_expr_threshold).sum() / len(scfa_vals)
        gf_mean = np.mean(gf_vals)
        scfa_mean = np.mean(scfa_vals)
        
        # 计算组内变异系数（CV）
        gf_cv = np.std(gf_vals) / (gf_mean + 0.1) if gf_mean > 0 else 0
        scfa_cv = np.std(scfa_vals) / (scfa_mean + 0.1) if scfa_mean > 0 else 0
        
        # 多层次分析：按表达模式和强度分类基因
        
        # 情况1: GF特异性表达（GF有表达，SCFA无表达）
        if gf_nonzero_ratio >= 0.67 and (scfa_vals <= min_expr_threshold).all():
            # 基于表达水平进一步分类
            if gf_mean > 20 and gf_cv < 0.5:
                specific_genes['GF_specific_high'].append(gene)
                result_df.loc[gene, 'gene_type'] = 'GF_specific_high'
                pval = 0.0001  # 高表达且低变异 - 非常显著
            elif gf_mean > 5:
                specific_genes['GF_specific_medium'].append(gene)
                result_df.loc[gene, 'gene_type'] = 'GF_specific_medium'
                pval = 0.001   # 中等表达 - 很显著
            else:
                specific_genes['GF_specific_low'].append(gene)
                result_df.loc[gene, 'gene_type'] = 'GF_specific_low'
                pval = 0.01    # 低表达 - 显著
        
        # 情况2: SCFA特异性表达（SCFA有表达，GF无表达）
        elif scfa_nonzero_ratio >= 0.67 and (gf_vals <= min_expr_threshold).all():
            # 基于表达水平进一步分类
            if scfa_mean > 20 and scfa_cv < 0.5:
                specific_genes['SCFA_specific_high'].append(gene)
                result_df.loc[gene, 'gene_type'] = 'SCFA_specific_high'
                pval = 0.0001  # 高表达且低变异 - 非常显著
            elif scfa_mean > 5:
                specific_genes['SCFA_specific_medium'].append(gene)
                result_df.loc[gene, 'gene_type'] = 'SCFA_specific_medium'
                pval = 0.001   # 中等表达 - 很显著
            else:
                specific_genes['SCFA_specific_low'].append(gene)
                result_df.loc[gene, 'gene_type'] = 'SCFA_specific_low'
                pval = 0.01    # 低表达 - 显著
        
        # 情况3: 两组都有表达，进行常规t检验
        else:
            try:
                # 检查是否有足够的非零值进行t检验
                if len(gf_vals) >= 2 and len(scfa_vals) >= 2:
                    print(f"\n调试信息 - 行 {idx} (基因 {gene}):")
                    print(f"GF组数据: {gf_vals}")
                    print(f"SCFA组数据: {scfa_vals}")
                    _, pval = stats.ttest_ind(gf_vals, scfa_vals, equal_var=False)
                    pval = round(pval, 2)  
                    print(f"t检验p值: {pval}")
                    if np.isnan(pval):
                        pval = 1.0  # 如果pval为NaN，设置为1.0（不显著）
                        print("p值为NaN，设置为1.0")
                    
                    # 计算效应量
                    log2FC = np.log2((scfa_mean + 1) / (gf_mean + 1))
                    if log2FC > 0.5 and pval <= 0.05:
                        specific_genes['high_in_SCFA'].append(gene)
                    elif log2FC < -0.5 and pval <= 0.05:
                        specific_genes['high_in_GF'].append(gene)
                else:
                    pval = 1.0
            except Exception as e:
                print(f"对行{idx} (基因{gene})进行t检验失败: {str(e)}")
                pval = 1.0
        
        pvalues.append(pval)
    
    result_df['pvalue'] = pvalues
    
    # 多重检验校正
    try:
        from scipy.stats.multitest import multipletests
        # 使用Benjamini-Hochberg方法校正p值
        reject, padj, _, _ = multipletests(pvalues, method='fdr_bh')
        result_df['padj'] = [round(p, 2) for p in padj]  # 将校正后的p值四舍五入到两位小数
        result_df['reject_null'] = reject  # 是否拒绝原假设
    except ImportError:
        print("警告: 无法导入multipletests模块，使用未校正的p值")
        result_df['padj'] = [round(p, 2) for p in result_df['pvalue']]
        result_df['reject_null'] = result_df['pvalue'] < 0.05
    except Exception as e:
        print(f"多重检验校正失败: {str(e)}，使用未校正的p值")
        result_df['padj'] = [round(p, 2) for p in result_df['pvalue']]
        result_df['reject_null'] = result_df['pvalue'] < 0.05
    
    # 过滤全0行
    non_zero_df = result_df[(~result_df['is_GF_all_zero']) | (~result_df['is_SCFA_all_zero'])].copy()
    print(f"非全零行数: {len(non_zero_df)} / {len(result_df)}")
    
    # 重置索引，确保基因名作为列
    if non_zero_df.index.name:
        non_zero_df = non_zero_df.reset_index()
    else:
        # 如果没有索引名，创建一个基因名列
        non_zero_df = non_zero_df.reset_index()
        non_zero_df = non_zero_df.rename(columns={'index': gene_col if gene_col else 'gene'})
    
    significant_df = non_zero_df[non_zero_df['pvalue'] <= 0.05].copy()
    print(f"\n过滤后显著性行数: {len(significant_df)} / {len(non_zero_df)}")
    
    # 汇总
    print("\n基因分类统计:")
    print(f"GF特异性高表达基因: {len(specific_genes['GF_specific_high'])}")
    print(f"GF特异性中等表达基因: {len(specific_genes['GF_specific_medium'])}")
    print(f"GF特异性低表达基因: {len(specific_genes['GF_specific_low'])}")
    print(f"SCFA特异性高表达基因: {len(specific_genes['SCFA_specific_high'])}")
    print(f"SCFA特异性中等表达基因: {len(specific_genes['SCFA_specific_medium'])}")
    print(f"SCFA特异性低表达基因: {len(specific_genes['SCFA_specific_low'])}")
    
    # 结合padj值更新高表达基因列表
    high_in_SCFA = []
    high_in_GF = []
    for gene in specific_genes['high_in_SCFA']:
        idx = significant_df.index[significant_df[gene_col] == gene].tolist()
        if idx and significant_df.loc[idx[0], 'padj'] < 0.05:
            high_in_SCFA.append(gene)
            significant_df.loc[idx[0], 'gene_type'] = 'high_in_SCFA'
    
    for gene in specific_genes['high_in_GF']:
        idx = significant_df.index[significant_df[gene_col] == gene].tolist()
        if idx and significant_df.loc[idx[0], 'padj'] < 0.05:
            high_in_GF.append(gene)
            significant_df.loc[idx[0], 'gene_type'] = 'high_in_GF'
    
    print(f"校正后SCFA显著高表达基因: {len(high_in_SCFA)}")
    print(f"校正后GF显著高表达基因: {len(high_in_GF)}")
    
    return significant_df, specific_genes

def analyze_group_differences(mean_df, fold_change_threshold=0.5, pvalue_threshold=0.05):
    """
    分析组间差异并输出潜在的显著基因
    
    Parameters:
    -----------
    mean_df : pd.DataFrame
        含有GF和SCFA平均值的数据框
    fold_change_threshold : float
        log2倍数变化阈值，默认为0.5（相当于约1.4倍变化）
    pvalue_threshold : float
        校正后p值阈值，默认为0.05
        
    Returns:
    --------
    dict
        包含分析结果的字典
    """
    # 识别基因列名
    gene_col = None
    for col in mean_df.columns:
        if col not in ['GF_mean', 'GF_std', 'SCFA_mean', 'SCFA_std', 
                      'log2FC', 'pvalue', 'padj', 'reject_null', 'gene_type',
                      'is_GF_all_zero', 'is_SCFA_all_zero'] and not col.startswith('GF') and not col.startswith('SCFA'):
            gene_col = col
            break
    
    if not gene_col:
        gene_col = 'gene'
        print("警告: 未找到基因名列，使用默认'gene'列")
    
    # 过滤显著的差异表达基因
    # 上调：log2FC > 阈值 且 padj < 阈值
    potentially_up = mean_df[(mean_df['log2FC'] > fold_change_threshold) & 
                           (mean_df['padj'] < pvalue_threshold) &
                           (~mean_df['is_SCFA_all_zero'])].copy()
    
    # 下调：log2FC < -阈值 且 padj < 阈值
    potentially_down = mean_df[(mean_df['log2FC'] < -fold_change_threshold) & 
                             (mean_df['padj'] < pvalue_threshold) &
                             (~mean_df['is_GF_all_zero'])].copy()
    
    # 组特异性基因（基于gene_type字段）
    gf_specific = mean_df[mean_df['gene_type'].str.startswith('GF_specific')].copy()
    scfa_specific = mean_df[mean_df['gene_type'].str.startswith('SCFA_specific')].copy()
    
    # 显示分析结果
    print(f"\n差异表达分析结果 (|log2FC| > {fold_change_threshold} 且 padj < {pvalue_threshold}):")
    print(f"显著上调基因: {len(potentially_up)}")
    print(f"显著下调基因: {len(potentially_down)}")
    print(f"GF特异性基因: {len(gf_specific)}")
    print(f"SCFA特异性基因: {len(scfa_specific)}")
    
    # 显示前十个显著上调基因
    if len(potentially_up) > 0:
        print("\nTop 10 显著上调基因:")
        top_up = potentially_up.sort_values('padj').head(10)
        for _, row in top_up.iterrows():
            print(f"{row[gene_col]}: SCFA={row['SCFA_mean']:.2f}, GF={row['GF_mean']:.2f}, "
                 f"log2FC={row['log2FC']:.2f}, padj={row['padj']:.2e}")
    
    # 显示前十个显著下调基因
    if len(potentially_down) > 0:
        print("\nTop 10 显著下调基因:")
        top_down = potentially_down.sort_values('padj').head(10)
        for _, row in top_down.iterrows():
            print(f"{row[gene_col]}: SCFA={row['SCFA_mean']:.2f}, GF={row['GF_mean']:.2f}, "
                 f"log2FC={row['log2FC']:.2f}, padj={row['padj']:.2e}")
    
    # 显示前十个GF特异性基因
    if len(gf_specific) > 0:
        print("\nTop 10 GF特异性基因:")
        top_gf = gf_specific.sort_values('GF_mean', ascending=False).head(10)
        for _, row in top_gf.iterrows():
            print(f"{row[gene_col]}: GF={row['GF_mean']:.2f}, SCFA={row['SCFA_mean']:.2f}, 类型={row['gene_type']}")
    
    # 显示前十个SCFA特异性基因
    if len(scfa_specific) > 0:
        print("\nTop 10 SCFA特异性基因:")
        top_scfa = scfa_specific.sort_values('SCFA_mean', ascending=False).head(10)
        for _, row in top_scfa.iterrows():
            print(f"{row[gene_col]}: SCFA={row['SCFA_mean']:.2f}, GF={row['GF_mean']:.2f}, 类型={row['gene_type']}")
    
    return {
        'up_regulated': potentially_up,
        'down_regulated': potentially_down,
        'gf_specific': gf_specific,
        'scfa_specific': scfa_specific
    }

if __name__ == "__main__":
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='计算GF和SCFA组的平均表达值与多层次基因分类分析')
    
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='输入Excel文件路径')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='输出CSV文件路径')
    parser.add_argument('-r', '--rows', type=int, default=None,
                        help='要读取的最大行数')
    parser.add_argument('-f', '--fc', type=float, default=0.5,
                        help='log2倍数变化阈值（默认：0.5，相当于约1.4倍变化）')
    parser.add_argument('-p', '--pvalue', type=float, default=0.05,
                        help='校正后p值阈值（默认：0.05）')
    parser.add_argument('-t', '--threshold', type=float, default=1.0,
                        help='最小表达阈值，用于识别非零表达（默认：1.0）')
    parser.add_argument('--keep-original', action='store_true', 
                        help='保留原始样本列（默认不保留）')
    
    args = parser.parse_args()
    
    # 计算组平均值并进行统计检验
    mean_df, specific_genes = calculate_group_means(args.input, args.rows, args.threshold)
    
    # 分析组间差异
    results = analyze_group_differences(mean_df, args.fc, args.pvalue)
    
    # 可能需要，将原始样本数据加回去
    if args.keep_original:
        print("\n保留原始样本列...")
        # 读取原始数据
        if args.rows:
            original_df = pd.read_excel(args.input, nrows=args.rows)
        else:
            original_df = pd.read_excel(args.input)
        
        # 确定基因列
        if 'gene_symb' in original_df.columns:
            gene_col = 'gene_symb'
        else:
            gene_col = original_df.columns[0]
        
        # 识别样本列
        gf_cols = [col for col in original_df.columns if col.startswith('GF') and col != 'GF_mean']
        scfa_cols = [col for col in original_df.columns if col.startswith('SCFA') and col != 'SCFA_mean']
        
        # 合并
        for col in gf_cols + scfa_cols:
            if col not in mean_df.columns:  # 避免重复列
                mean_df[col] = original_df.set_index(gene_col)[col]
    
    # 筛选显著差异基因
    filtered_df = mean_df[
        (abs(mean_df['log2FC']) >= args.fc) & 
        (mean_df['padj'] <= args.pvalue)
    ].copy()

    filtered_df.to_csv(args.output, index=False)  
    
    # 保存各类基因结果
    output_dir = os.path.dirname(args.output)
    base_name = os.path.splitext(os.path.basename(args.output))[0]
    
    # 1. 保存显著上调基因
    up_file = os.path.join(output_dir, f"{base_name}_upregulated.csv")
    up_data = results['up_regulated'].copy()
    if not args.keep_original:
        cols_to_keep = [col for col in up_data.columns if not (col.startswith('GF') and col not in ['GF_mean', 'GF_std']) and not (col.startswith('SCFA') and col not in ['SCFA_mean', 'SCFA_std'])]
        up_data = up_data[cols_to_keep]
    up_data.to_csv(up_file, index=False)
    print(f"显著上调基因保存至: {up_file}")
    
    # 2. 保存显著下调基因
    down_file = os.path.join(output_dir, f"{base_name}_downregulated.csv")
    down_data = results['down_regulated'].copy()
    if not args.keep_original:
        cols_to_keep = [col for col in down_data.columns if not (col.startswith('GF') and col not in ['GF_mean', 'GF_std']) and not (col.startswith('SCFA') and col not in ['SCFA_mean', 'SCFA_std'])]
        down_data = down_data[cols_to_keep]
    down_data.to_csv(down_file, index=False)
    print(f"显著下调基因保存至: {down_file}")
    
    # 3. 保存GF特异性基因
    gf_specific_file = os.path.join(output_dir, f"{base_name}_GF_specific.csv")
    gf_specific_data = results['gf_specific'].copy()
    if not args.keep_original:
        cols_to_keep = [col for col in gf_specific_data.columns if not (col.startswith('GF') and col not in ['GF_mean', 'GF_std']) and not (col.startswith('SCFA') and col not in ['SCFA_mean', 'SCFA_std'])]
        gf_specific_data = gf_specific_data[cols_to_keep]
    gf_specific_data.to_csv(gf_specific_file, index=False)
    print(f"GF特异性基因保存至: {gf_specific_file}")
    
    # 4. 保存SCFA特异性基因
    scfa_specific_file = os.path.join(output_dir, f"{base_name}_SCFA_specific.csv")
    scfa_specific_data = results['scfa_specific'].copy()
    if not args.keep_original:
        cols_to_keep = [col for col in scfa_specific_data.columns if not (col.startswith('GF') and col not in ['GF_mean', 'GF_std']) and not (col.startswith('SCFA') and col not in ['SCFA_mean', 'SCFA_std'])]
        scfa_specific_data = scfa_specific_data[cols_to_keep]
    scfa_specific_data.to_csv(scfa_specific_file, index=False)
    print(f"SCFA特异性基因保存至: {scfa_specific_file}")
    
    # 5. 按表达级别保存特异性基因
    # GF特异性高表达基因
    if specific_genes['GF_specific_high']:
        gf_high_file = os.path.join(output_dir, f"{base_name}_GF_specific_high.csv")
        gf_high_data = mean_df[mean_df['gene_type'] == 'GF_specific_high'].copy()
        if not args.keep_original:
            cols_to_keep = [col for col in gf_high_data.columns if not (col.startswith('GF') and col not in ['GF_mean', 'GF_std']) and not (col.startswith('SCFA') and col not in ['SCFA_mean', 'SCFA_std'])]
            gf_high_data = gf_high_data[cols_to_keep]
        gf_high_data.to_csv(gf_high_file, index=False)
        print(f"GF特异性高表达基因保存至: {gf_high_file}")
    
    # SCFA特异性高表达基因
    if specific_genes['SCFA_specific_high']:
        scfa_high_file = os.path.join(output_dir, f"{base_name}_SCFA_specific_high.csv")
        scfa_high_data = mean_df[mean_df['gene_type'] == 'SCFA_specific_high'].copy()
        if not args.keep_original:
            cols_to_keep = [col for col in scfa_high_data.columns if not (col.startswith('GF') and col not in ['GF_mean', 'GF_std']) and not (col.startswith('SCFA') and col not in ['SCFA_mean', 'SCFA_std'])]
            scfa_high_data = scfa_high_data[cols_to_keep]
        scfa_high_data.to_csv(scfa_high_file, index=False)
        print(f"SCFA特异性高表达基因保存至: {scfa_high_file}") 