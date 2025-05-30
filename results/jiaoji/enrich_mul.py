#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import os
import numpy as np

def read_input_file(file_path, column='gene_symb'):
    """
    读取输入文件，支持多种格式（CSV, TSV, TXT）
    
    参数:
    file_path: 输入文件路径
    column: 基因符号列名
    
    返回:
    pandas DataFrame
    """
    try:
        # 根据文件扩展名自动判断分隔符
        if file_path.endswith('.csv'):
            data = pd.read_csv(file_path)
        elif file_path.endswith('.tsv') or file_path.endswith('.txt'):
            data = pd.read_csv(file_path, sep='\t')
        else:
            # 尝试自动检测分隔符
            with open(file_path, 'r') as f:
                first_line = f.readline()
                if ',' in first_line:
                    data = pd.read_csv(file_path)
                else:
                    data = pd.read_csv(file_path, sep='\t')
        
        print(f"成功读取文件: {file_path}")
        print(f"文件包含 {len(data)} 行数据")
        
        # 检查列名
        if column not in data.columns:
            # 如果指定的列名不存在，但只有一列，直接使用该列
            if len(data.columns) == 1:
                data = data.rename(columns={data.columns[0]: column})
                print(f"使用唯一列作为基因符号列: {data.columns[0]}")
            else:
                available_columns = ', '.join(data.columns)
                raise ValueError(f"未找到列 '{column}'。可用的列: {available_columns}")
        
        return data
        
    except Exception as e:
        print(f"读取文件错误: {e}")
        raise

def main():
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='基因集富集分析脚本')
    
    # 添加基本参数
    parser.add_argument('-l', '--list_libraries', action='store_true', help='列出所有可用的基因集库')
    parser.add_argument('-s', '--species', default='Mouse', help='物种 (默认: Mouse)')
    parser.add_argument('-n', '--number', type=int, default=20, help='可视化展示的通路数量 (默认: 20)')
    parser.add_argument('-c', '--column', default='gene_symb', help='包含基因符号的列名 (默认: gene_symb)')
    parser.add_argument('-k', '--kegg_library', default='', help='指定基因集库名称 (如果不指定则使用默认)')
    parser.add_argument('-p', '--pvalue', type=float, default=0.05, help='原始P值显著性阈值 (默认: 0.05)')
    parser.add_argument('-d', '--debug', action='store_true', help='启用调试模式，显示更多信息')
    parser.add_argument('--use_original_pval', action='store_true', help='使用原始P值而非调整后P值进行筛选')
    
    # 只有当不是列出库的模式时，才需要输入和输出参数
    input_output_group = parser.add_argument_group('输入输出参数')
    input_output_group.add_argument('-i', '--input', help='输入文件路径 (支持 CSV, TSV, TXT 格式)')
    input_output_group.add_argument('-o', '--output', help='输出结果文件夹路径')
    
    args = parser.parse_args()
    
    # 列出所有可用的基因集库
    if args.list_libraries:
        print("可用的基因集库:")
        libraries = gp.get_library_name()
        for lib in libraries:
            print(f"- {lib}")
        return
    
    # 检查是否提供了必需的输入和输出参数
    if not args.input or not args.output:
        parser.error("当不使用 -l/--list_libraries 时，必须提供 -i/--input 和 -o/--output 参数")
    
    # 确保输出目录存在
    os.makedirs(args.output, exist_ok=True)
    
    # 读取输入文件
    try:
        data = read_input_file(args.input, args.column)
    except Exception as e:
        print(f"读取文件失败: {e}")
        return
    
    # 提取基因符号列表，确保所有值都是字符串
    if args.column in data.columns:
        # 将基因符号列转换为字符串类型，避免float对象没有strip方法的错误
        data[args.column] = data[args.column].astype(str)
        
        # 移除空格，过滤掉空值和NA值
        data[args.column] = data[args.column].str.strip()
        data = data[~data[args.column].isna() & (data[args.column] != '') & (data[args.column].str.lower() != 'nan')]
        
        gene_list = data[args.column].tolist()
        print(f"共提取到 {len(gene_list)} 个基因符号")
        
        # 调试模式下显示部分基因符号
        if args.debug:
            print("基因列表前10个样本:")
            print(gene_list[:10])
    else:
        print(f"错误: 输入文件中没有找到列 '{args.column}'")
        return
    
    # 先获取所有可用的库名称
    all_libraries = gp.get_library_name()
    print("可用的基因集库:")
    kegg_libraries = [lib for lib in all_libraries if 'KEGG' in lib.upper()]
    for lib in kegg_libraries:
        print(f"- {lib}")
    
    # 选择合适的基因集库
    if args.kegg_library:
        # 使用用户指定的库
        gene_sets = args.kegg_library
    else:
        # 自动选择合适的库
        species_key = args.species.lower()
        # 尝试查找匹配物种的KEGG库
        matched_libs = [lib for lib in kegg_libraries if species_key in lib.lower()]
        
        if matched_libs:
            gene_sets = matched_libs[0]  # 使用第一个匹配的库
        else:
            # 如果没有匹配物种的库，使用任何KEGG库
            gene_sets = kegg_libraries[0] if kegg_libraries else 'KEGG_2016'
    
    print(f"使用数据库: {gene_sets}")
    
    # 确定库类型（用于输出消息）
    if 'KEGG' in gene_sets:
        library_type = 'KEGG通路'
        english_library_type = 'KEGG Pathways'
    elif 'GO_Biological_Process' in gene_sets:
        library_type = 'GO生物学过程'
        english_library_type = 'GO Biological Processes'
    elif 'GO_Molecular_Function' in gene_sets:
        library_type = 'GO分子功能'
        english_library_type = 'GO Molecular Functions'
    elif 'GO_Cellular_Component' in gene_sets:
        library_type = 'GO细胞组分'
        english_library_type = 'GO Cellular Components'
    elif 'Reactome' in gene_sets:
        library_type = 'Reactome通路'
        english_library_type = 'Reactome Pathways'
    elif 'WikiPathway' in gene_sets:
        library_type = 'WikiPathway通路'
        english_library_type = 'WikiPathways'
    else:
        library_type = '基因集'
        english_library_type = 'Gene Sets'
    
    # 进行基因集富集分析
    try:
        print(f"正在进行{library_type}富集分析...原始P值阈值: {args.pvalue}")
        
        # 如果基因列表特别长，可能需要进行分批处理或采样
        if len(gene_list) > 3000 and args.debug:
            print(f"警告: 基因列表包含 {len(gene_list)} 个基因，这可能导致分析速度变慢或内存问题")
            
        # 确保所有基因符号都是有效的字符串
        valid_genes = [str(gene).strip() for gene in gene_list if gene and str(gene).strip().lower() != 'nan']
        if len(valid_genes) < len(gene_list):
            print(f"注意: 从原始列表中移除了 {len(gene_list) - len(valid_genes)} 个无效基因符号")
            gene_list = valid_genes
        
        # 为确保获取所有结果，先使用宽松的P值运行分析    
        enr = gp.enrichr(gene_list=gene_list, 
                         gene_sets=gene_sets,
                         organism=args.species,
                         outdir=args.output,
                         cutoff=1.0)  # 使用宽松P值以获取全部结果
        
        # 获取结果并按原始P值排序
        result = enr.results
        result = result.sort_values('P-value')  # 按原始P值排序
        
        # 根据用户选择的P值阈值筛选结果
        filtered_result = result[result['P-value'] <= args.pvalue]
        
        # 保存全部结果
        all_result_file = os.path.join(args.output, f'{library_type.replace(" ", "_")}_all_results.csv')
        result.to_csv(all_result_file, index=False)
        print(f"所有富集分析结果已保存至: {all_result_file}")
        
        # 保存筛选后的结果
        if not filtered_result.empty:
            filtered_file = os.path.join(args.output, f'{library_type.replace(" ", "_")}_filtered_pval_{args.pvalue}.csv')
            filtered_result.to_csv(filtered_file, index=False)
            print(f"按原始P值 {args.pvalue} 筛选后的结果已保存至: {filtered_file}")
        
        # 检查是否有结果
        if filtered_result.empty:
            print(f"警告: 在原始P值阈值 {args.pvalue} 下没有找到显著富集的{library_type}。")
            print(f"尝试提高P值阈值，例如使用 -p 0.1 或更高，或查看完整结果文件。")
            
            # 如果没有显著结果，但原始结果不为空，使用原始结果进行可视化
            if not result.empty:
                print(f"使用所有结果进行可视化...")
                # 使用前10个结果进行展示
                top_results = result.head(10)
                print(f"\n前10个{library_type}(无P值过滤):")
                print(top_results[['Term', 'Overlap', 'P-value', 'Odds Ratio']])
            else:
                print(f"没有找到任何{library_type}结果，可能是基因列表与数据库匹配问题")
                
                # 调试信息
                if args.debug:
                    debug_file = os.path.join(args.output, 'gene_list_debug.txt')
                    with open(debug_file, 'w') as f:
                        f.write('\n'.join(gene_list[:1000]))  # 保存前1000个基因作为示例
                    print(f"已将部分基因列表保存到: {debug_file} 以供调试")
                    
                return
        
        # 选择用于可视化的数据
        vis_data = filtered_result if not filtered_result.empty else result.head(min(args.number, len(result)))
        
        # 创建可视化图表 - 使用自定义函数，不依赖enr.plot
        print("生成可视化图表...")
        try:
            # 获取要可视化的结果数据
            n_terms = min(args.number, len(vis_data))
            plot_data = vis_data.head(n_terms)
            
            # 1. 创建条形图 - 不使用enr.plot
            plt.figure(figsize=(12, max(8, n_terms * 0.4)))  # 根据条目数量调整高度
            
            # 提取数据
            terms = plot_data['Term'].tolist()
            pvals = -np.log10(plot_data['P-value'].values)
            
            # 创建Y坐标位置
            y_pos = np.arange(len(terms))
            
            # 绘制横向条形图
            bars = plt.barh(y_pos, pvals, align='center', alpha=0.7)
            plt.yticks(y_pos, terms)
            plt.xlabel('-log10(P-value)')
            plt.title(f'Top {n_terms} Enriched {english_library_type} (P-value)')
            
            # 添加P值标签
            for i, bar in enumerate(bars):
                plt.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2, 
                         f'P={plot_data["P-value"].iloc[i]:.4f}', 
                         va='center')
            
            plt.tight_layout()
            
            # 保存条形图
            barplot_file = os.path.join(args.output, f'{library_type.replace(" ", "_")}_barplot.png')
            plt.savefig(barplot_file, dpi=300)
            plt.close()
            print(f"条形图已保存至: {barplot_file}")
            
            # 2. 创建点图
            try:
                # 提取更多数据 - 与上面相同但逆序
                plot_data_rev = plot_data.iloc[::-1]
                sizes = []
                labels = []
                
                for i, overlap in enumerate(plot_data_rev['Overlap']):
                    if isinstance(overlap, str) and '/' in overlap:
                        num, total = overlap.split('/')
                        num = int(num)
                        total = int(total)
                        sizes.append(num * 30)  # 点大小与基因数成比例
                        ratio = num / total * 100
                        labels.append(f"{num}/{total} ({ratio:.1f}%)")
                    else:
                        sizes.append(100)  # 默认大小
                        labels.append("N/A")
                
                # 创建散点图
                plt.figure(figsize=(12, max(8, n_terms * 0.4)))
                
                # 绘制点图
                scatter = plt.scatter(pvals, y_pos, s=sizes, alpha=0.7, c=pvals, cmap='viridis')
                
                plt.yticks(y_pos, terms)
                plt.xlabel('-log10(P-value)')
                plt.title(f'Top {n_terms} Enriched {english_library_type} - Dot Plot (P-value)')
                plt.grid(axis='x', linestyle='--', alpha=0.7)
                
                # 添加注释：基因重叠数
                for i, (x, y) in enumerate(zip(pvals, y_pos)):
                    plt.annotate(labels[i], 
                                 xy=(x, y), 
                                 xytext=(10, 0), 
                                 textcoords="offset points",
                                 va='center')
                
                plt.tight_layout()
                
                # 保存点图
                dotplot_file = os.path.join(args.output, f'{library_type.replace(" ", "_")}_dotplot.png')
                plt.savefig(dotplot_file, dpi=300)
                plt.close()
                print(f"点图已保存至: {dotplot_file}")
            except Exception as e:
                print(f"生成点图时出错: {e}")
                import traceback
                traceback.print_exc()
        except Exception as e:
            print(f"生成可视化图表时出错: {e}")
            import traceback
            traceback.print_exc()
        
        # 显示富集结果的基本信息
        if not filtered_result.empty:
            top_filtered = filtered_result.head(10)
            print(f"\n原始P值 <= {args.pvalue} 的通路富集结果:")
            print(top_filtered[['Term', 'Overlap', 'P-value', 'Odds Ratio']])
        else:
            top_original = result.head(10)
            print(f"\n前10个{library_type}富集结果(按原始P值排序):")
            print(top_original[['Term', 'Overlap', 'P-value', 'Odds Ratio']])
        
    except Exception as e:
        print(f"富集分析过程中出错: {e}")
        import traceback
        traceback.print_exc()  # 打印详细的错误堆栈，有助于调试
        return
    
    print("\n分析完成!")

if __name__ == "__main__":
    main()
