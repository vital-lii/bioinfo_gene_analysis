import requests
import pandas as pd
import time
import xml.etree.ElementTree as ET
from urllib.parse import quote
import sys
import os
import re

# 配置信息
email = "vitalii@163.com"
api_key = "aa6c3591c99ff84bfd76c068a0984d38a208"

def read_p_values_csv(file_path):
    """读取包含p值的CSV文件"""
    try:
        df = pd.read_csv(file_path)
        # 检查文件是否包含必要的列
        if 'gene_id' not in df.columns or 'p_value' not in df.columns:
            print(f"警告: 文件缺少必要的列。当前列名: {df.columns.tolist()}")
            print("继续处理，但可能会丢失一些信息...")
        
        # 处理多个基因ID的情况（通过分号分隔）
        df_expanded = df.copy()
        
        # 对gene_id列进行处理，处理多个ID的情况
        if 'gene_id' in df.columns:
            # 将nan值和'-'转换为空字符串
            df_expanded['gene_id'] = df_expanded['gene_id'].fillna('').astype(str)
            df_expanded['gene_id'] = df_expanded['gene_id'].replace('-', '')
            
            # 提取所有非空的基因ID
            valid_gene_ids = []
            for idx, row in df_expanded.iterrows():
                gene_id = row['gene_id']
                if gene_id and gene_id != 'nan':
                    # 处理多个基因ID的情况（通过分号分隔）
                    if ';' in gene_id:
                        ids = [id.strip() for id in gene_id.split(';')]
                        valid_gene_ids.extend(ids)
                    else:
                        valid_gene_ids.append(gene_id)
            
            # 移除重复的基因ID
            valid_gene_ids = list(set([id for id in valid_gene_ids if id and id != 'nan']))
            
            print(f"从CSV文件中提取了 {len(valid_gene_ids)} 个唯一基因ID")
            return df_expanded, valid_gene_ids
        else:
            print("文件中没有找到gene_id列")
            return df_expanded, []
            
    except Exception as e:
        print(f"读取文件出错: {e}")
        sys.exit(1)

def fetch_gene_symbols_batch(gene_ids, start_idx, batch_size):
    """批量获取基因信息"""
    # 正确的NCBI E-utilities基本URL
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    batch = gene_ids[start_idx:start_idx+batch_size]
    results = []
    
    try:
        # 直接使用ESummary查询
        id_str = ",".join(map(str, batch))
        summary_url = f"{base_url}esummary.fcgi?db=gene&id={id_str}&api_key={api_key}&email={quote(email)}&retmode=xml"
        
        response = requests.get(summary_url, timeout=30)
        response.raise_for_status()
        
        # 解析结果XML
        root = ET.fromstring(response.content)
        
        for doc in root.findall(".//DocumentSummary"):
            entrez_id = doc.get("uid")
            
            # 提取基因符号和全名
            symbol_elem = doc.find(".//Name")
            description_elem = doc.find(".//Description")
            
            symbol = symbol_elem.text if symbol_elem is not None and hasattr(symbol_elem, 'text') else "Not Found"
            description = description_elem.text if description_elem is not None and hasattr(description_elem, 'text') else "Not Found"
            
            # 提取物种信息
            organism = "Not Found"
            organism_elem = doc.find(".//Organism")
            if organism_elem is not None:
                if hasattr(organism_elem, 'text'):
                    organism = organism_elem.text
            
            results.append({
                "gene_id": entrez_id,
                "gene_symbol": symbol,
                "description": description,
                "organism": organism
            })
        
        return results, None
    
    except requests.exceptions.RequestException as e:
        return [], f"请求错误: {e}"
    except ET.ParseError as e:
        return [], f"XML解析错误: {e}"
    except Exception as e:
        return [], f"未知错误: {e}"

def fetch_gene_symbols(gene_ids, batch_size=50, max_retries=3, retry_delay=5):
    """获取基因符号主函数，包含重试逻辑"""
    all_results = []
    
    for i in range(0, len(gene_ids), batch_size):
        retry_count = 0
        success = False
        
        while not success and retry_count < max_retries:
            print(f"处理基因ID {i+1}-{min(i+batch_size, len(gene_ids))} (共{len(gene_ids)}个)")
            
            results, error = fetch_gene_symbols_batch(gene_ids, i, batch_size)
            
            if error is None:
                all_results.extend(results)
                success = True
                # 遵守NCBI API使用规则：每秒不超过3个请求
                time.sleep(0.34)
            else:
                retry_count += 1
                print(f"错误: {error}")
                print(f"重试 ({retry_count}/{max_retries})...")
                time.sleep(retry_delay)
        
        if not success:
            print(f"警告: 跳过基因ID {i+1}-{min(i+batch_size, len(gene_ids))}，达到最大重试次数")
    
    return pd.DataFrame(all_results)

def fetch_single_gene(gene_id):
    """单个基因查询备选方案"""
    # 正确的NCBI E-utilities基本URL
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    try:
        # 使用ESummary查询单个基因
        summary_url = f"{base_url}esummary.fcgi?db=gene&id={gene_id}&api_key={api_key}&email={quote(email)}&retmode=xml"
        
        response = requests.get(summary_url, timeout=30)
        response.raise_for_status()
        
        # 解析结果XML
        root = ET.fromstring(response.content)
        
        doc = root.find(".//DocumentSummary")
        if doc is None:
            return {
                "gene_id": gene_id,
                "gene_symbol": "Not Found",
                "description": "Not Found",
                "organism": "Not Found"
            }
            
        # 提取基因符号和全名
        symbol_elem = doc.find(".//Name")
        description_elem = doc.find(".//Description")
        
        symbol = symbol_elem.text if symbol_elem is not None and hasattr(symbol_elem, 'text') else "Not Found"
        description = description_elem.text if description_elem is not None and hasattr(description_elem, 'text') else "Not Found"
        
        # 提取物种信息
        organism = "Not Found"
        organism_elem = doc.find(".//Organism")
        if organism_elem is not None and hasattr(organism_elem, 'text'):
            organism = organism_elem.text
        
        return {
            "gene_id": gene_id,
            "gene_symbol": symbol,
            "description": description,
            "organism": organism
        }
    
    except Exception as e:
        print(f"查询基因 {gene_id} 出错: {e}")
        return {
            "gene_id": gene_id,
            "gene_symbol": "Error",
            "description": "Error",
            "organism": "Error"
        }

def backup_method(gene_ids):
    """备用方法：尝试使用mygene包"""
    print("尝试使用mygene包查询基因...")
    try:
        import mygene
        mg = mygene.MyGeneInfo()
        results = mg.querymany(gene_ids, scopes='entrezgene', 
                             fields=['symbol', 'name', 'taxid', 'organism'],
                             species='all', returnall=True)
        
        output = []
        for hit in results['out']:
            gene_id = hit.get('query', '')
            gene_symbol = hit.get('symbol', 'Not Found')
            description = hit.get('name', 'Not Found')
            organism = hit.get('organism', 'Not Found')
            
            output.append({
                'gene_id': gene_id,
                'gene_symbol': gene_symbol,
                'description': description,
                'organism': organism
            })
        
        return pd.DataFrame(output)
    except ImportError:
        print("mygene包未安装，使用NCBI API逐个查询...")
        return backup_method_single(gene_ids)

def backup_method_single(gene_ids):
    """逐个查询基因的备用方法"""
    print("使用备用方法逐个查询基因...")
    results = []
    
    for i, gene_id in enumerate(gene_ids):
        print(f"处理基因 {i+1}/{len(gene_ids)}: {gene_id}")
        result = fetch_single_gene(gene_id)
        results.append(result)
        time.sleep(0.34)  # 遵守NCBI API使用规则
    
    return pd.DataFrame(results)

def merge_gene_info_with_original_data(original_df, gene_info_df):
    """将基因信息合并到原始数据中"""
    print("将基因符号信息合并到原始数据中...")
    
    # 如果没有gene_id列，直接返回原始数据
    if 'gene_id' not in original_df.columns:
        print("原始数据中没有gene_id列，无法合并信息")
        return original_df
    
    # 创建基因ID到符号的映射字典
    gene_symbol_dict = dict(zip(gene_info_df['gene_id'], gene_info_df['gene_symbol']))
    
    # 创建一个存储基因符号的新列
    original_df['mapped_gene_symbol'] = ''
    
    # 遍历原始数据
    for idx, row in original_df.iterrows():
        gene_id = str(row['gene_id'])
        if not gene_id or gene_id == 'nan' or gene_id == '-':
            continue
            
        # 处理多个基因ID的情况
        if ';' in gene_id:
            symbols = []
            for gid in gene_id.split(';'):
                gid = gid.strip()
                if gid in gene_symbol_dict:
                    symbols.append(gene_symbol_dict[gid])
            original_df.at[idx, 'mapped_gene_symbol'] = '; '.join(symbols) if symbols else ''
        else:
            if gene_id in gene_symbol_dict:
                original_df.at[idx, 'mapped_gene_symbol'] = gene_symbol_dict[gene_id]
    
    return original_df

def main():
    # 文件路径
    input_csv = "acu_asth_vs_asth_with_geneid_p_less_than_0.05.csv"  # 默认使用相对路径
    output_file = "enriched_gene_symbols.csv"  # 输出文件名
    
    # 命令行参数处理
    if len(sys.argv) > 1:
        input_csv = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 读取CSV文件
    print(f"读取CSV文件: {input_csv}")
    try:
        # 尝试读取文件
        original_df, gene_ids = read_p_values_csv(input_csv)
    except FileNotFoundError:
        print(f"文件未找到: {input_csv}")
        absolute_path = os.path.abspath(input_csv)
        print(f"尝试绝对路径: {absolute_path}")
        try:
            original_df, gene_ids = read_p_values_csv(absolute_path)
        except FileNotFoundError:
            print(f"请输入文件的完整路径:")
            input_path = input().strip('"\'')
            original_df, gene_ids = read_p_values_csv(input_path)
    
    # 如果没有有效的基因ID，直接返回原始数据
    if not gene_ids:
        print("没有找到有效的基因ID，将直接保存原始数据")
        original_df.to_csv(output_file, index=False)
        print(f"原始数据已保存到 {output_file}")
        return
    
    print(f"共提取 {len(gene_ids)} 个唯一基因ID")
    
    # 测试NCBI连接
    print("测试NCBI API连接...")
    try:
        # 使用正确的NCBI E-utilities基本URL
        test_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?api_key={api_key}"
        response = requests.get(test_url, timeout=10)
        response.raise_for_status()
        print("NCBI API连接正常")
    except Exception as e:
        print(f"NCBI API连接测试失败: {e}")
        print("请检查您的API密钥和互联网连接")
        print("将尝试使用备用方法...")
    
    # 获取基因符号
    print("从NCBI获取基因符号...")
    gene_info_df = None
    
    try:
        # 尝试安装mygene包
        try:
            import mygene
            print("使用mygene包查询...")
            gene_info_df = backup_method(gene_ids)
        except ImportError:
            print("mygene包未安装，尝试使用NCBI API...")
            gene_info_df = fetch_gene_symbols(gene_ids)
            
            # 如果结果为空，尝试备用方法
            if len(gene_info_df) == 0:
                print("批量查询未返回结果，尝试单个查询方法...")
                gene_info_df = backup_method_single(gene_ids)
    except Exception as e:
        print(f"错误: {e}")
        print("尝试备用方法...")
        gene_info_df = backup_method_single(gene_ids)
    
    # 合并基因信息到原始数据
    if gene_info_df is not None and not gene_info_df.empty:
        enriched_df = merge_gene_info_with_original_data(original_df, gene_info_df)
        
        # 保存结果
        print(f"保存结果到 {output_file}...")
        enriched_df.to_csv(output_file, index=False)
        
        print("完成！")
        print(f"共映射 {len(gene_info_df)} 个基因ID")
        print(f"结果已保存到 {output_file}")
    else:
        print("未能获取基因符号信息，将保存原始数据")
        original_df.to_csv(output_file, index=False)
        print(f"原始数据已保存到 {output_file}")

if __name__ == "__main__":
    main()
