#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
增强版基因互作网络分析 - 提高基因识别率
使用STRING API并添加一级互作邻居
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import requests
import time
import json

def parse_args():
    parser = argparse.ArgumentParser(description='Enhanced gene interaction network analysis')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file with gene_symb column')
    parser.add_argument('-o', '--output', required=True, help='Output prefix for files')
    parser.add_argument('-s', '--species', required=True, choices=['human', 'mouse', 'rat'], 
                        help='Species (human, mouse, or rat)')
    parser.add_argument('--score', type=int, default=400, 
                        help='STRING interaction score threshold (0-1000, default: 400)')
    parser.add_argument('--include_neighbors', action='store_true',
                        help='Include first-degree neighbor interactions')
    parser.add_argument('--neighborhood_score', type=int, default=700,
                        help='Score threshold for neighbor interactions (0-1000, default: 700)')
    return parser.parse_args()

# 物种名到NCBI分类ID的映射
SPECIES_MAP = {
    'human': 9606,
    'mouse': 10090,
    'rat': 10116
}

def get_string_interactions(gene_list, organism, min_score=400):
    """从STRING数据库获取蛋白质互作关系"""
    string_api_url = "https://string-db.org/api/json/network"
    
    params = {
        "identifiers": "\n".join(gene_list),  # 基因名称列表
        "species": organism,                  # 物种ID
        "caller_identity": "pypath_alt_tool", # 调用者标识
        "required_score": min_score,          # 最小互作分数(0-1000)
        "add_nodes": 0                       # 不添加额外节点
    }
    
    try:
        response = requests.post(string_api_url, data=params)
        data = response.json()
        return data
    except Exception as e:
        print(f"Error fetching STRING interactions: {e}")
        return []

def get_string_id_mapping(gene_list, organism):
    """获取基因符号到STRING ID的映射"""
    string_api_url = "https://string-db.org/api/json/get_string_ids"
    
    params = {
        "identifiers": "\n".join(gene_list),  # 基因名称列表
        "species": organism,                  # 物种ID
        "limit": 1,                          # 每个名称返回的匹配数量
        "caller_identity": "pypath_alt_tool", # 调用者标识
        "echo_query": 1                      # 在结果中包含查询标识符
    }
    
    try:
        response = requests.post(string_api_url, data=params)
        data = response.json()
        
        # 创建映射字典
        mapping = {}
        for item in data:
            if 'preferredName' in item and 'queryItem' in item:
                mapping[item['queryItem']] = item['preferredName']
        
        return mapping, data
    except Exception as e:
        print(f"Error fetching STRING ID mapping: {e}")
        return {}, []

def get_string_interactions_by_id(string_ids, organism, min_score=400, add_nodes=0):
    """使用STRING ID获取互作关系"""
    string_api_url = "https://string-db.org/api/json/network"
    
    # 将STRING ID转换为不带版本号的格式
    cleaned_ids = [string_id.split('.')[0] for string_id in string_ids if string_id]
    
    params = {
        "identifiers": "\n".join(cleaned_ids),  # STRING ID列表
        "species": organism,                    # 物种ID
        "caller_identity": "pypath_alt_tool",   # 调用者标识
        "required_score": min_score,            # 最小互作分数
        "add_nodes": add_nodes,                # 是否添加额外节点
        "network_type": "functional"           # 功能网络(包括间接互作)
    }
    
    try:
        response = requests.post(string_api_url, data=params)
        data = response.json()
        return data
    except Exception as e:
        print(f"Error fetching STRING interactions by ID: {e}")
        return []

def main():
    # 解析命令行参数
    args = parse_args()
    
    # 设置物种
    organism = SPECIES_MAP.get(args.species.lower())
    if not organism:
        print(f"Error: Unsupported species '{args.species}'. Use 'human', 'mouse', or 'rat'.")
        sys.exit(1)
    
    # 创建输出目录
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # 读取输入文件
    try:
        print(f"Reading input file: {args.input}")
        data = pd.read_csv(args.input)
        if 'gene_symb' not in data.columns:
            print("Error: Input file must contain 'gene_symb' column")
            sys.exit(1)
            
        gene_symbols = data['gene_symb'].tolist()
        # 过滤掉非字符串或空值
        gene_symbols = [g for g in gene_symbols if isinstance(g, str) and g.strip()]
        print(f"Loaded {len(gene_symbols)} valid genes from input file")
    except Exception as e:
        print(f"Error loading input file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # 第1步：获取基因符号到STRING ID的映射
    print("Step 1: Mapping gene symbols to STRING IDs...")
    gene_to_string, mapping_data = get_string_id_mapping(gene_symbols, organism)
    
    # 保存映射数据
    with open(f"{args.output}_string_mapping.json", 'w') as f:
        json.dump(mapping_data, f, indent=2)
    
    recognized_genes = [gene for gene in gene_symbols if gene in gene_to_string]
    print(f"Recognized {len(recognized_genes)}/{len(gene_symbols)} genes ({len(recognized_genes)/len(gene_symbols)*100:.2f}%)")
    
    # 保存识别和未识别的基因列表
    with open(f"{args.output}_recognized_genes.txt", 'w') as f:
        f.write("\n".join(recognized_genes))
    
    unrecognized_genes = [gene for gene in gene_symbols if gene not in gene_to_string]
    with open(f"{args.output}_unrecognized_genes.txt", 'w') as f:
        f.write("\n".join(unrecognized_genes))
    
    # 如果识别的基因太少，尝试其他方法
    if len(recognized_genes) < 0.2 * len(gene_symbols):
        print("Warning: Less than 20% of genes were recognized. Trying alternative methods...")
        # 尝试转换为大写
        print("Trying uppercase gene symbols...")
        upper_genes = [gene.upper() for gene in unrecognized_genes]
        upper_to_string, _ = get_string_id_mapping(upper_genes, organism)
        
        # 更新映射
        for i, gene in enumerate(unrecognized_genes):
            if upper_genes[i] in upper_to_string:
                gene_to_string[gene] = upper_to_string[upper_genes[i]]
                recognized_genes.append(gene)
        
        print(f"After uppercase conversion: Recognized {len(recognized_genes)}/{len(gene_symbols)} genes ({len(recognized_genes)/len(gene_symbols)*100:.2f}%)")
    
    # 第2步：使用STRING ID获取互作关系
    print(f"Step 2: Retrieving interactions with score threshold {args.score}...")
    string_ids = [gene_to_string[gene] for gene in recognized_genes if gene in gene_to_string]
    
    # 分批处理，STRING API有限制
    batch_size = 100
    all_interactions = []
    
    for i in range(0, len(string_ids), batch_size):
        batch = string_ids[i:i+batch_size]
        print(f"  Processing batch {i//batch_size + 1}/{(len(string_ids)-1)//batch_size + 1} ({len(batch)} genes)...")
        
        interactions = get_string_interactions_by_id(batch, organism, min_score=args.score, add_nodes=0)
        if interactions:
            all_interactions.extend(interactions)
            print(f"    Retrieved {len(interactions)} interactions")
        else:
            print("    No interactions found for this batch")
        
        if i + batch_size < len(string_ids):  # 不是最后一批
            print("    Waiting 1 second to avoid API rate limits...")
            time.sleep(1)
    
    # 保存原始互作数据
    with open(f"{args.output}_raw_interactions.json", 'w') as f:
        json.dump(all_interactions, f, indent=2)
    
    # 第3步（可选）：获取一级邻居互作
    if args.include_neighbors and len(recognized_genes) > 0:
        print(f"Step 3: Retrieving first-degree neighbor interactions with score threshold {args.neighborhood_score}...")
        
        neighbor_interactions = []
        
        # 选择顶部基因以减少API调用
        # 如果识别的基因太多，只使用前60个
        top_genes = recognized_genes[:min(60, len(recognized_genes))]
        top_string_ids = [gene_to_string[gene] for gene in top_genes if gene in gene_to_string]
        
        for i, string_id in enumerate(top_string_ids):
            print(f"  Retrieving neighbors for gene {i+1}/{len(top_string_ids)}: {string_id}...")
            
            # 为每个基因获取高可信度邻居(添加10个)
            interactions = get_string_interactions_by_id([string_id], organism, 
                                                       min_score=args.neighborhood_score, 
                                                       add_nodes=10)
            if interactions:
                neighbor_interactions.extend(interactions)
                print(f"    Retrieved {len(interactions)} neighbor interactions")
            else:
                print("    No neighbor interactions found")
            
            if i < len(top_string_ids) - 1:  # 不是最后一个基因
                print("    Waiting 1 second to avoid API rate limits...")
                time.sleep(1)
        
        # 将邻居互作添加到所有互作中
        all_interactions.extend(neighbor_interactions)
        
        # 保存包含邻居的原始互作数据
        with open(f"{args.output}_with_neighbors_raw.json", 'w') as f:
            json.dump(all_interactions, f, indent=2)
    
    # 第4步：构建网络
    print("Step 4: Building interaction network...")
    network = nx.Graph()
    edges_added = set()  # 跟踪已添加的边以避免重复
    
    # 创建从STRING ID到基因符号的逆映射
    string_to_gene = {v: k for k, v in gene_to_string.items()}
    
    for interaction in all_interactions:
        source_id = interaction.get('stringId_A', '')
        target_id = interaction.get('stringId_B', '')
        source_name = interaction.get('preferredName_A', source_id)
        target_name = interaction.get('preferredName_B', target_id)
        score = float(interaction.get('score', 0)) / 1000.0  # 转换为0-1
        
        # 检查这两个基因是否已经是边的一部分
        edge_key = tuple(sorted([source_name, target_name]))
        if edge_key in edges_added:
            continue
        
        # 将名称映射回原始基因符号(如果可能)
        source_gene = string_to_gene.get(source_name, source_name)
        target_gene = string_to_gene.get(target_name, target_name)
        
        # 添加节点和边
        network.add_node(source_gene, string_id=source_id, is_input=source_gene in gene_symbols)
        network.add_node(target_gene, string_id=target_id, is_input=target_gene in gene_symbols)
        network.add_edge(source_gene, target_gene, weight=score)
        
        edges_added.add(edge_key)
    
    # 计算有多少原始输入基因在网络中
    input_genes_in_network = [n for n, d in network.nodes(data=True) if d.get('is_input', False)]
    print(f"Included {len(input_genes_in_network)}/{len(gene_symbols)} input genes in the network ({len(input_genes_in_network)/len(gene_symbols)*100:.2f}%)")
    
    # 保存网络数据
    network_data = {
        'source': [],
        'target': [],
        'score': [],
        'source_in_input': [],
        'target_in_input': []
    }
    
    for u, v, data in network.edges(data=True):
        network_data['source'].append(u)
        network_data['target'].append(v)
        network_data['score'].append(data.get('weight', 0.0))
        network_data['source_in_input'].append(u in gene_symbols)
        network_data['target_in_input'].append(v in gene_symbols)
    
    edges_df = pd.DataFrame(network_data)
    edges_df.to_csv(f"{args.output}_edges.csv", index=False)
    
    # 节点数据
    node_data = {
        'gene': list(network.nodes()),
        'in_input': [n in gene_symbols for n in network.nodes()],
        'degree': [network.degree(n) for n in network.nodes()],
        'string_id': [network.nodes[n].get('string_id', '') for n in network.nodes()]
    }
    
    # 计算中心性指标
    if len(network.nodes()) > 0:
        print("Calculating centrality metrics...")
        try:
            betweenness = nx.betweenness_centrality(network)
            node_data['betweenness'] = [betweenness.get(n, 0) for n in network.nodes()]
        except Exception as e:
            print(f"Could not calculate betweenness centrality: {e}")
            node_data['betweenness'] = [0] * len(network.nodes())
            
        try:
            closeness = nx.closeness_centrality(network)
            node_data['closeness'] = [closeness.get(n, 0) for n in network.nodes()]
        except Exception as e:
            print(f"Could not calculate closeness centrality: {e}")
            node_data['closeness'] = [0] * len(network.nodes())
    
    nodes_df = pd.DataFrame(node_data)
    nodes_df.to_csv(f"{args.output}_nodes.csv", index=False)
    
    print(f"Network data saved to {args.output}_edges.csv and {args.output}_nodes.csv")
    
    # 生成可视化
    if network.number_of_nodes() > 0:
        print("Generating visualizations...")
        
        # 只保留原始输入列表中的基因
        if args.include_neighbors:
            input_only_network = network.copy()
            non_input_nodes = [n for n in network.nodes() if n not in gene_symbols]
            print(f"Creating input-only visualization (removing {len(non_input_nodes)} neighbor nodes)...")
            # 不直接删除，而是创建子图
            input_only_network = network.subgraph([n for n in network.nodes() if n in gene_symbols])
            print(f"Input-only network has {input_only_network.number_of_nodes()} nodes and {input_only_network.number_of_edges()} edges")
        
        # 1. 基本网络图
        plt.figure(figsize=(14, 12))
        
        # 根据节点类型设置颜色和大小
        node_color = []
        node_size = []
        
        for n in network.nodes():
            if n in gene_symbols:
                # 输入基因 - 使用度数设置大小
                node_size.append(min(network.degree(n) * 50 + 100, 2000))
                node_color.append('red')  # 输入基因为红色
            else:
                # 非输入基因(邻居) - 固定大小
                node_size.append(100)
                node_color.append('lightblue')  # 邻居基因为浅蓝色
        
        # 设置布局
        pos = nx.spring_layout(network, k=0.3, iterations=50)
        
        # 绘制网络
        nx.draw_networkx(
            network, pos,
            node_color=node_color,
            node_size=node_size,
            font_size=8 if network.number_of_nodes() < 50 else 0,
            edge_color='gray',   # 改为灰色提高可见度
            width=[min(data.get('weight', 0.4) * 1000, 2.0) for u, v, data in network.edges(data=True)],  # 大幅增加边的宽度
            alpha=0.8
        )
        
        plt.title(f'Gene Interaction Network ({args.species}) - {network.number_of_nodes()} nodes, {network.number_of_edges()} edges')
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(f"{args.output}_network.png", dpi=300)
        
        # 2. 顶级枢纽基因子网络图
        if network.number_of_nodes() > 5:
            # 找出顶级枢纽基因，只考虑输入基因
            degrees = {n: d for n, d in network.degree() if n in gene_symbols}
            top_count = min(15, len(degrees))
            top_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:top_count]
            top_genes = [gene for gene, _ in top_nodes]
            
            # 创建顶级枢纽基因的子网络
            hub_network = network.subgraph(top_genes)
            
            plt.figure(figsize=(12, 10))
            hub_pos = nx.spring_layout(hub_network, k=0.4, iterations=100)
            
            node_size = [hub_network.degree(n) * 100 + 300 for n in hub_network.nodes()]
            
            nx.draw_networkx(
                hub_network, hub_pos,
                node_color='red',
                node_size=node_size,
                font_size=9,
                font_weight='bold',
                edge_color='gray',
                width=[min(hub_network[u][v].get('weight', 0.4) * 1000, 2.0) for u, v in hub_network.edges()],  # 使用与主网络相同的放大比例
                alpha=0.9
            )
            
            plt.title(f'Top {top_count} Hub Genes Network ({args.species})')
            plt.axis('off')
            plt.tight_layout()
            plt.savefig(f"{args.output}_hub_network.png", dpi=300)
        
        # 3. 度分布图
        plt.figure(figsize=(10, 6))
        # 只包括输入基因的度分布
        input_degrees = [d for n, d in network.degree() if n in gene_symbols]
        plt.hist(input_degrees, bins=min(20, len(set(input_degrees))), alpha=0.7, color='red', label='Input genes')
        
        if args.include_neighbors:
            # 包括邻居基因的度分布
            neighbor_degrees = [d for n, d in network.degree() if n not in gene_symbols]
            if neighbor_degrees:
                plt.hist(neighbor_degrees, bins=min(20, len(set(neighbor_degrees))), 
                         alpha=0.5, color='blue', label='Neighbor genes')
                plt.legend()
        
        plt.title("Node Degree Distribution")
        plt.xlabel("Degree")
        plt.ylabel("Frequency")
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.savefig(f"{args.output}_degree_dist.png", dpi=300)
        
        # 输出重要统计信息到文本文件
        with open(f"{args.output}_stats.txt", 'w') as f:
            f.write(f"Network Statistics for {args.output}\n")
            f.write("=" * 50 + "\n")
            f.write(f"Species: {args.species} (NCBI Taxonomy ID: {organism})\n")
            f.write(f"Total genes in input: {len(gene_symbols)}\n")
            f.write(f"Genes recognized by STRING: {len(recognized_genes)} ({len(recognized_genes)/len(gene_symbols)*100:.2f}%)\n")
            f.write(f"Input genes in network: {len(input_genes_in_network)} ({len(input_genes_in_network)/len(gene_symbols)*100:.2f}%)\n")
            f.write(f"Total nodes in network: {network.number_of_nodes()}\n")
            f.write(f"Total edges in network: {network.number_of_edges()}\n")
            f.write(f"Interaction score threshold: {args.score}\n")
            
            if args.include_neighbors:
                f.write(f"Neighbor nodes included: {len(network.nodes()) - len(input_genes_in_network)}\n")
                f.write(f"Neighbor score threshold: {args.neighborhood_score}\n")
            
            f.write(f"Network density: {nx.density(network):.4f}\n")
            
            try:
                f.write(f"Average clustering coefficient: {nx.average_clustering(network):.4f}\n")
            except Exception as e:
                f.write(f"Average clustering coefficient: N/A (Error: {e})\n")
            
            f.write("\nTop 20 Hub Genes (from input list):\n")
            f.write("-" * 30 + "\n")
            for gene, degree in sorted([(n, d) for n, d in network.degree() if n in gene_symbols], 
                                     key=lambda x: x[1], reverse=True)[:20]:
                f.write(f"{gene}: {degree} connections\n")
        
        print(f"Network statistics saved to {args.output}_stats.txt")
    else:
        print("Warning: Network is empty. No genes with interactions were found.")
    
    print("Analysis completed successfully!")

if __name__ == "__main__":
    main()
