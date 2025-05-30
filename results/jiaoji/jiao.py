import pandas as pd
import argparse

# 创建参数解析器
parser = argparse.ArgumentParser(description='查找两列中重复的基因')
parser.add_argument('-i', '--input', required=True, help='输入文件路径')
parser.add_argument('-o', '--output', help='输出文件路径')

args = parser.parse_args()

# 尝试用 utf-8-sig 读取，自动去除BOM
df = pd.read_csv(args.input, sep=',', encoding='utf-8-sig')

# 打印实际列名，检查是否有隐藏字符
print(df.columns.tolist())

# 去除首尾空格后再试
df.columns = df.columns.str.strip()
print(df.columns.tolist())

# 获取两列的基因名
col1_genes = set(df['gene_symb'].dropna())
col2_genes = set(df['gene_symb_1'].dropna())

# 找出重复的基因
common_genes = col1_genes.intersection(col2_genes)

# 创建一个新的DataFrame来存储结果
result_df = df.copy()
result_df['common_genes'] = ''

# 将重复的基因添加到第三列
if common_genes:
    print(f"找到{len(common_genes)}个重复的基因：")
    print(sorted(list(common_genes)))
    
    # 在原始位置标记重复的基因
    for gene in common_genes:
        mask = (df['gene_symb'] == gene) | (df['gene_symb_1'] == gene)
        result_df.loc[mask, 'common_genes'] = gene

# 设置默认输出路径（如果未指定）
if not args.output:
    args.output = args.input.rsplit('.', 1)[0] + '_with_common.csv'

# 保存结果
result_df.to_csv(args.output, sep=',', index=False, encoding='utf-8')
print(f"\n结果已保存至：{args.output}")