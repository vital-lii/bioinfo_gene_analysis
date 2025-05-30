#需建立.env文件，并设置email和api_key
from Bio import Entrez
import pandas as pd
import time
import os
import sys
import argparse
from datetime import datetime
import logging
from dotenv import load_dotenv
from utils.path_manager import PathManager

# 加载环境变量
load_dotenv()

def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='PubMed文献检索工具')
    parser.add_argument('-o', '--output', help='指定输出目录路径')
    return parser.parse_args()

class PubMedSearcher:
    """PubMed文献检索类"""
    
    def __init__(self, email=None, api_key=None, output_dir=None):
        """
        初始化PubMed检索器
        
        Args:
            email (str, optional): NCBI账户邮箱
            api_key (str, optional): NCBI API密钥
            output_dir (str, optional): 自定义输出目录
        """
        # 设置路径管理器
        self.path_manager = PathManager(output_dir)
        self.path_manager.ensure_dirs()
        
        # 设置日志
        self.setup_logging()
        
        # 设置Entrez (优先使用参数值，否则使用环境变量)
        self.setup_entrez(
            email or os.getenv('email'),
            api_key or os.getenv('api_key')
        )
        
    def setup_entrez(self, email, api_key):
        """配置Entrez"""
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
            
    def setup_logging(self):
        """配置日志"""
        log_dir = self.path_manager.get_path('logs')
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = os.path.join(log_dir, f"pubmed_search_{timestamp}.log")
        
        logging.basicConfig(
            filename=log_file,
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        
    def search_pubmed(self, term, retmax=100, sort="relevance"):               #在此修改篇数
        """
        执行PubMed检索
        
        Args:
            term (str): 检索词
            retmax (int): 最大返回结果数
            sort (str): 排序方式 ("relevance" or "date")
            
        Returns:
            list: PMIDs列表
        """
        try:
            logging.info(f"开始检索: {term}")
            handle = Entrez.esearch(
                db="pubmed",
                term=term,
                retmax=retmax,
                sort=sort
            )
            record = Entrez.read(handle)
            pmids = record["IdList"]
            logging.info(f"检索到 {len(pmids)} 篇文献")
            return pmids
        except Exception as e:
            logging.error(f"检索出错: {str(e)}")
            return []
            
    def fetch_article_details(self, pmid, max_retries=3, timeout=10):
        """
        获取单篇文献详细信息，添加重试机制
        
        Args:
            pmid (str): PubMed ID
            max_retries (int): 最大重试次数
            timeout (int): 超时时间（秒）
            
        Returns:
            dict: 文献信息字典
        """
        for attempt in range(max_retries):
            try:
                time.sleep(0.5)  # 避免API限制
                handle = Entrez.efetch(
                    db="pubmed", 
                    id=pmid, 
                    retmode="xml",
                    timeout=timeout
                )
                article_data = Entrez.read(handle)
                article = article_data["PubmedArticle"][0]["MedlineCitation"]["Article"]
                
                # 提取摘要
                abstract = article.get("Abstract", {}).get("AbstractText", ["No abstract"])[0]
                
                # 提取MeSH词
                mesh_terms = article_data["PubmedArticle"][0]["MedlineCitation"].get("MeshHeadingList", [])
                keywords = ", ".join([mesh.get("DescriptorName", "") for mesh in mesh_terms])
                
                # 提取年份
                pub_date = article["Journal"]["JournalIssue"]["PubDate"]
                year = pub_date.get("Year", "Unknown")
                
                return {
                    "PMID": pmid,
                    "Title": article["ArticleTitle"],
                    "Abstract": abstract,
                    "Keywords": keywords,
                    "Year": year,
                    "Journal": article["Journal"]["Title"]
                }
                
            except Exception as e:
                if attempt < max_retries - 1:
                    print(f"\n获取PMID:{pmid}失败，正在重试({attempt+2}/{max_retries})...")
                    time.sleep(2)  # 重试前等待时间
                    continue
                else:
                    logging.error(f"获取文献 {pmid} 详情时出错: {str(e)}")
                    return None
            
    def batch_fetch_articles(self, search_terms):
        """批量检索并保存文献信息"""
        output_dir = self.path_manager.get_path('output')
        print(f"文献将保存在: {output_dir}")
        
        # 计算总主题数
        total_topics = len(search_terms)
        current_topic = 0
        
        for category, term in search_terms.items():
            current_topic += 1
            print(f"\n[主题 {current_topic}/{total_topics}] 处理类别: {category}")
            print(f"检索词: {term}")
            
            # 检索文献
            pmids = self.search_pubmed(term)
            print(f"找到 {len(pmids)} 篇文献")
            
            if not pmids:
                print("未找到相关文献，跳过")
                continue
            
            # 计算预估时间
            estimated_time = len(pmids) * 0.5 / 60  # 转换为分钟
            print(f"正在获取文献详细信息... (预计需要 {estimated_time:.1f} 分钟)")
            
            # 获取详细信息
            articles = []
            start_time = time.time()
            
            for i, pmid in enumerate(pmids, 1):
                # 计算进度和剩余时间
                elapsed_time = time.time() - start_time
                progress = i / len(pmids)
                if i > 1:  # 至少处理一篇后才计算剩余时间
                    remaining_time = (elapsed_time / (i-1)) * (len(pmids) - i+1) / 60
                    print(f"进度: {progress*100:.1f}% | "
                          f"处理第 {i}/{len(pmids)} 篇文献 (PMID: {pmid}) | "
                          f"预计还需 {remaining_time:.1f} 分钟", end='\r')
                else:
                    print(f"进度: {progress*100:.1f}% | "
                          f"处理第 {i}/{len(pmids)} 篇文献 (PMID: {pmid})", end='\r')
                
                article_info = self.fetch_article_details(pmid)
                if article_info:
                    articles.append(article_info)
            
            # 保存为DataFrame
            if articles:
                df = pd.DataFrame(articles)
                output_file = os.path.join(output_dir, f"{category}_literature.csv")
                df.to_csv(output_file, index=False, encoding='utf-8')
                print(f"\n已保存 {len(articles)} 篇文献到: {output_file}")
                
                # 显示成功率
                success_rate = len(articles) / len(pmids) * 100
                print(f"文献获取成功率: {success_rate:.1f}%")
            else:
                print("\n未能获取任何文献的详细信息")
            
        print("\n所有主题处理完成！")
        print(f"结果保存在: {output_dir}")

def main():
    try:
        print("开始运行PubMed文献检索...")
        
        # 解析命令行参数
        args = parse_args()
        
        # 配置搜索词 可进行修改或添加
        search_terms = {
            "Autophagy_Treg_plus": (
                "(Autophagy OR 'autophagy') AND "
                "(Treg OR 'regulatory T cell') "
                "OR (HDAC3 OR 'histone deacetylase 3')"
            )
        }
        
        print(f"\n将检索以下主题的文献:")
        for category, term in search_terms.items():
            print(f"- {category}: {term}")
        
        # 初始化检索器
        print("\n初始化PubMed检索器...")
        searcher = PubMedSearcher(output_dir=args.output)
        
        # 执行批量检索
        print("\n开始批量检索文献...")
        searcher.batch_fetch_articles(search_terms)
        
    except Exception as e:
        print(f"程序运行出错: {str(e)}")
        import traceback
        print(traceback.format_exc())

if __name__ == "__main__":
    main()
