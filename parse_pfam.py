import pandas as pd
from collections import defaultdict
import argparse

def main(file_path):
    # 读取文件
    df = pd.read_csv(file_path, sep='\t')

    # 创建一个字典来存储基因家族和对应的PFAM结构域及其出现次数
    gene_family_pfam = defaultdict(lambda: defaultdict(int))
    gene_family_genes = defaultdict(set)  # 存储每个基因家族的基因

    # 遍历每一行，统计每个基因家族中PFAM结构域的出现次数
    for index, row in df.iterrows():
        gene_family = row['GeneFamliy']
        gene_id = row['GeneID']
        hmm_accession = row['HMM_Accession']
        
        gene_family_pfam[gene_family][hmm_accession] += 1
        gene_family_genes[gene_family].add(gene_id)

    # 输出每个基因家族的PFAM结构域比例
    for gene_family, pfam_counts in gene_family_pfam.items():
        total_genes = len(gene_family_genes[gene_family])  # 该基因家族中的基因总数
        print(f"Gene Family: {gene_family}")
        print(f"  Total Genes: {total_genes}")
        
        for pfam_domain, count in pfam_counts.items():
            proportion = count / total_genes  # 计算该PFAM结构域的比例
            print(f"    PFAM Domain: {pfam_domain}, Proportion: {proportion:.2f}")
        print("\n")

    # 将结果按PFAM结构域数量分类
    results_by_domain_count = defaultdict(list)

    for gene_family, pfam_domains in gene_family_pfam.items():
        pfam_list = sorted(list(pfam_domains.keys()))  # 对PFAM结构域进行排序
        domain_count = len(pfam_list)  # PFAM结构域的数量
        results_by_domain_count[domain_count].append([gene_family] + pfam_list)

    # 保存结果到不同的文件
    for domain_count, results in results_by_domain_count.items():
        # 创建列名
        columns = ['Gene Family'] + [f'PFAM Domain{i+1}' for i in range(domain_count)]
        
        # 将结果转换为DataFrame
        result_df = pd.DataFrame(results, columns=columns)
        
        # 保存到文件
        output_file = f'PFAM_Domain{domain_count}.tab'
        result_df.to_csv(output_file, sep='\t', index=False)
        print(f"Saved {len(results)} gene families with {domain_count} PFAM domains to {output_file}")

if __name__ == "__main__":
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description="Process a PFAM file and classify gene families by the number of PFAM domains.")
    parser.add_argument('file_path', type=str, help='Path to the input PFAM file (tab-separated).')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 调用主函数
    main(args.file_path)
