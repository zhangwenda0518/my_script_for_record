import os
import pandas as pd
import argparse
from collections import defaultdict

# 读取PFAM_Domain文件
def read_pfam_domains(pfam_domain_files):
    pfam_domains = defaultdict(list)
    
    for file in pfam_domain_files:
        # 读取PFAM_Domain文件
        pfam_domain = pd.read_csv(file, sep='\t')
        
        # 遍历每一行
        for _, row in pfam_domain.iterrows():
            gene_family = row['Gene Family']
            domains = [row[col] for col in row.index if col.startswith("PFAM Domain")]
            
            # 将每个文件中的 PFAM Domain 条件存储为独立的规则
            pfam_domains[gene_family].append({"domains": domains, "require_all": True})
    
    return pfam_domains

# 筛选modified.blast.pfam.tab文件
def filter_pfam_results(pfam_domains, pfam_results_folder):
    # 创建输出文件夹（如果不存在）
    output_folder = os.path.join(pfam_results_folder, "filtered_results")
    os.makedirs(output_folder, exist_ok=True)
    
    for gene_family, domain_rules in pfam_domains.items():
        pfam_result_file = os.path.join(pfam_results_folder, f"{gene_family}_modified.blast.pfam.tab")
        if not os.path.exists(pfam_result_file):
            print(f"Warning: File not found: {pfam_result_file}")
            continue
        
        try:
            # 读取pfam结果文件
            pfam_results = pd.read_csv(pfam_result_file, sep='\t')
            
            # 检查是否包含必要的列
            required_columns = ['GeneID', 'HMM_Accession', 'HMM_Name', 'E-value', 'Score', 'Coverage']
            if not all(col in pfam_results.columns for col in required_columns):
                print(f"Error: Required columns not found in {pfam_result_file}. Expected columns: {required_columns}")
                continue
            
            # 初始化筛选结果
            final_results = pd.DataFrame()
            non_compliant_results = pd.DataFrame()
            
            # 遍历每个 PFAM Domain 规则
            for rule in domain_rules:
                domains = rule["domains"]
                require_all = rule["require_all"]
                
                # 筛选符合条件的行
                if require_all:
                    # 需要同时满足所有结构域
                    filtered_results = pfam_results[pfam_results['HMM_Accession'].isin(domains)]
                    # 检查是否同时包含所有需要的结构域
                    gene_ids = filtered_results['GeneID'].unique()
                    for gene_id in gene_ids:
                        gene_results = filtered_results[filtered_results['GeneID'] == gene_id]
                        if set(domains).issubset(set(gene_results['HMM_Accession'])):
                            final_results = pd.concat([final_results, gene_results])
                        else:
                            non_compliant_results = pd.concat([non_compliant_results, gene_results])
                else:
                    # 只需要满足其中一个结构域
                    filtered_results = pfam_results[pfam_results['HMM_Accession'].isin(domains)]
                    final_results = pd.concat([final_results, filtered_results])
            
            # 统计筛选保存的基因数目和不符合要求的基因数目
            num_filtered_genes = final_results['GeneID'].nunique()
            num_non_compliant_genes = non_compliant_results['GeneID'].nunique()
            
            # 输出统计信息
            print(f"Gene family: {gene_family}")
            print(f"Number of genes saved: {num_filtered_genes}")
            print(f"Number of genes that do not meet the domain requirements: {num_non_compliant_genes}")
            
            # 输出筛选结果
            if not final_results.empty:
                output_file = os.path.join(output_folder, f"{gene_family}_filtered.pfam.tab")
                final_results.drop_duplicates().to_csv(output_file, sep='\t', index=False)
                print(f"Filtered results saved to: {output_file}")
            else:
                print(f"No results found for gene family: {gene_family}")
            
            # 输出不符合要求的基因信息
            if not non_compliant_results.empty:
                unfiltered_file = os.path.join(output_folder, f"{gene_family}_unfiltered.pfam.tab")
                non_compliant_results.drop_duplicates().to_csv(unfiltered_file, sep='\t', index=False)
                print(f"Unfiltered results saved to: {unfiltered_file}")
            
            print("-" * 50)  # 分隔线
        
        except Exception as e:
            print(f"Error processing {pfam_result_file}: {e}")

# 主函数
def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="Filter PFAM results based on PFAM domain files.")
    parser.add_argument("pfam_domain_files", nargs="+", help="Paths to the PFAM_Domain files.")
    parser.add_argument("pfam_results_folder", help="Path to the folder containing modified.blast.pfam.tab files.")
    args = parser.parse_args()

    # 读取PFAM_Domain文件
    pfam_domains = read_pfam_domains(args.pfam_domain_files)
    
    # 筛选pfam结果文件
    filter_pfam_results(pfam_domains, args.pfam_results_folder)

if __name__ == "__main__":
    main()
