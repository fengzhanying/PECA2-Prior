import pandas as pd

def extract_tss_from_gtf(gtf_file, output_file):
    tss_data = []
    
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            
            gene_name = None
            for attr in attributes.split(";"):
                if "gene_name" in attr:
                    gene_name = attr.split('"')[1]
                    break
            gene_id = None
            for attr in attributes.split(";"):
                if "gene_id" in attr:
                    gene_id = attr.split('"')[1]
                    break
            gene_type = None
            for attr in attributes.split(";"):
                if "gene_biotype" in attr:
                    gene_type = attr.split('"')[1]
                    break
            if gene_type != "protein_coding":
                continue
            if not gene_name:
                continue
            if feature_type != "gene" and gene_name not in []:
                continue
            tss = start if strand == "+" else end
            tss_data.append([gene_id, gene_name, chrom, tss, strand])
    

    tss_df = pd.DataFrame(tss_data, columns=["GeneID", "GeneName", "Chromosome", "TSS",  "Strand"])
    tss_df.to_csv(output_file, sep="\t", index=False)
    print(f"TSS信息已保存到 {output_file}")

# 示例使用
gtf_file = "/lustre/home/zhangfy/fetalmonkey/Mul10-ensemble/Macaca_mulatta.Mmul_10.114.chr.gtf"
output_file = "Monkey_TSS.txt"
extract_tss_from_gtf(gtf_file, output_file)
