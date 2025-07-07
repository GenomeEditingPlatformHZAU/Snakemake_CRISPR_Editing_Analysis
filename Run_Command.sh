bsub -q interactive -Is bash
export PATH=/public/home/zpxu/miniconda2/condabin:$PATH
snakemake -s workflow/Snakefile -np
snakemake -s workflow/Snakefile --cores 10 --use-conda
snakemake --report report.html --report-stylesheet config/custom-stylesheet.css

## 结果打包分享
cd /public/home/zpxu/Bohou/CRISPR/Snakemake_CRISPR_Barcode_HiTom_Analysis/results/ZZN_Z1_Barcode_pair_Results
tar -Ipigz -cf Cas9.tar.gz Cas9
