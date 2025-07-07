cd /public/home/zpxu/Bohou/CRISPR/Snakemake_CRISPR_Barcode_HiTom_Analysis
mkdir -p Share_Software/Snakemake_CRISPR_Barcode_HiTom_Analysis
cp -r config/ Share_Software/Snakemake_CRISPR_Barcode_HiTom_Analysis
cp -r schemas/ Share_Software/Snakemake_CRISPR_Barcode_HiTom_Analysis
cp -r workflow/ Share_Software/Snakemake_CRISPR_Barcode_HiTom_Analysis
mkdir -p Share_Software/Snakemake_CRISPR_Barcode_HiTom_Analysis/logs
mkdir -p Share_Software/Snakemake_CRISPR_Barcode_HiTom_Analysis/results
mkdir -p Share_Software/Snakemake_CRISPR_Barcode_HiTom_Analysis/Resources/Barcode
mkdir -p Share_Software/Snakemake_CRISPR_Barcode_HiTom_Analysis/Resources/HiTom
cp README.md Share_Software/Snakemake_CRISPR_Barcode_HiTom_Analysis
cp LICENSE Share_Software/Snakemake_CRISPR_Barcode_HiTom_Analysis

cd Share_Software/
tar -Ipigz -cf Snakemake_CRISPR_Barcode_HiTom_Analysis.tar.gz Snakemake_CRISPR_Barcode_HiTom_Analysis
mv Snakemake_CRISPR_Barcode_HiTom_Analysis.tar.gz ../

cd ../
rm -rf Share_Software/

echo "Success. Can be share file [/public/home/zpxu/Bohou/CRISPR/Snakemake_CRISPR_Barcode_HiTom_Analysis/Snakemake_CRISPR_Barcode_HiTom_Analysis.tar.gz]."
