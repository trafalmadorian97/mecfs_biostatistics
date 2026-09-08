---
tags:
- MAGMA
hide:
- toc
---
# MAGMA GTEx
I applied [MAGMA](../../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] to [DecodeME](../../../Data_Sources/DecodeME.md), partially reproducing analysis from the DecodeME preprint[@genetics2025initial].  


## MAGMA Gene Analysis

I applied MAGMA's SNP-wise-mean model to the summary statistics from DecodeME GWAS-1. 


In this step:

- Data from the 1000 genomes projects was downloaded from the [MAGMA website](https://cncr.nl/research/magma/) and used as a linkage disequilibrium reference.
- [Build 151 of dbSNP](https://hgw2.soe.ucsc.edu/cgi-bin/hgTables?hgsid=2912494930_cRufLdpdc1ynRc2sCM3g1WGAWAgH&hgta_doSchemaDb=hg19&hgta_doSchemaTable=snp151Flagged
  ) was used to assign RSIDs to SNPs.
- Magma's default proximity-based rules were used to assign SNPs to genes.

MAGMA produces a table of genes, effect sizes, and p values.  Filtering these genes via the Benjamini-Hochberg procedure[@benjamini1995controlling] at a false discovery rate of 0.01, and joining with a database of gene descriptions from [Ensembl Biomart](https://jun2026.archive.ensembl.org/info/data/biomart/index.html) produces the following table:

{{ markdown_table("docs/_figs/decode_me_gwas_1_magma_gene_table_markdown.mdx", title="MAGMA Genes") }}



It is also useful get a genome-wide view of the patterns of gene significance induced by the SNP-wise-mean model. The plot below serves this purpose:



{{ plotly_embed("docs/_figs/decode_me_magma_gene_plot.html", id="magma-decode-me-gene-manhattan-window") }}




The MAGMA SNP-wise-mean model is sensitive to the rules used for allocation of variants to genes.  The results above use a conservative assignment strategy: only coding variants are assigned to genes.  It is interesting to compare these with results produced by less-restrictive assumptions.  In the plot below, variants between 35 kilobases upstream and 10 kilobases downstream of a given gene are associated with that gene.



{{ plotly_embed("docs/_figs/decode_me_magma_gene_plot_with_window.html", id="magma-decode-me-gene-manhattan") }}


## MAGMA Gene Property Analysis

I next applied [MAGMA's](../../../Bioinformatics_Concepts/MAGMA_Overview.md) gene property analysis module to [DecodeME](../../../Data_Sources/DecodeME.md).  This step combined the conservative gene analysis results above with tissue-specific RNA expression values from the [GTEx project](../../../Data_Sources/GTEx_RNAseq_Data.md)[@gtex2020gtex].  The aim was to identify tissues enriched for genes associated with ME/CFS.  The results are plotted below:

{{ plotly_embed("docs/_figs/decode_me_gwas_1_build_37_magma_ensemble_specific_tissue_gene_covar_analysis_bar_plot/magma_gene_set_plot.html", id="decode-me-magma-gtex-tissue-bar-plot") }}
In this plot, the y-axis corresponds to $-\log_{10}(p)$ values, the x-axis corresponds to tissue type (only tissues with the smallest p values are shown), and bars are colored according to whether their p value meets Bonferroni-corrected significance threshold,



These results unambiguously point to the nervous system as a major site of ME/CFS gene activity.



## Reproducing
To reproduce this analysis, run the {{ api_link("DecodeME Initial Analysis Script", "mecfs_bio.analysis.decode_me_initial_analysis") }}.

## Follow-Up Questions
1. Do other approaches to identify significant tissues from GWAS-summary statistics produce concordant or discordant results?
2. How reliable is the GTEx-based MAGMA gene-set-analysis approach to identifying significant tissues?  In other words: for diseases with well-understood pathological processes, does it produce results consistent with these processes?