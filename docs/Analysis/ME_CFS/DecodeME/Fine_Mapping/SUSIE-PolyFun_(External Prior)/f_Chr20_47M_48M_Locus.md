---
tags:
  - SuSiE
hide:
  - toc
---
# Chr20 47M-48.2M
I applied PolyFun[@weissbrod2020functionally] SUSIE[@wang2020simple] fine-mapping to the DecodeME[@genetics2025initial] GWAS-1 signal on Chromosome 20, using the same methodology I previously applied to the [chromosome 1 locus](a_Polyfun_Chr1_173M_174M_Locus.md).  



### Comparison of configurations
The UpSetPlots below illustrate respectively

- The overlap across the four SUSIE runs of all variants found in credible sets, and
- The overlap across the four SUSIE runs of the minimal set of variants required to achieve a total PIP of 50%.



{{
png_embed("docs/_figs/decode_me_polyfun_explainchr20_47000000_48200000_palindromes_keep_polyfun_upset_all_cs_variants.png",
alt="upset plot for chrom 20")
}}


{{
png_embed("docs/_figs/decode_me_polyfun_explainchr20_47000000_48200000_palindromes_keep_polyfun_upset_cs50_variants.png",
alt="50 PIP upset plot for chrom 17")
}}

Unlike the previous loci, but consistent with my [uniform-prior runs at the chromosome 20 locus](../SUSIE/f_Chr20_47M_48M_Locus.md), here the results produced by SUSIE depend on the chosen configuration: The $L=10$ and $L=2$ runs are similar, and both produce two credible sets, while the $L=1$ and strict $L=10$ runs are also similar, and produce a single credible set.


### Results ($L=10$)

The plot and table below show the results for the $L=10$ run, which is similar to the $L=2 run$



{{
susie_polyfun_explain_plot("docs/_figs/decode_me_polyfun_explainchr20_47000000_48200000_palindromes_keep_l10_explain_plot_svg.svg")
}}


{{
susie_polyfun_data_table(src="docs/_figs/decode_me_polyfun_explainchr20_47000000_48200000_palindromes_keep_l10_explain_detailed_table.parquet",
id="chr20_polyfun_susie_table_l10")
}}


### Results ($L=1$)



The plot and table results below show the results for the $L=1$ run, which is similar to the strict $L=10$ run.

{{
susie_polyfun_explain_plot("docs/_figs/decode_me_polyfun_explainchr20_47000000_48200000_palindromes_keep_l1_explain_plot_svg.svg")
}}


{{
susie_polyfun_data_table(src="docs/_figs/decode_me_polyfun_explainchr20_47000000_48200000_palindromes_keep_l1_explain_detailed_table.parquet",
id="chr20_polyfun_susie_table_l1")
}}
