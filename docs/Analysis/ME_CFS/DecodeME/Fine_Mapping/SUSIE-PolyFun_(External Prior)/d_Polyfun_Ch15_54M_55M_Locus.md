---
tags:
  - SuSiE
---
# Polyfun Chr15 54M-55M

I applied PolyFun[@weissbrod2020functionally] SUSIE[@wang2020simple] fine-mapping to the DecodeME GWAS-1 signal[@genetics2025initial] on Chromosome 15, using the same methodology as I applied to the [chromosome 1 locus](a_Polyfun_Chr1_173M_174M_Locus.md).  


## Results

### Comparison of configurations

The first UpSetPlot below shows that all 4 SUSIE configurations found the same set of variants.

{{
png_embed("docs/_figs/decode_me_polyfun_explainchr15_54500000_55500000_palindromes_keep_polyfun_upset_all_cs_variants.png",
alt="upset plot for chrom 15")
}}

The second UpSetPlot compares the minimal set of variants needed to achieve 50% PIP.  Again, this set of variants is identical across all 4 SUSIE configurations.

{{
png_embed("docs/_figs/decode_me_polyfun_explainchr15_54500000_55500000_palindromes_keep_polyfun_upset_cs50_variants.png",
alt="50 PIP upset plot for chrom 15")
}}


### Detailed Fine mapping results

In contrast to the [earlier chromosome 1 run](a_Polyfun_Chr1_173M_174M_Locus.md) at the chromosome 15 locus, the uniform-prior and PolyFun-prior runs produce very similar results. The plot below provides an overview of the results


{{
susie_polyfun_explain_plot("docs/_figs/decode_me_polyfun_explainchr15_54500000_55500000_palindromes_keep_l10_explain_plot_svg.svg")
}}

The table below provides detailed information on $L=10$ SUSIE credible-set variants with and without the polyfun prior.  The _lift_ column shows that at this chromosome 15 locus, no genetic variant is significantly affected by the PolyFun prior.

{{
susie_polyfun_data_table(src="docs/_figs/decode_me_polyfun_explainchr15_54500000_55500000_palindromes_keep_l10_explain_detailed_table.parquet",
id="chr15_polyfun_susie_table")
}}


