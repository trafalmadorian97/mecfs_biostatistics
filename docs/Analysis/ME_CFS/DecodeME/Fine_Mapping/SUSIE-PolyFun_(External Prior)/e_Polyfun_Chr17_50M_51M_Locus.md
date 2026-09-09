---
tags:
  - SuSiE
---
# Polyfun Chr17 50M-51M

I applied PolyFun[@weissbrod2020functionally] SUSIE[@wang2020simple] fine-mapping to the DecodeME[@genetics2025initial] GWAS-1 signal on Chromosome 17, using the same methodology I previously applied to the [chromosome 1 locus](a_Polyfun_Chr1_173M_174M_Locus.md).  


### Comparison of configurations

The UpSetPlots below illustrate respectively 

- The overlap across the four SUSIE runs of all variants found in credible sets, and 
- The overlap across the four SUSIE runs of the minimal set of variants required to achieve a total PIP of 50%.


These plots show that in general, the variant sets selected by the four runs are very similar, though the $L=1$ run differs slightly from the others.



{{
png_embed("docs/_figs/decode_me_polyfun_explainchr17_50000000_51000000_palindromes_keep_polyfun_upset_all_cs_variants.png",
alt="upset plot for chrom 17")
}}


{{
png_embed("docs/_figs/decode_me_polyfun_explainchr17_50000000_51000000_palindromes_keep_polyfun_upset_cs50_variants.png",
alt="50 PIP upset plot for chrom 17")
}}


### Detailed Fine mapping results


The plot and table below illustrates the result of the $L=10$ PolyFun-prior SUSIE run.


{{
susie_polyfun_explain_plot("docs/_figs/decode_me_polyfun_explainchr17_50000000_51000000_palindromes_keep_l10_explain_plot_svg.svg")
}}


{{
susie_polyfun_data_table(src="docs/_figs/decode_me_polyfun_explainchr17_50000000_51000000_palindromes_keep_l10_explain_detailed_table.parquet",
id="chr17_polyfun_susie_table")
}}



While the uniform-prior SUSIE run produces a large and diffuse credible set, the PolyFun-prior SUSIE credible set is more concentrated, focusing especially on the two evolutionarily conserved variants **17:50291040:C:T**  and **17:50299079:G:A**.


The next table shows the functional annotations for these two top variants in full detail.


{{susie_polyfun_variant_detail_table(src="docs/_figs/decode_me_polyfun_explainchr17_50000000_51000000_palindromes_keep_l10_explain_per_variant_annotation_table.parquet" ,id="chr17_polyfun_susie_characterization")}}

Consistent with the above, both variants have high-weight evolutionarily conservation annotations (_Conserved_Primate_phastCons46way_common, Conserved_Mammal_phastCons46way_common, Conserved_LindbladToh_common_).  These annotations explain the greater PIP assigned to **17:50291040:C:T**  and **17:50299079:G:A** in the PolyFun-prior SUSIE runs.