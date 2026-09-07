---
tags:
- SuSiE
hide:
- toc
---
# Polyfun Chr1 173.5M-174.5M


## Methodology

To extend the results from my [earlier](../SUSIE/a_Chr1_173M_174M_Locus.md) SUSIE[@wang2020simple] fine-mapping of the DecodeME GWAS-1 signal[@genetics2025initial], I applied SUSIE again, but this time used a Bayesian prior derived from functional genomic annotations, instead of a uniform prior.


As a linkage disequilibrium reference, I used a [UK Biobank LD matrix hosted on AWS Open Data](https://registry.opendata.aws/ukbb-ld/).  Because this LD reference uses GRCh37 coordinates, I used GWASLab to liftover the DecodeME GWAS-1 summary statistics to GRCh37.

As before, to assess robustness and sensitivity to configuration, I ran SUSIE four times

- Once with $L=10$,
- Once with $L=2$,
- Once with $L=1$,
- Once with $L=10$ and strict variant filtering.

As before, in my SUSIE runs, I retained palindromic SNPs whose strand orientation GWASLAB was able to determine from allele frequencies in the Thousand Genomes Project, and discarded other palindromic SNPs.

### Prior Construction

For this analysis, I used the precomputed prior provided by the authors of PolyFun[@weissbrod2020functionally][^prior_note]. The PolyFun authors created this prior as follows:

1.  Select 15 [UK biobank](../../../../../Data_Sources/UKBB.md) traits with mutual [genetic correlations](../../../../../Bioinformatics_Concepts/Genetic_Correlation.md) less than 0.2.
2.  Run l2-penalized [stratified linkage disequilibrium score regression](../../../../../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md) on these traits to estimate allocation of [heritability](../../../../../Bioinformatics_Concepts/Heritability.md) enrichment of functional annotations[^annotation_note] for each trait.
3.  Average these heritability enrichments across the 15 traits to produce cross-trait heritability enrichment for each functional annotation.
4.  Use these heritability enrichments to define a Bayesian prior: the prior probability that a variant is causal is proportional to the cross-trait heritability enrichment implied by its functional annotations.
5. For robustness, modify the prior by limiting its dynamic range and binning the prior probabilities.

The resulting prior upweights genetic variants with functional annotations that are associated with high heritability across a range of diverse traits in the UK Biobank

## Results

### Comparison across runs

I begin by comparing the credible set variants across the $L=1$, $L=2$, $L=10$, and strict $L=10$ runs. The results are plotted in the UpsetPlot below:

{{
png_embed("docs/_figs/decode_me_polyfun_explainchr1_173500000_174500000_palindromes_keep_polyfun_upset_all_cs_variants.png",
alt="upset plot for chrom 1")
}}

Restricting to the minimal set of variants constituting a total PIP exceeding 50% produces the UpsetPlot:

{{png_embed("docs/_figs/decode_me_polyfun_explainchr1_173500000_174500000_palindromes_keep_polyfun_upset_cs50_variants.png",
alt="cs50 upset plot for chrom 1")
}}

These results show that at the chromosome 1 locus, SUSIE is insensitive to configuration: we get the same variant set regardless.


### Detailed Fine mapping results


The plot below illustrates the results of $L=10$ SUSIE fine mapping with and without the PolyFun prior


{{
susie_polyfun_explain_plot("docs/_figs/decode_me_polyfun_explainchr1_173500000_174500000_palindromes_keep_l10_explain_plot_svg.svg")
}}

The table below provides detailed information on $L=10$ SUSIE credible-set variants with and without the polyfun prior.

{{
susie_polyfun_data_table(src="docs/_figs/decode_me_polyfun_explainchr1_173500000_174500000_palindromes_keep_l10_explain_detailed_table.parquet",
id="chr1_polyfun_susie_table")
}}


These results show that while with a uniform prior SUSIE produces a rather diffuse credible set, with the PolyFun prior, the credible set is significantly more concentrated on a small number of variants.  The top PolyFun SUSIE variant is the SNP **chr1:173855298:A:T**.  PolyFun appears to have assigned a high prior weight to this variant due to annotation in the _coding_  and _conserved_ families.  The second PolyFun SUSIE variant is the insertion **chr1:173838788:T:TG**.  PolyFun appears to have assigned high prior weight to this variant due to annotations in the **coding** and **promoter_or_enhancer** families.


To investigate further, we can look in detail at the full set of annotations for the top variants.

{{susie_polyfun_variant_detail_table(src="docs/_figs/decode_me_polyfun_explainchr1_173500000_174500000_palindromes_keep_l10_explain_per_variant_annotation_table.parquet" ,id="chr1_polyfun_susie_characterization")}}

The above table reveals that:

-  **chr1:173855298:A:T** is variant in the 3' untranslated region (_UTR_3_UCSC_common_) of the gene [ZBTB37](https://www.genecards.org/card/ZBTB37) that is strongly evolutionarily conserved in primates and mammals (_Conserved_Primate_phastCons46way_common, Conserved_Mammal_phastCons46way_common)_.
- **chr1:173838788:T:TG** is a variant in an evolutionarily ancient promoter (_Ancient_Sequence_Age_Human_Promoter_common_) for the ZBTB37 gene.



[^prior_note]: That is, I used the first approach [listed on the PolyFun wiki](https://github.com/omerwe/polyfun/wiki/1.-Computing-prior-causal-probabilities-with-PolyFun).

[^annotation_note]: The functional annotations used here come from the baseline model first described in Finucane et al. 2015[@finucane2015partitioning] and which have been extended by the Broad institute since then.  The version of the baseline model used by Polyfun author included 187 functional annotations[@weissbrod2020functionally], which cover domains as diverse as evolutionarily conserved regions, qtls, [epigenetic marks](../../../../../Bioinformatics_Concepts/Epigenetics.md), non-synomous regions, promoters and enhancers, and more.