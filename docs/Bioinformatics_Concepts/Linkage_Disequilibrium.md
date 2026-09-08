# Linkage Disequilibrium
Linkage disequilibrium (LD) refers to statistical dependence between genetic variants. LD is central to statistical genomics[^handbook_note][^name_note].

## Measures

- When we are interested in patterns of LD across a genomic region with $n$ variants, it is common to report the LD matrix $R\in\mathbb{R}^{n\times n}$, whose $(i,j)$ component is $r_{i,j}$ the [Pearson correlation](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) between variant $i$ and variant $j$.
- Note, however, that this matrix $R$ reflects only pairwise dependence, and so is not a complete characterization of LD.  In particular, there are many possible higher-order dependence structures consistent with any given $R$ matrix[^corr_example].

## Drivers

There are two main drivers of LD: mutation and recombination.

### Mutation

For simplicity, first consider LD in the absence of recombination, as occurs in mitochondrial DNA and certain regions of the Y chromosome.  In such recombination-free regions, the distance between two variants is irrelevant to their LD. Instead, LD is a function of historical mutations, and the fates of populations containing these mutations.  

Figure 7 from the Hapmap paper[@international2005haplotype] illustrates the concept: variants in non-recombining regions tends to be highly correlated if they arose on the same branch of a genealogical tree.


![hapmap-mutation-fig](https://github.com/user-attachments/assets/d56eb383-5edf-4601-b4c5-a991acd25931)



### Recombination


Besides mutation, the other major driver of LD in the eukaryotic genome is recombination. In regions subject to recombination, LD decays as the distance between variants increases, because the odds of an intervening [recombination event](https://en.wikipedia.org/wiki/Genetic_recombination) correspondingly increase. However, due to the complex structure of eukaryotic DNA, the odds of recombination events are non-uniform across a chromosome. Thus, the rate of LD decay with genomic distance is not constant. Instead, LD displays a block-like structure, with block boundaries determined by recombination hotspots.


As an illustrative example, here is a plot of the absolute value of the correlation between genetic variants in a region of chromosome 1.  This plot was generated from the [UK Biobank LD matrices stored on AWS OpenData](https://registry.opendata.aws/ukbb-ld/).  In the plot, the x and y axes correspond genomic position, while color indicates absolute correlation.


![ld_example_plot](https://github.com/user-attachments/assets/a05681d5-91f3-4b89-8023-d3d50a22b8bd)

Consistent with the above, we observe irregularly spaced LD blocks. 

### Other Factors

Besides mutation and recombination, LD is affected broadly by historical patterns of migration, population isolation, and natural selection.


## Genomic Distance

It is frequently useful to measure distance along the genome not in number of base pairs, but in recombination frequency. For this purpose, the preferred unit is the [centimorgan](https://en.wikipedia.org/wiki/Centimorgan). Two genomic positions are one centimorgan apart if there is 1% chance of a recombination event between them per generation[^time_love].



## Genotyping

The strong linkage disequilibrium between variants in close proximity means that it is possible to capture a large proportion of common human genetic variation without sequencing all variants.  Instead, researchers can genotype a set of carefully chosen common variants called "tag variants".  Un-genotyped common variants can be imputed with high accuracy from their LD relationships with tag variants.



[^handbook_note]: For an overview of linkage disequilibrium with reference to population genomic models like the famous [coalescent](https://en.wikipedia.org/wiki/Coalescent_theory) see chapter 2 of the Handbook of Statistical Genomics[@balding2019handbook].

[^name_note]: As noted in chapter 2 of the Handbook of Statistical Genomics, the word "disequilibrium" in "linkage disequilibrium" is somewhat inaccurate, in the sense that statistical dependence between variants occurs in population-genetic equilibrium.

[^corr_example]: The following example is illustrative.  Consider 3 variants A, B and C. If all 3 variants are independent, the correlation matrix $R$ will equal the identity matrix.  If A and B are independent, but C is the [exclusive or](https://en.wikipedia.org/wiki/Exclusive_or) of A and B, their correlation matrix is also the identity matrix.

[^time_love]: For a readable popular-science account of the early history of genetics and the development of the centimorgan from the studies of fruit-fly mutants, see _Time, Love, and Memory_[@weiner2000time].