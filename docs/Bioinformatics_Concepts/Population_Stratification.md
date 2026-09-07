# Population Stratification

## Causal privilege

Suppose we run an epidemiological study to understand the effect of variable A on variable B.  We detect an association between A and B. There are three main possibilities:


- **Causality**: A causes B.


``` mermaid
graph LR
A(A) --> B(B);

classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;
class A,B normal;
```


- **Reverse causality**: B causes A.


``` mermaid
graph LR
B(B) --> A(A);

classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;
class A,B normal;
```


- **Confounding**: A and B are both caused by a third variable C.


``` mermaid
graph LR
C(C) --> A(A);
C --> B(B);

classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;
class A,B,C normal;
```


Determining how much of the association to attribute to each of these possibilities is a challenge that frustrates much traditional epidemiological research[@hernan2010causal].


In contrast to traditional epidemiological research, genetic association studies are causally privileged.  Let A be a person's genotype and B be a phenotype of interest observed after birth.  In general:

- A person's genotype is fixed at conception, so reverse causality can be ruled out.
- Most kinds of environmental effects do not affect a person's genotype, so environmental confounding can be ruled out.


Thus a genotype-phenotype association is much more likely to be causal than a general epidemiological association.

This causal privilege is a significant advantage, but it does not mean that genetic studies are immune to causal inference complications.  One such complication is population stratification[@dattani2022clarifying].


## Types of stratification


### Genetic population stratification

Genetic population stratification occurs when the population under study contains multiple subpopulations, and mating within subpopulations has historically been more common than mating across subpopulations.  Normally, [linkage disequilibrium](Linkage_Disequilibrium.md) in humans decays to zero at a distance of a few megabases, and does not cross chromosomal boundaries.  Genetic population stratification changes this.  For example, SNP P on chromosome 1 and SNP Q on chromosome 2 may both be more common in a subpopulation than in the general population due to historical non-random mating.  Thus having  P increases your odds of being a member of the subpopulation, which increases your odds of having Q. P and Q are therefore correlated, despite being on different chromosomes. The is illustrated in the causal diagram below.



``` mermaid
graph LR
C(Subpopulation) --> A(SNP P);
C --> B(SNP Q);


classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;
class A,B,C normal;
```



Suppose now that P has a true causal effect on the phenotype of interest but Q does not.  The long-range correlation between P and Q will produce a GWAS association of Q with the phenotype, creating the false impression of causal GWAS hit in the vicinity of Q.  See below.


``` mermaid
graph LR
C(Subpopulation) --> A(SNP P);
C --> B(SNP Q);
A --> D(Phenotype)


classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;
class A,B,C,D normal;
```

The non-causal association is induced by the backdoor path[^backdoor_note]:

$$
Q \gets \text{Subpopulation} \to P \to \text{Phenotype}.
$$


### Environmental population stratification

It is common for different subpopulations to be exposed to different environments.  These different environments  may differentially affect the phenotype of interest.  This phenomenon is called environmental population stratification.  On its own, environmental population stratification does not confound GWAS.  


However, if both genetic and environmental population stratification are present,  environmental stratification can combine with genetic stratification to induce non-causal GWAS associations.  Having both genetic and environmental stratification is common: genetically distinct people often inhabit distinct environments. An instance of combined environmental and genetic stratification is illustrated below:


``` mermaid
graph LR
C(Subpopulation) --> D(Environment);
C --> E(SNP)
D --> B(Phenotype)


classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;
class B,C,D,E normal;
```

In the scenario illustrated by the diagram, a non-causal association will be induced between the SNP and the phenotype due to the backdoor path: 


$$
\text{SNP}\gets\text{Subpopulation}\to \text{Environment} \to \text{Phenotype}.
$$

Even in the extreme case where the phenotype is entirely environmental and does not depend on genetics at all, combined genetic and environmental population stratification can induce widespread genotype-phenotype associations.


## Adjusting for stratification

We have established that although associations in genetic studies are causally privileged, they can still be confounded by population stratification.  So what can be done?  There are a variety of techniques to mitigate the effects of population stratification.


### Controlling for PCs

The classical causal inference strategy to remove confounding is to adjust for the confounder, which in our case is subpopulation membership.  This strategy is illustrated in the causal diagram below, where (as is traditional in the causal inference literature) we draw a box around conditioned variables.


``` mermaid
graph LR
C(Subpopulation) --> D(Environment);
C --> E(SNP)
D --> B(Phenotype)


classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;
class B,D,E normal;
class C conditioned;
```

By conditioning on the subpopulation, we break the non-causal association between the phenotype and the SNP.



Unfortunately, human population structure is sufficiently complex that it is impossible to mathematical describe it in full detail.  Thus, we must use a proxy.  Population genetics research indicates a person's subpopulation membership can be well-approximated by the allocation of their genotype to genetic principal components[@price2006principal]. Thus, a strategy to approximately adjust for confounding due to population stratification is to condition on genetic principal components.  In the following derivation, we follow Hoffman (2013)[@hoffman2013correcting].


- Let $N$ denote the number of study participants.
- Let $M$ denote the number of genetic variants.
- Let  $X\in\mathbb{R}^{N\times M}$ denote the genotype matrix, normalized to have column means of 0.
- Let $y\in\mathbb{R}^N$ be the phenotype vector, normalized to have mean 0.
- Let $x_j\in\mathbb{R}^N$ be the $j$th column of $X$.
- Let $\hat\beta_j\in\mathbb{R}$ be the scalar marginal regression coefficient of the $j$th genetic variant. This is the quantity we will report in our GWAS summary statistics file for variant $j$.
- Let $\epsilon\in\mathbb{R}^N$ be the random vector of residual environmental and genetic effects.
- Let $X=USV^T$ be the singular value decomposition of $X$.  Thus $U,V^T \in\mathbb{R}^{N\times N}$ are orthogonal  and $S\in\mathbb{R}^{N\times N}$ is diagonal.
- Let $q\in\mathbb{Z}_{++}$ be the number of principal components we retain.  Let $U_{1:q}\in\mathbb{R}^{n\times q}$ be the matrix formed from the first $q$ columns of $U$.
- Let $\sigma^2_e>0$ denote the scale of the residual effects.

The PC-controlled marginal GWAS regression for genetic variant $j$ is:

$$
\begin{align}
y &= \hat\beta_j x_j + U_{1:q} \omega + \epsilon\\
\epsilon  &\sim \mathcal{N}(0, \sigma^2_e I).
\end{align}
$$

We estimate $\hat\beta_j$ by maximum likelihood:

$$
\begin{align}
\hat\beta_j,\omega &= \operatorname*{argmax}_{\hat\beta_j,\omega} \mathcal{N}(y|\hat\beta x_j + U_{1:q} \omega,\sigma^2_e I)
\end{align}
$$

where $\mathcal{N}(y|\mu,\Sigma)$ denotes the multivariate normal density with mean $\mu$ and covariance $\Sigma$ evaluated at $y$.

$\hat\beta_j$ is retained as the marginal GWAS-effect estimate, while $\omega$ is discarded as a nuisance parameter. In this way, we estimate the marginal GWAS association of variant $j$ while controlling for genetic principal components and thus approximately controlling for population stratification. 



### LMMs

Linear mixed models (LMMs) are another popular method to control for population stratification. Here, I will explain them following Hoffman's derivation[@hoffman2013correcting], which clarifies their connection to the PC-control approach. Using the same notation as above, consider the following model for the marginal gwas effect of variant $j$.

$$
\begin{align}
y &= \hat\beta x_j + R\gamma + \epsilon \label{gamma_form}\\
\epsilon  &\sim \mathcal{N}(0, \sigma^2_e I)\\
\gamma & \sim \mathcal{N}(0, \sigma^2_\gamma I ),
\end{align}
$$

where
 

- $R :=U_{1:q}S_{1:q,1:q} \in \mathbb{R}^{N \times q}$ and $S_{1:q,1:q} \in \mathbb{R}^{q\times q}$ is the first $q$ rows and columns of $S$.
- $\gamma \in \mathbb{R}^q$ is the represents the population structure effect.   In contrast to $\omega$ above,  $\gamma$ has a Bayesian prior.


We fit this model via maximum likelihood:

$$
\begin{align}
\hat\beta_j &= \operatorname*{argmax}_{\hat\beta_j} \int \mathcal{N}(y | \hat\beta_j x_j + R\gamma, \sigma^2_e I ) \mathcal{N}(\gamma|0, \sigma^2_\gamma)  \, \mathrm{d}\gamma
\end{align}
$$


Comparing this model to the PC-control model of the previous section, the following points are salient:

- In both cases, we use principal components to control for subpopulation membership.
- With the direct PC-control model, we are limited in the number of principal components we can include.  Including too many may result in a model where the number fit parameters approaches or exceeds $N$, the number of study participants, resulting in poor conditioning or non-uniqueness.  With the LMM we face no such restriction.
- With the LMM, more variable principal components can have a larger effect on the phenotype.  In contrast, in the direct PC-control model all components equal.


While the formulation $(\ref{gamma_form})$ is useful for revealing the connection between LMMs and direct PC control, it is not how LMMs are typically written.  To convert $(\ref{gamma_form})$ to standard LMM form, pick $q=M$ and define $\alpha:= R \gamma$.  By the properties the multivariate normal distribution[^mvnormal_note], 

$$
\begin{align}
\alpha &\sim \mathcal{N}(0, \sigma^2_\gamma RR^T)\\
&= \mathcal{N}(0, \sigma^2_\gamma U S^2 U^T )\\
&= \mathcal{N}(0, \sigma^2_\gamma U SV^T V S U^T )\\\
&=\mathcal{N}(0, \sigma^2_\gamma X X^T)\\\
\end{align}
$$

Thus the restated LMM becomes:

$$
\begin{align}
y &= \hat\beta_j x_j + \alpha + \epsilon\\
\alpha &\sim \mathcal{N}(0,\sigma^2_\gamma K)\\
\epsilon &\sim\mathcal{N}(0, \sigma^2_e I)
\end{align}
$$

Where $K:=XX^T\in\mathbb{R}^{n\times n}$ is called the "genetic relatedness matrix" whose $(i,j)$ entry measures genetic similarity between study participants $i$ and $j$.  This standard is usually used in presentations of LMMs. 



### LOCO 

While LMMs are effective at controlling for population stratification, if care is not taken they can unduly reduce GWAS statistical power.  This reduction in statistical power can occur for two separate reasons: proximal contamination and ascertainment bias.

We begin with proximal contamination. Consider again the equation $(\ref{gamma_form})$. Note that for large values of $q$, we can easily have


$$
\begin{align}
\mathrm{span}(x_j) \approx \mathrm{Range}(R)
\end{align}
$$

where $\mathrm{Range}(R)$ denotes the subspace spanned by the columns of $R$.  This causes a statistical problem resembling multicolinearity: the part of the phenotype vector $y$ that lies in $\mathrm{span}(x_j)$ can be predicted either using $\hat\beta_j$ or $\gamma$.  Simulations[@yang2014advantages] show that this ambiguity can significant reduce statistical power.


The above-described problem is called "proximal contamination" because it results from the inclusion of variant $j$ and nearby variants in close LD with it in the matrix $X$ from which $R$ is constructed.


The standard solution to proximal contamination is the exclusion of these problematic proximal variants when $R$ is constructed.  One approach is called LOCO (leave one chromosome out) LMM. In this approach is $R$ is replaced by $R_{-\mathrm{chr}(j)}$, in which we construct $R$ without the chromosome including variant $j$.  Thus $(\ref{gamma_form})$ is replaced by


$$
\begin{align}
y&= \hat\beta_j x_j R_{-\mathrm{chr}(j)}\gamma + \epsilon\\
\epsilon  &\sim \mathcal{N}(0, \sigma^2_e I)\\
\gamma & \sim \mathcal{N}(0, \sigma^2_\gamma I ),
\end{align}
$$

Thus the LMM term $R_{-\mathrm{chr}(j)}\gamma $ is different for variants on different chromosomes.  Simulations and theoretical analysis[@yang2014advantages] suggest that this approach recovers the power lost by the standard LMM method.  Moreover, since population stratification produced multi-chromosme effects, the exclusion of one chromosome does not affect the ability of the LMM to control for population stratification


### Ascertainment Bias


### REGENIE 

todo

### LDSC


todo


[^backdoor_note]: See _Chapter 7: Confounding_ in Hernan and Robins[@hernan2010causal] for a discussion of backdoor paths.

[^mvnormal_note]: See _Section 4.9: Multivariate normal distribution_ in Grimmet and Stirzaker[@grimmett2020probability].

[//]: # (A key advantage of genetic studies over non-genetic epidemiological studies is that genetic studies are causally privileged.  Specifically, genetic studies benefit from the following advantages:)

[//]: # ()
[//]: # (1.  A person's genes are fixed at birth. Therefore, when we detect an association of a genotype with a phenotype observed later in life, we can be confident there is no reverse causation.  That is, the direction of causality is from the genotype to the phenotype, and not the reverse.)

[//]: # (2. )

[//]: # ()
[//]: # (``` mermaid)

[//]: # (graph LR)

[//]: # (A[Genotype] --> C[Phenotype];)

[//]: # (B[Environment] --> C;)

[//]: # (```)
