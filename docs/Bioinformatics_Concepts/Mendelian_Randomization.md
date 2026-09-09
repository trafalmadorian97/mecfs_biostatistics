

# Mendelian Randomization

[//]: # (To discuss: Formalism of Causal inferece.  How the three instrumental variables assumptions dont let us do causal inference without a `parametric` model.  Wald ratio.  Wald ratio error.  Other methods.)
## Categories of Data Science
Hernan et al.[@hernan2019second] assert that the core activities of the data scientist are, in order of increasing complexity:


1. **Description**: Computing summary statistics, to enable large datasets to be easily comprehended. An example would be producing a chart summarizing the characteristics of heart disease patients.
2. **Prediction**: Estimating the distribution of one variable conditional on other variables.  An example would be calculation of a person's odds of heart disease, conditional on their demographics, genetics, and [LDL](../Analysis/Low_Density_Lipoprotein/Verma_et_al/LDL_Overview.md) level.
3. **Causal Inference**:   Reasoning about counterfactuals.   An example would be evaluating the evidence in support of a statement like: "If all patients with high LDL were put on statins, occurrence of heart disease would be reduced by $X$ percent".  


Here, we focus on causal inference.


[//]: # (Since causal effects are the core of scientific knowledge, there is understandable interest in specifying the criteria according to which valid causal inferences can be made.)

## Methods of causal inference

The randomized controlled trial (RCT) is the gold standard causal inference method. RCTs have the advantage of requiring few assumptions, but the disadvantages of being costly and time-consuming.  This motivates alternative causal inference techniques, which require more assumptions, but are cheap and fast.  Mendelian Randomization (MR) is such a technique.



## Mendelian Randomization

As an illustrative example, consider the question of how LDL affects heart disease. We might notice from epidemiological studies that people with heart disease tend to have higher LDL than people without it. We recall, however, that correlation is not causation. Thus we apply MR, in hopes of determining whether there is a true causal effect.


### The three main MR assumptions

The figure below from Hartley et al.[@hartley2022guide] summarizes MR:

![mr-diagram](https://github.com/user-attachments/assets/afcab2e5-b9d0-423f-a43b-0b4f47b172e4)

We apply MR to estimate the causal effect of an exposure (e.g. LDL) on an outcome (e.g. heart disease).  MR requires a genetic variant ("Genetic Instrument") associated with the exposure.  

The validity of MR depends on three main assumptions:

- **IV1**: The instrument must be associated with the exposure.
- **IV2**: The instrument cannot share a cause with the exposure or the outcome.
- **IV3**: The instrument cannot be causally associated with the outcome, except via the exposure.

In MR, natural variation in the genetic instrument plays a role analogous to the random treatment assignment in an RCT.  The analogy between MR and an RCT is illustrated in the diagram below, from Zuber et al.[@zuber2022combining]:

![mr-analogy](https://github.com/user-attachments/assets/39aac59f-c0c6-4bce-91a9-80fc910b2753)



### The Fourth MR Assumption

Many texts on MR emphasize IV1, IV2, and IV3. However, these conditions are not sufficient to uniquely determine the causal effect of the exposure on the outcome (See Hernan and Robins Chapter 16[@hernan2010causal]).  Another assumption is required.  There are numerous possible variants forms of this fourth assumption, which can largely be grouped into 

- Homogeneity assumptions: roughly, the causal effect of the exposure on the outcome is the same for all individuals.
- Monotonicity assumptions: roughly, the instrument moves the exposure in the same direction for all individuals.


MR texts tend to ignore the need for a fourth assumption because most practical applications of MR assume that the outcome $Y$ is a linear function of the exposure $X$:

$$
\begin{align}
Y &= \beta_{Y,X}X +\alpha_1 + \epsilon &\text{ where }\epsilon \sim N(0,\sigma_1) 
\end{align}
$$

Under this linear model, homogeneity holds automatically.




## Pitfalls in Mendelian Randomization

### Horizontal Pleiotropy

A common source of error in Mendelian Randomization occurs when the genetic instrument affects the outcome through a causal pathway that does not involve the exposure.  This scenario, which is known as "horizontal pleiotropy", results in a violation of IV3.  If the strength of the causal effect through the alternative pathway is large, the results of Mendelian Randomization can be misleading.

The diagram below illustrates horizontal pleiotropy:

``` mermaid
graph LR
A[Genetic Instrument] --> B[Protein 1];
B --> C[Exposure];
C --> D[Outcome];
A --> E[Protein 2];
E ----> D;


classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;

class A,B,C,D,E normal;
```

A genetic variant affects the levels of two proteins.  One protein affects the outcome through the exposure, while the other affects the outcome independently of the exposure. If we believe the [Omnigenic Model](Omnigenic_Model.md), we should expect this kind of horizontal pleiotropy to be relatively common.


The risk of horizontal pleiotropy is magnified when the connection between the genetic variant and the exposure is complex and indirect.  An example would be a study in which the genetic instrument affects neurodevelopment, and the exposure is tobacco use.  The risk is reduced when the connection is straightforward and direct.  An example would be an MR study in which the exposure is the plasma level of a a protein, and the genetic instrument is a cis-regulatory variant for that protein (a cis-pQTL).


## Methods of Mendelian Randomization

### The Wald Ratio


todo

## Links

[Talk on MR by Dr. Jean Morrison](https://www.youtube.com/watch?v=kxhmy3EWH8Q&t=34m13s)

