# Liability Threshold Model

## Overview

Modern statistical genetics treats the genetic contribution to common disease as quantitative and continuous. In contrast, many common disease phenotypes are binary: a person either has the disease or does not.  How can we bridge this gap? The liability threshold model is simple, elegant solution.


Under the liability threshold model, the genetic contribution to disease $G$ combines additively with the environmental contribution $E$ to produce the latent disease liability $L$. When $L$ exceeds the threshold $T$, the patient expresses the disease phenotype $Y=1$.  Otherwise, the patient expresses the normal phenotype $Y=0$.  In equations, we have


$$
\begin{align}
L&=G+E\\
Y&= 1_{L>T}.
\end{align}
$$


It is typical to model $G$ and $E$ as independent normal random variables.  This allows the liability threshold model to analyzed via the theory of the [truncated normal distribution](https://en.wikipedia.org/wiki/Truncated_normal_distribution).


A key advantage of the liability threshold model is that it allows many statistical techniques originally developed for quantitative phenotypes to be applied to the binary phenotypes. One just applies such a technique to the underlying liability $L$ instead of the observed phenotype $Y$.