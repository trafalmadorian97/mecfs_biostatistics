# Liability Threshold Model

## Overview

Modern statistical genetics treats the genetic contribution to common disease as quantitative and continuous. In contrast, many common disease phenotypes are most naturally modeled as binary variables: a person either has the disease or does not.  How can we bridge the gap between binary phenotypes our continuous theory of genetic risk? The liability threshold model is simple, elegant approach.


Under the liability threshold model, the genetic contribution to disease risk $G$ combines additivity with environmental contributions $E$ to produce the latent disease liability variable $L$. When $L$ exceeds the threshold $T$, the patient has the disease phenotype $Y=1$.  Otherwise, the patient has the normal phenotype $Y=0$.  In equations, we have


$$
\begin{align}
L&=G+E\\
Y&= 1_{L>T}.
\end{align}
$$


It is typical to model $G$ and $E$ as independent normal random variables.  This allows the liability threshold model to analyzed via the theory of the [truncated normal distribution](https://en.wikipedia.org/wiki/Truncated_normal_distribution).


A key advantage of the liability threshold model is that it allows many statistical techniques originally developed for quantitative phenotypes to be applied to the binary phenotypes. In essence, one just applies such a technique to the underlying liability $L$.