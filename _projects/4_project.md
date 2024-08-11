---
layout: page
title: Percolation thresholds under coverings.
description: Summer 2024, TIFR.
img:
importance: 3
category: Projects
---
This project was carried out as part of the Visiting Students Research Programme, at the Tata Institute of Fundamental Research under the guidance of Prof. Subhajit Goswami. I gave a presentation at TIFR, which could be found [here](https://ishaan44.github.io/assets/pdf/VSRP_Presentation.pdf).

It is a well known fact that for a d-regular graph $$G, p_c(G) \geq \frac{1}{d-1}$$. Now if for a d-regular graph we set $$p_c = \frac{1}{d-1}$$ does it imply that the graph is a tree? We answer this question here. 

### Quasi-transitive $$d$$-regular graphs with $$p_c = \frac{1}{d-1}$$ are trees.
$$\textbf{Main Result.}$$ Let $$G$$ be a $$d$$-regular quasi-transitive graph then if $$G$$ is not isomorphic to a tree. Then,  $$p_c(G) > \frac{1}{d-1}.$$

### Non quasi-transitive counter-example.
We start with following results,

$$\textbf{Result 1.}$$ Let G be a graph of degree d. Then $$p_c(G) \geq \frac{1}{d-1}.$$

$$Proof:$$ The proof is simple and follows from the first moment method. Let $$p_c(G) < \frac{1}{d-1}$$ and fix a $$p_c < p < \frac{1}{d-1}$$. Now let $$X_n = \Sigma_{\gamma \in \Gamma_n}\mathbb{1}(\gamma \ \text{is open})$$, then $$ \mathbb{P}_p(0 \leftrightarrow\partial \Lambda_n) = \mathbb{P}_p(X_n > 0). $$ 

But, $$\mathbb{P}_p(X_n > 0) \leq \mathbb{E}(X_n) = \Sigma_{\gamma \in \Gamma_n}\mathbb{P}(\gamma \ \text{is open}) = p^n |\Gamma_n|.$$ Now $$|\Gamma_n| \leq d(d-1)^{n-1.}$$ 

Therefore,
$$\mathbb{P}_p(X_n > 0) \leq p^n d(d-1)^{n-1} = (p(d-1))^{n-1} pd \xrightarrow{} 0$$ as $$n \rightarrow{} \infty.$$ But for $$p>p_c$$ this limit should be stricly positive, this is a contradiction. Thus $$p_c(G) \geq \frac{1}{d-1}.$$



$$\textbf{Result 2.}$$ Let G be a d-regular tree. Then $$p_c(G) = \frac{1}{d-1}$$


