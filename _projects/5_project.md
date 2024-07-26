---
layout: page
title: Percolation of words
description: A study of percolation of arbitrary words for site percolation on non-amenable graphs. July 2024 - present. 
img: 
importance: 3
category: Projects
---
This project is along with Ritvik Radhakrishnan who is a Ph.D. student in the group of Prof. Vincent Tassion at ETH Zurich. The problem that we will be exploring was suggested by Prof. Tassion. 

### Setup and Notation
The setup is Bernoulli site percolation a transitive (maybe quasi-transitive) graph $$G$$. We start by some notation first, let $$W = \{0,1\}^{\mathbb{N}}$$ be the set of all infinite binary sequences. We call this set the set of all words. Now we say that a word $$w$$ is $$\textbf{seen}$$ at a vertex $$v$$ if there exists a path $$v, v_1, v_2 \cdots $$ such that $$w(v_i) = X(v_i)$$, where $$X(v)$$ is the state of the site $$v$$. We let $$W(v) = \{w \in W: w \ \text{is seen at v}\}$$ and let $$W_{\infty} = \{w: w \ \text{is seen somewhere}\}.$$ Now the primary question is that are all words seen somewhere, i.e. is $$W_{\infty} = W$$ a.s.? This question was first asked by Harry Kesten and Itai Benjamini for the Euclidean and Triangular lattice. Following the notation of Nolin, et.al we call a word $$w \ \textbf{exceptional}$$ if it is not seen anywhere.



### Previous Work

### Problem Statement
Consider Bernoulli site percolation on a non-amenable transitive graph with $$p_c^{site}<1/2$$. Show that there are no exceptional words at $$p=1/2$$.


### Comments
$$\textbf{What leads to the problem?}$$ The motivating factor is that since the result was fairly easy to show in high dimensions, it should be simple for a graph with infinite dimension (i.e a non-amenable graph). Although it is interesting to note that this is might not always be case, for instance the mean-field criticality in complete generality is not known for non-amenable graphs. However for $$\mathbb{Z}^d, d \geq 11$$ we know mean-field criticality due to the breakthrough work of Hara and Slade. 

$$\textbf{Hermon and Hutchcroft}$$. Since we are in the supercritical regime, it is important to understand the series of papers by Hermon, Hutchcroft detailing the nature of supercritical percolation on non-amenable graphs. One of the key quantities that we are concerned with is the so-called intrinsic radius, in particular if we have a sharp enough control on the tail of the intrinsic radius we will be through. The series of paper by Hutchcroft do indeed show that an exponential decay of the tail of the intrinsic radius happens, however the problem is that the constants are extremely weak and so we don't have a sharp enough decay. 

$$\textbf{Perturbative result}.$$ A paper by Schramm and Benjamini resolves this problem for non-amenable graphs with vertex Cheeger constant $$> 2$$. Schramm and Benjamini show that such graphs have embeddings of a tree with the root having degree = 2 and all other vertices having degree 4. Since we are supercritical for that tree as well the result then follows by the Kesten-Benjamini result for trees. A similar argument can be used to get other perturbative results for small enough spectral radius, etc. 

$$\textbf{Division into two classes.}$$ Famously the study of percolation on non-amenable setting has been often divided into unimodular and non-unimodular graphs. Popular examples are the BLPS result of $$\theta{(p_c)} = 0$$ for non-amenable unimodular graphs, Hutchcroft's paper showing $$p_c < p_u$$ for non-unimodular graphs. Hence, perhaps a similar division would help with this problem as well. Although it is not clear which would be the easier case. I expect the unimodular case to be simpler. 