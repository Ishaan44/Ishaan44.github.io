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

<h2>Using the Critical Threshold to Characterize a Tree</h2>
<p><strong>Abstract:</strong> In this article, we study the critical threshold for an arbitrary quasi-transitive d-regular graph.
We show that for such graphs G, p<sub>c</sub>(G) = 1/(d-1) if and only if G is a tree. This gives a neat probabilistic characterization 
for regular trees. We also give counterexamples for the non-quasi-transitive case.</p>

<h3>1 Introduction</h3>
<p>We consider independent bond percolation on a graph G, i.e., we retain edges with probability p and throw them away with probability 1-p.
By a simple first moment, a graph with maximal degree Δ &lt; ∞ satisfies p<sub>c</sub>(G) ≥ 1/(Δ-1). For trees of degree Δ, by a dual 
second moment upper bound, we can get p<sub>c</sub> ≤ 1/(Δ-1), thus implying that p<sub>c</sub> = 1/(Δ-1).</p>

<p>The question we ask is the following: Consider percolation on a Δ-regular graph G, are trees the only graphs with p<sub>c</sub> = 1/(Δ-1)?
For the rest of this article, we restrict ourselves to percolation on regular graphs. In light of our question, we show that if one assumes
quasi-transitivity, then for any Δ-regular graph G which is not a tree, p<sub>c</sub> &gt; 1/(Δ-1).</p>

<h3>Theorem 1</h3>
<p>Let G be a quasi-transitive regular graph of degree Δ. Then p<sub>c</sub>(G) ≥ 1/(Δ-1), and equality holds if and only if G is a tree.
It is important to note that the requirement of being a quasi-transitive graph is essential.</p>

<h3>2 Percolation on Trees</h3>
<h4>2.1 Branching number and the critical point for trees</h4>
<p>Suppose T = (V,E) is an infinite locally finite tree with root O. We imagine the tree T as growing downward from the root O. For x, y ∈ V,
we write x ≤ y if x is on the shortest path from O to y; and T<sub>x</sub> for the subtree of T containing all the vertices y with y ≥ x.</p>

<p>For a vertex x ∈ V, we denote by d(x,O) the graph distance from O to x. We want to understand the critical point for a tree, motivated 
by the comparison from Galton-Watson branching processes. This naturally leads us to the study of the average number of branches coming out 
of a vertex, which is called the branching number of a tree. To rigorously define this, we use conductances and flows on trees. For each edge e, 
we define the conductance of an edge to be c(e) := λ<sup>-|e|</sup>, where |e| denotes the distance of the edge e from the root O. It is 
natural to define conductances decreasing exponentially with the distance since trees grow exponentially.</p>

<p>If λ is very small, then due to large conductances, there is a non-zero flow on the tree satisfying 0 ≤ θ(e) ≤ λ<sup>-|e|</sup>. While 
increasing the value of λ, we observe a critical value λ<sub>c</sub> above which such a flow does not exist. This is precisely the branching number.</p>

<h4>Theorem 2. (R. Lyons, [Lyo90])</h4>
<p>Let T be a tree, then p<sub>c</sub>(T) = 1/br(T), where br(T) is the branching number of the tree.</p>

<h3>3 Non-quasi-transitive Counterexamples</h3>
<p>Let T = (V,E) be the graph of black edges shown in the figure, define G := T∪{e} = (V,E∪{e}) to be the graph obtained after adding the red edge e. 
In the figure, T is the underlying tree rooted at O and G is a 3-regular graph. We claim that p<sub>c</sub>(G) = 1/(Δ-1), for Δ = 3. Similar counters 
can be constructed for d ≥ 4. Now for the tree T, |T<sub>1</sub>| = 3, |T<sub>2</sub>| = 6, |T<sub>3</sub>| = 10, after this point, every point has 2 branches 
coming out, so |T<sub>3+n</sub>| = 10× 2<sup>n</sup>. Therefore, grT = 2.</p>

<p>T is obviously subperiodic since for all x from Level 3 onwards, T<sub>x</sub> is exactly T<sub>A</sub>, so we can define the function f(v) = ϕ(v), where ϕ is the isomorphism 
between T<sub>x</sub> and T<sub>A</sub>. For the rest of the vertices, we can simply choose the identity map. Thus, T is 2-subperiodic. Now by Theorem 3, 
brT = grT = 2. Thus, p<sub>c</sub>(T) = 1/2, but p<sub>c</sub>(G) ≤ p<sub>c</sub>(T) = 1/2. However, by a first moment bound, p<sub>c</sub>(G) ≥ 1/2. Thus, 
p<sub>c</sub>(G) = 1/2.</p>

<h3>4 Proof of the Theorem</h3>
<p>We now show that if G is a quasi-transitive Δ-regular graph, then p<sub>c</sub>(G) > 1/(Δ-1). The idea is to cover every Δ-regular graph 
by a Δ-regular tree. The proof of this does not require quasi-transitivity. Next, we use the strict monotonicity (under coverings) result 
of Martineau and Severo [MS19]. To apply this result, we want our covering map to satisfy certain properties.</p>

<h3>References</h3>
<ul>
    <li>R. Lyons, [Lyo90], Random walks and percolation on trees. The Annals of Probability, 18(3):931–958, 1990.</li>
    <li>F. Severo, S. Martineau, [MS19], Strict monotonicity of percolation thresholds under covering maps, 2019.</li>
</ul>

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


