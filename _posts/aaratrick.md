---
layout: post
title: Every group of order 120 has a subgroup of index 3 or 5
date: 2024-07-15
# description: A solution to Problem 1C.4 from Isaacs's _Finite Group Theory_
categories: problems 
tags: group-theory
giscus_comments: true
# thumbnail: assets/img/9.jpg
---

I have been reviewing basic algebra over this summer as preparation for my PhD quals, and took this opportunity to solve some problems from Isaacs's _Finite Group Theory_. All of the exercises in section 1C, which deals with the Sylow theorems, were reasonably easy, **except** for 1C4, which asks us to prove the statement given in the title. I felt quite good after actually solving this problem, and so I thought it would be apt to post my solution as the inaugural post for this blog. 

Acknowledgement: I had basically covered some trivial cases one night, and called up my good friend Saraswata to cross-check if I had made any silly errors till then. I went to bed afterwards without making any serious progress, but somehow when I sat down with the problem again the next morning, I was able to essentially bypass the case where I had gotten stuck; thanks to Saraswata for listening to me ramble at 2 am :laughing:

#### Notation:

For a finite group $$G$$ and a prime $$p$$, we denote the set of $$p-$$Sylow subgroups of $$G$$ by $$\mathrm{Syl}_p(G)$$ and also set $$n_p(G) = \# \mathrm{Syl}_p(G)$$. We omit the mention of $$G$$ if it's clear from context. We will use the following generalization of the Sylow subgroup counting theorem, proved in Isaacs as Theorem 1.16:

> Suppose that $$G$$ is a finite group such that $$n_p(G) > 1$$. 
> Let $$S$$ and $$T$$ be Sylow-$$p$$ subgroups of $$G$$ such that $$\vert S \cap T \vert$$ is maximal.
> Then, 
> $$ n_p(G) \equiv 1 \pmod{[S:S \cap T]}$$.
{: .block-tip }

A nice example application of this result is discussed in the book right after the statement, and I encourage everyone to take a look at that. Essentially, we use the result to get stronger bounds on the index $$[S:S \cap T]$$, which in turn enables us to give arguments about normality of Sylow subgroups and so on. Let's now dive into the actual solution.

#### The solution:
Let $$G$$ be a group of order 120. Note that $$120 = 2^3 \cdot 3 \cdot 5$$. By the Sylow theorems, we thus have

- $$ n_2 \equiv 1 \mod 2 $$ and $$ n_2 \mid 15 $$, so that $$ n_2 \in \{1,3,5,15\} $$.
- $$ n_3 \equiv 1 \mod 3 $$ and $$ n_3 \mid 40 $$, so that $$ n_3 \in \{1,4,10,40\} $$.
- $$ n_5 \equiv 1 \mod 5 $$ and $$ n_5 \mid 24 $$, so that $$ n_5 \in \{1,6\} $$.

We will deal with the 4 possible values of $$n_2$$ separately. Suppose first that $$n_2=1$$, and let $$N \vartriangleleft G$$ be the Sylow-2 subgroup. Then $$\vert G/N \vert = 15$$. But every group of order $$pq$$, where $$p < q$$ are distinct primes such that $$q \not\equiv 1 \mod p$$, is known to be cyclic, and so we get $$G/N \simeq \mathbb{Z}/ 15\mathbb{Z}$$. We also know that this has a subgroup $$K/N$$ of order 3, and hence, we get a subgroup $$K \leq G$$ of index 5.

We assume $$n_2 > 1$$ henceforth. If either one of $$n_3,n_5$$ is 1, we pick this $$3-$$ or $$5-$$Sylow subgroup $$P$$ and consider the product $$PQ$$, where $$Q \in \mathrm{Syl}_2$$: this is a subgroup because of the normality of $$P$$, and has order either 24 or 40 (because $$P \cap Q = \{e\}$$), and thus index either 5 or 3. In either case, we are done; we are now left the cases where $$n_5=6$$ and $$n_3>1$$. 

Consider a *non-maximal* $$P \in \mathrm{Syl}_2$$, i.e, there is a subgroup $$H$$ such that $$P < H < G$$. But note that $$\vert P \vert = 8 $$ must divide $$\vert H \vert$$, which gives us the only possibility is that $$\vert H \vert $$ is either 24 or 40, and so we are done in this case as well. Therefore, we can assume all $$2-$$Sylow subgroups are maximal, without any loss of generality.

By maximality, $$N_G(S)$$ is either $$S$$ or $$G$$ for all $$S \in \mathrm{Syl}_2$$. As $$n_2 > 1$$, we in fact get $$N_G(S) = S$$ for all $$2-$$Sylow subgroups $$S$$. But, the Sylow theorems then give us

$$n_2 = [G: N_G(S)] = [G:S] = 15$$

Let's take stock of where we are: we have shown the existence of the required subgroup in all but 3 cases, namely when $$n_2 = 15, n_5 = 6$$ and $$n_3>1$$. The time is now ripe to use the result mentioned previously. We pick $$S,T \in \mathrm{Syl}_2$$ such that $$D := S \cap T$$ is of maximal order. The result then gives us $$n_2 = 15 \equiv 1 \mod [S:D] $$. Note that $$[S:D] = \frac{8}{\vert D \vert}$$. As $$D$$ must be properly contained in both $$S$$ and $$T$$, it has order less than 8. This leads to the only possibility that $$[S:D] = 2$$, i.e, $$\vert D \vert = 4$$.

As $$D < S$$ is a subgroup of index 2, it is normal. By the symmetry of our arguments in $$S$$ and $$T$$ so far, we actually have $$D$$ is normal in both of them, i.e, $$S,T \subseteq N_G(D)$$. Note that equality here would imply $$S = T$$, which then means $$\vert D \vert = 8$$, a contradiction! Hence, by the maximality of $$2-$$Sylow subgroups, $$N_G(D) = G$$. Thus, we get $$D \vartriangleleft G$$.

We now consider the group $$\bar{D} = G/D$$, which has order $$30 = 2 \cdot 3 \cdot 5$$. *Suppose* we have managed to find a subgroup $$H/D \leq \bar D$$ of order 6 or 10: in either case, we get a subgroup $$H \leq G$$ whose index is

$$[G:H] = \frac{[G:D]}{[H:D]} \in \{3,5\} $$

Therefore, should we manage to prove the existence of such a subgroup of $$\bar D$$, we would finish the proof. Applying the Sylow theorems to $$\bar D$$, we get $$n_3(\bar D) \in \{1,10\}$$ and $$n_5(\bar D) \in \{1,6\}$$. Now (repeating the argument given in the beginning) if either one of $$n_3, n_5$$ is 1, we can again take the product of this normal Sylow subgroup with any $$2-$$Sylow subgroup of $$\bar D$$ to get a subgroup of order 6 or 10. Hence, it remains to show that at least one of $$n_3(\bar D)$$ and $$n_5(\bar D)$$ must be 1.

We argue this last bit by contradiction: suppose $$n_3(\bar D) = 10$$ and $$n_5(\bar D) = 6$$. Because the exponent of both 3 and 5 in $$\vert \bar D \vert$$ is 1, we can count the number of elements having orders 3 and 5:

- For order 3: $$\# \{x \in \bar D \mid x^3 = \bar e, x \neq \bar e\} = n_3(\bar D) \cdot 2 = 20$$
- For order 5: $$\# \{x \in \bar D \mid x^5 = \bar e, x \neq \bar e\} = n_5(\bar D) \cdot 5 = 24$$

Adding these two counts gives 44, which exceeds the order of $$\bar D$$! Hence, this situation cannot happen, and so we are done. $$\blacksquare$$

***

Do let me know if you find any mistakes or typos by commenting below, or even by [emailing](mailto:ukg7ef@virginia.edu) me. Hope you found this post useful!