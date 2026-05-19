#set text(size: 12pt)

#align(center)[
  #text(24pt, weight: "bold")[DFCP Notes]
]

= 1.1
- Allele: gene variant at a locus on a chromosome
- Diploid: humans have 2 sets of chromosomes, 1 from mom and 1 from dad
- Haplotype: a sequence on a chromosome, genotype = 2 haplotypes
- Phased data: separates parental haplotypes
- biallelic marker: 2 possible alleles
  - diploid, minor can occur 0 (homozygous major allele), 1 (heterozygous allele), 2 (homozygous minor allele) times
  - Single Nucleotide Polymorphism (SNP): 2 possible base pairs at a locus
- Imputation: biallelic markers so binary (major / minor) $N times L$ matrix

= 1.2
- Kingman's Coalescent
  - time going backwards
  - all pairs merge at rate $1 / N_e$, where $N_e$ is the population size $PP["choose same parent"] = 1/N_e$
  - time until next merge with $k$ lineages $~ "Exp"(lambda_k = 1/N_e binom(k, 2))$
  - unions genetic material
- coalescent with recombination
  - kingman's coalescent + recombination
  - recombination: $~ "Exp"(lambda_k = rho k \/ 2)$, split lineage in half
  - splits genetic material
- spatial repr of coalescent with recombination
  - instead of time axis, move along chromosome
  - recombination event = change tree
  - not markov
- mutations: Poisson process with rate $theta k(t)$, uniform locations $l ~ "Uniform"(0, 1)$
- sequentially markov coalescent (SMC): markov approx of spation coalescent
  - start of chromosome: genealogy sampled from kingman's coalescent
  - markov jump process
    - jumps happen with rate $rho L(l) \/ 2$, total branch length $L$ (sum of all edges) and location $l$
    - choose point uniformly on geneaology: pick branch (edge) based on length, uniformly on branch
    - regraft: choose time of reattachment by rate $rho / 2 k(t)$, extending cut branch back, then choose random lineage to attach to
  - prune and regraft moving to right on chromosome
- PAC (Product of Approximate Conditionals) model
  - HMM, states $Z_(i,l) in {1..i-1}$, copying from previous sequence $Z_(i,l)$
  - for sequence $i$: $i-1$ states (previous sequence)
  - transition with probability $c_l$ to uniform random state
  - noisy copying: copy with probability $theta$, otherwise random
  - first sequence is just drawn uniform randomly from $2^L$ binary haplotypes
  - problem: not exchangeable, depends on order of sequences considered
    - can average over random orderings
- fastPHASE (location dependent HMM)
  - at each position, copies one of $K$ latent haplotype instead of previous sequence
  - transition with prob $c_l$
  - $pi_(l,k)$ to transition to latent haplotype $k$ at location $l$ (if transition occurs)
  - learned allele emission frequencies
  - exchangeable
  - I THINK THIS ONE IS GOOD
- STRUCTURE (sequence dependent HMM)
  - transition prob based only on sequence idx

= CRP
partiton $cal(R) ~ "CRP"(R, alpha)$

$n$th customer:
- $PP["joins new table"] = alpha / (n - 1 + alpha)$
- $PP["joins existing table with " \#T "people"] = (\#T) / (n-1 + alpha)$ \

$
A &= {{1}, {2, 3, 7}, {4, 5}, {6}} \

PP[A] &= alpha / (1-1+alpha) dot alpha / (2-1+alpha) dot 1 / (3-1+alpha) dot alpha / (4-1+alpha) dot 1 / (5-1+alpha)dot alpha / (6-1+alpha) dot 2 / (7-1+alpha) \

&= alpha^(\#A) dot 1 / ((alpha) (alpha + 1) ... (alpha + \#R - 1)) dot product_(a in A) (\#a-1)!

$


$
cal(R) ~ "CRP"(R, alpha) \

PP[cal(R) = A] =
Gamma(alpha) / Gamma(alpha + \#R) alpha^(\#A) product_(a in A) Gamma(\#a)
$

CRP is exchangeable. Probability does not depend on order of arrival, just on table sizes.

== Discounted CRP
$n$th customer
- $PP["joins new table"] = (alpha + d K) / (n-1 + alpha)$

- $PP["joins existing table"] = (\#a - d) / (n-1 + alpha)$

where $K$ is the number of existing tables, $d$ is the discount parameter, $alpha$ is the concentration parameter


Kramp's symbol: $[x]^n_d = (x)(x+d)(x+2d)... (x+(n-1)d)$

denom: $(1-1+alpha)(2-1+alpha)...(\#R-1+alpha) = (alpha)(alpha+1)...(alpha+\#R-1) = [alpha]_1^(\#R)$

new tables: $alpha^(\#A) -> alpha(alpha+d)(alpha+2d)...(alpha+(\#A-1)d) = [alpha]_d^(\#A)$

join existing $a$: $(1)(2)...(\#a - 1) -> (1-d)(2-d)...(\#a-1-d) = [1-d]_1^(\#a-1)$

$
PP[A] &= [alpha]_d^(\#A) / [alpha]_1^(\#R) product_(a in A) [1-d]_1^(\#a-1) \
&= [alpha+d]_d^(\#A-1) / [alpha+1]_1^(\#R-1) product_(a in A) [1-d]_1^(\#a-1)
$

Prior partition block growth
- undiscounted: $\#R = cal(O)(alpha log n)$
- discounted: $\#R = cal(O)(n^d)$

= Fragmentation
$cal(Q) ~ "Frag"(cal(R), alpha=0, d)$

For each cluster in $cal(R)$, CRP to partition that cluster into further fragments.

CRP: $PP(A) = ([alpha+d]_d^(\#A-1)) / ([alpha+1]_1^(\#R-1)) product_(a in A) [1-d]_1^(\#a-1)$

$
PP(cal(Q)) &= product_(a in cal(R)) [([d]_d^(\#F_a-1)) / ([1]_1^(\#a-1)) product_(b in F_a) [1-d]_1^(\#b-1)] \

&= product_(a in cal(R)) [(Gamma(\#F_a) d^(\#F_a-1)) / (Gamma(\#a) Gamma(1-d)^(\#F_a)) product_(b in F_a) Gamma(\#b - d)] \

&= (d^(\#cal(Q) - \#cal(R))) / (Gamma(1-d)^(\#cal(Q))) (product_(a in cal(R)) Gamma(\#F_a) / Gamma(\#a)) (product_(b in cal(Q)) Gamma(\#b - d))
$

$F_a = {b in cal(Q): b subset.eq a}$ are fragments of $a$.

= Coagulation
$cal(R) ~ "Coag"(cal(Q), alpha \/ d, 0)$

CRP to partition $cal(Q)$, join subsets of $cal(Q)$ to make a coarser partition.

CRP: $PP[A] = Gamma(alpha) / Gamma(alpha + \#R) alpha^(\#A) product_(a in A) Gamma(\#a)$

$
PP[cal(R)] &= Gamma(alpha\/d) / Gamma(alpha\/d + \#cal(Q)) (alpha\/d)^(\#A) product_(a in cal(R)) Gamma(\#C_a)
$

$C_a = {b in cal(Q): b subset.eq a}$ are all subsets coagulated into $a$.

= Leave-one-out conditionals
$a_i, b_i$ are the respective cluster assignments of $i$ in $cal(R), cal(Q)$.

== CRP
$
cal(R) ~ "CRP"(R, alpha, d=0) \

PP[a_i=a | cal(R)^(-i)] = cases(
  (\#a) / (\#R -1 + alpha) "if" a in cal(R)^(-i) "(join existing cluster)",
  alpha / (\#R -1 + alpha) "if" a = emptyset "(join new cluster)",
)
$

== Frag
$
cal(Q) ~ "Frag"(cal(R), alpha=0, d) \

PP[b_i = b | a_i = a, cal(R)^(-i), cal(Q)^(-i)] = cases(
  (\#b - d) / (\#a) wide &"if" a in cal(R^(-i)) and b in cal(R^(-i)) wide &"(join existing cluster)",
  (\#F_a d) / (\#a) wide &"if" a in cal(R^(-i)) and b = emptyset wide &"(join new cluster)",
  1         wide &"if" a = emptyset and b = emptyset wide &"(stay in singleton cluster)"
)
$

== Coag
$
cal(R) = "Coag"(cal(Q), alpha/d, 0) \

PP[a_i=a|b_i=b, cal(R)^(-i), cal(Q)^(-i)] = cases(
  1 wide & "if" a in cal(R)^(-i) and b in C_a wide & "b coaged into a, i must be in a",
  (d \#C_a) / (alpha + d \#Q^(-i)) wide & "if" a in cal(R)^(-i) and b = emptyset wide & "singleton to cluster" PP prop \#C_a,
  alpha / (alpha + d \#Q^(-i)) wide & "if" a = emptyset = b wide & "singleton to singleton" PP prop alpha / d
)
$

= DFCP
$cal(Q) ~ "CRP"(cal(R), 0, d)$

$cal(R) ~ "CRP"(cal(Q), alpha \/ d, 0)$

$d -> 0: cal(R_l) = cal(R_(l+1))$ b/c Frag chooses existing, Coag chooses new

Generative model:
- initial partition: $cal(R_1) ~ "CRP"(R, alpha, 0)$
- fragment: $cal(Q_l) ~ "CRP"(R_l, 0, d_l)$
- coagulate: $cal(R_(l+1)) ~ "CRP"(Q_l, alpha \/ d_l, 0)$
- $beta_l ~ "Beta"(gamma_l / 2, gamma_l / 2)$
- for each cluster $a$: $theta_(a l) ~ "Bernoulli"(beta_l)$
- for all observations (SNPs) in cluster at location: $x_(i l) = theta_(a l)$

other priors:
- $log alpha ~ N(m, v)$
- $log d_l ~ "Uniform"(log d_min, 0)$
- $log gamma_l ~ "Uniform"(log gamma_min, 0)$

$
PP[cal(R)_(1...L), cal(Q)_(1...L-1)] = &

PP[cal(R_1)] (product_(l=1)^(L-1) PP[cal(Q)_l|cal(R)_(l-1)]) (product_(l=1)^(L-1) PP[cal(R)_(l+1)|cal(Q)_l]) \

= & Gamma(alpha) / Gamma(alpha + N) alpha^(\#cal(R)_1) product_(a in cal(R)_1) Gamma(\#a) \

& dot product_(l=1)^(L-1) [(d_l^(\#cal(Q)_l - \#cal(R)_l)) / (Gamma(1-d_l)^(\#cal(Q)_l)) (product_(a in cal(R)_l) Gamma(\#F_a) / Gamma(\#a)) (product_(b in cal(Q)_l) Gamma(\#b - d_l))] \

& dot product_(l=1)^(L-1) [Gamma(alpha\/d_l) / Gamma(alpha\/d_l + \#cal(Q)_l) (alpha\/d_l)^(\#cal(R)_(l+1)) product_(a in cal(R)_(l+1)) Gamma(\#C_a)] \

= & Gamma(alpha) / Gamma(alpha + N) med alpha^(sum_(l=1)^L \# cal(R)_l) \

& dot product_(a in cal(R)_1) Gamma(\#F_a) med dot product_(a in cal(R)_L) Gamma(\#C_a) \

& dot product_(l=2)^(L-1) Gamma(\#F_a) / Gamma(\#a) Gamma(\#C_a) \

& dot product_(l=1)^(L-1) (d_l^(\#cal(Q)_l - \#cal(R)_l - \#cal(R)_(l+1)) Gamma(alpha \/ d_l)) / (Gamma(1-d_l)^(\#cal(Q)_l) Gamma(alpha \/ d_l  + \#cal(Q)_l))

(product_(b in cal(Q)_l) Gamma(\#b - d_l))
$


#let todo(body) = box(
  fill: rgb("#ff6666"),
  stroke: rgb("#cc4444"),
  inset: (x: 6pt, y: 4pt),
  [*TODO:* #body],
)
#todo[Question: 4.10 eqs dont' include prior probs?]
