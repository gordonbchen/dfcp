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

$a_l, b_l$ are the clusters for sequence $i$ at location $l$ in $cal(R), cal(Q)$.

$
cal(R)_1 ~ "CRP"(R, alpha, 0) \

PP[a_1 = a | cal(R)^(-i)_1] =& cases(
  (\#a) / (N-1 + alpha)      wide & a in cal(R)^(-i)    wide & "join existing",
  (alpha) / (N-1 + alpha)    wide & a = emptyset        wide & "join new"
) \ \


cal(Q)_l ~ "CRP"(cal(R)_l, 0, d_l) \
PP[b_l = b | a_l = a, cal(R)^(-i)_l, cal(Q)^(-i)_l] =& cases(
  (\#b - d_l) / (\#a)            wide & a in cal(R)^(-i)_l and b in cal(Q)^(-i)_l  wide & "cluster to cluster",
  (\#F_l (a) med d_l) / (\#a)    wide & a in cal(R)^(-i)_l and b = emptyset  wide & "cluster to singleton",
  1                      wide & a = emptyset = b      wide & "singleton to singleton"
) \ \


cal(R)_(l+1) ~ "CRP"(cal(Q)_l, alpha \/ d, 0) \

PP[a_(l+1) = a | b_l = b, cal(R)^(-i)_(l+1), cal(Q)^(-i)_l] =& cases(
  (d_l \#C_l (a)) / (alpha + d_l \#Q^(-i)_l) wide & a in cal(R)^(-i)_(l+1) and b = emptyset       wide & "singleton to cluster",
  alpha / (alpha + d_l \#Q^(-i)_l)  wide & a = emptyset = b                              wide & "singleton to singleton",
  1          wide & a in cal(R)^(-i)_(l+1) and b in C_l (a)        wide & "follows into cluster"
)
$

$F_l (a) = {b in cal(Q)_l : b subset.eq a}$ are fragements of $a in cal(R)_l$.

$C_l (a) = {b in cal(Q)_l : b subset.eq a}$ are fragments that get coagulated into $a in cal(R)_(l+1)$.

== Messages
messages: prob of observations to the right after taking sequence $i$ out

$
m_C^l (a) =& PP[x_(i, l+1:L) | a_l = a, cal(R)^(-i)_(l:L), cal(Q)^(-i)_(l:L-1)] \

=& sum_(b in cal(Q)^(-i)_l union {emptyset}) PP[b_l = b | a_l = a, cal(R)^(-i)_l, cal(Q)^(-i)_l] med  m^l_F (b) \

=& cases(
  m_F^l (b=emptyset) wide a = emptyset,

  1 / (\#a) [sum_(b in F_l (a)) (\#b -d_l) m_F^l (b) + \#F_l(a) d_l m_F^l (b=emptyset)]  wide a in cal(R)^(-i)_l
)
\ \ \


m_F^l (b) =& PP[x_(i, l+1:L) | b_l = b, cal(R)^(-i)_(l:L), cal(Q)^(-i)_(l:L-1)] \

=& sum_(a in cal(R)^(-i)_(l+1) union {emptyset}) PP[a_(l+1) = a | b_l = b, cal(R)^(-i)_l, cal(Q)^(-i)_l] med Lambda(x_(i,l+1) | a) med m^(l+1)_C (a) \

=& cases(
  Lambda(x_(i,l+1) | a) med m^(l+1)_C (a)  wide  b in cal(Q)^(-i)_l  wide  a "st" b in C_l (a),

  1 / (alpha + d_l \#Q^(-i)_l) [alpha Lambda(x_(i,l+1) | emptyset) med m^(l+1)_C (emptyset) 
  + sum_(a in cal(R)^(-i)_(l+1)) d_l \#C_l (a) Lambda(x_(i,l+1) | a) med m^(l+1)_C (a)]
  wide b = emptyset
)
$

== Likelihood $Lambda$
$a != emptyset$, match cluster allele: $Lambda(x | a) = delta(x_(i,l) = theta_(a,l))$ \ \

$a = emptyset$, $beta_l ~ "Beta"(gamma_l / 2, gamma_l / 2)$, $theta_(a,l) ~ "Bernoulli"(beta_l)$

$
Lambda(x | a=emptyset) = cases(
  (gamma_l\/2 + n_(1,l)) / (gamma_l + n_(0,l) + n_(1, l)) wide x = 1,
  (gamma_l\/2 + n_(0,l)) / (gamma_l + n_(0,l) + n_(1, l)) wide x = 0,
)
$

$n_(1,l) = \#{a in cal(R)^(-i)_l : theta_(a, l) = 1}$ is the number of clusters that emit the major allele at location $l$.

== Posteriors
$
PP[a_1 = a|x_i, cal(R)^(-i)_(1:L), cal(Q)^(-i)_(1:L-1)]

prop & PP[a_1=a|cal(R)^(-i)_1] PP[x_(i 1)|a_1 = a] PP[x_(i, 2:L)|a_1=a, cal(R)^(-i)_(1:L), cal(Q)^(-i)_(1:L-1)] \

prop & PP[a_1=a|cal(R)^(-i)_1] dot Lambda(x_(i 1)|a) dot m_C^1 (a) \ \ \


PP[b_l = b|x_i, cal(R)^(-i)_(1:L), cal(Q)^(-i)_(1:L-1)]

prop & PP[b_l = b|a_l = a, cal(R)^(-i)_l, cal(Q)^(-i)_l] PP[x_(i, l+1:L) | b_l = b, cal(R)^(-i)_(l:L), cal(Q)^(-i)_(l:L-1)] \

prop & PP[b_l = b|a_l = a, cal(R)^(-i)_l, cal(Q)^(-i)_l] m_F^l (b) \ \ \


PP[a_l = a|x_i, cal(R)^(-i)_(1:L), cal(Q)^(-i)_(1:L-1)]

prop & PP[a_l = a|b_(l-1) = b, cal(R)^(-i)_l, cal(Q)^(-i)_(l-1)] PP[x_(i, l)|a] PP[x_(i, l+1:L) | a_l = a, cal(R)^(-i)_(l:L), cal(Q)^(-i)_(l:L-1)] \

prop & PP[a_l = a|b_(l-1) = b, cal(R)^(-i)_l, cal(Q)^(-i)_(l-1)] dot Lambda(x_(i, l) | a) dot m_C^l (a) \ \ \
$

== Slice sampling $alpha, d_l, gamma_l$
TODO: many mistakes here
- 4.9 R instead of Q in denom (not critical)
- 4.10 doesn't have priors
- 4.20 and 4.21: extra -L in alpha exp and +1 in d_l exponent


$
PP[alpha|cal(R), cal(Q), d] prop

& PP[alpha]

dot Gamma(alpha) / Gamma(alpha + N) med alpha^(sum_(l=1)^L \# cal(R)_l)

dot product_(l=1)^(L-1) Gamma(alpha \/ d_l) / Gamma(alpha \/ d_l  + \#cal(Q)_l)) \



PP[d_l|cal(R), cal(Q), alpha] prop

& PP[d_l]

dot (d_l^(\#cal(Q)_l - \#cal(R)_l - \#cal(R)_(l+1)) Gamma(alpha \/ d_l)) / (Gamma(1-d_l)^(\#cal(Q)_l) Gamma(alpha \/ d_l  + \#cal(Q)_l))

dot product_(b in cal(Q)_l) Gamma(\#b - d_l) \ \


PP[theta_(a l)|beta_l] =& beta_l^(theta_(a l)) (1-beta_l)^(1-theta_(a l)) \

PP[theta_l|beta_l] =& beta_l^(n_(1 l)) (1-beta_l)^(n_(0 l)) wide n_(1 l) = \#{a in cal(R)_l : theta_(a l) = 1} \

PP[beta_l|gamma_l] =& Gamma(gamma_l) / Gamma(gamma_l / 2)^2 beta_l^(gamma_l / 2 - 1) (1-beta_l)^(gamma_l / 2 - 1)  \


PP[theta_l|gamma_l] =& integral_0^1 PP[theta_l|beta_l]PP[beta_l|gamma_l] d beta_l \

=& integral_0^1 Gamma(gamma_l) / Gamma(gamma_l / 2)^2 beta_l^(gamma_l / 2 + n_(1 l) - 1) (1-beta_l)^(gamma_l / 2 + n_(0 l) - 1)  d beta_l \

=& (Gamma(gamma_l) Gamma(gamma_l/2+n_(1 l)) Gamma(gamma_l/2+n_(0 1))) / 
(Gamma(gamma_l/2)^2 Gamma(gamma_l + n_(1 l) +n_(0 l))) \


PP[gamma_l|theta_l, cal(R)_l] prop& PP[gamma_l]
  dot (Gamma(gamma_l) Gamma(gamma_l/2+n_(1 l)) Gamma(gamma_l/2+n_(0 1))) /
  (Gamma(gamma_l/2)^2 Gamma(gamma_l + n_(1 l) +n_(0 l))) \
$

QUESTION: sequences in a cluster must agree for every single emission, otherwise likelihood is 0. this seems too hard a constraint. why not softer?


= Slice Sampling

$
u ~& "Uniform"(0, f(x)) \

A_u =& {x' : f(x') >= u} \

x' ~& "Uniform"(A_u)
$

Computing $A_u$ via stepping out and shrinkage

$
r ~ "Uniform"(0, 1) \

L = x - r w \

R = L + w \
$

Step out
- while $f(L) > u$, $L <- L - w$
- while $f(R) > u$, $R <- R + w$

Shrinkage
- $x' ~ "Uniform(L, R)"$
- if $f(x') >= u$: accept
- else
  - if $x' > x: R = x'$
  - if $x' < x: L = x'$

== Log space
$
u ~& "Uniform"(0, f(x)) \

u =& U f(x) wide U ~ "Uniform"(0, 1) \

log u =& log U + log f(x) \ \ \

PP[-log U <= x] =& PP[log U >= -x] \
=& PP[U >= e^(-x)] \
=& 1 - e^(-x) wide "so" -log U ~ "Exp"(lambda = 1) \ \ \

E ~& "Exp"(lambda=1) \

log u =& log f(x) - E
$
