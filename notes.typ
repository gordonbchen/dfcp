#set text(size: 12pt)

#let todo(body) = box(
  fill: rgb("#ff6666"),
  stroke: rgb("#cc4444"),
  inset: (x: 6pt, y: 4pt),
  [*TODO:* #body],
)


#align(center)[
  #text(24pt, weight: "bold")[DFCP Notes]
]



= CRP
partiton $cal(R) ~ "CRP"(R, alpha)$

$n$th customer:
- $PP["joins new table"] = alpha / (n - 1 + alpha)$
- $PP["joins existing table with " \#T "people"] = (\#T) / (n-1 + alpha)$ \

$
A &= {{1}, {2, 3, 7}, {4, 5}, {6}} \

PP[A] &= alpha / (1-1+alpha) dot alpha / (2-1+alpha) dot 1 / (3-1+alpha) dot alpha / (4-1+alpha) dot 1 / (5-1+alpha)dot alpha / (6-1+alpha) dot 2 / (7-1+alpha) \

&= alpha^(\#A) dot 1 / ((alpha) (alpha + 1) ... (alpha + \#R - 1)) dot product_(a in A) (\#a-1)! \ \
$

$
cal(R) ~& "CRP"(R, alpha) \

PP[cal(R) = A] =& Gamma(alpha) / Gamma(alpha + \#R) alpha^(\#A) product_(a in A) Gamma(\#a)
$

CRP is exchangeable. Probability does not depend on order of arrival, just on table sizes.



= Discounted CRP
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



= Messages
messages: prob of observations to the right after taking sequence $i$ out

$
m_C^l (a) =& PP[x_(i, l+1:L) | a_l = a, cal(R)^(-i)_(l:L), cal(Q)^(-i)_(l:L-1)] \

=& sum_(b in cal(Q)^(-i)_l union {emptyset}) PP[b_l = b | a_l = a, cal(R)^(-i)_l, cal(Q)^(-i)_l] med  m^l_F (b) \

=& cases(
  m_F^l (b=emptyset) wide a = emptyset,

  1 / (\#a) [sum_(b in F_l (a)) (\#b -d_l) m_F^l (b) + \#F_l (a) d_l m_F^l (b=emptyset)]  wide a in cal(R)^(-i)_l
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



= Likelihood $Lambda$
$a != emptyset$, match cluster allele: $Lambda(x | a) = delta(x_(i,l) = theta_(a,l))$ \ \

$a = emptyset$, $beta_l ~ "Beta"(gamma_l / 2, gamma_l / 2)$, $theta_(a,l) ~ "Bernoulli"(beta_l)$

$
Lambda(x | a=emptyset) = cases(
  (gamma_l\/2 + n_(1,l)) / (gamma_l + n_(0,l) + n_(1, l)) wide x = 1,
  (gamma_l\/2 + n_(0,l)) / (gamma_l + n_(0,l) + n_(1, l)) wide x = 0,
)
$

$n_(1,l) = \#{a in cal(R)^(-i)_l : theta_(a, l) = 1}$ is the number of clusters that emit the major allele at location $l$.



= Posteriors
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



= Slice sampling $alpha, d_l, gamma_l$
many mistakes here
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



#pagebreak()
= ME paper

== CRP
$
cal(R)_l ~& "CRP"(R, alpha, 0) \

PP[cal(R)_l=A] =& Gamma(alpha) / Gamma(alpha + N)  alpha^(\#A)  product_(a in A) Gamma(\#a)
$

== Frag
$
cal(Q)_l|cal(R)_l ~& "Frag"(cal(R)_l, 0, d) \

PP[cal(Q)_l|cal(R)_l]

=& product_(a in cal(R)_l) [(Gamma(\#F_a) d^(\#F_a-1))/(Gamma(\#a) Gamma(1-d)^(\#F_a)) product_(b in F_a) Gamma(\#b - d) ] \

=& d^(\#cal(Q)_l - \#cal(R)_l) / Gamma(1-d)^(\#cal(Q)_l) [product_(a in cal(R)_l) Gamma(\#F_a)/Gamma(\#a)] [product_(b in cal(Q)_l) Gamma(\#b - d)] \

F_a =& {b in cal(Q)_l: b subset.eq a}
$

== Coag
$
cal(R)_(l+1)|cal(Q)_l ~& "Coag"(cal(Q)_l, alpha\/d, 0) \

PP[cal(R)_(l+1)|cal(Q)_l] =& Gamma(alpha\/d) / Gamma(alpha\/d + \#cal(Q))  (alpha\/d)^(\#cal(R)_(l+1))  product_(a in cal(R)_(l+1)) Gamma(\#C_a) \

C_a =& {b in cal(Q)_l: b subset.eq a}
$

== Generative model
$
alpha ~& "Gamma"(tau_1, tau_2) \

d_l ~& "Beta"(v_1, v_2) \ \


cal(R)_l ~& "CRP"(R, alpha, 0) \

cal(Q)_l|cal(R)_l ~& "Frag"(cal(R)_l, 0, d_l) \

cal(R)_(l+1)|cal(Q)_l ~& "Coag"(cal(Q)_l, alpha\/d_l, 0) \ \


gamma_l ~& "Gamma"(phi.alt_1, phi.alt_2) \

beta_l|gamma_l ~& "Dirichlet"(gamma_(l,1) ... gamma_(l,2)) \ 

theta_(l a)|beta_l ~& "Categorical"(beta_l) \ \


x_(i l)|a_(i l) =& theta_(l a_(i l))
$




#pagebreak()
= Maximization Expectation
- maximization: $C_i^* = "argmax"_(C_i) EE_q [log p(C, theta|cal(D))]$
- expectation: $log q_i (theta_i) prop EE_(q_(-i))[ log p(theta_i |C^*, theta_(-i), cal(D)) ]$


Mean field approx: $q(alpha, gamma, d) = q(alpha) product_(l=1)^L p(gamma_l) product_(l=1)^(L-1) p(d_l)$


$
theta^* =& "argmin"_theta [ "KL"(q_theta (z) || p(z|x) )] \

=& "argmin"_theta [ EE_(q_theta) [log q_theta (z)] - EE_(q_theta) [log p(z|x) ] ]\

=& "argmax"_theta [ EE_(q_theta) [log p(z,x) ] - EE_(q_theta) [log q_theta (z)] ]\ \

"ELBO"(q_theta) =& EE_(q_theta) [log p(z,x) ] - EE_(q_theta) [log q_theta (z)] \ \
$


$
log p(C, gamma, alpha, d, X)

=& log( product_(l=1)^L Gamma(4 gamma_l) / Gamma(4 gamma_l + \#cal(R)_l) product_(k in {A,T,C,G}) Gamma(gamma_l + n_k) / Gamma(gamma_l) ) \

+& log[ Gamma(alpha) / Gamma(alpha + N) alpha^(\#cal(R)_1) product_(a in cal(R)_1) Gamma(\#a) ]\

+& log[ product_(l=1)^(L-1) ((d_l^(\#cal(Q)_l - \#cal(R)_l)) / (Gamma(1-d_l)^(\#cal(Q)_l)) (product_(a in cal(R)_l) Gamma(\#F_a) / Gamma(\#a)) (product_(b in cal(Q)_l) Gamma(\#b - d_l))) ]\

+& log [product_(l=1)^(L-1) (Gamma(alpha\/d_l) / Gamma(alpha\/d_l + \#cal(Q)_l) (alpha\/d_l)^(\#cal(R)_(l+1)) product_(a in cal(R)_(l+1)) Gamma(\#C_a)) ]\

+& log "Gamma"(tau_1, tau_2) + sum_(gamma_l) log "Gamma"(phi.alt_1, phi.alt_2) + sum_(d_l) log "Beta"(v_1, v_2)
$

$
=& log Gamma(alpha) / Gamma(alpha + N) + (sum_(l=1)^L \#cal(R)_l + tau_1 - 1) log alpha - tau_2 alpha 

+ tau_1 log tau_2 - log Gamma(tau_1) \

+& sum_(l=1)^(L-1) [
  (\#cal(Q)_l - \#cal(R)_l - \#cal(R)_(l+1) + v_1 - 1) log d_l

  + (v_2 - 1) log (1 - d_l) - log B(v_1, v_2) \

  &- \#cal(Q)_l log Gamma(1-d_l)

  + sum_(b in cal(Q)_l) log Gamma(\#b - d_l)
] \

+& sum_(l=1)^(L-1) log Gamma(alpha\/d_l) / Gamma(alpha\/d_l + \#cal(Q)_l) \

+& sum_(l=1)^L [
  phi.alt_1 log phi.alt_2 - log Gamma(phi.alt_1)

  + (phi.alt_1 - 1) log gamma_l - phi.alt_2 gamma_l

  + log Gamma(K gamma_l) / Gamma(K gamma_l + \#cal(R)_l)

  + sum_(k in {A,T,C,G}) log Gamma(gamma_l + n_k) / Gamma(gamma_l)
] \

+& sum_(a in cal(R)_1) log Gamma(\#a)

+ sum_(l=1)^(L-1) ( sum_(a in cal(R)_l) [ log Gamma(\#F_a) - log Gamma(\#a) ] + sum_(a in cal(R)_(l+1)) log Gamma(\#C_a) ) \
\ \ \
$




#pagebreak()
= Bivariate delta approx
$
f(x,y) approx& f(mu_x, mu_y) \

+& (x-mu_x) f_x (mu_x, mu_y) + (y-mu_y) f_y (mu_x, mu_y) \

+& 1/2 (x-mu_x)^2 f_(x x)(mu_x, mu_y) \

+& (x-mu_x)(y-mu_y) f_(x y)(mu_x, mu_y) \

+& 1/2 (y-mu_y)^2 f_(y y)(mu_x, mu_y) \ \

EE_(x y) [f(x,y)] approx& f(mu_x, mu_y) + 1/2 sigma^2_x f_(x x)(mu_x, mu_y) + "Cov"(x, y) f_(x y)(mu_x, mu_y) + 1/2 sigma^2_y f_(y y)(mu_x, mu_y) \ \ \
$

By mean field approx, $"Cov"(x, y) = 0$.
$
EE_(x y) [f(x,y)] approx& f(mu_x, mu_y) + 1/2 sigma^2_x f_(x x)(mu_x, mu_y) + 1/2 sigma^2_y f_(y y)(mu_x, mu_y)
$


$
(d)/(d d_l) log Gamma(alpha/d_l) =& psi(alpha/d_l) dot (-alpha)/d_l^2 \

(d^2)/(d d_l^2) log Gamma(alpha/d_l) =& psi(alpha/d_l) dot (2 alpha)/d_l^3 + psi_1(alpha/d_l) dot alpha^2 / d_l^4
$




#pagebreak()
= Variational update for $gamma_l$
$
log q^*_(gamma_l) (gamma_l)

prop& log p(gamma_l) + log p(X|C, gamma_l) \

prop& log( gamma_l^(phi.alt_1-1) e^(-phi.alt_2 gamma_l) )

+ log( Gamma(K gamma_l) / Gamma(K gamma_l + \#cal(R)_l) product_(k in {A,T,C,G}) Gamma(gamma_l + n_k) / Gamma(gamma_l) ) \


prop& (phi.alt_1 - 1) log gamma_l - phi.alt_2 gamma_l

+ log Gamma(K gamma_l) / Gamma(K gamma_l + \#cal(R)_l) + sum_(k in {A,T,C,G}) log Gamma(gamma_l + n_k) / Gamma(gamma_l) \
$

laplace approx in log space: same as $alpha$
- $eta = log gamma_l$

- solve $hat(eta) = "argmax"_eta [ hat(h)(eta) ]$ where $hat(h)(eta) = h(e^eta) + eta$

- $q_eta (eta) approx cal(N)(hat(eta), sigma_eta^2 = -1/(hat(h)''(hat(eta))))$

  $hat(h)''(hat(eta)) = hat(gamma)_l^2 h''(hat(gamma)) - 1$

- $q_(gamma_l) (gamma_l) approx "LogNormal"(hat(eta), sigma_eta^2)$ \ \

$
h''(gamma_l) =& (1 - phi.alt_1) / gamma_l^2

+ K^2 [psi_1(K gamma_l) - psi_1(K gamma_l + \#cal(R)_l)] + sum_(k in {A,T,C,G}) [ psi_1(gamma_l + n_k) - psi_1(gamma_l) ]
$




#pagebreak()
= Variational update for $alpha$
$
log q^*_alpha (alpha)

prop& log Gamma(alpha) / Gamma(alpha + N) + (sum_(l=1)^L \#cal(R)_l + tau_1 - 1) log alpha - tau_2 alpha \

+& sum_(l=1)^(L-1) EE_(d_l)[ log Gamma(alpha\/d_l) / Gamma(alpha\/d_l + \#cal(Q)_l) ] \ \ \


prop& log Gamma(alpha) - log Gamma(alpha + N) + (sum_(l=1)^L \#cal(R)_l + tau_1 - 1) log alpha - tau_2 alpha \

+& sum_(l=1)^(L-1) (EE_(d_l)[ log Gamma(alpha\/d_l) ] - EE_(d_l) [log Gamma(alpha\/d_l + \#cal(Q)_l) ] )\ \ \
$


== Delta method for $EE_x [log Gamma(a/x + b)]$
$
\
EE[f(x)] approx& f(mu) + 1/2 sigma^2 f''(mu) \

d/(d x) [log Gamma(a/x + b)] =& psi(a/x + b) dot -a/x^2 \


d^2/(d x^2) [log Gamma(a/x + b)] =& psi_1(a/x + b) dot a^2/x^4 + psi(a/x+b) dot (2a)/x^3 \


EE[log Gamma(a/x + b)] approx& log Gamma(a/mu + b) + 1/2 sigma^2 f''(mu) \
$


== Laplace approx
$
EE[log Gamma(a/x + b)] approx log Gamma(a/mu + b) + 1/2 sigma^2 [psi_1(a/mu + b) dot a^2/mu^4 + psi(a/mu+b) dot (2a)/mu^3] \ \

d/(d a) = 1/mu psi(a/mu + b) + 1/2 sigma^2 [
  psi_2(a/mu + b) dot a^2/mu^5 + psi_1(a/mu+b) dot (4a)/mu^4 + psi(a/mu+b) dot 2/mu^3
] \ \


d^2/(d a^2) = 1/mu^2 psi_1(a/mu + b) + 1/2 sigma^2 [
  psi_3(a/mu + b) dot a^2/mu^6 + psi_2(a/mu + b) dot (6a)/mu^5 + psi_1(a/mu+b) dot 6/mu^4
] \ \
$


$
h''(alpha) =& psi_1(alpha) - psi_1(alpha+N) + (1 - tau_1 - sum_(l=1)^(L) \#cal(R)_l) / alpha^2 \

+& sum_(l=1)^(L-1) (d^2/(d alpha^2) EE_(d_l)[log Gamma(alpha/d_l)] - d^2/(d alpha^2) EE_(d_l)[log Gamma(alpha/d_l + \#cal(Q)_l)])\
$



#pagebreak()
= Laplace approx
- $log p(z) prop h(z)$

- find mode $hat(z) = "argmax"_z h(z)$

- taylor expansion: $h(z) approx h(hat(z)) + 1/2 h''(hat(z)) (z-hat(z))^2$

  - $h(z) approx h(hat(z)) + h'(hat(z)) (z - hat(z)) + 1/2 h''(hat(z)) (z-hat(z))^2$
  - $h'(hat(z)) = 0$ b/c $hat(z)$ is the mode

- let $sigma^2 = - 1/(h''(hat(z)))$

  - $h''(hat(z)) < 0$ since $hat(z)$ is a maximum

- $h(z) approx h(hat(z)) - 1/(2 sigma^2)(z-hat(z))^2$

- $p(z) prop exp(h(z)) approx exp(h(hat(z))) exp(-1/(2 sigma^2) (z-hat(z))^2) ~ cal(N)(hat(z), -1/(h''(hat(z))))$



== Laplace approx in log space
- $alpha > 0$. let $eta = log alpha$

- $q_eta (eta) = q_alpha (e^eta) |(d alpha) / (d eta)| = q_alpha (e^eta) e^eta$

- $q_alpha (alpha) prop exp(h(alpha))$

  $q_eta (eta) prop exp(h(e^eta)) e^eta$

- $log q_eta (eta) prop hat(h) (eta) = h(e^eta) + eta$

- laplace approx
  - solve $hat(eta) = "argmax"_eta hat(h)(eta)$

  - $q_eta (eta) approx cal(N) (hat(eta), sigma_eta^2 = -1/ (hat(h)''(hat(eta))))$

- $hat(h)'(eta) = d/(d eta) [h(e^eta) + eta] = e^eta h'(e^eta) + 1 = alpha h'(alpha) + 1$

  $hat(h)''(eta) = d/(d eta) [h'(e^eta) e^eta + 1] = e^(2eta) h''(e^eta) + e^eta h'(e^eta) =  alpha^2 h''(alpha) + alpha h'(alpha)$ \ \

  $hat(h)'(hat(eta)) = hat(alpha) h'(hat(alpha)) + 1 = 0$, so $h'(hat(alpha)) = -1/hat(alpha)$

  $hat(h)''(hat(eta)) = hat(alpha)^2 h''(hat(alpha)) -1$ \ \

- $q_alpha (alpha) approx "LogNormal"(hat(eta), sigma_eta^2)$

  $mu_alpha = exp(hat(eta) + sigma_eta^2 / 2)$

  $sigma_alpha = (exp(sigma_eta^2) - 1) exp(2 hat(eta) + sigma_eta^2)$




#pagebreak()
= Variational update for $d_l$
$
log q_(d_l)^*(d_l)

prop& (\#cal(Q)_l - \#cal(R)_l - \#cal(R)_(l+1) + v_1 - 1) log d_l + (v_2 - 1) log (1 - d_l) \

-& \#cal(Q)_l log Gamma(1-d_l) + sum_(b in cal(Q)_l) log Gamma(\#b - d_l) \

+& EE_alpha [ log Gamma(alpha\/d_l) / Gamma(alpha\/d_l + \#cal(Q)_l) ] \ \ \


prop& (\#cal(Q)_l - \#cal(R)_l - \#cal(R)_(l+1) + v_1 - 1) log d_l + (v_2 - 1) log (1 - d_l) \

-& \#cal(Q)_l log Gamma(1-d_l) + sum_(b in cal(Q)_l) log Gamma(\#b - d_l) \

+& EE_alpha [ log Gamma(alpha\/d_l) ] - EE_alpha [ log Gamma(alpha\/d_l + \#cal(Q)_l) ] \ \ \
$

Paper correction: $EE_(alpha) [log (alpha\/d_l)^(\#cal(R)_(l+1))] -> -\#cal(R)_(l+1) log d_l$



== Delta method
$
EE_x [ f(x) ] approx& f(mu) + 1/2 sigma^2 f''(mu) \

d^2/(d x^2) log Gamma(a x + b) =& a^2 psi_1(a x + b) \

EE_x [ log Gamma(a x + b) ] approx& log Gamma(a mu + b) + 1/2 sigma^2 a^2 psi_1(a mu + b)
$


== derivs
$
EE_alpha [log Gamma(mu_alpha / d_l + b)] approx log Gamma(mu_alpha / d_l + b) + 1/2 sigma^2 psi_1(mu_alpha / d_l + b) 1/d_l^2 \
$

$
d/(d d_l)

=& psi(mu_alpha / d_l + b) dot (-mu_alpha)/d_l^2 + 1/2 sigma^2 [
  psi_2(mu_alpha / d_l + b) (-mu_alpha)/d_l^4 + psi_1(mu_alpha / d_l + b) (-2)/d_l^3
] \

d^2/(d d_l^2)

=& psi_1(mu_alpha / d_l + b) dot mu_alpha^2 / d_l^4 + psi(mu_alpha / d_l + b) dot (2 mu_alpha)/d_l^3 \

+& 1/2 sigma^2 [
  psi_3(mu_alpha / d_l + b) (mu_alpha^2)/d_l^6

  + psi_2(mu_alpha / d_l + b) (6mu_alpha)/d_l^5

  + psi_1(mu_alpha / d_l + b) 6/d_l^4
]
\ \ \
$

$
h''(d_l)

=& (\#cal(R)_l - \#cal(Q)_l + \#cal(R)_(l+1) + 1 - v_1) / d_l^2 - (v_2-1)/ (1-d_l)^2 \

-& \#cal(Q)_l dot psi_1(1-d_l) + sum_(b in cal(Q)_l) psi_1(\#b - d_l) \

+& d^2/(d d_l^2) EE_alpha [ log Gamma(alpha / d_l) ] - d^2/(d d_l^2) EE_alpha [ log Gamma(alpha / d_l + \#cal(Q)_l) ]
$



== Laplace approximation in logit space
- $log q_d (d) prop h(d)$

- let $d = sigma(psi) = 1 / (1+e^(-psi))$. then $psi = log d / (1-d)$

- $q_psi (psi) = q_d (sigma(psi)) |(d d)/(d psi)| = q_d (d) dot d(1-d)$

- $q_d (d) prop exp(h(d))$

  $q_psi (psi) prop exp(h(sigma(psi))) dot d (1-d)$

- $log q_psi (psi) prop hat(h)(psi) = h(d) + log d + log (1-d)$

- laplace approx
  - solve $psi_0 = "argmax"_psi hat(h)(psi)$

  - $q_psi (psi) approx cal(N)(psi_0, sigma^2_psi = -1/(hat(h)''(psi_0)))$

- $
  hat(h)'(psi) =& d/(d psi) [h(d) + log d + log (1-d)] \
  =& d/(d d) [h(d) + log d + log (1-d)] (d d)/(d psi) \
  =& (h'(d) + 1/d - 1/(1-d)) dot d (1-d) \
  =& d (1-d) dot h'(d) + 1-2d \
  $

  $
  hat(h)''(psi) =& d/(d psi) [d (1-d) dot h'(d) + 1 - 2d] \
  =& d/(d d) [d (1-d) dot h'(d) + 1 - 2d] (d d)/(d psi) \
  =& [(1-2d)h'(d) + d(1-d) h''(d) - 2] dot d(1-d) \
  $

- $
  hat(h)'(psi_0) =& d_0 (1-d_0) dot h'(d_0) + 1-2d_0 = 0 \

  h'(d_0) =& (2d_0-1) / (d_0(1-d_0))
  $

  $
  hat(h)''(psi_0) =& [(1-2d_0)h'(d_0) + d_0(1-d_0) h''(d_0) - 2] dot d_0(1-d_0) \
  =& -(1-2d_0)^2 + [d_0(1-d_0)]^2 h''(d_0) - 2 d_0 (1-d_0) \
  $

- $q_psi (psi) approx cal(N)(psi_0, sigma_psi^2)$

  $q_d (d) approx "LogitNormal"(psi_0, sigma_psi^2)$

- 2nd order Delta approx for mean
  $
  EE[d] = EE[sigma(psi)] approx& sigma(psi_0) + 1/2 sigma_psi^2 dot sigma''(psi_0) \
  approx& d_0 + 1/2 sigma_psi^2 (1-2d_0) dot d_0 (1-d_0) \ \ \

  sigma'(psi_0) =& d_0 (1-d_0) \
  sigma''(psi_0) =& (1-2d_0) dot d_0(1-d_0)
  $

- 1st order delta approx for variance
  $
  d =& sigma(psi) approx sigma(psi_0) + sigma'(psi_0) (psi - psi_0) \ \

  "Var"[d] =& "Var"[sigma(psi)] \

  approx& "Var"[sigma(psi_0) + sigma'(psi_0) (psi - psi_0)] \

  approx& [sigma'(psi_0)]^2 "Var"[psi] \

  approx& d_0^2 (1-d_0)^2 sigma_psi^2
  $

- 2nd order delta approx for variance

$
EE[d^2] = EE[sigma(psi)^2]

approx& sigma(psi_0)^2 + 1/2 sigma_psi^2 dot lr(d^2/(d psi^2) [sigma(psi)^2] |)_(psi=psi_0) \

approx& d_0^2 + 1/2 sigma_psi^2 dot (4 d_0 - 6d_0^2) dot d_0 (1-d_0) \ \ \


lr(d/(d psi) [sigma(psi)^2] |)_(psi=psi_0) =& 2 d_0^2 (1-d_0) \

lr(d^2/(d psi^2) [sigma(psi)^2] |)_(psi=psi_0) =& (4 d_0 - 6d_0^2) dot d_0 (1-d_0) \
$

$
"Var"[d] = EE[d^2] - EE[d]^2
approx& d_0^2 + 1/2 sigma_psi^2 dot (4 d_0 - 6d_0^2) dot d_0 (1-d_0) - EE[d]^2 \
$



#pagebreak()
= Maximization
$C_i^* = "argmax"_(C_i) EE_q(Theta) [log p(C, Theta|cal(D))]$

Paper correction: take out and put back 1 sequence at a time.


== Likelihood
$a != emptyset$
$
Lambda(x_(i l) | a) = cases(
  1 wide& x_(i l) = theta_(a l),
  0 wide& "otherwise"
)
$

$a = emptyset$

$
beta_l ~ "Dirichlet"(gamma_l) \

theta_(a l) ~ "Categorical"(beta_l)
$

$
Lambda(x_(i l) = k | a=emptyset) = (gamma_l + n_(k l)) / (K gamma_l + \#cal(R)_l^(-i))
$

$n_(k l) = \#{a in cal(R)^(-i)_l : theta_(a l) = k}$ is the number of clusters that emit the allele $k$ at location $l$.

$K$ = \# alleles.

$
EE_q log Lambda(x_(i l) | a) = cases(
  cases(
    0 wide& x_(i l) = theta_(a l),
    -infinity wide& "otherwise"
  )
  wide& a != emptyset,

  EE_(gamma_l)[ log (gamma_l + n_(k l)) ] - EE_(gamma_l)[ log (K gamma_l + \#cal(R)_l^(-i)) ]
  wide& a = emptyset
)
$


== Viterbi
$
m_F^l (b)

=& max_(a in cal(R)_(l+1)^(-i) union {emptyset}) EE_q [
  log PP[a_(l+1)=a | b_l = b, cal(R)_(l+1)^(-i), cal(Q)_l^(-i)]
  + m_C^(l+1) (a) ] \

=& cases(
  - EE_q [log (alpha + d_l \#Q^(-i)_l)]
  + max_(a in cal(R)_(l+1)^(-i) union {emptyset}) cases(
    EE_alpha [log alpha] + m^(l+1)_C (emptyset) wide& a = emptyset,
    EE_(d_l)[log d_l] + log \#C_l (a) + m^(l+1)_C (a) wide& a in cal(R)^(-i)_(l+1)
  ) wide& b = emptyset,

  m_C^(l+1)(a) wide wide a "st" b in C_l (a) wide& b in cal(Q)_l^(-i)
)


\ \ \
m_C^l (a)

=& EE_q [ log Lambda(x_(i, l)|a) ]
  + max_(b in cal(Q)_l^(-i) union {emptyset}) EE_q [
  log PP[b_l=b | a_l=a, cal(R)_l^(-i), cal(Q)_l^(-i)]
  + m_F^l (b) ]\

=& EE_q [ log Lambda(x_(i, l)|a) ]

+ cases(
  m_F^l (emptyset) wide& a = emptyset,

  -log \#a
  + max_(b in F_l (a) union {emptyset}) cases(
    EE_(d_l)[log (\#b -d_l)] + m_F^l (b) wide& b in F_l (a),
    log \#F_l (a) + EE_(d_l)[log d_l] + m_F^l (emptyset) wide& b = emptyset
  )
  wide& a in cal(R)_l^(-i)
)

\ \ \
m_C^L (a) =& EE_q [ log Lambda(x_(i, L)|a) ]
$


== $EE_alpha [log alpha]$
$
alpha ~& "Gamma"(tau_1, tau_2) \

E_alpha [log alpha] =& psi(tau_1) - ln tau_2 wide wide psi "is digamma function" \ \ \


alpha ~& "LogNormal"(mu_alpha, sigma^2_alpha) \
E_alpha [log alpha] =& mu_alpha
$

== $EE_(d_l)[log d_l]$
$
d_l ~& "Beta"(v_1, v_2) \

EE_(d_l)[log d_l] =& psi(v_1) - psi(v_1 + v_2) \ \ \
$

== Not closed form
$
EE_x [f(x)] approx& f(mu) + 1/2 sigma^2 f''(mu) \

EE_x [log(a x + b)] approx& log(a mu + b) + 1/2 sigma^2 f''(mu) \
approx& log(a mu + b) - 1/2 sigma^2 a^2 / (a mu + b)^2 \

d / (d x) log(a x + b) =&  a / (a x + b) \
d^2 / (d x^2) log(a x + b) =& (-a^2) / (a x + b)^2
$


- $EE_(gamma_l)[ log (gamma_l + n_(k l)) ]$
- $EE_(gamma_l)[ log (K gamma_l + \#cal(R)_l^(-i)) ]$ 
- $EE_(d_l)[log( \#b - d_l)]$
- $EE_(d_l)[log d_l]$


$EE_q [log (alpha + d_l \#Q^(-i)_l)]$
- Let $Y = alpha + d_l \#Q^(-i)_l$

- $
  mu_Y =& mu_alpha + \#Q^(-i)_l mu_(d_l) \
  sigma_Y^2 =& sigma^2_alpha + (\#Q^(-i)_l)^2 sigma^2_(d_l)
  $

  note that $alpha, d_l$ are independent by the mean field assumption.

- $EE_q [log (alpha + d_l \#Q^(-i)_l)] = EE_Y [log Y]$





#pagebreak()
= Predictive checks
- minor allele frequency distribution across variant sites

== Linkage disequilibrium
- 2 biallelic loci: A/a at locus 1, B/b at locus 2

- independent: $p_(A B) = p_(A) p_(B)$

- linkage disequilibrium (not indep): $p_(A B) != p_(A) p_(B)$

  $D = p_(A B) - p_(A) p_(B)$

  $D > 0: wide p_(A B) > p_(A) p_(B) wide$ AB more frequently than expected

- $D = p_(A B) p_(a b) - p_(a B) p_(B a)$

  determinant = 0: rank 1 matrix, outer product of marginals, decomposable = indep

- $r^2 = D^2 / (p_(A) p_(a) p_(B) p_(b))$

  $r^2 = 0 -> D = 0$, uncorrelated

- LD decay curve: $1/(|cal(P)_d|) sum_((i, j) in cal(P)_d) r^2_(i,j)$ for all pairs $cal(P)_d$ in a distance bin $d$




#pagebreak()
= Soft DFCP

== DFCP
$
alpha ~ "Gamma"(tau_1, tau_2) \

d_l ~ "Beta"(v_1, v_2) \ \ \


cal(R)_0 ~ "CRP"(R, alpha, 0) \

cal(Q)_l ~ "Frag"(cal(R)_l, 0, d_l) \

cal(R)_(l+1) ~ "Coag"(cal(Q)_l, alpha\/d_l, 0) \ \ \


gamma_l ~ "Gamma"(phi.alt_1, phi.alt_2) \

beta_l ~ "Dirichlet"(gamma_l) \

theta_(a l) ~ "Categorical"(beta_l) \ \ \


x_(i l) = theta_(a_(i l) l)
$


== Soft DFCP
$
gamma_l ~ "Gamma"(phi.alt_1, phi.alt_2) \

beta_(a l) ~ "Dirichlet"(gamma_l) \ \

x_(i l) ~ "Categorical"(beta_(a_(i l) l))
$ \ \


$
log p(C, gamma, alpha, d, X)


=& log p(cal(R)_1|alpha)

+ sum_(l=1)^(L-1) [ log p(cal(Q)_l|cal(R)_l, d_l) + log p(cal(R)_(l+1)|cal(Q)_l, alpha, d_l)] \ 

+& sum_(i=1)^N sum_(l=1)^L log p(x_(i l)|gamma_l) \

+& log p(alpha|tau_1, tau_2) + sum_(l=1)^(L-1) log p(d_l|v_1, v_2) + sum_(l=1)^L log p(gamma_l|phi.alt_1, phi.alt_2)
$ \

$
sum_(i=1)^N sum_(l=1)^L log p(x_(i l)|gamma_l) = sum_(l=1)^L sum_(a in cal(R)_l)

  log integral_(beta_(a l)) p(beta_(a l)|gamma_l) product_(i in a) p(x_(i l)|beta_(a l))
$ \

$
p(beta_(a l)|gamma_l) =& Gamma(K gamma_l) / Gamma(gamma_l)^K product_(k=1)^K (beta_(a l))_k^(gamma_l-1) \


product_(i in a) p(x_(i l)|beta_(a l))

  =& product_(i in a) product_(k=1)^K (beta_(a l))_k ^ (bb(1){x_(i l) = k})

  = product_(k=1)^K (beta_(a l))_k ^ n_(a l k) \


p(beta_(a l)|gamma_l) product_(i in a) p(x_(i l)|beta_(a l))

=& Gamma(K gamma_l) / Gamma(gamma_l)^K product_(k=1)^K (beta_(a l))_k^(gamma_l + n_(a l k) - 1) \


integral_(beta_(a l)) p(beta_(a l)|gamma_l) product_(i in a) p(x_(i l)|beta_(a l))

=& integral_(beta_(a l)) Gamma(K gamma_l) / Gamma(gamma_l)^K product_(k=1)^K (beta_(a l))_k^(gamma_l + n_(a l k) - 1) \

=& Gamma(K gamma_l) / Gamma(gamma_l)^K  integral_(beta_(a l)) product_(k=1)^K (beta_(a l))_k^(gamma_l + n_(a l k) - 1) \

=& Gamma(K gamma_l) / Gamma(gamma_l)^K (product_(k=1)^K Gamma(gamma_l + n_(a l k))) / Gamma(K gamma_l + \#a) \

=& Gamma(K gamma_l) / (Gamma(gamma_l)^K  Gamma(K gamma_l + \#a))product_(k=1)^K Gamma(gamma_l + n_(a l k)) 
$

$
sum_(i=1)^N sum_(l=1)^L log p(x_(i l)|gamma_l)

=& sum_(l=1)^L sum_(a in cal(R)_l) log integral_(beta_(a l)) p(beta_(a l)|gamma_l) product_(i in a) p(x_(i l)|beta_(a l)) \

=& sum_(l=1)^L sum_(a in cal(R)_l) [
  log Gamma(K gamma_l) / (Gamma(gamma_l)^K  Gamma(K gamma_l + \#a)) 

  + sum_(k=1)^K log Gamma(gamma_l + n_(a l k))
]
$ \


#pagebreak()
= Variational update for $gamma_l$
$
log q_(gamma_l)^* (gamma_l) prop& log p(gamma_l) + log p(X_(a l)|C, gamma_l) \

prop& (phi.alt_1 - 1) log gamma_l - phi.alt_2 gamma_l + sum_(a in cal(R)_l) [
  log Gamma(K gamma_l) / (Gamma(gamma_l)^K  Gamma(K gamma_l + \#a)) 

  + sum_(k=1)^K log Gamma(gamma_l + n_(a l k))
]
$ \

$
h''(gamma_l) =& (1-phi.alt_1)/gamma_l^2 + sum_(a in cal(R)_l) [
  K^2 [ psi_1(K gamma_l) - psi_1(K gamma_l + \#a) ]

  - K psi_1(gamma_l)

  + sum_(k=1)^K psi_1(gamma_l + n_(a l k))
]
$


= Maximization Likelihood
$
Lambda(x_(i l) = k|a_(-i), gamma_l, x_(-i l))

=& integral_(beta_(a l)) p(beta_(a l)|gamma_l, x_(-i l)) p(x_(i l)|beta_(a l)) \

=& integral_(beta_(a l)) Gamma(K gamma_l + \#a_(-i)) / (product_(k=1)^K Gamma(gamma_l + n_(a l k))) product_(k=1)^K (beta_(a l))_k^(gamma_l + n_(a l k) + bb(1){x_(i l) = k} - 1) \

=& Gamma(K gamma_l + \#a_(-i)) / (product_(k=1)^K Gamma(gamma_l + n_(a l k)))

(product_(k=1)^K Gamma(gamma_l + n_(a l k) + bb(1){x_(i l) = k})) / Gamma(K gamma_l + \#a_(-i) + 1) \

=& (gamma_l + n_(a l k)) / (K gamma_l + \#a_(-i)) \ \ \



Lambda(x_(i l) |a = emptyset, gamma_l, x_(-i l)) = 1/K
$



#pagebreak()
= Hyperparams

== $gamma_l$ 
Cluster has $n$ observed alleles, all A.

$
p(x_(i l) = k|a, x^(-i)) =& (gamma_l + n_(a l k)^(-i)) / (K gamma_l + n_(a l)^(-i)) \

p(x_(i l) != A|n_A = n) =& ((K-1) gamma_l) / (K gamma_l + n_(a l)^(-i))
$

prob of mismatch should be $approx$ bit flip ratio?

try $gamma_l ~ "Gamma"(1, 20)$

$
p =& ((K-1) gamma_l) / (K gamma_l + n_(a l)^(-i)) \

p K gamma_l + p n =& (K-1) gamma_l \

p n =& ((1-p) K -1) gamma_l \

(p n) / ((1-p) K -1) &= gamma_l \
$


== $alpha$
$
E[\#cal(R)_l | alpha] approx& alpha log (1 + N / alpha)
$




#pagebreak()
= Trees
- fastsimcoal gives coalescent trees at each position
- given a tree, is the dfcp marginal clustering in it? metric for how far away

== Parsimony score
- given a tree and clusters, how can I assign cluster labels to internal ancestral nodes to explain the
  observed cluster labels with as few changes along edges as possible

- $s_T (cal(X)) = min_(hat(cal(X)) : V(T) -> cal(K)) \#{(u, v) in E(T) : hat(cal(X))(u) != hat(cal(X))(v)}$

  subject to $cal(X)(x) = hat(cal(X))(x)$ on leaves $x$

- optimal parsimony score for $k$ clusters = $k-1$

  excess parsimony score: $s_T (cal(X)) - (k-1)$

- at optimality: clusters represent subtrees

== Fitch's Algorithm
- let $P(v)$ be the parsimony score of a subtree rooted at $v$
- let $S(v)$ be the set of cluster labels for the vertex $v$ that can achieve $P(v)$

base case: at edges, $P(v) = 0$, $S(v) = {C(v)}$

recursion: vertex $v$ has children $a, b$.
- if $S(a) inter S(b) = emptyset$
  - there is no cluster label that can be assigned to $v$ to not increase parsimony
  - $P(v) = P(a) + P(b) + 1$
  - $S(v) = S(a) union S(b)$
- else $S(a) inter S(b) != emptyset$
  - $P(v) = P(a) + P(b)$
  - $S(v) = S(a) inter S(b)$


== Other metrics
+ % of MRCA (Most Recent Common Ancestor): |a| / (\# leaves under MRCA)
  - weighted avg: $sum_a (|a|) / N (|a|) / m(a)$

+ cluster stability: jaccard index $J(C_l, C_(l+1)) = (|B_l inter B_(l+1)|) / (|B_l union B_(l+1)|)$


=== Bad
+ average within cluster coalescent time: but penalizes model for not splitting tree farther
  even if no mutation

+ fixed time clustering: dfcp doesn't have to split the tree at the same height,
  mutations can appear at different times on different branches.



#pagebreak()
= Robust and Interpretable Statistical Genetic Modelling

== statistics independent of reference genealogy
- average \# marginal clusters: $1/L sum_l |cal(R)_l|$

- cluster entropy: $H(C) = -sum_(i in U) (|x_i|)/(|C|) log (|x_i|)/(|C|)$
  - $U$ unique haplotypes in a cluster
  - entropy of haplotypes in a cluster

- partition entropy: $H(cal(R)_l) = -sum_(C in cal(R)_l) (|C|)/N log (|C|)/N$

- $"purity"(cal(R)_l) = 1/N sum_k max_j |C_k inter x_j|$
  - normalized sum of cluster purity, how many in cluster are the same as a representative haplotype

- mutual information: how much more information you get about X by knowing Y.

  $
  I(X, Y) =& H(X) - H(X|Y) \
  =& - sum_x p(x) log p(x) + sum_(x y) p(x y) log p(x|y) \
  =& - sum_x p(x) log p(x) + sum_(x y) p(x y) log (p(x y))/(p(y)) \
  =& - sum_y sum_x p(x y) log p(x) + sum_(x y) p(x y) log (p(x y))/(p(y)) \
  =& sum_(x y) p(x y) log (p(x y))/(p(x) p(y))
  $

  mutual information b/t partition and classes: how much info you get about the cluster by knwoing the haplotype.

  $
  "MI"(cal(R)_l, U) =& sum_(k=1)^K sum_(j=1)^U p(C_k x_j) log (p(C_k x_j)) / (p(C_k) p(x_j)) \
  =& sum_(k=1)^K sum_(j=1)^U (|C_k inter x_j|)/N log N (|C_k inter x_j|) / (|C_k| |x_j|) \
  $

  mutual information maximized for degenerate each haplotype gets its own cluster

  $"NMI" =& ("MI"(cal(R), U)) / sqrt(H(cal(R)) H(U))$
  or
  $"NMI" =& (2"MI"(cal(R), U)) / (H(cal(R)) + H(U))$

== statistics conditional on reference geneaology
- mutual information b/t reference partition and dfcp partition

  $
  "MI"(cal(R)_1, cal(R)_2) =& sum_(k_1) sum_(k_2) p(k_1, k_2) log (p(k_1, k_2)) / (p(k_1) p(k_2)) \
  =& sum_(k_1) sum_(k_2) (|k_1 inter k_2|)/N log N (|k_1 inter k_2|) / (|k_1| |k_2|) \
  $

- adjusted mutual information

  $
  "AMI"(cal(R)_1, cal(R)_2) =& ("MI"(R_1, R_2) - EE["MI"(R_1, R_2)])
    / (max(H(R_1), H(R_2)) - EE["MI"(R_1, R_2)])
  $

  where $EE["MI"(R_1, R_2)]$ is the expected mutual information for random clusters of the same size and \# haplotypes.

- prune and regraft distance, NP hard in general

- marginal tree recovery: $sum_l sum_(C in R) bold(1){C in "Clades"(T_l)}$

- jaccard distance b/t dfcp clusters and marginal trees:
  $sum_l sum_(C in R) max_(S in "Clades"(T_l)) J(S, C)$

- importance scores: downweight degenerate clusters

  $
  f_(l,r)(x|alpha, beta) = Gamma(alpha+beta)/(Gamma(alpha) Gamma(beta))
    ((x-l) / (r-l))^(alpha-1) (1 - (x-l) / (r-l))^(beta-1)
  $

  then using $x = |c|, c in cal(R)$ and $l=1, r=N$, we can downweight degenerate clusters.

  scaled jaccard distance: using $alpha=beta = 1, 2, 5$ penalize degenerate clusters.

  $
  (sum_l sum_(C in R) f_(1,N)(C) max_(S in "Clades"(T_l)) J(S, C))
    / (sum_l sum_(C in R) f_(1,N)(C))
  $

- adjusted rand index

  2 clusterings $A, B$

  $n_(q,r) = |A_q inter B_r|$, shared haplotypes in clusters.

  $binom(n_(q,r), 2)$ is the number of haplotype pairs in both $A_q inter B_r$

  $S_(A B) = sum_(q,r) binom(n_(q,r),2)$, is the number of haplotype pairs in the same cluster between all clusters in $A, B$.

  $T = binom(N, 2)$ is the total number of hapltoype cluster pairs.

  Expected shared pairs is $E = (S_A S_B) / T$

  $"ARI"(A, B) = (S_(A,B) - E) / (1/2 (S_A + S_B) - E)$


== Useful stats
unsupervised
- average \# marginal clusters
- cluster purity
- partition entropy
- NMI: guard against degenerate N clusters.

supervised
- importance score weighted jaccard distance from dfcp cluster to tree



#pagebreak()
#show link: it => underline(text(fill: blue)[#it])
== #link("https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0060123&type=printable")[
  The Effect of Single Recombination Events on Tree Height and Shape
]

- spatial coalescent is not markov?
- 


#pagebreak()
= Noisy DFCP

== DFCP
$
alpha ~& "Gamma"(tau_1, tau_2) \

d_l ~& "Beta"(v_1, v_2) \ \


R_0 ~& "CRP"(R, alpha, 0) \

Q_l ~& "Frag"(R_l, 0, d_l) \

R_(l+1) ~& "Coag"(Q_l, alpha\/d_l, 0) \ \ \


gamma_l ~& "Gamma"(phi.alt_1, phi.alt_2) \

beta_l ~& "Dirichlet"(gamma_l) \

theta_(a l) ~& "Categorical"(beta_l) \ \

x_(i l) =& theta_(a_(i l), l)
$

== Noisy
$
gamma_l ~& "Gamma"(phi.alt_1, phi.alt_2) \

beta_l ~& "Dirichlet"(gamma_l) \

theta_(a l) ~& "Categorical"(beta_l) \ \


epsilon.alt ~& "Beta"(lambda_1, lambda_2) \ \

PP[x_(i l) = k] =& cases(
  k = theta_(a l) wide&    1 - epsilon.alt,
  "otherwise"      wide&    epsilon.alt \/ (K-1)
) \ \
$

$
log p(C, gamma, alpha, d, X)

=& log p(R_1|alpha)

+ sum_(l=1)^(L-1) [ log p(Q_l|R_l, d_l) + log p(R_(l+1)|Q_l, alpha, d_l)] \ 

+& log p(alpha|tau_1, tau_2) + sum_(l=1)^(L-1) log p(d_l|v_1, v_2) + sum_(l=1)^L log p(gamma_l|phi.alt_1, phi.alt_2) \

+& log p(epsilon.alt|lambda_1, lambda_2) \

+& sum_(i=1)^N sum_(l=1)^L log p(x_(i l)|gamma_l, epsilon.alt) \ \
$

$
sum_(i=1)^N sum_(l=1)^L log p(x_(i l)|gamma_l, epsilon.alt) =&

log Gamma(K gamma_l) / Gamma(K gamma_l + \#cal(R)_l)
  + sum_(k in {A,T,C,G}) log Gamma(gamma_l + n_k) / Gamma(gamma_l) \

  +& m_l log (1-epsilon.alt) + (N-m_l) log epsilon.alt / (K-1) \ \

m_l =& sum_(i=1)^N bold(1){x_(i l) = theta_(a_(i l) l)}
$

Note that variational update for $alpha, d_l, gamma_l$ is the same. \ \


== Variational update for $epsilon.alt$
$
log q^*(epsilon.alt)

prop& log p(X|C, gamma_l, epsilon.alt) + log p(epsilon.alt|lambda_1, lambda_2) \ \


prop& m_l log (1-epsilon.alt) + (N-m_l) log epsilon.alt / (K-1)

+ (lambda_1-1) log epsilon.alt + (lambda_2-1) log (1-epsilon.alt) \ \


prop& (lambda_1+N-m_l - 1) log epsilon.alt + (lambda_2+m_l - 1) log (1-epsilon.alt) \ \ \


epsilon.alt ~& "Beta"(lambda_1 + N - m_l, lambda_2 + m_l) \ \
$

#todo[interesting that despite the $1\/K-1$ factor, this is still a beta bernoulli update.]




#pagebreak()
= TODO

== initialization
- durbin (biological sequence analysis): positional bwt
- shapeit clustering

== hp tuning
- botorch with pareto (https://botorch.org/docs/multi_objective)

== max clade iou
- denom should be based on clustering of true marginal tree
- cut tree w/ same \# of clusters as dfcp, then use that normalization weighting

== dfcp w/ epsilon switch
- dfcp with mismatch = fixed constant epsilon likelihood instead of 0

== ARG model
- read The Effect of Single Recombination Events on Coalescent Tree Height and Shape
- need probabilistic model tying marginal trees together, w/ kingman's coalescent distributed marginal trees

== next possible
- eval internal repr of shapeit / beagle
- multimodal: adding likelihood terms

