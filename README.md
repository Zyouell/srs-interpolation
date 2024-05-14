# srs-interpolation

## Background

So you have a few lists of $d= 2^{n}$ evaluations for polynomials of degree $d - 1$ defined over some finite field $\mathbb{F}_{r}$ with 2-adicity at least $d$ that also happens to be the scalar field of a pairing friendly elliptic curve. You wish to commit these polynomial using the [KZG](https://www.iacr.org/archive/asiacrypt2010/6477178/6477178.pdf) polynomial commitment scheme but you don't want to perform an inverse fast Fourier transform each time. Well you're in luck! You can perform an FFT-like operation just once on your KZG structured reference string. This will produce a list of commitments to the Lagrange basis of size $d$, then you can simply perform an MSM on each of your evaluation lists to recover the KZG commitment.

The code in this repository is based off of techniques described in [this](https://eprint.iacr.org/2017/602.pdf) paper.

## High-Level Overview

The SRS of the KZG commitment scheme is generated using a trusted setup to generate powers of a scalar $\tau\in \mathbb{F}_{r}$ and multiply the generator point $g \in \mathrm{E}(\mathbb{F}_{q})$ by these powers. This results in a list of points $$\begin{align}
g, \tau\cdot g, \tau^{2}\cdot g, \dots, \tau^{d-1}\cdot g.
\end{align}$$

The Lagrange basis polynomials if the evaluation points are primitive-$2^{n}$ ($\omega\in\mathbb{F}_{r}$ with $\omega^{2^{n}} = 1$) roots of unity are defined as

$$
\begin{align} w_{j} =& \ \prod_{j\neq m}(\omega^{j} - \omega^{m})^{-1} \\ l(x) =& \ \prod_{i = 1}^{d} (x - \omega^{i}) \\
\mathrm{L}_{i}(x) =& \ l(x)\cdot \frac{w_{i}}{x-\omega^{i}}. \end{align}
$$

From this we see that $\mathrm{L}_{i}(\omega^{i})=1$ and $\mathrm{L}_{i}(\omega^{j})=0$ if $i \neq j$. Naively we could produce commitments to the $\mathrm{L}_{i}(x)$ by multiplying everything out and calculating coefficients, but this is horribly inefficient. Instead consider the polynomial
$$ P(x, Y) := 1 + x\cdot Y + x^{2}\cdot Y^{2} + \dots + x^{d-1} \cdot Y^{d-1}. $$
By substituting $Y = \omega^{-i}$ the above becomes
$$ P(x, \omega^{-i}) = 1 + \left ( \frac{x}{\omega^{i}}\right) + \left ( \frac{x}{\omega^{i}}\right)^{2}+\dots+\left ( \frac{x}{\omega^{i}}\right)^{d-1} $$
which has the familiar property $P(\omega^{i}, \omega^{-i}) = 1$ and $P(\omega^{j}, \omega^{-i}) = 0$ if $i\neq j$, hence $P(x,\omega^{-i}) = \mathrm{L}_{i}(x)$ (by the uniqueness of Lagrange interpolation).

Now we can use the normal FFT trick of splitting $P$ into its even and odd coefficients $$\begin{align*}P_{0}(x, \omega^{-i}) =& \sum_{j=0}^{d/2 -1} \left ( \frac{1}{\omega^{i}}\right)^{2j}\cdot x^{j} \\ P_{1}(x, \omega^{-i}) =& \sum_{j=0}^{d/2 -1} \left ( \frac{1}{\omega^{i}}\right)^{2j + 1}\cdot x^{j} \\ P(x, \omega^{-i}) =& P_{0}(x^{2}, \omega^{-i}) + x\cdot P_{1}(x^{2}, \omega^{-i})\end{align*}$$
and repeating this procedure recursively until $P_{0}$ and $P_{1}$ are both constant polynomials. We then recombine them using the butterfly trick.

### A Simple Example

We have the polynomial $$A(x) = A_{0} + A_{1}\cdot x + A_{2}\cdot x^{2} + A_{3}\cdot x^{3}$$
and wish to find its evaluations on the set $\{1, \omega, \omega^{2}, \omega^{3}\}$ where $\omega^{4}=1$. First we split $A$ into lists of odd and even coefficients $$\begin{align*} A_{\mathrm{even}}:=& \ \{A_{0}, A_{2}\}\\ A_{\mathrm{odd}}:=& \ \{A_{1}, A_{3}\} \end{align*}$$
and now we perform the same FFT procedure on $A_{\mathrm{even}}$ and $A_{\mathrm{odd}}$ to get $Y_{\mathrm{even}}$ and $Y_{\mathrm{odd}}$, the evaluations on the set $\{1, \omega^{2}\}$ (since $\omega^{4}=1$ we have $\omega^{2}= -1$). Observe $$\begin{align*} Y_{\mathrm{even}} =& \ \{A_{0}+A_{2}, A_{0}-A_{2}\} \\ Y_{\mathrm{odd}} =& \ \{A_{1}+A_{3},A_{1}-A_{3}\}\end{align*}$$
thus $$\begin{align*} Y=\{&A_{0}+A_{1}+A_{2}+A_{3},\\&A_{0}+\omega\cdot A_{1} - A_{2} -\omega\cdot A_{3} ,\\ &A_{0}-A_{1}+A_{2}-A_{3},\\ &A_{0}-\omega\cdot A_{1} - A_{2}+\omega\cdot A_{3}\} \end{align*}$$
which are precisely the evaluations of $A(x)$ at $\{1, \omega, \omega^{2}, \omega^{3}\}$.
