# retrodesign
[![Build Status](https://travis-ci.org/andytimm/retrodesign.svg?branch=master)](https://travis-ci.org/andytimm/retrodesign)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/retrodesign)](https://cran.r-project.org/package=retrodesign)

Tools for calculating and working with Gelman et al's Type S (sign), and M (magnitude or exaggeration ratio) errors for analyzing hypothesis tests in R.

Features:

* Functions for calculating Power and Type S/M errors across variety of effect sizes, building on code provided in Gelman and Carlin (2014, see below).

* Graphics function for visualizing these errors with or without ggplot

* Implementation of Lu et al's (2018, see below) closed form solution for Type M error, providing a speed-up.

# Usage

You can install retrodesign with:

```
install.packages("retrodesign")
```

You can find an online version of retrodesign's vignette [on my website](https://andytimm.github.io/2019/02/05/Intro_To_retrodesign.html), which provides an introduction to both `retrodesign` and Type S/M errors.

# More Reading on Type S/M Errors

1. **Gelman and Tuerlinckx's [Type S error rates for classical and Bayesian single and multiple comparisons
procedures](http://www.stat.columbia.edu/~gelman/research/published/francis8.pdf) (2000)**: A comparison of the properties of Type S errors of frequentist and Bayesian confidence statements. Useful for how this all plays out in a Bayesian context. Bayesian confidence statements have the desirable property of being more conservative than frequentist ones.
2. **Gelman and Carlin's [Assessing Type S and Type M Errors](http://www.stat.columbia.edu/~gelman/research/published/retropower20.pdf)
 (2014)**: Gelman and Carlin compare their suggested design analysis, as we've written about above,
 to more traditional design analysis, through several examples, and discuss the
 desirable properties it has in more depth than I do here. It is also the source of the original
 retrodesign() function, which I re-use in the package with permission.
3. **Lu et al's [A note on Type S/M errors in hypothesis testing](https://onlinelibrary.wiley.com/doi/full/10.1111/bmsp.12132)
(2018)**: Lu and coauthors go further into the mathematical properties of Type S/M errors, and prove the closed form solutions implemented in `retrodesign`.
4. **McShane et al's [Abandon Statistical Significance](https://arxiv.org/abs/1709.07588) (2017)**: If you
want a starting point on the challenges with NHST that have led many statisticians
to argue for abandoning NHST all together, and starting points for alternative
ways of doing science.
