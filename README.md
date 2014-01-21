Bayesian Entropy estimator for binary vector observations
=========================================================

Estimating Shannon's entropy from data is difficult, especially when you have less data compared to all possible symbols.
However, if your data are in the form of binary vectors, you're in luck!
BDentropy provides two state-of-the-art Bayesian entropy estimators for binary vectors.
It was primiarily developed for estimating entropy of neural spike trains, however, if your data has a similar structure that our prior assumes, it could work very well for you too.

The MATLAB code is a reference implementation for the results described in the following paper:

- Evan Archer, Il Memming Park, Jonathan W. Pillow. [Bayesian entropy estimation for binary spike train data using parametric prior knowledge. Neural Information Processing Systems](http://papers.nips.cc/paper/4873-bayesian-entropy-estimation-for-binary-spike-train-data-using-parametric-prior-knowledge) [(NIPS) 2013](http://books.nips.cc/nips26.html)

Installation
------------
Downloading the source code:

    $ git clone 
    $ git submodule update --init # to pull the PYMentropy submodule

To compiling the MEX code, run the following in `src` folder.

    >> makeMex

Use the `startup.m` script in `src` to add the relevant paths to MATLAB.

Quick example
-------------
Let's say that your observation is in the following matrix form:

![](doc/figs/binary_vector_observations.png)

Here, each column of your matrix corresponds to an `m` dimensional observation.
To estimate entropy

    >> words = [...
	[0 1 0 0 0 0 1]', ...
	[1 0 0 0 0 0 1]', ...
	[0 0 0 0 0 1 0]', ...
	[0 1 0 0 0 0 1]', ...
	[0 0 0 0 0 1 0]', ...
    ];
    >> [nn ocnts] = words2nnOcnts(words); % compact representation
    >> H = computeH_BD(nn,ocnts,size(words,1)) % estimate entropy

    H =
	2.6644

`computeH_BD` returns the estimated entropy in unit of bits.

Tutorial
========
The code corresponding to the following tutorial is in [tutorial.m](src/tutorial.m)

Compact representation
----------------------
To estimate the entropy, we only need to know the number of unique words with TODO

    >> [nn ocnts] = words2nnOcnts(words); % compact representation

Variance
--------
If more than one return value are requested to `computeH_BD`, it samples from the posterior and returns the variance, confidence intervals, and samples. Due to the sampling, it is significantly slower than just computing the mean.

TODO: put the syntax and examples


External links
==============
A closely related sister-paper also appeared at the same conference:

* Il Memming Park, Evan W. Archer, Kenneth Latimer, Jonathan W. Pillow. [Universal models for binary spike patterns using centered Dirichlet processes](http://papers.nips.cc/paper/5050-universal-models-for-binary-spike-patterns-using-centered-dirichlet-processes) NIPS 2013

This code shares some code from the [PYMentropy](https://github.com/pillowlab/PYMentropy) project. The PYM entropy estimator is a generic discrete entropy estimator, not restricted to binary vector obsevations.
