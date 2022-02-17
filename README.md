# epee
This R package wraps a small C++ library that uses expectation propagation to calculate a variational approximation to the likelihood of the phylogenetic threshold model. The advantage of the approximation is that it scales to a nearly arbitrary number of traits and tips.

For examples of use, see the functions `testcase_continuous()` (for continuous traits), `testcase_categorical()` (for binary traits), and `testcase_mixed()` (for mixed continuous-binary traits).

Slides decribing the method and implementation are found [here](https://github.com/nspope/epee/raw/master/inst/pope_ep_evolution17.pdf).
