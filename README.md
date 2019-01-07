# Reduced Word Enumeration, Complexity, and Randomization
Cara Monical, Benjamin Pankow, and Alexander Yong

Usage Instructions
----------

This repository contains three estimators for counting reduced words:
* `hecke_estimator.py` estimates the number of Hecke words of a given permutation `P` and length `e` using the random descent algorithm
* `transition_estimator.py` estimates the number of reduced words of a permutation `P` using the transition algorithm
* `combined_hecke_estimator.py` estimates the number of Hecke words of a permutation `P` and length `e` using both the random descent algorithm and transition algorithm

To use, execute the chosen estimator using Python 3: `python3 hecke_estimator.py`. Example usage is as follows:

```
Permutation (comma separated): 5,4,3,2,1
e: 7 
Number of trials T: 10000
Executed in 2 seconds

Estimated number of Hecke words:
3.31863 x 10^5
```
