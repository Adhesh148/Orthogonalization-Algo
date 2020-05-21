# Orthogonalization of Matrix-FS_183 6 
Plotting the loss of Orthogonality of a matrix at each iteration step due to four different methods of Orthogonalization 
This is a mini project in which we observe the numerical behaviour of various orthogonalization algorithms. We focus on the orthogonality of the computed vectors which may be lost in the classical or modified Gram-Shmidt algorithm, while the Gram-Schmidt algorithm with reorthogonalization has been shown to compute vectors which are orthogonal to machine precision level. 
All the algorithms have been implemented in the python notebook and the loss of orthogonality have also been plotted. 
###### Note:
All of these algorithms have been implemented and studied in the paper (Luc Giraud et al). This is a mere recreation of the algorithms in the paper using python for individual edification. Please refer to:
http://www.cerfacs.fr/algor/reports/2003/TR_PA_03_25.pdf

The matrix used in the python code is available in the matrix market under the name FS 183 6. The link for the same is provided here:
https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/facsimile/fs_183_6.html

### What is Orthogonalization
In simple terms, orthogonalization is the process of making two vectors orthogonal (perpendicular) to each other. Two vectors are said to be orthogonal to each other if their inner product is zero.
Refer to the following to get a better understading of orthogonal vectors:

### Why Orthogonalize ?

Orthogonality also makes a difference in how statistical tests are run. Orthogonal models only have one way to estimate model parameters and to run statistical tests. Non-orthogonal models have several ways to do this, which means that the results can be more complicated to interpret. In general, more correlation between independent variables means that you should interpret result more cautiously.

### What are the common methods of orthogonalization ?
* Classical Graham Scmidt (CGS)
* Modified Graham Scmidt (MGS)
* CGS with Reorthogonaliztion
* MGS with Reorthogonaliztion

### Which method works better ?
Keeping loss of orthogonality over each iteration as parameter to define the quality of an algorithm it is quite evident from the plots that the MGS method with reorthogonaliztion works the best (very closely followed by CGS with orthogonalization) where as the Classical Graham Scmidt suffers from great loss of orthogonality. We quite comfortably say that MGS with reorthogonaliztion is a better algorithm.
