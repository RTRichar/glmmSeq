# glmmSeq
The taxonomic classification of DNA sequences is an important component of most metabarcoding and metagenomic applications. Various classifiers have been put forward to satisfy this need, however, none have been shown to be both probabilistically accurate and robust to gaps in reference sequence data (i.e. the model overfitting problem). For example, the sintax classifier routinely scores among the top of the field in terms of accuracy, yet even this classifier exhibits remarkable misclassification and overclassifcation rates (). 


Sequence classification with generalized linear mixed models

GLMM clearly applicable to issue of sequence classification
problem of overfitting with other classifiers

desirable framework with respect to limitting the potential of overfitting and explicitely modelling log odds of classification sucess

CV binomial outcomes ~ predicted confidence
test robustness to overfitting (acc ~ db size)
compare to SINTAX
