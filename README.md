### glmmSeq: Taxonomic sequence classification with generalized linear mixed models

Contact: Rodney Richardson -- rodney.richardson {at} umces.edu

### License
GNU General Public License v3.0

### Dependencies

Unix/Linux/OSX  
Python 3.6.8  
vsearch 2.14.1  
R 3.6.1  
MuMIn R package  
glmmTMB R package  

### What glmmSeq Does
glmmSeq is designed for the classification of DNA sequences produced in metagenetic and metagenomic applications. The software can be trained to classify sequences from any genetic locus and cross-validation results indicate that it is highly generalizable across markers with variable genetic characteristics (e.g. sequences from plant ITS2, plant rbcL and arthropod COI). The current version of the code is undergoing beta testing.

### Motivation
A number of methods exist for the taxonomic classification of DNA sequences, such as those produced in metabarcoding, metagenetic and metagomic applications. Accurate estimation of classification confidence is among the most desirable features of any such classification procedure. Unfortunately, classifiers have typically been designed for specific use cases, such as classification of microbial 16S data, and with little follow-up regarding the generalization of accurate classification performance to other systems and markers (e.g. arthropod COI or plant *rbcL*, *trnL* and ITS2). For example, SINTAX (https://www.biorxiv.org/content/10.1101/074161v1.abstract), one of the highest performing classifiers (https://peerj.com/articles/4652/), appears to generalize poorly to some commonly used plant metagenetic markers (Figure 1), with classification confidence values that strongly underestimate actual confidence and are thus difficult for end users to evaluate and use. Ideally, accurate classification confidence estimates would closely correspond to the probability of correct classification (Figure 2).

**Figure 1:** 10-fold family-level cross-validation outcomes (1=misclassification, 0=correct classification) regressed against SINTAX classificaition confidence values for plant *rbcL* data. A LOESS best-fit regression line is shown in red and an ideal relationship between confidence values and classification outcomes is shown with a black line. 
  
<img src="https://github.com/RTRichar/glmmSeq/blob/master/Figures/Sintax.png" width="350" height="300">

**Figure 2:** 10-fold family-level cross-validation outcomes (1=misclassification, 0=correct classification) regressed against SINTAX classificaition confidence values for plant *rbcL* data. A LOESS best-fit regression line is shown in red and an ideal relationship between confidence values and classification outcomes is shown with a black line. 
  
<img src="https://github.com/RTRichar/glmmSeq/blob/master/Figures/glmmSeq.png" width="350" height="300">

### How glmmSeq Works
With glmmSeq, classifications are performed based on glmmTMB binomial generalized linear mixed models (https://www.biorxiv.org/content/10.1101/132753v1), a method which has not been considered in previous work to our knowledge. For these models, the categorical taxonomic identities of alignment subject sequences are specified as random intercept terms to account for random variance in the relationship between the log-odds of classification success and alignment percent identity (Figure 3). Briefly, glmmSeq can be trained on reference sequences from any marker or taxonomic study system. During training, simple alignment-based classifications are performed using pairwise VSEARCH (https://peerj.com/articles/2584/) alignment, yielding classification outcomes (0 for correct classification, 1 for incorrect classification) and alignment percent identities. This alignment-based classification procedure is conducted within a k-fold cross-validation framework, and classification outcomes are then regressed against alignment percent identity values with additional optional fixed and random effect terms. For classifications at the family, genus and species levels, the default classification model is shown below and this method typically results in highly accurate classification confidence estimates (Figure 2).
  
Classification Outcome ~ Alignment Percent Identity + (1|Subject Class/Subject Order/Subject Family)
  
**Figure 3:** Illustration of the best-fit relationship between the probability of misclassification and pairwise alignment percent identity for a sample of arthropod COI sequences. The black line illustrates the global relationship for the data, while red, green and blue lines are fit for specific taxonomic groups to illustrate the considerable variance in these relationships.

<img src="https://github.com/RTRichar/glmmSeq/blob/master/Figures/MotivationalPlot.png" width="400" height="300">

### Getting started

TRAIN_glmmSeq.py --help  
CLASSIFY_glmmSeq.py --help
