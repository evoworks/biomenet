
## BiomeNET

This is a repository for a forked version of our code that implements **a novel hierarchical Bayesian model, called BiomeNET (Bayesian inference of metabolic networks), for inferring differential prevalence of metabolic subnetworks among microbial communities**. To infer the structure of community-level metabolic interactions, BiomeNET applies a mixed-membership modelling framework to enzyme abundance information. The basic idea is that the mixture components of the model (metabolic reactions, subnetworks, and networks) are shared across all groups (microbiome samples), but the mixture proportions vary from group to group. Through this framework, the model can capture nested structures within the data. BiomeNET is unique in modeling each metagenome sample as a mixture of complex metabolic systems (metabosystems). The metabosystems are composed of mixtures of tightly connected metabolic subnetworks. BiomeNET differs from other unsupervised methods by allowing researchers to discriminate groups of samples through the metabolic patterns it discovers in the data, and by providing a framework for interpreting them. We employ a collapsed Gibbs sampler for inference of the mixture weights under BiomeNET, and we use simulation to validate the inference algorithm.


The paper describing the method is:

>Shafiei, M., Dunn, K. A., Chipman, H., Gu, H., & Bielawski, J. P. (2014). BiomeNet: a Bayesian model for inference of metabolic divergence among microbial communities. PLoS computational biology, 10(11), e1003918.

Please cite this paper if you use this repository or the original code.

The original code was written by Mahdi Shafiei (Bielawski Group postdoc). The repository for the original version is here: https://sourceforge.net/projects/biomenet/
