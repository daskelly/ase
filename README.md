# Allele-specific expression

This repository contains code to implement the model of allele-specific expression
(ASE) described in Skelly et al. 2011 Genome Research 
[dx.doi.org/10.1101/gr.119784.110](link).

The directory [biased_SNPs](biased_SNPs) provides code to implement the model
for detecting biased SNPs from DNA data. These SNPs are then filtered out for all
subsequent analyses. This model is described in section S1.3.1 of the 
supplementary material. See the README [here](biased_SNPs/README.md).

The directory [DNA_model](DNA_model) provides code to implement the model
for genomic DNA read counts that estimates overdispersion in this "null" data where 
no genes should show ASE. This model is described in section 1.3.2 of the 
supplementary material. See the README [here](DNA_model/README.md).

The directory [RNA_model](RNA_model) provides code to implement the model
to detect ASE in read counts derived from RNA. This model is described in section 1.3.3 
of the supplementary material. See the README [here](RNA_model/README.md).
