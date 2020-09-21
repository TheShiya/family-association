# Family-based association study design

An empirical study of several flavors of family-based vs population-based design to examine each of their advantages.

## Overview of project goals

We use simulation studies to examine family-based vs population-based design. The gain in power of analyzing families comes from the fact that the effect variants should be enriched within the family. 
However, since family members are correlated the effective sample size is smaller. Also controls in families are likely to be carriers of disease alleles which can hurt power, compared to using population controls. 

We therefore look into different strategies to assess the trade-off:

1. Analyze families using LMM including both affected and unaffected individuals. We can either analyze all families, or analyze them separately and combine the outcome later.
2. Analyze unrelated case control design of the same sample size and phenotypic model, and see how much larger the effect size have to be in families (in step 1), to offset the loss of power due to using related individuals, compared to unrelated design.
3. Analyze affected individuals in the families but replaced the unaffected members with related controls from other families of the same structure, where no one is affected in their family. That is, the controls are still related, but perhaps no longer carrying a risk allele.
4. Analyze affected family members and replaced the unaffected members with population-based controls. Here the controls are unrelated.
5. Analyze a sample of unrelated individuals with the same number of affected and unaffected individuals as in analysis 1-3.
