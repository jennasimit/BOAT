# BOAT
BOAT (Bayesian Overlap Analysis Tool) identifies variants that are associated with two traits and tests for enrichment (i.e. whether there are more shared associated variants than expected by chance). Only summary statistics are required to implement BOAT.

Diseases often co-occur in individuals more often than expected by chance, and may be explained by shared underlying genetic etiology. When summary statistics are available, p-values are often used to assess association at a selected p-value threshold for both traits. However, p-values do not account for differences in power, whereas Bayes' factors (BFs) do, and may be approximated using summary statistics. Consequently, in overlap analyses, the use of BFs tend to result in a lower type I error rate than when p-values are used.

BOAT can be used to identify overlap variants between two traits by comparing ABFs, as well as by comparing p-values. It uses McNemar's mid-P test to assess overlap enrichment.

Asimit JL, Panoutsopoulou K, Wheeler E, Berndt SI, GIANT consortium et al. 2015. A Bayesian Approach to the Overlap Analysis of Epidemiologically Linked Traits. Genetic epidemiology. 39;8;624-34.
PUBMED: 26411566; PMC: 4832282; DOI: 10.1002/gepi.21919
