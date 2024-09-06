# STUDY

Our research investigates the genetic mechanisms underlying adaptations in Isophya rizeensis, a univoltine bush cricket species that inhabits a broad altitudinal range from sea level to 2500 meters. We aimed to understand how this species adapts to varying environmental conditions, particularly temperature and precipitation fluctuations associated with altitude.

To achieve this, we analyzed genome-wide polymorphisms in 71 individuals and 11 subpopulations of I. rizeensis using RAD sequencing data. Our study employed various analytical methods, including Principal Component Analysis (PCA), Admixture, and Snmf analyses, to explore genetic differentiation among populations. Additionally, Discriminant Analysis of Principal Components (DAPC) was used to identify specific genetic loci associated with altitudinal variations.

We sought to elucidate the genetic diversity and adaptive responses of I. rizeensis to environmental changes. This research contributes to a deeper understanding of the evolutionary dynamics and adaptive mechanisms of species in response to altitude-related environmental factors, enriching the field of biodiversity research.

## Data Manipulation 

First of all, we need to create our list files for individuals and populations. You can use the following commands:

```bash
# Create a list of individuals
ls *.fastq.gz | cut -d'_' -f1-2 > isphya.list

# Extract unique population identifiers from the list
cut -c3-5 isphya.list | uniq > pop.list

