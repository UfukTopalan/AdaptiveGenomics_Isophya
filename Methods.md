## STUDY

Our research investigates the genetic mechanisms underlying adaptations in Isophya rizeensis, a univoltine bush cricket species that inhabits a broad altitudinal range from sea level to 2500 meters. We aimed to understand how this species adapts to varying environmental conditions, particularly temperature and precipitation fluctuations associated with altitude.

To achieve this, we analyzed genome-wide polymorphisms in 71 individuals and 11 subpopulations of I. rizeensis using RAD sequencing data. Our study employed various analytical methods, including Principal Component Analysis (PCA), Admixture, and Snmf, Fst analyses, to explore genetic differentiation among populations. Additionally, Discriminant Analysis of Principal Components (DAPC) was used to identify specific genetic loci associated with altitudinal variations. Because we do not have a reference genome for this species or from any other closely related species we will use a denovo set of RAD contigs discovered using the CAYMY population of I. rizeensis as both an outgroup and a reference genome for alignment. 

We sought to elucidate the genetic diversity and adaptive responses of I. rizeensis to environmental changes. This research contributes to a deeper understanding of the evolutionary dynamics and adaptive mechanisms of species in response to altitude-related environmental factors, enriching the field of biodiversity research.

## Data Manipulation 

First, we need to create list files for both individuals and populations. You can use the following commands:

```bash
# Create a list of individuals
ls *.fastq.gz | cut -d'_' -f1-2 > isphya.list

# Extract unique population identifiers from the list
cut -c3-5 isphya.list | uniq > pop11_list
```
## Alignment with BWA
### Alignment and Processing
In our analysis, we utilize a series of tools to align sequencing reads and process the resulting BAM files:

- BWA-MEM: The BWA-MEM algorithm is employed to align sequencing reads to a reference genome. It is well-suited for handling long reads and provides accurate alignments by performing local alignment. This algorithm uses a seed-and-extend approach, which enhances its sensitivity to correctly align reads that span across large regions of the genome.

- Samtools: Post-alignment, Samtools is used to process the BAM files. It performs tasks such as sorting the aligned reads and filtering to retain only properly paired reads. This ensures that the data is of high quality and ready for further analysis.

- Picard: The Picard tool is used to mark duplicate reads in the BAM files. Identifying duplicates is important for reducing biases that may arise during sequencing and ensuring the accuracy of subsequent analyses.

Before running the alignment script, we must index our reference genome.
```bash
~/bin/bwa/bwa index -a bwtsw ~/isophya/references/isophya_contigs_CAYMY.fasta
```

We use this [script](scripts_folder/alignment.sh) for aligning sequencing reads to a reference genome using BWA and processing the resulting BAM files with Samtools and Picard. To run the script, use the following command:

```bash
sbatch align_pe_reads.sh ~/isophya/raw/isphya.list ~/isophya/references/isophya_contigs_CAYMY.fasta
```
Now that we have 71 BAM files at low/medium depth, we can create a BAM list for the dataset as a whole, and also for each subpopulation to use in downstream analysis.

```bash
ls ~/isophya/*_sorted_mdup.bam > isophya71.bamlist
for i in cat pop11_list; do grep $i isophya71.bamlist > ${i}.bamlist; done
```

For understanding the samtools flags, you can refer to their official documentation at [Samtools Format](https://www.samformat.info/sam-format-flag#google_vignette) ; Let’s check one of the files using the `samtools flagstat` command:
```bash
samtools flagstat CKLYU_010_sorted_mdup.bam
24236 + 0 in total (QC-passed reads + QC-failed reads)
4412 + 0 duplicates
15408 + 0 mapped (63.57%:-nan%)
24236 + 0 paired in sequencing
12080 + 0 read1
12156 + 0 read2
13148 + 0 properly paired (54.25%:-nan%)
14714 + 0 with itself and mate mapped
694 + 0 singletons (2.86%:-nan%)
1554 + 0 with mate mapped to a different chr
826 + 0 with mate mapped to a different chr (mapQ>=5)
```

Now, to check the stats for all BAM files, you can use a simple script available [here](scripts_folder/count_no_align.sh).
Run the script with the following command:
```bash
sbatch count_no_align.sh
```
With the output, you can review the mapped reads in the third column and filter out any individuals or reads with low quality.


## Eliminating Paralogs

After individual elimination, we can perform filtering for paralog sites. Before doing so, we need to create a position file that contains the positions for all populations. This can be achieved by genotyping all subpopulations and creating a `.mafs` file, then extracting the positions from it.

To genotype all groups, you can use the following [script](scripts_folder/genotyping_groups.sh) You'll need to provide a population list file, which contains all subpopulations, and create a `.bamlist` file for each subpopulation. You will also need to provide a reference FASTA file.

To run the script:

```bash
sbatch genotyping_groups.sh pop11_list ~/references/isophya_contigs_CAYMY.fasta
```
After creating `.mafs` files for each subpopulation, extract the positions and create a `.pos` file with the following code:

```bash
zless CANCK.mafs.gz | cut -f1-2 | sed 1d > CANCK.pos
```
Next, we can filter the paralog sites using the script [here](scripts_folder/get_paralog.sh). Note that the code provided is an example, and you'll need to organize it for your use case. In this example, I have performed it separately for each subpopulation.

To run the script:
```bash
#pop11 list contains my subpopulations names with .bamlist extensions and in that file they have their bam files paths.
sbatch get_paralog.sh pop11_list
```

This process will generate a file with a .paralogs extension. By conducting a Bonferroni correction (0.05 confidence level) on our SNP number, we can compute a chi-square value using [this website](http://courses.atlas.illinois.edu/spring2016/STAT/STAT200/pchisq.html) and set the threshold based on that (df = 1). Paralogs are eliminated if their values are higher than this threshold. 

For example for the CANCK subpopulation:
```bash
less CANCK.paralogs | awk '$5>35.51  {print $1, $2}'> CANCK_paralog.list
awk 'NR==FNR{a[$1,$2]; next} !(($1,$2) in a)' CANCK_paralog.list CANCK.pos > CANCK_sites
```
Now that we have created a sites file for each population, excluding any paralog sites, we can proceed to create a new sites file. In this new file, we will extract the sites present across all populations. Below is an example of how to do this:
```bash
awk 'NR==FNR{a[$0];next} $0 in a' CANCK_sites PLVYL_sites > CANCK_PLVYL.sites
```
To extract the sites across all groups, use a simple [loop](scripts_folder/extract_common_sites.sh):

To run the loop:
```bash
sbatch extract_common_sites.sh
```

We now have a `nonparalog_sites` file.
At this point, we have filtered out paralog sites, but the process isn't finished. We still need to perform further filtering on the data.

### Filtering Only High Coverage Sites
Now, for obtaining and continue only with high coverage sites which this time we'll use 6x coverage, we can use the same pop11_list and create a mafs file for each of the populations as above but this time with this [script](scripts_folder/high_coverage.sh).

To run the script:
```bash
sbatch high_coverage.sh pop11_list
```

Now we're doing the same work as above, extracting the sites from mafs file, and than running the loop (don't forget to update it for a different output) for obtaining the sites that are present across all populations. Then we must do one more step, we will continue with the sites that appear on both 6x coverage sites file and non_paralogs sites files. 

For that, let's crerate our final sites file then also index it for genotyping.

```bash
awk 'NR==FNR{a[$0];next} $0 in a' 6x_coverage_sites nonparalog_sites > 6cov_nonparalog.sites
```

Now for genotyping and variant calling proccess, we have everything that we need.

## Genotyping and Variant Calling
Now, variant calling and genotyping are performed using ANGSD, a toolkit designed for analyzing next-generation sequencing data. The process involves identifying genetic variants such as single nucleotide polymorphisms (SNPs) and determining the genotype of individuals at each site across the genome.

ANGSD allows for flexible analysis by calculating genotype likelihoods rather than directly calling genotypes, which is ideal for working with low to medium-depth sequencing data. This method accounts for uncertainty in base calling and sequence quality, providing more robust results when sequencing coverage is limited. The toolkit also calculates allele frequencies across populations, helping to identify polymorphic sites.

Several output formats, such as Beagle, VCF, and PLINK, are supported to enable further analysis using a variety of downstream tools. This flexibility allows for comprehensive population genetics analyses, ensuring accurate and reliable detection of genetic variants within the dataset.

We'll use several filtering options for obtaining only high quality data. To understand and apply the filtering options for your data, you can check [angsd filter section](https://www.popgen.dk/angsd/index.php/Filters).

Now, for calling genotypes on our bamfiles, we must give a bamlist which we've created earlier. Also again, we need to specify the references sequence too. 
The options used for variant calling in ANGSD include `-GL 1` for genotype likelihoods using the SAMtools model, `-doMajorMinor 1` for inferring the major and minor alleles, `-doMaf 1` to calculate minor allele frequencies, `-doGlf 2` for outputting genotype likelihoods in Beagle format, `-doGeno 5` for generating genotype probabilities, `-doBcf 1` to create a BCF file, `-doPost 1` for posterior probability calculations with a cutoff of `-postCutoff 0.95`, filtering by mapping quality (`-minMapQ 10`), base quality (`-minQ 20`), a minimum number of individuals (`-minInd $mInd`), SNP p-value threshold (`-SNP_pval 1e-12`), minor allele frequency (`-minMaf 0.05`), and limiting the analysis to sites from the file `group_coverage/6cov_nonparalog.sites` which we obtain earlier by removing the paralogs and including only has 6x coverage sites.
We can use this [script](scripts_folder/genotype.sh) to obtain all the files above.
To run this code:
```bash
sbatch genotype.sh isophya71.bamlist ~/references/isophya_contigs_CAYMY.fasta
```
We created the beagle file, mafs file, bcf file, geno file and tped file which we'll use all this file in our analysis.

## Population structure (PCA)

To conduct Principal Component Analysis (PCA) with `pcangsd`, we work directly with genotype likelihoods in Beagle format, without calling genotypes. This approach accounts for genotype uncertainty, making it especially suitable for low-depth sequencing data. The command specifies a Beagle file (`${pop}.beagle.gz`) as input and saves the PCA results in the `results_pca` directory. `pcangsd` efficiently computes covariance matrices based on genotype likelihoods and outputs principal components, which can be used to explore genetic structure among populations.
We can use this [script](scripts_folder/get_PCA.sh).
To run it
```bash
sbatch get_PCA.sh isophya71
```
Running the `get_PCA.sh` script generates an output file called `isophya71.cov`, which contains the covariance matrix. This matrix will be used for eigendecomposition to identify and visualize the main axes of genetic variation. The first two principal components (PC axes) can be plotted using an R script [here](scripts_folder/plotPCA.R), and a simple cluster file can be created to label populations with distinct shapes and for dark and pale population using two colors for clear visualization in the plot.

Our `isophya.firtina.valley.tsv` file contains information about individuals' colors, altitudes, body size, and temperature seasonality. We can use this file to create a simple `.clst` file for labeling populations.

**Example of the TSV file:**

| Indv       | Color | Alt   | Bsize | TS   |
|------------|-------|-------|-------|------|
| CANCK_007  | Pale  | 1200  | 0.66  | 77.8 |
| CANCK_00D  | Pale  | 1200  | 0.66  | 77.8 |
| CANCK_010  | Pale  | 1200  | 0.66  | 77.8 |
| CANCK_011  | Pale  | 1200  | 0.66  | 77.8 |

**Creating the `.clst` file:**

```bash
cut -f1 isophya.firtina.valley.tsv | sed 1d | sed 1iFID > FID
cut -f1 isophya.firtina.valley.tsv | sed 1d | cut -c1-5 | sed 1iIID > IID
less isophya/isophya.firtina.valley.tsv | grep -f taken.list | sort -nk3 | cut -f2 > CLUSTER  # Include only selected individuals
paste -d' ' FID IID CLUSTER > isophya71.clst

# To simplify the plot, we can replace IID with numbers for subpopulation labeling:
sed -i 's/ IST06/ 1/g' isophya.firtina71.clst
...
```
Running the R script for PCA plot:
```bash
Rscript plot.PCA.R -i isophya71.cov -c1-2 -a isophya71.clst -o isophya71_pca.pdf
```



