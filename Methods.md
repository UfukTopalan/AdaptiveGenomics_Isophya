##Data Ordering 

First of all, we need to create our list files for individuals and populations. You can use the following commands:

```bash
# Create a list of individuals
ls *.fastq.gz | cut -d'_' -f1-2 > isphya.list

# Extract unique population identifiers from the list
cut -c3-5 isphya.list | uniq > pop.list
