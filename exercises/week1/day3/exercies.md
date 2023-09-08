# Questions:
## PART 1
- How large (in base pairs) is the bacterial reference genome?
```cat mleprae_reference_genome.fasta.gz  | tail -n +2 | wc -m```
Ther result is 3,226,951 base pairs.

- Use bwa aln (and samse) to align the two fastq files to the bacterial reference (Hint: look at the exercises from https://github.com/ANGSD/adv_binf_2023_week1/tree/main/day2).
Align first 
```bwa aln mleprae_reference_genome.fasta.gz l10.fastq.gz > l10_SA.sai```
```bwa aln mleprae_reference_genome.fasta.gz l30.fastq.gz > l30_SA.sai```
Produce the SAM files: 
``` bwa samse mleprae_reference_genome.fasta.gz l10_SA.sai l10.fastq.gz > output_l10_aln.sam ```
``` bwa samse mleprae_reference_genome.fasta.gz l30_SA.sai l30.fastq.gz > output_l30_aln.sam ```

- Sort and index the resulting files using samtools, and make sure they are saved in bam format.
Sort first and output the .bam files withouth any additional arguments so that it is compatible with the index command. 
```samtools sort output_l10_aln.sam -o output_sorted_l10_aln.bam```
```samtools sort output_l30_aln.sam -o output_sorted_l30_aln.bam```
check: samtools quickcheck -vvv output_sorted_l10_aln.bam 
Then, indexing the two binary files: 
```samtools index output_sorted_l10_aln.bam```
```samtools index output_sorted_l30_aln.bam```

- Filter out the unaligned reads and create new bam files that contain only mapped reads. Identify which flag to filter out on (https://broadinstitute.github.io/picard/explain-flags.html)
Filtering the two sorted binary files with Flag 4: 
```samtools view -F 4 output_sorted_l10_aln.bam -o output_sorted_filtered_l10_aln.bam```
```samtools view -F 4 output_sorted_l30_aln.bam -o output_sorted_filtered_l30_aln.bam```

From the files we can see that l10 aligned map does not contain any unmapped whereas the l30 contains 1 unmapped read.

- By looking at the reads aligning and the chromosomal positions, we can calculate how many reads are mapped correctly with a 0 nucleotide difference between the origin in the read ID and the mapping coordinate. The script /TEACHING/BIOINF22/adv_binf_2023_week1/day3/get_stats.sh takes a bam file and outputs two numbers - firstly the number of reads that map correctly (i.e. those where the start position in the bam file matches the read ID), and then the number that do not map correctly.

Use this script to find out how many reads map correctly and incorrectly in the two bam files

copy first 
```cp /TEACHING/BIOINF23/adv_binf_2023_week1/day3/get_stats.sh .```

Run the command: 
```./get_stats.sh output_sorted_filtered_l10_aln.bam```
It has 15171 correctly mapped and 84829 incorrectly mapped 


```./get_stats.sh output_sorted_filtered_l30_aln.bam ```
It has 99699 correctly mapped and 300 incorrectly mapped


- Now create two new bam files where you filter the reads so we only retain reads with a mapping quality of greater than or equal to 1. (Hint: look at the exercises from https://github.com/ANGSD/adv_binf_2023_week1/tree/main/day2).
Repeat the above steps exercise with the two new filtered bam files. How do the numbers differ?

```samtools view -F 4 -q 1 output_sorted_l10_aln.bam -o output_sorted_filtered_map_qual_l10_aln.bam```
```./get_stats.sh output_sorted_filtered_map_qual_l10_aln.bam```
Not it has 1681 correctly mapped and 53 incorrectly mapped 

```samtools view -F 4 -q 1 output_sorted_l30_aln.bam -o output_sorted_filtered_map_qual_l30_aln.bam```
```./get_stats.sh output_sorted_filtered_map_qual_l30_aln.bam```
All of the reads are correctly mapped (99453)


## PART 2 - Data analysis:
In this part you are given two compressed FASTQ files in the directory /TEACHING/BIOINF23/day3/

One of these files is ancient, and one is modern. They are both reads from the pathogen Pseudomonas aeruginosa (reference genome fasta: paeruginosa.fasta.gz)

jessica.fastq.gz
fletcher.fastq.gz 
The goal is to ascertain which of the files (Jessica or Fletcher) belongs to the modern or ancient sample. For this exercise, you can look at the computer exercises from day1 and day2 for inspirations (https://github.com/ANGSD/adv_binf_2023_week1/tree/main/day1 and https://github.com/ANGSD/adv_binf_2023_week1/tree/main/day2)

Questions:
- Perform, using fastp, adapter trimming on the provided sample sequence files. Make sure to discard reads under 30 bp (-l parameter).
First, I would like to know the average length of the reads of each fasta file 
`zcat file | awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}'`
It turns out that Jessica and Fletcher have 87.8 and 125 average length of reads respectively. 
Next we want to trim the adapter of those two files.
``` fastp -l 30 --detect_adapter_for_pe -i jessica.fastq.gz -o jessica_trimmed_filtered.fastq.gz```
Detecting adapter sequence for read1...
>Illumina TruSeq Adapter Read 1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

Read1 before filtering:
total reads: 800000
total bases: 70223027
Q20 bases: 67578862(96.2346%)
Q30 bases: 65815231(93.7231%)

Read1 after filtering:
total reads: 477112
total bases: 28207281
Q20 bases: 27296289(96.7704%)
Q30 bases: 26679081(94.5822%)

Filtering result:
reads passed filter: 477112
reads failed due to low quality: 2056
reads failed due to too many N: 634
reads failed due to too short: 320198
reads with adapter trimmed: 705878
bases trimmed due to adapters: 34460849

Duplication rate (may be overestimated since this is SE data): 0.714589%



```fastp -l 30 --detect_adapter_for_pe -i fletcher.fastq.gz -o fletcher_trimmed_filtered.fastq.gz```
>Illumina TruSeq Adapter Read 1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

Read1 before filtering:
total reads: 800000
total bases: 100000000
Q20 bases: 95946778(95.9468%)
Q30 bases: 92582372(92.5824%)

Read1 after filtering:
total reads: 800000
total bases: 71977242
Q20 bases: 69399016(96.418%)
Q30 bases: 67273806(93.4654%)

Filtering result:
reads passed filter: 800000
reads failed due to low quality: 0
reads failed due to too many N: 0
reads failed due to too short: 0
reads with adapter trimmed: 798732
bases trimmed due to adapters: 28022758

Duplication rate (may be overestimated since this is SE data): 0.33228%

- How many reads contained adapters in both datasets?
Jessica contains 705878 reads with adapters 
Fletcher contains 798732 reads with adapters

- What is the mean length of the reads before and after trimming?
We can use awk command for it. 
For Jessica after the filtering: 
```zcat jessica_trimmed_filtered.fastq.gz | awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}'```
59.1209 (before it is 87.8)
For Fletcher after the filtering: 
```zcat fletcher_trimmed_filtered.fastq.gz | awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}'```
89.9716 (before it is 125)


- Perform bwa alignment using aln and samse. For each sample, sort the sam file, save it as a bam, and index it. Remember that the reference is NOT the same as the previous exercise.
Align first
```bwa aln paeruginosa.fasta.gz jessica_trimmed_filtered.fastq.gz > jessica_trimmed_filtered_SA.sai```
```bwa aln paeruginosa.fasta.gz fletcher_trimmed_filtered.fastq.gz > fletcher_trimmed_filtered_SA.sai```
Output the sam
```bwa samse paeruginosa.fasta.gz jessica_trimmed_filtered_SA.sai jessica_trimmed_filtered.fastq.gz > output_jessica_aln.sam```
```bwa samse paeruginosa.fasta.gz fletcher_trimmed_filtered_SA.sai fletcher_trimmed_filtered.fastq.gz > output_fletcher_aln.sam```
Sort the sam and output as bam
```samtools sort output_jessica_aln.sam -o output_sorted_jessica_aln.bam```
```samtools sort output_fletcher_aln.sam -o output_sorted_fletcher_aln.bam```
index the bam 
```samtools index output_sorted_jessica_aln.bam```
```samtools index output_sorted_fletcher_aln.bam```


- For each sample, create new bam files with only the aligned reads. What proportion of reads remain?
```samtools view -F 4 output_sorted_jessica_aln.bam -o output_sorted_filtered_jessica_aln.bam```
Before filtered: 477112 reads 
After filtered: 175601 reads

```samtools view -F 4 output_sorted_fletcher_aln.bam -o output_sorted_filtered_fletcher_aln.bam```
Before filtered: 800000 reads
After filtered: 717393 reads

- What is the mean depth of coverage for each sample? (use samtools depth test.bam|datamash mean 3)
```samtools depth output_sorted_filtered_jessica_aln.bam|datamash mean 3```
jessica: 2.388 

```samtools depth output_sorted_filtered_fletcher_aln.bam|datamash mean 3```
fletcher: 10.775

- Use mapDamage to identify the nucleotide misincorporation and fragmentation patterns for each sample
```mapDamage -i ../../output_sorted_filtered_jessica_aln.bam -r ../../paeruginosa.fasta.gz```
```mapDamage -i ../../output_sorted_filtered_fletcher_aln.bam -r ../../paeruginosa.fasta.gz --no-stats```

Jessica is the ancient