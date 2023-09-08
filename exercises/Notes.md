# Advance Bioinformatics Notes

## Sequence 
You might see capital and lowercase in reference genome (.fasta file): 
Capital Letters
Functional Regions: Sometimes, uppercase letters are used to indicate regions of known function, such as exons in a gene sequence.

Lowercase Letters
Non-functional Regions: In some annotated sequences, lowercase letters might be used to indicate non-coding or non-functional regions, like introns or intergenic regions.

## Important General Commands: 
1. tr -> tr command is a Linux command-line utility that translates or deletes characters from standard input (stdin). often used with pipes (|) and redirects (>>)
2. cut -> select particular fields (cut -f 1,2 will select column 1 and 2 only) 
3. cat -> Print all of the text from a text file (not recommended)
4. gunzip -> unzipping compressed (zipped) files: 'gunzip -k zip_file.zip' [Dont forget to '-k' for keeping the original files, otherwise it will be overwritten]
5. bgzip -> bcftools requires that files are compressed with bgzip. So you have to: gunzip -k chr21.fa.gz & bgzip chr21.fa
6. zcat -> Print all of the text from a zipped file (not recommended)
7. less -> Looking at text file (press "q" to close)
8. 'zless' or 'zmore' -> Opening zipped file like using "less" (press "q" to close)
9. wc -> word counting. To count for the lines: 'wc -l' characters: 'wc -m'
10. grep -> finding pattern from data (return the line(s) wich contain the patern. We can also count the number of lines it found). Examples:  
   finding how mane sequence from a fasta file: grep -c '^>' fasta_file.fa
11. sed -> Stream editor. It's like grep but it has ability to edit the file. looping, and so on. 
   It can execute several command with: 
   1. sed -e 'command1' -e 'command2' target_files.txt 
   2. sed -f file_contains_commands target_files.txt
   sed can be used to see e.g. the zipped file of FASTQ file. : zcat data.fastq.gz | sed -n '2~4p'  | head -5
   It means that sed (with slient -n) print start from the second line and reprint every 4th line: 2, 6, 10 ..... and print only five seqeunces.  
12. awk -> Pattern finding and processing. More complex than "sed".
   command: "awk options 'selection_criteria {action }' input-file > output-file"
   examples: 
   1. printing second and thrid column of file awk: '{print $2 "\t" $3}' file.txt
   2. printing all columns and linse: awk ' {print $0}' file.txt
   3. printing the 3rd and 4th columns where the letter ‘a’ appears in either: awk '/a/ {print $3 "\t" $4}' file.txt
   4. Printing the sequences from a fasta file :  zcat DATA_L001_R1.fastq.gz | awk 'NR%4==2{print $0}'|head -5    
   (NR%4==2 is a condition that checks if the line number modulo 4 is equal to 2. It is a selection_criteria before doing action which is print.)
   5. Printing the average length of the reads: zcat DATA_L001_R1.fastq.gz | awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}'
   more resources: https://www.grymoire.com/Unix/Awk.html
13. datamash -> It is a useful command for doing statistical thing i.e. max, min, avg on a file just like tidyverse in R.
   The file scores.txt contains tests scores of college students of different majors (Arts, Business, Health and Medicine, Life Sciences, Engineering, Social Sciences).
   The files has three columns: Name, Major, Score:
   Shawn     Arts  65
   Marques   Arts  58
   Fernando  Arts  78
   Paul      Arts  63
   Walter    Arts  75
   ...

   Using datamash, find the lowest (min) and highest (max) score for each College Major (Major is in column 2, the score values are in column 3):
   $ datamash -g 2 min 3 max 3 < scores.txt
   Arts            46  88
   Business        79  94
   Health-Medicine 72  100
   Social-Sciences 27  90
   Life-Sciences   14  91
   Engineering     39  99 

   It can also do the average length of the read as we did on awk: zcat DATA_L001_R1.fastq.gz | awk 'NR%4==2{print length($0)}' | datamash mean 1
   more resources: https://www.gnu.org/software/datamash/examples/#example_files

## Quality Control and Preprocessing 
### FASTQC
provide a comprehensive quality assessment of raw sequencing data. It's a diagnostic tool to understand the quality of your data before and after preprocessing.
- User Interface: Provides a graphical user interface (GUI) as well as a command-line interface.
- Output: Generates a detailed report that includes multiple visualizations, such as base quality distribution, GC content, sequence length distribution, sequence duplication levels, and more.

### FASTP 
- Multifunctionality: fastp is an all-in-one tool designed for both quality control and preprocessing. It can perform tasks like adapter trimming, quality filtering, error correction, and more, all within a single command.
- Auto-Detection: fastp can automatically detect and remove adapter sequences without requiring the user to specify them, although you can also manually specify if needed.
- Performance: It is optimized for high performance and can process data very quickly.
- Quality Reports: One of fastp's standout features is its capability to generate HTML and JSON reports summarizing data quality before and after processing.
- Support for Paired-End Data: Designed to handle both single-end and paired-end reads efficiently.
- Ease of Use: Due to its auto-detection capabilities and the comprehensive reports it generates, fastp is relatively easy to use, especially for those new to sequencing data preprocessing.
- Quality Control: It provides both read-level and base-level quality control. You can filter or trim reads based on their overall quality scores, specific base positions, or other criteria.
more resources on it: https://github.com/OpenGene/fastp
example: 
fastp --in1 ../../data/fastq/DATA_L001_R1.fastq.gz \ 
--in2 ../../data/fastq/DATA_L001_R2.fastq.gz \ 
--out1 DATA_L001_R1_trimmed.fastq.gz \ 
--out2 DATA_L001_R2_trimmed.fastq.gz \ 
--merge \ 
--merged_out DATA_L001_merged_trimmed.fastq.gz \ 
--unpaired1 DATA_L001_unpaired_R1.fastq.gz \
--unpaired2 DATA_L001_unpaired_R2.fastq.gz \ 
--length_required 30 \ 
--detect_adapter_for_pe


## BWA
### Index 
he bwa index command is used in bioinformatics to create an index for a reference genome. BWA (Burrows-Wheeler Aligner) is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. Creating an index of the reference genome is a prerequisite for performing efficient mapping or alignment of sequencing reads against that genome.
After running bwa index, you'll get several files that collectively make up the index for the reference genome. The files typically have the following extensions:

.amb: Contains information regarding the length of each sequence, as well as the number of sequences.
.ann: Contains information about sequence names and lengths.
.bwt: The Burrows-Wheeler Transform of the reference sequences, a crucial data structure for efficient alignment.
.pac: The reference sequences in a compressed, binary format.
.sa: The suffix array corresponding to the Burrows-Wheeler Transform, another crucial data structure for alignment.
These files are generated in the same directory where the original reference FASTA file resides, and they usually have the same prefix as the original FASTA file (i.e., if the FASTA file is reference.fasta, the index files will start with reference).

### aln, samse, and sampe
The bwa aln and bwa samse (for single-end reads) or bwa sampe (for paired-end reads) commands work in tandem to align sequencing reads to a reference genome and produce an aligned SAM (Sequence Alignment/Map) file.
```bwa aln -t 2 reference_file readfile_1 > readfile1_SA.sai```
Then we will generate alignments in the SAM format given the single-end read by searching in the SA intervals
```bwa samse reference_file readfile1_SA.sai readfile_1 > output_aln.sam```


### MEM 
The bwa mem command is another alignment algorithm offered by BWA, and it is often recommended for high-throughput sequencing data, especially for longer reads (usually 100bp or longer). The mem in bwa mem stands for Maximal Exact Matches, which is the algorithm's foundational concept. Unlike bwa aln followed by bwa samse or bwa sampe, bwa mem performs the alignment and generates the SAM file in a single step, making it more convenient and usually faster for most modern sequencing data types.
``` bwa mem -t 2 reference_file readfile_1 > output_mem.sam ```


## SAMTOOLS
### quickcheck
quick check whether it is corrupted before doing bam analysis
```samtools quickcheck DATA.bam```

### faidx 
Slicing faster into FASTA files. It creates file.fa.fai for its index. Afterward, we can slice e.g. with 
``` samtools faidx file.fa seq_name:start-stop```

### view 
The samtools view command is mainly used to manipulate and filter SAM and BAM files, which are formats for storing sequence alignment data.

- Input/Output: Can convert between SAM,CRAM, and BAM formats.
- Filtering: Filters reads based on group, target, tag, mapping quality, flags, or other attributes.
Use Case for Filtering Unaligned Reads ("view -F 4"):
- Reducing Noise: Unaligned reads don't map to your reference genome, which means they might just add noise to downstream analyses like variant calling.
- Specialized Analyses: Sometimes you may be interested in the reads that do not align to the reference genome, perhaps because you're looking for novel sequences or contaminants. In such cases, you'd filter to keep the unaligned reads for separate analysis.
- Reducing File Size: Unaligned reads take up space. 
- Alignment Quality: In certain contexts, a focus on aligned reads can improve the robustness and reliability of downstream analyses such as differential expression in RNA-seq or variant identification in DNA-seq.

Count the number of unmapped reads

```samtools view -f 4 -c DATA.bam```

### sort
Use Case for Sorting BAM Files:
- Variant Calling: Tools like GATK or FreeBayes usually require sorted BAM files to efficiently call variants.
- Genome Browsers: Programs like IGV work more efficiently with sorted and indexed BAM files, as this allows for quicker data retrieval.
- Merging BAM Files: If you're combining multiple BAM files into one, having each sorted in the same way is crucial for the final merged file to be valid and useful.
- Coverage Calculations: Tools that calculate read depth or coverage statistics usually require sorted BAM files to operate efficiently.
- Efficiency: Sorted files are generally more efficient for data storage and retrieval, which can be crucial when dealing with large genomic datasets.
Finally, often you can also have your aligner write directly to samtools sort:
bwa mem genome.fa reads.fastq | samtools sort -o myfile_sorted.bam

### flagstat
used to knowo the statistics of each flag status. good to know i.e. for percentage QC pass or mapped reads

### stat
contains comprehensive statistics of sam files

### tview 
text alignment viewer (based on the ncurses library). In the viewer, press `?' for help and press `g' to check the alignment start from a region in the format like `chr10:10,000,000' or `=10,000,000' when viewing the same reference sequence.
The top line shows the reference sequence, or 'N's if unknown. Underneath this is the consensus, derived from the sequence alignments. Below the consensus the sequence alignment records are shown. Uppercase and lowercase is used to distinguish the sequence strand, with uppercase being the top/forward strand.
When the reference is known, both consensus and alignment record sequences are displayed in a dot-notation where a matching character is shown as '.' (forward strand) or ',' (reverse strand) and only mismatching bases and missing bases are shown. This mode can be toggled with the "." command.
```samtools tview DATA.bam -p chr21:14338386-14338400```

### mpileup
Use for summarizing at a base level with the following columns: 
- Chromosome name.
- 1-based position on the chromosome.
- Reference base at this position (this will be “N” on all lines if -f/--fasta-ref has not been used).
- Number of reads covering this position.
Read bases. This encodes information on matches, mismatches, indels, strand, mapping quality, and starts and ends of reads.
- Base qualities
- Map qualities
command for obtaining the depth distribution:
```samtools mpileup DATA.bam | cut -f4 | sort -n | uniq -c  > mpileup_depth.txt```

## Summary on Mapping Tools
### BWA (Burrows-Wheeler Aligner)
- Purpose: DNA sequence alignment.
- Algorithm: Utilizes the Burrows-Wheeler Transform to create a compressed index of the reference genome.
- File Types: Works primarily with DNA sequences in FASTA and FASTQ formats.
- Features: Includes different algorithms like BWA-backtrack, BWA-SW, and BWA-MEM for various types of alignment tasks. BWA-MEM is generally recommended for high-quality reads over 100 bp.
- Use-Cases: Genomic resequencing, variant calling, structural variant discovery.

### Bowtie2
- Purpose: DNA sequence alignment.
- Algorithm: Also uses the Burrows-Wheeler Transform for indexing but offers different speed and sensitivity settings.
- File Types: Works with DNA sequences in FASTA and FASTQ formats.
- Features: Local and end-to-end alignment, multiple mapping strategies.
- Use-Cases: Often used for transcriptome alignments, genomic alignments, ChIP-Seq.

### HISAT2
Purpose: DNA and RNA sequence alignment.
Algorithm: Extension of Bowtie2 designed for efficient RNA-Seq alignment. Utilizes a hierarchical indexing scheme.
- File Types: DNA and RNA sequences in FASTA and FASTQ formats.
- Features: Spliced alignments, suitable for RNA-Seq data. Detects exon-intron boundaries.
- Use-Cases: Specifically optimized for RNA-Seq data, gene expression quantification.

### Salmon
- Purpose: Transcript quantification.
- Algorithm: Does not perform alignment in the traditional sense. Instead, it uses quasi-mapping or lightweight alignment to quantify abundances of transcripts.
- File Types: RNA sequences in FASTA and FASTQ formats. Can also take pre-aligned reads in SAM/BAM format.
- Features: Extremely fast and memory-efficient. Provides measures like Transcripts Per Million (TPM) and estimated counts.
- Use-Cases: RNA-Seq data analysis, differential gene expression.

## Variant Calling 
### bcftools 
```bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf```
The first mpileup part generates genotype likelihoods at each genomic position with coverage. The second call part makes the actual calls. The -m switch tells the program to use the default calling method, the -v option asks to output only variant sites, finally the -O option selects the output format. In this example we chosen binary compressed BCF, which is the optimal starting format for further processing, such as filtering.
