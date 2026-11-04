# leeHom: Bayesian reconstruction of ancient DNA fragments

QUESTIONS: gabriel [dot] reno [ at sign ] gmail.com

## About
leeHom is a program to trim adapters and merge overlapping sequencing using a Bayesian reconstruction. This is particularly useful for ancient DNA and or DNA shorter fragment sizes like cell-free DNA.


## Downloading:
Go to https://github.com/grenaud/leeHom and do one of the following:

1. Clone the repository: `git clone --recursive https://github.com/grenaud/leeHom.git`
2. [Download the ZIP archive](https://github.com/grenaud/leeHom/archive/master.zip) and then manually download the submodules

## Installation:
Make sure you have `autoconf` and zlib (`zlib1g-dev`) installed, e.g. `apt-get install autoconf zlib1g-dev` on Ubuntu. The BAM files are read and written with [htslib](https://github.com/samtools/htslib), which is built from the submodule.

1. Build the submodules and main code by typing :
   
    make
   
To use an htslib that is already installed, type `make HTSLIBINC=/path/to/include-parent HTSLIBLIB=/path/to/lib`. `make test` runs the tests.

## Running the program:
To launch the program simply type:
```
src/leeHom
```




## Running faster
- `-t [threads]` uses several cores. The output is the same, and in the same order as the input, whatever the number of threads (use `--seed` to also fix the random choice between equally good bases).
- `--quick` makes the most expensive step faster, at the price of a very small loss in sensitivity. leeHom compares every possible position of the adapter, and for each one it adds the log likelihood of every pair of overlapping bases. With `--quick`, once a position has compared its first `--quickbases` bases (default 12) it is only computed in full if the likelihood reached so far is clearly better than what random sequences would give. Otherwise the remaining terms are added as their expected value under random sequences (3/4 mismatches, 1/4 matches), so the number of terms of the likelihood is unchanged. `--quickfraction [0-1]` sets how far between the likelihood of random sequences (0) and of a perfect overlap (1) a position must be at that point (default 0.4).
  `--quick` only speeds up the likelihood computation, and the gain is modest: on simulated 75 bp reads (1M pairs, one thread, output not compressed) the defaults are 5 to 7% faster, `--quickbases 12 --quickfraction 0.5` about 10% faster (about 0.1% of the merged pairs are kept as pairs instead), and the most aggressive settings (about 25%) lose over 10% of the merges. A `--quickfraction` of 0, the literal "below what random sequences give" rule, is slower than not using `--quick`. Larger `--quickbases` or smaller `--quickfraction` are more cautious and give less speedup.
- `--gzlevel [0-9]` sets the compression level of the output (fastq or BAM). When most reads are written compressed, compression is often more than half of the CPU time: `--gzlevel 1` is much faster than the default (6 for fastq) for files about 20% bigger.

## Fastq input, BAM output
With fastq input, `-o [file]` (or `-u` for an uncompressed BAM on standard output, good for a pipe) writes an unaligned BAM instead of the fastq files, with the mates of a pair consecutive (flags 77 and 141) and single reads with flag 4. It cannot be combined with `-fqo`/`-fqt`. The reads are the same as in the fastq output and in the same order as the input. `--phred64` is supported. Example:
```
src/leeHom --ancientdna -fq1 reads_1.fq.gz -fq2 reads_2.fq.gz -t 4 -u | samtools sort -o sorted.bam -
```

## Interpreting the output using the fastq mode
If the input is fastq, the output will also be fastq. Normally you will get different files. Here is a table explaining their meaning/content. The **[PREFIX]** is the output prefix specified using -fqo.


### Single-end mode:

| File suffix                 | Meaning                                                                                                                                                                                                                |
|-----------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **[PREFIX].fq.gz**            | Sequences that were either trimmed or untrimmed confidently by leeHom                      |
| **[PREFIX].fail.fq.gz**       | Sequences where the most likely original sequence length was **not** several fold more likely than the second most likely sequence length. The most likely sequence length is are nonetheless found in the file.                      |


### Paired-end mode:

| File suffix                 | Meaning                                                                                                                                                                                                                |
|-----------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **[PREFIX].fq.gz**            | Sequences that were trimmed and merged confidently by leeHom                      |
| **[PREFIX]_r1.fq.gz**         | Forward reads that were neither trimmed nor merged confidently by leeHom                      |
| **[PREFIX]_r2.fq.gz**         | Reverse reads that were neither trimmed nor merged confidently by leeHom                      |
| **[PREFIX].fail.fq.gz**       | Sequences that were trimmed and merged by leeHom and where the most likely original sequence length was **not** several fold more likely than the second most likely sequence length                      |
| **[PREFIX]_r1.fail.fq.gz**    | Forward reads that were neither trimmed nor merged by leeHom and where the most likely original sequence length was **not** several fold more likely than the second most likely sequence length                      |
| **[PREFIX]_r2.fail.fq.gz**    | Reverse reads that were neither trimmed nor merged by leeHom and where the most likely original sequence length was **not** several fold more likely than the second most likely sequence length                      |





## Interpreting the log messages
By default, leehom prints a log messsage to stderr. leehom counts 2 paired-end reads or one single-end read as one cluster. This is how to interpret the log:

| Field                       | Meaning                                                                                                                                                                                                                |
|-----------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Total**                   | Total number of clusters analyzed by leehom                                                                                                                                                                            |
| **Merged (trimming)**       | Number of clusters where the adapter was found and, for paired-end reads, a likely overlap between the forward and reverse reads was found                                                                             |
| **Merged (overlap)**        | When using a prior or `--ancientdna` and paired-end data, the number of clusters where an likely overlap between the ends of paired-end reads was found and reads were merged as a single molecule                     |
| **Kept PE/SR**              | Cluster where no sufficient evidence for merging paired-end or trimming single-end reads was found                                                                                                                     |
| **Trimmed SR**              | When using single-end reads, the number of reads where a match to the adapter was found                                                                                                                                |
| **Adapter dimers/chimeras** | Number of cluster where a chimeric sequence (two adapters merged together without insert) was found and were flagged as QC failed                                                                                      |
| **Failed Key**              | If the user provides a key (sequences must start with this sequences is it was part of the adapter and not covered by the sequencing primer) the number of clusters where the key was not found will be reported here. |


## How / When to set a prior?
If you have previous data and have a representative sample of reasonable size for your library, using a prior will give you the best reconstruction possible. leehom models the distribution as a log-normal one. To infer the parametes of the log-normal distribution, you can use src/approxDist.R to get the parameters and use those as prior in leehom. To create this data from aligned modern DNA or (preferably aligned) pre-trimmed ancient, you can launch the following command on the BAM file:

```bash
samtools view aligned.bam | gawk 'and($2,5)==0 { print length($10) }; and($2,2)==2 && $9>0 { print $9 }' | gzip  > sizedata.dat.gz
```

(We recommend the use of this command with BWA especially as it uses the properly paired flag).

Given that the resulting file has one molecule length per line. You can use the following R script to get the maximum likelihood parameters for the log-normal distribution:

```
src/approxDist.R testData/sizedata.dat.gz
Loading required package: survival
Loading required package: splines
     meanlog         sdlog
  5.1444983185   0.2885764766
 (0.0009125589) (0.0006452766)
```

The first number (`5.1444983185`) represents the "Location" for log-normal distribution whereas the second number (`0.2885764766`) represents the "Scale".
Those parameters can now be used in leehom as `--loc 5.1444983185 --scale 0.2885764766`.

If you have no previous data, you can launch leehom using a uniform prior (default parameters). If you have modern DNA with relatively long insert sizes, you can use leehom as is. If you have ancient DNA and/or suspect that molecules are short, use the `--ancientdna` which will merge reads that share only a partial overlap if they have sufficient probabilistic support.


## What if the adapters are stored in fasta files?
You can use a file descriptor for each. For example if they are stored as FASTA:

```bash
src/leeHom -f `cat /mydata/forward.fa | tail -n+2 | tr -d "\n"` -s `cat /mydata/reverse.fa | tail -n+2 | tr -d "\n"`
```

## How do I get the sequence of the adapters?
Contact your sequencing center or core and ask them to send you the sequence of the sequencing adapters. They should be known in advance since, currently, leehom cannot infer the sequence but uses this information as a prior.

## Does my BAM file have to be in a specific order?
For single-end, no. For paired-end, mate have to be consecutive in the file (i.e. `First mate 1`, `second mate 1`, `First mate 2`, `second mate 2`, _etc.._ ). This can be done by either using raw data from the sequencer or sorting by name.

## Does my BAM file have to be unmapped?
BAM files can also be used to store unmapped reads. leeHom should be used prior to mapping since reads will map to a more accurate location if the adapters have been removed and overlapping stretches have been merged.


## What do the `FF` tags in the BAM file mean?
Here is a table with the meaning for each FF flag:

| FF            |                     Meaning                        |
|-------------- |----------------------------------------------------|
| 1             | single-end reads trimmed                           |
| 2             | paired-end merged overlapping portions             |
| 3             | paired-end trimmed and merged overlapping portions |


## Test data:

A small test set of 50,000 ancient DNA sequences from the Altai Neandertal can be found here: `testData/rawAncientDNA.bam`

To launch the reconstruction program for this test data set, use the following command:
```bash
src/leeHom -f AGATCGGAAGAGCACACGTCTGAACTCCAG -s GGAAGAGCGTCGTGTAGGGAAAGAGTGTAG --ancientdna -o testData/reconsAncientDNA.bam testData/rawAncientDNA.bam
```

The reconstructed ancient DNA fragments will be in `testData/reconsAncientDNA.bam`

For fastq data:
```bash
src/leeHom -f AGATCGGAAGAGCACACGTCTGAACTCCAG -s GGAAGAGCGTCGTGTAGGGAAAGAGTGTAG --ancientdna -fq1 testData/rawAncientDNA.f1.gz -fq2 testData/rawAncientDNA.f2.gz -fqo testData/outfq
```

The output files will have the `testData/outfq*` prefix.

Alternatively, you can specify the names of all 6 files:

```bash
src/leeHom -f AGATCGGAAGAGCACACGTCTGAACTCCAG -s GGAAGAGCGTCGTGTAGGGAAAGAGTGTAG --ancientdna -fq1 testData/rawAncientDNA.f1.gz -fq2 testData/rawAncientDNA.f2.gz -fqose testData/outfq.fq.gz -fqosef testData/outfq.fail.fq.gz -fqope1 testData/outfq_r1.fq.gz -fqope1f testData/outfq_r1.fail.fq.gz -fqope2 testData/outfq_r2.fq.gz -fqope2f testData/outfq_r2.fail.fq.gz
```


