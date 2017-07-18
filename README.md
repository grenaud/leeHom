# leeHom: Bayesian reconstruction of ancient DNA fragments

_QUESTIONS:_ gabriel [dot] reno [ at sign ] gmail.com

## About
leeHom is a program for the Bayesian reconstruction of ancient DNA (the binary is mergeTrimReadsBAM).


## Downloading:
Go to https://github.com/grenaud/leeHom and do one of the following:

1. Clone the repository: `git clone --recursive https://github.com/grenaud/leeHom.git`
2. [Download the ZIP archive](https://github.com/grenaud/leeHom/archive/master.zip) and then manually download the submodules

## Installation:
Make sure you have `cmake` installed. If required, install it by
typing `apt-get install cmake` for Ubuntu for example.

1. Build Bamtools first:
    ```bash
    cd bamtools/
    mkdir build/
    cd build/
    cmake ..
    make
    cd ../..
    ```

2. Build the submodules and main code by typing :
    ```bash
    make
    ```

## Running the program:
To launch the program simply type:
```
src/leeHom
```

or for if you have a multi-core machine, you can use them using the `-t` option in the following multi-threaded program:
```
src/leeHomMulti
```


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
BAM files can also be used to store unmapped reads. Ideally, leehom should be used prior to mapping since reads will map to a more accurate location if the adapters have been removed and overlapping stretches have been merged.


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
