# ibd-cluster

The **ibd-cluster** program infers multi-individual identity by descent (IBD)
from phased sequence data. The **ibd-cluster** program is able to analyze
biobank data with hundreds of thousands of individuals.

Version 0.2 of **ibd-cluster** incorporates a probabilistic model and
is more accurate than version 0.1

Last updated: May 9, 2025
Version: 0.2.0

## Contents

* [Installation](#installation)
* [Running ibd-cluster](#running-ibd-cluster)
  * [Required Parameters](#required-parameters)
  * [Optional Parameters](#optional-parameters)
* [Output files](#output-files)
* [License](#license)
* [Citation](#citation)
* [References](#references)

## Installation

You can download the most recent version of the **ibd-cluster** program
with the command:

    wget https://faculty.washington.edu/browning/ibd-cluster.jar

or you can download and compile the source files and create the
ibd-cluster.jar file with the commands:

    git clone https://github.com/browning-lab/ibd-cluster.git
    javac -cp ibd-cluster/src/ ibd-cluster/src/ibdcluster/ClustMain.java
    jar cfe ibd-cluster.jar ibdcluster/ClustMain -C ibd-cluster/src/ ./

[Back to Contents](#contents)

## Running ibd-cluster

The **ibd-cluster** program requires Java version 1.11 (or a later version).
Use of an earlier Java version may produce an "Unsupported Class Version" error.

The command:

    java -jar ibd-cluster.jar

prints a summary of the command line arguments.

To run **ibd-cluster**, enter the following command:

    java -Xmx[GB]g -jar ibd-cluster.jar [arguments]

where [GB] is the number of gigabytes of memory available for the analysis, and
[arguments] is a space-separated list of parameters. Each parameter
has the form **name=value**.

The shell script
[run.ibd-cluster.test](https://raw.githubusercontent.com/browning-lab/ibd-cluster/master/test/run.ibd-cluster.test)
will run an **ibd-cluster** test analysis.

[Back to Contents](#contents)

### Required Parameters

The **ibd-cluster** program has three required parameters:

* **gt=[file]** where **[file]** is a
[Variant Call Format](https://faculty.washington.edu/browning/intro-to-vcf.html)
(VCF) file containing a GT FORMAT subfield.  All genotypes must be phased, have
the phased allele separator ('|'), and have no missing alleles. If your data
is unphased, you can phase your data with the
[Beagle](https://faculty.washington.edu/browning/beagle/beagle.html) program.
A VCF record may have multiple ALT alleles. Any input VCF records that are
outside the genetic map are excluded from the analysis
(see the **map** parameter). If the VCF filename ends in ".gz" or ".bgz", the
VCF file is assumed to be compressed with gzip. If the VCF filename ends
in ".bref3", the VCF file is assumed to be compressed with
[bref3](https://faculty.washington.edu/browning/beagle/beagle.html).

* **map=[file]** where **[file]** is a
[PLINK format genetic map](http://zzz.bwh.harvard.edu/plink/data.shtml#map)
with cM units for the chromosome or chromosomes in the input VCF file.
We recommend the deCODE genetic map
[[1]](#references) for **ibd-cluster** analyses of human data.
The **ibd-cluster** program uses linear interpolation to estimate the genetic
position of loci between map positions.
Any input VCF records that are outside the genetic map are excluded from the
analysis. Each chromosome identifier in the input VCF file
must match a chromosome identifier in the map file. If the map filename ends
in ".gz", the map file is assumed to be compressed with gzip.

* **out=[string]** where **[string]** is the output filename prefix.

[Back to Contents](#contents)

### Optional Parameters

#### Data Parameters

The following three optional parameters restrict the data that will be
analyzed:

* **chrom=[chrom]:[start]‑[end]** specifies the chromosome or chromosome interval
to be analyzed: **[chrom]** is the CHROM identifier in the
input VCF file, **[start]** is the first base pair position in the interval, and
**[end]** is the last base pair position in the interval.
An entire chromosome, the beginning of a chromosome, or the end of a
chromosome may be specified with "**chrom=[chrom]**", "**chrom=[chrom]:‑[end]**",
and "**chrom=[chrom]:[start]‑**" respectively.  If no **chrom** parameter
is specified, all records in the input VCF file will be analyzed.

* **excludesamples=[file]** where **[file]** is a text file containing samples
(one sample per line) to be excluded from the analysis.

* **excludemarkers=[file]** where **[file]** is a text file containing markers
(one marker identifier per line) to be excluded from the analysis.
A marker identifier can be an identifier from the VCF record ID field, or it
can be a VCF record's CHROM and POS fields separated by a colon
(i.e. "CHROM:POS").

[Back to Contents](#contents)

#### Algorithm Parameters

All algorithm parameters are optional and have sensible default values.

* **min-maf=[0.0 < number ≤ 0.5]** specifies the minimum minor allele
frequency. Markers with minor allele frequency less than the minimum
frequency will be excluded from the analysis. For multi-allelic markers,
the minor allele frequency is defined to be the second-largest allele
frequency (**default: min-maf=0.1**).

* **min-ibs-cm=[number > 0.0]** specifies the minimum centimorgan length
of an identity-by-state (IBS) segment from which the endpoints of an
identity-by-descent (IBD) segment will be estimated (**default: min-ibs-cm=1.0**).

* **min-ibd-cm=[number > 0.0]** specifies the minimum centimorgan length
of an identity-by-descent (IBD) segment (**default: min-ibd-cm=1.0**).

* **pbwt=[integer > 0]** specifies the number of interleaved
Positional Burrows-Wheeler Transform analyses that will be performed
(**default: pbwt=4**).

* **trim=[number ≥ 0.0]** specifies the centimorgan length that will be
trimmed from each end of an IBD segment (**default: trim=0.5**).

* **discord=[0.0 < number < 1.0]** specifies the probability that two
alleles in a pairwise IBD segment are discordant (**default: discord=0.0005**).

* **out-cm=[number > 0.0]** specifies the centimorgan distance between
consecutive output  positions (**default: out-cm=0.02**).

* **nthreads=[integer ≥ 1]** specifies the number of computational threads.
The default **nthreads** parameter is the number of CPU cores.

[Back to Contents](#contents)

## Output files
The **ibd-cluster** program produces two output files: a **log** file, and
an **ibdclust** file.

The **log** file (.log) contains a summary of the analysis. The summary
includes:
* the analysis parameters
* the number of haplotypes after sample filtering
* the number of VCF records after marker filtering
* the number of positions in the **ibdclust** output file
* the mean number of IBD clusters per output position
* the allele discordance rate within trimmed IBD segments

The **ibdclust** file (.ibdclust.gz) reports the IBD cluster indices for
each haplotype at a sequence of equally-spaced genomic positions.
Each line of the **ibdclust** file is tab-delimited. The first line is a header
line which describes the data in each column. The fields of the header line
are "CHROM", "POS", "CM" and the sample identifiers from the input VCF file.
Each subsequent line contains the chromosome identifier, base coordinate, and
centimorgan coordinate of a genomic locus and two IBD cluster indices,
separated by '|', for each sample (one IBD cluster index for each haplotype).
The order of a sample's two IBD cluster indices at a locus matches the order
of the sample's two haplotypes in the input VCF file.


[Back to Contents](#contents)

## License
The **ibd-cluster** program is licensed under the Apache License, Version 2.0
(the License). You may obtain a copy of the License from
[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0).

[Back to Contents](#contents)

## Citation

If you use **ibd-cluster** in a published analysis, please report the program
version printed in the **log** file and cite the article describing
the **ibd-cluster** method:

> S R Browning, B L Browning (2025). Estimating gene conversion rates from
population data using multi-individual identity by descent.
bioRxiv 2025.02.22.639693;
doi: [https://doi.org/10.1101/2025.02.22.639693](https://doi.org/10.1101/2025.02.22.639693)

This study used the python script
[find.gcs.py](https://raw.githubusercontent.com/browning-lab/ibd-cluster/master/src/find.gc.py)
to detect alleles changed by gene conversion.

## References

**[1]** Halldorsson BV, Palsson G, Stefansson OA, Jonsson H, Hardarson MT,
Eggertsson HP, Gunnarsson B, Oddsson A, Halldorsson GH, Zink F, Gudjonsson SA,
Frigge ML, Thorleifsson G, Sigurdsson A, Stacey SN, Sulem P, Masson G,
Helgason A, Gudbjartsson DF, Thorsteinsdottir U, Stefansson K.
Characterizing mutagenic effects of recombination through a sequence-level
genetic map. Science. 2019 Jan 25;363(6425):eaau1043.
doi: [10.1126/science.aau1043](https://doi.org/10.1126/science.aau1043).
PMID: [30679340](https://pubmed.ncbi.nlm.nih.gov/30679340/).

[Contents](#contents)



