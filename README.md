# ibd-cluster

The **ibd-cluster** program uses multi-individual identity by descent to
cluster haplotypes. Haplotypes in the same cluster at a genomic position
are identical by descent at that position.

The **ibd-cluster** program can analyze biobank-scale sequence data.

Last updated: November 6, 2023  
Current version: 0.1.0

## Contents

* [Installation](#installation)
* [Running ibd-cluster](#running-ibd-cluster)
  * [Required parameters](#required-parameters)
  * [Optional parameters](#optional-parameters)
* [Output files](#output-files)
* [License](#license)
* [Citation](#citation)

## Installation

You can download the latest executable file,
[ibd-cluster.jar](https://faculty.washington.edu/browning/ibd-cluster.jar),
with the command:

    wget https://faculty.washington.edu/browning/ibd-cluster.jar

or you can download the source files and create the executable file
with the commands:

    git clone https://github.com/browning-lab/ibd-cluster.git
    javac -cp ibd-cluster/src/ ibd-cluster/src/ibdcluster/ClustMain.java
    jar cfe ibd-cluster.jar ibdcluster/ClustMain -C ibd-cluster/src/ ./
    jar -i ibd-cluster.jar

[Contents](#contents)

## Running ibd-cluster

The **ibd-cluster** program requires Java 1.8 or a later version. Use of an
earlier Java version will produce an "Unsupported Class Version" error.

The command:

    java -jar ibd-cluster.jar

prints a summary of the command line arguments.

You can run an **ibd-cluster** analysis with the command:

    java -Xmx[GB]g -jar ibd-cluster.jar [arguments]

where **[GB]** is the maximum number of gigabytes of memory to use, and
**[arguments]** is a space-separated list of parameter values, each expressed as
**parameter=value**.  If a Java out of memory error occurs, you probably will
need to increase the Java -Xmx parameter.

The shell script
[run.ibd-cluster.test](https://raw.githubusercontent.com/browning-lab/ibd-cluster/master/test/run.ibd-cluster.test)
runs a test **ibd-cluster** analysis.

[Contents](#contents)

### Required parameters

The **ibd-cluster** program has three required parameters. All other parameters
are optional and have reasonable default values.  Input files having a name
ending in ".gz" are assumed to be gzip-compressed.

* **gt=[file]** where **[file]** is a
[Variant Call Format](https://faculty.washington.edu/browning/intro-to-vcf.html)
(VCF) that contains phased, nonmissing genotype data for each individual.
If the VCF file has unphased or sporadic missing genotypes, you can phase and
impute the genotypes using the
[Beagle](https://faculty.washington.edu/browning/beagle/beagle.html) program
before running **ibd-cluster**.  If the VCF filename ends in ".bref3", it is
assumed to be bref3-compressed.  Software for bref3 compression and
decompression can be downloaded from the
[Beagle web site](https://faculty.washington.edu/browning/beagle/beagle.html).

* **map=[file]** where **[file]** is a
[PLINK format genetic map](https://zzz.bwh.harvard.edu/plink/data.shtml#map)
with cM units. Positions of markers that are between genetic map positions are
estimated using linear interpolation. The chromosome identifiers
in the genetic map and the input VCF files must match. Analysis is restricted
to the region within the genetic map. HapMap genetic maps
in cM units are available for
[GRCh36](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/),
[GRCh37](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/), and
[GRCh38](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/).

* **out=[string]** where **[string]** is the output filename prefix.

[Contents](#contents)

### Optional parameters

* **chrom=[chrom]:[start]‑[end]** specifies the chromosome or chromosome interval
to be analyzed: **[chrom]** is the chromosome identifier in the
input VCF and map files, **[start]** is the first base pair coordinate, and
**[end]** is the last base pair coordinate.
An entire chromosome, the beginning of a chromosome, or the end of a
chromosome may be specified with "**chrom=[chrom]**", "**chrom=[chrom]:‑[end]**",
and "**chrom=[chrom]:[start]‑**" respectively. If a **chrom** parameter is not
specified, all chromosomes in the input VCF file will be analyzed.

* **length=[number > 0.0]** specifies the minimum cM length of an
identity-by-state (IBS) segment that will be used to cluster haplotypes
(**default: length=2.0**). The IBS segment will be trimmed prior to
clustering. The **length** parameter must be at least twice as large as the
**trim** parameter. See the **trim** parameter for more information.

* **trim=[number > 0.0]** specifies the cM length that will be trimmed
from each end of an identity-by-state (IBS) segment (**default: trim=1.0**).
The **length** parameter must be at least twice as large as the **trim**
parameter. See the **length** parameter for more information.

* **min-maf=[number < 0.5]** specifies the minimum minor allele frequency
(**default: min-maf=0.1**). Markers with minor allele frequency less than
**min-maf** will be ignored. For multi-allelic markers, the minor allele
frequency is the second-largest allele frequency.

* **aggregate=[number > 0.0]** specifies the maximum cM length of a
haplotype segment that will be treated as a single allele
(**default: aggregate=0.005**).

* **nthreads=[integer ≥ 1]** specifies the number of computational threads to
use for the analysis. The default **nthreads** parameter is the number of
CPU cores.

* **excludesamples=[file]** where **[file]** is a text file containing samples
(one sample per line) that are to be excluded from the analysis.

* **excludemarkers=[file]** where [file] is a text file containing markers
(one marker identifier per line) that are to be excluded from the analysis.
A marker identifier can be an identifier from the VCF record ID field, or it
can be a VCF record's CHROM and POS fields separated by a colon
(i.e. "CHROM:POS").

[Contents](#contents)

## Output files
The **ibd-cluster** program produces two output files: a **log** file and an
**IBD cluster** file.
* The **log** file (.log) contains a summary of the analysis.
* The **IBD cluster** file (.ibdclust.gz) is a tab-delimited file that
contains the haplotype clusters at each position.

The first line of the **IBD cluster** file is a header line that describes
the data in each column.  The first three fields of the header line are
"CHROM", "POS", "CM".  These three columns contain the chromosome (CHROM),
base position (POS), and cM position (CM) of the positions at which IBD
clustering is performed.  The remaining fields of the header line are the
sample identifiers. These columns contain two haplotype clusters for each
individual at each position.

An individual's two haplotype clusters are separated by a vertical bar ('|').
The first and second haplotype clusters respectively contain the individual's
first and second haplotypes in the input VCF file.
The haplotype clusters at each position are indexed by consecutive integers
starting with 0. Haplotypes with the same cluster index at a position are
identical by descent at that position. The cluster indices at a position
apply only to that position. The same cluster index may correspond to
different clusters at different positions.

[Contents](#contents)

## License
The **ibd-cluster** program is licensed under the Apache License,
Version 2.0 (the License). You may obtain a copy of the License from
[https://www.apache.org/licenses/LICENSE-2.0](https://www.apache.org/licenses/LICENSE-2.0)

[Contents](#contents)

## Citation

If you use **ibd-cluster** in a published analysis, please report the program
version printed in the **log** file and cite the article describing
the **ibd-cluster** method:

> S R Browning, B L Browning (2023). Biobank-scale inference of multi-individual
identity by descent and gene conversion.  bioRxiv 2023.11.03.565574;
doi: https://doi.org/10.1101/2023.11.03.565574

[Sharon Browning](https://sites.uw.edu/sguy/) developed the **ibd-cluster** clustering method.  
[Brian Browning](https://faculty.washington.edu/browning) developed the **ibd-cluster** algorithms and software.

[Contents](#contents)
