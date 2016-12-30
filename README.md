===============
The molecular landscape of cutaneous neurofibromas
===============
This project is focused on analyzing the molecular landscape of cutaneous NF and is
funded by the Children's Tumor Foundation.  All available data and key analysis
steps uploaded to http://www.synapse.org/cutaneousNF

To access the data you will need to register with Synapse and click through to
get approval to use the data for your own research. One you are authenticated,
there are numerous scripts provided here to facilitate analysis of the data.
These scripts access the data through the [Synapse R and Python
Clients](http://docs.synapse.org/articles/getting_started.html#installing-synapse-clients).


## Data available for download

This project contains diverse types of high-throughput data. Links to the
analysis of each type of data can be found here:
* Copy Number Alteration data from SNP OMNI Arrays
* Whole Genome Sequencing
* RNA-Sequencing
* iTRAQ Proteomics

Each dataset is described in more detail below.

### SNP Arrays
SNP Array data collection is described [on the
wiki](https://www.synapse.org/#!Synapse:syn4984604/wiki/400306) with files
available as output from GenomeStudio.  Files can be downloaded in their
original form using the `cnv_unprocessed_files` function in `bin/dermalNFData.R`
script available [here](./bin/dermalNFData.R).

### Whole Genome Sequencing (WGS)
The WGS data was processed in numerous ways to compare and contrast various
methods of calling variants and mutation calls. More details are available [on
the wiki](https://www.synapse.org/#!Synapse:syn4984604/wiki/400307).

TODO: add

### RNA-Seq
RNA-Seq collection is described [on the
wiki]() with BAM files and RNA quantities available for download.

Gene matrices can be downloaded via various functions in `bin/dermalNFData.R`
script available [here](./bin/dermalNFData.R).

## Structure of this repository
Binaries used to process this data can be found in [the binaries directory](./bin/) while specific
analytical tasks can be found in [the analysis directory](./analysis).

TODO: fill in with more details
