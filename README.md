# MixSIH

MixSIH: a mixture model for single individual haplotyping.

## Features
MixSIH solves the haplotype assembly problem with mixture model. 

## Reference

Matsumoto, H., & Kiryu, H. (2013). MixSIH: a mixture model for single individual haplotyping. BMC Genomics, 14(Suppl 2), S5. [<a href="http://www.ncbi.nlm.nih.gov/pubmed/23445519">Pubmed</a>].

## How to build

```
git clone https://github.com/hmatsu1226/MixSIH
cd MixSIH
make
```

Or download from "Download ZIP" button and unzip it.

## Running MixSIH
##### Usage
```
./MixSIH <Option> <Input_file> <Output_file1> <Output_file2>
```

##### Options

* -a DOUBLE : Error rate. By default, it is set to 0.1.

##### Example of running MixSIH
```
./MixSIH -a 0.05 frag_sample.txt profile.txt haplotype.txt
```

##### Format of Input_file

* First line : Describes the number of lines of Input_file (which corresponds to the number of SNP fragments - 1).
* After first line : Each line describes a SNP fragment as below.
	- \segment_num \fragment_name \start_site1 \sequence1 \start_site2 \sequence2.....
	- \segment_num : The number of the segments which don't have gaps in the segments.
	- \fragment_name : The name of the fragment.
	- \start_allele'i' : The first site's position of i-th segment (1-origin).
	- \sequence'i' : The sequence of i-th segment.

The fragments must be sorted by the value of the third column and these can be sorted as follows:
```
sort -n -k 3 frag.txt > frag_sorted.txt
```

##### Example of Input_file
```
16
1 frag1 1 000
1 frag2 1 0001
2 frag3 1 11 4 111
2 frag4 2 11 5 11
1 frag5 3 000
1 frag6 4 000000
1 frag7 4 1111
3 frag8 5 11 8 1 10 1
1 frag9 6 000
1 frag10 7 0000
1 frag11 7 1011
1 frag12 7 111
1 frag13 8 010
1 frag14 9 11
1 frag15 9 00
```

##### Format of Output_file
##### Format of Output_file1

The Output_file1 is composed of the list of the blocks.
Each block contains a header which is composed of the relative position of the first site,
the number of the　sites in the block and the number of the sites which can be phased.
Each block consists of columns as follows:

* Col1: relative position of the site
* Col2: probability that the phase is (0,1)
* Col3: probability that the phase is (1,0)
* Col4: exp(connectivity) of the site

##### Example of Output_file1

```
BLOCK: offset: 1 len: 10 phased: 10
1    0.127    0.873    0.000
2    0.102    0.898    4.424
3    0.102    0.898    7.351
4    0.224    0.776    6.891
5    0.072    0.928    10.566
6    0.072    0.928    11.954
7    0.073    0.927    11.126
8    0.166    0.834    12.211
9    0.149    0.851    8.469
10    0.078    0.922    7.296
```

##### Format of Output_file2

The Output_file2 is composed of the list of the blocks.
Each block contains a header which is composed of the relative position of the first site,
the number of the sites in the block, the number of the sites which can be phased.
Each block consists of columns as follows:

* Col1: relative position of the site
* Col2: allele in the first haplotype
* Col3: allele in the second haplotype

##### Example of Output_file2
```
BLOCK: offset: 1 len: 10 phased: 10
1    0    1
2    0    1
3    0    1
4    0    1
5    0    1
6    0    1
7    0    1
8    0    1
9    0    1
10    0    1
```


## Extract reliable regions
extract_reliable_region.rb divides the haplotypes so that MC of the divided regions are higher than threshold.

##### Usage
```
ruby extract_reliable_region.rb <Input_file1> <Input_file2> <Output_file1> <Output_file2> <threshold>
```

##### Example
```
ruby extract_reliable_region.rb profile.txt haplotype.txt profile_6.txt haplotype_6.txt 6.0
```
Format of Input_file1 and Output_file1:
It is the same format of Output_file1 of MixSIH.

Format of Input_file2 and Output_file2:
It is the same format of Output_file2 of MixSIH.


## Dataset
For real data, we used the SNP fragments of Duitama's group [1] and it was downloaded from <a href="http://owww.molgen.mpg.de/~genetic-variation/SIH/data/">http://owww.molgen.mpg.de/~genetic-variation/SIH/data/</a>.

As the correct haplotypes, we used the haplotypes which are determined by pedigree genotypes and we downloaded these from 1000 Genomes Project <a href="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2010_07/trio/snps/CEU.trio.2010_03.genotypes.vcf.gz">(ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2010_07/trio/snps/CEU.trio.2010_03.genotypes.vcf.gz)</a>.

1. Duitama, Jorge, et al. "Fosmid-based whole genome haplotyping of a HapMap trio child: evaluation of Single Individual Haplotyping techniques." Nucleic acids research 40.5 (2012): 2041-2053.

## Slides
<iframe src="//www.slideshare.net/slideshow/embed_code/key/rOShgwymSEEHxJ" width="595" height="485" frameborder="0" marginwidth="0" marginheight="0" scrolling="no" style="border:1px solid #CCC; border-width:1px; margin-bottom:5px; max-width: 100%;" allowfullscreen> </iframe> <div style="margin-bottom:5px"> <strong> <a href="//www.slideshare.net/HirotakaMatsumoto/iscbasiasccg2012" title="MixSIH: a mixture model for single individual haplotyping" target="_blank">MixSIH: a mixture model for single individual haplotyping</a> </strong> from <strong><a href="//www.slideshare.net/HirotakaMatsumoto" target="_blank">Hirotaka Matsumoto</a></strong> </div>
