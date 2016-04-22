# RG (ReliableGenome)

ReliableGenome (RG) is a method for partitioning genomes into high and low concordance regions with respect to a set of surveyed VC pipelines. RG integrates variant call sets created by multiple pipelines from arbitrary numbers of input datasets and interpolates expected concordance for genomic regions without data resulting in a genome-wide concordance score.
Ultimately, genomic regions of high/low concordance are calculated from this genome-wide signal.	

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

Read more about RG in [this poster](docs/2016_RG_poster.pdf) that was presented at the [NGS16](https://www.iscb.org/ngs2016) conference.

## RG genomic partitions

* A genomic partition calculated from 219 deep WGS alignments can be [downloaded here](data/public/20160422_RG-win1000-score1_-3-RELIABLE-above0.5.bed.gz)([.tbi](data/public/20160422_RG-win1000-score1_-3-RELIABLE-above0.5.bed.gz.tbi)). The partition contains 2,209,778 concordance intervals located on chromosomes 1,..,22,X,Y.
  * Reference sequence: hs37d5 (GRCh37 + decoy)
  * Read mapper: bwa + stampy
  * Variant callers: GATK HaplotypeCaller, platypus, samtools
  * Params: w<sup>c</sup> = 1; w<sup>d</sup> = -3; t<sub>c</sub>=t<sub>d</sub>=0.5; window size x=1000
* The same partition with LCR and HD regions removed ([RG-LCR-HD](data/public/20160422_RG-win1000-score1_-3-RELIABLE-above0.5.bed.gz-min-LCR-min-HD.bed.gz)([.tbi](data/public/20160422_RG-win1000-score1_-3-RELIABLE-above0.5.bed.gz-min-LCR-min-HD.bed.gz.tbi))), see paper.
* The same partition with LCR100 and HD regions removed ([RG-LCR100-HD](data/public/20160422_RG-win1000-score1_-3-RELIABLE-above0.5.bed.gz-min-LCR100-min-HD.bed.gz)([.tbi](data/public/20160422_RG-win1000-score1_-3-RELIABLE-above0.5.bed.gz-min-LCR100-min-HD.bed.gz.tbi))), see paper.

---

## Building and Running RG 

### Prerequisites
* maven3 
* jdk 1.7+

### Building 
* Clone source code from github
* Run the bash script [build.sh](build.sh) 
* If everything worked you will find a standalone executable JAR in the bin directory and a library jar (containing only the classes from the RG source tree) in the target directory.

### Running RG

You can run RG via `java -jar bin/wtchg-rg-1.0.jar` which will print basic usage information. 
Use `java -Xmx12g -jar ...` to run RG with more dedicated heap space (recommended).

#### RG JOIN

To join VCF files from different variant callers, run:
`java -jar bin/wtchg-rg-1.0.jar CalcreliabilitySignals join`
to get detailed usage information.

Usage example:
`java -Xmx12g -jar bin/wtchg-rg-1.0.jar join 
    -d <GATK.VCF> -dl GATK 
    -d <PLAT.VCF> -dl PLAT 
    -d <SAMT.VCF> -dl SAMT 
    -o <SNV.out.vcf> 
    -oi <INDEL.out.vcf> 
    -dontCheckSort -dropAllFiltered -indelMergeWin 5`

#### RG CALC

To calculate the genome-wide concordance score signal, run:
`java -jar bin/wtchg-rg-1.0.jar CalcreliabilitySignals calc`
to get detailed usage information.

Usage example:
`java -Xmx12g -jar bin/wtchg-rg-1.0.jar calc 
    -o <OUTDIR> 
    -w 1000 
    -scoringSchema 1,-3 
    -thresholds 0.5,0.5 
    -dontCheckSort -v`

---

## Test data

Find some test VCF files that are ready to JOIN in [data/public/vcf/](data/public/vcf/).
Usage example (please modify paths to vcf/jar files as required): 
`java -Xmx12g -jar wtchg-rg-1.0.jar calc -o results -w 1000 -scoringSchema 1,-3 -thresholds 0.5,0.5 -createWigs -dontCheckSort -v -d vcf/AW_CRS_1631.DP+MDI.vcf.gz -d vcf/AW_CRS_1632.DP+MDI.vcf.gz -d vcf/AW_CRS_1806.DP+MDI.vcf.gz -d vcf/AW_CRS_1807.DP+MDI.vcf.gz -d vcf/AW_CRS_4103.DP+MDI.vcf.gz -d vcf/AW_CRS_4917.DP+MDI.vcf.gz -d vcf/AW_SC_4654.DP+MDI.vcf.gz -d vcf/AW_SC_4655.DP+MDI.vcf.gz -d vcf/AW_SC_4659.DP+MDI.vcf.gz
`
Please note that the "-createWigs" switch results in the creation of [WIG files](https://genome.ucsc.edu/goldenpath/help/wiggle.html) containing the genome-wide (interpolated) score signal and a signal showing the number of contributing datasets per position ("power signal"). The produced WIG files are too large to load them into a genome browser directly and should be converted, e.g., to the [BigWig format](https://genome.ucsc.edu/goldenpath/help/bigWig.html) using the following commandline
`wigToBigWig <WIG> <CHRSIZES> <BIGWIG>`.
(a chromosome sizes file [is provided here](data/public/hs37d5.fa) for convenience).
