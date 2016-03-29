# Welcome to RG (ReliableGenome)

ReliableGenome/RG is a method for partitioning genomes into high and low concordance regions with respect to a set of surveyed VC pipelines. RG integrates variant call sets created by multiple pipelines from arbitrary numbers of input datasets and interpolates expected concordance for genomic regions without data resulting in a genome-wide concordance score.
Ultimately, genomic regions of high/low concordance are calculated from this genome-wide signal.	

## RG genomic partitions

* A genomic partition calculated from 219 deep WGS alignments can be [data/public/20160329_RG-win1000-score1_-3-RELIABLE-above0.5.bed.gz](downloaded here). The partition contains 2,209,778 concordance intervals located on chromosomes 1,..,22,X,Y.
** Reference sequence: hs37d5 (GRCh37 + decoy)
** Read mapper: bwa + stampy
** Variant callers: GATK HaplotypeCaller, platypus, samtools
** Params: w<sup>c</sup> = 1; w<sup>d</sup> = -3; t<sub>c</sub>=t<sub>d</sub>=0.5; window size x=1000



## prerequisites

* maven3 
* jdk 1.7+

## building and running RG 

* Run the bash script build.sh
* If everything works fine you will find a standalone executable JAR in the bin directory and a library jar (containing only the classes from the RG source tree) in the target directory


## Running RG

You can run RG via `java -jar bin/wtchg-rg-1.0.jar` which will print basic usage information. 
Use `java -Xmx12g -jar ...` to run RG with more dedicated heap space (recommended).

### RG JOIN

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

### RG CALC

To calculate the genome-wide concordance score signal, run:
`java -jar bin/wtchg-rg-1.0.jar CalcreliabilitySignals calc`
to get detailed usage information.

Usage example:
`java -Xmx12g -jar bin/wtchg-rg-1.0.jar calc 
    -o <OUTDIR> 
    -w 1000 \
    -scoringSchema 1,-3 
    -thresholds 0.5,0.5 
    -dontCheckSort -v`

## Test data

Find some test VCF files that are ready to JOIN in [data/public/vcf/](data/public/vcf/)
