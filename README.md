# Welcome to RG

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
