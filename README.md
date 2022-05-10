# GATK germline CNV caller (DNAnexus Platform App)

Dx wrapper to run the GATK germlineCNVcaller, generate gCNV.bed file for copy ratio visualisation and summarise CNV calls across samples on a run.
This is the source code for an app that runs on the DNAnexus Platform.

## What does this app do?
Calls copy number variants (CNV) from a minimum of 30 samples based on variation in read depth of exome sequencing data.

Generates a gCNV.bed file for the visualisation of CNV calls in IGV.js.

Generates run-level statistics and summary of the calls in all samples.

## What are typical use cases for this app?
This app may be executed as standalone or as part of an analysis pipeline.

Input files (preprocessed and intervals) come from the prepIntervals app and input bam and bai files come from the Sentieon app of the NGS data processing workflow (eg Dias single).

## What data are required for this app to run?
This app requires:
* at least 30 sample.bam and sample.bai files from the same sequencing run (array of files)
* preprocessed.interval_list specifying the intervals where CNVs should be called (-iinterval_list)
* a corresponding annotation.tsv that has GC content and optionally other information about each of the intervals in the list (-iannotation_tsv)
* OPTIONAL to provide a prior probability.tsv for CNV calling (applet has a built-in default)

Parameters to calculate coverage, filter out low quality calls and call variants:
* low read depth threshold - discard intervals with low coverage across the majority of samples
* low percentage threshold - if interval is covered below threshold in samples above this percentage then filter out that interval
* minimum and maxinmum GC content of interval (-iminGC, -imaxGC)
* (minimum and maximum mappability of interval) - depending if annotation.tsv has this information
* (minimum and maximum segmental duplication of interval) - depending if annotation.tsv has this information

## What does this app output?
* a pair of _intervals.vcf (genotype of every interval in the input bed) and _segments.vcf (consecutive intervals with the same CNV status are merged) for each sample
* a summary and statistics of all CNV calls in all samples on the run
* gCNV.bed file for copy ratio visualisation in IGV.js

## Dependencies
The app requires a tar.gz of the Docker image for GATK 4.2+ as an input. Htslib and bedtools are bundled with the app.
Deb files for [parallel](https://ftp.gnu.org/gnu/parallel/) and its dependency [sysstat](http://sebastien.godard.pagesperso-orange.fr/download.html) are bundled with the app in /resources.

### This app was made by East GLH
