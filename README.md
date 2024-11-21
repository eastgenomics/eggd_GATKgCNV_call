# GATK Germline CNV Caller (DNAnexus Platform App)

DNAnexus app to run the GATK GermlineCNVCaller and generate per sample gCNV.bed file for copy ratio visualisation and summarise CNV calls across samples on a run.
This is the source code for an app that runs on the DNAnexus Platform.

## What does this app do?
Calls copy number variants (CNV) from a minimum of 30 samples based on variation in read depth of exome sequencing data.

Generates a gCNV.bed file for the visualisation of CNV calls in IGV.js.

## What are typical use cases for this app?
This app may be executed as standalone or as part of an analysis pipeline.

Input files (preprocessed and annotated target intervals) come from the [GATKgCNV_prep app](https://github.com/eastgenomics/eggd_GATKgCNV_prep) and input bam and bai files come from the Sentieon app of the NGS data processing workflow (eg Dias single).

## What data are required for this app to run?
**Required**
* `-iGATK_docker` (`file`) - Docker image of GATK
* `-ibambais` (`array:file`) - pairs of bam / bai files for each sample (minimum: 30 samples)
* `-iinterval_list` (`file`) - preprocessed.interval_list specifying the intervals where CNVs should be called
* `-iannotation_tsv` (`file`) - corresponding annotation.tsv that has GC content and optionally other information about each of the target intervals
* `-irun_name` (`str`) - name of run, used for naming run level output files


**Optional**
* `-iprior_prob` (`file`): prior probabilities for the copy number of the chromosomes, default bundled with app in [resources](https://github.com/eastgenomics/eggd_GATKgCNV_call/blob/main/resources/home/dnanexus/prior_prob.tsv)
* `-ipost_process_proc_per_sample` (`int`): number of CPU cores to allow per sample for PostProcessGermlineCNVCalls step. For larger captures it appears more efficient to allow more CPU cores per sample and process fewer in parallel (i.e setting this input to 2-4)
* `-iscatter_by_chromosome` (`bool`): controls if to split CNV calling across sub jobs by chromosome (i.e for larger captures)
* `-imax_sub_job_instance` (`str`): instance size to use for sub jobs where interval count >15,000 (default: `mem1_ssd1_v2_x36`)
* `-imid_sub_job_instance` (`str`): instance size to use for sub jobs where interval count <15,000 and >10,000 (default: `mem1_ssd1_v2_x16`)
* `-imin_sub_job_instance` (`str`): instance size to use for sub jobs where interval count <10,000 (default: `mem1_ssd2_v2_x8`)
* `-iCollectReadCounts_args` (`str`): optional command line arguments for CollectReadCounts
* `-iFilterIntervals_args` (`str`): optional command line arguments for FilterIntervals
* `-iDetermineGermlineContigPloidy_args` (`str`): optional command line arguments for DetermineGermlineContigPloidy
* `-iGermlineCNVCaller_args` (`str`): optional command line arguments for GermlineCNVCaller
* `-iPostprocessGermlineCNVCalls_args` (`str`): optional command line arguments for PostprocessGermlineCNVCalls
* `-ikeep_all_sample_traces` (`bool`): controls whether to keep all sample traces in per sample visualisation bed files, if true all other samples will be anonymised grey traces
* `-idebug_fail_start` (`bool`): automatically fail the job after inputs have been downloaded
* `-idebug_fail_end` (`bool`): automatically fail the job after all commands have finished


## What does this app output?
* `{sample_name}_intervals.vcf`: vcf with genotype of every interval in the input bed for the sample
* `{sample_name}_segments.vcf`: vcf of merged consecutive intervals with the same CNV status
* `{sample_name}_copy_ratios.gcnv.bed.gz`: bed file for copy ratio visualisation in igv.js
* `{sample_name}_copy_ratios.gcnv.bed.gz.tbi`: index for the copy ratio bed file
* `{run_name}_copy_ratios.gcnv.bed.gz`: bed file for the entire run copy ratio visualisation in igv.js
* `{run_name}_copy_ratios.gcnv.bed.gz.tbi`: index for the entire run copy ratio bed file
* `{run_name}_excluded_intervals.bed`: list of intervals excluded from CNV calling for the run (0-based)


**How to run this app**:

```bash
dx run app-eggd_GATKgCNV_call/ \
  -iGATK_docker="<GATK_docker.tar.gz>" \
  -iinterval_list="<output from prep app>" \
  -iannotation_tsv="<output from prep app>" \
  -irun_name="<name of run>" \
  -ibambais="<array of paired sample bam and index file IDs, specified once per file>"
```

> [!TIP]
> The array input for bam/bai files provided to `-ibambais` may be generated with the following:
>
> `$ dx api system findDataObjects '{"scope":{"project": "<project>","folder":"<folder>"}, "name": {"regexp": ".*bam$|.*bai$"}, "limit": 1000}' | jq -r '.results[].id' | sed 's/^/-ibambais=/g'`
>
> where `<project>` and `<folder>` are the project and path to the bam files, respectively.
>
> This will generate the following:
>```-ibambais=file-GPg03x84f1q7F1KZJvy9zVgY
> -ibambais=file-GPg054Q4YYz7Y4fzYYGYKz31
> -ibambais=file-GPg000Q4fFQx748FQ8jZbbBB
> -ibambais=file-GPfzzK049JXy1JKpGvJ5vJ7Q
> ...


## Notes

Parameters used by GATK to calculate coverage, filter out low quality intervals and call variants:
* low read depth threshold - discard intervals with low coverage across the majority of samples
* low percentage threshold - if interval is covered below threshold in samples above this percentage then filter out that interval
* minimum and maximum GC content of interval
* minimum and maximum mappability of interval - depending if annotation.tsv has this information
* minimum and maximum segmental duplication of interval - depending on whether annotation.tsv has this information


## Dependencies
The app requires a tar.gz of the Docker image for GATK 4.2+ as an input. Htslib and bedtools are bundled with the app as app assets.
Deb files for [parallel](https://ftp.gnu.org/gnu/parallel/) and its dependency [sysstat](http://sebastien.godard.pagesperso-orange.fr/download.html) are bundled with the app in resources/home/dnanexus/.

### This app was made by East GLH
