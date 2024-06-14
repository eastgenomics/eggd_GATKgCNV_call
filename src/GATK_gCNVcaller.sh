#!/bin/bash
# GATKgCNV_call
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.

# Exit at any point if there is any error and output each line as it is executed (for debugging)
set -e -x -o pipefail

main() {

    mark-section "Installing packages"
    echo "Installing packages"
    sudo dpkg -i sysstat*.deb
    sudo dpkg -i parallel*.deb
    cd packages
    pip install -q pytz-* python_dateutil-* pysam-* numpy-* pandas-* pybedtools-* PyVCF-*
    cd ..

    # Load the GATK docker image
    mark-section "Loading GATK Docker image"
    dx download "$GATK_docker" -o GATK.tar.gz
    docker load -i GATK.tar.gz

    # Parse the image ID from the list of docker images
    # need to export variables so they're available to parallel
    export GATK_image=$(docker images --format="{{.Repository}} {{.ID}}" | grep "^broad" | cut -d' ' -f2)
    export $CollectReadCounts_args
    export $PostprocessGermlineCNVCalls_args

    ## Create folder to collect input files:
    mkdir inputs

    # Prior probabilities tsv
    # file can be provided as input or a default is used bundled with the app
    if [[ ! -z $prior_prob ]]
    then
        echo "Prior prob file is provided as '$prior_prob'"
        dx download "$prior_prob" -o inputs/prior_prob.tsv
    else
        mv prior_prob.tsv inputs/
    fi

    cd inputs
    mkdir beds
    mark-section "Downloading interval files"
    # Intervals file (preprocessed bed from GATK_prep)
    dx download "$interval_list" -o beds/preprocessed.interval_list
    # Annotation tsv (from GATK_prep)
    dx download "$annotation_tsv" -o beds/annotated_intervals.tsv

    mkdir bams
    ## Download all input bam and bai files
    mark-section "Downloading input bam & bai files"
    echo ${bambais[@]} | jq -r '.["$dnanexus_link"]' | xargs -n1 -P$(nproc --all) dx download --no-progress -o bams/

    cd /home/dnanexus
    echo "All input files have downloaded to inputs/"

    # Optional to hold job after downloading all input files
    if [ "$debug_fail_start" == 'true' ]; then exit 1; fi

    # 1. Run CollectReadCounts:
    # takes one bam (and its index) file at a time along with the targets.interval_list
    echo "Running CollectReadCounts for all input bams"
    mark-section "CollectReadCounts"
    mkdir inputs/base_counts
    find inputs/bams/ -name "*.bam" | parallel -I filename --max-args 1 --jobs 8 \
    'sample_file=$( basename filename ); \
    sample_name="${sample_file%.bam}"; \
    echo $sample_name; \
    /usr/bin/time -v sudo docker run -v /home/dnanexus/inputs:/data $GATK_image gatk CollectReadCounts \
    -I /data/bams/${sample_file} \
    -L /data/beds/preprocessed.interval_list -imr OVERLAPPING_ONLY \
    ${CollectReadCounts_args} \
    -O /data/base_counts/${sample_name}_basecount.hdf5'

    # prepare a batch_input string that has all sample_basecount.tsv file as an input
    # and another that has the fileIDs for cnv_call subjobs
    batch_input=""
    #basecount_ids=""
    for base_count in inputs/base_counts/*_basecount.hdf5; do
        sample_file=$( basename $base_count )
        batch_input+="--input /data/base_counts/${sample_file} "
        #basecount_id=$( dx upload --brief $base_count )
        #basecount_ids+="-ibasecounts=$basecount_id "
    done

    # 2. Run FilterIntervals:
    # filters out low coverage or not uniquely mappable regions
    echo "Running FilterIntervals for the preprocessed intervals with sample basecounts"
    mark-section "FilterIntervals"
    /usr/bin/time -v docker run -v /home/dnanexus/inputs:/data $GATK_image gatk FilterIntervals \
        -L /data/beds/preprocessed.interval_list -imr OVERLAPPING_ONLY \
        --annotated-intervals /data/beds/annotated_intervals.tsv \
        $batch_input  $FilterIntervals_args \
        -O /data/beds/filtered.interval_list

    echo "Identifying excluded intervals from CNV calling on this run"
    # Convert interval_list to BED files
    docker run -v /home/dnanexus/inputs:/data $GATK_image gatk IntervalListToBed \
        -I /data/beds/preprocessed.interval_list \
        -O /data/beds/preprocessed.bed
    docker run -v /home/dnanexus/inputs:/data $GATK_image gatk IntervalListToBed \
        -I /data/beds/filtered.interval_list \
        -O /data/beds/filtered.bed

    # Identify regions that are in preprocessed but not in filtered, ie the excluded regions
    bedtools intersect -v -a inputs/beds/preprocessed.bed -b inputs/beds/filtered.bed > excluded_intervals.bed

    # 3. Run DetermineGermlineContigPloidy:
    # takes the base count tsv-s from the previous step, optional target_bed, and a contig plody priors tsv
    # outputs a ploidy model and ploidy-calls for each sample
    echo "Running DetermineGermlineContigPloidy for the calculated basecounts"
    mark-section "DetermineGermlineContigPloidy"
    mkdir inputs/ploidy-dir
    /usr/bin/time -v docker run -v /home/dnanexus/inputs:/data $GATK_image gatk DetermineGermlineContigPloidy \
        -L /data/beds/filtered.interval_list -imr OVERLAPPING_ONLY \
        $DetermineGermlineContigPloidy_args \
        $batch_input \
        --contig-ploidy-priors /data/prior_prob.tsv \
        --output-prefix ploidy \
        -O /data/ploidy-dir
    
    # Make batch list of ploidy fileIDs for cnv_call subjobs
    #ploidycall_ids=""
    #for ploidy_call in inputs/ploidy-dir/*; do
    #    ploidycall_id=$( dx upload --brief $ploidy_call )
    #    ploidycall_ids+="-iploidycalls=$ploidycall_id "
    #done

    dx upload -rp /home/dnanexus/inputs/ploidy-dir
    dx upload -rp /home/dnanexus/inputs/base_counts

    # 4. Run GermlineCNVCaller:
    # takes the base count tsv-s, target bed and contig ploidy calls from the previous steps
    # outputs a CNVcalling model and copy ratio files for each sample
    echo "Running GermlineCNVCaller for the calculated basecounts using the generated ploidy file"
    mark-section "GermlineCNVCaller"
    mkdir inputs/gCNV-dir
    tsv=/home/dnanexus/inputs/beds/annotated_intervals.tsv
    ints=/home/dnanexus/inputs/beds/filtered.interval_list
    if [ $split_by_chromosome == 'true' ]; then
        # make new intervals files (split by chromosome) & initiate cnv calling in subjobs
        chroms=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
        cnv_call_jobs=()
        for i in ${chroms[@]}; do
            echo "Chromosome $i"
            chr_tsv=/home/dnanexus/inputs/beds/chr"$i"_annotation.tsv
            chr_ints=/home/dnanexus/inputs/beds/chr"$i".interval_list
            grep ^@ $tsv > $chr_tsv; grep ^CONTIG $tsv >> $chr_tsv; grep -P "^$i\t" $tsv >> $chr_tsv
            grep ^@ $ints > $chr_ints; grep -P "^$i\t" $ints >> $chr_ints
            chr_tsv=$( dx upload --brief $chr_tsv )
            chr_ints=$( dx upload --brief $chr_ints )
            job_name='"$i"_cnv_call'
            command="dx-jobutil-new-job call_cnvs \
                -iannotation_tsv='$chr_tsv' -iinterval_list='$chr_ints' \
                -iGATK_docker='$GATK_docker' \
                -iGermlineCNVCaller_args='$GermlineCNVCaller_args' \
                --name $job_name"
            cnv_call_jobs+=($(eval $command))
        done
    else
        command="dx-jobutil-new-job call_cnvs -iannotation_tsv=$tsv -iinterval_list=$ints"
        cnv_call_jobs+=($(eval $command))
    fi

    # Wait for all subjobs to finish before grabbing outputs
    dx wait "${cnv_call_jobs[@]}"

    # 5. Run PostprocessGermlineCNVCalls:
    # takes CNV-model in, spits vcfs out
    echo "Running PostprocessGermlineCNVCalls"
    mark-section "PostprocessGermlineCNVCalls"
        # Required Arguments for 4.2.5.0: (4.2 onwards)
        # --calls-shard-path <File>     List of paths to GermlineCNVCaller call directories.  This argument must be specified atleast once. Required. 
        # --contig-ploidy-calls <File>  Path to contig-ploidy calls directory (output of DetermineGermlineContigPloidy). Required. 
        # --model-shard-path <File>     List of paths to GermlineCNVCaller model directories.  This argument must be specified atleast once. Required. 
        # --output-denoised-copy-ratios <File> Output denoised copy ratio file.  Required. 
        # --output-genotyped-intervals <File>  Output intervals VCF file.  Required. 
        # --output-genotyped-segments <File> Output segments VCF file.  Required. 

    mkdir inputs/vcfs
    # command finds sample data based on an arbitrary index which needs to be passed to parallel
    # index is created based on the number of input bams
    # triple colon at the end is the parallel way to provide an array of integers
    sample_num=$(ls inputs/bams/*.bam | wc -l)
    index=$(expr $sample_num - 1)
    parallel --jobs 8 '/usr/bin/time -v docker run -v /home/dnanexus/inputs:/data $GATK_image \
        gatk PostprocessGermlineCNVCalls \
        --sample-index {} \
        ${PostprocessGermlineCNVCalls_args} \
        --autosomal-ref-copy-number 2 \
        --allosomal-contig X \
        --allosomal-contig Y \
        --contig-ploidy-calls /data/ploidy-dir/ploidy-calls \
        --calls-shard-path /data/gCNV-dir/CNV-calls \
        --model-shard-path /data/gCNV-dir/CNV-model \
        --output-genotyped-intervals /data/vcfs/sample_{}_intervals.vcf \
        --output-genotyped-segments /data/vcfs/sample_{}_segments.vcf \
        --output-denoised-copy-ratios /data/vcfs/sample_{}_denoised_copy_ratios.tsv \
    ' ::: $(seq 0 1 $index)

    # 6. Rename output vcf files based on the sample they contain information about
    find inputs/vcfs/ -name "*_segments.vcf" | parallel -I{} --max-args 1 --jobs 8 ' \
        sample_file=$( basename {} ); file_name="${sample_file%_segments.vcf}"; \
        sample_name=$(bcftools view {} -h | tail -n 1 | cut -f 10 ); \
        mv inputs/vcfs/$file_name"_denoised_copy_ratios.tsv" inputs/vcfs/$sample_name"_denoised_copy_ratios.tsv"; \
        mv inputs/vcfs/$file_name"_segments.vcf" inputs/vcfs/$sample_name"_segments.vcf"; \
        mv inputs/vcfs/$file_name"_intervals.vcf" inputs/vcfs/$sample_name"_intervals.vcf" \
    '

    echo "CNV calling finished successfully"

    ## Create output directories
    vcf_dir=out/result_files/CNV_vcfs && mkdir -p ${vcf_dir}
    summary_dir=out/result_files/CNV_summary && mkdir -p ${summary_dir}
    vis_dir=out/result_files/CNV_visualisation && mkdir -p ${vis_dir}

    # and move result files into outdir to be uploaded
    mv inputs/vcfs/*.vcf ${vcf_dir}/
    mv excluded_intervals.bed ${summary_dir}/$run_name"_excluded_intervals.bed"

    mark-section "Creating copy ratio visualisation files"
    # 7. Generate gcnv bed files from copy ratios for visualisation in IGV
    echo "Generating gcnv bed files for all sample copy ratios"
    denoised_copy_ratio_files=$(find inputs/vcfs/ -name "*_denoised_copy_ratios.tsv")
    python3 generate_gcnv_bed.py --copy_ratios "$denoised_copy_ratio_files" -s \
    --run "$run_name"

    mv ./"$run_name"*.gcnv.bed.gz* "${summary_dir}"/
    mv ./*.gcnv.bed.gz* "${vis_dir}"/

    echo "All scripts finished successfully, uploading output files to dx"
    mark-section "Uploading outputs"
    if [ "$debug_fail_end" == 'true' ]; then exit 1; fi

    # Upload output files
    dx-upload-all-outputs --parallel

}

call_cnvs() {
    dx-download-all-inputs

    mkdir -p /home/dnanexus/in/basecounts
    dx download $( dx find data --project $DX_WORKSPACE_ID --name *hdf5 --brief )
    mv *hdf5 /home/dnanexus/in/basecounts/

    # Make basecount batch string
    batch_input=""
    for base_count in /home/dnanexus/in/basecounts/*_basecount.hdf5; do
        sample_file=$( basename $base_count )
        batch_input+="--input /data/basecounts/${sample_file} "
    done

    # Load the GATK docker image
    mark-section "Loading GATK Docker image"
    docker load -i /home/dnanexus/in/GATK_docker/GATK*.tar.gz

    # Declare ENV variable for command
    export GATK_image=$(docker images --format="{{.Repository}} {{.ID}}" | grep "^broad" | cut -d' ' -f2)

    # Run CNV caller
    /usr/bin/time -v docker run -v /home/dnanexus/in/:/data/ $GATK_image gatk GermlineCNVCaller \
        -L /data/interval_list/*.interval_list -imr OVERLAPPING_ONLY \
        --annotated-intervals /data/annotation_tsv/*_annotation.tsv \
        --run-mode COHORT \
        $GermlineCNVCaller_args \
        $batch_input \
        --contig-ploidy-calls /home/dnanexus/in/ploidy-dir/ploidy-calls/ \
        --output-prefix CNV \
        -O /data/gCNV-dir
    
    # Upload outputs back to parent
    dx-upload-all-outputs
}
