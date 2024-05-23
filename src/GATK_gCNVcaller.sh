#!/bin/bash
# GATKgCNV_call
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.

# Exit at any point if there is any error and output each line as it is executed (for debugging)
set -e -x -o pipefail

run_cnv_calling() {

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

    mkdir inputs/beds
    mark-section "Downloading interval files"
    # Intervals file (preprocessed bed from GATK_prep)
    dx download "$interval_list" -o inputs/beds/preprocessed.interval_list
    # Annotation tsv (from GATK_prep)
    dx download "$annotation_tsv" -o inputs/beds/annotated_intervals.tsv

    mkdir inputs/bams
    ## Download all input bam and bai files
    mark-section "Downloading input bam & bai files"
    echo ${bambais[@]} | jq -r '.["$dnanexus_link"]' | xargs -n1 -P$(nproc --all) dx download --no-progress -o inputs/bams/

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
    batch_input=""
    for base_count in inputs/base_counts/*_basecount.hdf5; do
        sample_file=$( basename $base_count )
        batch_input+="--input /data/base_counts/${sample_file} "
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

    # 4. Run GermlineCNVCaller:
    # takes the base count tsv-s, target bed and contig ploidy calls from the previous steps
    # outputs a CNVcalling model and copy ratio files for each sample
    echo "Running GermlineCNVCaller for the calculated basecounts using the generated ploidy file"
    mark-section "GermlineCNVCaller"
    mkdir inputs/gCNV-dir
    /usr/bin/time -v docker run -v /home/dnanexus/inputs:/data $GATK_image gatk GermlineCNVCaller \
        -L /data/beds/filtered.interval_list -imr OVERLAPPING_ONLY \
        --annotated-intervals /data/beds/annotated_intervals.tsv \
        --run-mode COHORT \
        $GermlineCNVCaller_args \
        $batch_input \
        --contig-ploidy-calls /data/ploidy-dir/ploidy-calls/ \
        --output-prefix CNV \
        -O /data/gCNV-dir

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

main() {

    # Make Intervals & Annotation files arrays to allow fileIDs to be paired up as inputs to subjobs
    declare -A interval_file_array
    for i in ${!interval_list[@]}; do
        prefix=$( dx describe --name "${interval_list[$i]}" | cut -d "." -f1 )
        id=$( dx describe --json "${interval_list[$i]}" | jq -r '.id' )
        interval_file_array+=( ["$prefix"]=$id )
    done

    declare -A annotation_tsv_array
    for i in ${!annotation_tsv[@]}; do
        prefix=$( dx describe --name "${annotation_tsv[$i]}" | cut -d "_" -f1 )
        id=$( dx describe --json "${annotation_tsv[$i]}" | jq -r '.id' )
        annotation_tsv_array+=( ["$prefix"]=$id )
    done

    # Make bambais string
    bambais_str=$( echo ${bambais[@]} | jq -r '.["$dnanexus_link"]' | sed s/file/\ -ibambais\=file/g )
    
    # FileID for docker image
    GATK_docker=$( echo ${GATK_docker[@]} | jq -r '.["$dnanexus_link"]' )

    # Set off subjobs per intervals file (usually per chromosome if split)
    cnv_call_jobs=()
    for i in ${!interval_file_array[@]}; do
        interval_list="${interval_file_array[$i]}"
        annotation_tsv="${annotation_tsv_array[$i]}"
        job_name='"$i"_cnv_call'
        command="dx-jobutil-new-job run_cnv_calling -iinterval_list=$interval_list \
                    -iannotation_tsv=$annotation_tsv $bambais_str -iGATK_docker=$GATK_docker \
                    -irun_name=$run_name --instance-type mem2_ssd1_v2_x16 --name $job_name"
        cnv_call_jobs+=($(eval $command))
    done

    # Wait for all subjobs to finish before grabbing outputs
    dx wait "${cnv_call_jobs[@]}"

    # Get the output from the cnv_call jobs
    echo "cnv_call jobs:"
    echo "${cnv_call_jobs[@]}"

    # Download all output files to head job instance (here)
    for job in ${cnv_call_jobs[@]}; do
        mkdir -p outputs/${job}
        for i in $(dx describe ${job}:result_files --json --multi | jq -r '.[] | .id'); do
            dx download $i -o outputs/${job}/
        done
    done

    # Merge all the files so there's one per chromosome per sample
    # Interval & Segments VCFs
    for vcf in $(find outputs/ -name *.vcf); do
        bcftools view $vcf -Oz -o "$vcf".gz
        bcftools index "$vcf".gz
    done
    mkdir -p out/result_files
    for unique_vcf_name in $( basename -s .vcf.gz $(find outputs/ -name "*.vcf.gz") | sort | uniq ); do
        bcftools merge outputs/*/"$unique_vcf_name".vcf.gz -Ov -o out/result_files/"$unique_vcf_name".vcf
    done

    # Copy Ratio BED files
    # Concatenate the files
    for unique_bed_name in $( basename -s .gcnv.bed.gz $(find outputs/ -name "*.gcnv.bed.gz") | sort | uniq ); do
        for bedfile in $( find outputs/ -name "$unique_bed_name".gcnv.bed.gz ); do
            echo "$( zcat $bedfile | head -n 2 )" >> "$unique_bed_name"_headers.txt
            zcat $bedfile | tail -n +3 >> "$unique_bed_name"_merged.bed
        done
        # Check header is in the same order in all files
        if (( $(cat "$unique_bed_name"_headers.txt | sort | uniq | wc -l) != 2 )); then 
            echo "HEADERS DO NOT MATCH"
            exit 1
        fi
        # Add header & copy to result directory
        cat "$unique_bed_name"_headers.txt| sort | uniq | tail -n 1 > out/result_files/"$unique_bed_name".gcnv.bed.gz
        cat "$unique_bed_name"_headers.txt| sort | uniq | head -n 1 >> out/result_files/"$unique_bed_name".gcnv.bed.gz
        cat "$unique_bed_name"_merged.bed >> out/result_files/"$unique_bed_name".gcnv.bed.gz
        tabix out/result_files/"$unique_bed_name".gcnv.bed.gz
    done

    # Excluded Interbals BED
    excluded_intervals=$( cat outputs/*/"$run_name"_excluded_intervals.bed > out/result_files/"$run_name"_excluded_intervals.bed)

    echo "Specifying output files"
    dx-upload-all-outputs --parallel
    #for job in ${cnv_call_jobs[@]}; do
    #    dx-jobutil-add-output result_files ${job}:result_files --class=array:jobref
    #    dx-jobutil-add-output stats_txts ${job}:stats_txt --class=array:jobref
    #done

}
