#!/bin/bash

# Exit at any point if there is any error and output each line as it is executed (for debugging)
set -exo pipefail

# prefixes all lines of commands written to stdout with datetime
PS4='\000[$(date)]\011'
export TZ=Europe/London

# set frequency of instance usage in logs to 10 seconds
kill $(ps aux | grep pcp-dstat | head -n1 | awk '{print $2}')
/usr/bin/dx-dstat 10

# control how many operations to open in parallel for download / upload, set one per CPU core
PROCESSES=$(nproc --all)

main() {
    mark-section "Installing packages and loading GATK Docker image"
    sudo dpkg -i sysstat*.deb
    sudo dpkg -i parallel*.deb
    sudo -H python3 -m pip install --no-index --no-deps packages/*

    dx download "$GATK_docker" -o GATK.tar.gz
    docker load -i GATK.tar.gz

    # create valid empty JSON file for job output, fixes https://github.com/eastgenomics/eggd_tso500/issues/19
    echo "{}" > job_output.json

    # Parse the image ID from the list of docker images
    # need to export variables (if set) so they're available to parallel
    export GATK_image=$(docker images --format="{{.Repository}} {{.ID}}" | grep "^broad" | cut -d' ' -f2)
    if [[ -n "$CollectReadCounts_args" ]]; then export $CollectReadCounts_args; fi
    if [[ -n  "$PostprocessGermlineCNVCalls_args" ]]; then export $CollectReadCounts_args; fi

    _download_parent_job_inputs

    # Optional to hold job after downloading all input files
    if [ "$debug_fail_start" == 'true' ]; then exit 1; fi

    _call_GATK_CollectReadCounts
    _call_GATK_FilterIntervals
    _call_GATK_IntervalListToBed

    # Identify regions that are in preprocessed but not in filtered, ie the excluded regions
    bedtools intersect -v -a inputs/beds/preprocessed.bed -b inputs/beds/filtered.bed > excluded_intervals.bed

    _call_GATK_DetermineGermlineContigPloidy
    _call_GATK_GermlineCNVCaller
    _call_GATKPostProcessGermlineCNVCalls

    _call_generate_gcnv_bed
    _format_output_files_for_upload

    if [ "$debug_fail_end" == 'true' ]; then exit 1; fi

    _upload_final_output
}


_download_parent_job_inputs() {
    : '''
    Download all required files for the parent job
    '''
    mark-section "Downloading inputs"
    mkdir -p inputs/beds inputs/bams

    # Prior probabilities tsv
    # file can be provided as input or a default is used bundled with the app
    if [[ -n $prior_prob ]]
    then
        echo "Prior probabilities file is provided as ${prior_prob}"
        dx download "$prior_prob" -o inputs/prior_prob.tsv
    else
        mv prior_prob.tsv inputs/
    fi

    # Intervals file (preprocessed bed from GATK_prep) and annotation tsv (from GATK_prep)
    dx download "$interval_list" -o inputs/beds/preprocessed.interval_list
    dx download "$annotation_tsv" -o inputs/beds/annotated_intervals.tsv

    SECONDS=0
    echo "${bambais[@]}" | jq -r '.["$dnanexus_link"]' \
        | xargs -n1 -P $PROCESSES dx download --no-progress -o inputs/bams/

    duration=$SECONDS
    total_size=$(du -sh /home/dnanexus/inputs/bams | cut -f1)
    total_files=$(find inputs/bams -type f | wc -l)

    echo "Downloaded ${total_files} bam/bai files (${total_size}) in $(($duration / 60))m$(($duration % 60))s"
    echo "All input files have downloaded to inputs/"
}

_format_output_files_for_upload() {
    : '''
    Rename and move around final output files for upload
    '''
    mark-section "Formatting output files for upload"

    vcf_dir=out/result_files/CNV_vcfs
    summary_dir=out/result_files/CNV_summary
    vis_dir=out/result_files/CNV_visualisation

    mkdir -p $vcf_dir $summary_dir $vis_dir

    mv inputs/vcfs/*.vcf ${vcf_dir}/
    mv excluded_intervals.bed ${summary_dir}/$run_name"_excluded_intervals.bed"
    mv ./"$run_name"*.gcnv.bed.gz* "${summary_dir}"/
    mv ./*.gcnv.bed.gz* "${vis_dir}"/

    total_files=$(find out/ -type f | wc -l)

    echo "Files moved ready to upload, ${total_files} to upload"
}

_upload_final_output() {
    : '''
    Upload the final app output files

    This is fairly quick, so using built in parallel upload
    '''
    mark-section "Uploading final output"

    SECONDS=0
    dx-upload-all-outputs --parallel
    duration=$SECONDS

    echo "Uploaded files in $(($duration / 60))m$(($duration % 60))s"
}

_upload_single_file() {
    : '''
    Uploads single file with dx upload and associates uploaded
    file ID to specified output field

    Arguments
    ---------
        1 : str
            path and file to upload
        2 : str
            app output field to link the uploaded file to
        3 : bool
            (optional) controls if to link output file to job output spec
    '''
    local file=$1
    local field=$2
    local link=$3

    # remove any parts of /home, /dnanexus, /out and /inputs from the upload path
    local remote_path=$(sed 's/\/home//;s/\/dnanexus//;s/\/out//;s/\/inputs//' <<< "$file")

    file_id=$(dx upload "$file" --path "$remote_path" --parents --brief)

    if [[ "$link" == true ]]; then
        dx-jobutil-add-output "$field" "$file_id" --array
    fi
}

_call_GATK_CollectReadCounts() {
    : '''
    Call GATK CollectReadCounts

    Called in parallel, providing each BAM and its index along with the provided targets.interval_list.

    Outputs files to /home/dnanexus/inputs/basecounts
    '''
    mark-section "Running CollectReadCounts for all input bams"
    mkdir inputs/basecounts

    SECONDS=0
    find inputs/bams/ -name "*.bam" | parallel -I filename --max-args 1 --jobs $PROCESSES \
        'sample_file=$( basename filename ); \
        sample_name="${sample_file%.bam}"; \
        sudo docker run -v /home/dnanexus/inputs:/data \
        "$GATK_image" gatk CollectReadCounts \
        -verbosity WARNING \
        -I /data/bams/${sample_file} \
        -L /data/beds/preprocessed.interval_list \
        -imr OVERLAPPING_ONLY \
        '"$CollectReadCounts_args"' \
        -O /data/basecounts/${sample_name}_basecount.hdf5'

    duration=$SECONDS
    echo "CollectReadCounts completed in $(($duration / 60))m$(($duration % 60))s"
}

_call_GATK_FilterIntervals() {
    : '''
    Call GATK FilterIntervals to filter out low coverage / non unique mapped regions
    '''
    mark-section "Running FilterIntervals for the preprocessed intervals with sample basecounts"

    # prepare a batch_input string that has all sample_basecount.tsv file as an input
    local batch_input
    batch_input=$(find inputs/basecounts/ -type f -name '*_basecount.hdf5'  -exec basename {} \; | sed 's/^/--input /g')
    # batch_input=""
    # for base_count in inputs/basecounts/*_basecount.hdf5; do
    #     sample_file=$( basename $base_count )
    #     batch_input+="--input /data/basecounts/${sample_file} "
    # done


    SECONDS=0
    /usr/bin/time -v docker run -v /home/dnanexus/inputs:/data \
        "$GATK_image" gatk FilterIntervals \
        -L /data/beds/preprocessed.interval_list \
        -imr OVERLAPPING_ONLY \
        --annotated-intervals /data/beds/annotated_intervals.tsv \
        $batch_input  $FilterIntervals_args \
        -O /data/beds/filtered.interval_list \
        --verbosity WARNING

    duration=$SECONDS
    echo "FilterIntervals completed in $(($duration / 60))m$(($duration % 60))s"
}

_call_GATK_IntervalListToBed() {
    : '''
    Call GATK IntervalListToBed to generate bed files
    '''
    mark-section "Running IntervalListToBed to identify excluded intervals from CNV calling on this run"
    SECONDS=0
    docker run \
        -v /home/dnanexus/inputs:/data \
        "$GATK_image" gatk IntervalListToBed \
        -I /data/beds/preprocessed.interval_list \
        -O /data/beds/preprocessed.bed \
        --VERBOSITY WARNING

    docker run \
        -v /home/dnanexus/inputs:/data \
        "$GATK_image" gatk IntervalListToBed \
        -I /data/beds/filtered.interval_list \
        -O /data/beds/filtered.bed \
        --VERBOSITY WARNING

    duration=$SECONDS
    echo "IntervalListToBed for preprocessed and filtered lists completed in $(($duration / 60))m$(($duration % 60))s"
}

_call_GATK_DetermineGermlineContigPloidy() {
    : '''
    Call GATK DetermineGernlineContigPloidy to generate the ploidy model and calls for each sample

    This is the most CPU/memory intensive and longest step
    '''
    mark-section "Running DetermineGermlineContigPloidy for the calculated basecounts"
    mkdir inputs/ploidy_dir

    local batch_input
    batch_input=$(find inputs/basecounts/ -type f -name '*_basecount.hdf5'  -exec basename {} \; | sed 's/^/--input /g')
    # batch_input=""
    # for base_count in inputs/basecounts/*_basecount.hdf5; do
    #     sample_file=$( basename $base_count )
    #     batch_input+="--input /data/basecounts/${sample_file} "
    # done

    SECONDS=0
    /usr/bin/time -v docker run \
        -v /home/dnanexus/inputs:/data \
        "$GATK_image" gatk DetermineGermlineContigPloidy \
        -L /data/beds/filtered.interval_list \
        -imr OVERLAPPING_ONLY \
        $DetermineGermlineContigPloidy_args \
        $batch_input \
        --contig-ploidy-priors /data/prior_prob.tsv \
        --output-prefix ploidy \
        -O /data/ploidy_dir \
        --verbosity WARNING

    duration=$SECONDS
    echo "DetermineGermlineContigPloidy completed in $(($duration / 60))m$(($duration % 60))s"
}

_call_GATK_GermlineCNVCaller() {
    : '''
    Calls GATK GermlineCNVCaller to call the actual CNVs.

    This may run entirely in the parent job, or may be split to sub jobs by chromosome
    or by interval, dependent on if -iscatter_by_interval_count or -iscatter_by_chromosome
    are specified, respectively.

    If splitting to sub jobs, _set_off_sub_jobs is called, which holds the parent job until
    all sub jobs have completed.
    '''
    mark-section "Running GermlineCNVCaller for the calculated basecounts using the generated ploidy file"
    mkdir inputs/gCNV-dir

    local total_intervals
    total_intervals=$(grep -v ^@ /home/dnanexus/inputs/beds/filtered.interval_list | wc -l)

    local batch_input
    batch_input=$(find inputs/basecounts/ -type f -name '*_basecount.hdf5'  -exec basename {} \; | sed 's/^/--input /g')

    if [ "$scatter_by_interval_count" == 'true' -a "$scatter_count" -lt "$total_intervals" ]; then
        mark-section "Scattering intervals into sublists of approximately $scatter_count"

        /usr/bin/time -v docker run -v /home/dnanexus/inputs:/data \
            "$GATK_image" gatk IntervalListTools \
            --INPUT /data/beds/filtered.interval_list \
            --SUBDIVISION_MODE INTERVAL_COUNT \
            --SCATTER_CONTENT "$scatter_count" \
            --OUTPUT /data/scatter-dir \
            --VERBOSITY WARNING

        # Set off subjobs, will be held here until all complete
        _set_off_subjobs
    elif [ "$scatter_by_chromosome" == "true" ]; then
        echo "Scattering intervals by chromosome"
        local chromosomes
        local ints

        chromosomes=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
        ints=/home/dnanexus/inputs/beds/filtered.interval_list

        for i in "${chromosomes[@]}"; do
            echo "Checking chromosome $i"
            mkdir -p /home/dnanexus/inputs/scatter-dir/chr"$i"

            local chr_ints
            chr_ints=/home/dnanexus/inputs/scatter-dir/chr"$i"/scattered.interval_list

            # Skip chromosome if no intervals present
            if [[ -z $(grep -P '^'$i'\t' $ints | head -n1) ]]; then
                echo "No intervals found for Chromosome $i, skipping..."
            else
                # Collect header & relevant lines for current chromosome
                grep ^@ $ints > $chr_ints
                grep -P '^'$i'\t' $ints >> $chr_ints
            fi

        done

        echo "Completed collecting intervals for all chromosomes"

        # Set off subjobs, will be held here until all complete
        _set_off_subjobs
    else
        # Set off cnv_calling together in the parent job
        /usr/bin/time -v docker run -v /home/dnanexus/inputs:/data \
            "$GATK_image" gatk GermlineCNVCaller \
            -L /data/beds/filtered.interval_list \
            -imr OVERLAPPING_ONLY \
            --annotated-intervals /data/beds/annotated_intervals.tsv \
            --run-mode COHORT \
            $GermlineCNVCaller_args \
            $batch_input \
            --contig-ploidy-calls /data/ploidy_dir/ploidy-calls/ \
            --output-prefix CNV \
            -O /data/gCNV-dir \
            --verbosity WARNING
    fi
}

_call_GATKPostProcessGermlineCNVCalls() {
    : '''
    Call GATK PostProcessGermlineCNVCalls to generate the output VCFs of CNVs from the model

    Required Arguments for 4.2.5.0: (4.2 onwards)
    --calls-shard-path <File>     List of paths to GermlineCNVCaller call directories. This argument must be specified at least once. Required.
    --contig-ploidy-calls <File>  Path to contig-ploidy calls directory (output of DetermineGermlineContigPloidy). Required.
    --model-shard-path <File>     List of paths to GermlineCNVCaller model directories. This argument must be specified at least once. Required.
    --output-denoised-copy-ratios <File> Output denoised copy ratio file.  Required.
    --output-genotyped-intervals <File>  Output intervals VCF file.  Required.
    --output-genotyped-segments <File> Output segments VCF file.  Required.
    '''
    mark-section "Running PostprocessGermlineCNVCalls"
    mkdir inputs/vcfs

    # Make batch input for model & calls shard paths
    local batch_input_postprocess
    local models

    batch_input_postprocess=""
    models=$(find inputs/gCNV -type d -name "*-model" -exec basename {} \; | cut -d'-' -f1)
    batch_input_postprocess+=$(sed 's/^/ --model-shard-path \/data\/gCNV-dir\//g; s/$/-model/g' <<< $models)
    batch_input_postprocess+=$(sed 's/^/ --calls-shard-path \/data\/gCNV-dir\//g; s/$/-calls/g' <<< $models)

    # command finds sample data based on an arbitrary index which needs to be passed to parallel
    # index is created based on the number of input bams
    # triple colon at the end is the parallel way to provide an array of integers
    local index
    index=$(expr $(find inputs/bams -type f -name '*.bam' | wc -l) - 1)

    SECONDS=0
    parallel --jobs $PROCESSES '/usr/bin/time -v docker run -v /home/dnanexus/inputs:/data $GATK_image \
        gatk PostprocessGermlineCNVCalls \
        --sample-index {} \
        '"$PostprocessGermlineCNVCalls_args"' \
        --autosomal-ref-copy-number 2 \
        --allosomal-contig X \
        --allosomal-contig Y \
        --contig-ploidy-calls /data/ploidy_dir/ploidy-calls \
        '"$batch_input_postprocess"' \
        --output-genotyped-intervals /data/vcfs/sample_{}_intervals.vcf \
        --output-genotyped-segments /data/vcfs/sample_{}_segments.vcf \
        --output-denoised-copy-ratios /data/vcfs/sample_{}_denoised_copy_ratios.tsv \
        --verbosity WARNING \
    ' ::: $(seq 0 1 $index)

    duration=$SECONDS
    echo "PostprocessGermlineCNVCalls completed in $(($duration / 60))m$(($duration % 60))s"
}

_call_generate_gcnv_bed() {
    : '''
    Call the generate_gcnv_bed.py script to generate CNV visualisation beds from the copy ratio tsvs
    '''
    mark-section "Generating gCNV copy ratio visualisation files"
    local denoised_copy_ratio_files
    denoised_copy_ratio_files=$(find inputs/vcfs/ -name "*_denoised_copy_ratios.tsv")

    python3 generate_gcnv_bed.py \
        --copy_ratios "$denoised_copy_ratio_files" \
        -s --run "$run_name"

    echo "Completed generating gCNV copy ratio visualisation files"
}

_launch_sub_jobs() {
    : '''
    Wrapper function to launch all sub jobs per chromosome / interval by calling
    dx-jobutil-new-job with the _sub_job function.

    This function will hold the parent job until all sub jobs complete.
    '''
    mark-section "Launching subjobs to run GATK GermlineCNVCaller"

    # Upload input files to workspace container so subjobs can access them
    tsv=$( dx upload --brief /home/dnanexus/inputs/beds/annotated_intervals.tsv )

    SECONDS=0
    echo "Uploading polidy and base counts for sub jobs"
    export -f _upload_single_file  # required to be accessible to xargs sub shell
    find /home/dnanexus/inputs/ploidy_dir /home/dnanexus/inputs/basecounts -type f \
        | xargs -P $(nproc --all) -n1 -I{} bash -c "_upload_single_file {} _ false"

    duration=$SECONDS
    echo "Uploaded ploidy and base count files in $(($duration / 60))m$(($duration % 60))s"

    for i in $( find /home/dnanexus/inputs/scatter-dir -name "scattered.interval_list" ); do
        job_name=$( echo $i | rev | cut -d '/' -f 2 | rev )
        ints=$( dx upload --brief $i )

        # Bump instance type up for large interval lists
        interval_num=$(grep -v ^@ $i | wc -l)
        if [ $interval_num -gt 15000 ]; then
            instance="mem1_ssd1_v2_x72"
        elif [ $interval_num -gt 10000 ]; then
            instance="mem1_ssd1_v2_x36"
        else
            instance="mem1_ssd1_v2_x16"
        fi

        dx-jobutil-new-job _call_cnvs \
            -iannotation_tsv="$tsv" \
            -iinterval_list="$ints" \
            -iGATK_docker="$GATK_docker" \
            -iGermlineCNVCaller_args="$GermlineCNVCaller_args" \
            --instance-type "$instance" \
            --extra-args='{"priority": "high"}' \
            --name "$job_name" >> job_ids
    done

    # Wait for all subjobs to finish before grabbing outputs
    SECONDS=0
    echo "$(wc -l job_ids) jobs launched, holding job until all to complete..."
    dx wait --from-file job_ids

    duration=$SECONDS
    echo "All subjobs for GermlineCNVCaller completed in $(($duration / 60))m$(($duration % 60))s"
}

_get_sub_job_output() {
    : '''
    Get the output of GerminelineCNVCaller run in each subjob and download back to the parent job.

    Should only be called once all sub jobs have completed from the parent job.
    '''
    mark-section "Downloading gCNV job output"

    echo "Waiting 30 seconds to ensure all files are hopefully in a closed state..."
    sleep 30

    SECONDS=0
    set +x  # suppress this going to the logs as its long

    # files from sub jobs will be in the container- project context of the
    # current job ($DX_WORKSPACE-id) => search here for all the files
    gCNV_files=$(dx find data --json --verbose --path "$DX_WORKSPACE_ID:/gCNV-dir")

    # turn describe output into id:/path/to/file to download with dir structure
    files=$(jq -r '.[] | .id + ":" + .describe.folder + "/" + .describe.name'  <<< $gCNV_files)

    # build aggregated directory structure and download all files
    cmds=$(for f in  $files; do \
        id=${f%:*}; path=${f##*:}; dir=$(dirname "$path"); \
        echo "'mkdir -p inputs/$dir && dx download --no-progress $id -o inputs/$path'"; done)

    echo $cmds | xargs -n1 -P${PROCESSES} bash -c

    set -x

    total=$(du -sh /home/dnanexus/inputs/gCNV-dir/ | cut -f1)
    duration=$SECONDS
    echo "Downloaded $(wc -w <<< ${files}) files (${total}) in $(($duration / 60))m$(($duration % 60))s"
}

_sub_job() {
    : '''
    App code for each sub job to run GATK GermlineCNVCaller per chromosome / interval

    The annotation tsv, interval list and GATK Docker image are provided as inputs to the app,
    but both the basecount and ploidy files are uploaded from the parent job and are available
    in the container, so are downloaded separately.
    '''
    # prefixes all lines of commands written to stdout with datetime
    PS4='\000[$(date)]\011'
    export TZ=Europe/London

    # set frequency of instance usage in logs to 10 seconds
    kill $(ps aux | grep pcp-dstat | head -n1 | awk '{print $2}')
    /usr/bin/dx-dstat 10

    # create valid empty JSON file for job output, fixes https://github.com/eastgenomics/eggd_tso500/issues/19
    echo "{}" > job_output.json

    # control how many operations to open in parallel for download / upload, set one per CPU core
    PROCESSES=$(nproc --all)

    # get chromosome name for output prefix
    name=$( cat dnanexus-job.json | jq -r '.name' )

    export -f _upload_single_file  # required to be accessible to xargs sub shell

    _sub_job_download_inputs

    # Make basecount batch string
    set +x
    batch_input=$(find in/basecounts/ -type f -name '*_basecount.hdf5'  -exec basename {} \; | sed 's/^/--input /g')
    set -x

    # Load the GATK docker image
    mark-section "Loading GATK Docker image"
    docker load -i /home/dnanexus/in/GATK_docker/GATK*.tar.gz
    export GATK_image=$(docker images --format="{{.Repository}} {{.ID}}" | grep "^broad" | cut -d' ' -f2)

    # Run CNV caller
    SECONDS=0
    /usr/bin/time -v docker run -v /home/dnanexus/in/:/data/ \
        "$GATK_image" gatk GermlineCNVCaller \
        -L /data/interval_list/$interval_list \
        -imr OVERLAPPING_ONLY \
        --annotated-intervals /data/annotation_tsv/$annotated_intervals \
        --run-mode COHORT \
        $GermlineCNVCaller_args \
        $batch_input \
        --contig-ploidy-calls /data/ploidy_dir/ploidy-calls/ \
        --output-prefix $name \
        -O /data/gCNV-dir \
        --verbosity WARNING

    duration=$SECONDS
    echo "GermlineCNVCaller completed in $(($duration / 60))m$(($duration % 60))s"

    _sub_job_upload_outputs
}

_sub_job_download_inputs() {
    : '''
    Util function for _sub_job().

    Downloads all inputs to the sub job, including:
        - interval lists
        - basecount files
        - ploidy files
    '''
    SECONDS=0
    dx-download-all-inputs --parallel

    interval_list=$( basename $( find /home/dnanexus/in/ -name '*.interval_list' ))
    annotated_intervals=$( basename $( find /home/dnanexus/in/ -name 'annotated_intervals.tsv' ))

    # Get basecount and ploidy files uploaded into container by parent job
    basecount_files=$(dx find data --json --verbose --path "$DX_WORKSPACE_ID:/basecounts")
    ploidy_files=$(dx find data --json --verbose --path "$DX_WORKSPACE_ID:/ploidy_dir")

    # turn describe output into id:/path/to/file to download with dir structure,
    # build aggregated directory structure and download all files to given directory path
    basecount_files=$(jq -r '.[] | .id + ":" + .describe.folder + "/" + .describe.name'  <<< $basecount_files)
    ploidy_files=$(jq -r '.[] | .id + ":" + .describe.folder + "/" + .describe.name'  <<< $ploidy_files)

    cmds=$(for f in $basecount_files $ploidy_files; do \
        id=${f%:*}; path=${f##*:}; dir=$(dirname "$path"); \
        echo "'mkdir -p in/$dir && dx download --no-progress $id -o in/$path'"; done)

    echo $cmds | xargs -n1 -P $PROCESSES bash -c

    duration=$SECONDS
    total_files=$(find in/ -type f | wc -l)
    total_size=$(du -sh in/ | cut -f 1)
    echo "Downloaded $total_files files ($total_size) in $(($duration / 60))m$(($duration % 60))s"
}

_sub_job_upload_outputs() {
    : '''
    Util function for _sub_job().

    Uploads required output back to container to be downloaded back to parent job for post processing
    '''
    mark-section "Uploading sub job output"
    mkdir -p out/gCNV-dir
    mv /home/dnanexus/in/gCNV-dir/$name-calls out/gCNV-dir/
    mv /home/dnanexus/in/gCNV-dir/$name-model out//gCNV-dir/

    total_files=$(find out/ -type f | wc -l)
    total_size=$(du -sh out/ | cut -f 1)

    SECONDS=0
    echo "Uploading sample output"
    find /home/dnanexus/out/ -type f | xargs -P $PROCESSES -n1 -I{} bash -c \
        "_upload_single_file {} _ false"

    duration=$SECONDS
    echo "Uploaded ${total_files} files (${total_size}) in $(($duration / 60))m$(($duration % 60))s"
}