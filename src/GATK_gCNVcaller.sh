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

    # control how many operations to open in parallel for download / upload
    THREADS=$(nproc --all)

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
    find inputs/bams/ -name "*.bam" | parallel -I filename --max-args 1 --jobs $THREADS \
    'sample_file=$( basename filename ); \
    sample_name="${sample_file%.bam}"; \
    echo $sample_name; \
    /usr/bin/time -v sudo docker run -v /home/dnanexus/inputs:/data $GATK_image gatk CollectReadCounts \
    -I /data/bams/${sample_file} \
    -L /data/beds/preprocessed.interval_list -imr OVERLAPPING_ONLY \
    '"$CollectReadCounts_args"' \
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
    
    # Upload preprocess file to workspace container so subjobs can access them
    dx upload -rp /home/dnanexus/inputs/ploidy-dir
    dx upload -rp /home/dnanexus/inputs/base_counts

    # OPTIONAL - split intervals via scatter
    if [ $split_by_chromosome == 'true' ]; then
        /usr/bin/time -v docker run -v /home/dnanexus/inputs:/data $GATK_image gatk IntervalListTools \
        --INPUT /data/beds/filtered.interval_list \
        --SUBDIVISION_MODE INTERVAL_COUNT \
        --SCATTER_CONTENT 5000 \
        --OUTPUT /data/scatter-dir
    fi

    # 4. Run GermlineCNVCaller:
    # takes the base count tsv-s, target bed and contig ploidy calls from the previous steps
    # outputs a CNVcalling model and copy ratio files for each sample
    echo "Running GermlineCNVCaller for the calculated basecounts using the generated ploidy file"
    mark-section "GermlineCNVCaller"
    mkdir inputs/gCNV-dir
    tsv=$( dx upload --brief /home/dnanexus/inputs/beds/annotated_intervals.tsv )
    if [ $split_by_chromosome == 'true' ]; then
        cnv_call_jobs=()
        for i in $( find /home/dnanexus/inputs/scatter-dir -name "scattered.interval_list" ); do
            job_name=$( echo $i | rev | cut -d '/' -f 2 | rev )
            ints=$( dx upload --brief $i )
            command="dx-jobutil-new-job call_cnvs \
                -iannotation_tsv='$tsv' -iinterval_list='$ints' \
                -iGATK_docker='$GATK_docker' \
                -iGermlineCNVCaller_args='$GermlineCNVCaller_args' \
                --instance-type mem2_ssd1_v2_x8 \
                --name $job_name"
            cnv_call_jobs+=($(eval $command))
        done
    else
        # TODO add normal gatk command within parent job
        ints=$( dx upload --brief /home/dnanexus/inputs/beds/filtered.interval_list )
        command="dx-jobutil-new-job call_cnvs -iannotation_tsv=$tsv -iinterval_list=$ints"
        cnv_call_jobs+=($(eval $command))
    fi

    # Wait for all subjobs to finish before grabbing outputs
    dx wait "${cnv_call_jobs[@]}"

    # download all scatter output
    _get_gCNV_job_outputs

    # Make batch input for model & calls shard paths
    batch_input_postprocess=""
    for shard_dir in inputs/gCNV-dir/*-model; do
        prefix=$( basename $shard_dir | cut -d '-' -f1 )
        batch_input_postprocess+="--calls-shard-path /data/gCNV-dir/$prefix-calls "
        batch_input_postprocess+="--model-shard-path /data/gCNV-dir/$prefix-model "
    done

    # 5. Run PostprocessGermlineCNVCalls:
    # takes CNV-model in, spits vcfs out
    echo "Running PostprocessGermlineCNVCalls"
    mark-section "PostprocessGermlineCNVCalls"
        # Required Arguments for 4.2.5.0: (4.2 onwards)
        # --calls-shard-path <File>     List of paths to GermlineCNVCaller call directories.  This argument must be specified at least once. Required. 
        # --contig-ploidy-calls <File>  Path to contig-ploidy calls directory (output of DetermineGermlineContigPloidy). Required. 
        # --model-shard-path <File>     List of paths to GermlineCNVCaller model directories.  This argument must be specified at least once. Required. 
        # --output-denoised-copy-ratios <File> Output denoised copy ratio file.  Required. 
        # --output-genotyped-intervals <File>  Output intervals VCF file.  Required. 
        # --output-genotyped-segments <File> Output segments VCF file.  Required. 

    mkdir inputs/vcfs
    # command finds sample data based on an arbitrary index which needs to be passed to parallel
    # index is created based on the number of input bams
    # triple colon at the end is the parallel way to provide an array of integers
    sample_num=$(ls inputs/bams/*.bam | wc -l)
    index=$(expr $sample_num - 1)
    parallel --jobs $THREADS '/usr/bin/time -v docker run -v /home/dnanexus/inputs:/data $GATK_image \
        gatk PostprocessGermlineCNVCalls \
        --sample-index {} \
        '"$PostprocessGermlineCNVCalls_args"' \
        --autosomal-ref-copy-number 2 \
        --allosomal-contig X \
        --allosomal-contig Y \
        --contig-ploidy-calls /data/ploidy-dir/ploidy-calls \
        '"$batch_input_postprocess"' \
        --output-genotyped-intervals /data/vcfs/sample_{}_intervals.vcf \
        --output-genotyped-segments /data/vcfs/sample_{}_segments.vcf \
        --output-denoised-copy-ratios /data/vcfs/sample_{}_denoised_copy_ratios.tsv \
    ' ::: $(seq 0 1 $index)

    # 6. Rename output vcf files based on the sample they contain information about
    find inputs/vcfs/ -name "*_segments.vcf" | parallel -I{} --max-args 1 --jobs $THREADS ' \
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

    dx-download-all-inputs --parallel

    interval_list=$( basename $( find /home/dnanexus/in/ -name '*.interval_list' ))
    annotated_intervals=$( basename $( find /home/dnanexus/in/ -name 'annotated_intervals.tsv' ))

    # get chromosome name for output prefix
    name=$( cat dnanexus-job.json | jq -r '.name' )

    # Get basecounts
    mkdir -p /home/dnanexus/in/basecounts
    dx download $( dx find data --project $DX_WORKSPACE_ID --name *hdf5 --brief )
    mv *hdf5 /home/dnanexus/in/basecounts/

    # Get ploidy calls
    mkdir ploidy-dir
    dx download -r $DX_WORKSPACE_ID:ploidy-dir/ploidy-calls -o ploidy-dir/
    mv ploidy-dir/ /home/dnanexus/in/

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
        -L /data/interval_list/$interval_list -imr OVERLAPPING_ONLY \
        --annotated-intervals /data/annotation_tsv/$annotated_intervals \
        --run-mode COHORT \
        $GermlineCNVCaller_args \
        $batch_input \
        --contig-ploidy-calls /data/ploidy-dir/ploidy-calls/ \
        --output-prefix $name \
        -O /data/gCNV-dir
    
    # Upload outputs back to parent
    mkdir -p out/GermlineCNVCaller/gCNV-dir
    mv /home/dnanexus/in/gCNV-dir/$name-calls out/GermlineCNVCaller/gCNV-dir/
    mv /home/dnanexus/in/gCNV-dir/$name-model out/GermlineCNVCaller/gCNV-dir/
    dx-upload-all-outputs --parallel

}

_get_gCNV_job_outputs() {

    # Download output from all gCNV jobs to run final gather step
    # Modified from Jethro's _get_scatter_job_outputs function in eggd_tso500

    echo "Downloading gCNV job output"

    echo "Waiting 60 seconds to ensure all files are hopefully in a closed state..."
    sleep 60

    SECONDS=0
    set +x  # suppress this going to the logs as its long

    # files from sub jobs will be in the container- project context of the
    # current job ($DX_WORKSPACE-id) => search here for  all the files
    gCNV_files=$(dx find data --json --verbose --path "$DX_WORKSPACE_ID:/gCNV-dir")

    # turn describe output into id:/path/to/file to download with dir structure
    files=$(jq -r '.[] | .id + ":" + .describe.folder + "/" + .describe.name'  <<< $gCNV_files)

    # build aggregated directory structure and download all files
    cmds=$(for f in  $files; do \
        id=${f%:*}; path=${f##*:}; dir=$(dirname "$path"); \
        echo "'mkdir -p inputs/$dir && dx download --no-progress $id -o inputs/$path'"; done)

    echo $cmds | xargs -n1 -P${THREADS} bash -c

    set -x

    total=$(du -sh /home/dnanexus/inputs/gCNV-dir/ | cut -f1)
    duration=$SECONDS
    echo "Downloaded $(wc -w <<< ${files}) files (${total}) in $(($duration / 60))m$(($duration % 60))s"
    
}