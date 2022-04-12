#!/bin/Python3.8
"""
    Script to parse CNV vcf file and annotate each call with exon information.

    Expected inputs:
        - folder of .vcf files
        - filename of the exons.bed for annotation with expected columns
            ["chrom", "start", "end", "gene", "transcript", "exon_num"]

    Sophie Ratkai 220317
"""
import os
import sys

import pandas as pd
import vcf
import pybedtools as bedtools


def vcf2df(vcf_file):
    CNVcalls = pd.DataFrame()

    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        CNV_call = {}
        call_sample = record.samples[0]  # vcf.model._Call
        CNV_call['sample'] = call_sample.sample
        CNV_call['chrom'] = record.CHROM  # str
        CNV_call['start'] = record.POS  # int
        CNV_call['end'] = record.INFO["END"]
        CNV_call['ID'] = record.ID  # str
        if record.ALT[0] is None:
            continue
        CNV_call['CNV'] = str(record.ALT[0]).strip("<>")

        format_fields = record.FORMAT.split(':')  # GT:CN:NP:QA:QS:QSE:QSS
        for i, field in enumerate(format_fields):
            CNV_call[field] = call_sample.data[i]

        CNVcalls = CNVcalls.append(CNV_call, ignore_index=True)
    return CNV_call['sample'], CNVcalls


def ID2pos(row):
    pos = row.ID.split('_')
    return pos[1], pos[2], pos[3], row.ID


def annotate(calls_df, exons_df):
    """
        Annotate each CNV call with overlapping exon information
    """
    calls_bed = bedtools.BedTool.from_dataframe(calls_df)
    exon_bed = bedtools.BedTool.from_dataframe(exons_df)

    # Intersecting calls bed file with transcript_exon information
    # requires 100% overlap on exon within call coordinates
    calls_w_exons = calls_bed.intersect(exon_bed, loj=True)
    # convert pybedtools object to dataframe
    annotated_calls_df = calls_w_exons.to_dataframe(index_col=False,
        names=["call_chrom", "call_start", "call_end", "call_ID",
        "exon_chrom", "exon_start", "exon_end", "gene", "transcript", "exon"])
    annotated_calls_df['exon_length'] = \
        annotated_calls_df['exon_end'] - annotated_calls_df['exon_start']
    # print("annotated_calls_df")
    # print(annotated_calls_df)

    calls = annotated_calls_df["call_ID"].unique().tolist()
    genes = []  # list that will become the genes column in the dfexon
    transcripts = []
    exons = []  # list that will become the exons column in the df
    lengths = []

    for call in calls:  # Depending on the number of exon annotation and
                        # whether there's annotation available at all:
        annotation_df = annotated_calls_df[annotated_calls_df["call_ID"] == call]
        annotation_df = annotation_df.reset_index(drop=True)
        # print("annotation_df for sample {}, call {}".format(
        #     calls_df["sample"][0], call))
        if len(annotation_df) == 1:  # call is a single exon
            call_gene = annotation_df["gene"][0]
            call_transcript = annotation_df["transcript"][0]
            call_exon = annotation_df["exon"][0]
            call_length = annotation_df["exon_length"][0]
        else:  # this call covers multiple exons
            call_genes = annotation_df["gene"].unique().tolist()
            if len(call_genes) == 1:  # call covers multiple exons in 1 gene
                call_gene = call_genes[0]
                call_transcript = annotation_df["transcript"].tolist()[0]
                exon_nums = sorted(annotation_df["exon"].unique().tolist())
                start_exon = str(exon_nums[0])
                end_exon = str(exon_nums[-1])
                call_exon = "-".join([start_exon, end_exon])
                call_length = sum(annotation_df['exon_length'].tolist())
            else:  # this call covers multiple exons in multiple genes
                call_gene = ', '.join(call_genes)
                call_transcripts = annotation_df["transcript"].unique().tolist()
                call_transcript = ', '.join(call_transcripts)
                call_exons = []
                exon_lengths = []
                for transcript in call_transcripts:
                    exon_annot_df = annotation_df[
                        annotation_df['transcript'] == transcript][
                            ['exon', 'exon_length']]
                    exon_nums = sorted(exon_annot_df["exon"].unique().tolist())
                    if len(exon_nums) == 1:
                        call_exons.append(str(exon_nums[0]))
                    else:
                        start_exon = str(exon_nums[0])
                        end_exon = str(exon_nums[-1])
                        call_exons.append("-".join([start_exon, end_exon]))
                    exon_lengths.append(exon_annot_df['exon_length'].sum())
                call_exon = ', '.join(call_exons)
                call_length = sum(exon_lengths)

        genes.append(call_gene)
        transcripts.append(call_transcript)
        exons.append(call_exon)
        lengths.append(call_length)

    # create annotation dataframe
    CNV_annotation = pd.DataFrame.from_dict({"ID": calls, "gene": genes,
                "transcript": transcripts, "exon": exons, "length": lengths})
    return CNV_annotation


if __name__ == "__main__":

    # Parse command inputs: folder of vcfs and an exons.tsv
    folder = sys.argv[1]
    exons_tsv = sys.argv[2]

    # Ensure that folder has vcf files, otherwise exit
    num_vcf = 0
    for file in os.listdir(folder):
        if file.endswith("_segments.vcf"):
            num_vcf += 1
    if num_vcf == 0:
        print("Directory has no segments vcfs, nothing to do, exiting now.")
        sys.exit(0)
    print("Analysis will begin of {} samples in {} folder".format(
        num_vcf, folder))

    # Load the exon annotation bed file into a dataframe
    exons_df = pd.read_csv(exons_tsv, sep='\t', names=["chrom", "start", "end",
        "gene", "transcript", "exon_num"])
    # print("exons dataframe head:")
    # print(exons_df.head())

    # Initialise dictionary that will hold CNV counts
    CNV_counts = {}
    # and a DataFrame that will hold CNV calls from all samples
    CNV_calls = pd.DataFrame()
    unique_CNVs = pd.DataFrame()
    # empty list for samples with no CNV calls
    samples_no_CNV = []

    # Loop over files in dir, read vcf in and annotate calls
    for counter, filename in enumerate(os.listdir(folder)):
        if filename.endswith("_segments.vcf"):
            # print("Annotating sample {}/{}".format(counter + 1, num_vcf))
            vcf_file = os.path.join(folder, filename)
            sample_name, sampleCNVcalls = vcf2df(vcf_file)
            if len(sampleCNVcalls) == 0:
                samples_no_CNV.append(sample_name)
                continue
            else:
                calls = sampleCNVcalls[["ID", "CNV"]].apply(tuple, axis=1)
                for call in calls:
                    if call not in CNV_counts.keys():
                        CNV_counts[call] = {}
                        CNV_counts[call]['count'] = 1
                        CNV_counts[call]['samples'] = sample_name
                    else:
                        CNV_counts[call]['count'] += 1
                        CNV_counts[call]['samples'] += ", " + sample_name
            CNV_calls = CNV_calls.append(sampleCNVcalls,
                ignore_index=True, sort=False).astype(
                    {'start': 'int64', 'end': 'int64', 'CN': 'int64'})
        else:
            pass

    CNV_counts_df = pd.DataFrame.from_dict(CNV_counts).transpose()
    CNV_counts_df.reset_index(inplace=True)
    CNV_counts_df.rename(
        columns={'level_0': 'ID', 'level_1': 'CNV'}, inplace=True
    )
    unique_CNVs[["chrom", "start", "end", "ID"]] = CNV_counts_df.apply(
        lambda row: ID2pos(row), axis=1, result_type='expand'
    )

    # Add annotation: "gene","transcript","exon","length" columns
    CNV_annotation = annotate(unique_CNVs, exons_df)
    annotated_CNVs = CNV_calls.merge(CNV_annotation, on='ID')
    annotated_CNVs.sort_values(by=['sample', 'chrom', 'start'],
                                inplace=True, ignore_index=True)

    CNV_counts_df = CNV_counts_df.merge(CNV_annotation, on='ID')
    CNV_counts_df.sort_values(by='ID', inplace=True, ignore_index=True)

    # Write annotated calls and counts to files
    annotated_CNVs.to_csv("Annotated_CNV_summary.tsv", sep='\t', index=False,
                    encoding='utf-8', header=True)
    CNV_counts_df.to_csv("CNV_counts.tsv", sep='\t', encoding='utf-8',
        index=False)
    with open("CNV_counts.tsv", mode='a') as fh:
        fh.write(f"no CNV detected \t {len(samples_no_CNV)} \t {', '.join(samples_no_CNV)}")
