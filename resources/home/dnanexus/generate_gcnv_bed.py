"""
    Script to generate bed file from a run of cnv calls and intervals list.

    Expected inputs:
        - denoised copy ratio tsv files of samples to visualise
        - run name prefix
    Optional input:
        - use '-s' to generate per sample bed files in addition to the run one

    Requires bgzip and tabix to be installed and on path.

    created: 211119 by Jethro Rainford
    last modified: 220325 by Sophie Ratkai
"""
import argparse
from pathlib import Path
import subprocess

import pandas as pd


def parse_args():
    """
        Parse command line arguments
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--copy_ratios', nargs='+',
        help='all copy ratio files to generate bed file from'
    )
    parser.add_argument(
        '--run', help='name of run to prefix output per run bed file with'
    )
    parser.add_argument(
        '-s', '--per_sample', action='store_true',
        help=(
            'generate per sample bed file, with just the sample highlighted '
            'and other sample names removed'
        )
    )

    args = parser.parse_args()

    # check if copy ratios are being passed from output of 'find'
    # where it is a list with all files as one string
    if len(args.copy_ratios) == 1:
        args.copy_ratios = args.copy_ratios[0].split('\n')

    # ensure interval list isn't accidentally passed as copy ratio file
    args.copy_ratios = [x for x in args.copy_ratios if 'interval_list' not in x]

    return args


def generate_bed(args):
    """
        Generates bed dataframe from multiple copy ratio files

        Returns:
        - copy_ratio_df (df): df of all samples and intervals
        - samples (list): list of sample names, read from header in files
    """
    # get intervals from the first file to use for comparing all against
    copy_ratio_df = pd.read_csv(
        args.copy_ratios[0], sep='\t', comment='@', header=0
    )
    copy_ratio_df = copy_ratio_df.iloc[:, 0:3]  # keep just chr, start and end
    copy_ratio_df.columns = ['chr', 'start', 'end']

    # read all files in and add to copy_ratio_df
    for file in args.copy_ratios:
        with open(file, 'r') as fh:
            # add sample name from first line as separate col name to
            # identify rows
            lines = fh.readlines()
            lines = [x for x in lines if x.startswith('@RG')]
            sample_name = lines[0].split('SM:')[1].replace('\n', '')

            # check we haven't read in twice
            try:
                assert sample_name not in copy_ratio_df.columns, (
                    f"Sample {sample_name} already in sample list!"
                    f"Current sample list: {copy_ratio_df.columns}"
                )
            except:
                print("Sample is already loaded, skipping it")
                continue

        # read copy ratio to df
        file_df = pd.read_csv(
            file, sep='\t', comment="@", header=0,
            names=['chr', 'start', 'end', sample_name]
        )

        # sense check regions for new file the same as what is in df
        # before adding
        assert copy_ratio_df.iloc[:, :3].equals(file_df.iloc[:, :3]), (
            f"Copy ratio file for {sample_name} has different intervals "
            "Length of intervals: {}\n".format(len(copy_ratio_df)),
            "Length of file: {}".format(len(file_df))
        )

        # add sample copy ratio to file
        copy_ratio_df[sample_name] = file_df[sample_name]

    return copy_ratio_df


def generate_per_sample_beds(copy_ratio_df):
    """
        Generates one dataframe / bed file per sample, with other samples
        anonymised

        Args:
            - copy_ratio_df (df): df of all copy ratios of all samples

        Returns
            - per_sample_dfs (list): list of tuples, with (name, df) per sample
    """
    samples = copy_ratio_df.columns.tolist()[3:]

    # calculate mean and std dev across all samples (rows)
    mean = copy_ratio_df[samples].mean(axis=1)
    std = copy_ratio_df[samples].std(ddof=1, axis=1)

    plus_std = mean + std
    minus_std = mean - std
    plus_std2 = mean + (std * 2)
    minus_std2 = mean - (std * 2)

    mean_std_df = pd.DataFrame.from_dict({'mean': mean,
        'mean_plus_std': plus_std, 'mean_plus_std2': plus_std2,
        'mean_minus_std': minus_std, 'mean_minus_std2': minus_std2
    })


    per_sample_dfs = []

    for sample in samples:

        sample_df = pd.concat(
            [copy_ratio_df[['chr', 'start', 'end', sample]].copy(),
                mean_std_df], axis='columns')

        per_sample_dfs.append((sample, sample_df))

    return per_sample_dfs


def write_outfile(copy_ratio_df, prefix, per_sample):
    """
        Write output bed files and compress with bgzip.
        Bed file with highlight_samples has random colouring of each sample
        trace to improve visibility.

        Args:
            - copy_ratio_df (df): df of all copy ratios to write
            - prefix (str): prefix for naming output file
            - per_sample (bool): controls writing header per sample or per run
    """
    outfile = "{}_copy_ratios.gcnv.bed".format(prefix)

    with open(outfile, 'w') as fh:
        # this line is needed at the beginning to tell igv.js that it is
        # a gcnv bed as it can't automatically set this track type
        if per_sample:
            # colour mapping for tracks
            # keys have to match colummn names in sample_df
            colours = {
                prefix: 'red',
                'mean': 'blue',
                'mean_plus_std': '#0B2559',
                'mean_minus_std': '#0B2559',
                'mean_plus_std2': '#36BFB1',
                'mean_minus_std2': '#36BFB1',
            }

            highlight = ' '.join([f"highlight={x};{y}" for x, y in colours.items()])
            fh.write(f"track  type=gcnv height=500 {highlight} \n")
        else:
            fh.write("track type=gcnv height=500 clickToHighlight=any \n")

    copy_ratio_df.to_csv(outfile, sep='\t', index=False, mode='a')

    # compress & index
    subprocess.call("bgzip {}".format(Path(outfile)), shell=True)
    subprocess.call("tabix -f -S 2 -b 2 -e 3 {}.gz".format(Path(outfile)), shell=True)


def main():

    args = parse_args()
    copy_ratio_df = generate_bed(args)

    # write output bed file
    write_outfile(copy_ratio_df, args.run, per_sample=False)

    if args.per_sample:
        # generating per sample bed files
        per_sample_dfs = generate_per_sample_beds(copy_ratio_df)

        for sample in per_sample_dfs:
            sample_name = sample[0]
            sample_df = sample[1]
            write_outfile(sample_df, sample_name, per_sample=True)

if __name__ == "__main__":
    main()
