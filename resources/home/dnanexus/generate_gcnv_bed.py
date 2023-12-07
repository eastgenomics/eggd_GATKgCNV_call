"""
Script to generate bed file from a run of cnv calls and intervals list.

Requires bgzip and tabix to be installed and on path.
"""
import argparse
from pathlib import Path
import subprocess

import pandas as pd


def parse_args():
    """Parse command line arguments"""
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
            'and run level mean and std deviations added'
        )
    )
    parser.add_argument(
        '--keep_all_samples', action='store_true', default=True,
        help=(
            'controls if to keep all sample traces in per sample files, '
            'if true all other samples will be coloured grey'
        )
    )
    parser.add_argument(
        '--hide_sample_names', action='store_true', default=True,
        help=(
            'controls if to hide all other sample labels in per sample files '
            'if --keep_all_samples has been set to true'
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


def generate_copy_ratio_df(args):
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

    # Start coordinates are currently 1-based, need to be offset by 1
    # to create a 0-based BED file for IGV visualisation
    copy_ratio_df['start'] = copy_ratio_df['start'].apply(lambda x: x-1)

    return copy_ratio_df


def generate_per_sample_copy_ratio_dfs(copy_ratio_df, keep_all_samples):
    """
        Generates one dataframe / bed file per sample, with run level mean
        and std deviations added

        Args:
            - copy_ratio_df (df): df of all copy ratios of all samples
            - keep_all_samples (bool): determines if to keep all sample
                traces in per sample files

        Returns
            - pd.DataFrame|list: either single DataFrame if keep_all_samples
                is True, or list of tuples, with (name, df) per sample
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

    if keep_all_samples:
        copy_ratio_df = pd.concat([copy_ratio_df, mean_std_df], axis='columns')
        return copy_ratio_df
    else:
        per_sample_dfs = []

        for sample in samples:

            sample_df = pd.concat(
                [copy_ratio_df[['chr', 'start', 'end', sample]].copy(),
                    mean_std_df], axis='columns')

            per_sample_dfs.append((sample, sample_df))

        return per_sample_dfs


def write_outfile(copy_ratio_df, prefix, per_sample, hide_sample_names):
    """
        Write output bed files and compress with bgzip.
        Bed file with highlight_samples has random colouring of each sample
        trace to improve visibility.

        Args:
            - copy_ratio_df (df): df of all copy ratios to write
            - prefix (str): prefix for naming output file
            - per_sample (bool): controls writing bed file header, if per sample the
                mean tracks are added, else if per run clickToHighlight is enabled
            - hide_sample_names (bool): controls if to hide other sample labels
                in per sample file where all sample traces are being added
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
            fh.write(f"track type=gcnv height=500 {highlight} \n")
        else:
            fh.write("track type=gcnv height=500 clickToHighlight=any \n")

        if not hide_sample_names:
            copy_ratio_df.to_csv(outfile, sep='\t', index=False, mode='a')
        else:
            # written all sample traces as grey lines to per sample file
            # => hide names of all but our sample
            to_keep = [
                'chr', 'start', 'end', 'mean', 'mean_plus_std', 'mean_plus_std2',
                'mean_minus_std', 'mean_minus_std2', prefix
            ]

            header = '\t'.join([
                '⠀' if x not in to_keep else x
                for x in copy_ratio_df.columns.tolist()
            ])

            fh.write(header)
            copy_ratio_df.to_csv(
                outfile, sep='\t', header=False, index=False, mode='a'
            )

    # compress & index
    subprocess.call("bgzip {}".format(Path(outfile)), shell=True)
    subprocess.call("tabix -f -S 2 -b 2 -e 3 {}.gz".format(Path(outfile)), shell=True)


def main():

    args = parse_args()
    copy_ratio_df = generate_copy_ratio_df(args)
    samples = copy_ratio_df.columns.tolist()[3:]

    # write output bed file
    write_outfile(
        copy_ratio_df, args.run, per_sample=False, hide_sample_names=False)

    if args.per_sample:
        # generating per sample bed files
        per_sample_dfs = generate_per_sample_copy_ratio_dfs(
            copy_ratio_df, args.keep_all_samples
        )

        if isinstance(per_sample_dfs, list):
            # writing per sample df with just sample trace + mean + std dev
            for sample in per_sample_dfs:
                sample_name = sample[0]
                sample_df = sample[1]
                write_outfile(
                    sample_df, sample_name, per_sample=True,
                    hide_sample_names=args.hide_sample_names
                )
        else:
            # writing per sample df with all samples, colouring the sample
            # trace red and leaving all other sample traces in grey
            for sample in samples:
                write_outfile(
                    per_sample_dfs, sample, per_sample=True,
                    hide_sample_names=args.hide_sample_names
                )


if __name__ == "__main__":
    main()
