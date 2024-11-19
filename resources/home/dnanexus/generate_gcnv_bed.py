"""
Script to generate bed file from a run of cnv calls and intervals list.

Requires bgzip and tabix to be installed and on path.
"""

import argparse
from concurrent.futures import as_completed, ProcessPoolExecutor
from os import cpu_count
from pathlib import Path
import subprocess
from timeit import default_timer as timer
from typing import Tuple

import pandas as pd


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--copy_ratios",
        nargs="+",
        help="all copy ratio files to generate bed file from",
    )
    parser.add_argument(
        "--run", help="name of run to prefix output per run bed file with"
    )
    parser.add_argument(
        "-s",
        "--per_sample",
        action="store_true",
        help=(
            "generate per sample bed file, with just the sample highlighted "
            "and run level mean and std deviations added"
        ),
    )
    parser.add_argument(
        "--keep_all_samples",
        action="store_true",
        default=True,
        help=(
            "controls if to keep all sample traces in per sample visualisation"
            " bed files, if true all other samples will be anonymised grey"
            " traces"
        ),
    )

    args = parser.parse_args()

    # ensure interval list isn't accidentally passed as copy ratio file
    args.copy_ratios = [
        x for x in args.copy_ratios if "interval_list" not in x
    ]

    return args


def read_single_copy_ratio_file(copy_ratio_file) -> Tuple[str, pd.DataFrame]:
    """
    Read single copy ratio file in to DataFrame, parsing the sample name
    from the `@RG` header line

    Parameters
    ----------
    copy_ratio_file : str
        copy ratio filename to read in

    Returns
    -------
    str
        name of sample
    pd.DataFrame
        intervals and copy ratio values read in from file
    """
    with open(copy_ratio_file, encoding="utf-8", mode="r") as fh:
        # read through file until last line of header with sample name
        while True:
            line = fh.readline()
            if line.startswith("@RG"):
                sample_name = line.split("SM:")[1].strip()
                break

    file_df = pd.read_csv(
        copy_ratio_file,
        sep="\t",
        comment="@",
        header=0,
        names=["chr", "start", "end", sample_name],
    )

    return sample_name, file_df


def read_all_copy_ratio_files(copy_ratio_files) -> pd.DataFrame:
    """
    Read in all copy ratio files to single dataframe.

    DataFrame will have the following structure:

    | chr | start    | end      | sample_1 | sample_2 | ...
    |-----|----------|----------|----------|----------|----
    | 1   | 10566155 | 10566275 | 1.9836   | 2.0889   | ...
    | 1   | 17345175 | 17345329 | 2.1342   | 2.6353   | ...
    | 1   | 17345329 | 17345483 | 1.9385   | 1.8833   | ...
    | ... | ...      | ...      | ...      | ...      | ...

    Parameters
    ----------
    copy_ratio_files : list
        list of copy ratio files for the run to read in

    Returns
    -------
    pd.DataFrame
        DataFrame of all intervals for all samples
    """
    print("\nReading in copy ratio tsv files")
    start = timer()

    # get intervals from the first file to use for comparing all against
    _, copy_ratio_df = read_single_copy_ratio_file(
        copy_ratio_file=copy_ratio_files[0]
    )
    copy_ratio_df = copy_ratio_df[["chr", "start", "end"]]

    for copy_ratio_file in copy_ratio_files:
        sample, sample_df = read_single_copy_ratio_file(copy_ratio_file)

        assert copy_ratio_df[["chr", "start", "end"]].equals(
            sample_df[["chr", "start", "end"]]
        ), f"Copy ratio file for {sample} has different intervals "

        copy_ratio_df[sample] = sample_df[sample]

    print(
        f"Completed reading {len(copy_ratio_files)} copy ratio "
        f"files in {round(timer() - start, 2)}s"
    )

    return copy_ratio_df


def calculate_mean_and_std_dev(copy_ratio_df) -> pd.DataFrame:
    """
    Calculate the mean and 2 standard deviations for each interval
    across all samples

    Parameters
    ----------
    copy_ratio_df : pd.DataFrame
        DataFrame of all intervals for all samples

    Returns
    -------
    pd.DataFrame
        DataFrame of all intervals for all samples with mean and std devs
    """
    print("\nCalculating mean values for copy ratios")

    samples = copy_ratio_df.columns.tolist()[3:]
    start = timer()

    # calculate mean and std dev across all samples (rows)
    mean = copy_ratio_df[samples].mean(axis=1)
    std = copy_ratio_df[samples].std(ddof=1, axis=1)

    plus_std = mean + std
    minus_std = mean - std
    plus_std2 = mean + (std * 2)
    minus_std2 = mean - (std * 2)

    mean_std_df = pd.DataFrame.from_dict(
        {
            "mean": mean,
            "mean_plus_std": plus_std,
            "mean_plus_std2": plus_std2,
            "mean_minus_std": minus_std,
            "mean_minus_std2": minus_std2,
        }
    )

    copy_ratio_df = pd.concat([copy_ratio_df, mean_std_df], axis="columns")

    print(f"Completed generating means in {round(timer() - start, 2)}s")
    return copy_ratio_df


def write_run_level_bed_file(copy_ratio_df, prefix) -> None:
    """
    Writes the run level bed file containing all samples with their
    labels, and clickToHighlight enabled for all

    Parameters
    ----------
    copy_ratio_df : pd.DataFrame
        DataFrame of copy ratios for all samples
    prefix : str
        name of run to prefix output file name with
    """
    outfile = f"{prefix}_copy_ratios.gcnv.bed"

    with open(outfile, "w") as fh:
        # this line is needed at the beginning to tell igv.js that it is
        # a gcnv bed as it can't automatically set this track type
        fh.write("track type=gcnv height=500 clickToHighlight=any \n")

    copy_ratio_df.to_csv(outfile, sep="\t", header=True, index=False, mode="a")
    compress_and_index_bed_file(outfile)


def write_sample_bed_file(
    copy_ratio_df, sample_name, keep_all_samples
) -> None:
    """
    Writes individual sample bed file.

    If `--keep_all_samples` has been specified, then all other sample
    values will be kept as grey traces with labels removed (i.e so that
    they are anonymised and only the given sample is labelled)

    Parameters
    ----------
    copy_ratio_df : pd.DataFrame
        DataFrame of all intervals for all samples
    sample_name : str
        name of sample
    keep_all_samples : bool
        controls if to keep all sample values in the output file
    """
    print(f"Creating output bed file for {sample_name}")

    outfile = f"{sample_name}_copy_ratios.gcnv.bed"

    # colour mapping for tracks, keys have to match column names in dataframe
    colours = {
        sample_name: "red",
        "mean": "blue",
        "mean_plus_std": "#0B2559",
        "mean_minus_std": "#0B2559",
        "mean_plus_std2": "#36BFB1",
        "mean_minus_std2": "#36BFB1",
    }

    minimum_columns = [
        "chr",
        "start",
        "end",
        "mean",
        "mean_plus_std",
        "mean_plus_std2",
        "mean_minus_std",
        "mean_minus_std2",
        sample_name,
    ]

    if not keep_all_samples:
        # not keeping all other samples as grey traces => remove the columns
        copy_ratio_df = copy_ratio_df[minimum_columns]
        header = "\t".join(copy_ratio_df.columns.tolist())
    else:
        # keeping other samples => remove their column labels from header
        # by setting them to a blank space
        header = "\t".join(
            [
                "⠀" if x not in minimum_columns else x
                for x in copy_ratio_df.columns.tolist()
            ]
        )

    with open(outfile, encoding="utf-8", mode="w") as fh:
        # this line is needed at the beginning to tell igv.js that it is
        # a gcnv bed as it can't automatically set this track type
        # onlyHandleClicksForHighlightedSamples set so that the other
        # sample traces in grey with no labels can't be clicked on
        highlight = " ".join(
            [f"highlight={x};{y}" for x, y in colours.items()]
        )
        fh.write(
            "track type=gcnv height=500"
            f" onlyHandleClicksForHighlightedSamples=true {highlight} \n"
        )

        fh.write(f"{header}\n")

    copy_ratio_df.to_csv(
        outfile, sep="\t", header=False, index=False, mode="a"
    )

    compress_and_index_bed_file(outfile)


def compress_and_index_bed_file(bed_file) -> None:
    """
    Call bgzip and tabix on output bed file to compress and index

    Parameters
    ----------
    bed_file : str
        bed file to compress and index
    """
    subprocess.run(f"bgzip {Path(bed_file)}", shell=True, check=True)
    subprocess.run(
        f"tabix -f -S 2 -b 2 -e 3 {Path(bed_file)}.gz",
        shell=True,
        check=True,
    )


def main() -> None:
    args = parse_args()

    copy_ratio_df = read_all_copy_ratio_files(
        copy_ratio_files=args.copy_ratios
    )
    write_run_level_bed_file(copy_ratio_df=copy_ratio_df, prefix=args.run)

    if args.per_sample:
        copy_ratio_df = calculate_mean_and_std_dev(copy_ratio_df)

        samples = [
            Path(x).name.replace("_denoised_copy_ratios.tsv", "")
            for x in args.copy_ratios
        ]

        with ProcessPoolExecutor(max_workers=cpu_count()) as executor:
            concurrent_jobs = {
                executor.submit(
                    write_sample_bed_file,
                    sample_name=sample,
                    copy_ratio_df=copy_ratio_df,
                    keep_all_samples=args.keep_all_samples,
                ): sample
                for sample in samples
            }

            for future in as_completed(concurrent_jobs):
                # access returned output as each is returned in any order
                try:
                    future.result()
                except Exception as error:
                    print(
                        "Error in writing output file for"
                        f" {concurrent_jobs[future]}"
                    )
                    raise error


if __name__ == "__main__":
    main()
