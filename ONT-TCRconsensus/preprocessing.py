import subprocess
import os
from typing import Union
import ray


@ray.remote
def dorado_trim(
    fastq: Union[str, os.PathLike[str]],
    fastq_out_dir: Union[str, os.PathLike[str]],
    logs_dir: Union[str, os.PathLike[str]],
    dorado_threads: int,
    dorado_excutable: Union[str, os.PathLike[str]],
    nanopore_tcr_seq_primers_fasta: Union[
        str, os.PathLike[str]
    ],
    max_reads: int = None,
):
    fastq_trimmed = os.path.join(fastq_out_dir, os.path.basename(fastq).split(".")[0] + "_trimmed.fastq")

    log_file = os.path.join(logs_dir, os.path.basename(fastq).split(".")[0] + "_dorado_trim")

    with open(log_file + ".err", "w") as ferr, open(fastq_trimmed, "w") as fastq_out:
        if not max_reads:
            subprocess.run(
                [
                    dorado_excutable,
                    "trim",
                    "--primer-sequences",
                    nanopore_tcr_seq_primers_fasta,
                    "--emit-fastq",
                    "-v",
                    "--threads",
                    str(dorado_threads),
                    fastq,
                ],
                stderr=ferr,
                stdout=fastq_out,
            )
        else:
            subprocess.run(
                [
                    dorado_excutable,
                    "trim",
                    "--primer-sequences",
                    nanopore_tcr_seq_primers_fasta,
                    "--max-reads",
                    str(max_reads),
                    "--emit-fastq",
                    "-v",
                    "--threads",
                    str(dorado_threads),
                    fastq,
                ],
                stderr=ferr,
                stdout=fastq_out,
            )

    return fastq_trimmed


@ray.remote
def seqkit_fastq_filtering(
    fastq: Union[str, os.PathLike[str]],
    library_fastq_out_dir_dict: dict,
    library_logs_dir_dict: dict,
    minimal_average_qual: int,
    minimal_length: int,
    subsample_proportion: float = None,
):
    library = os.path.basename(os.path.dirname(fastq))
    fastq_out_dir = library_fastq_out_dir_dict[library]
    logs_dir = library_logs_dir_dict[library]
    filtered_fastq = os.path.join(fastq_out_dir, os.path.basename(fastq).split(".")[0] + "_filtered.fastq")
    filtered_fastq_sample = os.path.join(
        fastq_out_dir, os.path.basename(fastq).split(".")[0] + "_filtered_sample.fastq"
    )

    log_file = os.path.join(logs_dir, os.path.basename(fastq).split(".")[0] + "_seqkit")

    with open(log_file + "_stats_before_filtering.txt", "w") as txt_out:
        subprocess.run(["seqkit", "stat", "-a", fastq], stdout=txt_out)

    with open(log_file + ".err", "w") as ferr:
        subprocess.run(
            ["seqkit", "seq", "-m", str(minimal_length), "-Q", str(minimal_average_qual), fastq, "-o", filtered_fastq],
            stderr=ferr,
        )

        if subsample_proportion:
            subprocess.run(
                ["seqkit", "sample", "-p", str(subsample_proportion), filtered_fastq, "-o", filtered_fastq_sample],
                stderr=ferr,
            )
            subprocess.run(["rm", filtered_fastq])
            subprocess.run(["mv", filtered_fastq_sample, filtered_fastq])

    with open(log_file + "_stats_after_filtering.txt", "w") as txt_out:
        subprocess.run(["seqkit", "stat", "-a", filtered_fastq], stdout=txt_out)

    return filtered_fastq


@ray.remote(num_cpus=1)  # vsearch filtering does not allow more than 1 cpu...
def vsearch_fastq_filtering(
    fastq: Union[str, os.PathLike[str]],
    library_fastq_out_dir_dict: dict,
    library_logs_dir_dict: dict,
    max_ee_rate_base: float,
    minimal_length: int,
):
    library = os.path.basename(os.path.dirname(fastq))
    fastq_out_dir = library_fastq_out_dir_dict[library]
    logs_dir = library_logs_dir_dict[library]
    filtered_fastq = os.path.join(fastq_out_dir, os.path.basename(fastq).split(".")[0] + "_filtered.fastq")

    log_file = os.path.join(logs_dir, os.path.basename(fastq).split(".")[0] + "_vsearch_fastq_filtering.log")
    stats_file = os.path.join(logs_dir, os.path.basename(fastq).split(".")[0])

    # subprocess.run(['vsearch', '--fastq_stats', fastq,
    #                 '--fastq_ascii', '33',
    #                 '--fastq_qmax', '93',
    #                 '--quiet',
    #                 '--log', stats_file + '_stats_before_vsearch_filtering.log'])

    with open(stats_file + "_stats_before_vsearch_filtering.txt", "w") as txt_out:
        subprocess.run(["seqkit", "stat", "-a", fastq], stdout=txt_out)

    subprocess.run(
        [
            "vsearch",
            "--fastq_filter",
            fastq,
            "--fastq_ascii",
            "33",
            "--fastq_qmax",
            "93",
            "--fastq_maxee_rate",
            str(max_ee_rate_base),
            "--fastq_minlen",
            str(minimal_length),
            "--quiet",
            "--log",
            log_file,
            "-fastqout",
            filtered_fastq,
        ]
    )

    # subprocess.run(['vsearch', '--fastq_stats', filtered_fastq,
    #                 '--fastq_ascii', '33',
    #                 '--fastq_qmax', '93',
    #                 '--quiet',
    #                 '--log', stats_file + '_stats_after_vsearch_filtering.log'])

    with open(stats_file + "_stats_after_vsearch_filtering.txt", "w") as txt_out:
        subprocess.run(["seqkit", "stat", "-a", filtered_fastq], stdout=txt_out)

    return filtered_fastq
