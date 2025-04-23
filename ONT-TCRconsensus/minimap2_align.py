import collections
import os
import subprocess
from typing import Union

import numpy as np
import pandas as pd
import pysam

from nanopore_tcr_consensus.region_split import generate_region_length_dict_from_ref_fa


def calculate_blast_id(entry: pysam.libcalignedsegment.AlignedSegment):
    number_of_alignment_columns = sum(
        [cigar_tuple[1] for cigar_tuple in entry.cigartuples if cigar_tuple[0] in [0, 1, 2]]
    )
    number_of_matches = number_of_alignment_columns - entry.get_tag("NM")
    return number_of_matches / number_of_alignment_columns


def count_cs_tags(bam_file: Union[str, os.PathLike[str]]):
    cs_tag_counter = collections.Counter()
    cs_tag_region_counter = collections.defaultdict(lambda: collections.Counter())
    cs_tag_blast_id_counter = collections.defaultdict(lambda: collections.Counter())

    with pysam.AlignmentFile(bam_file, "rb") as bam_in:
        for entry in bam_in:
            if entry.is_unmapped or entry.is_secondary or entry.is_supplementary:
                continue

            if entry.has_tag("cs"):
                cs_tag = entry.get_tag("cs")
                cs_tag_counter[cs_tag] += 1
                cs_tag_region_counter[cs_tag][entry.reference_name] += 1
                cs_tag_blast_id_counter[cs_tag][calculate_blast_id(entry=entry)] += 1

    return cs_tag_counter, cs_tag_region_counter, cs_tag_blast_id_counter


def minimap2_ref_self_homology_map(reference: Union[str, os.PathLike[str]], minimap2_threads: int):
    log_file = os.path.join(
        os.path.dirname(reference), os.path.basename(reference).split(".")[0] + "_minimap2_ref_self_homology_map"
    )
    ref_homology_paf = os.path.join(
        os.path.dirname(reference), os.path.basename(reference).split(".")[0] + "_ref_homology_out.paf"
    )

    with open(ref_homology_paf, "w") as paf_out, open(log_file + ".err", "w") as ferr:
        minimap2_process = subprocess.Popen(
            [
                "minimap2",
                "-DP",
                "-c",
                "-k19",
                "-w 19",
                "-U50,500",
                "-g10k",
                "-m200",
                "--end-bonus=10",
                "-t",
                str(minimap2_threads),
                "--cap-kalloc",
                "500m",  # New --cap-kalloc default setting from 2.28 in 2.24
                reference,
                reference,
            ],
            stdout=paf_out,
            stderr=ferr,
        )

        minimap2_process.wait()

    return ref_homology_paf


def minimap2_ont_align(
    fastx: Union[str, os.PathLike[str]],
    reference: Union[str, os.PathLike[str]],
    bam_out_dir: Union[str, os.PathLike[str]],
    logs_dir: Union[str, os.PathLike[str]],
    minimap2_threads: int,
    minimap2_params: str = None,
):
    bam_file = os.path.join(bam_out_dir, os.path.basename(fastx).split(".")[0] + ".bam")
    log_file = os.path.join(logs_dir, os.path.basename(fastx).split(".")[0] + "_minimap2")
    flagstat_file = os.path.join(logs_dir, os.path.basename(bam_file).split(".")[0] + "_flagstat.out")

    with open(log_file + ".err", "w") as ferr:
        if minimap2_params:
            minimap2_process = subprocess.Popen(
                [
                    "minimap2",
                    "-ax",
                    "map-ont",
                    "-k19",  # Manually implement lr:hq preset from 2.27 in 2.24
                    "-w 19",
                    "-U50,500",
                    "-g10k",
                    "--end-bonus=10",
                    "-t",
                    str(minimap2_threads),
                    "--cap-kalloc",
                    "500m",  # New --cap-kalloc default setting from 2.28 in 2.24
                    minimap2_params,
                    reference,
                    fastx,
                ],
                stdout=subprocess.PIPE,
                stderr=ferr,
            )
        else:
            minimap2_process = subprocess.Popen(
                [
                    "minimap2",
                    "-ax",
                    "map-ont",
                    "-k19",  # Manually implement lr:hq preset from 2.27 in 2.24
                    "-w 19",
                    "-U50,500",
                    "-g10k",
                    "--end-bonus=10",
                    "-t",
                    str(minimap2_threads),
                    "--cap-kalloc",
                    "500m",  # New --cap-kalloc default setting from 2.28 in 2.24
                    "--cs",
                    reference,
                    fastx,
                ],
                stdout=subprocess.PIPE,
                stderr=ferr,
            )
        subprocess.run(
            ["samtools", "sort", "-@", str(minimap2_threads), "-o", bam_file],
            stdin=minimap2_process.stdout,
            stderr=ferr,
        )
        subprocess.run(["samtools", "index", "-@", str(minimap2_threads), bam_file])

        cs_tag_counter, cs_tag_region_counter, cs_tag_blast_id_counter = count_cs_tags(bam_file=bam_file)
        ferr.write("\nTop 40 most common cs tags:\n")
        top_40_cs_tags = cs_tag_counter.most_common(40)
        for cs_tag_tuple in top_40_cs_tags:
            ferr.write(str(cs_tag_tuple) + "\n")
        ferr.write("\nTop 4 most common regions counted for each of the top 40 most common cs tags:\n")
        for cs_tag_tuple in top_40_cs_tags:
            ferr.write(str(cs_tag_tuple[0]) + " " + str(cs_tag_region_counter[cs_tag_tuple[0]].most_common(4)) + "\n")
        ferr.write("\nTop 4 most common blast identities counted for each of the top 40 most common cs tags:\n")
        for cs_tag_tuple in top_40_cs_tags:
            ferr.write(str(cs_tag_tuple[0]) + " " + str(cs_tag_blast_id_counter[cs_tag_tuple[0]].most_common(4)) + "\n")

    with open(flagstat_file, "w") as flag_out:
        subprocess.run(["samtools", "flagstat", "-@", str(minimap2_threads), bam_file], stdout=flag_out)

    return bam_file


def filter_consensus_alignments(
    consensus_bam_file: Union[str, os.PathLike[str]],
    reference: Union[str, os.PathLike[str]],
    logs_dir: Union[str, os.PathLike[str]],
    blast_id_threshold: float,
    minimal_region_overlap: float = 0.9975,
    max_softclip_5_end: int = 73,
    max_softclip_3_end: int = 68,
):
    filtered_consensus_bam_file = os.path.splitext(consensus_bam_file)[0] + "_filtered.bam"
    log_file = os.path.join(logs_dir, os.path.basename(consensus_bam_file).split(".")[0] + "_bam_filter.log")
    nt_too_short_csv = os.path.join(logs_dir, os.path.basename(consensus_bam_file).split(".")[0] + "_nt_too_short.csv")
    region_nt_too_short_csv = os.path.join(
        logs_dir, os.path.basename(consensus_bam_file).split(".")[0] + "_region_nt_too_short.csv"
    )
    nt_too_long_csv = os.path.join(logs_dir, os.path.basename(consensus_bam_file).split(".")[0] + "_nt_too_long.csv")
    region_nt_too_long_csv = os.path.join(
        logs_dir, os.path.basename(consensus_bam_file).split(".")[0] + "_region_nt_too_long.csv"
    )
    blast_id_csv = os.path.join(logs_dir, os.path.basename(consensus_bam_file).split(".")[0] + "_blast_id.csv")
    region_blast_id_csv = os.path.join(
        logs_dir, os.path.basename(consensus_bam_file).split(".")[0] + "_region_blast_id.csv"
    )
    num_subreads_blast_id_csv = os.path.join(
        logs_dir, os.path.basename(consensus_bam_file).split(".")[0] + "_number_of_subreads_blast_id.csv"
    )

    allowed_nt_too_short_list = []
    allowed_nt_too_long_list = []
    allowed_nt_diff_similarity_list = []
    region_length_dict = generate_region_length_dict_from_ref_fa(reference=reference)

    for region_len in region_length_dict.values():
        allowed_nt_too_short_list.append(region_len - (region_len * minimal_region_overlap))
        allowed_nt_too_long_list.append((region_len * (2 - minimal_region_overlap)) - region_len)
        allowed_nt_diff_similarity_list.append(region_len - (region_len * blast_id_threshold))

    n_short = 0
    n_long = 0
    n_written = 0
    n_primary_mapped = 0
    n_primary_mapped_correct_len = 0
    region_blast_id_list_dict = collections.defaultdict(list)
    region_nt_too_short_list_dict = collections.defaultdict(list)
    region_nt_too_long_list_dict = collections.defaultdict(list)
    num_subreads_blast_id_list_dict = collections.defaultdict(list)

    with (
        pysam.AlignmentFile(consensus_bam_file, "rb") as bam_in,
        pysam.AlignmentFile(filtered_consensus_bam_file, "wb", template=bam_in) as bam_out,
    ):
        for entry in bam_in:
            if entry.is_unmapped:
                continue
            if entry.is_secondary or entry.is_supplementary:
                continue
            n_primary_mapped += 1

            if entry.reference_length < (region_length_dict[entry.reference_name] * minimal_region_overlap):
                n_short += 1
                # nt_too_short_list.append((region_length_dict[entry.reference_name] * minimal_region_overlap) - entry.reference_length)
                region_nt_too_short_list_dict[entry.reference_name].append(
                    (region_length_dict[entry.reference_name] * minimal_region_overlap) - entry.reference_length
                )
                continue
            if entry.query_length > (
                region_length_dict[entry.reference_name] * (2 - minimal_region_overlap)
                + (max_softclip_5_end + max_softclip_3_end)
            ):
                n_long += 1
                # nt_too_long_list.append(entry.query_length - (region_length_dict[entry.reference_name] * (2 - minimal_region_overlap) + (max_softclip_5_end + max_softclip_3_end)))
                region_nt_too_long_list_dict[entry.reference_name].append(
                    entry.query_length
                    - (
                        region_length_dict[entry.reference_name] * (2 - minimal_region_overlap)
                        + (max_softclip_5_end + max_softclip_3_end)
                    )
                )
                continue
            n_primary_mapped_correct_len += 1

            blast_id = calculate_blast_id(entry=entry)
            region_blast_id_list_dict[entry.reference_name].append(blast_id)
            num_subreads = entry.query_name.split("_")[-1]
            num_subreads_blast_id_list_dict[num_subreads].append(blast_id)
            if blast_id > blast_id_threshold:
                n_written += 1
                bam_out.write(entry)

    subprocess.run(["samtools", "index", filtered_consensus_bam_file])

    nt_too_short_list = [nt for nt_list in region_nt_too_short_list_dict.values() for nt in nt_list]
    nt_too_long_list = [nt for nt_list in region_nt_too_long_list_dict.values() for nt in nt_list]
    blast_id_list = [blast_id for blast_id_list in region_blast_id_list_dict.values() for blast_id in blast_id_list]

    pd.DataFrame(
        [
            (num_subreads, blast_id)
            for num_subreads, blast_id_list in num_subreads_blast_id_list_dict.items()
            for blast_id in blast_id_list
        ],
        columns=["number_of_subreads", "blast_id"],
    ).to_csv(num_subreads_blast_id_csv, index=False)
    pd.DataFrame(
        [(region, nt) for region, nt_list in region_nt_too_short_list_dict.items() for nt in nt_list],
        columns=["region", "number_of_nt"],
    ).to_csv(region_nt_too_short_csv, index=False)
    pd.DataFrame(
        [(region, nt) for region, nt_list in region_nt_too_long_list_dict.items() for nt in nt_list],
        columns=["region", "number_of_nt"],
    ).to_csv(region_nt_too_long_csv, index=False)
    pd.DataFrame(
        [
            (region, blast_id)
            for region, blast_id_list in region_blast_id_list_dict.items()
            for blast_id in blast_id_list
        ],
        columns=["region", "blast_id"],
    ).to_csv(region_blast_id_csv, index=False)
    pd.DataFrame(nt_too_short_list, columns=["number_of_nt"]).to_csv(nt_too_short_csv, index=False)
    pd.DataFrame(nt_too_long_list, columns=["number_of_nt"]).to_csv(nt_too_long_csv, index=False)
    pd.DataFrame(blast_id_list, columns=["blast_id"]).to_csv(blast_id_csv, index=False)

    with open(log_file, "w") as log_out:
        log_out.write("Consensus alignment filtering performed with the following parameters:\n")
        log_out.write("- minimal region overlap: " + str(minimal_region_overlap) + "\n")
        log_out.write("- minimal blast identity with reference: " + str(blast_id_threshold) + "\n")
        log_out.write("From these parameters follows:\n")
        log_out.write("- Minimal Phred Q = " + str(round(-10 * np.log10(1 - blast_id_threshold), 2)) + "\n")
        log_out.write("- Median region nucleotide length: " + str(np.median(list(region_length_dict.values()))) + "\n")
        log_out.write(
            "- Median allowed too few nucleotides/region: " + str(round(np.median(allowed_nt_too_short_list), 2)) + "\n"
        )
        log_out.write(
            "- Median allowed too many nucleotides/region: " + str(round(np.median(allowed_nt_too_long_list), 2)) + "\n"
        )
        log_out.write(
            "- Median allowed variants (i.e., indels and mismatches)/region: "
            + str(round(np.median(allowed_nt_diff_similarity_list), 2))
            + "\n\n"
        )
        log_out.write("Total # primary alignments in bam file before filtering: " + str(n_primary_mapped) + "\n")
        log_out.write(
            "% of primary alignments that have shorter overlap than minimal region overlap: "
            + str(round(100 * n_short / n_primary_mapped, 2))
            + "\n"
        )
        log_out.write(
            "% of primary alignments that have too long reads: " + str(round(100 * n_long / n_primary_mapped, 2)) + "\n"
        )
        log_out.write(
            "Median # of too few nucleotides of alignments that have shorter than minimal region overlap: "
            + str(round(np.median(nt_too_short_list), 2))
            + "\n"
        )
        log_out.write(
            "Median # of too many nucleotides of alignments that have longer than minimal region overlap: "
            + str(round(np.median(nt_too_long_list), 2))
            + "\n"
        )
        log_out.write("Median blast identity of alignments: " + str(round(np.median(blast_id_list), 6)) + "\n")
        log_out.write(
            "15th percentile blast identity of alignments: " + str(round(np.quantile(blast_id_list, 0.15), 6)) + "\n"
        )
        log_out.write(
            "25th percentile blast identity of alignments: " + str(round(np.quantile(blast_id_list, 0.25), 6)) + "\n"
        )
        log_out.write(
            "75th percentile blast identity of alignments: " + str(round(np.quantile(blast_id_list, 0.75), 6)) + "\n"
        )
        log_out.write("Average blast identity of alignments: " + str(round(np.mean(blast_id_list), 6)) + "\n")
        log_out.write("Total # primary alignments written to filtered bam file: " + str(n_written) + "\n")
        log_out.write(
            "% of primary alignments written to filtered bam file: "
            + str(round(100 * n_written / n_primary_mapped, 2))
            + "\n"
        )
        log_out.write(
            "% primary alignments of primary alignments with minimal overlap written to filtered bam file: "
            + str(round(100 * n_written / n_primary_mapped_correct_len, 2))
            + "\n"
        )

        cs_tag_counter, cs_tag_region_counter, cs_tag_blast_id_counter = count_cs_tags(
            bam_file=filtered_consensus_bam_file
        )
        log_out.write("\nTop 40 most common cs tags:\n")
        top_40_cs_tags = cs_tag_counter.most_common(40)
        for cs_tag_tuple in top_40_cs_tags:
            log_out.write(str(cs_tag_tuple) + "\n")
        log_out.write("\nTop 4 most common regions counted for each of the top 40 most common cs tags:\n")
        for cs_tag_tuple in top_40_cs_tags:
            log_out.write(
                str(cs_tag_tuple[0]) + " " + str(cs_tag_region_counter[cs_tag_tuple[0]].most_common(4)) + "\n"
            )
        log_out.write("\nTop 4 most common blast identities counted for each of the top 40 most common cs tags:\n")
        for cs_tag_tuple in top_40_cs_tags:
            log_out.write(
                str(cs_tag_tuple[0]) + " " + str(cs_tag_blast_id_counter[cs_tag_tuple[0]].most_common(4)) + "\n"
            )

    return filtered_consensus_bam_file


def estimate_precision_at_num_subreads(
    consensus_bam_file: Union[str, os.PathLike[str]],
    reference: Union[str, os.PathLike[str]],
    logs_dir: Union[str, os.PathLike[str]],
    minimal_region_overlap: float = 1.0,
    blast_id_threshold: float = 1.0,
    max_softclip_5_end: int = 73,
    max_softclip_3_end: int = 68,
    num_subreads_threshold: int = 6,
):
    precision_at_num_subreads_csv = os.path.join(
        logs_dir, os.path.basename(consensus_bam_file).split(".")[0] + "_precision_at_different_num_subreads.csv"
    )

    region_length_dict = generate_region_length_dict_from_ref_fa(reference=reference)
    num_subreads_primary_mapped_counter = collections.defaultdict(int)
    num_subreads_written_counter = collections.defaultdict(int)
    n_primary_mapped = 0
    n_primary_mapped_above_num_subreads_threshold = 0
    n_written_above_num_subreads_threshold = 0

    with pysam.AlignmentFile(consensus_bam_file, "rb") as bam_in:
        for entry in bam_in:
            if entry.is_unmapped:
                continue
            if entry.is_secondary or entry.is_supplementary:
                continue

            n_primary_mapped += 1
            num_subreads = int(entry.query_name.split("_")[-1])
            num_subreads_primary_mapped_counter[num_subreads] += 1
            if num_subreads >= num_subreads_threshold:
                n_primary_mapped_above_num_subreads_threshold += 1

            if entry.reference_length < (region_length_dict[entry.reference_name] * minimal_region_overlap):
                continue
            if entry.query_length > (
                region_length_dict[entry.reference_name] * (2 - minimal_region_overlap)
                + (max_softclip_5_end + max_softclip_3_end)
            ):
                continue

            blast_id = calculate_blast_id(entry=entry)
            if blast_id >= blast_id_threshold:
                num_subreads_written_counter[num_subreads] += 1
                if num_subreads >= num_subreads_threshold:
                    n_written_above_num_subreads_threshold += 1

    precision_at_num_subreads_series = (
        (pd.Series(num_subreads_written_counter) / pd.Series(num_subreads_primary_mapped_counter))
        .fillna(0)
        .sort_index()
    )
    precision_at_num_subreads_series.name = "precision"
    precision_at_num_subreads_series.to_csv(precision_at_num_subreads_csv, header=True, index_label="num_subreads")
    print(n_written_above_num_subreads_threshold, n_primary_mapped_above_num_subreads_threshold)
    print(
        "median # of subreads:",
        np.median(
            [num_subreads for num_subreads, count in num_subreads_primary_mapped_counter.items() for _ in range(count)]
        ),
    )
    print(
        "mean # of subreads:",
        np.mean(
            [num_subreads for num_subreads, count in num_subreads_primary_mapped_counter.items() for _ in range(count)]
        ),
    )
    print("n primary mapped:", n_primary_mapped)
    print(
        "Overall precision above " + str(num_subreads_threshold) + " is:",
        n_written_above_num_subreads_threshold / n_primary_mapped_above_num_subreads_threshold,
    )
    # (pd.Series(num_subreads_written_counter) / pd.Series(num_subreads_primary_mapped_counter)).fillna(0).sort_index().to_csv(precision_at_num_subreads_csv, header='precision', index_label='num_subreads')
