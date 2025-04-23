import itertools
import os
from typing import Union

import edlib
import pysam
import ray


def reverse_complement(seq: str):
    tab = str.maketrans("ACTG", "TGAC")
    return seq.translate(tab)[::-1]


# def reverse_qual(qual: str):
#     return qual[::-1]


def extract_umi(
    query_seq: str,
    # query_qual: str,
    pattern: str,
    max_edit_dist: int,
):
    # umi_qual = None
    equalities = [
        ("M", "A"),
        ("M", "C"),
        ("R", "A"),
        ("R", "G"),
        ("W", "A"),
        ("W", "T"),
        ("S", "C"),
        ("S", "G"),
        ("Y", "C"),
        ("Y", "T"),
        ("K", "G"),
        ("K", "T"),
        ("V", "A"),
        ("V", "C"),
        ("V", "G"),
        ("H", "A"),
        ("H", "C"),
        ("H", "T"),
        ("D", "A"),
        ("D", "G"),
        ("D", "T"),
        ("B", "C"),
        ("B", "G"),
        ("B", "T"),
        ("N", "A"),
        ("N", "C"),
        ("N", "G"),
        ("N", "T"),
        ("m", "a"),
        ("m", "c"),
        ("r", "a"),
        ("r", "g"),
        ("w", "a"),
        ("w", "t"),
        ("s", "c"),
        ("s", "g"),
        ("y", "c"),
        ("y", "t"),
        ("k", "g"),
        ("k", "t"),
        ("v", "a"),
        ("v", "c"),
        ("v", "g"),
        ("h", "a"),
        ("h", "c"),
        ("h", "t"),
        ("d", "a"),
        ("d", "g"),
        ("d", "t"),
        ("b", "c"),
        ("b", "g"),
        ("b", "t"),
        ("n", "a"),
        ("n", "c"),
        ("n", "g"),
        ("n", "t"),
        ("a", "A"),
        ("c", "C"),
        ("t", "T"),
        ("g", "G"),
    ]

    result = edlib.align(
        pattern,
        query_seq,
        task="path",
        mode="HW",
        k=max_edit_dist,
        additionalEqualities=equalities,
    )
    if result["editDistance"] == -1:
        return None, None

    edit_dist = result["editDistance"]
    locs = result["locations"][0]
    umi_start_pos = locs[0]
    umi_end_pos = locs[1] + 1
    umi = query_seq[umi_start_pos:umi_end_pos]
    # umi_qual = query_qual[locs[0]:locs[1]+1]

    return edit_dist, umi


def extract_adapters(
    entry: pysam.libcfaidx.FastxRecord,
    adapter_length_5_end: int,
    adapter_length_3_end: int,
):
    read_5p_seq = None
    # read_5p_qual = None
    read_3p_seq = None
    # read_3p_qual = None

    read_5p_seq = entry.sequence[:adapter_length_5_end]
    read_3p_seq = entry.sequence[-adapter_length_3_end:]

    # read_5p_qual = entry.quality[:adapter_length_5_end]
    # read_3p_qual = entry.quality[-adapter_length_3_end:]

    return read_5p_seq, read_3p_seq


def get_read_name(entry: pysam.libcfaidx.FastxRecord):
    return entry.name.split(";")[0]


def get_read_strand(entry: pysam.libcfaidx.FastxRecord):
    strand = entry.name.split("strand=")
    if not len(strand) > 1:
        raise Exception("Read strand not annotated!")
    return strand[1]


def combine_umis_fasta(
    seq_5p: str,
    seq_3p: str,
    #    qual_5p: str,
    #    qual_3p: str,
    strand: str,
):
    if strand == "+":
        return seq_5p + seq_3p

    else:
        return reverse_complement(seq_3p) + reverse_complement(seq_5p)


def write_fasta(
    entry: pysam.libcfaidx.FastxRecord,
    strand: str,
    result_5p_fwd_umi_dist: int,
    result_3p_rev_umi_dist: int,
    result_5p_fwd_umi_seq: str,
    result_3p_rev_umi_seq: str,
    # result_5p_fwd_umi_qual: str,
    # result_3p_rev_umi_qual: str,
    out_fasta: Union[str, os.PathLike[str]],
):
    combined_umi_seq = combine_umis_fasta(
        seq_5p=result_5p_fwd_umi_seq,
        seq_3p=result_3p_rev_umi_seq,
        # qual_5p = result_5p_fwd_umi_qual,
        # qual_3p = result_3p_rev_umi_qual,
        strand=strand,
    )

    print(
        ">{};strand={};umi_fwd_dist={};umi_rev_dist={};umi_fwd_seq={};umi_rev_seq={};seq={}".format(
            entry.name,
            strand,
            result_5p_fwd_umi_dist,
            result_3p_rev_umi_dist,
            result_5p_fwd_umi_seq,
            result_3p_rev_umi_seq,
            entry.sequence,
        ),
        file=out_fasta,
    )

    print(combined_umi_seq, file=out_fasta)


@ray.remote(num_cpus=1)
def extract_umis(
    fastx_file: Union[str, os.PathLike[str]],
    umi_fasta_out_dir: Union[str, os.PathLike[str]],
    write_region: bool,
    # logs_dir: Union[str, os.PathLike[str]],
    adapter_length_5_end: int = 73,
    adapter_length_3_end: int = 68,
    max_pattern_dist: int = 3,
    umi_fwd: str = "TTTVVVVTTVVVVTTVVVVTTVVVVTTT",
    umi_rev: str = "AAABBBBAABBBBAABBBBAABBBBAAA",
):
    if write_region:
        region = os.path.basename(fastx_file).split(".")[0]
        umi_fasta_out_dir = os.path.join(umi_fasta_out_dir, region + "_detected_umis.fasta")
    else:
        umi_fasta_out_dir = os.path.join(umi_fasta_out_dir, "_detected_umis.fasta")

    # log_file = os.path.join(logs_dir, os.path.basename(region_fastq_file).split('.')[0] + '_extract_umis.err')

    # n_total = 0
    n_both_umi = 0
    # strand_stats = {"+": 0, "-": 0}
    with (
        pysam.FastxFile(fastx_file) as fastq_in,
        open(umi_fasta_out_dir, "w") as fasta_out,
    ):
        for entry in fastq_in:
            strand = get_read_strand(entry)
            entry.name = get_read_name(entry)
            # n_total += 1

            read_5p_seq, read_3p_seq = extract_adapters(
                entry=entry,
                adapter_length_5_end=adapter_length_5_end,
                adapter_length_3_end=adapter_length_3_end,
            )

            if read_5p_seq is None or read_3p_seq is None:
                continue

            # strand_stats[strand] += 1

            # Extract fwd UMI
            result_5p_fwd_umi_dist, result_5p_fwd_umi_seq = extract_umi(
                query_seq=read_5p_seq,
                # query_qual = read_5p_qual,
                pattern=umi_fwd,
                max_edit_dist=max_pattern_dist,
            )
            # Extract rev UMI
            result_3p_rev_umi_dist, result_3p_rev_umi_seq = extract_umi(
                query_seq=read_3p_seq,
                # query_qual = read_3p_qual,
                pattern=umi_rev,
                max_edit_dist=max_pattern_dist,
            )

            if not result_5p_fwd_umi_seq or not result_3p_rev_umi_seq:
                continue

            n_both_umi += 1

            write_fasta(
                entry=entry,
                strand=strand,
                result_5p_fwd_umi_dist=result_5p_fwd_umi_dist,
                result_3p_rev_umi_dist=result_3p_rev_umi_dist,
                result_5p_fwd_umi_seq=result_5p_fwd_umi_seq,
                result_3p_rev_umi_seq=result_3p_rev_umi_seq,
                # result_5p_fwd_umi_qual,
                # result_3p_rev_umi_qual,
                out_fasta=fasta_out,
            )

    if n_both_umi:
        return umi_fasta_out_dir
    else:
        return None


@ray.remote(num_cpus=1)
def count_single_umi_overlaps(umi_1_seq: str, umi_fasta_region_2_seqs: list, overlapping_umi_edit_threshold: int):
    # equalities = [("M", "A"), ("M", "C"), ("R", "A"), ("R", "G"), ("W", "A"), ("W", "T"), ("S", "C"), ("S", "G"), ("Y", "C"), ("Y", "T"), ("K", "G"), ("K", "T"), ("V", "A"), ("V", "C"),
    #           ("V", "G"), ("H", "A"), ("H", "C"), ("H", "T"), ("D", "A"), ("D", "G"), ("D", "T"), ("B", "C"), ("B", "G"), ("B", "T"), ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T"),
    #           ("m", "a"), ("m", "c"), ("r", "a"), ("r", "g"), ("w", "a"), ("w", "t"), ("s", "c"), ("s", "g"), ("y", "c"), ("y", "t"), ("k", "g"), ("k", "t"), ("v", "a"), ("v", "c"),
    #           ("v", "g"), ("h", "a"), ("h", "c"), ("h", "t"), ("d", "a"), ("d", "g"), ("d", "t"), ("b", "c"), ("b", "g"), ("b", "t"), ("n", "a"), ("n", "c"), ("n", "g"), ("n", "t"),
    #           ("a", "A"), ("c", "C"), ("t", "T"), ("g", "G")]

    umi_overlap_count = 0

    for umi_2_seq in umi_fasta_region_2_seqs:
        # if edlib.align(umi_1_seq,
        #                umi_2_seq,
        #                task="path",
        #                mode="HW",
        #                additionalEqualities=equalities)['editDistance'] <= overlapping_umi_edit_threshold:
        #     umi_overlap_count += 1
        if umi_1_seq == umi_2_seq:
            umi_overlap_count += 1

    return umi_overlap_count


@ray.remote(
    num_cpus=2
)  # TODO: test whether num_cpus can be changed to 1 without giving memory errors that kill workers
def count_overlapping_umis_between_2_regions(
    region_1_dir: Union[str, os.PathLike[str]],
    region_2_dir: Union[str, os.PathLike[str]],
    regions_w_overlapping_umis_tsv: Union[str, os.PathLike[str]],
    overlapping_umi_edit_threshold: int,
):
    region_1 = os.path.basename(region_1_dir)
    region_2 = os.path.basename(region_2_dir)
    logs_dir = os.path.dirname(regions_w_overlapping_umis_tsv)

    umi_fasta_region_2_seqs = []
    with pysam.FastxFile(os.path.join(region_2_dir, "umi_clusters_consensus.fasta")) as umi_fasta_region_2:
        for entry in umi_fasta_region_2:
            umi_fasta_region_2_seqs.append(entry.sequence)

    compare_edit_dist_umis_futures = []
    with pysam.FastxFile(os.path.join(region_1_dir, "umi_clusters_consensus.fasta")) as umi_fasta_region_1:
        for entry in umi_fasta_region_1:
            compare_edit_dist_umis_futures.append(
                count_single_umi_overlaps.remote(
                    umi_1_seq=entry.sequence,
                    umi_fasta_region_2_seqs=umi_fasta_region_2_seqs,
                    overlapping_umi_edit_threshold=overlapping_umi_edit_threshold,
                )
            )

    umi_overlap_count_list = ray.get(compare_edit_dist_umis_futures)

    if max(umi_overlap_count_list) > 1:
        with open(os.path.join(logs_dir, "region_region_umi_comparison.stderr"), "a") as ferr:
            print(
                "WARNING: there are UMIs from",
                region_1,
                "that match more than 1 UMI within",
                region_2,
                file=ferr,
            )

    umi_overlap_count = sum(umi_overlap_count_list)
    if umi_overlap_count:
        with open(os.path.join(regions_w_overlapping_umis_tsv), "a") as tsv_out:
            print(region_1, region_2, str(umi_overlap_count), sep="\t", file=tsv_out)

    if umi_overlap_count:
        return True
    else:
        return False


def count_overlapping_umis_between_all_regions(
    smolecule_filtered_fa_list: list,
    overlapping_umi_edit_threshold: int,
    logs_dir: Union[str, os.PathLike[str]],
):
    region_dirs = [os.path.dirname(smolecule_fa) for smolecule_fa in smolecule_filtered_fa_list]

    regions_w_overlapping_umis_tsv = os.path.join(logs_dir, "regions_w_overlapping_umis.tsv")
    with open(regions_w_overlapping_umis_tsv, "a") as tsv_out:
        print("region_1", "region_2", "umi_overlap_count", sep="\t", file=tsv_out)

    compare_umis_region_clusters_futures = []
    for region_1_dir, region_2_dir in itertools.combinations(region_dirs, 2):
        compare_umis_region_clusters_futures.append(
            count_overlapping_umis_between_2_regions.remote(
                region_1_dir=region_1_dir,
                region_2_dir=region_2_dir,
                regions_w_overlapping_umis_tsv=regions_w_overlapping_umis_tsv,
                overlapping_umi_edit_threshold=overlapping_umi_edit_threshold,
            )
        )

    region_comparisons_w_overlapping_umis_bool_list = ray.get(compare_umis_region_clusters_futures)

    return region_comparisons_w_overlapping_umis_bool_list
