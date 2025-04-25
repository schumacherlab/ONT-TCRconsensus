import os
from typing import Union
import pysam
import collections
import pandas as pd
import numpy as np
import json


def parse_bed(bed_regions):
    with open(bed_regions) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                raise Exception("Bed file is broken because line contains less than 4 columns!")

            region = {
                "chr": cols[0],
                "start": int(cols[1]),
                "end": int(cols[2]),
                "name": cols[3],
            }
            yield region


def generate_regions_set_from_ref_fa(
    reference: Union[str, os.PathLike[str]],
    remove_negative_control_regions: bool = False,
    negative_control_regions_suffix_str_tuple: tuple = None,
):
    if remove_negative_control_regions:
        regions_set = set()
        with pysam.FastxFile(reference) as fasta_in:
            for entry in fasta_in:
                if entry.name.endswith(negative_control_regions_suffix_str_tuple):
                    continue
                regions_set.add(entry.name)

        return regions_set

    regions_set = set()
    with pysam.FastxFile(reference) as fasta_in:
        for entry in fasta_in:
            regions_set.add(entry.name)

    return regions_set


def generate_region_length_dict_from_ref_fa(reference: Union[str, os.PathLike[str]]):
    region_length_dict = {}
    with pysam.FastxFile(reference) as fasta_in:
        for entry in fasta_in:
            region_length_dict[entry.name] = len(entry.sequence)

    return region_length_dict


def greedy_minimap2_most_similar_region_clustering(
    minimap2_most_similar_region_tuple_list: list, similarity_threshold: float
):
    sorted_data = sorted(minimap2_most_similar_region_tuple_list, key=lambda x: x[2], reverse=True)

    similar_region_clusters = []
    seen = set()

    for string1, string2, similarity in sorted_data:
        if string1 not in seen and string2 not in seen:
            if similarity >= similarity_threshold:
                similar_region_clusters.append({string1, string2})
                seen.update([string1, string2])
        elif string1 in seen or string2 in seen:
            for cluster in similar_region_clusters:
                if string1 in cluster or string2 in cluster:
                    if similarity >= similarity_threshold:
                        cluster.update([string1, string2])
                        seen.update([string1, string2])  # Mark both as seen
                    break

    return similar_region_clusters


def generate_region_split_dict_from_minimap2_homology_paf(
    ref_homology_paf: Union[str, os.PathLike[str]],
    reference: Union[str, os.PathLike[str]],
    id_cluster_threshold: float = 0.98,
):
    # minimap2 -DP -c -k19 -w19 -U50,500 -g10k -m200 -t 99 hc_cd8_cd4_combined_nanopore.fa hc_cd8_cd4_combined_nanopore.fa > ref_homology_out.paf
    log_file = os.path.join(
        os.path.dirname(ref_homology_paf),
        os.path.basename(ref_homology_paf).split(".")[0] + "_generate_region_split_dict.log",
    )
    region_cluster_dict_json = os.path.join(
        os.path.dirname(ref_homology_paf), os.path.basename(ref_homology_paf).split(".")[0] + "_region_split_dict.json"
    )

    minimap2_most_similar_region_dict_json = os.path.join(
        os.path.dirname(ref_homology_paf),
        os.path.basename(ref_homology_paf).split(".")[0] + "_most_similar_region_dict.json",
    )
    # region_cluster_list_txt = os.path.join(os.path.dirname(ref_homology_paf),
    #                                        os.path.basename(ref_homology_paf).split('.')[0] + '_region_split_dict.txt')

    with open(log_file, "w") as log_out:
        ref_homology_out_df = pd.read_csv(ref_homology_paf, sep="\t", header=None)
        ref_homology_out_df.loc[:, "min_tcr_length"] = ref_homology_out_df.loc[:, [1, 6]].min(axis=1)
        print(
            "Homology .paf shape before minimal alignment overlap filtering:",
            ref_homology_out_df.shape[0],
            file=log_out,
        )
        ref_homology_out_df = ref_homology_out_df.loc[
            ref_homology_out_df[10] > ref_homology_out_df.loc[:, "min_tcr_length"] * 0.99, :
        ]
        ref_homology_out_df.reset_index(inplace=True, drop=True)
        print(
            "Homology .paf shape after minimal alignment overlap filtering:", ref_homology_out_df.shape[0], file=log_out
        )
        ref_homology_out_df["sorted_query_target_pair"] = ref_homology_out_df.apply(
            lambda row: tuple(sorted([row[0], row[5]])), axis=1
        )
        ref_homology_out_df.drop_duplicates("sorted_query_target_pair", inplace=True)
        print(
            "Homology .paf shape after query-target symmetric pairs filtering:",
            ref_homology_out_df.shape[0],
            file=log_out,
        )
        ref_homology_out_df["blast_identity"] = ref_homology_out_df[9] / ref_homology_out_df[10]

        minimap2_most_similar_region_tuple_list = []
        percent_id_idxmax = ref_homology_out_df.groupby(0)["blast_identity"].idxmax()
        for tcr_1, idx in zip(percent_id_idxmax.index, percent_id_idxmax):
            minimap2_most_similar_region_tuple_list.append(
                (tcr_1, ref_homology_out_df.loc[idx, 5], ref_homology_out_df.loc[idx, "blast_identity"])
            )
        print(
            "Median blast identity of most similar regions:",
            str(np.median([tuple[2] for tuple in minimap2_most_similar_region_tuple_list])),
            file=log_out,
        )
        print(
            "0.925 quantile blast identity of most similar regions:",
            str(np.quantile([tuple[2] for tuple in minimap2_most_similar_region_tuple_list], 0.925)),
            file=log_out,
        )
        print(
            "0.950 quantile blast identity of most similar regions:",
            str(np.quantile([tuple[2] for tuple in minimap2_most_similar_region_tuple_list], 0.950)),
            file=log_out,
        )
        print(
            "0.975 quantile blast identity of most similar regions:",
            str(np.quantile([tuple[2] for tuple in minimap2_most_similar_region_tuple_list], 0.975)),
            file=log_out,
        )
        print(
            "0.990 quantile blast identity of most similar regions:",
            str(np.quantile([tuple[2] for tuple in minimap2_most_similar_region_tuple_list], 0.990)),
            file=log_out,
        )
        print(
            "Maximal blast identity of most similar regions:",
            str(np.max([tuple[2] for tuple in minimap2_most_similar_region_tuple_list])),
            file=log_out,
        )

    # TODO: minimap2_most_similar_region_dict now does not contain blast_id distances for all TCR pairs,
    # as some TCRs cannot be aligneed with alignment length of min_tcr_length*0.99
    # should we remove min_tcr_length filtering to allow all pairs to be stored in minimap2_most_similar_region_dict?
    # Not sure...
    minimap2_most_similar_region_dict = collections.defaultdict(list)
    for region1, region2, blast_id in minimap2_most_similar_region_tuple_list:
        minimap2_most_similar_region_dict[region1].append(blast_id)
        minimap2_most_similar_region_dict[region2].append(blast_id)

    with open(minimap2_most_similar_region_dict_json, "w") as json_out:
        json.dump(minimap2_most_similar_region_dict, json_out, indent=4)

    similar_region_clusters_list = greedy_minimap2_most_similar_region_clustering(
        minimap2_most_similar_region_tuple_list, similarity_threshold=id_cluster_threshold
    )
    # print(len(similar_region_clusters_list))
    # with open(region_cluster_list_txt, 'w') as txt_out:
    #     for cluster_set in similar_region_clusters_list:
    #         txt_out.write(str(cluster_set) + '\n')

    region_cluster_dict = {}

    cluster_idx = 0
    if similar_region_clusters_list:
        for cluster_set in similar_region_clusters_list:
            for region in cluster_set:
                region_cluster_dict[region] = cluster_idx
            cluster_idx += 1

    regions_set = generate_regions_set_from_ref_fa(reference=reference)

    for region in regions_set:
        if region in region_cluster_dict:
            continue
        else:
            region_cluster_dict[region] = cluster_idx
            cluster_idx += 1

    if len(regions_set.intersection(set(region_cluster_dict.keys()))) != len(regions_set):
        raise Exception(
            "Regions are missing from the final_region_cluster_split_list!"
            "Reference region minimap2 self-homology map clustering is broken!"
        )

    with open(region_cluster_dict_json, "w") as json_out:
        json.dump(region_cluster_dict, json_out, indent=4)

    return region_cluster_dict_json, np.max([tuple[2] for tuple in minimap2_most_similar_region_tuple_list])


def filter_and_split_reads_by_region_cluster(
    bam_file: Union[str, os.PathLike[str]],
    region_cluster_dict_json: Union[str, os.PathLike[str]],
    reference: Union[str, os.PathLike[str]],
    logs_dir: Union[str, os.PathLike[str]],
    region_fasta_out_dir: Union[str, os.PathLike[str]],
    minimal_region_overlap: float = 0.95,
    max_softclip_5_end: int = 73,
    max_softclip_3_end: int = 68,
):
    n_unmapped = 0
    n_short = 0
    n_long = 0
    n_primary_mapped = 0

    log_file = os.path.join(
        logs_dir, os.path.basename(bam_file).split(".")[0] + "_filter_and_split_reads_by_region_cluster.err"
    )
    # flagstat_file = os.path.join(logs_dir, os.path.basename(bam_file).split('.')[0] + '_flagstat.out')

    # with open(flagstat_file, 'r') as flag_in:
    #     lines = flag_in.readlines()
    #     n_primary = int(lines[1].split(' +')[0])
    #     n_primary_mapped = int(lines[7].split(' +')[0])

    with open(region_cluster_dict_json, "r") as json_in:
        region_cluster_dict = json.load(json_in)

    region_length_dict = generate_region_length_dict_from_ref_fa(reference=reference)

    region_cluster_fasta_set = set()
    detected_regions_set = set()
    n_reads_region_cluster_counter = collections.defaultdict(int)
    with pysam.AlignmentFile(bam_file, "rb") as bam_in:
        for entry in bam_in:
            if entry.is_unmapped:
                n_unmapped += 1
                continue
            if entry.is_secondary or entry.is_supplementary:
                continue
            n_primary_mapped += 1

            if entry.reference_length < (region_length_dict[entry.reference_name] * minimal_region_overlap):
                n_short += 1
                continue
            if entry.query_length > (
                region_length_dict[entry.reference_name] * (2 - minimal_region_overlap)
                + (max_softclip_5_end + max_softclip_3_end)
            ):
                n_long += 1
                continue

            region_cluster = region_cluster_dict[entry.reference_name]
            n_reads_region_cluster_counter[region_cluster] += 1
            region_cluster_fasta = os.path.join(region_fasta_out_dir, "region_cluster{}.fasta".format(region_cluster))
            with open(region_cluster_fasta, "a") as fasta_out:
                if entry.is_reverse:
                    print(">{};strand={}".format(entry.query_name, "-"), file=fasta_out)
                    print(entry.get_forward_sequence(), file=fasta_out)
                else:
                    print(">{};strand={}".format(entry.query_name, "+"), file=fasta_out)
                    print(entry.query_sequence, file=fasta_out)

            region_cluster_fasta_set.add(region_cluster_fasta)
            detected_regions_set.add(entry.reference_name)

    logging_str = "Total # primary alignments in bam file: " + str(n_primary_mapped) + "\n"
    # logging_str += '% unmapped reads of all primary reads: ' + str(round(n_unmapped/n_primary, 3)) + '\n'
    logging_str += (
        "median # of primary alignments in region clusters that have minimal region overlap and are not too long: "
        + str(round(np.median(list(n_reads_region_cluster_counter.values())), 3))
        + "\n"
    )
    logging_str += (
        "% of primary alignments that have shorter overlap than minimal region overlap: "
        + str(round(100 * n_short / n_primary_mapped, 2))
        + "\n"
    )
    logging_str += (
        "% of primary alignments that have too long reads: " + str(round(100 * n_long / n_primary_mapped, 2)) + "\n"
    )
    # logging_str += '% of primary alignments that are concatamers: ' + str(round(100 * n_concatamer / n_primary_mapped, 2)) + '\n'

    regions_set = generate_regions_set_from_ref_fa(
        reference=reference,
        remove_negative_control_regions=True,
        negative_control_regions_suffix_str_tuple=("_v_n", "cdr3j_n", "full_n"),
    )

    detected_regions_set = set(
        [region for region in detected_regions_set if not region.endswith(("_v_n", "cdr3j_n", "full_n"))]
    )
    fraction_regions_detected = len(regions_set.intersection(detected_regions_set)) / len(regions_set)
    number_missing_regions = len(regions_set.difference(detected_regions_set))
    missing_regions = regions_set.difference(detected_regions_set)
    logging_str += (
        "fraction detected regions of total regions in reference in initial non-polished read alignments: "
        + str(round(fraction_regions_detected, 4))
        + "\n"
    )
    logging_str += (
        "# of missing regions from reference in initial non-polished read alignments: "
        + str(number_missing_regions)
        + "\n"
    )
    logging_str += (
        "missing/non-detected regions from reference in initial non-polished read alignments: "
        + str(missing_regions)
        + "\n"
    )

    with open(log_file, "w") as ferr:
        ferr.write(logging_str)

    return list(region_cluster_fasta_set)


def filter_and_split_reads_by_region(
    bam_file: Union[str, os.PathLike[str]],
    reference: Union[str, os.PathLike[str]],
    logs_dir: Union[str, os.PathLike[str]],
    region_fasta_out_dir: Union[str, os.PathLike[str]],
    minimal_region_overlap: float = 0.95,
    max_softclip_5_end: int = 73,
    max_softclip_3_end: int = 68,
):
    n_unmapped = 0
    n_short = 0
    n_long = 0
    n_primary_mapped = 0

    log_file = os.path.join(
        logs_dir, os.path.basename(bam_file).split(".")[0] + "_filter_and_split_reads_by_region.err"
    )
    # flagstat_file = os.path.join(logs_dir, os.path.basename(bam_file).split('.')[0] + '_flagstat.out')

    # with open(flagstat_file, 'r') as flag_in:
    #     lines = flag_in.readlines()
    #     n_primary = int(lines[1].split(' +')[0])
    #     n_primary_mapped = int(lines[7].split(' +')[0])

    region_length_dict = generate_region_length_dict_from_ref_fa(reference=reference)

    region_fasta_set = set()
    detected_regions_set = set()
    n_reads_region_counter = collections.defaultdict(int)
    with pysam.AlignmentFile(bam_file, "rb") as bam_in:
        for entry in bam_in:
            if entry.is_unmapped:
                n_unmapped += 1
                continue
            if entry.is_secondary or entry.is_supplementary:
                continue
            n_primary_mapped += 1

            if entry.reference_length < (region_length_dict[entry.reference_name] * minimal_region_overlap):
                n_short += 1
                continue
            if entry.query_length > (
                region_length_dict[entry.reference_name] * (2 - minimal_region_overlap)
                + (max_softclip_5_end + max_softclip_3_end)
            ):
                n_long += 1
                continue

            n_reads_region_counter[entry.reference_name] += 1
            region_fasta = os.path.join(region_fasta_out_dir, "region_{}.fasta".format(entry.reference_name))
            with open(region_fasta, "a") as fasta_out:
                if entry.is_reverse:
                    print(">{};strand={}".format(entry.query_name, "-"), file=fasta_out)
                    print(entry.get_forward_sequence(), file=fasta_out)
                else:
                    print(">{};strand={}".format(entry.query_name, "+"), file=fasta_out)
                    print(entry.query_sequence, file=fasta_out)

            region_fasta_set.add(region_fasta)
            detected_regions_set.add(entry.reference_name)

    logging_str = "Total # primary alignments in bam file: " + str(n_primary_mapped) + "\n"
    # logging_str += '% unmapped reads of all primary reads: ' + str(round(n_unmapped/n_primary, 3)) + '\n'
    logging_str += (
        "median # of primary alignments in regions that have minimal region overlap and are not too long: "
        + str(round(np.median(list(n_reads_region_counter.values())), 3))
        + "\n"
    )
    logging_str += (
        "% of primary alignments that have shorter overlap than minimal region overlap: "
        + str(round(100 * n_short / n_primary_mapped, 2))
        + "\n"
    )
    logging_str += (
        "% of primary alignments that have too long reads: " + str(round(100 * n_long / n_primary_mapped, 2)) + "\n"
    )
    # logging_str += '% of primary alignments that are concatamers: ' + str(round(100 * n_concatamer / n_primary_mapped, 2)) + '\n'

    regions_set = generate_regions_set_from_ref_fa(
        reference=reference,
        remove_negative_control_regions=True,
        negative_control_regions_suffix_str_tuple=("_v_n", "cdr3j_n", "full_n"),
    )

    detected_regions_set = set(
        [region for region in detected_regions_set if not region.endswith(("_v_n", "cdr3j_n", "full_n"))]
    )
    fraction_regions_detected = len(regions_set.intersection(detected_regions_set)) / len(regions_set)
    number_missing_regions = len(regions_set.difference(detected_regions_set))
    missing_regions = regions_set.difference(detected_regions_set)
    logging_str += (
        "fraction detected regions of total regions in reference: " + str(round(fraction_regions_detected, 4)) + "\n"
    )
    logging_str += "# of missing regions from reference: " + str(number_missing_regions) + "\n"
    logging_str += "missing/non-detected regions from reference: " + str(missing_regions) + "\n"

    with open(log_file, "w") as ferr:
        ferr.write(logging_str)

    return list(region_fasta_set)


def filter_and_split_reads_by_region_old(
    bam_file: Union[str, os.PathLike[str]],
    reference: Union[str, os.PathLike[str]],
    bed_ref_file: Union[str, os.PathLike[str]],
    logs_dir: Union[str, os.PathLike[str]],
    region_fastq_out_dir: Union[str, os.PathLike[str]],
    minimal_region_overlap: float = 0.95,
    max_softclip_5_end: int = 73,
    max_softclip_3_end: int = 68,
):
    n_short = 0
    n_long = 0

    log_file = os.path.join(
        logs_dir, os.path.basename(bam_file).split(".")[0] + "_filter_and_split_reads_by_region.err"
    )
    flagstat_file = os.path.join(logs_dir, os.path.basename(bam_file).split(".")[0] + "_flagstat.out")

    with open(flagstat_file, "r") as flag_in:
        lines = flag_in.readlines()
        n_primary_mapped = int(lines[7].split(" +")[0])

    logging_str = "Total # primary alignments in bam file: " + str(n_primary_mapped) + "\n"

    region_fastq_list = []
    detected_region_list = []
    n_reads_region_list = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for region in parse_bed(bed_ref_file):
            logging_str += "Processing region: " + region["name"] + "\n"
            n_reads_region = 0
            region_fastq = os.path.join(region_fastq_out_dir, "{}.fastq".format(region["name"]))
            with open(region_fastq, "w") as out:
                region_length = region["end"] - region["start"]
                for read in bam.fetch(contig=region["chr"], start=region["start"], stop=region["end"]):
                    if read.is_secondary:
                        continue
                    if read.is_supplementary:
                        continue

                    # if read.query_alignment_length < (read.query_length - (max_softclip_5_end + max_softclip_3_end)):
                    #     n_concatamer += 1
                    #     continue

                    if read.reference_length < (region_length * minimal_region_overlap):
                        n_short += 1
                        continue

                    if read.query_length > (
                        region_length * (2 - minimal_region_overlap) + (max_softclip_5_end + max_softclip_3_end)
                    ):
                        n_long += 1
                        continue

                    n_reads_region += 1
                    if read.is_reverse:
                        print("@{};strand={}".format(read.query_name, "-"), file=out)
                        print(read.get_forward_sequence(), file=out)
                        print("+", file=out)
                        print(pysam.qualities_to_qualitystring(read.get_forward_qualities()), file=out)
                    else:
                        print("@{};strand={}".format(read.query_name, "+"), file=out)
                        print(read.query_sequence, file=out)
                        print("+", file=out)
                        print(pysam.qualities_to_qualitystring(read.query_qualities), file=out)

            if n_reads_region > 0:
                region_fastq_list.append(region_fastq)
                detected_region_list.append(region["name"])
                n_reads_region_list.append(n_reads_region)
            #     logging_str += "# primary alignments with above minimal region overlap length: {}".format(n_reads_region) + '\n'
            # else:
            #     logging_str += "region has no primary alignments with above minimal region overlap length! It will be skipped in further processing!\n"

    logging_str += (
        "median # of primary alignments that have minimal region overlap and are not too long: "
        + str(round(np.median(n_reads_region_list), 3))
        + "\n"
    )
    logging_str += (
        "% of primary alignments that have shorter overlap than minimal region overlap: "
        + str(round(100 * n_short / n_primary_mapped, 2))
        + "\n"
    )
    logging_str += (
        "% of primary alignments that have too long reads: " + str(round(100 * n_long / n_primary_mapped, 2)) + "\n"
    )
    # logging_str += '% of primary alignments that are concatamers: ' + str(round(100 * n_concatamer / n_primary_mapped, 2)) + '\n'

    ref_names_set = set()
    with pysam.FastxFile(reference, "r") as fasta_in:
        for entry in fasta_in:
            ref_names_set.add(entry.name)

    detected_regions_set = set(detected_region_list)
    fraction_regions_detected = len(ref_names_set.intersection(detected_regions_set)) / len(ref_names_set)
    number_missing_regions = len(ref_names_set.difference(detected_regions_set))
    missing_regions = ref_names_set.difference(detected_regions_set)
    logging_str += "fraction detected regions of total regions in reference: " + str(fraction_regions_detected) + "\n"
    logging_str += "# of missing regions from reference: " + str(number_missing_regions) + "\n"
    logging_str += "missing/non-detected regions from reference: " + str(missing_regions) + "\n"

    with open(log_file, "w") as ferr:
        ferr.write(logging_str)

    return region_fastq_list
