import collections
import json
import os
from typing import TextIO, Union

import pysam
import ray


def polish_cluster(
    id_cluster: int,
    clustering_out_dir: Union[str, os.PathLike[str]],
    polish_cluster_out_dir: Union[str, os.PathLike[str]],
    stat_out: TextIO,
    smolecule_out: TextIO,
    logging_str: str,
    min_reads_per_cluster: int = 20,
    max_reads_per_cluster: int = 60,
    balance_strands: bool = False,
    cons_umi: str = None,
):
    reads_found = 0

    n_fwd = 0
    n_rev = 0
    # reads_skipped = 0

    cluster_full_seq_fasta = os.path.join(polish_cluster_out_dir, "cluster{}.fasta".format(id_cluster))

    reads_fwd = {}
    reads_rev = {}

    required_fwd = max_reads_per_cluster
    required_rev = max_reads_per_cluster

    with pysam.FastxFile(os.path.join(clustering_out_dir, "cluster{}".format(id_cluster))) as fasta_in:
        for entry in fasta_in:
            cols = entry.name.split(";")
            if len(cols) != 7:
                raise Exception(
                    id_cluster,
                    "cluster fasta entry header has",
                    len(cols),
                    "cols while it should contain 7!",
                    entry.name,
                    cols,
                )
            strand = cols[1].split("strand=")[1]
            read_id = cols[0]

            # TODO: compute edit distance between cons_umi and umi. Ignore read if too large
            reads_found += 1

            if strand == "+":
                if required_fwd > 0:
                    reads_fwd[read_id] = entry
                    required_fwd -= 1
                n_fwd += 1
            elif strand == "-":
                if required_rev > 0:
                    reads_rev[read_id] = entry
                    required_rev -= 1
                n_rev += 1
            else:
                raise Exception("Strand annotation is", strand, "but only - or + are allowed!")

    if balance_strands:
        min_fwd = int(min_reads_per_cluster / 2)
        min_rev = int(min_reads_per_cluster / 2)
        max_reads_after_balance = min(n_fwd * 2, n_rev * 2, max_reads_per_cluster)
        max_fwd = int(max_reads_after_balance / 2)
        max_rev = int(max_reads_after_balance / 2)

    else:
        min_fwd = 0
        min_rev = 0

        if n_fwd > n_rev:
            max_rev = min(n_rev, int(max_reads_per_cluster / 2))
            max_fwd = min(max_reads_per_cluster - max_rev, n_fwd)
        else:
            max_fwd = min(n_fwd, int(max_reads_per_cluster / 2))
            max_rev = min(max_reads_per_cluster - max_fwd, n_rev)

    n_reads = max_fwd + max_rev
    if n_reads > max_reads_per_cluster:
        raise Exception("n_reads is higher than max_reads_per_cluster! max_fwd and max_rev calculation is incorrect!")

    logging_str += "Cluster: {} has {}/{} fwd and {}/{} rev reads\n".format(
        cluster_full_seq_fasta, n_fwd, max_fwd, n_rev, max_rev
    )

    if n_fwd >= min_fwd and n_rev >= min_rev and n_reads >= min_reads_per_cluster:
        entries_fwd = list(reads_fwd.values())[:max_fwd]
        entries_rev = list(reads_rev.values())[:max_rev]
        entries_to_write = entries_fwd + entries_rev
        entries_to_write = entries_to_write[:max_reads_per_cluster]  # Final cap

        reads_written_fwd = len(entries_fwd)
        reads_written_rev = len(entries_rev)
        reads_written_clusters = len(entries_to_write)
        cluster_written = 1

        with open(cluster_full_seq_fasta, "w") as cluster_full_seq_out:
            for entry in entries_to_write:
                cols = entry.name.split(";")
                strand = cols[1].split("strand=")[1]
                read_seq = cols[6].split("seq=")[1]
                read_id = cols[0]

                print(f">{read_id}", file=cluster_full_seq_out)
                print(read_seq, file=cluster_full_seq_out)

                if smolecule_out:
                    print(f">{id_cluster}", file=smolecule_out)
                    print(read_seq, file=smolecule_out)
    else:
        cluster_written = 0
        reads_written_fwd = 0
        reads_written_rev = 0
        reads_written_clusters = 0
        logging_str += "Cluster {} skipped\n".format(id_cluster)

    logging_str += "Cluster: {} has {} reads written: {} fwd - {} rev\n".format(
        cluster_full_seq_fasta, reads_written_clusters, reads_written_fwd, reads_written_rev
    )

    print(
        "cluster{}".format(id_cluster),
        n_fwd,
        n_rev,
        reads_written_fwd,
        reads_written_rev,
        reads_found,
        reads_written_clusters,
        cluster_written,
        sep="\t",
        file=stat_out,
    )
    return (cluster_written, reads_found, reads_written_clusters, logging_str)


@ray.remote
def parse_umi_clusters(
    consensus_umi_fasta: Union[str, os.PathLike[str]],
    regions_wo_clusters_txt: Union[str, os.PathLike[str]],
    min_reads_per_cluster: int = 20,
    max_reads_per_cluster: int = 60,
    region_cluster_dict_json: Union[str, os.PathLike[str]] = None,
    balance_strands: bool = False,
    max_clusters: int = None,
):
    n_clusters = 0
    n_written = 0
    reads_found = 0
    reads_written_clusters = 0

    if region_cluster_dict_json:
        with open(region_cluster_dict_json, "r") as json_in:
            region_cluster_dict = json.load(json_in)
        inverted_region_cluster_dict = collections.defaultdict(list)
        for region, region_cluster in region_cluster_dict.items():
            inverted_region_cluster_dict[region_cluster].append(region)

    region = os.path.basename(os.path.dirname(consensus_umi_fasta))
    region_library_clustering_dir = os.path.dirname(consensus_umi_fasta)
    polish_cluster_out_dir = os.path.join(region_library_clustering_dir, "clusters_fa")
    smolecule_fa = os.path.join(region_library_clustering_dir, "smolecule_clusters.fa")

    stat_file = os.path.join(region_library_clustering_dir, "vsearch_cluster_stats.tsv")
    log_file = os.path.join(region_library_clustering_dir, "parse_cluster.log")
    logging_str = ""

    if not os.path.exists(polish_cluster_out_dir):
        os.mkdir(polish_cluster_out_dir)
    else:
        raise Exception(polish_cluster_out_dir, "should not exist yet but does exist!")

    with pysam.FastxFile(consensus_umi_fasta) as fasta_in:
        for entry in fasta_in:
            n_clusters += 1

    with open(stat_file, "w") as stat_out, open(smolecule_fa, "w") as smolecule_out:
        print(
            "id_cluster",
            "n_fwd",
            "n_rev",
            "written_fwd",
            "written_rev",
            "n",
            "written",
            "cluster_written",
            sep="\t",
            file=stat_out,
        )

        with pysam.FastxFile(consensus_umi_fasta) as fasta_in:
            for entry in fasta_in:
                cols = entry.name.split(";")
                id_cluster = int(cols[-1].split("=")[1])
                cons_umi = None  # entry.sequence.replace("T", "")

                # if n_cluster < min_read_per_cluster:
                #     continue

                cluster_written, reads_found, reads_written_clusters, logging_str = polish_cluster(
                    id_cluster=id_cluster,
                    clustering_out_dir=region_library_clustering_dir,
                    polish_cluster_out_dir=polish_cluster_out_dir,
                    logging_str=logging_str,
                    min_reads_per_cluster=min_reads_per_cluster,
                    max_reads_per_cluster=max_reads_per_cluster,
                    stat_out=stat_out,
                    smolecule_out=smolecule_out,
                    balance_strands=balance_strands,
                    cons_umi=cons_umi,
                )

                n_written += cluster_written
                reads_found += reads_found
                reads_written_clusters += reads_written_clusters
                if max_clusters and n_written > max_clusters:
                    break

    if n_written == 0 or reads_found == 0:
        with open(regions_wo_clusters_txt, "a") as txt_out:
            if region_cluster_dict_json:
                print(region, inverted_region_cluster_dict[int(region.split("region_cluster")[1])], file=txt_out)
            else:
                print(region, file=txt_out)
        return None
        # raise RuntimeError("ERROR - did not find a single cluster passing the min_reads_per_cluster threshold!")

    logging_str += "Clusters: {}% written ({})\n".format(int(n_written * 100.0 / n_clusters), n_written)
    logging_str += "Reads: {} found\n".format(reads_found)
    logging_str += "Reads: {}% in written clusters\n".format(int(reads_written_clusters * 100.0 / reads_found))

    if logging_str:
        with open(log_file, "w") as log_out:
            log_out.write(logging_str)

    return smolecule_fa
