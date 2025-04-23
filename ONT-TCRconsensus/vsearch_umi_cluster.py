import os
import subprocess
from typing import Union

import ray


@ray.remote
def vsearch_cluster(
    umi_fasta: Union[str, os.PathLike[str]],
    clustering_out_dir: Union[str, os.PathLike[str]],
    threads: int,
    min_umi_length: int = 50,
    max_umi_length: int = 60,
    identity: float = 0.94,
):
    consensus_umi_fasta = os.path.join(clustering_out_dir, "umi_clusters_consensus.fasta")
    # vsearch_clusters_uc = os.path.join(clustering_out_dir, 'vsearch_clusters.tsv')
    log_file = os.path.join(clustering_out_dir, "vsearch_cluster.log")

    subprocess.run(
        [
            "vsearch",
            "--clusterout_id",
            "--clusters",
            clustering_out_dir + "/cluster",
            "--consout",
            consensus_umi_fasta,
            "--minseqlength",
            str(min_umi_length),
            "--maxseqlength",
            str(max_umi_length),
            "--threads",
            str(threads),
            "--cluster_fast",
            umi_fasta,
            "--strand",
            "both",
            #'--uc', vsearch_clusters_uc,
            "--log",
            log_file,
            "--quiet",
            "--no_progress",
            "--clusterout_sort",
            "--gapopen",
            "0E/40I",
            "--mismatch",
            "-40",
            "--match",
            "10",
            "--id",
            str(identity),
        ]
    )

    return consensus_umi_fasta


@ray.remote
def vsearch_cluster_consensus(
    umi_fasta: Union[str, os.PathLike[str]],
    clustering_consensus_out_dir: Union[str, os.PathLike[str]],
    threads: int,
    min_umi_length: int = 50,
    max_umi_length: int = 60,
    identity: float = 0.97,
):
    consensus_umi_fasta = os.path.join(clustering_consensus_out_dir, "umi_clusters_consensus.fasta")
    log_file = os.path.join(clustering_consensus_out_dir, "vsearch_cluster_consensus.log")

    subprocess.run(
        [
            "vsearch",
            "--clusterout_id",
            "--clusters",
            clustering_consensus_out_dir + "/cluster",
            "--consout",
            consensus_umi_fasta,
            "--minseqlength",
            str(min_umi_length),
            "--maxseqlength",
            str(max_umi_length),
            "--threads",
            str(threads),
            "--cluster_fast",
            umi_fasta,
            "--strand",
            "both",
            "--log",
            log_file,
            "--quiet",
            "--no_progress",
            "--clusterout_sort",
            "--id",
            str(identity),
        ]
    )

    return consensus_umi_fasta
