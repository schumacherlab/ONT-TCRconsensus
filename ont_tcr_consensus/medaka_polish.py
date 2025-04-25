import os
import subprocess
from typing import Union

import numpy as np
import pandas as pd
import pysam
import ray


def medaka_memory_task(
    smolecule_fa: Union[str, os.PathLike[str]],
    # memory_gb_per_umi_cluster: float = 0.0015,
    memory_gb_task_overhead: float = 2.75,
):
    vsearch_cluster_stats_tsv = os.path.join(os.path.dirname(smolecule_fa), "vsearch_cluster_stats.tsv")
    vsearch_cluster_stats_df = pd.read_csv(vsearch_cluster_stats_tsv, sep="\t")
    vsearch_cluster_stats_df = vsearch_cluster_stats_df.loc[vsearch_cluster_stats_df.loc[:, "cluster_written"] == 1, :]
    max_subreads_umi = vsearch_cluster_stats_df.loc[:, "written"].max()
    mem_gb_per_umi_cluster = 0.0143 * max_subreads_umi + 0.0286 # constants were optimized for Intel Xeon Silver

    with open(smolecule_fa, "r") as fasta_in:
        num_umis_per_smolecule_fa = len({line.split("_")[0] for line in fasta_in if line.startswith(">")})

    medaka_memory_bytes = (mem_gb_per_umi_cluster * num_umis_per_smolecule_fa + memory_gb_task_overhead) * 1024**3
    return medaka_memory_bytes, mem_gb_per_umi_cluster


@ray.remote
def generate_medaka_smolecule_fa_batches(
    smolecule_fa: Union[str, os.PathLike[str]],
    max_cap_medaka_memory_gb: float,
    # memory_gb_per_umi_cluster: float = 0.0015,
    memory_gb_task_overhead: float = 2.75,
):
    medaka_memory_bytes, memory_gb_per_umi_cluster = medaka_memory_task(
        smolecule_fa=smolecule_fa,
        # memory_gb_per_umi_cluster = memory_gb_per_umi_cluster,
        memory_gb_task_overhead=memory_gb_task_overhead,
    )

    if medaka_memory_bytes < max_cap_medaka_memory_gb * 1024**3:
        return [(smolecule_fa, medaka_memory_bytes)]

    clustering_dir = os.path.dirname(smolecule_fa)
    seen_cluster_id_list = []
    with pysam.FastxFile(smolecule_fa) as fasta_in:
        for entry in fasta_in:
            cluster_id = int(entry.name.split("_")[0])
            if cluster_id not in seen_cluster_id_list:
                seen_cluster_id_list.append(cluster_id)

    cum_memory_arr = np.cumsum(np.full(len(seen_cluster_id_list), memory_gb_per_umi_cluster))
    cluster_batches = (cum_memory_arr // (max_cap_medaka_memory_gb - memory_gb_task_overhead)).astype(int)
    memory_last_batch = cum_memory_arr[-1] - (
        (max_cap_medaka_memory_gb - memory_gb_task_overhead) * cluster_batches[-1]
    )
    cluster_batches_memory_gb = [max_cap_medaka_memory_gb] * cluster_batches.max()
    cluster_batches_memory_gb.append(memory_last_batch + memory_gb_task_overhead)
    cluster_batches_memory_bytes = [memory_gb * 1024**3 for memory_gb in cluster_batches_memory_gb]
    number_of_batches = cluster_batches.max() + 1
    cluster_batch_map = dict(zip(seen_cluster_id_list, cluster_batches))

    batch_file_handles = {
        batch_idx: open(os.path.join(clustering_dir, f"smolecule_clusters_{batch_idx}.fa"), "w")
        for batch_idx in range(number_of_batches)
    }

    with pysam.FastxFile(smolecule_fa) as fasta_in:
        for entry in fasta_in:
            cluster_id = int(entry.name.split("_")[0])
            batch_idx = cluster_batch_map[cluster_id]
            batch_file_handles[batch_idx].write(f">{entry.name}\n{entry.sequence}\n")

    for handle in batch_file_handles.values():
        handle.close()

    os.remove(smolecule_fa)

    return [
        (batch_smolecule_fa.name, required_memory_bytes)
        for batch_smolecule_fa, required_memory_bytes in zip(batch_file_handles.values(), cluster_batches_memory_bytes)
    ]


def bin_medaka_memory_list_entries(medaka_memory_list: list, num_bins: int = 75):
    medaka_memory_arr = np.array(medaka_memory_list)
    bin_edges = np.linspace(medaka_memory_arr.min(), medaka_memory_arr.max(), num_bins + 1)
    bin_indices = np.digitize(medaka_memory_arr, bin_edges, right=True)
    binned_medaka_memory_arr = bin_edges[bin_indices]

    return binned_medaka_memory_arr.tolist()


@ray.remote
def medaka_polish_clusters(
    smolecule_fa: Union[str, os.PathLike[str]],
    medaka_out_dir: Union[str, os.PathLike[str]],
    medaka_model: str,
    threads: int,
):
    region = os.path.basename(os.path.dirname(smolecule_fa))

    if len(os.path.basename(smolecule_fa).split("_")) > 2:
        region = region + "_" + os.path.basename(smolecule_fa).split("_")[-1].split(".")[0]

    medaka_tmp_out_dir = os.path.join(medaka_out_dir, region + "_consensus_tmp")
    log_file = os.path.join(medaka_out_dir, region + "_consensus_smolecule.log")
    cluster_stats_tsv = os.path.join(os.path.dirname(smolecule_fa), "vsearch_cluster_stats.tsv")

    with open(log_file, "w") as log:
        try:
            process = subprocess.run(
                [
                    "medaka",
                    "smolecule",
                    "--model",
                    medaka_model,
                    "--method",
                    "spoa",
                    "--depth",
                    "2",
                    "--length",
                    "50",
                    # '--quiet',
                    "--threads",
                    str(threads),
                    medaka_tmp_out_dir,
                    smolecule_fa,
                ],
                stderr=log,
                timeout=(60 * 60 * 3),
                check=False,
            )

            if process.returncode != 0:
                print(
                    f"Medaka process for region {region} failed with exit code {process.returncode}. Check the log for details."
                )
                return None

        except subprocess.TimeoutExpired:
            print(f"Medaka process for region {region} timed out", file=log)
            return None

    if os.path.isfile(os.path.join(medaka_tmp_out_dir, "consensus.fasta")):
        subprocess.run(
            [
                "mv",
                os.path.join(medaka_tmp_out_dir, "consensus.fasta"),
                os.path.join(medaka_out_dir, region + "_consensus.fasta"),
            ]
        )
        # subprocess.run(['mv', os.path.join(medaka_tmp_out_dir, 'subreads_to_spoa.bam'), os.path.join(medaka_out_dir, region + '_consensus.bam')])
        # subprocess.run(['mv', os.path.join(medaka_tmp_out_dir, 'subreads_to_spoa.bam.bai'), os.path.join(medaka_out_dir, region + '_consensus.bam.bai')])
        cluster_stats_df = pd.read_csv(cluster_stats_tsv, sep="\t")
        with (
            pysam.FastxFile(os.path.join(medaka_out_dir, region + "_consensus.fasta")) as fasta_in,
            open(os.path.join(medaka_out_dir, region + "_consensus_tmp.fasta"), "w") as fasta_out,
        ):
            for entry in fasta_in:
                cluster_id = entry.name.split("_")[0]
                name = (
                    region
                    + "_"
                    + entry.name
                    + "_"
                    + str(
                        cluster_stats_df.loc[
                            cluster_stats_df.loc[:, "id_cluster"].str.split("cluster").str[1] == cluster_id,
                            "written",
                        ].values[0]
                    )
                )
                print(">" + name, file=fasta_out)
                print(entry.sequence, file=fasta_out)
        os.replace(
            os.path.join(medaka_out_dir, region + "_consensus_tmp.fasta"),
            os.path.join(medaka_out_dir, region + "_consensus.fasta"),
        )

        subprocess.run(["rm", "-rf", medaka_tmp_out_dir])
        return os.path.join(medaka_out_dir, region + "_consensus.fasta")
    else:
        return None
