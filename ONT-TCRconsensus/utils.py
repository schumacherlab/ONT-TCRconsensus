import os
from typing import Union


def init_nanopore_analysis_dir(fastq: Union[str, os.PathLike[str]], nano_dir: Union[str, os.PathLike[str]]):
    library = os.path.basename(fastq).split(".")[0]
    library_dir = os.path.join(nano_dir, library)
    library_logs_dir = os.path.join(library_dir, "logs")
    library_align_dir = os.path.join(library_dir, "align")
    library_region_cluster_fasta_dir_dict = os.path.join(library_dir, "region_cluster_fasta")
    library_umi_fasta_dir = os.path.join(library_dir, "umi_fasta")
    library_clustering_dir = os.path.join(library_dir, "clustering")
    library_fasta_dir = os.path.join(library_dir, "fasta")
    library_clustering_consensus_dir = os.path.join(library_dir, "clustering_consensus")
    library_consensus_region_fasta_dir_dict = os.path.join(library_dir, "region_fasta")
    library_consensus_umi_fasta_dir = os.path.join(library_dir, "consensus_umi_fasta")
    library_umi_counts_dir = os.path.join(library_dir, "counts")
    outs_dir = os.path.join(library_dir, "outs")

    os.mkdir(library_dir)
    os.mkdir(library_logs_dir)
    os.mkdir(library_align_dir)
    os.mkdir(library_region_cluster_fasta_dir_dict)
    os.mkdir(library_umi_fasta_dir)
    os.mkdir(library_clustering_dir)
    os.mkdir(library_fasta_dir)
    os.mkdir(library_clustering_consensus_dir)
    os.mkdir(library_consensus_region_fasta_dir_dict)
    os.mkdir(library_consensus_umi_fasta_dir)
    os.mkdir(library_umi_counts_dir)
    os.mkdir(outs_dir)

    return (
        library_dir,
        library_logs_dir,
        library_align_dir,
        library_region_cluster_fasta_dir_dict,
        library_umi_fasta_dir,
        library_clustering_dir,
        library_fasta_dir,
        library_clustering_consensus_dir,
        library_consensus_region_fasta_dir_dict,
        library_consensus_umi_fasta_dir,
        library_umi_counts_dir,
    )


def dorado_num_cpus_task(total_num_cpus: int, fastq_list=list):
    # TODO: change this to ray.cluster_resources().get("CPU", 0) if this does correctly return the numer of cpus allocated in sbatch
    num_cpus = total_num_cpus // len(fastq_list)
    if num_cpus == 0:
        # TODO: print to general workflow log that it would be more suitable if more cores are allocated to sbatch!
        return 1
    else:
        return num_cpus


def vsearch_umis_num_cpus_task(total_num_cpus: int, futures_umi_extract: list):
    # TODO: change this to ray.cluster_resources().get("CPU", 0) if this does correctly return the numer of cpus allocated in sbatch
    num_cpus = total_num_cpus // len(futures_umi_extract)
    if num_cpus <= 25:
        # TODO: print to general workflow log that it would be more suitable if more cores are allocated to sbatch!
        return 25  # I have no idea whether 25 cores is a good idea...
    else:
        return num_cpus  # This is very unlikely to happen for large TCR library sizes
