import argparse
import glob
import json
import os
import shutil
import subprocess
import sys

import ray
from nanopore_tcr_consensus.count import write_region_umi_counts_to_csv
from nanopore_tcr_consensus.extract_umis import count_overlapping_umis_between_all_regions, extract_umis
from nanopore_tcr_consensus.medaka_polish import (
    bin_medaka_memory_list_entries,
    generate_medaka_smolecule_fa_batches,
    medaka_polish_clusters,
)
from nanopore_tcr_consensus.minimap2_align import (
    filter_consensus_alignments,
    minimap2_ont_align,
    minimap2_ref_self_homology_map,
)
from nanopore_tcr_consensus.parse_umi_clusters import parse_umi_clusters
from nanopore_tcr_consensus.preprocessing import dorado_trim, vsearch_fastq_filtering
from nanopore_tcr_consensus.region_split import (
    filter_and_split_reads_by_region,
    filter_and_split_reads_by_region_cluster,
    generate_region_split_dict_from_minimap2_homology_paf,
)
from nanopore_tcr_consensus.utils import dorado_num_cpus_task, init_nanopore_analysis_dir, vsearch_umis_num_cpus_task
from nanopore_tcr_consensus.vsearch_umi_cluster import vsearch_cluster, vsearch_cluster_consensus


def main():
    parser = argparse.ArgumentParser(description="Count unique TCR molecule nanopore consensus reads.")
    parser.add_argument("json_config_file", type=str, help="Path to analysis run JSON config file")
    args = parser.parse_args()

    with open(args.json_config_file, "r") as json_in:
        run_config = json.load(json_in)

    reference_file = run_config["reference_file"]
    fastq_pass_dir = run_config["fastq_pass_dir"]
    only_run_ref_self_homology_bool = run_config["only_run_reference_self_homology"]
    dorado_trim_subsample_fastq = run_config["dorado_trim_subsample_fastq"]
    minimal_length = run_config["minimal_length"]
    max_ee_rate_base = run_config["max_ee_rate_base"]
    minimal_region_overlap = run_config["minimal_region_overlap"]
    max_softclip_5_end = run_config["max_softclip_5_end"]
    max_softclip_3_end = run_config["max_softclip_3_end"]
    umi_fwd = run_config["umi_fwd"]
    umi_rev = run_config["umi_rev"]
    max_pattern_dist = run_config["max_pattern_dist"]
    min_umi_length = run_config["min_umi_length"]
    max_umi_length = run_config["max_umi_length"]
    vsearch_identity = run_config["vsearch_identity"]
    min_reads_per_cluster = run_config["min_reads_per_cluster"]
    max_reads_per_cluster = run_config["max_reads_per_cluster"]
    balance_strands = run_config["balance_strands"]
    compare_umi_overlap_between_regions = run_config["compare_umi_overlap_between_regions"]
    overlapping_umi_edit_threshold = run_config["overlapping_umi_edit_threshold"]
    medaka_model = run_config["medaka_model"]
    medaka_memory_gb_per_umi_cluster = run_config["medaka_memory_gb_per_umi_cluster"]
    medaka_memory_gb_task_overhead = run_config["medaka_memory_gb_task_overhead"]
    max_cap_medaka_memory_gb = run_config["max_cap_medaka_memory_gb"]
    vsearch_identity_consensus = run_config["vsearch_identity_consensus"]
    dorado_excutable = run_config["dorado_excutable"]
    nanopore_tcr_seq_primers_fasta = run_config["nanopore_tcr_seq_primers_fasta"]
    id_reference_self_homology_cluster_threshold = 1 - max_ee_rate_base
    minimal_region_overlap_consensus = run_config["minimal_region_overlap_consensus"]
    blast_id_threshold = run_config["blast_id_threshold"]
    delete_tmp_files = run_config["delete_tmp_files"]

    ray.init()

    slurm_job_result = subprocess.run(
        ["scontrol", "show", "job", os.getenv("SLURM_JOB_ID")], capture_output=True, text=True
    )
    total_num_cpus = int(slurm_job_result.stdout.splitlines()[16].split("AllocTRES=cpu=")[1].split(",mem=")[0])
    print("Available cpus:", total_num_cpus, file=sys.stderr)
    if not isinstance(total_num_cpus, int):
        raise Exception(total_num_cpus, "is not an int!")

    nano_dir = os.path.join(fastq_pass_dir, "nano_tcr")
    if os.path.isdir(nano_dir):
        raise Exception("Nanopore TCR UMI consensus dir was already initialized!")
    os.mkdir(nano_dir)

    fastq_list = glob.glob(fastq_pass_dir + "/barcode*/*fastq.gz")

    print("Estimating self-homology map reference", file=sys.stderr)
    ref_homology_paf = minimap2_ref_self_homology_map(reference=reference_file, minimap2_threads=total_num_cpus - 1)

    region_cluster_dict_json, max_blast_id = generate_region_split_dict_from_minimap2_homology_paf(
        ref_homology_paf=ref_homology_paf,
        reference=reference_file,
        id_cluster_threshold=id_reference_self_homology_cluster_threshold,
    )

    if not minimal_region_overlap_consensus:
        minimal_region_overlap_consensus = max_blast_id
    if not blast_id_threshold:
        blast_id_threshold = max_blast_id

    if only_run_ref_self_homology_bool:
        return

    library_dir_dict = {}
    library_logs_dir_dict = {}
    library_align_dir_dict = {}
    library_region_cluster_fasta_dir_dict = {}
    library_umi_fasta_dir_dict = {}
    library_clustering_dir_dict = {}
    library_fasta_dir_dict = {}
    library_clustering_consensus_dir_dict = {}
    library_consensus_region_fasta_dir_dict = {}
    library_consensus_umi_fasta_dir_dict = {}
    library_umi_counts_dir_dict = {}

    for fastq in fastq_list:
        library = os.path.basename(fastq).split(".")[0]

        dirs = init_nanopore_analysis_dir(fastq=fastq, nano_dir=nano_dir)

        library_dir_dict[library] = dirs[0]
        library_logs_dir_dict[library] = dirs[1]
        library_align_dir_dict[library] = dirs[2]
        library_region_cluster_fasta_dir_dict[library] = dirs[3]
        library_umi_fasta_dir_dict[library] = dirs[4]
        library_clustering_dir_dict[library] = dirs[5]
        library_fasta_dir_dict[library] = dirs[6]
        library_clustering_consensus_dir_dict[library] = dirs[7]
        library_consensus_region_fasta_dir_dict[library] = dirs[8]
        library_consensus_umi_fasta_dir_dict[library] = dirs[9]
        library_umi_counts_dir_dict[library] = dirs[10]

    library = None  # Just to be sure that this initialized by the first fastq for loop

    num_dorado_cpus = dorado_num_cpus_task(total_num_cpus=total_num_cpus, fastq_list=fastq_list)

    print("Trimming adapters from reads", file=sys.stderr)
    futures_dorado = []
    for fastq in fastq_list:
        library = os.path.basename(fastq).split(".")[0]
        futures_dorado.append(
            dorado_trim.options(num_cpus=num_dorado_cpus).remote(
                fastq=fastq,
                fastq_out_dir=library_dir_dict[library],
                logs_dir=library_logs_dir_dict[library],
                dorado_excutable=dorado_excutable,
                dorado_threads=num_dorado_cpus,
                nanopore_tcr_seq_primers_fasta=nanopore_tcr_seq_primers_fasta,
                max_reads=dorado_trim_subsample_fastq,
            )
        )

    print("Filtering reads on quality and minimal length", file=sys.stderr)
    futures_vsearch_filter = []
    for fastq_trimmed in futures_dorado:
        futures_vsearch_filter.append(
            vsearch_fastq_filtering.remote(
                fastq=fastq_trimmed,
                library_fastq_out_dir_dict=library_dir_dict,
                library_logs_dir_dict=library_logs_dir_dict,
                max_ee_rate_base=max_ee_rate_base,
                minimal_length=minimal_length,
            )
        )

    fastq_filtered_list = ray.get(futures_vsearch_filter)
    bam_files_list = []
    for filtered_fastq in fastq_filtered_list:
        library = os.path.basename(os.path.dirname(filtered_fastq))
        print("Aligning reads from:", library, file=sys.stderr)
        bam_files_list.append(
            minimap2_ont_align(
                fastx=filtered_fastq,
                minimap2_threads=total_num_cpus - 1,
                reference=reference_file,
                bam_out_dir=library_align_dir_dict[library],
                logs_dir=library_logs_dir_dict[library],
            )
        )
        os.remove(
            os.path.join(
                os.path.dirname(filtered_fastq), os.path.basename(filtered_fastq).rsplit("_filtered", 1)[0] + ".fastq"
            )
        )
        os.remove(filtered_fastq)

    for bam_file in bam_files_list:
        library = os.path.basename(os.path.dirname(os.path.dirname(bam_file)))
        print("Splitting alignments into region clusters:", library, file=sys.stderr)
        region_cluster_fasta_list = filter_and_split_reads_by_region_cluster(
            bam_file=bam_file,
            region_cluster_dict_json=region_cluster_dict_json,
            reference=reference_file,
            logs_dir=library_logs_dir_dict[library],
            region_fasta_out_dir=library_region_cluster_fasta_dir_dict[library],
            minimal_region_overlap=minimal_region_overlap,
            max_softclip_5_end=max_softclip_5_end,
            max_softclip_3_end=max_softclip_3_end,
        )
        os.remove(bam_file)
        os.remove(bam_file + ".bai")

        print("Extracting umis from alignments within region clusters:", library, file=sys.stderr)
        futures_umi_extract = []
        for region_fasta in region_cluster_fasta_list:
            futures_umi_extract.append(
                extract_umis.remote(
                    fastx_file=region_fasta,
                    umi_fasta_out_dir=library_umi_fasta_dir_dict[library],
                    write_region=True,
                    # logs_dir: Union[str, os.PathLike[str]],
                    adapter_length_5_end=max_softclip_5_end,
                    adapter_length_3_end=max_softclip_3_end,
                    max_pattern_dist=max_pattern_dist,
                    umi_fwd=umi_fwd,
                    umi_rev=umi_rev,
                )
            )

        umi_fasta_list = ray.get(futures_umi_extract)
        filtered_umi_fasta_list = [umi_fasta for umi_fasta in umi_fasta_list if umi_fasta]

        num_vsearch_cpus = vsearch_umis_num_cpus_task(
            total_num_cpus=total_num_cpus, futures_umi_extract=futures_umi_extract
        )

        print("Clustering umis within region clusters:", library, file=sys.stderr)
        futures_vsearch = []
        for umi_fasta in filtered_umi_fasta_list:
            region_cluster = os.path.basename(umi_fasta).split("_detected_umis.fasta")[0]
            region_library_clustering_dir = os.path.join(library_clustering_dir_dict[library], region_cluster)
            os.mkdir(region_library_clustering_dir)
            futures_vsearch.append(
                vsearch_cluster.options(num_cpus=num_vsearch_cpus).remote(
                    umi_fasta=umi_fasta,
                    clustering_out_dir=region_library_clustering_dir,
                    threads=num_vsearch_cpus,  #  A thread is not always the same as a cpu... Should I let vsearch figure out available threads itself.
                    min_umi_length=min_umi_length,
                    max_umi_length=max_umi_length,
                    identity=vsearch_identity,
                )
            )

        # consensus_umi_fasta_list = ray.get(futures_vsearch)
        region_clusters_wo_clusters_txt = os.path.join(
            library_logs_dir_dict[library], "region_clusters_wo_umi_clusters.txt"
        )

        print("Parsing umi clusters within region clusters:", library, file=sys.stderr)
        futures_parse_umi_clusters = []
        for consensus_umi_fasta in futures_vsearch:
            futures_parse_umi_clusters.append(
                parse_umi_clusters.options(num_cpus=num_vsearch_cpus).remote(
                    consensus_umi_fasta=consensus_umi_fasta,
                    regions_wo_clusters_txt=region_clusters_wo_clusters_txt,
                    min_reads_per_cluster=min_reads_per_cluster,
                    max_reads_per_cluster=max_reads_per_cluster,
                    region_cluster_dict_json=region_cluster_dict_json,
                    balance_strands=balance_strands,
                    max_clusters=None,
                )
            )

        smolecule_fa_list = ray.get(futures_parse_umi_clusters)
        smolecule_filtered_fa_list = [smolecule_fa for smolecule_fa in smolecule_fa_list if smolecule_fa]
        if not os.path.isfile(region_clusters_wo_clusters_txt):
            # write an empty region_clusters_wo_clusters_txt
            with open(region_clusters_wo_clusters_txt, "w") as txt_out:
                txt_out.write("")

        # if compare_umi_overlap_between_regions:
        #     region_comparisons_w_overlapping_umis_bool_list = count_overlapping_umis_between_all_regions(smolecule_filtered_fa_list = smolecule_filtered_fa_list,
        #                                                                                                  overlapping_umi_edit_threshold = overlapping_umi_edit_threshold,
        #                                                                                                  logs_dir = library_logs_dir_dict[library]
        #                                                                                                 )

        print("Calculating consensus TCR sequences of umi clusters:", library, file=sys.stderr)
        futures_medaka_batching = []
        for smolecule_fa in smolecule_filtered_fa_list:
            futures_medaka_batching.append(
                generate_medaka_smolecule_fa_batches.remote(
                    smolecule_fa=smolecule_fa,
                    max_cap_medaka_memory_gb=max_cap_medaka_memory_gb,
                    # memory_gb_per_umi_cluster = medaka_memory_gb_per_umi_cluster,
                    memory_gb_task_overhead=medaka_memory_gb_task_overhead,
                )
            )

        medaka_smolecule_fa_batch_list = ray.get(futures_medaka_batching)
        smolecule_fa_batch_list = [
            batch_tuple[0] for smolecule_fa_list in medaka_smolecule_fa_batch_list for batch_tuple in smolecule_fa_list
        ]
        medaka_memory_bytes_list = [
            batch_tuple[1] for smolecule_fa_list in medaka_smolecule_fa_batch_list for batch_tuple in smolecule_fa_list
        ]
        binnned_medaka_memory_bytes_list = bin_medaka_memory_list_entries(
            medaka_memory_list=medaka_memory_bytes_list, num_bins=75
        )
        # batch_count_seq_futures = []
        # for batch_smolecule_fa in smolecule_fa_batch_list:
        #     batch_count_seq_futures.append(count_sequences_in_fasta.remote(batch_smolecule_fa))
        # smolecule_fa_batch_entry_count_list = ray.get(batch_count_seq_futures)
        # medaka_memory_dict = {key: (memory, binned_memory, entry_count,
        #                             int(((memory/1024**3)- medaka_memory_gb_task_overhead) / medaka_memory_gb_per_umi_cluster))
        #                       for key, memory, binned_memory, entry_count
        #                       in zip(smolecule_fa_batch_list,
        #                              medaka_memory_bytes_list,
        #                              binnned_medaka_memory_bytes_list,
        #                              smolecule_fa_batch_entry_count_list
        #                              )}
        # with open(os.path.join(library_logs_dir_dict[library], 'medaka_memory_bytes_dict.json'), "w") as json_out:
        #     json.dump(medaka_memory_dict, json_out, indent=4)

        futures_medaka = []
        for batch_smolecule_fa, medaka_memory_bytes in zip(smolecule_fa_batch_list, binnned_medaka_memory_bytes_list):
            futures_medaka.append(
                medaka_polish_clusters.options(num_cpus=2, memory=medaka_memory_bytes).remote(
                    smolecule_fa=batch_smolecule_fa,
                    medaka_out_dir=library_fasta_dir_dict[library],
                    medaka_model=medaka_model,
                    threads=2,
                )
            )
        consensus_fasta_list = ray.get(futures_medaka)

        incomplete_smolecule_region_clusters = []
        for smolecule_fa, result in zip(smolecule_fa_batch_list, consensus_fasta_list):
            if result is None:
                incomplete_smolecule_region_clusters.append(os.path.basename(os.path.dirname(smolecule_fa)))

        if incomplete_smolecule_region_clusters:
            print(
                "There are:",
                str(len(incomplete_smolecule_region_clusters)),
                "medaka smolecule tasks that did not run to completion!",
                file=sys.stderr,
            )
            print(
                "The following region clusters files did not complete:",
                incomplete_smolecule_region_clusters,
                file=sys.stderr,
            )

        filtered_consensus_fasta_list = [consensus_fasta for consensus_fasta in consensus_fasta_list if consensus_fasta]

        merged_consensus_fasta = os.path.join(library_fasta_dir_dict[library], "merged_consensus.fasta")
        with open(merged_consensus_fasta, "a") as fasta_out:
            for consensus_fasta in filtered_consensus_fasta_list:
                with open(consensus_fasta, "r") as fasta_in:
                    fasta_out.write(fasta_in.read())

        print("Aligning unique molecule consensus TCR sequences:", library, file=sys.stderr)
        consensus_bam_file = minimap2_ont_align(
            fastx=merged_consensus_fasta,
            minimap2_threads=total_num_cpus - 1,
            reference=reference_file,
            bam_out_dir=library_align_dir_dict[library],
            logs_dir=library_logs_dir_dict[library],
        )

        filtered_consensus_bam_file = filter_consensus_alignments(
            consensus_bam_file=consensus_bam_file,
            reference=reference_file,
            logs_dir=library_logs_dir_dict[library],
            blast_id_threshold=blast_id_threshold,
            minimal_region_overlap=minimal_region_overlap_consensus,
            max_softclip_5_end=max_softclip_5_end,
            max_softclip_3_end=max_softclip_3_end,
        )

        # Below we will check whether initial self-homology clustering on the reference worked to bin all unique molecules in the same region_clusters for consensus polishing:
        print("Splitting consensus alignments into regions:", library, file=sys.stderr)
        region_fasta_list = filter_and_split_reads_by_region(
            bam_file=filtered_consensus_bam_file,
            reference=reference_file,
            logs_dir=library_logs_dir_dict[library],
            region_fasta_out_dir=library_consensus_region_fasta_dir_dict[library],
            minimal_region_overlap=minimal_region_overlap_consensus,
            max_softclip_5_end=max_softclip_5_end,
            max_softclip_3_end=max_softclip_3_end,
        )

        print("Extracting umis from consensus alignments within regions:", library, file=sys.stderr)
        futures_umi_extract = []
        for region_fasta in region_fasta_list:
            futures_umi_extract.append(
                extract_umis.remote(
                    fastx_file=region_fasta,
                    umi_fasta_out_dir=library_consensus_umi_fasta_dir_dict[library],
                    write_region=True,
                    # logs_dir: Union[str, os.PathLike[str]],
                    adapter_length_5_end=max_softclip_5_end,
                    adapter_length_3_end=max_softclip_3_end,
                    max_pattern_dist=max_pattern_dist,
                    umi_fwd=umi_fwd,
                    umi_rev=umi_rev,
                )
            )

        umi_fasta_list = ray.get(futures_umi_extract)
        filtered_umi_fasta_list = [umi_fasta for umi_fasta in umi_fasta_list if umi_fasta]

        num_vsearch_cpus = vsearch_umis_num_cpus_task(
            total_num_cpus=total_num_cpus, futures_umi_extract=futures_umi_extract
        )

        print("Clustering umis within regions:", library, file=sys.stderr)
        futures_vsearch = []
        for umi_fasta in filtered_umi_fasta_list:
            region = os.path.basename(umi_fasta).split("_detected_umis.fasta")[0]
            region_library_clustering_consensus_dir = os.path.join(
                library_clustering_consensus_dir_dict[library], region
            )
            os.mkdir(region_library_clustering_consensus_dir)
            futures_vsearch.append(
                vsearch_cluster_consensus.options(num_cpus=num_vsearch_cpus).remote(
                    umi_fasta=umi_fasta,
                    clustering_consensus_out_dir=region_library_clustering_consensus_dir,
                    threads=num_vsearch_cpus,  #  A thread is not always the same as a cpu... Should I let vsearch figure out available threads itself.
                    min_umi_length=min_umi_length,
                    max_umi_length=max_umi_length,
                    identity=vsearch_identity_consensus,
                )
            )

        # consensus_umi_fasta_list = ray.get(futures_vsearch)
        regions_wo_clusters_txt = os.path.join(library_logs_dir_dict[library], "regions_wo_umi_clusters.txt")

        print("Parsing umi clusters within regions:", library, file=sys.stderr)
        futures_parse_umi_clusters = []
        for consensus_umi_fasta in futures_vsearch:
            futures_parse_umi_clusters.append(
                parse_umi_clusters.options(num_cpus=num_vsearch_cpus).remote(
                    consensus_umi_fasta=consensus_umi_fasta,
                    regions_wo_clusters_txt=regions_wo_clusters_txt,
                    min_reads_per_cluster=1,
                    max_reads_per_cluster=max_reads_per_cluster,
                    balance_strands=False,
                    max_clusters=None,
                )
            )

        smolecule_fa_list = ray.get(futures_parse_umi_clusters)
        smolecule_filtered_fa_list = [smolecule_fa for smolecule_fa in smolecule_fa_list if smolecule_fa]
        if not os.path.isfile(regions_wo_clusters_txt):
            # write an empty regions_wo_clusters_txt
            with open(regions_wo_clusters_txt, "w") as txt_out:
                txt_out.write("")

        print("Counting consensus UMIs:", library, file=sys.stderr)
        write_region_umi_counts_to_csv(
            smolecule_filtered_fa_list=smolecule_filtered_fa_list,
            library_umi_counts_dir=library_umi_counts_dir_dict[library],
            region_name="TCR",
        )

        if compare_umi_overlap_between_regions:
            print("Testing for consensus umi matches between regions:", library, file=sys.stderr)
            region_comparisons_w_overlapping_umis_bool_list = count_overlapping_umis_between_all_regions(
                smolecule_filtered_fa_list=smolecule_filtered_fa_list,
                overlapping_umi_edit_threshold=overlapping_umi_edit_threshold,
                logs_dir=library_logs_dir_dict[library],
            )
        if delete_tmp_files:
            shutil.rmtree(library_region_cluster_fasta_dir_dict[library])
            shutil.rmtree(library_clustering_dir_dict[library])
            shutil.rmtree(library_umi_fasta_dir_dict[library])
            shutil.rmtree(library_fasta_dir_dict[library])
            shutil.rmtree(library_clustering_consensus_dir_dict[library])
            shutil.rmtree(library_consensus_region_fasta_dir_dict[library])
            shutil.rmtree(library_consensus_umi_fasta_dir_dict[library])
            print("\n\n", file=sys.stderr)

    print("Done running all barcodes!", file=sys.stderr)
    return


if __name__ == "__main__":
    main()
