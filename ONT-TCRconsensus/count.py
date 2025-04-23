import os
import subprocess
from typing import Union

import pandas as pd
import ray


@ray.remote(num_cpus=1)
def count_sequences_in_fasta(smolecule_fa):
    try:
        result = subprocess.run(["grep", "-c", "^>", smolecule_fa], capture_output=True, text=True)

        if result.returncode == 0:
            count = int(result.stdout.strip())
            return count
        else:
            raise RuntimeError(f"Error processing {smolecule_fa}: {result.stderr.strip()}")
    except Exception as e:
        raise RuntimeError(f"Exception occurred while processing {smolecule_fa}: {e}")


def umi_counter(smolecule_filtered_fa_list: list):
    region_umi_counter = {}

    count_sequences_in_fasta_futures = []

    for smolecule_fa in smolecule_filtered_fa_list:
        count_sequences_in_fasta_futures.append(count_sequences_in_fasta.remote(smolecule_fa))

    umi_count_list = ray.get(count_sequences_in_fasta_futures)

    for smolecule_fa, count in zip(smolecule_filtered_fa_list, umi_count_list):
        region_umi_counter[os.path.dirname(smolecule_fa).split("region_")[1]] = count

    return region_umi_counter


def write_region_umi_counts_to_csv(
    smolecule_filtered_fa_list: list,
    library_umi_counts_dir: Union[str, os.PathLike[str]],
    region_name: str = "TCR",
):
    umi_counts_csv = os.path.join(library_umi_counts_dir, "umi_consensus_counts.csv")
    region_umi_counter = umi_counter(smolecule_filtered_fa_list=smolecule_filtered_fa_list)

    counts_df = pd.DataFrame.from_dict(region_umi_counter, orient="index", columns=["Count"])
    counts_df.index.name = region_name
    counts_df.reset_index(inplace=True)

    counts_df.to_csv(umi_counts_csv, index=False)
