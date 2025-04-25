# ONT-TCRconsensus
ONT-TCRconsensus creates and counts high accuracy full-length unique TCR molecule consensus sequences. 


<!-- Installation -->

## Installation
### Conda from yml (currently only supported option)
1. Clone the repo:
   ```sh
   git clone https://github.com/schumacherlab/nanopore_tcr_consensus.git
   ```
2. Navigate to the project directory:
   ```sh
   cd ONT-TCRconsensus
   ```
3. Create conda environment:
   ```
   conda env create -f ont_tcr_consensus.yml
   conda activate ont_tcr_consensus_env
   ```
4. Install <a href="https://github.com/nanoporetech/dorado" target="_blank">dorado</a>


<!-- ~Performance -->
## Performance 

ONT-TCRconensus performance on a ~70 M reads Promethion R10.4.1 run:
| CPU model | # CPUs | Memory G | ~Run time h | 
|----------------|------------------|------------------|----------------|
| Intel Xeon Silver | 110    |  275G | 20-24 h |
| Intel Xeon Gold  | 128    | 800G  | 5-6 h| 

We refer to <a href="https://github.com/nanoporetech/dorado" target="_blank">dorado</a> for basecalling performance. 


<!-- GETTING STARTED -->
## Running

0. Run basecalling using dorado:
   ```
   sbatch run_basecall_pipeline_multi-gpu.sh
   ```

1. Generate a reference.fa file with the <a href="https://github.com/schumacherlab/TCRtoolbox" target="_blank">TCRtoolbox</a> package `generate_assembly_nanopore_nt_refs` function. 

2. Update run_config.json by inserting paths and run ONT-TCRconsensus:  
 
   ```
   # SLURM cluster (conda): 
   sbatch run_tcr_consensus_slurm.sh run_config.json

   # Any cluster (conda): 
   mkdir -p ./logs
   source ~/miniconda3/etc/profile.d/conda.sh
   conda activate ont_tcr_consensus
   tcr_consensus run_config.json
   ```
   
3. After ONT-TCRconsensus pipeline is done running, for each barcode run QC plots and a filtered UMI count .csv can be generated in a common `outs` directory generated in the nanopore run      directory with `analysis.ipynb`:
   1. Open `analysis.ipynb` and change `nanopore_project_dir =` to path to nanopore run directory
   2. Add a libraries.csv:
      ```
      # Without ref_library_name
      barcode,library_name,ref_library_name,log_umi_counts_filter_threshold
      barcode02,baseline_b,,1
      barcode05,ylq_cd69_pos_b,,1
      barcode08,glc_cd69_pos_b,,1
      barcode11,glc_cd69_neg_b,,1
      
      # With ref_library_name (TCR library identifier in reference.fa fasta header names when multiple libraries are multiplexed in a single ONT run)
      barcode,library_name,ref_library_name,log_umi_counts_filter_threshold
      barcode14,NSCLC57,NSCLC57,2.5
      barcode15,N03LAM397,N03LAM397,1.5
      barcode16,YWE,ywe,1.5
      barcode17,blos_c,blos_c,1
      barcode18,blos_p,blos_p,1
      barcode19,mlm_m,mlm_m,2
      barcode20,mlm_d,mlm_d,1
      barcode21,mlm_t,mlm_t,2
      barcode22,str_b,str_b,1
      ```
   3. Adjust the following plotting parameters if needed: 
      ```blast_id_ylim_max = 350 
         umi_count_xlim_max = 900 
         umi_count_bin_size = 10 
         umi_count_ylim_max = 0.1
         log_umi_count_xlim_max = 8
         log_umi_count_bin_size = 0.25
         log_umi_count_ylim_max = 1.0
         most_similar_blast_id_threshold = 0.99925 # for zooming 
         ```
    4. (Optional) add a custom reference name set to additionaly color in plots:
      ```
      gil_plate_1_duplicate_1_ref_set = set()
      pattern = r"^[A-H](?:[1-9]|1[0-9]|2[0-4])$"
      for ref in non_n_ref_names_set:
          if (ref.endswith("GILGFVFTL") and ref.startswith("1_1") and re.match(pattern, ref.split("_")[2])):
              gil_plate_1_duplicate_1_ref_set.add(ref)

      custom_tcr_set_dict = {"gil_plate_1_tcrs": gil_plate_1_duplicate_1_ref_set}
      ```
   6. Run QC plotting and count generation code. Results are written to `outs`.
  

<!-- Citing This Work -->
## Citing This Work


<!-- LICENSE -->
## License

Distributed under the Apache 2.0 License.


