import collections
import glob
import json
import os
import re
from typing import Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from matplotlib.cm import ScalarMappable, viridis
from matplotlib.colors import Normalize

viridis.set_bad("lightgrey")


def set_mpl_params(mpl):
    mylines = 0.75  # pt updated
    mpl.rcParams["axes.linewidth"] = mylines  # default 1
    mpl.rcParams["ytick.direction"] = "out"  # default 'in'
    mpl.rcParams["xtick.direction"] = "out"  # default 'in'
    mpl.rcParams["xtick.major.size"] = 3.5  # default 4
    mpl.rcParams["ytick.major.size"] = 3.5  # default 4
    mpl.rcParams["xtick.major.width"] = mylines  # default 0.5
    mpl.rcParams["ytick.major.width"] = mylines  # default 0.5
    mpl.rcParams["grid.linewidth"] = mylines / 1.5  # default 0.5
    mpl.rcParams["grid.color"] = "0.8"  # default 'k'
    mpl.rcParams["grid.linestyle"] = "solid"  # default ':'
    mpl.rcParams["legend.frameon"] = False  # default True
    mpl.rcParams["font.family"] = "Arial"  # added
    mpl.rcParams["figure.dpi"] = 300
    mpl.rc("savefig", dpi=300)
    mpl.rcParams["pdf.fonttype"] = 42  # embed fonts as editable text in illustrator
    mpl.rcParams["ps.fonttype"] = 42  # embed fonts as editable text in illustrator

    return mpl, mylines


def startfig(w=4, h=2, rows=1, columns=1, wrs=None, hrs=None, frameon=True, return_first_ax=True):
    ratio = 0.393701  # 1 cm in inch
    myfigsize = (w * ratio, h * ratio)
    fig = plt.figure(figsize=(myfigsize))
    gs = mpl.gridspec.GridSpec(rows, columns, width_ratios=wrs, height_ratios=hrs)
    if return_first_ax:
        a = fig.add_subplot(gs[0, 0], frameon=frameon)
        return a, fig, gs
    else:
        return fig, gs


def split_ref_names_set_into_subsets(ref_names_set: set):
    non_n_ref_names_set = set()
    full_n_ref_names_set = set()
    cdr3j_n_ref_names_set = set()
    v_n_ref_names_set = set()
    for ref in ref_names_set:
        if ref.endswith("full_n"):
            full_n_ref_names_set.add(ref)
        elif ref.endswith("cdr3j_n"):
            cdr3j_n_ref_names_set.add(ref)
        elif ref.endswith("v_n"):
            v_n_ref_names_set.add(ref)
        else:
            non_n_ref_names_set.add(ref)

    return (
        non_n_ref_names_set,
        full_n_ref_names_set,
        cdr3j_n_ref_names_set,
        v_n_ref_names_set,
    )


def parse_raw_nanopore_qual_from_fastq_stats(
    fastq_stats_before_filtering_txt: Union[str, os.PathLike[str]],
):
    stats_df = pd.read_csv(fastq_stats_before_filtering_txt, sep="\s+")
    raw_nanopore_qual = 1 - 10 ** (-stats_df.loc[:, "AvgQual"].values[0] / 10)
    return raw_nanopore_qual


def parse_merged_consensus_bam_filter_log(
    merged_consensus_bam_filter_log: Union[str, os.PathLike[str]],
):
    with open(merged_consensus_bam_filter_log, "r") as log_in:
        lines = log_in.readlines()
        minimal_blast_id = float(lines[2].split("- minimal blast identity with reference:")[1])
        max_too_many_nt = float(lines[6].split("- Median allowed too few nucleotides/region:")[1])
        max_too_few_nt = float(lines[7].split("- Median allowed too many nucleotides/region:")[1])
        percent_too_few_nt = float(
            lines[11].split("% of primary alignments that have shorter overlap than minimal region overlap:")[1]
        )
        percent_too_many_nt = float(lines[12].split("% of primary alignments that have too long reads:")[1])
        percent_correct_overlap_length = 100 - (percent_too_few_nt + percent_too_many_nt)

    return (
        minimal_blast_id,
        max_too_many_nt,
        max_too_few_nt,
        percent_too_few_nt,
        percent_too_many_nt,
        percent_correct_overlap_length,
    )


def parse_quantile_95_blast_id_from_self_homology_log(
    self_homology_log: Union[str, os.PathLike[str]],
):
    with open(self_homology_log, "r") as log_in:
        lines = log_in.readlines()
        quantile_95_blast_id = float(lines[5].split("0.950 quantile blast identity of most similar regions:")[1])
    return quantile_95_blast_id


def plot_blast_id_with_most_similar_tcr(
    nanopore_project_dir: Union[str, os.PathLike[str]],
    xlim_min: float,
    ylim_max: int,
    library_ref_id: str,
    minimal_blast_id: float = None,
    quantile_95_blast_id: float = None,
    raw_nanopore_qual: float = None,
    title: str = None,
    save_suffix: str = None,
    outs_dir: Union[str, os.PathLike[str]] = None,
):
    most_similar_region_blob = glob.glob(
        os.path.join(nanopore_project_dir, "*_ref_homology_out_most_similar_region_dict.json")
    )
    if not len(most_similar_region_blob) == 1:
        raise Exception(
            "There are multiple ref_homology_out_most_similar_region_dict.json! This should not be possible!"
        )
    with open(most_similar_region_blob[0], "r") as json_in:
        most_similar_region_dict = json.load(json_in)

    if library_ref_id:
        most_similar_region_dict_ref = {
            ref: blast_id_list for ref, blast_id_list in most_similar_region_dict.items() if library_ref_id in ref
        }
    else:
        most_similar_region_dict_ref = most_similar_region_dict

    if not minimal_blast_id:
        ax, fig, gs = startfig(8, 7)
    else:
        ax, fig, gs = startfig(13, 7)
    max_blast_id_list = [np.max(blast_id_list) for blast_id_list in most_similar_region_dict.values()]
    if library_ref_id:
        max_blast_id_list_ref = [np.max(blast_id_list) for blast_id_list in most_similar_region_dict_ref.values()]
        median_blast_id_ref = np.median(max_blast_id_list_ref)
    ax.hist(
        max_blast_id_list,
        color="black",
        alpha=0.25,
        edgecolor="none",
        bins=np.arange(0, 1.0 + 0.00015, 0.00015),
        density=True,
    )

    if library_ref_id:
        ax.hist(
            max_blast_id_list_ref,
            color="red",
            alpha=0.25,
            edgecolor="none",
            bins=np.arange(0, 1.0 + 0.00015, 0.00015),
            density=True,
        )
        ax.axvline(
            median_blast_id_ref,
            color="darkred",
            linewidth=0.5,
            label="Median BLAST identity library\nwith closest TCR = " + str(round(median_blast_id_ref, 6)),
        )

    if minimal_blast_id:
        ax.axvline(
            minimal_blast_id,
            color="red",
            linewidth=0.5,
            label="Minimal BLAST identity\nto distinguish all TCRs\n= " + str(round(minimal_blast_id, 6)),
        )
        # ax.axvline(quantile_95_blast_id, color = 'blue', linewidth = 0.5, label = 'Required minimal blast identity\nto distinguish 95% of all TCRs')
        ax.axvline(
            raw_nanopore_qual,
            color="orange",
            linewidth=0.5,
            label="Average raw nanopore\nsequencing quality",
        )
        # ax.axvline(0.99925, color = 'black', linewidth = 0.5)

    ax.set_xlim(xlim_min, 1.0)
    ax.axvline(
        np.median(max_blast_id_list),
        color="cornflowerblue",
        linewidth=0.5,
        label="Median BLAST identity with\nclosest TCR = " + str(round(np.median(max_blast_id_list), 6)),
    )
    ax.set_ylim(0, ylim_max)
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels([round(x, 4) for x in ax.get_xticks()], fontsize=7, rotation=90)
    ax.tick_params(axis="y", labelsize=7)
    ax.tick_params(top=False, right=False)
    ax.set_xlabel("BLAST identity with closest TCR", fontsize=7)
    ax.set_ylabel("Frequency", fontsize=7)
    ax.set_title(title, fontsize=8)
    if minimal_blast_id:
        ax.legend(fontsize=7, loc="center left", bbox_to_anchor=(1, 0.5))
        # ax.legend(fontsize=6, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=len(ax.get_legend_handles_labels()[1]))
    fig.tight_layout()

    if outs_dir:
        if save_suffix:
            fig.savefig(
                os.path.join(
                    outs_dir,
                    "blast_id_with_most_similar_tcr_hist_" + save_suffix + ".pdf",
                )
            )
            plt.close()
        else:
            fig.savefig(os.path.join(outs_dir, "blast_id_with_most_similar_tcr_hist.pdf"))
            plt.close()
    else:
        plt.show()


def plot_hist_nt_too_short(
    nt_too_short_csv: Union[str, os.PathLike[str]],
    tcr_refs_df: pd.DataFrame,
    max_too_few_nt: int,
    percent_too_few_nt: float,
    title: str = None,
    outs_dir: Union[str, os.PathLike[str]] = None,
):
    nt_too_short_df = pd.read_csv(nt_too_short_csv)

    ax, fig, gs = startfig(14, 7)

    ax.hist(
        nt_too_short_df.loc[:, "number_of_nt"],
        bins=np.arange(0, 100, 2),
        color="black",
        alpha=0.25,
        edgecolor="none",
    )
    ax.axvline(max_too_few_nt, color="red", linewidth=0.75, label="Max too few nt allowed")
    ax.axvline(
        tcr_refs_df.loc[:, "cdr3j_beta_nt_order_primers"].str.len().quantile(0.25) - 80 - 16 - 8,
        color="lightblue",
        linewidth=0.75,
        label="25th quantile CDR3b-Jb length",
    )
    ax.axvline(
        tcr_refs_df.loc[:, "cdr3j_beta_nt_order_primers"].str.len().quantile(0.50) - 80 - 16 - 8,
        color="blue",
        linewidth=0.75,
        label="Median CDR3b-Jb length",
    )
    ax.axvline(
        tcr_refs_df.loc[:, "cdr3j_beta_nt_order_primers"].str.len().quantile(0.75) - 80 - 16 - 8,
        color="darkblue",
        linewidth=0.75,
        label="75th quantile CDR3b-Jb length",
    )
    ax.axvline(
        tcr_refs_df.loc[:, "cdr3j_alpha_nt_order_primers"].str.len().quantile(0.50) - 80 - 16 - 8,
        color="orange",
        linewidth=0.75,
        label="Median length CDR3a-Ja length",
    )
    ax.tick_params(axis="both", labelsize=8)
    ax.tick_params(top=False, right=False)
    ax.set_xlabel("# nucleotides too short", fontsize=8)
    ax.set_ylabel("Frequency", fontsize=8)
    ax.set_title(
        title + "\n" + str(percent_too_few_nt) + "% of all TCR alignments\nhave less nt overlap than threshold",
        fontsize=8,
    )
    ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    if outs_dir:
        fig.savefig(os.path.join(outs_dir, "nucleotides_too_short_hist.pdf"))
        plt.close()
    else:
        plt.show()


def plot_hist_nt_too_long(
    nt_too_long_csv: Union[str, os.PathLike[str]],
    max_too_many_nt: int,
    percent_too_many_nt: float,
    title: str = None,
    outs_dir: Union[str, os.PathLike[str]] = None,
):
    nt_too_long_df = pd.read_csv(nt_too_long_csv)

    ax, fig, gs = startfig(11, 7)

    ax.hist(
        nt_too_long_df.loc[:, "number_of_nt"],
        bins=np.arange(0, 60, 2),
        color="black",
        alpha=0.25,
        edgecolor="none",
    )
    ax.tick_params(axis="both", labelsize=8)
    ax.tick_params(top=False, right=False)
    ax.axvline(max_too_many_nt, color="red", linewidth=0.75, label="Max too many nt\nallowed")
    ax.set_xlabel("# nucleotides too long", fontsize=8)
    ax.set_ylabel("Frequency", fontsize=8)
    ax.set_title(
        title + "\n" + str(percent_too_many_nt) + "% of all TCR alignments\nhave more nt overlap than threshold",
        fontsize=8,
    )
    ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    if outs_dir:
        fig.savefig(os.path.join(outs_dir, "nucleotides_too_long_hist.pdf"))
        plt.close()
    else:
        plt.show()


def plot_hist_percent_alignments_above_blast_id(
    blast_id_csv: Union[str, os.PathLike[str]],
    minimal_blast_id: float,
    quantile_95_blast_id: float,
    percent_correct_overlap_length: float,
    title: str = None,
    outs_dir: Union[str, os.PathLike[str]] = None,
):
    blast_id_df = pd.read_csv(blast_id_csv)

    ax, fig, gs = startfig(12, 7)

    blast_ids = blast_id_df.loc[:, "blast_id"]
    bins = np.arange(0.995, 1.0002, 0.0001)
    hist, bin_edges = np.histogram(blast_ids, bins=bins)

    total_values = len(blast_ids)
    hist_percentage = (hist / total_values) * 100  # Convert to percentage

    ax.bar(
        bin_edges[:-1] - (0.0001 / 2),
        hist_percentage,
        width=np.diff(bin_edges),
        color="black",
        alpha=0.25,
        edgecolor="none",
        align="edge",
    )

    ax.axvline(
        minimal_blast_id,
        color="red",
        linewidth=0.75,
        label="Required minimal blast identity\nto distinguish all TCRs",
    )
    ax.axvline(
        quantile_95_blast_id,
        color="blue",
        linewidth=0.75,
        label="Required minimal blast identity\nto distinguish 95% of all TCRs",
    )

    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels([round(x, 5) for x in ax.get_xticks()], fontsize=8, rotation=90)
    ax.tick_params(axis="y", labelsize=8)
    ax.tick_params(top=False, right=False)
    ax.set_xlabel("Blast identity with reference", fontsize=8)
    ax.set_ylabel("% of all TCR alignments\nwith correct overlap length", fontsize=8)
    ax.set_xlim(0.995, 1.001)
    ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_title(
        title
        + "\n"
        + str(round(percent_correct_overlap_length, 2))
        + "% of all TCR alignments\nhave correct overlap length",
        fontsize=8,
    )
    fig.tight_layout()
    if outs_dir:
        fig.savefig(os.path.join(outs_dir, "precision_blast_id_hist.pdf"))
        plt.close()
    else:
        plt.show()


def plot_hist_subreads_per_umi_cluster(
    number_subreads_blast_id_csv: Union[str, os.PathLike[str]],
    title: str = None,
    outs_dir: Union[str, os.PathLike[str]] = None,
):
    num_subreads_blast_id_df = pd.read_csv(number_subreads_blast_id_csv)

    if num_subreads_blast_id_df.loc[:, "number_of_subreads"].max() > 60:
        print(
            "Warning: max detected # of subreads/UMI cluster is:",
            num_subreads_blast_id_df.loc[:, "number_of_subreads"].max(),
            "! Did you allow more than 60 in run_config.json?",
        )
    ax, fig, gs = startfig(13, 7)
    counts, bins = np.histogram(
        num_subreads_blast_id_df.loc[:, "number_of_subreads"],
        bins=np.arange(0.5, 120.5, 1),
        density=False,
    )

    total_count = counts.sum()
    percentages = (counts / total_count) * 100
    bin_centers = (bins[:-1] + bins[1:]) / 2
    ax.bar(bin_centers, percentages, width=np.diff(bins), color="grey", edgecolor="none")

    ax.set_xlim(0.5, 121.5)
    ax.set_xticks(np.arange(1, 120, 1))
    ax.set_xticklabels(np.arange(1, 120, 1), fontsize=8, rotation=90)
    ax.set_ylabel("% UMI clusters", fontsize=8)
    ax.set_xlabel("# subreads in UMI cluster", fontsize=8)

    ax.set_title(title, fontsize=8)
    ax.tick_params("y", labelsize=8)
    ax.tick_params(top=False, right=False)
    fig.tight_layout()

    if outs_dir:
        fig.savefig(os.path.join(outs_dir, "number_subreads_umi_cluster_hist.pdf"))
        plt.close()
    else:
        plt.show()


def plot_box_blast_id_vs_subreads_umi_cluster(
    number_subreads_blast_id_csv: Union[str, os.PathLike[str]],
    minimal_blast_id: float,
    title: str = None,
    outs_dir: Union[str, os.PathLike[str]] = None,
):
    ax, fig, gs = startfig(20, 7)

    num_subreads_blast_id_df = pd.read_csv(number_subreads_blast_id_csv)
    vc = num_subreads_blast_id_df.loc[:, "number_of_subreads"].value_counts()
    num_subreads_blast_id_df = num_subreads_blast_id_df.loc[
        num_subreads_blast_id_df.loc[:, "number_of_subreads"].isin(vc.index), :
    ]
    unique_num_subreads = sorted(num_subreads_blast_id_df.loc[:, "number_of_subreads"].unique())

    boxplot_data = [
        list(
            num_subreads_blast_id_df.loc[
                num_subreads_blast_id_df["number_of_subreads"] == num_subread,
                "blast_id",
            ]
        )
        for num_subread in unique_num_subreads
    ]

    parts = ax.boxplot(
        boxplot_data,
        tick_labels=unique_num_subreads,
        whis=2.5,
        showmeans=True,
        meanline=True,
        meanprops={"color": "orange", "linestyle": "-", "linewidth": 1.0},
    )
    ax.set_xlabel("# subreads in UMI cluster", fontsize=8)
    ax.set_ylabel("Blast identity with reference", fontsize=8)
    ax.set_ylim(0.93, 1.0001)
    ax.axhline(minimal_blast_id, color="red", linewidth=0.5)

    linewidth = 0.5
    for i, box in enumerate(parts["boxes"]):
        box.set(color="black", linewidth=linewidth)
        parts["whiskers"][2 * i].set(color="black", linewidth=linewidth, linestyle="-")  # Lower whisker
        parts["whiskers"][2 * i + 1].set(color="black", linewidth=linewidth, linestyle="-")  # Upper whisker
        parts["caps"][2 * i].set(color="black", linewidth=linewidth)  # Lower cap
        parts["caps"][2 * i + 1].set(color="black", linewidth=linewidth)  # Upper cap
        parts["medians"][i].set(color="blue", linewidth=1.0)  # Median line
        parts["fliers"][i].set(markerfacecolor="black", marker="o", markersize=1)

    ax.tick_params(right=False)
    ax.tick_params("x", labelsize=8, rotation=90)
    ax.set_yticks(ax.get_yticks())
    ax.set_yticklabels([round(y, 5) for y in ax.get_yticks()], fontsize=8)

    ax.xaxis.set_ticks_position("both")
    secax = ax.secondary_xaxis("top")
    secax.set_xticks(ax.get_xticks())
    secax.set_xticklabels(vc.values, fontsize=8, rotation=90)
    secax.set_xlabel(title + "\n# UMI clusters", fontsize=8)
    fig.tight_layout()
    if outs_dir:
        fig.savefig(os.path.join(outs_dir, "number_subreads_umi_cluster_vs_blast_id_boxplot.pdf"))
        plt.close()
    else:
        plt.show()

    ax, fig, gs = startfig(30, 7)

    boxplot_data = [
        list(
            num_subreads_blast_id_df.loc[
                num_subreads_blast_id_df["number_of_subreads"] == num_subread,
                "blast_id",
            ]
        )
        for num_subread in unique_num_subreads
    ]

    parts = ax.boxplot(
        boxplot_data,
        tick_labels=unique_num_subreads,
        whis=2.5,
        showmeans=True,
        meanline=True,
        meanprops={"color": "orange", "linestyle": "-", "linewidth": 1.0},
    )
    ax.set_xlabel("# subreads in UMI cluster", fontsize=8)
    ax.set_ylabel("Blast identity with reference", fontsize=8)
    ax.set_ylim(0.995, 1.0001)
    ax.axhline(minimal_blast_id, color="red", linewidth=0.5)

    linewidth = 0.5
    for i, box in enumerate(parts["boxes"]):
        box.set(color="black", linewidth=linewidth)
        parts["whiskers"][2 * i].set(color="black", linewidth=linewidth, linestyle="-")  # Lower whisker
        parts["whiskers"][2 * i + 1].set(color="black", linewidth=linewidth, linestyle="-")  # Upper whisker
        parts["caps"][2 * i].set(color="black", linewidth=linewidth)  # Lower cap
        parts["caps"][2 * i + 1].set(color="black", linewidth=linewidth)  # Upper cap
        parts["medians"][i].set(color="blue", linewidth=1.0)  # Median line
        parts["fliers"][i].set(markerfacecolor="black", marker="o", markersize=1)

    ax.tick_params(right=False)
    ax.tick_params("x", labelsize=8, rotation=90)
    ax.set_yticks(ax.get_yticks())
    ax.set_yticklabels([round(y, 5) for y in ax.get_yticks()], fontsize=8)

    ax.xaxis.set_ticks_position("both")
    secax = ax.secondary_xaxis("top")
    secax.set_xticks(ax.get_xticks())
    secax.set_xticklabels(vc.values, fontsize=8, rotation=90)
    secax.set_xlabel(title + "\n# UMI clusters", fontsize=8)
    fig.tight_layout()
    if outs_dir:
        fig.savefig(
            os.path.join(
                outs_dir,
                "number_subreads_umi_cluster_vs_blast_id_boxplot_yaxis_zoom_in.pdf",
            )
        )
        plt.close()
    else:
        plt.show()


def log_normalize_umi_counts(counts_df: pd.DataFrame):
    counts_df.loc[:, "log_transformed_count"] = np.log(counts_df.loc[:, "Count"])
    return counts_df


def filter_counts_on_umi_quantile_threshold(counts_df: pd.DataFrame, quantile_umi_threshold: float = 0.05):
    counts_df = counts_df.loc[
        counts_df.loc[:, "Count"] > np.quantile(counts_df.loc[:, "Count"], quantile_umi_threshold),
        :,
    ].copy()
    return counts_df


def filter_counts_on_log_umi_count_threshold(counts_df: pd.DataFrame, log_umi_counts_filter_threshold: float):
    return counts_df.loc[counts_df.loc[:, "log_transformed_count"] > log_umi_counts_filter_threshold, :].copy()


def plot_umi_counts_hist(
    counts_df: pd.DataFrame,
    xlim_max: int,
    bin_size: int,
    ylim_max: float,
    title: str = None,
    plot_nb_dist_fit: bool = True,
    plot_percentiles: bool = True,
    save_suffix: str = None,
    outs_dir: Union[str, os.PathLike[str]] = None,
):
    counts_df = counts_df.copy()

    if xlim_max < counts_df["Count"].quantile(0.99):
        print(
            "Warning: xlim_max is smaller than the 99th percentile UMI count:",
            counts_df["Count"].quantile(0.99),
        )

    ax, fig, gs = startfig(12, 7)
    ax.hist(
        counts_df.loc[:, "Count"] - 1,
        bins=np.arange(0, xlim_max, bin_size),
        alpha=0.25,
        edgecolor="none",
        density=True,
        color="black",
        zorder=4,
    )

    if plot_percentiles:
        ax.axvline(np.median(counts_df["Count"] - 1), color="blue", zorder=6, label="median")
        ax.axvline(
            np.quantile(counts_df["Count"] - 1, 0.05),
            color="yellow",
            zorder=6,
            label="5th percentile",
        )
        ax.axvline(
            np.quantile(counts_df["Count"] - 1, 0.95),
            color="black",
            zorder=6,
            label="95th percentile",
        )
    ax.set_xlim(0, xlim_max)
    ax.set_ylim(0, ylim_max)
    ax.set_xticks(ax.get_xticks())
    ax.set_yticks(ax.get_yticks())
    ax.set_xticklabels([round(x, 2) for x in ax.get_xticks()], fontsize=8, rotation=90)
    ax.set_yticklabels([round(y, 2) for y in ax.get_yticks()], fontsize=8)
    ax.set_xlabel("UMI counts", fontsize=8)
    ax.set_ylabel("Density", fontsize=8)
    ax.tick_params(top=False, right=False)

    ax.set_title(title, fontsize=7)

    if plot_nb_dist_fit:
        mean, var = np.mean(counts_df["Count"]), np.var(counts_df["Count"])
        if var > mean:
            r = (mean**2) / (var - mean)
            p = mean / var
        else:
            raise ValueError("Variance must be greater than mean for negative binomial parameters.")
        ks_statistic, p_value = sp.stats.kstest(counts_df["Count"] - 1, "nbinom", args=(r, p))
        print(f"KS Statistic: {ks_statistic}, p-value: {p_value}")
        x = np.arange(0, np.max(counts_df["Count"]) + 1)
        pmf = sp.stats.nbinom.pmf(x - 1, r, p)
        ax.plot(x, pmf, "r-", label="Fitted Negative Binomial")

    ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()

    if outs_dir:
        if save_suffix:
            fig.savefig(os.path.join(outs_dir, "umi_counts_hist_" + save_suffix + ".pdf"))
            plt.close()
        else:
            fig.savefig(os.path.join(outs_dir, "umi_counts_hist.pdf"))
            plt.close()
    else:
        plt.show()


def plot_log_transformed_umi_counts_hist(
    counts_df: pd.DataFrame,
    xlim_max: float,
    bin_size: float,
    ylim_max: float,
    custom_tcr_set: set = None,
    label_custom_tcr_set: str = None,
    density: bool = True,
    plot_normal_dist_fit: bool = True,
    plot_percentiles: bool = True,
    log_umi_counts_filter_threshold: float = None,
    most_similar_blast_id_threshold: float = 0.99925,
    title: str = None,
    save_suffix: str = None,
    nanopore_project_dir: Union[str, os.PathLike[str]] = None,
    outs_dir: Union[str, os.PathLike[str]] = None,
):
    counts_df = counts_df.copy()

    if xlim_max < counts_df["log_transformed_count"].quantile(0.99):
        print(
            "Warning: xlim_max is smaller than the 99th percentile UMI count:",
            counts_df["log_transformed_count"].quantile(0.99),
        )

    ax, fig, gs = startfig(8, 5)
    ax.hist(
        counts_df.loc[:, "log_transformed_count"],
        bins=np.arange(0, xlim_max, bin_size),
        density=density,
        alpha=0.25,
        edgecolor="none",
        color="black",
        zorder=4,
        label="All TCRs",
    )

    if most_similar_blast_id_threshold:
        most_similar_region_blob = glob.glob(
            os.path.join(nanopore_project_dir, "*_ref_homology_out_most_similar_region_dict.json")
        )
        if not len(most_similar_region_blob) == 1:
            raise Exception(
                "There are multiple ref_homology_out_most_similar_region_dict.json! This should not be possible!"
            )
        with open(most_similar_region_blob[0], "r") as json_in:
            most_similar_region_dict = json.load(json_in)
        most_similar_counts_df = counts_df.loc[
            counts_df.loc[:, "TCR"].isin(
                [
                    region
                    for region, blast_id_list in most_similar_region_dict.items()
                    if np.max(blast_id_list) > most_similar_blast_id_threshold
                ]
            ),
            :,
        ]
        ax.hist(
            most_similar_counts_df.loc[:, "log_transformed_count"],
            density=density,
            alpha=0.25,
            edgecolor="none",
            bins=np.arange(0, xlim_max, bin_size),
            color="red",
            zorder=4,
            label="Most similar TCRs",
        )
    if custom_tcr_set:
        ax.hist(
            counts_df.loc[counts_df.loc[:, "TCR"].isin(custom_tcr_set), "log_transformed_count"],
            density=density,
            alpha=0.25,
            edgecolor="none",
            bins=np.arange(0, xlim_max, bin_size),
            color="blue",
            zorder=4,
            label=label_custom_tcr_set,
        )

    if log_umi_counts_filter_threshold:
        ax.axvline(
            log_umi_counts_filter_threshold,
            color="orange",
            zorder=6,
            label="Filter threshold",
        )

    if plot_percentiles:
        ax.axvline(
            np.quantile(counts_df["log_transformed_count"], 0.05),
            color="yellow",
            zorder=6,
            label="5th percentile",
        )
        ax.axvline(
            np.median(counts_df["log_transformed_count"]),
            color="blue",
            zorder=6,
            label="median",
        )
        ax.axvline(
            np.quantile(counts_df["log_transformed_count"], 0.95),
            color="black",
            zorder=6,
            label="95th percentile",
        )

    ax.set_xlim(0, xlim_max)
    ax.set_ylim(0, ylim_max)
    ax.set_xticks(ax.get_xticks())
    ax.set_yticks(ax.get_yticks())
    ax.set_xticklabels([round(x, 2) for x in ax.get_xticks()], fontsize=7, rotation=90)
    ax.set_yticklabels([round(y, 2) for y in ax.get_yticks()], fontsize=7)
    ax.set_xlabel("log(TCR UMI counts)", fontsize=7)
    ax.set_ylabel("Density", fontsize=7)
    ax.tick_params(top=False, right=False)

    log10_fold_diff_5th_95th_percentile = round(
        np.log10(counts_df.loc[:, "Count"].quantile(0.95)) - np.log10(counts_df.loc[:, "Count"].quantile(0.05)),
        2,
    )

    title = title + "\nlog10 diff. 95th vs 5th percentile = " + str(log10_fold_diff_5th_95th_percentile)
    ax.set_title(title, fontsize=7)

    if plot_normal_dist_fit:
        mean, std = (
            np.mean(counts_df["log_transformed_count"]),
            np.std(counts_df["log_transformed_count"]),
        )
        ks_statistic, p_value = sp.stats.kstest(counts_df["log_transformed_count"], "norm", args=(mean, std))
        print(f"KS Statistic: {ks_statistic}, p-value: {p_value}")
        x = np.linspace(
            np.min(counts_df["log_transformed_count"]),
            np.max(counts_df["log_transformed_count"]),
            100,
        )
        pdf = sp.stats.norm.pdf(x, mean, std)
        ax.plot(x, pdf, "r-", label="Fitted\nNormal Distribution")

    ax.legend(fontsize=7, loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()

    if outs_dir:
        if save_suffix:
            fig.savefig(os.path.join(outs_dir, "log_transformed_umi_counts_hist_" + save_suffix + ".pdf"))
            plt.close()
        else:
            fig.savefig(os.path.join(outs_dir, "log_transformed_umi_counts_hist.pdf"))
            plt.close()
    else:
        plt.show()


def write_results_summary_txt(
    counts_df: pd.DataFrame,
    merged_consensus_bam_filter_log: Union[str, os.PathLike[str]],
    region_name_col: str,
    library_name: str,
    non_n_ref_names_set: set,
    outs_dir: Union[str, os.PathLike[str]] = None,
):
    counts_df = counts_df.copy()
    detected_regions_set = set(counts_df.loc[:, region_name_col])
    max_plate_num = max([int(ref.split("_")[1]) for ref in non_n_ref_names_set])
    not_detected_tcr_set = non_n_ref_names_set.difference(detected_regions_set)
    num_missing_tcr_per_plate_series = pd.Series(
        [int(region.split("_")[1]) for region in not_detected_tcr_set]
    ).value_counts()

    log10_fold_diff_5th_95th_percentile = round(
        np.log10(counts_df.loc[:, "Count"].quantile(0.95)) - np.log10(counts_df.loc[:, "Count"].quantile(0.05)),
        2,
    )

    if outs_dir:
        with open(os.path.join(outs_dir, "results_summary.txt"), "w") as txt_out:
            print("Library:", library_name, file=txt_out)
            with open(merged_consensus_bam_filter_log) as log_in:
                lines = log_in.readlines()
                for line in lines[:24]:
                    print(line.strip("\n"), file=txt_out)

            print(
                "# total number " + region_name_col + " in reference:",
                len(non_n_ref_names_set.union(detected_regions_set)),
                file=txt_out,
            )
            print(
                "# detected " + region_name_col + ":",
                len(non_n_ref_names_set.intersection(detected_regions_set)),
                file=txt_out,
            )
            print(
                "# not detected " + region_name_col + ":",
                len(not_detected_tcr_set),
                file=txt_out,
            )
            print(
                "Sensitivity:",
                round(
                    len(non_n_ref_names_set.intersection(detected_regions_set)) / len(non_n_ref_names_set),
                    4,
                ),
                "\n",
                file=txt_out,
            )
            print("not detected:\n" + str(not_detected_tcr_set), "\n", file=txt_out)
            print("# of TCRs not detected per plate:", file=txt_out)
            for plate in np.arange(1, max_plate_num + 1):
                if plate in num_missing_tcr_per_plate_series.index:
                    print(
                        str(plate) + ": " + str(num_missing_tcr_per_plate_series.loc[plate]),
                        file=txt_out,
                    )
                else:
                    print(str(plate) + ": " + str(0), file=txt_out)
            print(
                "\nlog10 fold difference between 5th and 95th percentile:",
                log10_fold_diff_5th_95th_percentile,
                file=txt_out,
            )

    else:
        print("Library:", library_name)
        with open(merged_consensus_bam_filter_log) as log_in:
            lines = log_in.readlines()
            for line in lines[:24]:
                print(line.strip("\n"))

        print(
            "# total number " + region_name_col + " in reference:",
            len(non_n_ref_names_set.union(detected_regions_set)),
        )
        print(
            "# detected " + region_name_col + ":",
            len(non_n_ref_names_set.intersection(detected_regions_set)),
        )
        print("# not detected " + region_name_col + ":", len(not_detected_tcr_set))
        print(
            "Sensitivity:",
            round(
                len(non_n_ref_names_set.intersection(detected_regions_set)) / len(non_n_ref_names_set),
                4,
            ),
            "\n",
        )
        print("not detected: " + str(not_detected_tcr_set), "\n")
        print("# of TCRs not detected per plate:")
        for plate in np.arange(1, max_plate_num + 1):
            print(str(plate) + ": " + str(num_missing_tcr_per_plate_series.loc[plate]))
        print("\nlog10 fold difference between 5th and 95th percentile:")


def plot_plate_umi_counts(
    counts_df: pd.DataFrame,
    non_n_ref_names_set: set,
    region_name_col: str,
    outs_dir: Union[str, os.PathLike[str]] = None,
):
    counts_df.copy()
    ref_names_wells_set = set([ref.split("_")[1] + "_" + ref.split("_")[2] for ref in non_n_ref_names_set])
    counts_df = counts_df.loc[counts_df.loc[:, region_name_col].isin(non_n_ref_names_set), :].copy()
    counts_df.loc[:, "detected_wells"] = [
        ref.split("_")[1] + "_" + ref.split("_")[2] for ref in counts_df.loc[:, region_name_col]
    ]
    detected_wells_set = set(counts_df.loc[:, "detected_wells"])
    min_plate_num = min([int(ref.split("_")[1]) for ref in non_n_ref_names_set])
    max_plate_num = max([int(ref.split("_")[1]) for ref in non_n_ref_names_set])
    plate_well_list_dict = collections.defaultdict(list)
    for ref in ref_names_wells_set:
        plate_well_list_dict[int(ref.split("_")[0])].append(ref.split("_")[1])

    norm = Normalize(
        vmin=counts_df.loc[:, "log_transformed_count"].min(),
        vmax=counts_df.loc[:, "log_transformed_count"].quantile(0.99),
    )
    mappable = ScalarMappable(norm=norm, cmap=viridis)

    plot_plate_df_dict = collections.defaultdict(pd.DataFrame)
    letter_to_number_map = dict(
        zip(
            [letter for letter in "ABCDEFGHIJKLMNOP"],
            [number for number in np.arange(len("ABCDEFGHIJKLMNOP"))],
        )
    )
    well_list = [x + str(y) for x in "ABCDEFGHIJKLMNOP" for y in range(1, 25)]

    for plate_name in np.arange(min_plate_num, max_plate_num + 1):
        print("Plotting plate:", plate_name)
        plot_plate_df_dict[plate_name] = pd.DataFrame(
            np.zeros((len(letter_to_number_map), 24)),
            index=letter_to_number_map.keys(),
            columns=np.arange(24),
        )
        for well in well_list:
            plate_well_name = str(plate_name) + "_" + well
            i = int(letter_to_number_map[re.split("\d", well)[0]])
            j = int(re.split("[A-P]", well)[1]) - 1
            if well not in plate_well_list_dict[plate_name]:
                plot_plate_df_dict[plate_name].iloc[i, j] = np.nan
                continue

            if plate_well_name in detected_wells_set:
                log_count = counts_df.loc[
                    counts_df.loc[:, "detected_wells"] == plate_well_name,
                    "log_transformed_count",
                ]
                plot_plate_df_dict[plate_name].iloc[i, j] = log_count.values[0]

        ax, fig, gs = startfig(w=36, h=55)
        ax.matshow(plot_plate_df_dict[plate_name], cmap=viridis, norm=norm)
        ax.set_xticks(
            np.arange(plot_plate_df_dict[plate_name].shape[1]),
            labels=[str(int(column) + 1) for column in np.arange(plot_plate_df_dict[plate_name].shape[1])],
            fontsize=8,
        )
        ax.set_yticks(
            np.arange(plot_plate_df_dict[plate_name].shape[0]),
            labels=letter_to_number_map.keys(),
            fontsize=8,
        )
        ax.xaxis.set_ticks_position("top")
        ax.yaxis.set_ticks_position("left")
        ax.set_title("Plate: " + str(plate_name), pad=25)
        cb = fig.colorbar(mappable, ax=ax, orientation="vertical", fraction=0.02, pad=0.025)
        cb.set_label("Log transformed\nUMI count", labelpad=20)
        fig.tight_layout()

        if outs_dir:
            fig.savefig(os.path.join(outs_dir, "plate_" + str(plate_name) + "_umi_counts.pdf"))
            plt.close()
        else:
            plt.show()


def plot_v_gene_fold_change_over_input_barplot(
    counts_df: pd.DataFrame,
    tcr_refs_df: pd.DataFrame,
    trav_col: str = "TRAV_IMGT_allele_collapsed",
    trbv_col: str = "TRBV_IMGT_allele_collapsed",
    title: str = None,
    outs_dir: Union[str, os.PathLike[str]] = None,
):
    counts_df = counts_df.copy()
    counts_df.loc[:, "TRAV"] = [
        tcr_refs_df.loc[tcr_refs_df["name"] == tcr, trav_col].values[0] for tcr in counts_df.loc[:, "TCR"]
    ]
    counts_df.loc[:, "TRBV"] = [
        tcr_refs_df.loc[tcr_refs_df["name"] == tcr, trbv_col].values[0] for tcr in counts_df.loc[:, "TCR"]
    ]
    counts_df.loc[:, "fraction_count"] = counts_df.loc[:, "Count"].div(counts_df.loc[:, "Count"].sum())

    trav_value_counts_df = pd.DataFrame(tcr_refs_df[trav_col].value_counts())
    trav_value_counts_df["normalized_count"] = trav_value_counts_df["count"] / trav_value_counts_df["count"].sum()
    trav_collapsed_df = (
        counts_df.groupby("TRAV").sum("fraction_count").sort_values("fraction_count", ascending=False).copy()
    )

    trav_fold_change_list = []
    for trav in trav_value_counts_df.index:
        if trav in trav_collapsed_df.index:
            trav_fold_change_list.append(
                trav_collapsed_df.loc[trav_collapsed_df.index == trav, "fraction_count"].values[0]
                / trav_value_counts_df.loc[trav_value_counts_df.index == trav, "normalized_count"].values[0]
            )
        else:
            trav_fold_change_list.append(
                0 / trav_value_counts_df.loc[trav_value_counts_df.index == trav, "normalized_count"].values[0]
            )

    trav_value_counts_df["fold_change_over_input"] = trav_fold_change_list
    trav_value_counts_df = trav_value_counts_df.sort_values("fold_change_over_input", ascending=False)

    ax, fig, gs = startfig(20, 7)

    ax.bar(
        np.arange(0, trav_value_counts_df.shape[0]),
        trav_value_counts_df["fold_change_over_input"] / trav_value_counts_df["fold_change_over_input"].median(),
        edgecolor="Black",
        linewidth=0.5,
        color="lightblue",
    )

    ax.axhline(1, color="red", linewidth=0.75)
    ax.set_xlim(-0.5, trav_value_counts_df.shape[0] + 0.5)
    ax.set_xticks(np.arange(0, trav_value_counts_df.shape[0]))
    ax.set_xticklabels(trav_value_counts_df.index, rotation=90, fontsize=7)
    ax.set_ylabel("Fold change over input\n(normalized to median)", fontsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.tick_params(right=False, top=False)
    ax.set_title(title, fontsize=8)
    fig.tight_layout()

    if outs_dir:
        fig.savefig(os.path.join(outs_dir, "TRAV_fold_change_over_input_barplot.pdf"))
        plt.close()
    else:
        plt.show()

    trbv_value_counts_df = pd.DataFrame(tcr_refs_df[trbv_col].value_counts())
    trbv_value_counts_df["normalized_count"] = trbv_value_counts_df["count"] / trbv_value_counts_df["count"].sum()
    trbv_collapsed_df = (
        counts_df.groupby("TRBV").sum("fraction_count").sort_values("fraction_count", ascending=False).copy()
    )

    trbv_fold_change_list = []
    for trbv in trbv_value_counts_df.index:
        if trbv in trbv_collapsed_df.index:
            trbv_fold_change_list.append(
                trbv_collapsed_df.loc[trbv_collapsed_df.index == trbv, "fraction_count"].values[0]
                / trbv_value_counts_df.loc[trbv_value_counts_df.index == trbv, "normalized_count"].values[0]
            )
        else:
            trbv_fold_change_list.append(
                0 / trbv_value_counts_df.loc[trbv_value_counts_df.index == trbv, "normalized_count"].values[0]
            )

    trbv_value_counts_df["fold_change_over_input"] = trbv_fold_change_list
    trbv_value_counts_df = trbv_value_counts_df.sort_values("fold_change_over_input", ascending=False)

    ax, fig, gs = startfig(20, 7)

    ax.bar(
        np.arange(0, trbv_value_counts_df.shape[0]),
        trbv_value_counts_df["fold_change_over_input"] / trbv_value_counts_df["fold_change_over_input"].median(),
        edgecolor="Black",
        linewidth=0.5,
        color="lightblue",
    )

    ax.axhline(1, color="red", linewidth=0.75)
    ax.set_xlim(-0.5, trbv_value_counts_df.shape[0] + 0.5)
    ax.set_xticks(np.arange(0, trbv_value_counts_df.shape[0]))
    ax.set_xticklabels(trbv_value_counts_df.index, rotation=90, fontsize=7)
    ax.set_ylabel("Fold change over input\n(normalized to median)", fontsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.tick_params(right=False, top=False)
    ax.set_title(title, fontsize=8)
    fig.tight_layout()
    if outs_dir:
        fig.savefig(os.path.join(outs_dir, "TRBV_fold_change_over_input_barplot.pdf"))
        plt.close()
    else:
        plt.show()


def plot_v_gene_fraction_missing_tcrs(
    counts_df: pd.DataFrame,
    tcr_refs_df: pd.DataFrame,
    non_n_ref_names_set: set,
    region_name_col: str,
    trav_col: str = "TRAV_IMGT_allele_collapsed",
    trbv_col: str = "TRBV_IMGT_allele_collapsed",
    title: str = None,
    outs_dir: Union[str, os.PathLike[str]] = None,
):
    detected_regions_set = set(counts_df.loc[:, region_name_col])
    not_detected_tcr_set = non_n_ref_names_set.difference(detected_regions_set)

    if len(not_detected_tcr_set) == 0:
        return "No missing TCRs!"

    missing_trav_counter = collections.Counter()
    missing_trbv_counter = collections.Counter()

    for missing_tcr in not_detected_tcr_set:
        missing_trav_counter[tcr_refs_df.loc[tcr_refs_df["name"] == missing_tcr, trav_col].values[0]] += 1
        missing_trbv_counter[tcr_refs_df.loc[tcr_refs_df["name"] == missing_tcr, trbv_col].values[0]] += 1

    ax, fig, gs = startfig(
        len([trav_count_tuple[1] for trav_count_tuple in missing_trav_counter.most_common()]) / 1.5,
        7,
    )

    ax.bar(
        np.arange(
            0,
            len([trav_count_tuple[1] for trav_count_tuple in missing_trav_counter.most_common()]),
        ),
        [
            trav_count_tuple[1] / sum([trav_count_tuple[1] for trav_count_tuple in missing_trav_counter.most_common()])
            for trav_count_tuple in missing_trav_counter.most_common()
        ],
        edgecolor="Black",
        linewidth=0.5,
        color="lightblue",
    )

    ax.set_xlim(
        -0.5,
        len([trav_count_tuple[1] for trav_count_tuple in missing_trav_counter.most_common()]) + 0.5,
    )
    ax.set_ylim(0, 1)
    ax.set_xticks(
        np.arange(
            0,
            len([trav_count_tuple[1] for trav_count_tuple in missing_trav_counter.most_common()]),
        )
    )
    ax.set_xticklabels(
        [trav_count_tuple[0] for trav_count_tuple in missing_trav_counter.most_common()],
        rotation=90,
        fontsize=7,
    )
    ax.set_ylabel("Fraction of missing TCRs", fontsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.tick_params(right=False, top=False)
    ax.set_title(
        title
        + ", # missing TCRs = "
        + str(sum([trav_count_tuple[1] for trav_count_tuple in missing_trav_counter.most_common()])),
        fontsize=8,
    )
    fig.tight_layout()
    if outs_dir:
        fig.savefig(os.path.join(outs_dir, "TRAV_counter_missing_tcr_barplot.pdf"))
        plt.close()
    else:
        plt.show()
        plt.close()

    ax, fig, gs = startfig(
        len([trbv_count_tuple[1] for trbv_count_tuple in missing_trbv_counter.most_common()]) / 1.5,
        7,
    )

    ax.bar(
        np.arange(
            0,
            len([trbv_count_tuple[1] for trbv_count_tuple in missing_trbv_counter.most_common()]),
        ),
        [
            trbv_count_tuple[1] / sum([trbv_count_tuple[1] for trbv_count_tuple in missing_trbv_counter.most_common()])
            for trbv_count_tuple in missing_trbv_counter.most_common()
        ],
        edgecolor="Black",
        linewidth=0.5,
        color="lightblue",
    )

    ax.set_xlim(
        -0.5,
        len([trbv_count_tuple[1] for trbv_count_tuple in missing_trbv_counter.most_common()]) + 0.5,
    )
    ax.set_ylim(0, 1)
    ax.set_xticks(
        np.arange(
            0,
            len([trbv_count_tuple[1] for trbv_count_tuple in missing_trbv_counter.most_common()]),
        )
    )
    ax.set_xticklabels(
        [trbv_count_tuple[0] for trbv_count_tuple in missing_trbv_counter.most_common()],
        rotation=90,
        fontsize=7,
    )
    ax.set_ylabel("Fraction of missing TCRs", fontsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.tick_params(right=False, top=False)
    ax.set_title(
        title
        + ", # missing TCRs = "
        + str(sum([trbv_count_tuple[1] for trbv_count_tuple in missing_trbv_counter.most_common()])),
        fontsize=8,
    )
    fig.tight_layout()
    if outs_dir:
        fig.savefig(os.path.join(outs_dir, "TRBV_counter_missing_tcr_barplot.pdf"))
        plt.close()
    else:
        plt.show()
        plt.close()
