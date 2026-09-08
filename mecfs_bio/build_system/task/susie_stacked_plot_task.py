"""
A Task and associated helper methods to plot the results and context of SUSIE fine mapping using stacked panels.
"""

import gc
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns
import structlog
import xarray as xr
from attrs import frozen

# import pandas as pd
from matplotlib import gridspec
from matplotlib.colors import TABLEAU_COLORS, ListedColormap
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from polars import String

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.plot_file_meta import GWASPlotFileMeta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_mlog10p_pipe import (
    ComputeMlog10pIfNeededPipe,
)
from mecfs_bio.build_system.task.pipes.compute_p_from_beta_se import (
    ComputePFromBetaSEPipeIfNeeded,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    COMBINED_CS_FILENAME,
    CS_COLUMN,
    FILTERED_GWAS_FILENAME,
    FILTERED_LD_FILENAME,
    NO_CS_FOUND_FILENAME,
    PIP_COLUMN,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_MLOG10P_COL,
    GWASLAB_POS_COL,
)
from mecfs_bio.util.plotting.save_fig import write_plots_to_dir


@frozen
class BinOptions:
    num_bins: int


logger = structlog.get_logger()

GENE_INFO_START_COL = "gene_start"
GENE_INFO_END_COL = "gene_end"
GENE_INFO_NAME_COL = "gene_name"
GENE_INFO_STRAND_COL = "strand"
GENE_INFO_CHROM_COL = "chrom"

_gwas_pipe = CompositePipe(
    [ComputePFromBetaSEPipeIfNeeded(), ComputeMlog10pIfNeededPipe()]
)

seaborn_rocket_cmap = sns.color_palette("rocket", n_colors=256)

_matplotlib_rocket_cmap = ListedColormap(seaborn_rocket_cmap)


@frozen
class RegionSelectOverride:
    chrom: int
    start: int
    end: int


@frozen
class RegionSelectDefault:
    pass


RegionSelect = RegionSelectOverride | RegionSelectDefault


def get_region(mode: RegionSelect, susie_output_path: Path) -> tuple[int, int, int]:
    """
    Determine the plotting region based on RegionSelect options
    """
    if isinstance(mode, RegionSelectOverride):
        return mode.chrom, mode.start, mode.end
    df = pl.read_parquet(susie_output_path / FILTERED_GWAS_FILENAME)
    assert df[GWASLAB_CHROM_COL].n_unique() == 1
    chrom = df[GWASLAB_CHROM_COL][0]
    start = int(max(int(df[GWASLAB_POS_COL].to_numpy().min()) - 1, 0))
    end = int(df[GWASLAB_POS_COL].to_numpy().max() + 2)
    return (
        chrom,
        start,
        end,
    )


HeatMapPlotMode = Literal["ld2", "ld_abs"]


@frozen
class HeatmapOptions:
    heatmap_bin_options: BinOptions | None
    mode: HeatMapPlotMode
    cmap: str | ListedColormap = "plasma"


@frozen
class SusieStackPlotTask(Task):
    """
    Create a plot to illustrate the results of a SUSIE run on a given locus.
    The resulting plot is a stack of panels showing
    - LD structure
    - marginal associations (i.e. Manhattan plot)
    - SUSIE PIPs
    - Genes
    """

    meta: Meta
    susie_task: Task
    gene_info_task: Task
    gene_info_pipe: DataProcessingPipe
    region_mode: RegionSelect
    heatmap_options: HeatmapOptions

    @property
    def deps(self) -> list["Task"]:
        deps = [self.susie_task, self.gene_info_task]
        return deps

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        susie_asset = fetch(self.susie_task.asset_id)
        assert isinstance(susie_asset, DirectoryAsset)
        susie_dir = susie_asset.path
        if (susie_dir / NO_CS_FOUND_FILENAME).exists():
            logger.debug("No credible sets to plot. Skipping susie panel.")
            (scratch_dir / NO_CS_FOUND_FILENAME).write_text("No credible sets.")
            susie_df = None
        else:
            susie_df = pl.read_parquet(susie_dir / COMBINED_CS_FILENAME)

        gene_info_asset = fetch(self.gene_info_task.asset_id)
        assert isinstance(gene_info_asset, FileAsset)
        gene_info_df = (
            self.gene_info_pipe.process(
                scan_dataframe_asset(gene_info_asset, self.gene_info_task.meta)
            )
            .collect()
            .to_polars()
        )

        chrom, start, end = get_region(self.region_mode, susie_output_path=susie_dir)
        gene_info_df = gene_info_df.filter(
            pl.col(GENE_INFO_CHROM_COL).cast(String()) == pl.lit(chrom).cast(String())
        )
        loaded = pl.read_parquet(susie_dir / FILTERED_GWAS_FILENAME)

        fig = plot_locus_tracks_matplotlib(
            gwas_df=_gwas_pipe.process_eager_polars(loaded),
            susie_cs_df=susie_df,
            ld_np=np.load(susie_dir / FILTERED_LD_FILENAME),
            gene_df=gene_info_df,
            start_bp=start,
            end_bp=end,
            chrom=chrom,
            heatmap_options=self.heatmap_options,
        )
        out_path = scratch_dir
        write_plots_to_dir(out_path, {"plot": fig})
        del fig
        del loaded
        del gene_info_df
        gc.collect()
        return FileAsset(out_path / "plot.png")

    @classmethod
    def create(
        cls,
        asset_id: str,
        susie_task: Task,
        gene_info_task: Task,
        gene_info_pipe: DataProcessingPipe,
        region_mode: RegionSelect,
        heatmap_options: HeatmapOptions,
    ):
        src_meta = susie_task.meta
        assert isinstance(src_meta, ResultDirectoryMeta)
        meta = GWASPlotFileMeta(
            trait=src_meta.trait,
            project=src_meta.project,
            id=AssetId(asset_id),
            extension=".png",
        )
        return cls(
            meta=meta,
            susie_task=susie_task,
            gene_info_task=gene_info_task,
            gene_info_pipe=gene_info_pipe,
            region_mode=region_mode,
            heatmap_options=heatmap_options,
        )


def plot_locus_tracks_matplotlib(
    gwas_df: pl.DataFrame,
    susie_cs_df: pl.DataFrame | None,
    ld_np: np.ndarray,
    gene_df: pl.DataFrame,
    start_bp: int,
    end_bp: int,
    chrom: int,
    heatmap_options: HeatmapOptions,
    # heatmap_bin_options: BinOptions | None,
    *,
    gwas_pos_col: str = GWASLAB_POS_COL,
    gwas_mlog10p_col: str = GWASLAB_MLOG10P_COL,
    susie_pos_col: str = GWASLAB_POS_COL,
    susie_pip_col: str = PIP_COLUMN,
    susie_cs_col: str = CS_COLUMN,
    gene_start_col: str = GENE_INFO_START_COL,
    gene_end_col: str = GENE_INFO_END_COL,
    gene_name_col: str = GENE_INFO_NAME_COL,
    gene_strand_col: str = GENE_INFO_STRAND_COL,
    max_mlog10p: float = 200,
) -> Figure:
    """
    Helper function to create the matplotlib plot consisting of stacked panels

    """

    gwas_df = gwas_df.with_columns(
        pl.min_horizontal(pl.lit(max_mlog10p), pl.col(GWASLAB_MLOG10P_COL)).alias(
            GWASLAB_MLOG10P_COL
        ),
    )

    lead = int(np.argmax(gwas_df[gwas_mlog10p_col]))

    ld2_slice = ld_np[lead, :] ** 2

    # initialize figure with 4 by 2 grid.  Right column is used for legends and colorbars
    fig = plt.figure(figsize=(12, 9))
    gs = gridspec.GridSpec(
        nrows=4,
        ncols=2,
        height_ratios=[1.8, 2.2, 1.6, 1.4],
        width_ratios=[1.0, 0.1],
        hspace=0.08,
        wspace=0.05,
    )

    # colorbars
    ld_cax_container = fig.add_subplot(gs[0, 1])
    ld_cax_container.axis("off")
    manh_cax_container = fig.add_subplot(gs[1, 1])
    manh_cax_container.axis("off")
    ld_cax = inset_axes(
        ld_cax_container,
        width="50%",
        height="50%",
        loc="center left",  # Anchor it to the left side of the container
        borderpad=0,  # No padding between anchor and axis
    )
    manh_cax = inset_axes(
        manh_cax_container,
        width="50%",
        height="50%",
        loc="center left",  # Anchor it to the left side of the container
        borderpad=0,  # No padding between anchor and axis
    )

    # manhattan plot
    ax_manh = draw_manhattan_track(
        fig=fig,
        target_gridspec_cell=gs[1, 0],
        gwas_df=gwas_df,
        ld2_colors=ld2_slice,  # Pass the color array
        gwas_pos_col=gwas_pos_col,
        gwas_mlog10p_col=gwas_mlog10p_col,
        break_at=20.0,
        max_break_proportion=0.5,
        lead=lead,
        colorbar_axis=manh_cax,
    )

    # setup other axis
    ax_ld = fig.add_subplot(gs[0, 0], sharex=ax_manh)
    ax_pip = fig.add_subplot(gs[2, 0], sharex=ax_manh)
    ax_gene = fig.add_subplot(gs[3, 0], sharex=ax_manh)

    pip_legend_ax = fig.add_subplot(gs[2, 1])
    pip_legend_ax.axis("off")

    # 2: SUSIE track
    plot_susie_track(
        susie_cs_df=susie_cs_df,
        pip_legend_ax=pip_legend_ax,
        ax_pip=ax_pip,
        susie_pip_col=susie_pip_col,
        susie_pos_col=susie_pos_col,
        susie_cs_col=susie_cs_col,
    )
    # #3 ld heatmap track
    plot_ld_heatmap(
        ld_np=ld_np,
        gwas_df=gwas_df,
        options=heatmap_options,
        ax_ld=ax_ld,
        fig=fig,
        ld_cax=ld_cax,
        gwas_pos_col=gwas_pos_col,
    )
    # 4: Gene tracks
    plot_gene_tracks(
        ax=ax_gene,
        gene_df=gene_df,
        start_bp=start_bp,
        end_bp=end_bp,
        gene_start_col=gene_start_col,
        gene_end_col=gene_end_col,
        gene_name_col=gene_name_col,
        gene_strand_col=gene_strand_col,
    )

    # axis config
    ax_gene.set_xlim(start_bp - 1, end_bp + 1)
    ax_gene.set_xlabel(f"Chromosome {chrom} (bp)")

    for ax in [ax_manh, ax_pip, ax_ld]:
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.tick_params(
            axis="x",  # changes apply to the x-axis
            which="both",  # major and minor ticks
            bottom=False,  # ticks along the bottom edge are off
            top=False,  # ticks along the top edge are off
            labelbottom=False,  # labels along the bottom edge are off
        )

    for ax in [ax_manh, ax_pip, ax_ld, ax_gene]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if ax == ax_gene:
            ax.spines["left"].set_visible(False)

    return fig


def get_array_and_edges_for_ld_heatmap(
    ld_abs: np.ndarray,
    pos: np.ndarray,
    bin_options: BinOptions | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Use xarray to bin LD data to facilitate creation of an LD heatmap
    """
    da = xr.DataArray(
        ld_abs,
        coords={
            "x": pos,
            "y": pos,
        },
        dims=["x", "y"],
    )
    if bin_options is not None:
        da = (
            da.groupby_bins(
                "x",
                bins=bin_options.num_bins,
            )
            .mean()
            .groupby_bins("y", bins=bin_options.num_bins)
            .mean()
        )
        da = da.rename({"x_bins": "x", "y_bins": "y"})
        da.coords["x"] = [i.mid for i in da.coords["x"].values]
        da.coords["y"] = [i.mid for i in da.coords["y"].values]
    else:
        da = da.groupby("x").mean().groupby("y").mean()

    new_pos = da.coords["x"].to_numpy()
    mid = (new_pos[:-1] + new_pos[1:]) / 2.0
    first = new_pos[0] - (mid[0] - new_pos[0])
    last = new_pos[-1] + (new_pos[-1] - mid[-1])
    edges = np.concatenate([[first], mid, [last]])
    return da.to_numpy(), edges


def plot_gene_tracks(
    ax,
    gene_df: pl.DataFrame,
    start_bp: int,
    end_bp: int,
    gene_start_col: str,
    gene_end_col: str,
    gene_name_col: str,
    gene_strand_col: str,
    font_size: int = 9,
    min_dist_between_genes: float = 0.03,  # Slightly increased buffer
):
    """
    Generated With Gemini
    Plots gene tracks with smart label centering and collision avoidance.
    """
    # 1. Filter genes
    region_df = gene_df.filter(
        (pl.col(gene_end_col) >= start_bp) & (pl.col(gene_start_col) <= end_bp)
    ).sort(gene_start_col)

    if region_df.height == 0:
        ax.text(0.5, 0.5, "No genes in region", transform=ax.transAxes, ha="center")
        ax.set_yticks([])
        return

    # 2. Setup Packing Parameters
    total_bp = end_bp - start_bp
    # Heuristic: How many bp does one character take up?
    # You might need to tune '0.015' based on your specific figure width/dpi
    char_width_factor = 0.015
    bp_per_char = total_bp * char_width_factor
    min_gap = total_bp * min_dist_between_genes

    lanes: list = []
    gene_placements: list = []

    for row in region_df.iter_rows(named=True):
        g_start = max(start_bp, row[gene_start_col])
        g_end = min(end_bp, row[gene_end_col])
        g_name = row[gene_name_col]

        # --- NEW LOGIC: Calculate Text Dimensions & Position ---
        text_width_bp = len(g_name) * bp_per_char
        gene_width_bp = g_end - g_start

        # Decision: Center label or Left-align?
        if gene_width_bp > text_width_bp:
            # Case A: Gene is longer than text -> CENTER text
            text_x = (g_start + g_end) / 2
            ha = "center"
            # Visual end is just the gene end (text is inside)
            visual_end = g_end + min_gap
        else:
            # Case B: Gene is short -> LEFT-ALIGN text at start
            text_x = g_start
            ha = "left"
            # Visual end is the text end
            visual_end = g_start + text_width_bp + min_gap
            # Ensure we don't accidentally clip the gene if they are nearly same size
            visual_end = max(visual_end, g_end + min_gap)

        # --- Packing Algorithm ---
        y_level = -1
        for i, lane_end in enumerate(lanes):
            if lane_end < g_start:
                y_level = i
                lanes[i] = visual_end
                break

        if y_level == -1:
            lanes.append(visual_end)
            y_level = len(lanes) - 1

        gene_placements.append(
            {
                "start": row[gene_start_col],
                "end": row[gene_end_col],
                "name": g_name,
                "strand": row[gene_strand_col],
                "y": y_level,
                "text_x": text_x,
                "ha": ha,
            }
        )

    # 3. Plotting
    n_lanes = len(lanes)
    ax.set_ylim(-0.5, n_lanes + 0.5)

    if n_lanes > 10:
        font_size = max(6, font_size - 2)

    for g in gene_placements:
        y = g["y"]
        # Gene Body
        ax.plot([g["start"], g["end"]], [y, y], color="navy", lw=1.5)

        # Exon/End Ticks
        tick_h = 0.2
        ax.plot([g["start"], g["start"]], [y - tick_h, y + tick_h], color="navy", lw=1)
        ax.plot([g["end"], g["end"]], [y - tick_h, y + tick_h], color="navy", lw=1)

        # Strand Arrow
        mid = (max(start_bp, g["start"]) + min(end_bp, g["end"])) / 2
        marker = ">" if g["strand"] == "+" else "<"
        ax.scatter([mid], [y], marker=marker, color="navy", s=30, zorder=10)

        # Label with Background Box
        # We use a white 'bbox' with alpha=0.8 to hide any lines crossing behind the text
        ax.text(
            g["text_x"],
            y + 0.35,
            g["name"],
            ha=g["ha"],
            va="center",
            fontsize=font_size,
            style="italic",
            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8),
            zorder=20,  # Ensure text sits on top of lines
        )

    ax.set_yticks([])
    ax.invert_yaxis()


def draw_manhattan_track(
    fig,
    target_gridspec_cell,  # e.g. gs[0, 0]
    colorbar_axis,
    gwas_df: pl.DataFrame,
    ld2_colors: np.ndarray,  # Pre-calculated colors for the scatter points
    gwas_pos_col: str,
    lead: int,
    gwas_mlog10p_col: str,
    min_y: int = 2,
    break_at: float = 20.0,
    max_break_proportion: float = 0.5,  # Cap the top section at 50% height
    saturation_point: float = 100.0,  # Point where we reach max height proportion
    significance_threshold: float = 7.8239087,
):
    """
    Generated with Gemini and then modified.

    Draws a Manhattan plot into the provided GridSpec cell.
    Automatically handles axis breaking if values exceed 'break_at'.
    """

    # 1. Determine Dynamic Ratios
    max_val = gwas_df[gwas_mlog10p_col].to_numpy().max()

    # Logic: Calculate what % of height the top part gets.
    # If max_val < 20, top gets 0%.
    # If max_val > 100, top gets 50%.
    # In between, scale linearly.
    if max_val <= break_at:
        top_fraction = 0
    else:
        # Linear interpolation
        slope = max_break_proportion / (saturation_point - break_at)
        top_fraction = slope * (max_val - break_at)
        # Clamp between a minimum visibility (e.g. 0.2) and max_break_proportion
        # If it's barely over 20, we still need enough space to draw ticks, so we set a floor.
        top_fraction = max(0.25, min(max_break_proportion, top_fraction))

    # 2. Case A: Standard Plot (No Break)
    if top_fraction == 0:
        ax = fig.add_subplot(target_gridspec_cell)
        sc = ax.scatter(
            gwas_df[gwas_pos_col],
            gwas_df[gwas_mlog10p_col],
            s=10,
            c=ld2_colors,
            linewidths=0,
            cmap="plasma",
        )
        ax.set_ylabel(r"$-\log_{10}(p)$")

        ax.scatter(
            gwas_df[gwas_pos_col][lead],
            gwas_df[gwas_mlog10p_col][lead],
            s=35,
            marker="^",
            c="black",
        )
        ax.axhline(
            y=significance_threshold,
            color="grey",
            linestyle="--",
            linewidth=1.5,
            label="Significance Threshold",
        )
        ax.set_ylim(min_y, (gwas_df[gwas_mlog10p_col]).to_numpy().max() * 1.05)

    else:
        # 3. Case B: Broken Axis (Nested GridSpec)
        # Create a nested grid INSIDE the target cell
        # height_ratios takes relative weights.
        # If top_fraction is 0.3, bottom is 0.7.
        gs_inner = gridspec.GridSpecFromSubplotSpec(
            nrows=2,
            ncols=1,
            subplot_spec=target_gridspec_cell,
            height_ratios=[top_fraction, 1.0 - top_fraction],
            hspace=0.05,  # Tiny gap for the break marks
        )

        ax_top = fig.add_subplot(gs_inner[0])
        ax_bottom = fig.add_subplot(gs_inner[1], sharex=ax_top)

        # Plot data on BOTH
        # (Optimization: You could filter data for ax_top, but scatter is fast enough usually)
        for ax in [ax_top, ax_bottom]:
            sc = ax.scatter(
                gwas_df[gwas_pos_col],
                gwas_df[gwas_mlog10p_col],
                s=10,
                c=ld2_colors,
                linewidths=0,
                cmap="plasma",
            )

        ax_top.scatter(
            gwas_df[gwas_pos_col][lead],
            gwas_df[gwas_mlog10p_col][lead],
            s=35,
            marker="^",
            c="black",
        )

        ax_bottom.axhline(
            y=significance_threshold,
            color="grey",
            linestyle="--",
            linewidth=1.5,
            label="Significance Threshold",
        )

        # Set Limits
        # Bottom: 0 to 20
        ax_bottom.set_ylim(min_y, break_at)
        # Top: 20 to Max
        ax_top.set_ylim(break_at, max_val * 1.05)

        # Visual Styling for the Break
        # Hide the spines between them
        ax_bottom.spines["top"].set_visible(False)
        ax_bottom.spines["right"].set_visible(False)
        ax_top.spines["bottom"].set_visible(False)
        ax_top.spines["right"].set_visible(False)

        # Remove x-ticks from top
        ax_top.tick_params(labelbottom=False, bottom=False)

        # Draw Diagonal "Cut" Lines
        d = 0.006
        kwargs = dict(transform=ax_top.transAxes, color="k", clip_on=False)
        ax_top.plot((-d, +d), (-d, +d), **kwargs)  # top-left
        ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right

        kwargs.update(transform=ax_bottom.transAxes)
        ax_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left
        ax_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right

        # Label: Center it on the "entire" side using figure coords is best,
        # or just place it on the top axis to be safe.
        ax_top.set_ylabel(r"$-\log_{10}(p)$")
        # Ensure bottom axis label doesn't conflict
        ax_bottom.set_ylabel("")
        ax = ax_bottom
        # Return the bottom axis so other tracks can link their x-axis to it

    cbar = fig.colorbar(sc, cax=colorbar_axis, shrink=0.4)
    cbar.set_label(r"$r^2$ with lead")
    return ax  # Return the single axis for x-linking later


def plot_susie_track(
    susie_cs_df: pl.DataFrame | None,
    ax_pip,
    pip_legend_ax,
    susie_cs_col: str = CS_COLUMN,
    susie_pip_col: str = PIP_COLUMN,
    susie_pos_col: str = GWASLAB_POS_COL,
):
    """
    Plot SUSIE pip using a bar graph colored by credible set
    """
    if susie_cs_df is None:
        return
    pip_traces = []
    palette = list(TABLEAU_COLORS.values())
    if len(susie_cs_df) > 0:
        for i, (cs, sub) in enumerate(
            susie_cs_df.group_by(susie_cs_col, maintain_order=True)
        ):
            color = palette[i % len(palette)]
            ax_pip.vlines(
                sub[susie_pos_col].to_numpy(),
                0.0,
                sub[susie_pip_col].to_numpy(),
                linewidth=1.5,
                label=f"CS {i + 1}",
                color=color,
            )
            handle = Line2D([0], [0], color=color, lw=2, label=f"CS {i + 1}")
            pip_traces.append(handle)
            # pip_trace_labels.append(f"CS {i+1}")
        # ax_pip.legend(loc="upper right", fontsize=8, frameon=False, ncols=2)
        pip_legend_ax.legend(
            handles=pip_traces,
            loc="center left",  # Vertically centered in the panel
            borderaxespad=0,  # Tight alignment to the left edge
            frameon=False,  # Clean look (no box)
            fontsize=9,
        )
    ax_pip.set_ylabel("PIP (SUSIE)")


def plot_ld_heatmap(
    ld_np: np.ndarray,
    gwas_df: pl.DataFrame,
    options: HeatmapOptions,
    ax_ld,
    fig,
    ld_cax,
    gwas_pos_col: str = GWASLAB_POS_COL,
):
    if options.mode == "ld_abs":
        to_plot: np.ndarray = abs(ld_np)
    elif options.mode == "ld2":
        to_plot = ld_np**2
    else:
        raise ValueError(f"Unknown mode: {options.mode}")
    ar, edges = get_array_and_edges_for_ld_heatmap(
        ld_abs=to_plot,
        pos=gwas_df[gwas_pos_col].to_numpy(),
        bin_options=options.heatmap_bin_options,
    )

    mesh = ax_ld.pcolormesh(
        edges, edges, ar, shading="auto", vmin=0, vmax=1, cmap=options.cmap
    )
    ax_ld.set_ylabel("bp")
    ax_ld.set_ylim(float(edges[0]), float(edges[-1]))

    cbar = fig.colorbar(mesh, cax=ld_cax, shrink=0.8)
    label = r"$|r|$" if options.mode == "ld_abs" else "$r^2$"
    cbar.set_label(label)
