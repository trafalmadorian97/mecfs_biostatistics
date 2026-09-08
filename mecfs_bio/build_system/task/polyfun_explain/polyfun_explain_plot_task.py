"""Stacked explainability figure for a polyfun-vs-uniform SUSIE result.

Panels, top to bottom, sharing the genomic-position x-axis:
  1. Manhattan (-log10 p), points colored by LD r^2 with the min-p lead
     variant, with local recombination rate (cM/Mb from the hg19 genetic map)
     on a secondary axis.
  2. PIP, uniform run (credible-set variants only).
  3. PIP, polyfun run (credible-set variants only). Shares its y-scale with the
     uniform panel so the two are directly comparable, and carries a callout on
     each credible set's prior-boosted variant naming its key annotation
     families (from the contrast task's callouts.parquet).
  4. Genes.

Writes both explain_plot.png and explain_plot.svg. Inspired by
SusieStackPlotTask but independent of it.

Figures are built through matplotlib's object-oriented API (a directly
constructed Figure), not pyplot. That keeps rendering off any global backend
state: the process-wide interactive backend is never selected or mutated, the
output format is chosen per file extension by savefig (so the .svg is a true
vector file and the .png raster), and no figure is registered in pyplot's
global manager, so there is nothing to close to reclaim memory.

The contrast_task supplies the per-credible-set callouts (callouts.parquet). The
annotation_parquet_task and ridge_weights_task deps are retained (the contrast
task consumes them upstream); the family-panel and prior-fold tracks that this
plot previously drew from them have been dropped.
"""

from pathlib import Path, PurePath

import narwhals as nw
import numpy as np
import polars as pl
import textalloc as ta
from attrs import frozen
from matplotlib.figure import Figure

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (
    scan_dataframe_asset,
)
from mecfs_bio.build_system.meta.result_directory_meta import (
    ResultDirectoryMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import (
    DataProcessingPipe,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    CALLOUT_LABEL_COL,
    CALLOUT_PIP_PF_COL,
    CALLOUTS_FILENAME,
    _load_run_variants,
)
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    COMBINED_CS_FILENAME,
    FILTERED_GWAS_FILENAME,
    FILTERED_LD_FILENAME,
    PIP_COLUMN,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
    GENE_INFO_STRAND_COL,
    plot_gene_tracks,
    plot_susie_track,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.genetic_map_constants import (
    GMAP_POS_COL,
    GMAP_RATE_COL,
)
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)

PLOT_PNG_FILENAME = "explain_plot.png"
PLOT_SVG_FILENAME = "explain_plot.svg"
# PIP units of empty space kept above the tallest stem for the callout labels.
_PIP_LABEL_HEADROOM = 0.4


@frozen
class PolyfunExplainPlotTask(Task):
    """Render the polyfun-vs-uniform explainability figure."""

    meta: Meta
    susie_uniform_task: Task
    susie_polyfun_task: Task
    contrast_task: Task
    annotation_parquet_task: Task
    gene_info_task: Task
    ridge_weights_task: Task
    genetic_map_task: Task
    genome_build: GenomeBuild = "19"
    gene_info_pipe: DataProcessingPipe = IdentityPipe()

    @property
    def deps(self) -> list["Task"]:
        return [
            self.susie_uniform_task,
            self.susie_polyfun_task,
            self.contrast_task,
            self.ridge_weights_task,
            self.annotation_parquet_task,
            self.gene_info_task,
            self.genetic_map_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        uni_dir = _dir(fetch, self.susie_uniform_task)
        pf_dir = _dir(fetch, self.susie_polyfun_task)

        pf_variants = _load_run_variants(pf_dir).sort(GWASLAB_POS_COL)
        # Credible-set membership per run: the PIP panels plot only these variants,
        # colored/legended by credible set (empty frame if the run found none).
        uni_cs = pl.read_parquet(uni_dir / COMBINED_CS_FILENAME)
        pf_cs = pl.read_parquet(pf_dir / COMBINED_CS_FILENAME)
        pf_full = pl.read_parquet(pf_dir / FILTERED_GWAS_FILENAME)
        ld = np.load(pf_dir / FILTERED_LD_FILENAME)
        contrast_dir = _dir(fetch, self.contrast_task)
        callouts = pl.read_parquet(contrast_dir / CALLOUTS_FILENAME)

        chrom = int(pf_variants[GWASLAB_CHROM_COL][0])
        bp_min = int(pf_variants[GWASLAB_POS_COL].to_numpy().min())
        bp_max = int(pf_variants[GWASLAB_POS_COL].to_numpy().max())

        genes = (
            self.gene_info_pipe.process(
                scan_dataframe_asset(
                    fetch(self.gene_info_task.asset_id), self.gene_info_task.meta
                )
            )
            .collect()
            .to_polars()
        )

        _render(
            scratch_dir=scratch_dir,
            uni_cs=uni_cs,
            pf_cs=pf_cs,
            pf_full=pf_full,
            callouts=callouts,
            ld=ld,
            genes=genes,
            fetch=fetch,
            genetic_map_task=self.genetic_map_task,
            chrom=chrom,
            bp_min=bp_min,
            bp_max=bp_max,
            genome_build=self.genome_build,
        )
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        susie_uniform_task: Task,
        susie_polyfun_task: Task,
        contrast_task: Task,
        annotation_parquet_task: Task,
        gene_info_task: Task,
        ridge_weights_task: Task,
        genetic_map_task: Task,
        genome_build: GenomeBuild = "19",
        gene_info_pipe: DataProcessingPipe = IdentityPipe(),
    ) -> "PolyfunExplainPlotTask":
        source_meta = susie_polyfun_task.meta
        if not isinstance(source_meta, ResultDirectoryMeta):
            raise ValueError(f"Unknown meta for polyfun susie task: {source_meta}")
        meta = ResultDirectoryMeta(
            id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=PurePath("analysis"),
        )
        return cls(
            meta=meta,
            susie_uniform_task=susie_uniform_task,
            susie_polyfun_task=susie_polyfun_task,
            contrast_task=contrast_task,
            annotation_parquet_task=annotation_parquet_task,
            gene_info_task=gene_info_task,
            ridge_weights_task=ridge_weights_task,
            genetic_map_task=genetic_map_task,
            genome_build=genome_build,
            gene_info_pipe=gene_info_pipe,
        )


def _dir(fetch: Fetch, task: Task) -> Path:
    asset = fetch(task.asset_id)
    assert isinstance(asset, DirectoryAsset)
    return asset.path


def _norm_sf(z: np.ndarray) -> np.ndarray:
    from scipy.stats import norm

    return norm.sf(z)


def _load_recomb(
    fetch: Fetch, task: Task, chrom: int, bp_min: int, bp_max: int
) -> pl.DataFrame:
    """Locus-windowed recombination rate (cM/Mb) from the hg19 genetic map."""
    return (
        scan_dataframe_asset(fetch(task.asset_id), task.meta)
        .filter(
            (nw.col(GWASLAB_CHROM_COL) == chrom)
            & (nw.col(GMAP_POS_COL) >= bp_min)
            & (nw.col(GMAP_POS_COL) <= bp_max)
        )
        .select(GWASLAB_CHROM_COL, GMAP_POS_COL, GMAP_RATE_COL)
        .collect()
        .to_polars()
        .sort(GMAP_POS_COL)
    )


def _plot_recomb(ax, recomb: pl.DataFrame) -> None:
    bp = recomb[GMAP_POS_COL].to_numpy().astype(float)
    rate = recomb[GMAP_RATE_COL].to_numpy().astype(float)
    if len(bp) >= 1:
        (line,) = ax.plot(bp, rate, color="tab:red", alpha=0.4, linewidth=1)
        line.set_rasterized(True)
    # The recomb line is red, but its axis label/ticks stay black to match the
    # other panels' axis labels (the color only carries meaning on the line).
    ax.set_ylabel("cM/Mb")


def _render(
    scratch_dir: Path,
    uni_cs: pl.DataFrame,
    pf_cs: pl.DataFrame,
    pf_full: pl.DataFrame,
    callouts: pl.DataFrame,
    ld: np.ndarray,
    genes: pl.DataFrame,
    fetch: Fetch,
    genetic_map_task: Task,
    chrom: int,
    bp_min: int,
    bp_max: int,
    genome_build: GenomeBuild,
) -> None:
    # manhattan + 2 pip + genes, one panel each.
    n_panels = 1 + 2 + 1
    fig = Figure(figsize=(12, 1.7 * n_panels))
    # Left column holds the tracks; the narrow right column is reserved for
    # legends/colorbars so nothing overlaps the data (mirrors SusieStackPlotTask).
    gs = fig.add_gridspec(
        nrows=n_panels,
        ncols=2,
        width_ratios=[1.0, 0.15],
        hspace=0.12,
        wspace=0.02,
    )
    ax0 = fig.add_subplot(gs[0, 0])
    axes = [ax0] + [fig.add_subplot(gs[i, 0], sharex=ax0) for i in range(1, n_panels)]
    x = pf_full[GWASLAB_POS_COL].to_numpy()

    # Panel 1: Manhattan colored by LD with lead (min-p) variant + recomb rate.
    z = (pf_full[GWASLAB_BETA_COL] / pf_full[GWASLAB_SE_COL]).to_numpy()
    neglogp = -np.log10(2.0 * _norm_sf(np.abs(z)))
    lead = int(np.argmax(np.abs(z)))
    r2 = ld[lead, :] ** 2
    sc = ax0.scatter(x, neglogp, c=r2, cmap="viridis", vmin=0, vmax=1, s=10)
    sc.set_rasterized(True)
    ax0.set_ylabel("-log10 p")
    # Right-column cell for panel 1: a shortened colorbar in the upper portion
    # (anchored right, clear of the recomb axis's ticks/label), and below it a
    # small legend telling the reader the red twin-axis line is recombination
    # rate.
    cbar_cell = fig.add_subplot(gs[0, 1])
    cbar_cell.axis("off")
    # Colorbar and recomb key share one aligned column in the RIGHT half of the
    # cell (the left half holds the cM/Mb twin-axis ticks and label): the colorbar
    # bar occupies [cbar_x0, cbar_x0+cbar_w] and the recomb line segment sits
    # directly beneath it over the same x-extent, so the two keys line up.
    cbar_x0, cbar_w = 0.52, 0.20
    cax = cbar_cell.inset_axes((cbar_x0, 0.42, cbar_w, 0.52))
    fig.colorbar(sc, cax=cax, label="r$^2$ w/ lead")

    recomb = _load_recomb(fetch, genetic_map_task, chrom, bp_min, bp_max)
    ax0b = ax0.twinx()
    _plot_recomb(ax0b, recomb)
    cbar_cell.plot(
        [cbar_x0, cbar_x0 + cbar_w],
        [0.24, 0.24],
        transform=cbar_cell.transAxes,
        color="tab:red",
        alpha=0.6,
        linewidth=1.2,
    )
    cbar_cell.text(
        cbar_x0 + cbar_w / 2.0,
        0.14,
        "recomb rate",
        transform=cbar_cell.transAxes,
        ha="center",
        va="top",
        fontsize=7,
    )

    # PIP uniform / polyfun as vertical stems restricted to credible-set variants,
    # colored per credible set with a legend in the right column (reuses the
    # stackplot's susie track). Empty frame -> no stems, no legend.
    _plot_pip_panel(fig, gs, 1, axes, uni_cs, "PIP (uniform)")
    _plot_pip_panel(fig, gs, 2, axes, pf_cs, "PIP (polyfun)")
    # Share one y-scale across both PIP panels so their stem heights are directly
    # comparable (PIP in [0, 1]; scale to the taller of the two, else full range).
    pip_top = _shared_pip_top(uni_cs, pf_cs)
    # Reserve a fixed band above the tallest stem for the callout labels; raise
    # BOTH PIP panels equally so their data scale stays shared (directly
    # comparable). A fixed (not stem-proportional) headroom keeps the band just
    # large enough for labels even when a stem already reaches PIP ~1.
    label_top = pip_top + _PIP_LABEL_HEADROOM
    axes[1].set_ylim(0.0, label_top)
    axes[2].set_ylim(0.0, label_top)
    # PIP is a probability, so never draw ticks/labels above 1.0 even though the
    # panel extends higher to fit the callouts.
    _cap_pip_ticks(axes[1])
    _cap_pip_ticks(axes[2])
    # Place the callout labels in the headroom above the stems, angled off to the
    # side so the leader line is clearly distinct from a vertical PIP stem.
    _place_callouts(
        axes[2],
        callouts,
        pf_cs,
        xlims=(float(bp_min), float(bp_max)),
        ylims=(0.0, label_top),
    )

    # Genes: reuse the stackplot's lane-packed gene track. Filter to this
    # chromosome first (the reference lists every chromosome; the helper windows
    # only by position).
    ax_gene = axes[-1]
    genes_chrom = genes.filter(
        pl.col(GENE_INFO_CHROM_COL).cast(pl.String) == str(chrom)
    )
    plot_gene_tracks(
        ax=ax_gene,
        gene_df=genes_chrom,
        start_bp=bp_min,
        end_bp=bp_max,
        gene_start_col=GENE_INFO_START_COL,
        gene_end_col=GENE_INFO_END_COL,
        gene_name_col=GENE_INFO_NAME_COL,
        gene_strand_col=GENE_INFO_STRAND_COL,
    )
    ax_gene.set_ylabel("genes")
    ax_gene.set_xlabel(f"hg{genome_build} chr{chrom} position (bp)")

    # Lock every panel to the locus window and tidy the shared x-axis: only the
    # bottom (genes) panel keeps tick labels; drop top+right spines throughout.
    # The gene panel additionally drops its left spine so it has no vertical
    # frame lines at all (matching the stackplot).
    ax0.set_xlim(bp_min, bp_max)
    for ax in axes[:-1]:
        ax.tick_params(axis="x", which="both", labelbottom=False, bottom=False)
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    ax_gene.spines["left"].set_visible(False)

    fig.savefig(scratch_dir / PLOT_PNG_FILENAME, dpi=150, bbox_inches="tight")
    fig.savefig(scratch_dir / PLOT_SVG_FILENAME, bbox_inches="tight")


def _shared_pip_top(uni_cs: pl.DataFrame, pf_cs: pl.DataFrame) -> float:
    """Common y-axis top for both PIP panels: a little above the tallest stem
    across the two runs, or the full [0, 1] range when neither has a credible
    set."""
    tops = [
        float(df[PIP_COLUMN].to_numpy().max())
        for df in (uni_cs, pf_cs)
        if df.height > 0 and PIP_COLUMN in df.columns
    ]
    if not tops:
        return 1.0
    return min(1.0, max(tops) * 1.05)


def _cap_pip_ticks(ax) -> None:
    """Keep the PIP axis's ticks/labels within [0, 1]: the panel extends above 1
    to make room for callout labels, but PIP is a probability so a tick above 1
    would misread. Drops any auto tick outside [0, 1] (a no-op for low-PIP panels
    whose ticks already stay under 1)."""
    ax.set_yticks([t for t in ax.get_yticks() if -1e-9 <= t <= 1.0 + 1e-9])


def _plot_pip_panel(
    fig: Figure,
    gs,
    row: int,
    axes: list,
    cs_df: pl.DataFrame,
    label: str,
) -> None:
    """Draw one PIP panel: credible-set-colored stems on the left-column axis and
    a per-credible-set legend in the reserved right-column cell."""
    ax_pip = axes[row]
    legend_ax = fig.add_subplot(gs[row, 1])
    legend_ax.axis("off")
    plot_susie_track(susie_cs_df=cs_df, ax_pip=ax_pip, pip_legend_ax=legend_ax)
    ax_pip.set_ylabel(label)


def _place_callouts(
    ax_pf,
    callouts: pl.DataFrame,
    stem_df: pl.DataFrame,
    xlims: tuple[float, float],
    ylims: tuple[float, float],
) -> None:
    """Annotate the polyfun PIP panel: one text label per callout row, anchored at
    (POS, pip_pf). textalloc treats every PIP stem (0 -> pip) as an obstacle line
    and places each label in the free space above them, angled north-east so the
    leader line meets the box corner and reads distinctly from the vertical stems.
    xlims/ylims are the panel's true axis limits (textalloc normalizes distances
    against them). Empty frame -> no-op."""
    if callouts.height == 0:
        return
    xs = callouts[GWASLAB_POS_COL].to_numpy().astype(float).tolist()
    ys = callouts[CALLOUT_PIP_PF_COL].to_numpy().astype(float).tolist()
    # A long (many-family) callout is wrapped one family per line at a smaller
    # font, keeping its box narrow so textalloc can seat it without crossing a
    # neighbouring stem; short callouts render unchanged (font 7, single line).
    wrapped = [_wrap_callout_label(t) for t in callouts[CALLOUT_LABEL_COL].to_list()]
    texts = [text for text, _ in wrapped]
    sizes = [size for _, size in wrapped]
    # Every stem (vertical line from 0 to its PIP) is an obstacle to route labels
    # and leader lines around.
    stem_x = stem_df[GWASLAB_POS_COL].to_numpy().astype(float)
    stem_pip = stem_df[PIP_COLUMN].to_numpy().astype(float)
    x_lines: list[np.ndarray | list[float]] = [[float(px), float(px)] for px in stem_x]
    y_lines: list[np.ndarray | list[float]] = [[0.0, float(pp)] for pp in stem_pip]
    # Seed so textalloc's candidate search is reproducible across builds (keeps
    # the committed SVG stable); textalloc draws from numpy's global RNG.
    np.random.seed(0)
    ta.allocate(
        ax_pf,
        xs,
        ys,
        texts,
        x_scatter=stem_x.tolist(),
        y_scatter=stem_pip.tolist(),
        x_lines=x_lines,
        y_lines=y_lines,
        textsize=sizes,
        linecolor="black",
        linewidth=0.6,
        direction="northeast",
        min_distance=0.03,
        max_distance=0.5,
        xlims=xlims,
        ylims=ylims,
        avoid_label_lines_overlap=True,
    )


def _wrap_callout_label(label: str, max_chars: int = 45) -> tuple[str, float]:
    """Lay out one callout label. A label longer than max_chars (a variant with
    several or verbose annotation families) is wrapped to one family per line at a
    smaller font, so its box is only as wide as the longest single family and can
    be placed clear of neighbouring stems. Shorter labels are returned unchanged.

    Parses the "pos:nea:ea (fam ++, fam +, ...)" form produced by the contrast
    task; anything not matching that shape is returned as-is."""
    if len(label) <= max_chars or " (" not in label or not label.endswith(")"):
        return label, 7.0
    head, inner = label.split(" (", 1)
    families = inner[:-1].split(", ")
    return "\n".join([head, *families]), 6.5
