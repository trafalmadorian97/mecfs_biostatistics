"""Asset generator for polyfun explainability fine-mapping.

Inner generator: given a locus's shared inputs and one run config, produce a
matched pair of SUSIE runs (polyfun precomputed prior + uniform) plus the
contrast and plot tasks that explain the pair in annotation terms. Outer
generator: build the per-locus shared inputs (LD interval, renamed LD labels,
harmonized sumstats), then call the inner generator for each of the four run
configs (L=1, L=2, L=10, L=10-strict) -> 8 SUSIE runs per locus.

The per-locus shared setup (LD interval lookup, label renaming, harmonization)
mirrors generate_assets_broad_ukbb_fine_map but is inlined here as a private
helper so the existing, separately tested fine-mapping generator is left
untouched.
"""

from pathlib import PurePath
from typing import Mapping

import structlog
from attrs import frozen

from mecfs_bio.asset_generator.ukbb_broad_ld_matrix_generator import (
    get_genomic_interval_stem_name,
    get_ld_labels_and_matrix_task_for_genomic_interval_build_37,
    get_optimal_ukbb_ld_interval,
)
from mecfs_bio.assets.reference_data.genetic_map.genetic_map_hg19 import (
    GENETIC_MAP_HG19,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.raw.magma_ensembl_gene_location_reference_data_build_37 import (
    MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
)
from mecfs_bio.assets.reference_data.polyfun.annotations.annotation_ridge_weights import (
    BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS,
)
from mecfs_bio.assets.reference_data.polyfun.annotations.baseline_lf_annotations import (
    BASELINE_LF_ANNOTATION_MATRIX,
)
from mecfs_bio.assets.reference_data.polyfun.precomputed_prior.polyfun_precomputed_prior import (
    COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS,
    POLYFUN_PRIOR_COL,
    create_prior_col_pipe,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.copy_file_from_directory_task import (
    CopyFileFromDirectoryTask,
)
from mecfs_bio.build_system.task.dataframe_output import (
    ParquetOutFormat,
)
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_chrom_pos_alleles import (
    ChromRange,
    HarmonizeGWASWithReferenceViaAlleles,
)
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_rsid import (
    PalindromeStrategy,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.concat_str_pipe import ConcatStrPipe
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.pipes.min_variants_for_cumulative_mass import (
    MinVariantsForCumulativeMass,
)
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    DETAILED_DISPLAY_TABLE_FILENAME,
    PER_VARIANT_ANNOTATION_TABLE_FILENAME,
    TOP_LINE_DISPLAY_TABLE_FILENAME,
    PolyfunExplainContrastTask,
    SecondaryPositionFromSnpid,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_plot_task import (
    PLOT_PNG_FILENAME,
    PLOT_SVG_FILENAME,
    PolyfunExplainPlotTask,
)
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    COMBINED_CS_FILENAME,
    CS_COLUMN,
    PIP_COLUMN,
    BroadInstituteFormatLDMatrix,
    PriorInfo,
    SusieRFinemapTask,
)
from mecfs_bio.build_system.task.upset_plot_task import (
    DirSetSource,
    UpSetPlotTask,
)
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)

logger = structlog.get_logger()


@frozen
class RunConfig:
    """One SUSIE parameterization applied to both members of a matched pair."""

    label: str
    max_credible_sets: int
    z_score_filtering_threshold: float = 2.0


RUN_CONFIGS: tuple[RunConfig, ...] = (
    RunConfig(label="l1", max_credible_sets=1),
    RunConfig(label="l2", max_credible_sets=2),
    RunConfig(label="l10", max_credible_sets=10),
    RunConfig(
        label="l10_strict", max_credible_sets=10, z_score_filtering_threshold=1.0
    ),
)


@frozen
class SharedFineMapInputs:
    """Per-locus inputs shared by every run config's matched pair."""

    base_name: str
    harmonized_sumstats_task: Task
    ld_labels_task: Task
    ld_matrix_task: Task
    gene_info_task: Task
    effective_sample_size: int
    genome_build: GenomeBuild = "19"
    q_factor: int = 100
    secondary_position_from_snpid: SecondaryPositionFromSnpid | None = None


@frozen
class PolyfunExplainGroup:
    """A matched uniform/polyfun SUSIE pair plus the tasks explaining it.

    The SUSIE fields are typed concretely (not the base Task) so consumers can
    reach run-specific attributes such as prior_info.
    """

    susie_uniform: SusieRFinemapTask
    susie_polyfun: SusieRFinemapTask
    contrast: Task
    plot: Task
    # The plot's png and svg copied out of its directory as standalone FileAssets,
    # so docs can include either figure format on its own.
    plot_png: Task
    plot_svg: Task
    # The two display tables copied out of the contrast directory as standalone
    # FileAssets, so docs can include either table without pulling in the
    # contrast task's other (detail) outputs.
    top_line_table: Task
    detailed_table: Task
    # The granular per-variant annotation characterization table, copied out of the
    # contrast directory as a standalone FileAsset.
    per_variant_annotation_table: Task
    label: str


@frozen
class PolyfunExplainOuterGroup:
    """All matched pairs (one per run config) for a single locus, plus the two
    UpSet plots comparing the polyfun runs' credible-set variants across the four
    run configs."""

    groups: list[PolyfunExplainGroup]
    upset_all_polyfun: Task
    upset_cs50_polyfun: Task

    def terminal_tasks(self) -> list[Task]:
        out: list[Task] = []
        for g in self.groups:
            out += [
                # g.susie_uniform,
                # g.susie_polyfun,
                # g.contrast,
                g.plot_png,
                g.plot_svg,
                g.top_line_table,
                g.detailed_table,
                g.per_variant_annotation_table,
            ]
        out += [self.upset_all_polyfun, self.upset_cs50_polyfun]
        return out

    @property
    def groups_by_label(self) -> Mapping[str, PolyfunExplainGroup]:
        return {group.label: group for group in self.groups}


def generate_polyfun_explain_group(
    shared: SharedFineMapInputs, config: RunConfig
) -> PolyfunExplainGroup:
    """Build one matched pair (uniform vs polyfun prior) under one run config,
    with the contrast and plot tasks explaining that pair."""
    stem = f"{shared.base_name}_{config.label}"
    prior_info = PriorInfo(
        prior_task=COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS,
        prior_pipe=create_prior_col_pipe(shared.q_factor),
        prior_col=POLYFUN_PRIOR_COL,
    )
    susie_uniform = SusieRFinemapTask.create(
        asset_id=f"{stem}_susie_uniform",
        gwas_data_task=shared.harmonized_sumstats_task,
        ld_labels_task=shared.ld_labels_task,
        ld_matrix_source=BroadInstituteFormatLDMatrix(shared.ld_matrix_task),
        effective_sample_size=shared.effective_sample_size,
        max_credible_sets=config.max_credible_sets,
        z_score_filtering_threshold=config.z_score_filtering_threshold,
        prior_info=None,
    )
    susie_polyfun = SusieRFinemapTask.create(
        asset_id=f"{stem}_susie_polyfun",
        gwas_data_task=shared.harmonized_sumstats_task,
        ld_labels_task=shared.ld_labels_task,
        ld_matrix_source=BroadInstituteFormatLDMatrix(shared.ld_matrix_task),
        effective_sample_size=shared.effective_sample_size,
        max_credible_sets=config.max_credible_sets,
        z_score_filtering_threshold=config.z_score_filtering_threshold,
        prior_info=prior_info,
    )
    contrast = PolyfunExplainContrastTask.create(
        asset_id=f"{stem}_explain_contrast",
        susie_uniform_task=susie_uniform,
        susie_polyfun_task=susie_polyfun,
        ridge_weights_task=BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS,
        annotation_parquet_task=BASELINE_LF_ANNOTATION_MATRIX,
        secondary_position=shared.secondary_position_from_snpid,
    )
    plot = PolyfunExplainPlotTask.create(
        asset_id=f"{stem}_explain_plot",
        susie_uniform_task=susie_uniform,
        susie_polyfun_task=susie_polyfun,
        contrast_task=contrast,
        annotation_parquet_task=BASELINE_LF_ANNOTATION_MATRIX,
        gene_info_task=shared.gene_info_task,
        ridge_weights_task=BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS,
        genetic_map_task=GENETIC_MAP_HG19,
        genome_build=shared.genome_build,
        gene_info_pipe=IdentityPipe(),
    )
    plot_png = CopyFileFromDirectoryTask.create_from_result_plot(
        asset_id=f"{stem}_explain_plot_png",
        source_directory_task=plot,
        path_inside_directory=PurePath(PLOT_PNG_FILENAME),
        extension=".png",
    )
    plot_svg = CopyFileFromDirectoryTask.create_from_result_plot(
        asset_id=f"{stem}_explain_plot_svg",
        source_directory_task=plot,
        path_inside_directory=PurePath(PLOT_SVG_FILENAME),
        extension=".svg",
    )
    top_line_table = CopyFileFromDirectoryTask.create_result_table(
        asset_id=f"{stem}_explain_top_line_table",
        source_directory_task=contrast,
        path_inside_directory=PurePath(TOP_LINE_DISPLAY_TABLE_FILENAME),
        extension=".parquet",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    )
    detailed_table = CopyFileFromDirectoryTask.create_result_table(
        asset_id=f"{stem}_explain_detailed_table",
        source_directory_task=contrast,
        path_inside_directory=PurePath(DETAILED_DISPLAY_TABLE_FILENAME),
        extension=".parquet",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    )
    per_variant_annotation_table = CopyFileFromDirectoryTask.create_result_table(
        asset_id=f"{stem}_explain_per_variant_annotation_table",
        source_directory_task=contrast,
        path_inside_directory=PurePath(PER_VARIANT_ANNOTATION_TABLE_FILENAME),
        extension=".parquet",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    )
    return PolyfunExplainGroup(
        susie_uniform=susie_uniform,
        susie_polyfun=susie_polyfun,
        contrast=contrast,
        plot=plot,
        plot_png=plot_png,
        plot_svg=plot_svg,
        top_line_table=top_line_table,
        detailed_table=detailed_table,
        per_variant_annotation_table=per_variant_annotation_table,
        label=config.label,
    )


def build_explainability_groups(
    shared: SharedFineMapInputs, configs: tuple[RunConfig, ...] = RUN_CONFIGS
) -> list[PolyfunExplainGroup]:
    """One matched pair (+ contrast + plot) per run config."""
    return [generate_polyfun_explain_group(shared, c) for c in configs]


_UPSET_VARIANT_ID_COL = "__variant_id"
_RUN_CONFIG_DISPLAY: dict[str, str] = {
    "l1": "L=1",
    "l2": "L=2",
    "l10": "L=10",
    "l10_strict": "L=10 strict",
}


def _polyfun_cs_variant_sources(
    groups: list[PolyfunExplainGroup],
    configs: tuple[RunConfig, ...],
    row_selection_pipe: DataProcessingPipe | None,
) -> list[DirSetSource]:
    """One UpSet set per run config, reading that config's polyfun combined_cs and
    exposing a per-variant id (chr__pos__ea__nea). row_selection_pipe, if given,
    trims the credible-set rows before the id is built (e.g. to the 50% credible
    set)."""
    assert len(groups) == len(configs)
    id_pipe = ConcatStrPipe(
        target_cols=[
            GWASLAB_CHROM_COL,
            GWASLAB_POS_COL,
            GWASLAB_EFFECT_ALLELE_COL,
            GWASLAB_NON_EFFECT_ALLELE_COL,
        ],
        sep="__",
        new_col_name=_UPSET_VARIANT_ID_COL,
    )
    pipe: DataProcessingPipe = (
        id_pipe
        if row_selection_pipe is None
        else CompositePipe([row_selection_pipe, id_pipe])
    )
    return [
        DirSetSource(
            name=_RUN_CONFIG_DISPLAY[cfg.label],
            task=g.susie_polyfun,
            file_in_dir=PurePath(COMBINED_CS_FILENAME),
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            pipe=pipe,
            col_name=_UPSET_VARIANT_ID_COL,
        )
        for cfg, g in zip(configs, groups)
    ]


def build_polyfun_upset_tasks(
    base_name: str,
    groups: list[PolyfunExplainGroup],
    configs: tuple[RunConfig, ...] = RUN_CONFIGS,
) -> tuple[Task, Task]:
    """Two UpSet plots over the four polyfun runs' credible-set variants:

    - all_cs: every variant in every credible set of each run.
    - cs50: the per-credible-set 50% credible set (the minimal highest-PIP
      variants whose cumulative PIP first reaches 0.5, unioned across sets),
      which trims the long low-PIP tail.
    """
    upset_all = UpSetPlotTask.create(
        asset_id=base_name + "_polyfun_upset_all_cs_variants",
        set_sources=_polyfun_cs_variant_sources(groups, configs, None),
    )
    upset_cs50 = UpSetPlotTask.create(
        asset_id=base_name + "_polyfun_upset_cs50_variants",
        set_sources=_polyfun_cs_variant_sources(
            groups,
            configs,
            MinVariantsForCumulativeMass(
                group_col=CS_COLUMN, value_col=PIP_COLUMN, threshold=0.5
            ),
        ),
    )
    return upset_all, upset_cs50


def _build_shared_locus_inputs(
    chrom: int,
    pos: int,
    build_37_sumstats_task: Task,
    base_name: str,
    sumstats_pipe: DataProcessingPipe,
    sample_size: int,
    gene_info_task: Task,
    q_factor: int,
    chrom_range: ChromRange | None,
    palindrome_strategy: PalindromeStrategy,
    genome_build: GenomeBuild,
    secondary_position_from_snpid: SecondaryPositionFromSnpid | None,
) -> SharedFineMapInputs:
    """Per-locus shared setup: LD interval lookup, LD-label renaming, and
    harmonization of the sumstats against the renamed labels. Mirrors the inline
    setup in generate_assets_broad_ukbb_fine_map (kept here rather than shared so
    that generator stays untouched)."""
    interval = get_optimal_ukbb_ld_interval(chrom=chrom, pos=pos)
    if chrom_range is not None:
        assert chrom == chrom_range.chrom
        assert pos >= chrom_range.start
        assert pos <= chrom_range.end
        assert interval.start <= chrom_range.start
        assert chrom_range.end <= interval.end
        base_name = (
            base_name + f"chr{chrom_range.chrom}_{chrom_range.start}_{chrom_range.end}"
        )
    else:
        base_name = base_name + "_" + get_genomic_interval_stem_name(interval)
    if palindrome_strategy != "drop":
        base_name = base_name + "_palindromes_" + palindrome_strategy

    logger.debug(
        f"To finemap position {pos} on chromosome {chrom}, interval {interval} was selected."
    )
    ld_labels_task, ld_matrix_task = (
        get_ld_labels_and_matrix_task_for_genomic_interval_build_37(interval=interval)
    )
    ld_labels_task_renamed = PipeDataFrameTask.create(
        source_task=ld_labels_task,
        asset_id=ld_labels_task.asset_id + "_renamed",
        out_format=ParquetOutFormat(),
        pipes=[
            RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
            RenameColPipe(old_name="chromosome", new_name=GWASLAB_CHROM_COL),
            RenameColPipe(old_name="position", new_name=GWASLAB_POS_COL),
            # allele1 is the non-effect allele in the Broad UKBB LD panel; see
            # https://github.com/omerwe/polyfun/issues/208#issuecomment-2563832487
            RenameColPipe(old_name="allele1", new_name=GWASLAB_NON_EFFECT_ALLELE_COL),
            RenameColPipe(old_name="allele2", new_name=GWASLAB_EFFECT_ALLELE_COL),
        ],
        backend="polars",
    )
    harmonized_sumstats_task = HarmonizeGWASWithReferenceViaAlleles.create(
        asset_id=base_name + "_gwas_harmonized_with_ref",
        gwas_data_task=build_37_sumstats_task,
        reference_task=ld_labels_task_renamed,
        palindrome_strategy=palindrome_strategy,
        gwas_pipe=CompositePipe(
            [
                sumstats_pipe,
                UniquePipe(
                    by=[
                        GWASLAB_CHROM_COL,
                        GWASLAB_POS_COL,
                        GWASLAB_EFFECT_ALLELE_COL,
                        GWASLAB_NON_EFFECT_ALLELE_COL,
                    ],
                    keep="none",
                    order_by=[
                        GWASLAB_CHROM_COL,
                        GWASLAB_POS_COL,
                        GWASLAB_EFFECT_ALLELE_COL,
                        GWASLAB_NON_EFFECT_ALLELE_COL,
                    ],
                ),
            ]
        ),
        chrom_range_filter=chrom_range,
    )
    return SharedFineMapInputs(
        base_name=base_name,
        harmonized_sumstats_task=harmonized_sumstats_task,
        ld_labels_task=ld_labels_task_renamed,
        ld_matrix_task=ld_matrix_task,
        gene_info_task=gene_info_task,
        effective_sample_size=sample_size,
        genome_build=genome_build,
        q_factor=q_factor,
        secondary_position_from_snpid=secondary_position_from_snpid,
    )


def generate_assets_polyfun_explain_fine_map(
    chrom: int,
    pos: int,
    build_37_sumstats_task: Task,
    base_name: str,
    sumstats_pipe: DataProcessingPipe,
    sample_size_or_effect_sample_size: int,
    gene_info_task: Task = MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
    q_factor: int = 100,
    chrom_range: ChromRange | None = None,
    palindrome_strategy: PalindromeStrategy = "drop",
    secondary_position_from_snpid: SecondaryPositionFromSnpid | None = None,
) -> PolyfunExplainOuterGroup:
    """Build the full explainability asset set for one locus: the per-locus
    shared inputs, then a matched uniform/polyfun SUSIE pair (+ contrast + plot)
    for each of the four run configs (8 SUSIE runs total).

    secondary_position_from_snpid, when given, adds a build-labelled secondary
    position column (e.g. pos_hg38) to the display tables, parsed from the
    variants' SNPID. Only correct when the SNPID position field is in the
    asserted build (true for gwaslab sumstats lifted over from that build)."""
    shared = _build_shared_locus_inputs(
        chrom=chrom,
        pos=pos,
        build_37_sumstats_task=build_37_sumstats_task,
        base_name=base_name,
        sumstats_pipe=sumstats_pipe,
        sample_size=sample_size_or_effect_sample_size,
        gene_info_task=gene_info_task,
        q_factor=q_factor,
        chrom_range=chrom_range,
        palindrome_strategy=palindrome_strategy,
        secondary_position_from_snpid=secondary_position_from_snpid,
        # Fixed to hg19: this generator runs on build-37 sumstats against the
        # Broad build-37 LD panel, the hg19 genetic map, and the hg19 baseline-LF
        # annotations. It drives the plot's x-axis coordinate-system label.
        genome_build="19",
    )
    groups = build_explainability_groups(shared)
    upset_all, upset_cs50 = build_polyfun_upset_tasks(shared.base_name, groups)
    return PolyfunExplainOuterGroup(
        groups=groups,
        upset_all_polyfun=upset_all,
        upset_cs50_polyfun=upset_cs50,
    )
