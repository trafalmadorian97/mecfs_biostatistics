"""Demonstrator: polyfun-vs-uniform explainability at the DecodeME chr6:26.2Mb
locus.

Same locus, sumstats, sample size, palindrome strategy, and chrom range as the
with_palindromes chr6_26 fine-mapping module, so the eight SUSIE runs here
operate on the same harmonized inputs.
"""

from mecfs_bio.asset_generator.polyfun_explain_fine_mapping_asset_generator import (
    generate_assets_polyfun_explain_fine_map,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_chrom_pos_alleles import (
    ChromRange,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    SecondaryPositionFromSnpid,
)

POLYFUN_EXPLAIN_CHR6_26 = generate_assets_polyfun_explain_fine_map(
    chrom=6,
    pos=26_215_000,
    build_37_sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.join_task,
    base_name="decode_me_polyfun_explain",
    sumstats_pipe=IdentityPipe(),
    sample_size_or_effect_sample_size=int(
        4 / (1 / 15_579 + 1 / 259_909)
    ),  # 4/(1/cases + 1/controls)
    palindrome_strategy="keep",
    chrom_range=ChromRange(6, 26_000_000, 27_000_000),
    secondary_position_from_snpid=SecondaryPositionFromSnpid(build_label="hg38"),
)
