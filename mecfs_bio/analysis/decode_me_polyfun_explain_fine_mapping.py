"""
Script to run the polyfun-vs-uniform explainability fine mapping on DecodeME
data, over the same loci as run_fine_mapping_decode_me_analysis.

Each locus produces an 8-run outer group (four SUSIE configs, each a matched
uniform/polyfun pair) plus a contrast and plot task per pair; terminal_tasks
covers all of them.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr1_174_128_548 import (
    POLYFUN_EXPLAIN_CHR1_174,
)


def run_polyfun_explain_fine_mapping_decode_me_analysis():
    """
    Function to run the polyfun explainability fine mapping across the main
    DecodeME loci: for each locus a matched uniform/polyfun SUSIE pair under
    each run config, with the contrast and callout plot explaining the pair.
    """
    DEFAULT_RUNNER.run(
        (
            POLYFUN_EXPLAIN_CHR1_174.terminal_tasks()
            # + POLYFUN_EXPLAIN_CHR6_26.terminal_tasks()
            # + POLYFUN_EXPLAIN_CHR6_97.terminal_tasks()
            # + POLYFUN_EXPLAIN_CHR15_54.terminal_tasks()
            # + POLYFUN_EXPLAIN_CHR17_50.terminal_tasks()
            # + POLYFUN_EXPLAIN_CHR20_47.terminal_tasks()
        ),
        incremental_save=True,
        must_rebuild_transitive=[
            group.contrast
            for chrom in [
                POLYFUN_EXPLAIN_CHR1_174,
                # POLYFUN_EXPLAIN_CHR6_26,
                # POLYFUN_EXPLAIN_CHR6_97,
                # POLYFUN_EXPLAIN_CHR15_54,
                # POLYFUN_EXPLAIN_CHR17_50,
                # POLYFUN_EXPLAIN_CHR20_47,
            ]
            for group in chrom.groups
        ],
    )


if __name__ == "__main__":
    run_polyfun_explain_fine_mapping_decode_me_analysis()
