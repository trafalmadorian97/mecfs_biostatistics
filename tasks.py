"""
Repo-level commands.  Run with

pixi r invoke <task>

Note that when running tasks this way, underscores (_) in task names should be replaced with dashes (-)
"""

import datetime as _dt
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml
from invoke import task

from mecfs_bio.util.link_check_streak import (
    dump_state,
    load_state,
    read_lychee_report,
    update_streaks,
)

NEW_UNIT_TEST_PATH = Path("test_mecfs_bio/unit")
SRC_PATH = Path("mecfs_bio")
DOCS_PATH = Path("docs")

USER_AGENT = '"Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"'

PULL_FIGURE_SCRIPT_PATH = Path("mecfs_bio/figures/key_scripts/pull_figures.py")
PUBLISH_FIGURES_SCRIPT_PATH = Path("mecfs_bio/figures/key_scripts/publish_figures.py")

FIGS_PATH = Path("docs/_figs/")
FIGS_PATTERN = "_figs"


GITHUB_TOKEN_CONFIG = Path("gh_token_config.yaml")


# Default location for the streak-state file. The CI workflow downloads this
# from the prior run's artifact and uploads the updated copy after each run.
LINK_CHECK_STATE_PATH = Path("link_check_state.json")
LINK_CHECK_REPORT_PATH = Path("lychee_report.json")
LINK_CHECK_DEFAULT_THRESHOLD = 7

# Lychee flags shared by every link-checking task. Keep this in one place so
# the strict and streak-aware tasks stay in sync.
#   - 403 in --accept tolerates anti-bot systems.
#   - --cache-exclude-status 400..=999 means failure results are NOT cached,
#     so a previously-failing URL is always re-checked next run (essential
#     for streak recovery detection).
LYCHEE_SHARED_ARGS = (
    f"--insecure --cache "
    f"--accept 100..=103,200..=299,403 "
    f"--cache-exclude-status 400..=999 "
    f"--exclude {FIGS_PATTERN} "
    f"--user-agent {USER_AGENT}"
)


# helper


def _set_gh_token() -> None:
    if not GITHUB_TOKEN_CONFIG.exists():
        print(f"No github token file found at {GITHUB_TOKEN_CONFIG}")
        return
    with open(GITHUB_TOKEN_CONFIG) as infile:
        loaded = yaml.load(infile, Loader=yaml.FullLoader)
    assert "GH_TOKEN" in loaded
    print("setting github token environment variable")
    os.environ["GH_TOKEN"] = loaded["GH_TOKEN"]


def _set_use_gh_token() -> None:
    os.environ["HAVE_GH_TOKEN"] = "True"


def _set_enable_api_autonav() -> None:
    os.environ["ENABLE_API_AUTONAV"] = "True"


# dev tasks
@task
def test(c):
    print(
        "Running unit and integration tests with pytest, while skipping un-needed tests with testmon..."
    )
    cmd = f" python   -m pytest  --testmon --typeguard-packages={SRC_PATH}  {NEW_UNIT_TEST_PATH}"
    print(cmd)
    c.run(cmd, pty=True)


@task
def test_debug(c):
    """
    Run tests with extra information printed.  Useful for debugging.
    """
    print("Running unit and integration tests with pytest...")
    cmd = f" python   -m pytest -s --typeguard-debug-instrumentation --typeguard-packages={SRC_PATH}  {NEW_UNIT_TEST_PATH}"
    print(cmd)
    c.run(cmd, pty=True)


@task
def format(c):
    """
    Format code
    """
    print("Formatting with ruff...")
    c.run(" ruff format", pty=True)


@task
def formatcheck(c):
    """
    Check for format errors
    """
    print("Checking format with black...")
    c.run(" ruff format --check .")


@task
def lintfix(c):
    print("linting and applying lint auto-fixes using ruff...")
    c.run("  ruff check --fix --unsafe-fixes")


@task
def lintcheck(c):
    print("linting using ruff...")
    c.run("ruff check")


@task
def typecheck(c):
    """
    Check for type errors
    """
    print("Typechecking with ty...")
    c.run(f"ty check {SRC_PATH} {NEW_UNIT_TEST_PATH}", pty=True)


@task
def spellcheck_docs(c):
    """
    check docs for spelling errors
    Edit _typos.toml to add exceptions
    """
    print(
        f"Checking documentation spelling using typos... (Add exceptions to {DOCS_PATH}/_typos.toml)"
    )
    c.run(f"typos {DOCS_PATH}", pty=True)


@task
def spellcheck_src(c):
    """
    check code for spelling errors
    Edit _typos.toml to add exceptions
    """
    print(
        f"Checking source spelling using typos... (Add exceptions to {SRC_PATH}/_typos.toml)"
    )
    c.run(f"typos {SRC_PATH}", pty=True)


@task
def checkimports(c):
    """
    Use import linter to enforce architectural constraints
    """
    print("Checking architectural constraints using import-linter...")
    c.run("lint-imports", pty=True)


def _missing_init_dirs() -> list[Path]:
    return [
        d
        for d in sorted(SRC_PATH.rglob("*"))
        if d.is_dir() and d.name != "__pycache__" and not (d / "__init__.py").exists()
    ]


@task
def check_init_files(c):
    """
    Verify that every directory under mecfs_bio/ contains an __init__.py.
    """
    print(f"Checking that every directory under {SRC_PATH}/ has an __init__.py...")
    missing = _missing_init_dirs()
    if missing:
        print("ERROR: The following directories are missing __init__.py:")
        for d in missing:
            print(f"  {d}")
        sys.exit(1)
    else:
        print("OK: all directories have __init__.py.")


@task
def fix_init_files(c):
    """
    Add a blank __init__.py to any directory under mecfs_bio/ that is missing one.
    """
    print(f"Adding missing __init__.py files under {SRC_PATH}/...")
    missing = _missing_init_dirs()
    if missing:
        for d in missing:
            (d / "__init__.py").touch()
            print(f"  created {d / '__init__.py'}")
        print(f"Created {len(missing)} file(s).")
    else:
        print("Nothing to fix.")


@task
def check_all_links(c):
    """
    Check all links with lychee.

    Shared flags live in ``LYCHEE_SHARED_ARGS``. 403 is accepted to tolerate
    anti-bot systems.
    """
    print("Checking links with lychee...")
    c.run(f"lychee {LYCHEE_SHARED_ARGS} {SRC_PATH} {DOCS_PATH}")


@task(help={"threshold": "Consecutive failing runs required to fail the task."})
def check_all_links_with_streak(c, threshold=LINK_CHECK_DEFAULT_THRESHOLD):
    """
    Run lychee and update the per-URL consecutive-failure streak state.

    Behaves like ``check-all-links`` but only fails when a URL has been broken
    for ``threshold`` runs in a row. Transient single-day outages no longer
    fail the daily workflow.

    State is read from / written to ``link_check_state.json`` at the repo
    root. The CI workflow round-trips this file through an artifact between
    runs (see ``.github/workflows/check_links.yml``).
    """

    threshold = int(threshold)
    state_path = LINK_CHECK_STATE_PATH
    report_path = LINK_CHECK_REPORT_PATH

    print(
        f"Checking links with lychee (streak-aware, threshold={threshold} consecutive runs)..."
    )

    # --verbose is required for lychee to populate success_map in the JSON
    # report; without it we cannot tell ok URLs from URLs that simply weren't
    # seen this run. --no-progress keeps the CI log readable.
    cmd = (
        f"lychee {LYCHEE_SHARED_ARGS} "
        f"--format json --verbose --no-progress "
        f"--output {report_path} "
        f"{SRC_PATH} {DOCS_PATH}"
    )
    # lychee exits non-zero when any URL fails; we want the report regardless.
    c.run(cmd, warn=True)

    if not report_path.exists():
        print(f"ERROR: lychee did not produce a JSON report at {report_path}.")
        sys.exit(1)

    # read_lychee_report raises LycheeReportSchemaError on schema drift. We
    # let that propagate so a lychee upgrade that changes the report shape
    # surfaces as a loud CI failure rather than a silently-empty check.
    report = read_lychee_report(report_path)
    previous_state = load_state(state_path)
    result = update_streaks(
        report,
        previous_state,
        run_date=_dt.date.today(),
        threshold=threshold,
    )
    dump_state(state_path, result.new_state)

    print(
        f"\nlychee summary: total={report.total}, errors={report.errors}, "
        f"timeouts={report.timeouts}"
    )
    print(f"State written to {state_path} ({len(result.new_state)} URLs tracked).")

    if result.cleared:
        print(f"\nRecovered URLs (streak reset to 0): {len(result.cleared)}")
        for url in result.cleared:
            print(f"  {url}")

    if result.new_failures:
        print(f"\nNew failures this run (streak=1): {len(result.new_failures)}")
        for url in result.new_failures:
            print(f"  {url}")

    if result.all_failing:
        print(f"\nAll currently failing URLs ({len(result.all_failing)}):")
        for entry in result.all_failing:
            print(f"  streak={entry.streak:>3}  {entry.url}")

    if result.persistent_failures:
        print(
            f"\nERROR: {len(result.persistent_failures)} URL(s) have been failing for "
            f">= {threshold} consecutive runs and are likely permanently broken:"
        )
        for entry in result.persistent_failures:
            print(f"  streak={entry.streak:>3}  {entry.url}")
        sys.exit(1)

    print(
        f"\nOK: no URLs have crossed the {threshold}-run persistent-failure threshold."
    )


@task
def check_table_trailing_newlines(c):
    """
    Added by Claude:
    Check that markdown docs ending with a table have at least one blank line after the table.

    The mkdocs-bibtex plugin breaks bibliography footnotes when a page ends with a
    markdown table and no trailing blank lines.
    """
    print(
        "Checking for markdown files ending with a table but missing trailing newlines..."
    )
    violations = []

    for md_file in sorted(DOCS_PATH.rglob("*.md")):
        content = md_file.read_text(encoding="utf-8")
        if not content:
            continue
        last_nonempty = next(
            (line for line in reversed(content.splitlines()) if line.strip()), None
        )
        if last_nonempty and last_nonempty.strip().startswith("|"):
            if not content.endswith("\n\n"):
                violations.append(md_file)

    if violations:
        print(
            "ERROR: The following markdown files end with a table but lack a trailing blank line.\n"
            "This breaks mkdocs-bibtex bibliography rendering. Add a blank line after the final table row."
        )
        for f in violations:
            print(f"  {f}")
        sys.exit(1)
    else:
        print("OK: all markdown files pass the table-trailing-newline check.")


@task
def fix_table_trailing_newlines(c):
    """
    Added by Claude:
    Autofix markdown files that end with a table but lack a trailing blank line.
    See check_table_trailing_newlines for why this is needed.
    """
    print("Fixing markdown files ending with a table but missing trailing newlines...")
    fixed = []

    for md_file in sorted(DOCS_PATH.rglob("*.md")):
        content = md_file.read_text(encoding="utf-8")
        if not content:
            continue
        last_nonempty = next(
            (line for line in reversed(content.splitlines()) if line.strip()), None
        )
        if last_nonempty and last_nonempty.strip().startswith("|"):
            if not content.endswith("\n\n"):
                md_file.write_text(content.rstrip("\n") + "\n\n", encoding="utf-8")
                fixed.append(md_file)

    if fixed:
        print(f"Fixed {len(fixed)} file(s):")
        for f in fixed:
            print(f"  {f}")
    else:
        print("No files needed fixing.")


@task
def check_local_links(c):
    """
    Check local links with lychee
    """
    print("Checking offline links with lychee...")
    cmd = f"lychee --exclude {FIGS_PATTERN}  --offline {SRC_PATH} {DOCS_PATH}"
    print(f"running {cmd}")
    c.run(cmd)


@task
def lint_actions(c):
    print("Linting github actions...")
    cmd = "actionlint"
    print(f"running {cmd}")
    c.run(cmd)


@task(
    pre=[
        lintfix,
        format,
        spellcheck_docs,
        spellcheck_src,
        check_local_links,
        fix_table_trailing_newlines,
        checkimports,
        # fix_init_files,
        typecheck,
        lint_actions,
        test,
    ]
)
def green(c):
    pass


@task
def install_r_packages(c):
    """
    Install R packages such as TwoSampleMR and LAVA
    """
    c.run("pixi r install-mr", pty=True)
    c.run("pixi r install-lava", pty=True)


### Figures and Documentation


@task()
def pfig(c):
    """
    Pull figures from github to populate the _figs directory
    """
    c.run(f"python {PULL_FIGURE_SCRIPT_PATH}")


@task()
def publish_figures(c):
    """
    Run the end-to-end figure-system workflow: generate any missing
    figures from ALL_FIGURE_TASKS, prune orphan manifest entries (after
    confirming no docs still reference them), then update the manifest
    and push new blobs to the GitHub release.
    """
    c.run(f"python {PUBLISH_FIGURES_SCRIPT_PATH}", pty=True)


@task
def serve_docs(
    c,
    strict: bool = True,
    include_authors: bool = False,
    enable_api_autonav: bool = False,
    live_reload:bool=False
):
    """
    Use mkdocs to serve documentation
    """
    if include_authors:
        _set_gh_token()
        _set_use_gh_token()
    if enable_api_autonav:
        _set_enable_api_autonav()
    cmd = " mkdocs serve "
    if strict:
        cmd += " --strict"
    if not live_reload:
        cmd += " --no-livereload"
    print("Serving documentation...")
    print(f"running {cmd}")
    c.run(cmd, pty=True)


@task(
    pre=[
        pfig,
    ]
)
def sdocs(c, include_authors: bool = False, enable_api_autonav: bool = False):
    """
    Retrieve figures, then serve docs
    """
    serve_docs(
        c, include_authors=include_authors, enable_api_autonav=enable_api_autonav
    )


@task
def build_docs(c):
    cmd = "mkdocs build --strict"
    print(f"running {cmd}")
    c.run(cmd)


### Asset store maintenance


# Metadata is deliberately not preserved: some asset stores live on mounts (for example
# Windows drives mounted in WSL) that reject utime and chmod, which makes a metadata
# preserving copy fail after the data has already been transferred.  Asset integrity is
# tracked by the build system's verifying trace, not by filesystem metadata.  --size-only
# then keeps a re-run cheap despite the missing timestamps.
_RSYNC_ARGS = (
    "-r --size-only --no-perms --no-times --no-group --no-owner "
    "--human-readable --info=progress2"
)


def _dir_size(path: Path) -> str:
    result = subprocess.run(
        ["du", "-sh", str(path)], capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        return "unknown"
    return result.stdout.split("\t")[0]


@task(
    help={
        "execute": "Actually move data.  Without this the task only reports what it would do.",
        "reverse": "Move data back from the remap roots to the default asset root.  Run this BEFORE deleting a rule from the config, since the task only knows about subtrees the config still declares.",
        "only": "Restrict the run to configured prefixes containing this substring, so a large migration can be staged one subtree at a time without editing the config.",
    }
)
def migrate_asset_store(
    c, execute: bool = False, reverse: bool = False, only: str | None = None
):
    """
    Reconcile the asset store on disk with the path_remap rules in
    default_runner_config.yaml.

    Remapping where an asset lives does not invalidate its verifying trace, because traces
    are content hashes and carry no location.  Changing the config without moving the data
    does, however, hide the existing assets from the build system, which then rebuilds them
    at the new location while the old copy keeps occupying the disk.  This task performs
    that move.

    It is a dry run by default, since a single invocation can move tens of gigabytes.  Pass
    --execute to do the work, --only to stage a large migration one subtree at a time, or
    --reverse --execute to bring a subtree back to the default root.

    A subtree that is present under BOTH roots is reported and skipped: the two copies may
    have diverged, and picking a winner is not something this task should decide.
    """
    # Imported here rather than at module scope because importing the default runner takes
    # ~1.5s and has side effects: it reads default_runner_config.yaml, logs it, and stats
    # the asset store to warn about unmigrated subtrees.  tasks.py is the entry point for
    # every repo command, so at module scope that noise and delay would land on unrelated
    # ones like `invoke format`.
    from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
    from mecfs_bio.build_system.runner.check_roots_available import (
        RemapRootUnavailableError,
        check_remap_roots_available,
    )

    asset_root = DEFAULT_RUNNER.asset_root
    rules = DEFAULT_RUNNER.path_remap
    if not rules:
        print("No path_remap rules configured; nothing to migrate.")
        return

    # Checked here as well as in the runner, and before rsync gets a chance to mkdir the
    # destination.  Copying to an absent remap root would materialize it on the local disk,
    # after which it exists and the runner's identical check passes forever, sending the
    # very assets this task moves back onto the disk being relieved.
    try:
        check_remap_roots_available(rules)
    except RemapRootUnavailableError as error:
        print(error)
        sys.exit(1)

    print(f"default asset root: {asset_root}")
    print(f"mode: {'REVERSE (remap root -> default)' if reverse else 'remap'}")
    if not execute:
        print("DRY RUN -- pass --execute to actually move data.\n")

    if only is not None:
        print(f"restricted to prefixes containing: {only}")

    planned: list[tuple[Path, Path]] = []
    for rule in rules:
        for prefix in rule.prefixes:
            if only is not None and only not in str(prefix):
                continue
            default_side = asset_root / prefix
            remap_side = rule.root / prefix
            src, dst = (
                (remap_side, default_side) if reverse else (default_side, remap_side)
            )
            if src.exists() and dst.exists():
                print(
                    f"  DIVERGED {prefix}: present under both {src} and {dst}.  "
                    f"Skipping -- resolve by hand."
                )
                continue
            if not src.exists():
                status = "already in place" if dst.exists() else "no data"
                print(f"  SKIP     {prefix}: {status}.")
                continue
            print(f"  MOVE     {prefix} ({_dir_size(src)}): {src} -> {dst}")
            planned.append((src, dst))

    if not planned:
        print("\nNothing to do.")
        return
    if not execute:
        print(f"\n{len(planned)} subtree(s) would be moved.  Re-run with --execute.")
        return

    for src, dst in planned:
        print(f"\nmoving {src} -> {dst}")
        dst.parent.mkdir(parents=True, exist_ok=True)
        c.run(f'rsync {_RSYNC_ARGS} "{src}/" "{dst}/"', pty=True)
        # Confirm the copy is complete before deleting anything.  A second pass that
        # reports no remaining transfers is the check.
        verify = c.run(
            f'rsync {_RSYNC_ARGS} --dry-run --itemize-changes "{src}/" "{dst}/"',
            hide=True,
        )
        remaining = [
            line for line in verify.stdout.splitlines() if line.startswith(("<", ">"))
        ]
        if remaining:
            print(
                f"REFUSING to delete {src}: {len(remaining)} file(s) still differ after "
                f"the copy.  First few: {remaining[:5]}"
            )
            continue
        shutil.rmtree(src)
        print(f"moved and removed source {src}")

    print("\nDone.  Re-run without --execute to confirm the store is consistent.")


# initialization
@task(pre=[install_r_packages, pfig, green])
def init(c):
    """
    Initial repo setup
    """
    pass
