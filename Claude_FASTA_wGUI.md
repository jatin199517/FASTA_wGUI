# FASTA Enhanced — GUI Development Context

> ## ⚠️ GUI Build — read first
>
> This file is the **GUI-development** context for FASTA Enhanced. The flow-analysis engine in **`ExtFunctions*.R`** (the base code documented below) is being used as the **base/backend for building a GUI for flow analysis**. The rest of this document (carried over from the engine `CLAUDE.md`) describes that backend so the GUI can call it correctly.
>
> **Ground rules for the GUI effort:**
> 1. **`ExtFunctions*.R` is a read-only, black-box dependency — DO NOT edit it going forward.** It may be **swapped for a newer version** at any time (the engine is developed/forked separately). Treat its functions (`apply_gating2`, `apply_compensation`, `apply_plotting`, `panel_QC`, `run_gating_qc`, the `.gate_*` plugins, etc.) as a stable API surface, not something to modify. If a behavior change is needed, it belongs in the engine project, not here.
> 2. **All GUI work lives in THIS folder** — GUI notes, design decisions, and progress go in this file (`Claude_FASTA_wGUI.md`); GUI source goes in new scripts kept in this folder.
> 3. **New GUI scripts must be independent of `ExtFunctions.R`** — they may *source/call* it as a dependency, but must not depend on edits to it and must not assume a specific engine version beyond the documented API. Keep GUI logic self-contained so a newer `ExtFunctions*.R` drops in without touching GUI code.
>
> The engine reference below is **descriptive** (how the backend behaves); it is **not** an invitation to modify the engine.

---

## Project

**FASTA Enhanced** is the shared, automated flow-cytometry analysis engine at **Lyell Immunopharma** — standardized gating, statistics, and layout plots from FCS files using R/Bioconductor (flowWorkspace, openCyto, ggcyto). One shared codebase (`ExtFunctions*.R`) is driven by a thin **orchestration script** — **now `FRISDA_0_73.R`** ("Flowcore's Really Insane Script for Deep immunophenotyping Analysis"; supersedes the legacy `FASTA_enhanced.R`), with PDF report `260608_Report_FRISDA_0.73.Rmd`. Downstream studies reuse the same engine via their own orchestration script — e.g. **IVANKA**, the in-vivo CAR PK adaptation, which is a **separate project/folder**.

**Owner:** Jatin Sharma (jsharma@lyell.com), Research & Science  
**Study context:** **LYL273** FDP Deep Phenotyping (Deep-IPT), unmixed. CAR-T IND/BLA preparation.  
**Branch:** `FASTA_0_5`.

> **Codebase round-trip (2026-06-13):** This CLAUDE.md was carried over to the **IVANKA** project so that every engine change made there to `ExtFunctions*.R` was captured in one place — **AutoSpill compensation** (`build_autospill_controls`/`run_autospill` + the `apply_compensation` `autospill` mode), the **`pop="+-"` 1D split gate**, the **`layout_dims`/`layout_gates`** plotting overrides, **FCS export args**, and the **2D-rect overlay fix**. We are now **back in FASTA Enhanced**, continuing from **IVANKA's last engine update (`ExtFunctions_0_72.R`), now bumped to `ExtFunctions_0_73.R`** here — which carries all of those features forward. The engine/compensation/plotting sections below describe that shared codebase and are valid here; **script-level wiring details that name IVANKA** (e.g. the `*_BackupCompMatrix_IVANKA.csv` fallback, the E44 QC plot) reflect IVANKA's *orchestration* and are not necessarily wired into `FASTA_enhanced.R`. IVANKA-only adaptation tasks moved with the IVANKA project (see Pending Work).

> **Working convention:** After successfully implementing or validating *any* change, **prompt the user to update this CLAUDE.md** (and the History / Pending Work sections) before moving on — don't silently skip it. CLAUDE.md is the single source of truth for this project.

> **Engine changes go in `ExtFunctions*.R` (it is forked into other projects).** ALL analysis/engine logic — gating, gate plugins, `run_gating_qc` checks, `apply_plotting` behavior, compensation, `panel_QC`, stats — belongs in `ExtFunctions*.R`, never duplicated into the orchestration script (`FRISDA_*.R`) or the Rmd report. The engine file is the shared skeleton forked into downstream projects (e.g. IVANKA), so a fix placed here propagates to every fork on the next merge; a fix placed in an Rmd/orchestration script would have to be redone per project. **Only genuinely report-/project-specific presentation or wiring stays out of the engine:** LaTeX preamble (`header-includes`), `kable`/`kable_styling` rendering options (`longtable`, `repeat_header`, caption escaping), data paths, `instructions_filename`, `PanelIteration`, and the locate-and-source bootstrap. When a change *could* live either place, prefer `ExtFunctions*.R`.

---

## Key Files

| File | Purpose |
|------|---------|
| `ExtFunctions_0_73.R` | **The entire codebase** — carried forward from IVANKA's last engine update (`0_72`), bumped to `0_73` here. Gate plugins, `apply_compensation`, **`build_autospill_controls`/`run_autospill`**, `apply_transformation`, **`apply_gating2`** (gating engine), `apply_plotting`, `panel_QC`, `run_gating_qc`, helpers |
| `FRISDA_0_73.R` | **Active main script** — bootstraps + sources `ExtFunctions*.R` via the `grepv("ExtFunctions", …)` glob (line ~50, so **no stale hardcoded re-source**), runs the **LYL273 Deep-IPT** pipeline. Orchestration only. Data dir: `…/Confidential LYL273 Translational Science…/Deep_IPT/260602_EXP26000696/…/Unmixed`. Live instructions `20260606_instructions_draft_LYL273.xlsx`; `PanelIteration <- "260527_LYL273_[JS]"` |
| `260608_Report_FRISDA_0.73.Rmd` | FRISDA PDF report (xelatex; Gating Scheme / Data Summary / Panel QC / Gating QC / Layouts / Session Info). **Note:** its setup chunk + header still grep `^FASCA_.*\.R$` and print "FASCA/FASTA" — stale naming that breaks `scriptName`/`scriptVer` extraction (no `FASCA_*.R` exists). Fix to `FRISDA` when convenient |
| `FASTA_enhanced.R` | **Legacy** main script (superseded by `FRISDA_0_73.R`). Still has the stale `source("ExtFunctions_0_70.R")` line — not used going forward |
| `20260606_instructions_draft_LYL273.xlsx` | **Live** instructions workbook (`instructions_filename` in the script): `metadata`, `markers`, `gating_template` sheets → `finalInstructions`. (Also in-folder: `20260531_…_BCA_panel.xlsx`, `20260604_instructions_cPARP_0.5.xlsx` — other panels, not the live one.) |
| `GatingTemplate.csv` | Legacy minimal draft template — vestigial (the live template is the xlsx `gating_template` sheet) |
| `Older versions/` | Prior-version backups (`0_5`, `0_59`, `0_60`, `0_63`, `0_69_wgtgating`, `0_70`, `0_71`). Safe because the ExtFunctions glob is **non-recursive** (top level only) |
| `validate_plot_stats.R` | Plot/stats validation helper |

**Filename is versioned** (`ExtFunctions_0_73.R`). The script locates it via `grepv("ExtFunctions", list.files(pattern=".*\\.R$"))` — a **non-recursive** glob, so keep **exactly one** `ExtFunctions*.R` at the **top level** (backups in the `Older versions/` subfolder are not matched and so are not sourced — fine; but two at top level would both source and silently double-define everything).

**Re-load gotchas (bit us repeatedly):**
- After editing `ExtFunctions_0_73.R`, you **must re-`source()`** it — re-running a pipeline call uses the stale in-memory function otherwise.
- After editing the xlsx, you **must re-read** `finalInstructions <- setNames(map(listOfSheets, ~read_excel(...)), listOfSheets)` — and re-execute downstream assignments (e.g. `PanelIteration`).
- A **gating** fix needs `gs <- apply_gating2(...)` re-run (re-sourcing alone doesn't recompute an existing `gs`); a **plotting-only** fix needs only `apply_plotting(gs, ...)`.

---

## Gating engine: `apply_gating2` (gt_gating-free)

**As of 2026-06-11, `apply_gating2` replaced the old `apply_gating`.** The old function and ALL its openCyto `gt_gating` machinery (`gt_gating`/`run_row`, `precompute_external_control`, the `.deferred_*` deferral logic, the `.logic_indices` restoration inside it) were **removed**. `apply_gating` no longer exists.

`apply_gating2(gs, template, mc.cores = 1, verbose = TRUE)` — sequential driver:
- Walks template rows **top-to-bottom**. Parents always exist before children → **no DAG resolution, no precompute cache, no deferrals**.
- Per row: resolve parent + channels → pick source data → compute a per-sample gate via `.fasta_compute_gate` → add via `gs_pop_add`.
- Clears `.fasta_gate_cache` at start so plugins compute fresh.

**`.fasta_compute_gate(fr, method, channels, filterId, args)`** — the one dispatcher: calls `.gate_<method>` (our plugins) or a builtin (`singletGate`, `gate_flowclust_2d`) **in `per_sample` mode** on whatever flowFrame the driver feeds it. (So the driver, not the plugin, decides which data/control to use.)

**Per gate type, how the driver adds it:**
- **1D rect** (mindensity / FMO / flowmeans / static): named per-sample list → `gs_pop_add`. `pop="-"` → flip gate to `[-Inf, thr]`. **`pop="+-"` or `"-+"` (equivalent) → 1D split:** from the one shared per-sample threshold, add BOTH `[thr, Inf]` as `<alias>+` and `[-Inf, thr]` as `<alias>-` (so set `alias` to the base name). Strictly 1D — a non-1D method with `pop="+-"` hard-errors. In `apply_plotting`, a `+-` row auto-overlays + labels **both** child pops (shared threshold line) via the same `.lg_override` path as `layout_gates`.
- **2D quad** (`gate_flowmeans_2d` / `gate_FMO_2d` / `gate_mindensity_2d`): plugin returns a `quadGate` → **decomposed into 4 `rectangleGate`s** (sidesteps the `gs_pop_add` quadGate per-sample collapse; counts identical to the old quadGate→polygon-children layout).
- **`gate_tmix_2d`**: returns a `filters` of 4 `polygonGate`s → added per sample as-is (curved contours preserved).
- **`gate_logic`**: **index-based** (NOT booleanFilter — see below). Dummy pass-through gate + per-sample AND/OR/NOT via `gh_pop_set_indices` (`NOT = parent_idx & !inner`). A gate_logic referencing **another** gate_logic restores the prior nodes' indices right after its own `recompute()` (else the dependency reads the wiped 100% dummy — this was the `CAR-` bug). Indices re-applied after the loop, and gate_logic **descendants are re-gated** against the restored parent (supports e.g. a `copy_gate` under `CAR_all`).
- **`copy_gate`**: copies gate(s) to a new parent — single-node sources, 2D-quad-children (tmix polygons or rect quads, matched via `name_x/name_y`), quadGate-node; **rejects** `gate_logic` sources (no geometry); skips inactive sources.

**Control modes:** `per_sample` (own data); `external_control` (driver feeds the **control sample's** data to the plugin — single `control_col/val`, per-axis `control_col_x/y`, and/or `groupBy` per-group control); `clusterFrom` (compute from the named population's data).

`.FASTA_QUAD_2D_METHODS = c("gate_flowmeans_2d","gate_FMO_2d","gate_tmix_2d","gate_mindensity_2d")` — used by `apply_plotting`/`run_gating_qc` for quad detection. (`apply_gating2` dispatches by the *returned gate type*, not this list.)

### Why `gate_logic` is index-based, not `booleanFilter`
`booleanFilter` is the "native" gt_gating-free path, but in this flowWorkspace build **its references never resolve at `recompute()`** (tried bare name, full path, root-relative — fail at every depth), **and every failed boolean recompute corrupts the shared cytoframe backing** (propagates back to `cs_comp`, needs an FCS reload to recover). The index approach (`gh_pop_set_indices`) is the proven mechanism the old engine used and never corrupts on failure. Don't reach for `booleanFilter` here.

---

## Custom Gating Methods (plugins in `ExtFunctions_0_73.R`)

| Method | Description |
|--------|-------------|
| `gate_FMO` | 1D FMO rectangleGate |
| `gate_FMO_2d` | 2D FMO quadGate. Per-axis controls via `control_col_x/y`, `control_val_x/y` |
| `gate_flowmeans` | 1D flowMeans gate |
| `gate_flowmeans_2d` | 2D flowMeans quadGate, per-axis control |
| `gate_mindensity` | 1D min-density (wraps openCyto). `upper_mult`/preset args honored both paths |
| `gate_mindensity_2d` | 2D min-density quadGate (in `.FASTA_QUAD_2D_METHODS`; `name_x/name_y`, combined-key precompute caching) |
| `gate_tmix_2d` | 2D quad wrapping openCyto `gate_quad_tmix` (joint t-mixture). Returns 4 `polygonGate`s; rectangular crosshair when well-resolved is expected openCyto behavior, not a bug. Args: `K`, `quantile1/3` (def 0.8), `usePrior`, `prior`, `trans`, `name_x/y`, `clusterFrom`, `groupBy` |
| `gate_singlet` | Singlet gate, optional external_control |
| `gate_static` | Fixed rectangleGate from `xmin/xmax/ymin/ymax` (or `min/max` 1D) |
| `gate_logic` | Boolean **AND / OR / NOT** of existing populations. `gates=` accepts comma-string `'A,B'` or `c('A','B')`. Index-based (see engine note) |
| `copy_gate` | Copies gate(s) from another population (incl. tmix/quad children, and under a gate_logic parent) |

**gt_gating-era leftovers:** the `register_plugins()` calls and the `add_boolean` helper were **removed (2026-06-11)**. The plugins' `method="external_control"` branches were **kept by decision** — inert for `apply_gating2` (always `per_sample`) but keep the `.gate_*` functions reusable as standalone dual-mode tools. `.fasta_gate_cache` is kept (plugins still write `quad_names` to it; `apply_gating2` clears it at start).

---

## Compensation (`apply_compensation`)

`apply_compensation(cs, metadata, autospill_matrix = NULL, export_fcs = FALSE, append_pdata_to_export = FALSE)` dispatches **per sample** on the metadata `Compensation_type` column:
- **`self`** — apply the acquisition-embedded `$SPILLOVER`/`$SPILL` matrix.
- **`unmixed` / NA / blank** — no compensation (data already unmixed/compensated).
- **`<file>.csv` / `<file>.xlsx`** — a square spillover matrix file in the working dir, applied via `compensate()`; files are loaded once and cached if reused.
- **`autospill`** — applies the **experiment-level** AutoSpill matrix passed via `autospill_matrix` (one matrix shared by all `autospill` samples; errors if `NULL`). All-`autospill` is fast-pathed via `compensate(cs, autospill_matrix)`. `apply_compensation` only **applies** the matrix — it never runs AutoSpill or touches files.
- **`autospectral`** — still a `stop()` **stub**.

**FCS export args** (one internal `do_export()` wraps every return path, so export is consistent across all modes):
- **`export_fcs`** (default FALSE) — when TRUE, writes each compensated frame to a **`Compensated/`** folder in `getwd()`, named `<original_basename>_<comptype>.fcs` (`<comptype>` = the sample's resolved `Compensation_type`, sanitized — e.g. `C1 117_autospill.fcs`, a `spill.csv` type → `_spill`). Exports whatever sample set was passed in (so a `cs[StainType != "Single Stain"]` subset writes just those).
- **`append_pdata_to_export`** (default FALSE) — when TRUE **and** `export_fcs` is TRUE, that sample's `pData` row is appended to the exported frame's FCS keywords (original acquisition keywords **+** one keyword per pData column; `NA`→`""`). The frame is detached via `cytoframe_to_flowFrame` first, so `cs_comp` itself is **not** mutated. No effect when `export_fcs` is FALSE.
- **Note:** set-level `pData(cs)` is *not* otherwise embedded in frames — plain `export_fcs` carries only the original acquisition keywords; use `append_pdata_to_export` to bake the layout metadata in (the old `FASTA_v0.851` did this via `annotate_changeParam_pdWrite_cs`).

### AutoSpill (`build_autospill_controls` + `run_autospill`, implemented 2026-06-12)
AutoSpill (`carlosproca/autospill` v0.2.0, GitHub-only) computes **one experiment-level** spillover matrix from single-color controls. Two engine helpers stage and compute it; the **script** wires them and writes the backup CSV (kept out of `apply_compensation`):
- **`build_autospill_controls(markers, controls, control_dir, ...)`** — the dye↔antigen pairing is **VETTED from the markers sheet** (`dye=dims`, `antigen=markers`); the method does **not** guess channels (an empirical max-signal peak is **unreliable** for spectral neighbours, e.g. `EF780`/`BV786`, `RB705`/`RB780`). Each control's **filename is auto-matched** from `controls` (a `filename`+`SampleID` df) by normalised fluorophore, with `alias = c(PB→"Pacific Blue", eFluor780→"Viability")` for the two SampleID tokens that diverge. **Hard-errors** on any unmatched vetted row.
  - **Beads used as-is:** SpectraComp/ViaComp beads each carry their own positive+negative populations in one file, so they're **symlinked unchanged** — AutoSpill gates pos/neg within each. No negative concatenation (AutoSpill has no native separate-negative input anyway). The `US SC`/`US VC` unstained files are therefore **not used**.
  - **GFP exception:** two files (`GFP+`/`GFP-`) → concatenated to `"GFP single stain pos_neg.fcs"` (so AutoSpill sees both pops); a **single** GFP file → assumed already pos/neg mixed, used as-is.
  - Symlinks + (optional) GFP concat go to a **temp dir** (nothing written to the data folder). **No `AutoF` row**. Prints the mapping for sign-off; returns `list(control_dir, control_def, mapping)` where **`control_def` is an in-memory data.frame** (the `fcs_control.csv` equivalent — not written to disk).
- **`run_autospill(control_dir, control_def, scatter_param = c("FS00-A","SS02-A"), param_set = "minimal")`** — wraps the AutoSpill `minimal` workflow (`read.flow.control`→`gate.flow.data`→`get.marker.spillover`→`refine.spillover`). **`minimal` writes nothing to disk.** `control_def` may be the **data.frame** (preferred) or a CSV path; AutoSpill's readers hardcode `read.csv()` on a path **and read it twice** (so a `textConnection` is exhausted and a dataframe can't be passed directly without forking the package) — `run_autospill` therefore serialises the dataframe to a **throwaway `tempfile()`** for the call only. Errors (not converges) propagate so the script can fall back. Returns `refined$spillover` (square, diagonal = 1, named by detector channel).
- **No native per-marker negative path:** `get.marker.spillover.posnegpop`/`process.posnegpop` derive pos/neg by **gating within** each control file (`flow.gate, flow.control, asp` only) — there is no separate-negative-file column. Hence the concat approach.
- **Scatter gotcha:** AutoSpill reads the **raw** FCS (pre-rename) and defaults scatter to `c("FSC-A","SSC-A")`; this instrument names scatter `FS00-A`/`SS02-A` (the orchestration renames `FS00→FSC`/`SS02→SSC` on the `cs` object only, via markers `dims`, line ~112) — so `run_autospill` overrides `asp$default.scatter.parameter` to the **raw** names, else `read.flow.control` aborts.
- **Channel names safe:** `-` is not a forbidden char, so `FL25-A`-style names are not mangled → matrix row/col names match `cs` channels for `compensate()` directly.
- **Script pattern + fallback:** `build_autospill_controls(markers, controls, control_dir)` → `compMatrix <- tryCatch(run_autospill(...), error → read latest local *_BackupCompMatrix_IVANKA.csv` and set `CompType <- "Local matrix")` (default `CompType <- "Autospill"`; the error is messaged for debugging) → `apply_compensation(cs, metadata, autospill_matrix = compMatrix)` → **on success only**, write `<YYYYMMDD>_<BenchlingID>_BackupCompMatrix_IVANKA.csv` (a fallback run does not rewrite the backup with itself) → **compensation QC plot**: the spike sample, CD3 (x) vs CD8 (y), faceted Uncompensated (`cs`) vs Compensated (`cs_comp`), via `flow_bin2d` on `asinh(value/cofac)` (cofac=150 display transform, since `cs`/`cs_comp` are untransformed; channels found by `grepl("CD3"/"CD8", markernames(cs))`).
- **Panel-resolution note:** AutoSpill can emit off-diagonal spillover **>1** for near-identical dyes (here `EF780`↔`BV786` etc.) — expected behaviour, a panel-design signal, not a bug.

**`autospectral`** is still an unimplemented stub.

---

## Gating Template Columns

Standard openCyto columns plus custom: `collapseDataForGating`, `groupBy`, `active_gate` (TRUE/FALSE — skip if FALSE), `layout_dims`, `layout_gates`, `layout_row`, `layout_col`, `layout_index`, `PlotType` ("Scatter"), `label_col`, `layout_labels` (which of gate_name/percent/count to show).

**`collapseDataForGating` / `groupBy` (wired into `apply_gating2` 2026-06-13):**
- **`collapseDataForGating`** (TRUE/FALSE/blank) — when TRUE, POOL the parent (or `clusterFrom`) events across the group, compute ONE gate via the row's method, and replicate it to every sample in the group. Reproducibly downsampled (per-sample cap = `collapse_max_events/n`, fixed `collapse_seed`; pools only the gating channels). **Mutually exclusive with `external_control`** (errors). Excluded from the identical-threshold QC checks (identical by design).
- **`groupBy`** — grouping for `collapseDataForGating` (and `external_control`'s per-group control): a pData column (or `":"/","` combo, e.g. `"Batch:Day"`), or a number `N` to group every N samples; blank = all samples in one group.
- `collapse_max_events` (def `5e5`) and `collapse_seed` (def `1234`) are **`apply_gating2()` arguments, NOT template `gating_args`** — set them in the orchestration call (`FRISDA_*.R`), not the template; they are not subject to the unrecognized-arg guard.

**`layout_dims` / `layout_gates` (wired into `apply_plotting` 2026-06-13):**
- **`layout_dims`** — comma-separated **markernames** (resolved to channels via `resolve_ch`); when non-blank, overrides the panel axes for that row (instead of the gating `dims`). Forces a `computePlotParams` recompute (bypasses the alias-keyed cache, which holds the default-dims limits/binwidth). Unresolvable marker → hard error.
- **`layout_gates`** — comma-separated **gate/population names** (resolved to full paths via `resolve_pop`, so no full path needed); when non-blank, **replaces** the overlay set — those gates' boundaries are drawn **only where their dims match the panel axes** (a non-matching listed gate is silently skipped, not an error), and the default primary-gate annotation + gate_logic crosshair are suppressed so ONLY the listed gates show. Unresolvable gate → hard error. **Each listed gate also gets its own `layout_labels` stats label** (gate_name/percent/count) — including on a `gate_logic` row, where the per-gate label is normally skipped (the overlay-loop label guard is relaxed under `layout_gates`). E.g. a `notBeads` panel with `layout_gates="Beads, notBeads"` shows the Beads box + `Beads %/count` label *and* the `notBeads %/count` centroid label.
- Both blocks sit **after** the final `overlay_pops`/channel computation (an earlier first attempt was silently clobbered by a later `overlay_pops <-` reassignment — watch for that).

**Layout QC:** `apply_plotting` stops with a descriptive error if two gates share the same `layout_index + layout_row + layout_col`.

---

## panel_QC

`panel_QC(gs, PanelMasterInventoryLocation, PanelIteration)` — checks GatingSet markernames/channels against a master panel inventory **iteration**. Returns a data.frame of mismatches (empty = clear).

- **`PanelIteration` MUST match the data's assay panel.** A wrong iteration flags everything. Deep-IPT data → e.g. `"260527_LYL273_[JS]"`, **not** `"LYL273 cPARP_1"` (the cPARP functional panel).
- **Antigen matching** is separator-insensitive via `.norm_ag` (strips `_ - . whitespace`): `PD-1`==`PD1`, `HLA-DR`==`HLADR`, `CD14_CD16`==`CD14CD16`. `/` and parens are **not** stripped → `TCF7/1`, `CD197(CCR7)` need a `Marker_Alias`.
- **Fluorochrome matching** via `.norm_fluor`: strips `-A`/`-H` channel suffix, spaces/separators, and maps `AF###`↔`Alexa Fluor ###`.
- `.clean` strips zero-width/invisible Unicode (the inventory is full of them).
- **Viability dyes are exempt both directions** (panel "Live Dead eFluor 780" ↔ GS marker "Viability").
- `Marker_Alias` column (inventory) is preferred over `Antigen` for matching — fill it for alias gaps.
- A **panel marker not run in a sample** (e.g. `CD45`) correctly flags — that's a true positive, not noise.

---

## Plotting Standards

**Always use `flow_bin2d()` and `computePlotParams()` — never `geom_hex()` or `ggcyto_par_set()`.** Panels are plain `ggplot` on extracted `exprs()` (columns renamed `x_var`/`y_var` for faceting) — *not* `ggcyto` objects, so `geom_gate()` does not apply here.

```r
params <- computePlotParams(gs, pop, xchannel, ychannel, adaptiveBins(gs, pop))
p <- ggplot(dat, aes(x=x_var, y=y_var)) +
  flow_bin2d(binwidth = params$binwidth) +
  coord_cartesian(xlim = params$xlim, ylim = params$ylim, clip = "on")
```

Auto gate annotation — note `apply_plotting` annotates the primary gate in **two passes** (the `.g_auto` primary-population block *and* the `overlay_pops` loop); both must stay consistent (they currently both draw every gate — overlapping identical draws). Each routes the `rectangleGate` case by how many plotted axes the gate constrains:
- **1D rect** (gate constrains only one of the plotted axes — e.g. `Live` = `[-Inf,thr]` on Viability) → a single **vline/hline** at the finite boundary (min **or** max, so `pop="-"` draws a line, not a box).
- **2D rect** (gate constrains **both** plotted axes — e.g. `notDebris` `[50000,Inf]×[21000,Inf]`) → draws **only the finite edges as `geom_segment`s** (open `±Inf` edges omitted entirely; open ends extended well past the view so `coord` clip trims them with no perpendicular cap) — i.e. the L-shaped corner. Avoids both the old full-panel **crosshair** (single-vline+single-hline) *and* a `geom_rect` box, whose top/right edges land at the panel boundary because `coord_cartesian`'s default expansion sits just inside the view.
- **quadGate children** → crosshairs.
- **tmix polygon quads** → polygon outline; labels at **quadrant centers** (a label-crosshair derived from each polygon's center-facing finite corner; falls back to event centroid). **No threshold VALUE label** — a polygon has no single axis-aligned threshold to print (by design, not a bug); the red "crosshair" seen on a well-resolved tmix is just the 4 polygon outlines meeting. Same for a `copy_gate` of a tmix.
- **gate_logic** → **no box** (the all-infinite dummy gate is skipped); a **crosshair derived from the component quadrant thresholds** (recursing one level into a referenced gate_logic, so `NOT(CAR_all)` finds the CAR thresholds); a **label at the population's event centroid**.

**Threshold VALUE labels (the red numbers): bake position into the data.frame, never compute it inside `aes()`.** `apply_plotting` builds ggplots in a loop and they're printed LATER (in the Rmd); an `aes()` expression referencing a loop-local like `.plot_xlim`/`.plot_ylim` is resolved at print time against the loop's FINAL value → label lands offscreen. Compute positions immediately (in the `data.frame(...)` call or as constants outside `aes`), map only baked columns in `aes` (`aes(y = yv, label = lbl)`). This was the 2026-06-14 y-label bug.

---

## Per-Sample Gate Application Pattern

`gs_pop_add(gs, namedListBySample, parent)` applies per-sample for rectangleGate/polygonGate. **quadGates collapse** (C++ `cpp_addGate` — last sample wins), which is why `apply_gating2` decomposes 2D quads into rectangleGates. For a manual per-sample fix: `for (sn in sampleNames(gs)) gh_pop_set_gate(gs[[sn]], node, per_sn_gate); recompute(gs)`.

---

## History / notable fixes

- **2026-06-14 — quad y-border value labels fixed (lazy-`aes` loop-var bug) + collapse caption.** The y-value threshold labels on quad gates were invisible because their position was computed **inside `aes()`** (`aes(y = pmin(yint + diff(.plot_ylim)..., .plot_ylim[2]...))`) — `aes()` resolves that lazily at **print** time, and `.plot_ylim` is an `apply_plotting` **loop variable** that, by the time the stored plots print (in the Rmd, after the loop), holds the LAST panel's value → label lands offscreen. The x-labels were accidentally immune (position baked into a data.frame column + constant `y=` outside `aes`). **Fix:** all 5 y-value `geom_text` now bake position + label into the data.frame (`data.frame(yv = pmin(<thr>+…, …), lbl = round(<thr>), …)`, `aes(y = yv, label = lbl)`). Swept the whole file — these were the only lazy-`.plot_*lim`-in-`aes` instances. Validated on LYL273: DUMP/CD3_CD19/CAR now show both x and y threshold values. **`gate_tmix_2d` (and its `copy_gate`) intentionally show NO threshold value label** — polygon quads have no single axis-aligned threshold (the red crosshair there is the 4 polygon outlines); confirmed by-design, not a bug. Also: the **plot caption now shows `collapseDataForGating: TRUE`** (next to Args/groupBy) when set. *Lesson: never reference loop-local vars inside `aes()` when ggplots are built in a loop and printed later — bake them into the data/constants (immediate eval).*

- **2026-06-13 — `collapseDataForGating` implemented + more QC/report fixes.** (a) **`collapseDataForGating`** column now honored by `apply_gating2` (mirrors openCyto): pool the parent/`clusterFrom` events across the group (groupBy: blank→all, numeric N→every N, pData col(s) `":"/","`→unique combo), compute ONE gate via any method, replicate to the group. **Reproducibly downsampled** — per-sample cap = `collapse_max_events/n` (new `apply_gating2` arg, default 5e5), fixed `collapse_seed` (default 1234), global RNG saved/restored so seeding is local. Pools **only the gating channels** (memory). **Mutually exclusive with external_control** (errors). Skipped in QC identical-threshold checks (1 & 2 — identical thresholds are by design). `collapse_max_events`/`collapse_seed` are **driver args, NOT gating_args** (don't go in the template, not subject to the arg guard). (b) **`threshold_outside_data_range` → Option A**: restricted to `per_sample` gates (external/FMO thresholds are control-derived, so checking against each sample's own P1–P99.95 is apples-to-oranges noise). (c) **Gating QC tables paginate** — `longtable=TRUE`+`repeat_header` in the Rmd + `\usepackage{longtable}` in the preamble (was `longtable=FALSE`+`scale_down`, which `\resizebox`'d onto one page and clipped long tables). (d) **Quad x-border value labels** fixed (`hjust=1` so the rotated label hangs into the panel, not clipped at the top). *y-border value labels are still missing — see Pending #4; the move-inside+clamp attempt did NOT fix it, so it's NOT edge-clipping.*

- **2026-06-13 — first-PDF review fixes (gating QC, mindensity arg guard, plot labels).** Reviewing the first FRISDA PDF surfaced several issues, fixed in `ExtFunctions_0_73.R`:
  - **`gate_logic_count_mismatch` false-flagged every `NOT` gate** — the check computed `expected` as `Reduce("&")` (= the component count) and compared it to the stored `NOT` count (`parent & !component`). Now handles `NOT` explicitly (reads parent indices). 
  - **`gate_logic_component_missing` false-flagged every comma-string `gates=`** — Check 9 used `lg_args$gates` raw, so `gates="A,B,C"` was treated as ONE bogus pop (the reported flag literally listed `GCC_CAR+CD19_CAR+,GCC_CAR+CD19_CAR-,GCC_CAR-CD19_CAR+` as one missing component). Now parsed via `.fasta_parse_pop_list()` like `apply_gating2` (this also un-blocked the OR count-check, which had `valid_comp` empty). 
  - **`threshold_outside_data_range`** band widened P1–P99 → **P1–P99.95** (`c(0.01, 0.9995)`) to cut noise from legitimately-high quantile/FMO gates.
  - **Unrecognized `gating_args` now error loudly** — `.fasta_compute_gate` compares each `.gate_*` plugin's `plugin_args` against `formals(fn)` and `stop()`s naming the offending arg (e.g. `x_min` instead of `min_x`, the DUMP bug). **`gate_static` is exempt** (it reads its bounds from `...` by design). Closes the silent-arg-drop class (`x_min`/`y_min`, `y_min, 5000`) for good.
  - **Quad x-border value labels were clipped** — rotated (`angle=90`) value labels anchored at `ylim[2]−5%` with `hjust=0.5` extended above the panel and were cut by `coord(clip="on")`. Changed all 5 rotated x-value `geom_text` to **`hjust=1`** (hang down into panel, axis-label convention). *y-border labels still clip for high `thr_y` — clamp pending (see Pending Work).*
  - **Strip text** `strip.text` size 7 → 6.
  - **DUMP `min`/`y_min` (Issue 4):** root cause = template used `x_min`/`y_min` but the plugin params are `min_x`/`min_y` (or `preset_min_x/y`) → silently swallowed by `...`, so `y_min=5000` did nothing and the y threshold sat at the density valley (~3800, floored only by `preset_min_y=3000`). Fix is template-side (`preset_min_x=4900, preset_min_y=5000`); the arg guard now makes this a loud error instead of a silent wrong gate.
  - **CD127+ mindensity parked at ~min (Issue 2): not a bug** — `gate_mindensity(min=2000)` only searches for a valley above 2000; that control pop is unimodal there, so openCyto returns the lower bound. *Lesson: parity/QC checks must mirror the engine's own parsing (`.fasta_parse_pop_list`) and logic (`NOT`) or they generate pure false positives.*

- **2026-06-13 — `clusterFrom` rejected on `gate_logic`/`copy_gate`.** Added guards at the top of both `apply_gating2` branches: a `clusterFrom`/`clusterFrom_x`/`clusterFrom_y` on a `gate_logic` or `copy_gate` row now `stop()`s with a clear message (neither does a data fit, so `clusterFrom` is meaningless). Kept for all other 1D/2D methods, where `clusterFrom` is driver-owned (resolved to `src_pop`, [[gt-gating-removal-apply-gating2]]).

- **2026-06-13 — FRISDA report renders end-to-end (LaTeX caption escaping fix).** First successful PDF render of `260608_Report_FRISDA_0.73.Rmd` (via `rmarkdown::render(ReportMarkdownFile)`, FRISDA line ~164), so the whole LYL273 pipeline — gating, plotting, `panel_QC`, `run_gating_qc`, summary plots — runs through to a PDF. **Bug fixed:** xelatex died with `! Missing $ inserted` at the Panel QC table because the caption `paste("Panel QC:", PanelIteration)` passed `260527_LYL273_[JS]` with **raw underscores** — `knitr::kable(escape=TRUE)` escapes *cell contents* but **NOT the `caption` argument**, so LaTeX read `_` as a math subscript. Fixed by wrapping it: `caption = paste("Panel QC:", knitr:::escape_latex(PanelIteration))` (matches the `knitr:::escape_latex(lbl)` already used in the Layouts `\subsection`). Diagnosis path: error → `tinytex` `latexmk` fail → read `…_FRISDA_0.73.log` for the `! …` + `l.<n>` line. **Env notes (Positron):** the interactive session finds pandoc 3.9.0.2 (`/opt/homebrew/bin`), xelatex (`/usr/local/bin`), `RSTUDIO_PANDOC` points at Positron's bundled quarto tools — render works fine from the console; this was a real LaTeX content bug, not an environment gap. *Lesson: `kable` `caption=` is NOT escaped — escape any dynamic caption with LaTeX specials (`_ # $ % & { }`).*

- **2026-06-13 — layout_gates+gate_logic overlay CRASHED Positron (reverted).** Attempted to make a `layout_gates` panel listing `gate_logic` gates (the `CAR-` panel: `layout_gates="CAR_all, CAR-"`, both gate_logic) draw each gate's component-derived crosshair + label instead of skipping the all-infinite dummy. The new overlay-loop path **hard-crashed the Positron R session during `apply_plotting`** (native segfault — no catchable R error; first execution of the path). Isolated by reverting just that path (the `.gl_xh` hoist was kept; it's inert). Suspected trigger: the per-gate **`gh_pop_get_data()` materialization on the index-based `gate_logic` dummy nodes** (`CAR_all`/`CAR-`; `CAR_all` is also the parent of the `CD4_CD8_CAR_all` copy_gate) that the removed `next` used to skip — i.e. materializing data on a gate_logic node in the overlay loop, not the geoms. **Pending reimplementation must avoid per-gate `gh_pop_get_data` on gate_logic nodes** — derive label counts from `gh_pop_get_count`/stored indices and label positions from the crosshair thresholds, not `exprs()` centroids. The `gate_logic` panel currently shows only its post-loop centroid label (no listed-gate overlays). *Lesson: data materialization on index-based gate_logic nodes inside the plot loop can segfault — read counts/gates, not data.*

- **2026-06-13 — `0_72` → `ExtFunctions_0_73.R`; `parse_args` hardened against malformed `gating_args`.** The LYL273 `apply_gating2` run died with the cryptic `argument is of length zero` (deep in the external_control `pick_in`/`pd_col` path). Root cause was a **template typo** in the `DUMP` row's `gating_args` — `y_min, 5000` instead of `y_min=5000` — which made `eval(parse(...))` throw, so `parse_args()` **silently returned `list()`**; the raw string still matched `grepl("external_control", …)`, so the row took the EC path with a `NULL` `control_col`. Fixed in the workbook (cell F6). Engine hardened (same class as the earlier `.gate_static`/`apply_transformation` guards): **`parse_args` now `stop()`s loudly** with the alias + offending string instead of swallowing the error, and **`pick_in` guards** a `NULL`/empty control column with a clear message. *Lesson (again): a silent `eval`-to-`list()` turns a one-char typo into a crash three calls deep.*

- **2026-06-13 — back in FASTA Enhanced; engine merged from IVANKA.** This CLAUDE.md round-tripped through the IVANKA project to capture all `ExtFunctions*.R` engine changes (AutoSpill, `pop="+-"` split gate, `layout_dims`/`layout_gates`, FCS export args, 2D-rect overlay fix). Returned to FASTA Enhanced carrying **`ExtFunctions_0_72.R`** (IVANKA's last engine state); re-pointed CLAUDE.md at the FASTA orchestration (`FASTA_enhanced.R`, LYL273 Deep-IPT, `20260606_instructions_draft_LYL273.xlsx`, `Older versions/` backups). The IVANKA/E44-validated engine features apply unchanged; re-validation on LYL273 data + the `0_70`→`0_73` re-source line (see Pending Work) are the open items.

- **2026-06-13 — 1D split gate (`pop="+-"`/`"-+"`).** A single-channel rectangle gate row can now emit BOTH polarity populations from one shared per-sample threshold: `<alias>+` = `[thr,Inf]`, `<alias>-` = `[-Inf,thr]`. `apply_gating2` adds both (errors if applied to a non-1D method); `apply_plotting` auto-overlays + labels both (reuses the `.lg_override` overlay path). Validated on E44 (`gate_mindensity` Viability split → `LD+ 85.8%` / `LD- 14.2%`, both labeled with the threshold line; `-+` equivalent; 2D method → clear error).

- **2026-06-13 — per-gate `layout_labels` for every `layout_gates` entry.** The overlay-loop stats label was guarded by `pt_row$gating_method != "gate_logic"`, so a gate_logic row (e.g. `notBeads`) labeled none of its listed gates (only its own centroid label, drawn post-loop). Relaxed the guard to `.lg_override || …` so each listed gate gets its `layout_labels` (gate_name/percent/count). Verified on E44: `notBeads` panel now shows both `Beads` and `notBeads` labels.
- **2026-06-13 — `layout_dims` / `layout_gates` wired into `apply_plotting`.** Previously scaffolded template columns that the engine never read. Now: `layout_dims` overrides panel axes (markernames→channels, forces `computePlotParams` recompute), `layout_gates` replaces the overlay set (resolve_pop'd, drawn only where dims match the axes, suppresses the default primary annotation), both hard-error on unresolvable names. **Debug note:** the override must go *after* the final `overlay_pops <-` computation (~line 2485) — an initial placement after an earlier `overlay_pops` assignment was silently clobbered by a later one; and the alias-keyed `plot_params_cache` had to be bypassed under `layout_dims` or the new axes got the old dims' limits/binwidth. Validated: axis override, gate overlay, dim-mismatch skip, both error paths, and default-unchanged.
- **2026-06-13 — 2D rect overlay fully fixed + visually verified (IVANKA).** The 2D `gate_static` (`notDebris`) overlay went crosshair → `geom_rect` box → still-a-box because **`apply_plotting` draws the primary gate in TWO passes**: path 1 = the `.g_auto` "primary population annotation" block (fixed first), path 2 = the `overlay_pops` loop (line ~2956). Both now draw **finite-edge `geom_segment`s** (open `±Inf` edges omitted). Rendered on E44 `C1 117.fcs` and confirmed: `notDebris` = clean L (`x≥50000` ∧ `y≥21000`, no top/right edges), `Single` = polygon, `Live` = single vline. **Known residual tech-debt:** both passes still draw every gate (harmless overlapping identical draws) — candidate for consolidation later.
- **2026-06-12 — First gating-layout QC + engine hardening (IVANKA).** Two instructions-xlsx errors surfaced in the `notDebris/Single/Live` layout: (a) `notDebris` `gating_args` typo `ymin21000` (missing `=`) → upstream `eval(parse(...))` silently dropped all bounds → fully-open gate → no overlay; (b) `flowjo` placed in `TransformArgs` instead of `TransformType` → `apply_transformation` did nothing → Viability stayed linear → Live population squished (`computePlotParams` was fine). Both fixed in the workbook. Hardening added: **`.gate_static` errors** if no recognised boundary arg survives parsing; **`apply_transformation` emits a QC message** listing marker channels left untransformed (TransformType NA/blank), on every return path. Plotting fix: **2D `rectangleGate` now draws only its finite edges as `geom_segment`s** (open `±Inf` edges omitted) instead of a full-panel vline+hline crosshair; an interim `geom_rect` was rejected because its open edges landed at the `coord` panel boundary. 1D gates still draw a single line — see Plotting Standards.
- **2026-06-12 — FCS export args on `apply_compensation`.** `export_fcs` (→ `Compensated/<basename>_<comptype>.fcs`) and `append_pdata_to_export` (append pData row to keywords via `cytoframe_to_flowFrame`, no mutation of `cs_comp`). See Compensation section.
- **2026-06-12 — AutoSpill compensation implemented + validated (IVANKA/E44).** Replaced the `apply_compensation` `autospill` stub with `build_autospill_controls` + `run_autospill` (see Compensation section). Key decisions through the session: dye↔antigen **vetted from the markers sheet** (not empirical peak — spectral neighbours like `EF780`/`BV786`, `RB705`/`RB780` mislead it); filename auto-matched from metadata `SampleID` with a 2-entry alias (`PB`↔Pacific Blue, `eFluor780`↔Viability); **beads used as-is** (they carry their own pos+neg; AutoSpill has no native separate-negative input — an earlier SC/VC concat design was dropped); **GFP** concatenated only if pos+neg both present; control definition kept as an **in-memory data.frame** → throwaway `tempfile()` (no persistent `fcs_control.csv`, no FCS exports); script-level **`CompType` fallback** to the latest local `*_BackupCompMatrix_IVANKA.csv` if AutoSpill fails. Validated end-to-end on E44 (15×15 matrix, converges, `apply_compensation` applies). Note: `compensate()` does not stamp the matrix into FCS keywords, so `spillover(cs_comp[[i]])` is `NULL` by design — the matrix lives in `compMatrix`.
- **2026-06-11/12 — gt_gating → `apply_gating2` cutover + follow-ups.** Removed `apply_gating` (1481 lines) + `register_plugins` calls + dead `add_boolean`. `gate_logic` index-based (booleanFilter abandoned — see engine note), gained AND/OR/**NOT** + comma-string `gates=`. `gate_mindensity_2d` completed as a 2D quad method. `copy_gate` fully ported (incl. quad/tmix children, under gate_logic). Full-template parity vs old engine confirmed (36 pops) before removal.
- **panel_QC** — root cause of mass mismatches was a **wrong `PanelIteration`** (cPARP panel vs Deep-IPT data), not code. Hardened anyway: `.norm_ag` separator-insensitive antigen match, `-A/-H` suffix strip, **symmetric** viability exemption, broadened invisible-char clean.
- **gate_logic / tmix plot annotation** reworked: no full-panel box, component-derived crosshair, centroid label; tmix quad labels centered.
- **`CAR-` = `NOT(CAR_all)` gating bug** — `recompute()` wiped `CAR_all` to 100% before the `NOT` read it → `CAR-`=0. Fixed by restoring prior gate_logic indices before reading components. **Was masked by the parity test** (the removed engine had the identical bug → both 0). *Lesson: parity proves engine agreement, not correctness.*
- **2026-06-07** (historical, in the now-removed `apply_gating`): per-sample quadGate collapse (→ why we decompose); `pop="-"` finite-boundary vline; layout duplicate-position QC.

---

## Pending Work

**FASTA Enhanced (active):**
1. **Re-validate the `0_73` engine on LYL273 Deep-IPT data** (via `FRISDA_0_73.R`) — ✅ *largely done 2026-06-13:* the full pipeline (gating `apply_gating2`, plotting `apply_plotting`, `run_gating_qc`, `panel_QC` @ `PanelIteration "260527_LYL273_[JS]"`, summary plots) runs end-to-end and renders the FRISDA PDF. Remaining: spot-check the gating/QC outputs for scientific correctness.
2. **`260608_Report_FRISDA_0.73.Rmd` stale naming** — setup chunk + header still grep `^FASCA_.*\.R$` / print "FASCA"/"FASTA"; `scriptName`/`scriptVer` extraction silently fails (no `FASCA_*.R`). Re-point at `^FRISDA.*\.R$`. (LaTeX caption-underscore crash already fixed — see History.)
3. **Re-do the `layout_gates`+`gate_logic` overlay (crash-safe).** The `CAR-` panel (`layout_gates="CAR_all, CAR-"`) should show the CAR-threshold crosshair + a label per listed gate, but the first attempt segfaulted Positron (see History) and was reverted. Reimplement WITHOUT per-gate `gh_pop_get_data()` on gate_logic nodes: derive the crosshair from `.gl_xh` (gate reads only), counts from `gh_pop_get_count`/indices, and label positions from threshold-based quadrant centers. Test incrementally.
4. ✅ **DONE 2026-06-14 — quad y-border value labels fixed.** Root cause was a lazy-`aes()` loop-variable bug (position computed inside `aes()` referencing the `.plot_ylim` loop var, resolved at print time → stale); fixed by baking position into the data.frame. See History. (tmix value labels remain absent by design — polygon quads have no single threshold.)
5. **`collapseDataForGating`** (implemented 2026-06-13) — ✅ *validated 2026-06-14:* exercised on the `CAR-_CD127+` row (`gate_mindensity`, `collapseDataForGating=TRUE`, `groupBy=FDP_lot`, `clusterFrom="CD3+CD19-"`) and confirmed working. **Key clarification:** collapse pools **within `groupBy` groups**, NOT across all samples. LYL273 `FDP_lot` = **5 lots × 7 samples** (C027.004/005/006/012/PBMC), so collapse yields **5 distinct gates** (one per lot), each replicated to its 7 samples. An initial "looks like per-sample mindensity" report was a **cross-lot comparison** (C027.012 vs C027.004 → correctly different gates); the real test is **within-lot identity** (two C027.012 samples must share the gate). To pool ALL samples into one gate, leave `groupBy` blank. `FRISDA_0_73.R` line ~118 calls `apply_gating2(..., collapse_max_events=5e5, collapse_seed=1234)`. Remaining: user's final within-lot identity spot-check.
6. *(legacy, not active)* `FASTA_enhanced.R` still hardcodes `source("ExtFunctions_0_70.R")` (`0_70` is gone). Superseded by `FRISDA_0_73.R`, which uses the glob bootstrap — only matters if anyone re-runs the old script.

**Engine-level (not active):**
7. **`gate_flowmeans` (1D)** — only gate type never exercised in parity. Validate if ever used.
8. **`gate_hinge_quad`** — FlowJo-style tilted quadrant gate, not implemented.
9. **`autospectral` compensation** — still an unimplemented `apply_compensation` stub. (`autospill` is implemented — see Compensation section.)

**IVANKA (separate project — not tracked here):** gating-template↔`markernames(cs)` alignment (`GCC_CAR`/`CD19_CAR` absent from the markers sheet), transformation check, cPARP-leftover repointing in stats/summary, and the `"InVivo CAR Panel_iteration1]"` PanelIteration fix moved with the IVANKA orchestration script.
