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
> The engine reference below is **descriptive** (how the backend behaves); it is **not** an invitation to modify the engine. **The engine in this folder is now `ExtFunctions_0_77.R`** (synced up from `0_73`; see the engine-sync note in History) — the reference below tracks that version.

---

## Project

**FASTA Enhanced** is the shared, automated flow-cytometry analysis engine at **Lyell Immunopharma** — standardized gating, statistics, and layout plots from FCS files using R/Bioconductor (flowWorkspace, openCyto, ggcyto). One shared codebase (`ExtFunctions*.R`) is driven by a thin **orchestration script** — here the **GUI driver `FASTA_wGUI_0_77.R`** (ClaudeR/httpgd MCP-enabled session for interactive GUI development). Downstream studies reuse the same engine via their own orchestration scripts — e.g. **FRISDA** (LYL273 Deep-IPT) and **IVANKA/IVNCA**, the in-vivo CAR PK adaptation, each a **separate project/folder**.

**Owner:** Jatin Sharma (jsharma@lyell.com), Research & Science  
**Study context:** GUI build over the shared engine; backend originally exercised on **LYL273** FDP Deep Phenotyping (Deep-IPT) and **E44** in-vivo CAR PK. CAR-T IND/BLA preparation.  
**Branch:** `main`.

> **Engine sync (2026-06-22): `ExtFunctions_0_73.R` → `ExtFunctions_0_77.R`.** The shared engine was pulled forward through the `0_74`/`0_75`/`0_76`/`0_77` work prototyped in the IVANKA/IVNCA and URSULA/LYL273 forks. New since `0_73` and now available to the GUI backend: **template pre-validation** (`prevalidate=TRUE`) in both `apply_gating2` and `apply_plotting`; **`gate_mixmethod_2d`** (per-axis fusion quad), **`gate_mindensityshoulder`** (1D) + **`gate_mindensityshoulder_2d`**, **`gate_flowclust_2d`** (ellipse), **`gate_quantile`** (1D interval) + **`gate_quantile_2d`** (2D box); `gate_mindensity` native `adjust`/`num_peaks`/`peaks` + `valley="simple"`; expanded `gate_singlet` args; `gate_static` made **data-free** (gates on 0-event parents); `pop="+-"` split `name_pos`/`name_neg`; `copy_gate` `name_suffix` + copy-of-copy hard-stop; `.fasta_check_numeric_args` quoted-numeric guard; **`collapseDataForGating`/`groupBy`** pooled gating; `control_val` `"x>y>z"` preference chains + ambiguity guard; **case-insensitive** channel/marker resolution; duplicate-child / forbidden-char / ambiguous-pop hard-stops; `init_error_log`; plotting: **`DotPlot`/`Pseudocolor` PlotTypes** + `DotPlot_color`, **visible-window binning** + `min_bins`/`max_bins`/`max_bins_per_axis`, threshold labels under `.lg_override`, empty-panel labels, full-path/full-gs scaling fixes; `panel_QC` combined-name fallback + non-fluor exemption. See History for the per-change detail.

> **Codebase round-trip (2026-06-13):** This CLAUDE.md was carried over to the **IVANKA** project so that every engine change made there to `ExtFunctions*.R` was captured in one place. The engine/compensation/plotting sections below describe that shared codebase and are valid here; **script-level wiring details that name IVANKA** (e.g. the `*_BackupCompMatrix_IVANKA.csv` fallback, the E44 QC plot) reflect IVANKA's *orchestration* and are not necessarily wired into the GUI driver. IVANKA-only adaptation tasks moved with the IVANKA project.

> **Working convention:** After successfully implementing or validating *any* change, **prompt the user to update this CLAUDE.md** (and the History / Pending Work sections) before moving on — don't silently skip it. CLAUDE.md is the single source of truth for this project.

> **Engine changes go in `ExtFunctions*.R` (it is forked into other projects).** ALL analysis/engine logic — gating, gate plugins, `run_gating_qc` checks, `apply_plotting` behavior, compensation, `panel_QC`, stats — belongs in `ExtFunctions*.R`, never duplicated into the orchestration/GUI script or a report. **For this GUI project specifically, the engine is a read-only swappable dependency (see GUI ground rules) — engine fixes belong in the engine project, not here.** Only genuinely GUI-/presentation-specific wiring stays in the GUI scripts: data paths, `instructions_filename`, `PanelIteration`, the locate-and-source bootstrap, and the GUI/MCP plumbing.

---

## Key Files

| File | Purpose |
|------|---------|
| `ExtFunctions_0_77.R` | **The entire backend engine** (read-only black-box dependency for the GUI). Gate plugins, `apply_compensation`, **`build_autospill_controls`/`run_autospill`**, `apply_transformation`, **`apply_gating2`** (gating engine), `apply_plotting`, `panel_QC`, `run_gating_qc`, helpers |
| `FASTA_wGUI_0_77.R` | **Active GUI driver** — installs/loads packages (incl. `ClaudeR`, `httpgd` via GitHub fallbacks), bootstraps + sources `ExtFunctions*.R` via the `grepv("ExtFunctions", …)` glob (line ~65), reads `instructions*.xlsx`, runs the pipeline. Orchestration/GUI only. `PanelIteration <- "Placeholder"` |
| `20260623_instructions.xlsx` | **Live** instructions workbook (matched by `list.files(pattern="instructions+\\.xlsx")`): `metadata`, `markers`, `gating_template` sheets → `finalInstructions` |
| `260608_Report_0.73.Rmd` / `.pdf` | Carried-over PDF report scaffold (stale `0.73` naming) |
| `GatingTemplate.csv` | Legacy minimal draft template — vestigial (the live template is the xlsx `gating_template` sheet) |
| `Versions/0_73/` | Prior-version backups (`ExtFunctions_0_73.R`, `FASTA_wGUI_0_73.R`). Safe because the ExtFunctions glob is **non-recursive** (top level only) |
| `TestSamples/` | Test FCS data for GUI development |

**Filename is versioned** (`ExtFunctions_0_77.R`). The script locates it via `grepv("ExtFunctions", list.files(pattern=".*\\.R$"))` — a **non-recursive** glob, so keep **exactly one** `ExtFunctions*.R` at the **top level** (backups in `Versions/` are not matched and so are not sourced — fine; but two at top level would both source and silently double-define everything).

**Re-load gotchas (bit us repeatedly):**
- After swapping `ExtFunctions_0_77.R`, you **must re-`source()`** it — re-running a pipeline call uses the stale in-memory function otherwise.
- After editing the xlsx, you **must re-read** `finalInstructions <- setNames(map(listOfSheets, ~read_excel(...)), listOfSheets)` — and re-execute downstream assignments (e.g. `PanelIteration`).
- A **gating** fix needs `gs <- apply_gating2(...)` re-run (re-sourcing alone doesn't recompute an existing `gs`); a **plotting-only** fix needs only `apply_plotting(gs, ...)`.

---

## Gating engine: `apply_gating2` (gt_gating-free)

**As of 2026-06-11, `apply_gating2` replaced the old `apply_gating`.** The old function and ALL its openCyto `gt_gating` machinery were **removed**. `apply_gating` no longer exists.

`apply_gating2(gs, template, mc.cores = 1, verbose = TRUE, prevalidate = TRUE, collapse_max_events = 5e5, collapse_seed = 1234)` — sequential driver:
- Walks template rows **top-to-bottom**. Parents always exist before children → **no DAG resolution, no precompute cache, no deferrals**.
- Per row: resolve parent + channels → pick source data → compute a per-sample gate via `.fasta_compute_gate` → add via `gs_pop_add`.
- Clears `.fasta_gate_cache` at start so plugins compute fresh.

**Pre-validation (`prevalidate = TRUE`):** a static precheck runs over ALL active template rows **before any gating** and aborts with a **table of every issue** (`row, alias, issue, detail`) — so a malformed template surfaces all problems at once. Checks: `gating_args` parses; method known; `pop` polarity; **forbidden chars in names** (`.FASTA_FORBIDDEN_NAME_CHARS = c("/", ",")` in the alias and quad `name_x`/`name_y`); dims count + channel resolution (`resolve_ch`); arg names (vs `formals`, `gate_static` exempt) + types (`.fasta_check_numeric_args`); per-method required/forbidden args (`gate_logic`→`gates`; `copy_gate`→`copy_from`; `gate_tmix_2d`→`K`; `gate_mixmethod_2d` per-axis exclusivity); `gate_singlet` `sidescatter` resolves; collapse×external_control and collapse×per-axis-mixmethod clashes; `control_col[_x/_y]` resolvable; **duplicate child population under a shared parent** (enumerates each row's produced child names — quad `name_x±name_y±`, 1D `alias`/`alias±`, gate_logic `alias`, copy_gate→source's children, all suffix-aware — and flags any `(parent-string, child)` made by >1 row; `gs_pop_add` can't create duplicate sibling nodes, so the 2nd silently fails → blank panels). Population + cache clearing happen **after** validation, so an abort leaves `gs` untouched. Escape hatch `prevalidate = FALSE`. **A prevalidate abort embeds a single-line issue summary in the `stop()` message** (`[row R 'alias' issue: detail]; …`) so the error log is self-contained, not just "see table above".

**`.fasta_compute_gate(fr, method, channels, filterId, args)`** — the one dispatcher: calls `.gate_<method>` (our plugins) or a builtin (`singletGate`) **in `per_sample` mode** on whatever flowFrame the driver feeds it. (So the driver, not the plugin, decides which data/control to use.)

**Per gate type, how the driver adds it:**
- **1D rect** (mindensity / FMO / flowmeans / static / shoulder / quantile): named per-sample list → `gs_pop_add`. `pop="-"` → flip gate to `[-Inf, thr]`. **`pop="+-"` or `"-+"` → 1D split:** from the one shared per-sample threshold, add BOTH `[thr, Inf]` as `<alias>+` and `[-Inf, thr]` as `<alias>-`. **`name_pos` / `name_neg`** (driver args) override the `+`/`-` child names **independently** (default `<alias>+`/`<alias>-`; set just one and the other keeps the default — e.g. `name_pos='Live', name_neg='Dead'`). Strictly 1D — a non-1D method with `pop="+-"` hard-errors. `name_pos`/`name_neg` are in `.FASTA_DRIVER_ARGS` (dropped before the 1D plugin), validated against `.FASTA_FORBIDDEN_NAME_CHARS`, counted by the duplicate-child check.
- **2D quad** (`gate_flowmeans_2d` / `gate_FMO_2d` / `gate_mindensity_2d` / `gate_mindensityshoulder_2d` / `gate_mixmethod_2d`): plugin returns a `quadGate` → **decomposed into 4 `rectangleGate`s** (sidesteps the `gs_pop_add` quadGate per-sample collapse).
- **`gate_mixmethod_2d`** (per-axis fusion): each axis is **independently** EITHER a fresh 1D method (`method_x`/`method_y` + `args_x`/`args_y`, per-axis sourcing) OR a copied boundary (`copy_x`/`copy_y`). No per-axis arg → driver runs the plugin through the normal per_sample/collapse path; per-axis sourcing → dedicated **`mix_per_axis`** branch assembles the per-sample `quadGate`.
- **`gate_tmix_2d`**: returns a `filters` of 4 `polygonGate`s → added per sample as-is (curved contours preserved).
- **`gate_flowclust_2d`**: returns ONE `fcEllipsoidGate` (ellipse) — kept as-is (NOT decomposed; not a quad).
- **`gate_logic`**: **index-based** (NOT booleanFilter — see below). Dummy pass-through gate + per-sample AND/OR/NOT via `gh_pop_set_indices` (`NOT = parent_idx & !inner`). A gate_logic referencing **another** gate_logic restores the prior nodes' indices right after its own `recompute()` (else the dependency reads the wiped 100% dummy — the `CAR-` bug). Indices re-applied after the loop; gate_logic **descendants are re-gated** against the restored parent.
- **`copy_gate`**: copies gate(s) to a new parent — single-node sources, 2D-quad-children (tmix polygons or rect quads, matched via `name_x/name_y`), quadGate-node; **rejects** `gate_logic` sources and **copy-of-a-copy** (hard-stop, names the traced original via `.copy_root`); skips inactive sources. `name_suffix` appends to every copied node name (the escape hatch for a deliberate same-parent copy).

**Control modes:** `per_sample` (own data); `external_control` (driver feeds the **control sample's** data to the plugin — single `control_col/val`, per-axis `control_col_x/y`, and/or `groupBy` per-group control); `clusterFrom` (compute from the named population's data).

**`control_val` preference chain `"x>y>z"`:** any `control_val` / `control_val_x` / `control_val_y` may encode a **`>`-separated ordered preference list** — the driver tries each token in order and uses the **first** with a matching sample in the pool (e.g. `control_val="FMO>FM2"` = "FMO, else FM2, …"). **Ambiguity guard:** if the chosen preference matches **>1** sample, the driver **stops and lists the candidates**, requiring an explicit `control_idx`/`control_idx_x`/`control_idx_y` (1-based, range-checked) to disambiguate.

`.FASTA_QUAD_2D_METHODS = c("gate_flowmeans_2d","gate_FMO_2d","gate_tmix_2d","gate_mindensity_2d","gate_mindensityshoulder_2d","gate_mixmethod_2d")` — used by `apply_plotting`/`run_gating_qc` for quad detection. (`apply_gating2` dispatches by the *returned gate type*, not this list.) `.FASTA_MIX_1D_METHODS = c("gate_mindensity","gate_flowmeans","gate_FMO","gate_static")` — the 1D methods `gate_mixmethod_2d` can fuse per axis. `.FASTA_DATA_FREE_METHODS = c("gate_static")` — exempt from the driver's 0-event-parent guards (a fixed gate never reads `exprs(fr)`).

### Why `gate_logic` is index-based, not `booleanFilter`
`booleanFilter` is the "native" gt_gating-free path, but in this flowWorkspace build **its references never resolve at `recompute()`**, **and every failed boolean recompute corrupts the shared cytoframe backing** (propagates back to `cs_comp`, needs an FCS reload). The index approach (`gh_pop_set_indices`) is the proven mechanism and never corrupts on failure. Don't reach for `booleanFilter` here.

---

## Custom Gating Methods (plugins in `ExtFunctions_0_77.R`)

| Method | Description |
|--------|-------------|
| `gate_FMO` | 1D FMO rectangleGate |
| `gate_FMO_2d` | 2D FMO quadGate. Per-axis controls via `control_col_x/y`, `control_val_x/y` |
| `gate_flowmeans` | 1D flowMeans gate |
| `gate_flowmeans_2d` | 2D flowMeans quadGate, per-axis control |
| `gate_mindensity` | 1D min-density (wraps openCyto). `upper_mult`/preset args honored. **Native openCyto knobs forwarded**: `adjust` (KDE bandwidth mult; default 2, SMALLER=0.5–1 resolves shallow valleys but non-monotonic — tune+inspect), `num_peaks` (force N-peak detection), `peaks` (explicit peak LOCATIONS). `min`/`max` are **DATA FILTERS** (drop events before density), NOT a search window — use `gate_range=c(lo,hi)` to constrain the cutpoint. **`valley=`**: `"KDE"` (default, openCyto valley-between-peaks) or **`"simple"`** = non-parametric histogram MIN-count within `gate_range` (no bimodality needed; tie broken by widest run). `"simple"` **requires `gate_range`** + uses **`window`** (bin width, default `diff(gate_range)/100`); `adjust`/`num_peaks`/`peaks` KDE-only, `window` simple-only |
| `gate_mindensityshoulder` | **1D "valley-else-shoulder"** min-density. Deterministic, no distribution modeling. Builds one KDE; **bimodal → delegates valley to `gate_mindensity`**; **else → SHOULDER** via `shoulder_method`: **`"curvature"`** (default, max-curvature elbow = foot of dominant peak) or **`"edge"`** (`peak_x + edge_k×right-HWHM` = negative's upper edge; for 2D-prominent / 1D-diffuse pops where curvature anchors too far left). Reuses all `gate_mindensity` args + `grid_n`/`bw`/`shoulder_method`/`edge_k`/`shoulder_floor`/`mode_min_height`/`mode_min_dip`. Returns `rectangleGate(c(thr,Inf))` — works with `pop="+-"` |
| `gate_mindensity_2d` | 2D min-density quadGate (in `.FASTA_QUAD_2D_METHODS`; `name_x/name_y`, combined-key precompute caching) |
| `gate_mindensityshoulder_2d` | **2D analogue** of `gate_mindensityshoulder` — independent valley-else-shoulder per axis → `quadGate`. Each axis computed by the **shared `.fasta_shoulder_threshold` helper** (same code path as 1D, can't drift). Args mirror `gate_mindensity_2d`'s `_x`/`_y` pattern + per-axis shoulder controls. In `.FASTA_QUAD_2D_METHODS` |
| `gate_tmix_2d` | 2D quad wrapping openCyto `gate_quad_tmix` (joint t-mixture). Returns 4 `polygonGate`s; rectangular crosshair when well-resolved is expected. Args: `K`, `quantile1/3` (def 0.8), `usePrior`, `prior`, `trans`, `name_x/y`, `clusterFrom`, `groupBy` |
| `gate_flowclust_2d` | **2D flowClust ellipse** — fits K-component bivariate t-mixture, returns ONE `fcEllipsoidGate` around the selected cluster (kept as ellipsoid, coerced to polygon at render only). Args: `K` (2), `target=c(x,y)` (**always set for reproducibility**), `quantile` (0.9, size knob), `transitional`/`translation`, `min`/`max` (**per-channel `c(x,y)`, length-2 — a SCALAR hard-errors; use `0` not `NA` for an open axis**), `usePrior`/`prior`/`trans`/`min.count`/etc. Driver owns sourcing. NOT in `.FASTA_QUAD_2D_METHODS` (single pop) |
| `gate_singlet` | Singlet gate (wraps `flowStats::singletGate`); needs 2 dims (area, height). Args (all `singletGate` pass-throughs): **`prediction_level`** (0.99, tightness), **`wider_gate`**, **`sidescatter='<chan\|marker>'`** (3-var fit; resolved against the frame, factor-coerced), **`subsample_pct`**, **`maxit`**. The inert `upper_mult`/`preset_min`/`preset_max` were **removed from formals** so the arg guard now **rejects** them here |
| `gate_static` | Fixed rectangleGate from `xmin/xmax/ymin/ymax` (or `min/max` 1D). **Data-free** (in `.FASTA_DATA_FREE_METHODS`): ignores parent events, placed even on a **0-event** parent. `pop="+-"` child names from `name_pos`/`name_neg` |
| `gate_quantile` | **1D quantile INTERVAL** — gate = `[Q(quantile_min), Q(quantile_max)]`. Defaults `quantile_min=0`/`quantile_max=1` (whole range; **0→-Inf, 1→+Inf** open ends). Narrow to trim: `0.99,1` = top 1%; `0.01,0.99` = middle 98%. 1D `rectangleGate`; use `pop="+"`. NOT in `.FASTA_MIX_1D_METHODS`. Validates `0<=min<max<=1` |
| `gate_quantile_2d` | **2D quantile BOX** `[Qx_min,Qx_max]×[Qy_min,Qy_max]` from per-axis quantiles (`quantile_min_x/max_x`, `quantile_min_y/max_y`; same 0→-Inf/1→+Inf rule). Returns a SINGLE box, **not** 4 quadrants — **not** in `.FASTA_QUAD_2D_METHODS` (drawn like a 2D `gate_static`) |
| `gate_logic` | Boolean **AND / OR / NOT** of existing populations. `gates=` accepts comma-string `'A,B'` or `c('A','B')`. Index-based (see engine note) |
| `copy_gate` | Copies gate(s) from another population (incl. tmix/quad children, under a gate_logic parent). **`copy_from` must point at the ORIGINAL gate, NOT another `copy_gate`** (copy-of-a-copy is rejected). **`name_suffix`** (default `""`) appended to every copied node name — the escape hatch for a deliberate same-parent copy. Validated against `.FASTA_FORBIDDEN_NAME_CHARS` |
| `gate_mixmethod_2d` | **2D per-axis fusion quad** — each axis independently = a fresh 1D method (`method_x`/`method_y` ∈ `.FASTA_MIX_1D_METHODS`, with `args_x`/`args_y` = nested `list()`) **or** a copied boundary (`copy_x`/`copy_y`). Per-axis sourcing: `control_col_x/y`+`control_val_x/y`, `clusterFrom_x/y`, or per_sample. `name_x/y` name the 4 children. Per axis, `method_`/`copy_` are mutually exclusive. `collapseDataForGating` is **blocked** with per-axis sourcing. Returns a `quadGate` → decomposed to 4 rects |

**Numeric-arg type guard (`.fasta_check_numeric_args`):** the unrecognized-arg guard catches typo'd arg *names* but not wrong *types* — a numeric arg passed as a **quoted string** (e.g. `min_x="4000"`) sails through into the gating math where R does **lexicographic** comparison (`10000 > "4000"` is `FALSE`) → silently wrong gate. `.fasta_check_numeric_args(args, ctx, allow_string=…)` now **stops loudly** naming the arg + the exact fix. Wired into `.fasta_compute_gate` (all fn-based methods incl. `gate_static`) and `.mix_axis_threshold`; `allow_string` exempts genuinely-character args.

**gt_gating-era leftovers:** the `register_plugins()` calls and the `add_boolean` helper were **removed (2026-06-11)**. The plugins' `method="external_control"` branches were **kept by decision** — inert for `apply_gating2` (always `per_sample`) but keep the `.gate_*` functions reusable as standalone dual-mode tools. `.fasta_gate_cache` is kept (plugins still write `quad_names` to it; `apply_gating2` clears it at start).

---

## Compensation (`apply_compensation`)

`apply_compensation(cs, metadata, autospill_matrix = NULL, export_fcs = FALSE, append_pdata_to_export = FALSE)` dispatches **per sample** on the metadata `Compensation_type` column:
- **`self`** — apply the acquisition-embedded `$SPILLOVER`/`$SPILL` matrix.
- **`unmixed` / NA / blank** — no compensation (data already unmixed/compensated).
- **`<file>.csv` / `<file>.xlsx`** — a square spillover matrix file in the working dir, applied via `compensate()`; files are loaded once and cached if reused.
- **`autospill`** — applies the **experiment-level** AutoSpill matrix passed via `autospill_matrix` (one matrix shared by all `autospill` samples; errors if `NULL`). All-`autospill` is fast-pathed via `compensate(cs, autospill_matrix)`. `apply_compensation` only **applies** the matrix — it never runs AutoSpill or touches files.
- **`autospectral`** — still a `stop()` **stub**.

**FCS export args** (one internal `do_export()` wraps every return path):
- **`export_fcs`** (default FALSE) — when TRUE, writes each compensated frame to a **`Compensated/`** folder in `getwd()`, named `<original_basename>_<comptype>.fcs` (`<comptype>` = the sample's resolved `Compensation_type`, sanitized). Exports whatever sample set was passed in.
- **`append_pdata_to_export`** (default FALSE) — when TRUE **and** `export_fcs` is TRUE, that sample's `pData` row is appended to the exported frame's FCS keywords. The frame is detached via `cytoframe_to_flowFrame` first, so `cs_comp` itself is **not** mutated.
- **Note:** set-level `pData(cs)` is *not* otherwise embedded in frames — plain `export_fcs` carries only the original acquisition keywords; use `append_pdata_to_export` to bake the layout metadata in.

### AutoSpill (`build_autospill_controls` + `run_autospill`, implemented 2026-06-12)
AutoSpill (`carlosproca/autospill` v0.2.0, GitHub-only) computes **one experiment-level** spillover matrix from single-color controls. Two engine helpers stage and compute it; the **orchestration script** wires them and writes the backup CSV (kept out of `apply_compensation`):
- **`build_autospill_controls(markers, controls, control_dir, ...)`** — the dye↔antigen pairing is **VETTED from the markers sheet** (`dye=dims`, `antigen=markers`); the method does **not** guess channels (an empirical max-signal peak is unreliable for spectral neighbours, e.g. `EF780`/`BV786`). Each control's **filename is auto-matched** from `controls` by normalised fluorophore, with `alias = c(PB→"Pacific Blue", eFluor780→"Viability")`. **Hard-errors** on any unmatched vetted row.
  - **Beads used as-is:** SpectraComp/ViaComp beads carry their own positive+negative populations in one file — **symlinked unchanged** (AutoSpill gates pos/neg within each). The `US SC`/`US VC` unstained files are therefore **not used**.
  - **GFP exception:** two files (`GFP+`/`GFP-`) → concatenated to `"GFP single stain pos_neg.fcs"`; a **single** GFP file → assumed already mixed.
  - Symlinks + (optional) GFP concat go to a **temp dir**. Returns `list(control_dir, control_def, mapping)` where **`control_def` is an in-memory data.frame**.
- **`run_autospill(control_dir, control_def, scatter_param = c("FS00-A","SS02-A"), param_set = "minimal")`** — wraps the AutoSpill `minimal` workflow. **`minimal` writes nothing to disk.** `control_def` may be the **data.frame** (preferred) or a CSV path; AutoSpill's readers hardcode `read.csv()` and read it twice, so `run_autospill` serialises the dataframe to a **throwaway `tempfile()`** for the call only. Errors propagate so the script can fall back. Returns `refined$spillover` (square, diagonal = 1, named by detector channel).
- **Scatter gotcha:** AutoSpill reads the **raw** FCS (pre-rename) and defaults scatter to `c("FSC-A","SSC-A")`; this instrument names scatter `FS00-A`/`SS02-A` — so `run_autospill` overrides `asp$default.scatter.parameter` to the **raw** names.
- **Channel names safe:** `-` is not a forbidden char, so `FL25-A`-style names are not mangled → matrix row/col names match `cs` channels for `compensate()` directly.
- **Script pattern + fallback:** `build_autospill_controls(...)` → `compMatrix <- tryCatch(run_autospill(...), error → read latest local *_BackupCompMatrix_*.csv` + set `CompType <- "Local matrix")` → `apply_compensation(cs, metadata, autospill_matrix = compMatrix)` → **on success only**, write `<YYYYMMDD>_<BenchlingID>_BackupCompMatrix_*.csv`.
- **Panel-resolution note:** AutoSpill can emit off-diagonal spillover **>1** for near-identical dyes — expected behaviour, a panel-design signal, not a bug.

**`autospectral`** is still an unimplemented stub.

---

## Gating Template Columns

Standard openCyto columns plus custom: `collapseDataForGating`, `groupBy`, `active_gate` (TRUE/FALSE — skip if FALSE), `layout_dims`, `layout_gates`, `layout_row`, `layout_col`, `layout_index`, `PlotType` (`Scatter`/`Pseudocolor`/`DotPlot`/`Histogram` — case-insensitive; `Scatter`==`Pseudocolor`, `DotPlot` = same panel via `geom_point`), `label_col`, `layout_labels` (which of gate_name/percent/count to show).

**`collapseDataForGating` / `groupBy`:**
- **`collapseDataForGating`** (TRUE/FALSE/blank) — when TRUE, POOL the parent (or `clusterFrom`) events across the group, compute ONE gate via the row's method, and replicate it to every sample in the group. Reproducibly downsampled (per-sample cap = `collapse_max_events/n`, fixed `collapse_seed`; pools only the gating channels). **Mutually exclusive with `external_control`** (errors). Excluded from the identical-threshold QC checks. Pools **within `groupBy` groups**, not across all samples.
- **`groupBy`** — grouping for `collapseDataForGating` (and `external_control`'s per-group control): a pData column (or `":"/","` combo, e.g. `"Batch:Day"`), or a number `N` to group every N samples; blank = all samples in one group.
- `collapse_max_events` (def `5e5`) and `collapse_seed` (def `1234`) are **`apply_gating2()` arguments, NOT template `gating_args`** — set them in the orchestration call.

**`layout_dims` / `layout_gates`:**
- **`layout_dims`** — comma-separated **markernames** (resolved to channels via `resolve_ch`); when non-blank, overrides the panel axes for that row. Forces a `computePlotParams` recompute **over the FULL `gs`** (all samples), so a `layout_dims` panel and a normal panel on the same parent are scaled/binned identically. Unresolvable marker → hard error.
- **`layout_gates`** — comma-separated **gate/population names** (resolved via `resolve_pop`); when non-blank, **replaces** the overlay set — those gates drawn **only where their dims match the panel axes** (non-matching listed gate silently skipped), and the default primary-gate annotation + gate_logic crosshair suppressed. **Each listed gate also gets its own `layout_labels` stats label** — including on a `gate_logic` row.
- Both blocks sit **after** the final `overlay_pops`/channel computation.

**Layout QC:** `apply_plotting` stops with a descriptive error if two gates share the same `layout_index + layout_row + layout_col`.

---

## panel_QC

`panel_QC(gs, PanelMasterInventoryLocation, PanelIteration)` — checks GatingSet markernames/channels against a master panel inventory **iteration**. Returns a data.frame of mismatches (empty = clear).

- **`PanelIteration` MUST match the data's assay panel.** A wrong iteration flags everything.
- **Antigen matching** is separator-insensitive via `.norm_ag` (strips `_ - . whitespace`): `PD-1`==`PD1`, `HLA-DR`==`HLADR`. `/` and parens are **not** stripped → `TCF7/1`, `CD197(CCR7)` need a `Marker_Alias`.
- **Fluorochrome matching** via `.norm_fluor`: strips `-A`/`-H` channel suffix, spaces/separators, and maps `AF###`↔`Alexa Fluor ###` (perl lookahead `af(?=[0-9])`).
- **Combined antigen+fluorophore fallback (`.norm_all`):** when a GS markername carries BOTH antigen and fluorophore in one string (e.g. desc `"CD3 BV421"`), a looser normalizer clears the row as a LAST RESORT if both the antigen and the fluorophore strings are substrings of the marker's combined (channel + desc) name.
- **Viability dyes are exempt both directions** (panel "Live Dead eFluor 780" ↔ GS marker "Viability"); both GS→panel checks use `.is_viability(desc, channel)` (a generically-named viability channel no longer leaks through).
- **`.is_nonfluor` exemption:** scatter (`FSC/SSC/spd`), time (`^T#`), event-info are not antibody conjugates → exempted from the GS→panel checks. **GFP stays flagged** (user adds it to the master inventory).
- `.clean` strips zero-width/invisible Unicode. `Marker_Alias` column (inventory) is preferred over `Antigen` for matching.
- A **panel marker not run in a sample** correctly flags — true positive, not noise.

---

## Plotting Standards

**Always use `flow_bin2d()` and `computePlotParams()` — never `geom_hex()` or `ggcyto_par_set()`.** Panels are plain `ggplot` on extracted `exprs()` (columns renamed `x_var`/`y_var`) — *not* `ggcyto` objects, so `geom_gate()` does not apply.

```r
params <- computePlotParams(gs, pop, xchannel, ychannel, adaptiveBins(gs, pop))
p <- ggplot(dat, aes(x=x_var, y=y_var)) +
  flow_bin2d(binwidth = params$binwidth) +
  coord_cartesian(xlim = params$xlim, ylim = params$ylim, clip = "on")
```

**PlotType `Scatter` / `Pseudocolor` / `DotPlot` / `Histogram`.** `.FASTA_PLOT_TYPES = c("scatter","pseudocolor","dotplot","histogram")` (case-insensitive). `Scatter` and `Pseudocolor` are **identical** (both `geom_bin2d` + log10 `scale_fill_gradientn` density). `DotPlot` renders the **same panel** (axes, gate overlays, crosshairs, threshold lines, stat labels, limits, faceting — all identical) but swaps the base layer to **`geom_point(size=0.6, alpha=0.5, colour=DotPlot_color)`** — a per-row choice for **low-event-count pops** where binning can't smooth. Point colour set via `apply_plotting(..., DotPlot_color = "black")` (default `"black"`).

**`apply_plotting` bins over the VISIBLE window; `max_bins_per_axis` is the OOM backstop.** `binwidth` is `dense(P5–P95)/bins`. `geom_bin2d` sizes its grid to the data range it is handed, so off-view tail events used to inflate the grid extent → grainy bins (and a degenerate channel could OOM-crash). **Fix (two coupled parts):** (1) `apply_plotting` **filters each panel's data to the plotted window** (`.plot_xlim`/`.plot_ylim` ± 1 binwidth) before `geom_bin2d` — nothing visible changes (those events are clipped anyway; counts/% labels come from `gh_pop` counts, not the binned data), but the grid spans only the visible range; (2) `computePlotParams`'s safety cap is sized to the **visible span** (`diff(xlim)/binwidth ≤ max_bins_per_axis`), so the fine dense binwidth is kept for normal panels. **`apply_plotting(..., max_bins_per_axis = 2000)`** (default 2000) is a pure **backstop** — RAISE for finer bins on a degenerate window (>5000 warns), LOWER for safer. When the cap bites it emits `[computePlotParams] binwidth capped … needs ~N bins/axis`. Separately, **`min_bins`/`max_bins`** (def 50/300) clamp `adaptiveBins`' `round(sqrt(median events))` — lowering `min_bins` coarsens only low-count panels (selective; high-count panels byte-identical).

Auto gate annotation — `apply_plotting` annotates the primary gate in **two passes** (the `.g_auto` primary block *and* the `overlay_pops` loop). The primary block runs only when **`!.lg_override`**; under `.lg_override` (set by `layout_gates`, or a `pop="+-"/"-+"` split) drawing routes through the `overlay_pops` loop. Each routes the `rectangleGate` case by how many plotted axes the gate constrains:
- **1D rect** (constrains one plotted axis — e.g. `Live` = `[-Inf,thr]`) → a single **vline/hline** at the finite boundary.
- **2D rect** (constrains **both** axes — e.g. `notDebris` `[50000,Inf]×[21000,Inf]`) → draws **only the finite edges as `geom_segment`s** (open `±Inf` edges omitted; extended past the view so `coord` clip trims with no cap) — the L-shaped corner.
- **quadGate children** → crosshairs.
- **tmix polygon quads** → polygon outline; labels at **quadrant centers**. **No threshold VALUE label** (a polygon has no single axis-aligned threshold — by design). Same for a `copy_gate` of a tmix.
- **ellipsoidGate** (`gate_flowclust_2d`) → outline coerced to polygon at render only; **no h/vline value** (tmix convention).
- **gate_logic** → **no box**; a **crosshair derived from the component quadrant thresholds** (recursing one level, so `NOT(CAR_all)` finds the CAR thresholds); a **label at the population's event centroid**.
- **quad children matched by EXACT name** (not grepl) — so a shorter sibling name (`CD45RA` vs `CD45RAt`) can't bleed children into the wrong panel.

**Pre-validation (`apply_plotting(..., prevalidate = TRUE)`):** static layout/plotting precheck mirroring `apply_gating2`'s — iterates the to-be-plotted rows, collects ALL issues into one table, aborts **before rendering any panel**. Checks PlotType ∈ `.FASTA_PLOT_TYPES`; numeric `layout_row`/`layout_col`; `gating_args` parse; `layout_dims`/`dims`/`parent`/`layout_gates` resolve; `label_col` ∈ pData; `layout_labels` ∈ `.FASTA_LAYOUT_LABELS`; cross-row refs; and **duplicate `layout_index+row+col`**. Static-only (does NOT predict empty/0-event panels). Escape hatch `prevalidate = FALSE`.

**Threshold VALUE labels (the red numbers): bake position into the data.frame, never compute it inside `aes()`.** `apply_plotting` builds ggplots in a loop printed LATER, so an `aes()` referencing a loop-local (`.plot_xlim`/`.plot_ylim`) resolves at print time against the loop's FINAL value → label offscreen. Compute positions immediately (in the `data.frame(...)` call) and map only baked columns in `aes`. (The 2026-06-14 y-label bug.)

**1D threshold labels under `.lg_override`.** The threshold VALUE label is now re-emitted by the `overlay_pops` loop's 1D vline/hline branches (mirroring the primary block's style), **deduped per panel** by rounded threshold — so 1D gates drawn via the overlay loop (`pop="+-"/"-+"` splits, `layout_gates` 1D gates) show the number, not just the line. 2D boxes in the overlay loop stay unlabeled. **Empty-panel (0-event) 1D gates** also get their threshold value label (positioned from `params$xlim/ylim`, finite even when empty). **Subtitle** shows up to **3 ancestry levels** (`great-grandparent/grandparent/parent | Method:`), built from the resolved `parent_path` (correct for aliases, quad children, explicit paths).

---

## Per-Sample Gate Application Pattern

`gs_pop_add(gs, namedListBySample, parent)` applies per-sample for rectangleGate/polygonGate. **quadGates collapse** (C++ `cpp_addGate` — last sample wins), which is why `apply_gating2` decomposes 2D quads into rectangleGates. For a manual per-sample fix: `for (sn in sampleNames(gs)) gh_pop_set_gate(gs[[sn]], node, per_sn_gate); recompute(gs)`. **`flowFrame(matrix)` is unsafe on 0-row input in this build** (writes `$TOT=0`/`$PnR=-Inf` → `read.FCS` rejects as corrupted) — the collapse `.pool_group` guards this (data-dependent → clean `stop`; data-free → 1-row placeholder).

---

## History / notable fixes

### Engine sync `0_73` → `0_77` (2026-06-22) — backend pulled forward from the IVANKA/IVNCA + URSULA/LYL273 forks

The GUI backend was bumped four versions. Grouped summary of what each version added (full per-change rationale lives in the engine project's CLAUDE.md):

- **`0_74` (engine sync from FASTA Enhanced master):** **`gate_mixmethod_2d`** (per-axis fusion quad + `.mix_axis_threshold` + `.FASTA_MIX_1D_METHODS`); **template pre-validation** `prevalidate=TRUE` in BOTH `apply_gating2` and `apply_plotting` (all issues → one table, abort before mutating `gs`/rendering); **`.fasta_check_numeric_args`** quoted-numeric guard; **`collapseDataForGating`/`groupBy`** pooled-per-group gating; plotting fixes (quad children **exact-name** match, **baked-position** threshold labels). The `-A` markername strip is commented out (markernames keep the `-A` suffix; template `dims` aligned to match).
- **`0_75`:** **`gate_mindensityshoulder`** (1D valley-else-shoulder, `curvature`/`edge`) + **`gate_mindensityshoulder_2d`** (+ shared `.fasta_shoulder_threshold` helper, 1D/2D can't drift); **low-event graceful fallback** (the three 1D threshold methods now `warning()`+conservative-gate instead of `stop()` killing the run); **`gate_static` data-free** (gates on 0-event parents; `.FASTA_DATA_FREE_METHODS`); **copy-of-a-copy hard-stop** (`.copy_root`); **ambiguous-pop hard-stop** (all four resolution paths list every match); **`name_pos`/`name_neg`** split-child names; **duplicate-child prevalidate** + `copy_gate` `name_suffix`; **forbidden-char guard** (`.FASTA_FORBIDDEN_NAME_CHARS = c("/", ",")`); root-safe child lookup (`.fasta_children_of`); subtitle grandparent + full-parent-path fixes; quad NA-limit fix (`na.rm` on xlim/ylim/dense); `.pool_group` empty-group `flowFrame` fix; **`init_error_log`** (global `options(error=)` handler appending to a log file in the dashboard format).
- **`0_76`:** **`gate_quantile`** (1D interval) + **`gate_quantile_2d`** (2D box); **`DotPlot`/`Pseudocolor` PlotTypes** + `apply_plotting(DotPlot_color=)`; **case-insensitive** channel/marker resolution (`resolve_ch` + `.rch`, ambiguity-guarded); `apply_plotting` **`min_bins`/`max_bins`** + **`max_bins_per_axis`** and **visible-window binning** (supersedes the `0_75` full-range OOM floor — fixes grainy high-count panels, no crash); `gate_mindensity` **`valley="simple"`** (non-parametric histogram min-count); prevalidate aborts embed the issue detail in the `stop()` message; `panel_QC` combined antigen+fluorophore fallback + both-direction viability + non-fluor exemption; `gate_singlet` real `flowStats::singletGate` args migrated (`prediction_level`/`sidescatter`/`subsample_pct`/`maxit`; inert `upper_mult`/preset args removed → now rejected); `gate_flowclust_2d` promoted to a full FASTA plugin (`fcEllipsoidGate`, render-time polygon, full driver sourcing); `gate_mindensity` native `adjust`/`num_peaks`/`peaks` forwarded; `computePlotParams` perf (single-pass extraction, `gh_pop_get_count`, lazy-memo); `control_val` **`"x>y>z"` preference chain** + ambiguity guard; 1D threshold label restored under `.lg_override`; OOM-guarded `computePlotParams` binwidth.
- **`0_77`:** routine version bump, **no functional changes**.

> Note: the IVANKA fork also renamed its orchestration script `IVANKA_* → IVNCA_*` at `0_76` — that's an **IVANKA-project** change, not an engine change, and does not apply to this GUI project.

### FASTA Enhanced / FRISDA-era fixes (pre-sync, retained)

- **2026-06-14 — quad y-border value labels fixed (lazy-`aes` loop-var bug) + collapse caption.** The y-value threshold labels on quad gates were invisible because their position was computed **inside `aes()`** referencing the `.plot_ylim` **loop variable**, resolved lazily at **print** time → label offscreen. Fix: bake position + label into the data.frame. `gate_tmix_2d` (and its `copy_gate`) intentionally show NO threshold value label. Caption now shows `collapseDataForGating: TRUE` when set. *Lesson: never reference loop-local vars inside `aes()` when ggplots are built in a loop and printed later.*
- **2026-06-13 — `collapseDataForGating` implemented + QC/report fixes.** Pool the parent/`clusterFrom` events across the group, compute ONE gate, replicate. Reproducibly downsampled; mutually exclusive with external_control; skipped in QC identical-threshold checks. `threshold_outside_data_range` → restricted to `per_sample` gates. Gating QC tables paginate (`longtable=TRUE`+`repeat_header`).
- **2026-06-13 — first-PDF review fixes.** `gate_logic_count_mismatch` false-flagged every `NOT` gate (now reads parent indices); `gate_logic_component_missing` false-flagged every comma-string `gates=` (now parsed via `.fasta_parse_pop_list()`); `threshold_outside_data_range` band widened P1–P99.95; **unrecognized `gating_args` now error loudly** (`gate_static` exempt); quad x-border value labels un-clipped (`hjust=1`); strip text 7→6.
- **2026-06-13 — `clusterFrom` rejected on `gate_logic`/`copy_gate`** (neither does a data fit).
- **2026-06-13 — `parse_args` hardened** against malformed `gating_args` (now `stop()`s loudly with alias + offending string instead of silently returning `list()`; `pick_in` guards a `NULL`/empty control column).
- **2026-06-13 — layout_gates+gate_logic overlay CRASHED Positron (reverted).** Per-gate `gh_pop_get_data()` on index-based `gate_logic` dummy nodes inside the plot loop segfaulted. *Pending reimplementation must avoid per-gate `gh_pop_get_data` on gate_logic nodes — read counts/gates, not data.*
- **2026-06-13 — 1D split gate (`pop="+-"`/`"-+"`).** A single-channel rect row emits BOTH polarity populations from one shared threshold. Validated on E44 (Viability split → `LD+ 85.8%` / `LD- 14.2%`).
- **2026-06-13 — `layout_dims` / `layout_gates` wired into `apply_plotting`** + per-gate `layout_labels` for every `layout_gates` entry.
- **2026-06-13/12 — 2D rect overlay fixed (finite-edge `geom_segment`s, L-corner)**; First gating-layout QC + engine hardening (`.gate_static` errors on no-boundary; `apply_transformation` QC message).
- **2026-06-12 — FCS export args on `apply_compensation`** (`export_fcs`, `append_pdata_to_export`).
- **2026-06-12 — AutoSpill compensation implemented + validated (E44).** See Compensation section.
- **2026-06-11/12 — gt_gating → `apply_gating2` cutover.** Removed `apply_gating` (1481 lines) + `register_plugins` + dead `add_boolean`. `gate_logic` index-based; AND/OR/**NOT** + comma-string `gates=`. `gate_mindensity_2d` completed; `copy_gate` fully ported. Full-template parity confirmed (36 pops) before removal.
- **panel_QC** — root cause of mass mismatches was a **wrong `PanelIteration`**, not code; hardened anyway.
- **`CAR-` = `NOT(CAR_all)` gating bug** — `recompute()` wiped `CAR_all` to 100% before the `NOT` read it. Fixed by restoring prior gate_logic indices before reading components. *Lesson: parity proves engine agreement, not correctness.*

---

## Pending Work

**GUI build (active):**
1. **Build the GUI over the `0_77` backend.** GUI work lives in this folder (`FASTA_wGUI_0_77.R` + new GUI scripts), calling `ExtFunctions_0_77.R` as a read-only swappable dependency. Keep GUI logic independent of engine internals.
2. **Wire `PanelIteration`** (currently `"Placeholder"` in `FASTA_wGUI_0_77.R`) and confirm the live `instructions*.xlsx` (`20260623_instructions.xlsx`) sheets resolve before running `panel_QC`.

**Backend re-validation (advisory — engine is owned by the engine project):**
3. **Spot-check the `0_77` engine through the GUI** on representative data (gating `apply_gating2`, plotting `apply_plotting`, `panel_QC`, `run_gating_qc`). New `0_74`–`0_76` features (prevalidate, `gate_quantile`, `gate_mindensityshoulder`, `DotPlot`, visible-window binning, case-insensitive resolution) ride along — verify they behave as documented in this folder's data.
4. **Re-do the `layout_gates`+`gate_logic` overlay (crash-safe)** if needed — derive the crosshair from gate reads only, counts from `gh_pop_get_count`/indices, NOT per-gate `gh_pop_get_data` (see 2026-06-13 segfault). *Engine-side change — belongs in the engine project.*

**Engine-level (not active, tracked in the engine project):**
5. **`gate_flowmeans` (1D)** — only gate type never exercised in parity. Validate if ever used.
6. **`gate_hinge_quad`** — FlowJo-style tilted quadrant gate, not implemented.
7. **`autospectral` compensation** — still an unimplemented `apply_compensation` stub. (`autospill` is implemented.)
8. ✅ **DONE (2026-06-21) — `min_bins`/`max_bins` exposed on `apply_plotting`** (de-grain low-event panels; defaults 50/300; selective). See engine-sync summary.
9. `gate_flowclust_2d` — prevalidate does not catch a scalar `min`/`max` length (a length-2 check is a candidate enhancement); `upper_x`/`upper_y` ellipse-size multipliers deferred.
