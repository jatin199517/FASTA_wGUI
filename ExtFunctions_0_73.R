
# =============================================================================
# ExtFunctions.R
# FASTA Enhanced -- External Helper Functions
# =============================================================================
# Functions are organized by pipeline stage. Internal helpers (not intended
# for direct use) are prefixed with a dot (e.g., .comp_self).
# Each function includes a header comment block for troubleshooting.
# =============================================================================


# -----------------------------------------------------------------------------
# SECTION 1: COMPENSATION
# Module-level cache for external control gate thresholds
# Populated by apply_gating() before gt_gating() runs.
.fasta_gate_cache <- new.env(parent = emptyenv())

# -----------------------------------------------------------------------------

#' Apply compensation to a cytoset based on Compensation_type in metadata
#'
#' Dispatches each sample to the appropriate compensation handler based on
#' the Compensation_type column in the metadata sheet of the instructions file.
#'
#' Allowed Compensation_type values (per sample):
#'   "self"         - apply the acquisition-embedded spillover matrix ($SPILLOVER
#'                    preferred over $SPILL if both are present)
#'   "unmixed"      - no compensation applied (data already unmixed/compensated)
#'   NA or blank    - treated identically to "unmixed"; no compensation applied
#'   <filename>     - a .xlsx or .csv file in the working directory;
#'                    imported and applied as the spillover matrix.
#'                    Different samples may reference different files.
#'                    Files are loaded once and cached if referenced by multiple samples.
#'   "autospill"    - experiment-level AutoSpill spillover matrix, supplied via
#'                    the `autospill_matrix` argument (one matrix shared by all
#'                    "autospill" samples). Compute it beforehand with
#'                    build_autospill_controls() + run_autospill(). This branch
#'                    only APPLIES the matrix -- it never reads/writes files.
#'   "autospectral" - not yet implemented (requires spectral acquisition setup)
#'
#' @param cs       cytoset loaded via load_cytoset_from_fcs()
#' @param metadata data.frame with columns "fcs_filename" and "Compensation_type"
#' @param autospill_matrix square spillover matrix (channels x channels) to apply
#'   to every sample whose Compensation_type is "autospill". Required if any
#'   sample resolves to "autospill"; ignored otherwise. Row/col names must match
#'   colnames(cs) for the spanned channels.
#'
#' @param export_fcs logical (default FALSE). When TRUE, writes each compensated
#'   frame to a "Compensated/" folder in the working directory, named
#'   "<original_basename>_<comptype>.fcs" where <comptype> is the sample's
#'   resolved Compensation_type (e.g. "C1 117_autospill.fcs").
#' @param append_pdata_to_export logical (default FALSE). When TRUE *and*
#'   export_fcs is TRUE, that sample's pData row is appended to the exported
#'   frame's FCS keywords (original acquisition keywords + one keyword per pData
#'   column). No effect when export_fcs is FALSE.
#'
#' @return cytoset with compensation applied per sample
#'
#' @note For "self": if all samples use "self", compensate(cs, "data") is used
#'   for efficiency. All-"autospill" is likewise fast-pathed via
#'   compensate(cs, autospill_matrix). Mixed types are handled per-sample via fsApply.
#'
#' @examples
#'   cs_comp <- apply_compensation(cs, metadata = finalInstructions$metadata)
apply_compensation <- function(cs, metadata, autospill_matrix = NULL,
                               export_fcs = FALSE,
                               append_pdata_to_export = FALSE) {

  # -- 1. Input validation --------------------------------------------------
  if (!inherits(cs, "cytoset")) {
    stop("[apply_compensation] cs must be a cytoset object.")
  }

  required_cols <- c("fcs_filename", "Compensation_type")
  missing_cols  <- setdiff(required_cols, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop("[apply_compensation] metadata is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }

  # -- 2. Resolve compensation type per sample ------------------------------
  resolve_comp_type <- function(raw_type, sn) {
    if (is.na(raw_type) || trimws(raw_type) == "") return("none")
    raw_lower <- tolower(trimws(raw_type))
    if (raw_lower == "self")         return("self")
    if (raw_lower == "unmixed")      return("none")
    if (raw_lower == "autospill")    return("autospill")
    if (raw_lower == "autospectral") return("autospectral")
    ext <- tolower(tools::file_ext(raw_type))
    if (!ext %in% c("xlsx", "xls", "csv")) {
      stop("[apply_compensation] Unrecognized Compensation_type for sample ", sn,
           ": ", raw_type,
           "\n  Expected: self, unmixed, NA, or a .xlsx/.csv filename.")
    }
    if (!file.exists(raw_type)) {
      stop("[apply_compensation] Spillover file not found: ", raw_type,
           "\n  Referenced by sample: ", sn,
           "\n  Working directory: ", getwd(),
           "\n  Ensure the file is in the project folder.")
    }
    return(raw_type)
  }

  resolved <- vapply(sampleNames(cs), function(sn) {
    row <- metadata[metadata$fcs_filename == sn, ]
    if (nrow(row) == 0) {
      stop("[apply_compensation] No metadata row found for sample: ", sn,
           "\n  Check fcs_filename values match sampleNames(cs).")
    }
    if (nrow(row) > 1) {
      warning("[apply_compensation] Multiple metadata rows for ", sn, " -- using the first.")
      row <- row[1, ]
    }
    resolve_comp_type(row$Compensation_type, sn)
  }, character(1))

  message("[apply_compensation] Compensation type(s) in this run: ",
          paste(unique(resolved), collapse = ", "))

  # -- Optional export of compensated FCS -----------------------------------
  # Identity passthrough with a side effect: when export_fcs, write each frame
  # to "Compensated/" suffixed by its resolved Compensation_type. Wraps every
  # return so all compensation paths export consistently.
  do_export <- function(cs_out) {
    if (!isTRUE(export_fcs)) return(cs_out)
    out_dir <- file.path(getwd(), "Compensated")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    pd_all <- if (isTRUE(append_pdata_to_export)) pData(cs_out) else NULL
    for (sn in sampleNames(cs_out)) {
      fr <- cs_out[[sn]]
      if (!is.null(pd_all)) {
        # detach to a standalone flowFrame so we don't mutate cs_out by reference,
        # then append this sample's pData row as keywords (original kw preserved)
        fr  <- cytoframe_to_flowFrame(fr)
        row <- pd_all[sn, , drop = FALSE]
        for (col in colnames(row)) {
          val <- as.character(row[[col]])
          keyword(fr)[[col]] <- if (length(val) == 0 || is.na(val)) "" else val
        }
      }
      suffix <- gsub("[^A-Za-z0-9]+", "-", tools::file_path_sans_ext(resolved[[sn]]))
      base   <- sub("\\.fcs$", "", sn, ignore.case = TRUE)
      flowCore::write.FCS(fr, file.path(out_dir, paste0(base, "_", suffix, ".fcs")))
    }
    message("[apply_compensation] Exported ", length(cs_out),
            " compensated FCS to: ", out_dir,
            if (isTRUE(append_pdata_to_export)) " (pData appended to keywords)" else "")
    cs_out
  }

  # -- 3. Fast path: all samples are "self" ---------------------------------
  if (all(resolved == "self")) {
    message("[apply_compensation] All samples: applying acquisition-defined matrix.")
    cs_comp <- tryCatch(
      compensate(cs, "data"),
      error = function(e) stop(
        "[apply_compensation] Failed to apply acquisition-defined compensation.\n",
        "  Error: ", conditionMessage(e), "\n",
        "  Ensure all FCS files contain a $SPILL or $SPILLOVER keyword.\n",
        "  Confirm in FlowJo: Edit > Instrument Settings."
      )
    )
    message("[apply_compensation] Done -- compensated ", length(cs_comp), " sample(s).")
    return(do_export(cs_comp))
  }

  # -- 4. Fast path: all samples are "none" ---------------------------------
  if (all(resolved == "none")) {
    message("[apply_compensation] All samples: no compensation applied.")
    return(do_export(cs))
  }

  # -- 4b. autospill requires the experiment-level matrix -------------------
  # AutoSpill is experiment-level (one matrix from the single-color controls),
  # not per-sample. Compute it beforehand with build_autospill_controls() +
  # run_autospill() and pass it in as autospill_matrix. This function only
  # APPLIES it -- it never runs AutoSpill or touches files.
  if (any(resolved == "autospill")) {
    if (is.null(autospill_matrix)) {
      stop("[apply_compensation] Compensation_type 'autospill' requires a matrix.\n",
           "  Compute it first, then pass it in:\n",
           "    ctrl   <- build_autospill_controls(control_dir, finalInstructions$markers)\n",
           "    cmat   <- run_autospill(ctrl$control_dir, ctrl$control_def_file)\n",
           "    cs_comp <- apply_compensation(cs, metadata, autospill_matrix = cmat)")
    }
    if (!is.matrix(autospill_matrix) || nrow(autospill_matrix) != ncol(autospill_matrix)) {
      stop("[apply_compensation] autospill_matrix must be a square channels x channels matrix.")
    }
    miss <- setdiff(rownames(autospill_matrix), colnames(cs[[1]]))
    if (length(miss) > 0) {
      stop("[apply_compensation] autospill_matrix channel(s) not in cs: ",
           paste(miss, collapse = ", "),
           "\n  The matrix uses detector channel names ($PnN, e.g. FL25-A);",
           "\n  ensure these are present and unrenamed in cs.")
    }
    message("[apply_compensation] autospill: applying ", nrow(autospill_matrix),
            "-channel experiment-level matrix.")
  }

  # -- 4c. Fast path: all samples are "autospill" ---------------------------
  if (all(resolved == "autospill")) {
    cs_comp <- compensate(cs, autospill_matrix)
    message("[apply_compensation] Done -- compensated ", length(cs_comp), " sample(s).")
    return(do_export(cs_comp))
  }

  # -- 5. Per-sample path: mixed or file-based types ------------------------
  file_cache <- list()

  load_spill_file <- function(filename) {
    if (!is.null(file_cache[[filename]])) return(file_cache[[filename]])
    ext <- tolower(tools::file_ext(filename))
    spill_df <- if (ext == "csv") {
      read.csv(filename, row.names = 1, check.names = FALSE)
    } else {
      raw <- as.data.frame(readxl::read_excel(filename, col_names = TRUE))
      rownames(raw) <- raw[[1]]
      raw[[1]] <- NULL
      raw
    }
    spill_mat <- as.matrix(spill_df)
    if (nrow(spill_mat) != ncol(spill_mat)) {
      stop("[apply_compensation] Spillover matrix in ", filename, " is not square (",
           nrow(spill_mat), " x ", ncol(spill_mat), ").\n",
           "  Row names and column names must list the same channel set.")
    }
    file_cache[[filename]] <<- spill_mat
    message("  Loaded spillover matrix: ", filename, " -- ", nrow(spill_mat), " channels")
    return(spill_mat)
  }

  comp_fs <- fsApply(cs, function(ff) {
    sn   <- identifier(ff)
    type <- resolved[sn]
    switch(type,
      "none" = {
        message("  [unmixed/NA] ", sn)
        ff
      },
      "self" = {
        message("  [self] ", sn)
        spill_list <- tryCatch(spillover(ff), error = function(e) NULL)
        spill <- spill_list[["$SPILLOVER"]] %||%
                 spill_list[["$SPILL"]]     %||%
                 spill_list[["SPILL"]]      %||%
                 Filter(Negate(is.null), spill_list)[[1]]
        if (is.null(spill)) {
          stop("[apply_compensation] [self] No spillover matrix in FCS header for: ", sn,
               "\n  Expected $SPILLOVER or $SPILL keyword.",
               "\n  Confirm in FlowJo: Edit > Instrument Settings.",
               "\n  Alternatively, supply a .xlsx/.csv filename in Compensation_type.")
        }
        compensate(ff, spill)
      },
      "autospill" = {
        message("  [autospill] ", sn)
        compensate(ff, autospill_matrix)
      },
      "autospectral" = stop(
        "[apply_compensation] [autospectral] Not yet implemented for: ", sn,
        "\n  Export unmixed data and use unmixed instead."
      ),
      { # Default: filename
        message("  [file: ", type, "] ", sn)
        spill_mat <- load_spill_file(type)
        tryCatch(
          compensate(ff, spill_mat),
          error = function(e) stop(
            "[apply_compensation] Failed to apply spillover from ", type,
            " to sample: ", sn,
            "\n  Error: ", conditionMessage(e),
            "\n  Channel names in the matrix must exactly match colnames(cs[[1]]).",
            "\n  Run: setdiff(rownames(spill_mat), colnames(cs[[1]])) to find mismatches."
          )
        )
      }
    )
  })

  cs_comp <- flowSet_to_cytoset(comp_fs)
  message("[apply_compensation] Done -- compensated ", length(cs_comp), " sample(s).")
  return(do_export(cs_comp))
}


# -- Null-coalescing operator ------------------------------------------------
`%||%` <- function(lhs, rhs) if (!is.null(lhs)) lhs else rhs


# -----------------------------------------------------------------------------
# SECTION 1b: AutoSpill compensation (experiment-level)
# -----------------------------------------------------------------------------
# AutoSpill (carlosproca/autospill) computes ONE experiment-level spillover
# matrix from the set of single-color controls. The two helpers below stage the
# controls (build_autospill_controls) and run the AutoSpill pipeline
# (run_autospill). The resulting matrix is then handed to apply_compensation()
# via its autospill_matrix argument. Writing the matrix to disk (the
# "BackupCompMatrix" CSV) is a script-level step, kept OUT of these functions.
# -----------------------------------------------------------------------------

#' Stage single-color controls for AutoSpill from a vetted dye/antigen table
#'
#' The dye<->antigen pairing is taken AS VETTED from the markers sheet
#' (dye = `dims`, antigen = `markers`) -- the method does NOT guess channels.
#' For each vetted control the FCS filename is auto-populated from `controls` by
#' matching the control's fluorophore to `controls$SampleID` (normalised; with
#' `alias` for the few tokens that differ, e.g. PB <-> Pacific Blue).
#'
#' E44 controls are SpectraComp/ViaComp beads, each carrying both a positive and
#' a negative bead population in one file, so they are used AS-IS (symlinked) --
#' AutoSpill gates pos/neg within each file. (AutoSpill has no native per-marker
#' separate-negative input, so there is nothing to concatenate for beads.)
#'
#' GFP is the exception: if BOTH a positive and negative GFP file exist they are
#' concatenated into one frame so AutoSpill sees both populations; if only ONE
#' GFP file is provided it is assumed to already contain both populations and
#' used as-is, like the beads.
#'
#' Symlinks + the (optional) GFP concat go in a temp dir; nothing is written to
#' the data folder. The control definition is RETURNED as a data.frame (the
#' fcs_control equivalent) -- run_autospill() materialises a throwaway tempfile
#' from it (AutoSpill's readers hardcode read.csv() on a path, read twice, so a
#' dataframe can't be passed directly without forking it).
#'
#' @param markers     markers sheet (finalInstructions$markers): `dims` (detector
#'                    channel = dye) + `markers` (fluorophore_antigen). The vetted
#'                    source of dye<->antigen.
#' @param controls    data.frame of single-stain controls with columns `filename`
#'                    (FCS basename in control_dir) and `SampleID` (fluorophore
#'                    token; also identifies the GFP files).
#' @param control_dir directory holding the raw control FCS files.
#' @param alias       named char vector: SampleID token -> markers fluorophore,
#'                    for tokens that differ between the two.
#' @param gfp_pattern regex identifying GFP control(s) by antigen / SampleID.
#' @param out_dir     temp staging dir (recreated fresh).
#' @param verbose     print the mapping table for sign-off.
#'
#' @return list(control_dir, control_def, mapping). `control_def` is the
#'   filename/dye/antigen/wavelength data.frame for run_autospill().
build_autospill_controls <- function(markers, controls, control_dir,
                                     alias = c("PB" = "Pacific Blue",
                                               "eFluor780" = "Viability"),
                                     gfp_pattern = "GFP",
                                     out_dir = file.path(tempdir(), "autospill_ctrl"),
                                     verbose = TRUE) {
  if (length(control_dir) != 1)
    stop("[build_autospill_controls] control_dir must be a single directory, got ",
         length(control_dir), ": ", paste(control_dir, collapse = ", "),
         "\n  (single-color controls span >1 folder, or a path is NA -- check metadata$path).")
  if (!dir.exists(control_dir))
    stop("[build_autospill_controls] control_dir not found: ", control_dir)
  for (cc in c("filename", "SampleID"))
    if (!cc %in% colnames(controls))
      stop("[build_autospill_controls] controls is missing column: ", cc)

  norm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

  # -- vetted dye <-> antigen from the markers sheet -------------------------
  mk <- as.data.frame(markers)
  mk <- mk[!is.na(mk$markers), c("dims", "markers")]
  vetted <- data.frame(dye = mk$dims,
                       antigen = sub("-A$", "", mk$markers),
                       stringsAsFactors = FALSE)
  vetted$fluor <- sub("_.*$", "", vetted$antigen)   # fluorophore = prefix before "_"

  # -- control SampleIDs, alias-normalised, for matching ---------------------
  ctl <- as.data.frame(controls)[, c("filename", "SampleID")]
  ctl$key <- norm(ifelse(ctl$SampleID %in% names(alias),
                         alias[ctl$SampleID], ctl$SampleID))

  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)
  dir.create(out_dir, recursive = TRUE)
  stage_asis <- function(f)
    file.symlink(normalizePath(file.path(control_dir, f)), file.path(out_dir, f))

  rows <- list()
  for (i in seq_len(nrow(vetted))) {
    dye <- vetted$dye[i]; antigen <- vetted$antigen[i]; nfluor <- norm(vetted$fluor[i])

    if (grepl(gfp_pattern, antigen, ignore.case = TRUE)) {
      gfp_rows <- ctl[grepl(gfp_pattern, ctl$SampleID, ignore.case = TRUE), ]
      if (nrow(gfp_rows) == 0)
        stop("[build_autospill_controls] no GFP control file found.")
      if (nrow(gfp_rows) >= 2) {              # pos + neg provided -> concatenate
        out_name <- "GFP single stain pos_neg.fcs"
        pos <- gfp_rows$filename[grepl("pos|\\+", gfp_rows$SampleID, ignore.case = TRUE)][1]
        neg <- gfp_rows$filename[grepl("neg|-",  gfp_rows$SampleID, ignore.case = TRUE)][1]
        if (is.na(pos) || is.na(neg)) { pos <- gfp_rows$filename[1]; neg <- gfp_rows$filename[2] }
        p <- flowCore::read.FCS(file.path(control_dir, pos), transformation = NULL)
        n <- flowCore::read.FCS(file.path(control_dir, neg), transformation = NULL)
        flowCore::exprs(p) <- rbind(flowCore::exprs(p), flowCore::exprs(n))
        suppressWarnings(flowCore::write.FCS(p, file.path(out_dir, out_name)))
        src <- paste0(pos, " + ", neg)
      } else {                                # single GFP file -> assume pos/neg mix, use as-is
        out_name <- gfp_rows$filename[1]
        stage_asis(out_name)
        src <- "single GFP file (as-is)"
      }
      rows[[length(rows) + 1]] <- data.frame(filename = out_name, dye = dye,
        antigen = antigen, source = src, stringsAsFactors = FALSE)
      next
    }

    # bead control: matched by fluorophore, used as-is (has its own pos+neg)
    pos <- ctl$filename[ctl$key == nfluor]
    if (length(pos) != 1)
      stop("[build_autospill_controls] vetted antigen '", antigen, "' (fluor '",
           vetted$fluor[i], "') matched ", length(pos), " controls by SampleID. ",
           "Add an alias if the SampleID token differs.")
    stage_asis(pos)
    rows[[length(rows) + 1]] <- data.frame(filename = pos, dye = dye,
      antigen = antigen, source = "as-is", stringsAsFactors = FALSE)
  }
  mapping <- do.call(rbind, rows)

  if (anyDuplicated(mapping$dye))
    stop("[build_autospill_controls] two controls map to the same dye: ",
         paste(mapping$dye[duplicated(mapping$dye)], collapse = ", "))

  control_def <- data.frame(
    filename   = mapping$filename,
    dye        = mapping$dye,
    antigen    = mapping$antigen,
    wavelength = NA,
    stringsAsFactors = FALSE)

  if (verbose) {
    message("[build_autospill_controls] ", nrow(mapping),
            " controls staged (beads as-is; GFP concatenated only if pos+neg)")
    message("  Control dir: ", out_dir)
    print(mapping)
  }
  invisible(list(control_dir = out_dir, control_def = control_def, mapping = mapping))
}


#' Run the AutoSpill pipeline and return the spillover matrix
#'
#' Thin wrapper over the AutoSpill `minimal` workflow (which writes no figures
#' or tables to disk). Overrides the scatter parameters to the acquisition's
#' raw channel names, since AutoSpill reads the raw FCS (pre-rename) and its
#' default c("FSC-A","SSC-A") will not match instruments that name scatter
#' FS00-A / SS02-A.
#'
#' @param control_dir   staging dir from build_autospill_controls()
#' @param control_def    control definition: either the data.frame returned by
#'   build_autospill_controls()$control_def (preferred -- the dataframe is the
#'   source of truth; a throwaway tempfile is written internally only because
#'   AutoSpill's readers hardcode read.csv() on a path), or a path to an existing
#'   fcs_control.csv.
#' @param scatter_param raw scatter channel names used for gating
#' @param param_set     AutoSpill parameter set (default "minimal")
#' @param verbose       pass-through to AutoSpill's iteration logging
#'
#' @return square spillover matrix (detector channels x detector channels),
#'   diagonal == 1; row/col names are detector channels ($PnN). Hand to
#'   apply_compensation(..., autospill_matrix = <this>).
run_autospill <- function(control_dir, control_def,
                          scatter_param = c("FS00-A", "SS02-A"),
                          param_set = "minimal", verbose = TRUE) {
  if (!requireNamespace("autospill", quietly = TRUE))
    stop("[run_autospill] package 'autospill' not installed ",
         "(devtools::install_github('carlosproca/autospill')).")

  # AutoSpill's read.marker()/read.flow.control() both read.csv() the control
  # definition (twice), so a dataframe can't be passed directly. Serialise it to
  # an ephemeral tempfile for this call only; never a persistent local CSV.
  if (is.data.frame(control_def)) {
    control_def_file <- tempfile("fcs_control_", fileext = ".csv")
    on.exit(unlink(control_def_file), add = TRUE)
    utils::write.csv(control_def, control_def_file, row.names = FALSE)
  } else {
    control_def_file <- control_def
    if (!file.exists(control_def_file))
      stop("[run_autospill] control_def is neither a data.frame nor an ",
           "existing file path: ", control_def_file)
  }

  asp <- autospill::get.autospill.param(param_set)
  asp$default.scatter.parameter <- scatter_param
  asp$verbose <- isTRUE(verbose)

  flow.control <- autospill::read.flow.control(control_dir, control_def_file, asp)
  flow.gate    <- autospill::gate.flow.data(flow.control, asp)
  sp.unco.untr <- autospill::get.marker.spillover(TRUE, flow.gate, flow.control, asp)
  refined      <- autospill::refine.spillover(sp.unco.untr, NULL, flow.gate,
                                              flow.control, asp)
  message("[run_autospill] converged -- ", nrow(refined$spillover),
          "-channel spillover matrix.")
  refined$spillover
}


# -----------------------------------------------------------------------------
# SECTION 2: TRANSFORMATION
# -----------------------------------------------------------------------------

#' Apply per-channel transformation to a compensated cytoset
#'
#' Reads TransformType and TransformArgs from the markers sheet and
#' applies the appropriate transform to each channel across all samples.
#'
#' Allowed TransformType values (per channel in markers sheet):
#'   NA or blank  - no transform
#'   "flowjo"     - flowjo biexponential via flowjo_biexp()
#'                  defaults: channelRange=9000, maxValue=522144,
#'                            pos=4.5, neg=-2.5, widthBasis=-50
#'   "asinh"      - arcsinhTransform(); defaults: a=0, b=0.05 (cofactor=20), c=0
#'   "biexp"      - data-driven via estimateLogicle() per sample; default: m=6
#'
#' TransformArgs (per channel, optional): semicolon-separated key=value pairs.
#'   Examples: "widthBasis=-100" | "a=0;b=0.002;c=0" | "m=5;t=262144"
#'
#' @param cs_comp  compensated cytoset (output of apply_compensation())
#' @param markers  data.frame with columns: dims, TransformType,
#'                 and optionally TransformArgs
#'
#' @return transformed cytoset
#'
#' @examples
#'   cs_comp_trans <- apply_transformation(cs_comp, markers = finalInstructions$markers)
apply_transformation <- function(cs_comp, markers, force = FALSE) {

  # -- 1. Input validation --------------------------------------------------
  if (!inherits(cs_comp, "cytoset")) {
    stop("[apply_transformation] cs_comp must be a cytoset object.")
  }

  # Guard against applying transformation multiple times.
  # Each run of apply_transformation tags the cytoset with an attribute.
  # Re-running on already-transformed data squishes the signal (biexp-on-biexp).
  if (!force && isTRUE(attr(cs_comp, "FASTA_transformed"))) {
    warning("[apply_transformation] This cytoset has already been transformed.
",
            "  Returning unchanged. Use force=TRUE to re-apply (not recommended).
",
            "  Tip: rebuild cs_comp_trans from cs_comp, not from cs_comp_trans.")
    return(cs_comp)
  }

  required_cols <- c("dims", "TransformType")
  missing_cols  <- setdiff(required_cols, colnames(markers))
  if (length(missing_cols) > 0) {
    stop("[apply_transformation] markers sheet is missing column(s): ",
         paste(missing_cols, collapse = ", "),
         "\n  Ensure dims and TransformType exist in the markers sheet.")
  }

  has_args_col <- "TransformArgs" %in% colnames(markers)

  # -- QC: marker channels left untransformed (no TransformType) -------------
  # Defined here, emitted before every return so it also fires on the early
  # "nothing to transform" path. Lists marker NAMES (falls back to channels if
  # the markers sheet has no `markers` column); scatter/time (no marker) excluded.
  emit_untransformed_qc <- function() {
    tt_blank <- is.na(markers$TransformType) | trimws(markers$TransformType) == ""
    mk_names <- if ("markers" %in% colnames(markers))
        markers$markers[tt_blank & !is.na(markers$markers)]
      else markers$dims[tt_blank]
    if (length(mk_names) > 0)
      message("[apply_transformation] QC -- ", length(mk_names),
              " marker channel(s) left UNTRANSFORMED (TransformType NA/blank): ",
              paste(mk_names, collapse = ", "),
              "\n  If unintended, check the TransformType column (e.g. 'flowjo' ",
              "placed in TransformArgs instead of TransformType).")
  }

  valid_types <- c("flowjo", "biexp", "asinh")
  raw_types   <- na.omit(unique(tolower(trimws(markers$TransformType))))
  raw_types   <- raw_types[raw_types != ""]
  bad_types   <- setdiff(raw_types, valid_types)
  if (length(bad_types) > 0) {
    stop("[apply_transformation] Unrecognized TransformType value(s): ",
         paste(bad_types, collapse = ", "),
         "\n  Valid options: flowjo, biexp, asinh, NA, or blank.")
  }

  # -- 2. Helpers -----------------------------------------------------------
  parse_args <- function(args_str) {
    if (is.na(args_str) || trimws(args_str) == "") return(list())
    pairs <- trimws(strsplit(args_str, ";")[[1]])
    pairs <- pairs[pairs != ""]
    parsed <- lapply(pairs, function(p) {
      kv <- strsplit(p, "=")[[1]]
      if (length(kv) != 2) {
        stop("[apply_transformation] Cannot parse TransformArgs entry: ", p,
             "\n  Expected format: key=value (e.g., m=6 or a=0;b=0.05;c=0)")
      }
      val <- suppressWarnings(as.numeric(trimws(kv[2])))
      if (is.na(val)) {
        stop("[apply_transformation] Non-numeric value in TransformArgs: ", p,
             "\n  All argument values must be numeric.")
      }
      setNames(list(val), trimws(kv[1]))
    })
    do.call(c, parsed)
  }

  get_channels <- function(type) {
    idx <- !is.na(markers$TransformType) &
           tolower(trimws(markers$TransformType)) == type
    markers$dims[idx]
  }

  get_args <- function(channel) {
    if (!has_args_col) return(list())
    parse_args(markers$TransformArgs[markers$dims == channel])
  }

  # -- 3. Group channels by type --------------------------------------------
  flowjo_channels <- get_channels("flowjo")
  asinh_channels  <- get_channels("asinh")
  biexp_channels  <- get_channels("biexp")

  if (length(c(flowjo_channels, asinh_channels, biexp_channels)) == 0) {
    emit_untransformed_qc()
    message("[apply_transformation] No channels marked for transformation. Returning unchanged.")
    return(cs_comp)
  }

  cs_trans <- cs_comp

  # -- 4. flowjo transform --------------------------------------------------
  # Defaults: channelRange=9000, maxValue=522144, pos=4.5, neg=-2.5, widthBasis=-50
  if (length(flowjo_channels) > 0) {
    message("[apply_transformation] Applying flowjo biexp to ",
            length(flowjo_channels), " channel(s): ",
            paste(flowjo_channels, collapse = ", "))
    flowjo_trans_list <- lapply(flowjo_channels, function(ch) {
      args     <- get_args(ch)
      defaults <- list(channelRange = 9000, maxValue = 522144,
                       pos = 4.5, neg = -2.5, widthBasis = -50)
      params   <- modifyList(defaults, args)
      tryCatch(
        do.call(flowjo_biexp, params),
        error = function(e) stop(
          "[apply_transformation] [flowjo] Failed for channel: ", ch,
          "\n  Error: ", conditionMessage(e),
          "\n  Check TransformArgs for this channel in the markers sheet."
        )
      )
    })
    flowjo_tl <- transformList(flowjo_channels, flowjo_trans_list)
    cs_trans  <- transform(cs_trans, flowjo_tl)
  }

  # -- 5. asinh transform ---------------------------------------------------
  # Default: a=0, b=0.05 (cofactor=20), c=0
  if (length(asinh_channels) > 0) {
    message("[apply_transformation] Applying asinh to ",
            length(asinh_channels), " channel(s): ",
            paste(asinh_channels, collapse = ", "))
    asinh_trans_list <- lapply(asinh_channels, function(ch) {
      args     <- get_args(ch)
      defaults <- list(a = 0, b = 0.05, c = 0)
      params   <- modifyList(defaults, args)
      tryCatch(
        arcsinhTransform(transformationId = paste0("asinh_", ch),
                         a = params$a, b = params$b, c = params$c),
        error = function(e) stop(
          "[apply_transformation] [asinh] Failed for channel: ", ch,
          "\n  Error: ", conditionMessage(e),
          "\n  Expected TransformArgs format: a=<num>;b=<num>;c=<num>"
        )
      )
    })
    asinh_tl <- transformList(asinh_channels, asinh_trans_list)
    cs_trans <- transform(cs_trans, asinh_tl)
  }

  # -- 6. biexp transform (estimateLogicle, per-sample) ---------------------
  # Default: m=6. Fails gracefully per sample -- untransformed sample is
  # returned and a summary warning is issued at the end.
  if (length(biexp_channels) > 0) {
    message("[apply_transformation] Applying estimateLogicle (biexp) to ",
            length(biexp_channels), " channel(s): ",
            paste(biexp_channels, collapse = ", "))

    biexp_params <- modifyList(list(m = 6), get_args(biexp_channels[1]))
    biexp_errors <- character(0)

    comp_fs <- fsApply(cs_trans, function(ff) {
      sn <- identifier(ff)
      tryCatch({
        trans <- do.call(estimateLogicle,
                         c(list(x = ff, channels = biexp_channels), biexp_params))
        transform(ff, trans)
      },
      error = function(e) {
        biexp_errors <<- c(biexp_errors, paste0("  - ", sn, ": ", conditionMessage(e)))
        message("[apply_transformation] [biexp] WARNING -- estimateLogicle failed for: ", sn,
                "\n  ", conditionMessage(e),
                "\n  Returning this sample untransformed.")
        ff
      })
    })

    if (length(biexp_errors) > 0) {
      warning(
        "[apply_transformation] [biexp] estimateLogicle failed for ",
        length(biexp_errors), " sample(s):\n",
        paste(biexp_errors, collapse = "\n"),
        "\n  Possible fixes:\n",
        "  1. Reduce m (e.g., m=4.5 or m=5) in TransformArgs.\n",
        "  2. Check samples for sparse or zero-skewed channel distributions.\n",
        "  3. Switch affected channels to flowjo or asinh.",
        call. = FALSE
      )
    }

    cs_trans <- flowSet_to_cytoset(comp_fs)
    pData(cs_trans) <- pData(cs_comp)[sampleNames(cs_trans), , drop=FALSE]
  }

  message("[apply_transformation] Done -- transformed ", length(cs_trans), " sample(s).")
  emit_untransformed_qc()
  # General pData guard — restore any columns lost through fsApply paths
  if (!identical(sort(colnames(pData(cs_comp))), sort(colnames(pData(cs_trans))))) {
    pData(cs_trans) <- pData(cs_comp)[sampleNames(cs_trans), , drop=FALSE]
  }
  attr(cs_trans, "FASTA_transformed") <- TRUE
  return(cs_trans)
}



# -----------------------------------------------------------------------------
# SECTION 3: PLOTTING UTILITIES
# -----------------------------------------------------------------------------

#' Compute axis limits and binwidth for a ggcyto plot from a GatingSet
#'
#' Calculates per-channel plot limits at the 0.5th/99.5th percentiles across
#' all samples in a GatingSet, adds 10% padding, and derives a square binwidth
#' suitable for geom_hex(). Intended to be called once per plot panel and the
#' result passed directly to coord_cartesian() and geom_hex().
#'
#' @param gs        GatingSet object
#' @param pop       character; gating population to extract data from
#'                  (e.g., "root", "Single", "CD3+"). Must exist in gs.
#' @param xchannel  character; channel name for the x-axis (must match colnames of the data)
#' @param ychannel  character; channel name for the y-axis
#' @param bins      numeric; adaptive bin count used to size the binwidth.
#'                  Typically the output of adaptiveBins(). Controls resolution
#'                  of geom_hex() -- higher values give finer bins.
#'
#' @return named list with three elements:
#'   $xlim     numeric(2) -- x-axis limits with 10% padding
#'   $ylim     numeric(2) -- y-axis limits with 10% padding
#'   $binwidth numeric(2) -- square binwidth (same value for x and y)
#'
#' @examples
#'   params <- computePlotParams(gs, pop = "root", xchannel = "FSC-A",
#'                               ychannel = "SSC-A", bins = 128)
#'   ggcyto(gs, aes(x = "FSC-A", y = "SSC-A"), subset = "root") +
#'     geom_hex(binwidth = params$binwidth) +
#'     coord_cartesian(xlim = params$xlim, ylim = params$ylim)
computePlotParams <- function(gs, pop, xchannel, ychannel, bins) {

  # -- 1. Input validation --------------------------------------------------
  if (!inherits(gs, "GatingSet")) {
    stop("[computePlotParams] gs must be a GatingSet object.")
  }

  valid_pops <- gs_get_pop_paths(gs)

  if (!pop %in% valid_pops) {
    matches <- valid_pops[basename(valid_pops) == pop]
    if (length(matches) == 0) {
      stop("[computePlotParams] Population not found: ", pop,
           "\n  Available populations:\n  ",
           paste(valid_pops, collapse = "\n  "))
    }
    if (length(matches) > 1) {
      warning("[computePlotParams] Multiple populations match short name ", pop,
              " -- using: ", matches[1])
    }
    pop <- matches[1]
  }

  # Check channels exist in at least the first sample
  sample_cols <- colnames(gh_pop_get_data(gs[[1]], pop))
  bad_channels <- setdiff(c(xchannel, ychannel), sample_cols)
  if (length(bad_channels) > 0) {
    stop("[computePlotParams] Channel(s) not found in population data: ",
         paste(bad_channels, collapse = ", "),
         "\n  Available channels: ", paste(sample_cols, collapse = ", "))
  }

  if (!is.numeric(bins) || length(bins) != 1 || bins <= 0) {
    stop("[computePlotParams] bins must be a single positive number.")
  }

  # -- 2. Compute per-channel limits across all samples ---------------------
  # Uses 0.5th / 99.5th percentile to exclude extreme outliers while keeping
  # biologically relevant signal in frame.
  get_range <- function(channel) {
    vals <- tryCatch(
      sapply(seq_along(gs), function(i)
        quantile(exprs(gh_pop_get_data(gs[[i]], pop))[, channel], c(0.005, 0.995))
      ),
      error = function(e) stop(
        "[computePlotParams] Failed to extract data for channel: ", channel,
        "\n  Population: ", pop,
        "\n  Error: ", conditionMessage(e)
      )
    )
    c(min(vals[1, ]), max(vals[2, ]))
  }

  xlim_raw <- get_range(xchannel)
  ylim_raw <- get_range(ychannel)

  # -- 3. Add 10% padding on each side -------------------------------------
  xlim <- xlim_raw + c(-1, 1) * 0.1 * diff(xlim_raw)
  ylim <- ylim_raw + c(-1, 1) * 0.1 * diff(ylim_raw)

  # -- 4. Derive square binwidth from the narrower axis --------------------
  # Binwidth: independent per axis using 5th-95th percentile range.
  # This ensures square bins sized to the actual data region, not to extreme outliers.
  get_dense_range <- function(channel) {
    vals <- sapply(seq_along(gs), function(i)
      quantile(exprs(gh_pop_get_data(gs[[i]], pop))[, channel], c(0.05, 0.95)))
    max(vals[2,]) - min(vals[1,])
  }
  dense_x <- get_dense_range(xchannel)
  dense_y <- get_dense_range(ychannel)
  binwidth <- c(dense_x / bins, dense_y / bins)

  list(xlim = xlim, ylim = ylim, binwidth = binwidth)
}


#' Flow cytometry colour-scaled 2D bin plot (drop-in for geom_hex)
#'
#' Returns a list of ggplot layers that adds a geom_bin2d with a log10-scaled
#' flow cytometry colour ramp (navy -> blue -> cyan -> green -> yellow -> red).
#' Empty bins are rendered as dark navy. The list adds cleanly to a ggplot or
#' ggcyto chain with +.
#'
#' Use bins= for most plots. For panels where x and y data ranges differ
#' substantially (e.g. a DUMP gate), pass binwidth= instead as a length-2
#' vector c(x_width, y_width) in data units -- this keeps bins visually square
#' and prevents horizontal stripe artefacts from a narrow y-range.
#'
#' @param bins      numeric scalar; number of bins in each dimension.
#'                  Typically from adaptiveBins(). Ignored if binwidth is set.
#' @param binwidth  numeric(2); bin width in data units c(x_width, y_width).
#'                  Use instead of bins when axis ranges differ substantially.
#'
#' @return list of ggplot layers (geom_bin2d + scale_fill_gradientn)
#'
#' @examples
#'   ggcyto(gs[i], aes(x = "FSC-A", y = "SSC-A"), subset = "root") +
#'     flow_bin2d(bins = adaptiveBins(gs, "root"))
#'
#'   # For asymmetric ranges (e.g. DUMP gate):
#'   ggcyto(gs[i], aes(x = "FL1-A", y = "SSC-A"), subset = "Single") +
#'     flow_bin2d(binwidth = params$binwidth)
flow_bin2d <- function(bins = NULL, binwidth = NULL) {

  if (is.null(bins) && is.null(binwidth)) {
    stop("[flow_bin2d] Provide either bins= (scalar) or binwidth= (length-2 vector).")
  }
  if (!is.null(bins) && !is.null(binwidth)) {
    stop("[flow_bin2d] Provide either bins= or binwidth=, not both.")
  }
  if (!is.null(bins) && (!is.numeric(bins) || length(bins) != 1 || bins <= 0)) {
    stop("[flow_bin2d] bins must be a single positive number.")
  }
  if (!is.null(binwidth) && (!is.numeric(binwidth) || length(binwidth) != 2)) {
    stop("[flow_bin2d] binwidth must be a numeric vector of length 2: c(x_width, y_width).")
  }

  bin_args <- if (!is.null(binwidth)) list(binwidth = binwidth) else list(bins = bins)

  list(
    do.call(geom_bin2d, bin_args),
    scale_fill_gradientn(
      colours  = c("navy", "#1E90FF", "cyan", "#32CD32", "yellow", "orange", "red"),
      trans    = "log10",
      na.value = "navy",
      name     = "Count"
    )
  )
}


#' Compute adaptive bin count from median event count across a GatingSet
#'
#' Estimates an appropriate bin count for geom_bin2d / geom_hex by taking the
#' square root of the median event count across all samples in a population.
#' The result is clamped to [min_bins, max_bins] to avoid degenerate plots on
#' very sparse or very large samples.
#'
#' @param gs        GatingSet object
#' @param pop       character; gating population to count events from
#' @param min_bins  numeric; minimum bin count (default 50)
#' @param max_bins  numeric; maximum bin count (default 300)
#'
#' @return single integer bin count, clamped to [min_bins, max_bins]
#'
#' @examples
#'   bins <- adaptiveBins(gs, pop = "Single")
#'   flow_bin2d(bins = bins)
#'   computePlotParams(gs, pop = "Single", xchannel = "FL1-A",
#'                     ychannel = "FL2-A", bins = bins)
adaptiveBins <- function(gs, pop, min_bins = 50, max_bins = 300) {

  # -- Input validation -----------------------------------------------------
  if (!inherits(gs, "GatingSet")) {
    stop("[adaptiveBins] gs must be a GatingSet object.")
  }

  valid_pops <- gs_get_pop_paths(gs)

  # Accept either full path ("/Cells/Singlets") or short name ("Singlets")
  if (!pop %in% valid_pops) {
    matches <- valid_pops[basename(valid_pops) == pop]
    if (length(matches) == 0) {
      stop("[adaptiveBins] Population not found: ", pop,
           "\n  Available populations:\n  ",
           paste(valid_pops, collapse = "\n  "))
    }
    if (length(matches) > 1) {
      warning("[adaptiveBins] Multiple populations match short name ", pop,
              " -- using: ", matches[1],
              "\n  Pass the full path to be explicit.")
    }
    pop <- matches[1]
  }

  if (!is.numeric(min_bins) || !is.numeric(max_bins) || min_bins <= 0 || max_bins <= 0) {
    stop("[adaptiveBins] min_bins and max_bins must be positive numbers.")
  }
  if (min_bins > max_bins) {
    stop("[adaptiveBins] min_bins (", min_bins, ") must be <= max_bins (", max_bins, ").")
  }

  # -- Compute ---------------------------------------------------------------
  n_events <- sapply(seq_along(gs), function(i) {
    nrow(exprs(gh_pop_get_data(gs[[i]], pop)))
  })

  if (all(n_events == 0)) {
    warning("[adaptiveBins] All samples have 0 events in population: ", pop,
            " -- returning min_bins (", min_bins, ").")
    return(as.integer(min_bins))
  }

  bins <- round(sqrt(median(n_events)))
  as.integer(max(min_bins, min(max_bins, bins)))
}


# -----------------------------------------------------------------------------
# SECTION 4: GATING
# -----------------------------------------------------------------------------

#' Apply automated gating to a GatingSet using an openCyto gating template
#'
#' Reads a gating template CSV file and applies it to a GatingSet using
#' openCyto's gt_gating(). The template drives all gate definitions --
#' gating method, parent population, channel(s), and optional arguments are
#' all read from the file, making the pipeline fully data-driven.
#'
#' Template columns (openCyto standard):
#'   alias       - name of the resulting population
#'   pop         - polarity: + (positive) or - (negative)
#'   parent      - parent population the gate is applied to
#'   dims        - channel or marker name(s); comma-separated for 2D gates
#'   gating_method     - openCyto gating function (e.g. gate_mindensity,
#'                       singletGate, gate_quantile). Custom methods can be
#'                       registered via openCyto::register_plugins().
#'   gating_args       - named args passed to the gating function as a
#'                       quoted string (e.g. "min=0,max=5000"). NA = defaults.
#'   collapseDataForGating - TRUE to POOL the parent (or clusterFrom) events
#'                           across the group (or all samples if groupBy blank),
#'                           compute ONE gate, and replicate it to every sample
#'                           in the group (useful for low-count populations).
#'                           Pooling is reproducibly downsampled (collapse_max_events
#'                           per group, collapse_seed). Mutually exclusive with
#'                           external_control. NA/FALSE = per-sample.
#'   groupBy           - grouping for collapse/external_control: a pData column
#'                       (or ":"/"," combo, e.g. "Batch:Day"), or a number N to
#'                       group every N samples. NA/blank = no grouping (all one).
#'   preprocessing_method - optional preprocessing function name. NA = none.
#'   preprocessing_args   - args for preprocessing function. NA = none.
#'
#' @param gs            GatingSet object (modified in place by openCyto)
#' @param template_path character; path to the gating template CSV file.
#'                      Defaults to "GatingTemplate.csv" in the working directory.
#' @param mc.cores      integer; number of cores for parallel gating (default 1).
#'                      Increase for large datasets if parallel backend is available.
#'
#' @return the GatingSet gs (invisibly), modified in place with gates added
#'
#' @note GatingSets are reference objects -- gs is modified in place regardless
#'   of whether the return value is captured. The return value is provided for
#'   use in pipelines and for explicit documentation of intent.
#'
#' @examples
#'   gs <- apply_gating(gs)
#'   gs <- apply_gating(gs, template_path = "GatingTemplate.csv", mc.cores = 4)
# Helper: walk up GS hierarchy to nearest non-empty ancestor flowFrame
walk_up_fr <- function(gh, pop) {
  fr <- gh_pop_get_data(gh, pop)
  if (nrow(exprs(fr)) > 0) return(fr)
  cur <- pop
  while (TRUE) {
    parent <- dirname(cur)
    if (parent == "." || parent == cur) break
    fr2 <- tryCatch(gh_pop_get_data(gh, parent), error=function(e) NULL)
    if (!is.null(fr2) && nrow(exprs(fr2)) > 0) return(fr2)
    cur <- parent
  }
  fr
}


# ── Section 5: openCyto plugins ─────────────────────────────────────────────
# =============================================================================
# recovered_plugins.R
# Recovered from session transcript (local_8c184fa4 / mcp-workspace-bash log)
# All functions as final-written to ExtFunctions.R
# =============================================================================

# -----------------------------------------------------------------------------
# smooth_polygon
# (used by .gate_debris; lives in SECTION 1 or earlier in ExtFunctions.R)
# -----------------------------------------------------------------------------

smooth_polygon <- function(pts, n_out = 300) {
  # Close polygon if not already closed
  if (!all(pts[1, ] == pts[nrow(pts), ])) pts <- rbind(pts, pts[1, , drop = FALSE])
  # Arc-length parameterization so spline spacing is uniform along the boundary
  dists <- cumsum(c(0, sqrt(rowSums(diff(pts)^2))))
  dists <- dists / max(dists)   # normalize to [0, 1]
  t_out <- seq(0, 1, length.out = n_out + 1)[seq_len(n_out)]
  x_s <- spline(dists, pts[, 1], xout = t_out, method = "periodic")$y
  y_s <- spline(dists, pts[, 2], xout = t_out, method = "periodic")$y
  # Return closed polygon
  rbind(cbind(x_s, y_s), c(x_s[1], y_s[1]))
}


# -----------------------------------------------------------------------------
# Module-level cache for external control gate thresholds
# Populated by apply_gating() before gt_gating() runs.
# -----------------------------------------------------------------------------

.fasta_gate_cache <- new.env(parent = emptyenv())


# -----------------------------------------------------------------------------
# SECTION 5: CUSTOM openCyto PLUGINS
# -----------------------------------------------------------------------------

# .gate_static ----------------------------------------------------------------

#' Static rectangle gate plugin for openCyto
#'
#' Defines a fixed rectangle gate from user-specified boundaries.
#' Use in the gating template by setting gating_method = "gate_static"
#' and providing boundaries in gating_args as a quoted key=value string.
#'
#' Supported gating_args:
#'   1D (one channel in dims):
#'     min=<value>   lower bound  (default: -Inf)
#'     max=<value>   upper bound  (default:  Inf)
#'
#'   2D (two channels in dims, comma-separated):
#'     xmin=<value>  lower bound on first channel   (default: -Inf)
#'     xmax=<value>  upper bound on first channel   (default:  Inf)
#'     ymin=<value>  lower bound on second channel  (default: -Inf)
#'     ymax=<value>  upper bound on second channel  (default:  Inf)
#'
#' @examples
#' # In GatingTemplate.csv:
#' #   alias  pop  parent  dims         gating_method  gating_args
#' #   Live    -   Single  Viability    gate_static    "min=0,max=3500"
#' #   CD3+    +   Live    FSC-A,SSC-A  gate_static    "xmin=5000,xmax=Inf,ymin=2000,ymax=Inf"
.gate_static <- function(fr, pp_res, channels, filterId = "static", ...) {
  args <- list(...)

  # -- Validate channel count -----------------------------------------------
  if (!length(channels) %in% c(1, 2)) {
    stop("[gate_static] Only 1 or 2 channels are supported. Got: ",
         length(channels), " (", paste(channels, collapse=", "), ")")
  }

  # -- Check all specified args are recognised -------------------------------
  if (length(channels) == 1) {
    valid_args <- c("min", "max")
  } else {
    valid_args <- c("xmin", "xmax", "ymin", "ymax")
  }
  unknown_args <- setdiff(names(args), valid_args)
  if (length(unknown_args) > 0) {
    warning("[gate_static] Unrecognised gating_args (ignored): ",
            paste(unknown_args, collapse=", "),
            "\n  Valid args for ", length(channels), "D gate: ",
            paste(valid_args, collapse=", "))
  }

  # -- Require at least one recognised boundary arg --------------------------
  # A static gate with no recognised bounds is almost always a malformed
  # gating_args string (e.g. a token missing '=', which upstream parsing
  # silently drops), yielding a fully-open [-Inf, Inf] gate that looks like a
  # missing gate overlay downstream. Error loudly instead of gating nothing.
  if (!any(valid_args %in% names(args))) {
    stop("[gate_static] No recognised boundary args parsed from gating_args.\n",
         "  Expected ", paste(valid_args, collapse = "/"),
         " but got: ", if (length(names(args))) paste(names(args), collapse = ", ") else "(none)",
         ".\n  Check the gating_args string for a token missing '=' ",
         "(e.g. 'ymin21000' should be 'ymin=21000').")
  }

  # -- Build gate -----------------------------------------------------------
  if (length(channels) == 1) {
    gate_min <- if (!is.null(args$min)) args$min else -Inf
    gate_max <- if (!is.null(args$max)) args$max else  Inf
    gate_list <- setNames(list(c(gate_min, gate_max)), channels)

  } else {
    xmin <- if (!is.null(args$xmin)) args$xmin else -Inf
    xmax <- if (!is.null(args$xmax)) args$xmax else  Inf
    ymin <- if (!is.null(args$ymin)) args$ymin else -Inf
    ymax <- if (!is.null(args$ymax)) args$ymax else  Inf
    gate_list <- setNames(list(c(xmin, xmax), c(ymin, ymax)), channels)
  }

  rectangleGate(.gate = gate_list, filterId = filterId)
}

# -- Plugin registration


# .gate_debris ----------------------------------------------------------------

#' Smooth concave-hull debris gate plugin for openCyto
#'
#' Per-sample debris detection in FSC-A/SSC-A space using:
#'   1. Log transform (handles negatives from compensation via pmax(...,1))
#'   2. flowMeans clustering (K=2) to separate debris from cells
#'   3. Cluster scoring by normalised FSC+SSC median (lowest = debris)
#'   4. Concave hull (concaveman) or convex hull around debris events
#'   5. Periodic spline smoothing for a clean gate boundary
#'   6. Origin-closing edges to fully enclose the debris region
#'
#' Use pop = "-" in the gating template to keep non-debris (cell) events.
#'
#' gating_args options:
#'   method=per_sample|external_control  (default per_sample)
#'   hull_quantile=<value>  quantile for trimming outliers before hull (default 0.995)
#'   concavity=<value>      concaveman concavity parameter (default 2; higher = tighter)
#'   n_edge=<value>         points along each closing edge (default 30)
#'   min_events=<value>     minimum events required per cluster (default 50)
#'
#' @examples
#' # In GatingTemplate.csv:
#' #   alias   pop  parent  dims         gating_method  gating_args
#' #   Debris  -    root    FSC-A,SSC-A  gate_debris
#' #   Debris  -    root    FSC-A,SSC-A  gate_debris    "hull_quantile=0.99,concavity=3"
.gate_debris <- function(fr, pp_res, channels, filterId = "debris",
                          method = "per_sample",
                          hull_quantile = 0.995, concavity = 2,
                          n_edge = 30, min_events = 50,
                          control_col = NULL, control_val = NULL, ...) {

  # External control: return cached gate
  if (method == "external_control") {
    cache_key <- paste(channels, collapse = ",")
    if (!exists(cache_key, envir = .fasta_gate_cache, inherits = FALSE))
      stop("[gate_debris] [external_control] No cached gate for: ", paste(channels, collapse=","))
    cached <- get(cache_key, envir = .fasta_gate_cache)
    cached@filterId <- filterId
    message("  [gate_debris] [external_control] ", identifier(fr), " | gate from cache.")
    return(cached)
  }

  sn <- identifier(fr)

  # -- 1. Validate ----------------------------------------------------------
  if (length(channels) != 2) {
    stop("[gate_debris] Exactly 2 channels required. Got: ",
         paste(channels, collapse = ", "))
  }

  # -- 2. Extract data and log-transform ------------------------------------
  x     <- exprs(fr)[, channels, drop = FALSE]
  x_log <- log(pmax(x, 1))

  # -- 3. Cluster with flowMeans (K=2, Mahalanobis) -------------------------
  fm <- tryCatch(
    flowMeans::flowMeans(x = x_log, varNames = channels,
                         MaxN = 2, NumC = 2, Mahalanobis = TRUE),
    error = function(e) stop(
      "[gate_debris] flowMeans failed for: ", sn,
      "\n  Error: ", conditionMessage(e)
    )
  )

  df_log <- as.data.frame(cbind(x_log, cluster = fm@Labels[[1]]))

  # -- 4. Score clusters: lowest normalised FSC+SSC median = debris ---------
  range_norm <- function(v) {
    rng <- range(v, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0.5, length(v)))
    (v - rng[1]) / diff(rng)
  }

  cluster_stats <- df_log %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(
      n       = dplyr::n(),
      fsc_med = median(.data[[channels[1]]]),
      ssc_med = median(.data[[channels[2]]]),
      .groups = "drop"
    ) %>%
    dplyr::filter(n >= min_events) %>%
    dplyr::mutate(score = range_norm(fsc_med) + range_norm(ssc_med))

  if (nrow(cluster_stats) == 0) {
    warning("[gate_debris] No cluster met min_events=", min_events,
            " for: ", sn, " -- returning NULL.")
    return(NULL)
  }

  debris_id   <- cluster_stats$cluster[which.min(cluster_stats$score)]
  debris_mask <- fm@Labels[[1]] == debris_id
  coords      <- x[debris_mask, channels, drop = FALSE]

  message("  [gate_debris] ", sn, " | debris cluster: ", sum(debris_mask),
          " events | id=", debris_id)

  # -- 5. Trim outliers before hull -----------------------------------------
  fsc_q <- quantile(coords[, channels[1]], c(1 - hull_quantile, hull_quantile))
  ssc_q <- quantile(coords[, channels[2]], c(1 - hull_quantile, hull_quantile))

  trim_mask   <- coords[, channels[1]] >= fsc_q[1] & coords[, channels[1]] <= fsc_q[2] &
                 coords[, channels[2]] >= ssc_q[1] & coords[, channels[2]] <= ssc_q[2]
  coords_trim <- coords[trim_mask, , drop = FALSE]

  if (nrow(coords_trim) < min_events) {
    warning("[gate_debris] Too few events after trimming for: ", sn, " -- returning NULL.")
    return(NULL)
  }

  # -- 6. Build hull --------------------------------------------------------
  if (requireNamespace("concaveman", quietly = TRUE)) {
    hull_pts <- concaveman::concaveman(as.matrix(coords_trim), concavity = concavity)
    colnames(hull_pts) <- channels
  } else {
    message("  [gate_debris] concaveman not available -- using convex hull.")
    idx      <- chull(coords_trim)
    hull_pts <- coords_trim[c(idx, idx[1]), , drop = FALSE]
  }

  hull_smooth           <- smooth_polygon(hull_pts, n_out = 300)
  colnames(hull_smooth) <- channels

  # -- 7. Extract upper arc (debris-cell boundary) --------------------------
  pts       <- hull_smooth[seq_len(nrow(hull_smooth) - 1), ]
  n_pts     <- nrow(pts)
  left_idx  <- which.min(pts[, channels[1]])
  right_idx <- which.max(pts[, channels[1]])

  if (left_idx < right_idx) {
    arc_a <- pts[left_idx:right_idx,              , drop = FALSE]
    arc_b <- pts[c(right_idx:n_pts, 1:left_idx),  , drop = FALSE]
  } else {
    arc_a <- pts[c(left_idx:n_pts, 1:right_idx),  , drop = FALSE]
    arc_b <- pts[right_idx:left_idx,               , drop = FALSE]
  }

  # Upper arc = larger mean SSC (debris-cell boundary)
  upper_arc <- if (mean(arc_a[, channels[2]]) >= mean(arc_b[, channels[2]])) arc_a else arc_b

  # -- 8. Close polygon with edges back to origin ----------------------------
  # Add n_edge points along each closing edge so smooth_polygon constrains the spline
  right_pt  <- upper_arc[nrow(upper_arc), , drop = FALSE]
  left_pt   <- upper_arc[1,               , drop = FALSE]
  origin_pt <- matrix(c(-Inf, -Inf), nrow = 1, dimnames = list(NULL, channels))

  right_to_origin <- matrix(
    c(seq(right_pt[, channels[1]], -Inf, length.out = n_edge),
      seq(right_pt[, channels[2]], -Inf, length.out = n_edge)),
    ncol = 2, dimnames = list(NULL, channels)
  )
  origin_to_left <- matrix(
    c(seq(-Inf, left_pt[, channels[1]], length.out = n_edge),
      seq(-Inf, left_pt[, channels[2]], length.out = n_edge)),
    ncol = 2, dimnames = list(NULL, channels)
  )

  closed_gate <- rbind(upper_arc, right_to_origin, origin_to_left, left_pt)
  polygonGate(.gate = closed_gate, filterId = filterId)
}

# Register plugin


# .gate_singlet ---------------------------------------------------------------

#' Singlet gate plugin with external_control support
#'
#' Wraps flowStats::singletGate with the same method='per_sample'|'external_control'
#' pattern used by gate_flowmeans and gate_debris.
#'
#' gating_args options:
#'   method=per_sample|external_control  (default: per_sample)
#'   wider_gate=TRUE|FALSE               passed to singletGate (default FALSE)
#'   control_col, control_val, control_idx  (for external_control only)
#'
#' @examples
#' # In GatingTemplate.csv:
#' #   Single | + | Cells | FSC-A,FSC-H | gate_singlet |
#' #   Single | + | Cells | FSC-A,FSC-H | gate_singlet | "method='external_control',control_col='Group',control_val='PBMC'"
.gate_singlet <- function(fr, pp_res, channels, filterId = "Single",
                           method = "per_sample", wider_gate = FALSE,
                           control_col = NULL, control_val = NULL,
                           groupBy = NULL,
                           upper_mult = 1, preset_min = NA, preset_max = NA, ...) {
  # Note: upper_mult/preset_min/preset_max not applied to polygon gate.

  if (length(channels) != 2)
    stop("[gate_singlet] Exactly 2 channels required (area, height). Got: ",
         paste(channels, collapse=", "))

  # External control: return cached gate
  if (method == "external_control") {
    sn         <- identifier(fr)
    sn_key    <- paste0(paste(channels, collapse=","), "|", sn)
    cache_key <- paste(channels, collapse=",")
    active_key <- if (exists(sn_key, envir=.fasta_gate_cache, inherits=FALSE)) sn_key else cache_key
    if (!exists(active_key, envir=.fasta_gate_cache, inherits=FALSE))
      stop("[gate_singlet] [external_control] No cached gate for: ", active_key)
    cached <- get(active_key, envir = .fasta_gate_cache)
    cached@filterId <- filterId
    message("  [gate_singlet] [external_control] ", identifier(fr),
            " | gate from cache.")
    return(cached)
  }

  # Per-sample: compute singlet gate from this sample's data
  tryCatch(
    flowStats::singletGate(fr, area = channels[1], height = channels[2],
                           wider_gate = wider_gate, filterId = filterId),
    error = function(e) stop(
      "[gate_singlet] singletGate failed for: ", identifier(fr),
      "\n  Error: ", conditionMessage(e)
    )
  )
}



# .gate_flowmeans -------------------------------------------------------------

#' 1D flowMeans cluster-based gate plugin for openCyto
#'
#' Finds K=2 clusters on a single channel using flowMeans, selects the target
#' cluster by its maximum value (following the LYL314 convention), and places a
#' gate threshold at the specified quantile of that cluster.
#'
#' Use pop="+" to keep events ABOVE the threshold (positive cells).
#' Use pop="-" to keep events BELOW the threshold (negative cells / live).
#'
#' Cluster selection logic (by max value, not median):
#'   cluster="pos" -- uses the HIGH-max cluster; threshold at its low quantile
#'                    (default quantile=0.004), placing the cut at the start of
#'                    the high-value (dead/positive) cluster.
#'   cluster="neg" -- uses the LOW-max cluster; threshold at its high quantile
#'                    (default quantile=0.999).
#'
#' gating_args options:
#'   cluster=pos|neg         which cluster drives the threshold (required)
#'   quantile=<value>        quantile within the selected cluster (see defaults)
#'   upper_mult=<value>      multiply threshold by this factor (default 1)
#'   clusterFrom=<pop>       compute threshold from this parent population instead
#'   preset_min=<value>      minimum floor for the threshold (default NA)
#'   preset_max=<value>      maximum ceiling for the threshold (default NA)
#'   method=per_sample|external_control
#'   control_col, control_val, control_idx  (for external_control only)
#'
#' @examples
#' # In GatingTemplate.csv:
#' #   alias    pop  parent  dims       gating_method    gating_args
#' #   Live      -   Single  Viability  gate_flowmeans   "cluster='pos',quantile=0.004"
#' #   T cells   +   Live    CD3        gate_flowmeans   "cluster='pos',quantile=0.004"
.gate_flowmeans <- function(fr, pp_res, channels, filterId = "gate",
                             method = "per_sample", cluster = "pos",
                             quantile = NULL, upper_mult = 1,
                             clusterFrom = NULL,
                             preset_min = NA, preset_max = NA,
                             control_col = NULL, control_val = NULL, ...) {

  sn <- identifier(fr)

  # -- 1. Validate ----------------------------------------------------------
  if (length(channels) != 1) {
    stop("[gate_flowmeans] Exactly 1 channel required. Got: ",
         paste(channels, collapse = ", "))
  }
  if (!cluster %in% c("pos", "neg")) {
    stop("[gate_flowmeans] cluster must be 'pos' or 'neg'. Got: ", cluster)
  }
  if (!method %in% c("per_sample", "external_control")) {
    stop("[gate_flowmeans] method must be 'per_sample' or 'external_control'. Got: ", method)
  }

  # -- External control path: read cached threshold from apply_gating ---------
  if (method == "external_control") {
    sn_key    <- paste(channels, sn, sep="|")  # per-sample key (clusterFrom)
    cache_key <- paste(channels, collapse = ",")             # single-threshold key (external_control)
    # Prefer sample-specific key (set when clusterFrom is used in per_sample mode)
    active_key <- if (exists(sn_key, envir=.fasta_gate_cache, inherits=FALSE)) sn_key else cache_key
    if (!exists(active_key, envir=.fasta_gate_cache, inherits=FALSE)) {
      stop("[gate_flowmeans] [external_control] No cached threshold for: ", filterId,
           "\n  Ensure apply_gating() pre-computed the threshold before gt_gating().")
    }
    threshold <- get(active_key, envir=.fasta_gate_cache) * upper_mult
    if (!is.na(preset_min)) threshold <- max(threshold, preset_min)
    if (!is.na(preset_max)) threshold <- min(threshold, preset_max)
    message("  [gate_flowmeans] [external_control] ", identifier(fr),
            " | threshold: ", round(threshold, 2), " (from cache, upper_mult=", upper_mult, ")")
    gate_list <- setNames(list(c(threshold, Inf)), channels)
    return(rectangleGate(.gate=gate_list, filterId=filterId))
  }

  # Check for per-sample cached threshold (set when clusterFrom is used)
  sn_cache_key <- paste(channels, identifier(fr), sep="|")
  if (exists(sn_cache_key, envir=.fasta_gate_cache, inherits=FALSE)) {
    threshold <- get(sn_cache_key, envir=.fasta_gate_cache) * upper_mult
    if (!is.na(preset_min)) threshold <- max(threshold, preset_min)
    if (!is.na(preset_max)) threshold <- min(threshold, preset_max)
    message("  [gate_flowmeans] ", sn, " | channel: ", channels,
            " | threshold: ", round(threshold, 2), " (from clusterFrom cache)")
    gate_list <- setNames(list(c(threshold, Inf)), channels)
    return(rectangleGate(.gate=gate_list, filterId=filterId))
  }

  # -- 2. Extract data ------------------------------------------------------
  x <- exprs(fr)[, channels, drop = FALSE]
  if (nrow(x) < 20) {
    stop("[gate_flowmeans] Too few events (", nrow(x), ") to cluster: ", sn)
  }

  # -- 3. Cluster with flowMeans K=2 ----------------------------------------
  fm <- tryCatch(
    flowMeans::flowMeans(x = x, varNames = channels, MaxN = 2, NumC = 2),
    error = function(e) stop(
      "[gate_flowmeans] flowMeans failed for: ", sn,
      "\n  Channel: ", channels,
      "\n  Error: ", conditionMessage(e)
    )
  )

  df <- as.data.frame(cbind(x, cluster_id = fm@Labels[[1]]))
  colnames(df) <- c("val", "cluster_id")

  # -- 4. Identify clusters by max value (LYL314 convention) ----------------
  clust_summary <- df %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarize(mx = max(val), .groups = "drop") %>%
    dplyr::arrange(mx)

  if (nrow(clust_summary) < 2) {
    warning("[gate_flowmeans] Only 1 cluster found for: ", sn,
            " on channel: ", channels)
  }

  pos_id <- clust_summary$cluster_id[which.max(clust_summary$mx)]  # highest max
  neg_id <- clust_summary$cluster_id[which.min(clust_summary$mx)]  # lowest max
  selected_id <- if (cluster == "pos") pos_id else neg_id

  cluster_vals <- df$val[df$cluster_id == selected_id]

  # -- 5. Default quantile --------------------------------------------------
  if (is.null(quantile)) {
    quantile <- if (cluster == "pos") 0.004 else 0.999
  }

  # -- 6. Compute threshold -------------------------------------------------
  threshold <- stats::quantile(cluster_vals, quantile, names = FALSE) * upper_mult
  if (!is.na(preset_min)) threshold <- max(threshold, preset_min)
  if (!is.na(preset_max)) threshold <- min(threshold, preset_max)

  message("  [gate_flowmeans] ", sn,
          " | channel: ", channels,
          " | cluster: ", cluster,
          " | threshold: ", round(threshold, 2))

  # -- 7. Return rectangleGate [threshold, Inf] -----------------------------
  gate_list <- setNames(list(c(threshold, Inf)), channels)
  rectangleGate(.gate = gate_list, filterId = filterId)
}



# .gate_flowmeans_2d ----------------------------------------------------------
# NOTE: registered as "gate_flowmeans_2d" (with underscore before 2d)

#' 2D flowMeans cluster-based gate plugin for openCyto
#'
#' Runs independent 1D flowMeans (K=2) on each of two channels, placing a
#' threshold on each axis. Returns a quadGate with thresholds at
#' (threshold_x, threshold_y), creating 4 populations.
#' All cluster/quantile/preset args are suffixed _x or _y to specify per-axis.
#'
#' gating_args options:
#'   cluster_x=pos|neg       cluster for x-axis (default pos)
#'   cluster_y=pos|neg       cluster for y-axis (default pos)
#'   quantile_x=<value>      quantile within x cluster (default: 0.004 for pos, 0.999 for neg)
#'   quantile_y=<value>      quantile within y cluster (default: 0.004 for pos, 0.999 for neg)
#'   upper_mult_x=<value>    multiply x threshold (default 1)
#'   upper_mult_y=<value>    multiply y threshold (default 1)
#'   clusterFrom=<pop>       compute both thresholds from this population
#'   preset_min_x=<value>    floor for x threshold (default NA)
#'   preset_max_x=<value>    ceiling for x threshold (default NA)
#'   preset_min_y=<value>    floor for y threshold (default NA)
#'   preset_max_y=<value>    ceiling for y threshold (default NA)
#'   name_x, name_y          labels for quad populations
#'   method=per_sample|external_control
#'   control_col, control_val, control_idx  (for external_control only)
#'
#' @examples
#' # In GatingTemplate.csv:
#' #   CD4+CD8+  +  T cells  CD4,CD8  gate_flowmeans_2d  "cluster_x='pos',cluster_y='pos'"
.gate_flowmeans_2d <- function(fr, pp_res, channels, filterId = "gate",
                               method      = "per_sample",
                               cluster_x   = "pos",  cluster_y   = "pos",
                               quantile_x  = NULL,   quantile_y  = NULL,
                               upper_mult_x = 1,    upper_mult_y = 1,
                               clusterFrom = NULL,
                               preset_min_x = NA,   preset_max_x = NA,
                               preset_min_y = NA,   preset_max_y = NA,
                               name_x = NULL, name_y = NULL,
                               control_col = NULL,  control_val = NULL, ...) {

  sn <- identifier(fr)

  # -- 1. Validate --------------------------------------------------------------
  if (length(channels) != 2)
    stop("[gate_flowmeans_2d] Exactly 2 channels required. Got: ",
         paste(channels, collapse=", "))
  for (cl in c(cluster_x, cluster_y))
    if (!cl %in% c("pos","neg"))
      stop("[gate_flowmeans_2d] cluster_x and cluster_y must be 'pos' or 'neg'. Got: ", cl)

  # -- 2. External control: return cached gate ----------------------------------
  if (method == "external_control") {
    cache_key <- paste(channels, collapse=",")
    if (!exists(cache_key, envir=.fasta_gate_cache, inherits=FALSE)) {
      stop("[gate_flowmeans_2d] [external_control] No cached gate for '", filterId,
           "' (channel: ", paste(channels, collapse=','), ").",
           "\n  Check that precompute_external_control ran for the control sample.",
           "\n  Or use clusterFrom=\"<population>\" to compute from a different parent.")
    } else {
      cached <- get(cache_key, envir=.fasta_gate_cache)
      cached@filterId <- filterId
      message("  [gate_flowmeans_2d] [external_control] ", sn, " | gate from cache.")
      return(cached)
    }
  }

  # Check clusterFrom per-sample cache before computing from fr
  sn_key_x <- paste(channels[1], sn, sep="|")
  sn_key_y <- paste(channels[2], sn, sep="|")
  if (exists(sn_key_x, envir=.fasta_gate_cache, inherits=FALSE) ||
      exists(sn_key_y, envir=.fasta_gate_cache, inherits=FALSE)) {
    thr_x2 <- if (exists(sn_key_x, envir=.fasta_gate_cache, inherits=FALSE))
      get(sn_key_x, envir=.fasta_gate_cache) * upper_mult_x else
      stats::quantile(exprs(fr)[, channels[1]], quantile_x, names=FALSE) * upper_mult_x
    thr_y2 <- if (exists(sn_key_y, envir=.fasta_gate_cache, inherits=FALSE))
      get(sn_key_y, envir=.fasta_gate_cache) * upper_mult_y else
      stats::quantile(exprs(fr)[, channels[2]], quantile_y, names=FALSE) * upper_mult_y
    if (!is.na(preset_min_x)) thr_x2 <- max(thr_x2, preset_min_x)
    if (!is.na(preset_max_x)) thr_x2 <- min(thr_x2, preset_max_x)
    if (!is.na(preset_min_y)) thr_y2 <- max(thr_y2, preset_min_y)
    if (!is.na(preset_max_y)) thr_y2 <- min(thr_y2, preset_max_y)
    message("  [gate_flowmeans_2d] [per_sample] ", sn, " (clusterFrom cache)",
            " | ", channels[1], "=", round(thr_x2,1), " | ", channels[2], "=", round(thr_y2,1))
    qx2 <- if (!is.null(name_x)) name_x else channels[1]
    qy2 <- if (!is.null(name_y)) name_y else channels[2]
    qn2 <- c(paste0(qx2,"-",qy2,"+"), paste0(qx2,"+",qy2,"+"),
             paste0(qx2,"+",qy2,"-"), paste0(qx2,"-",qy2,"-"))
    assign(paste0(filterId,":quad_names"), qn2, envir=.fasta_gate_cache)
    return(quadGate(.gate=setNames(list(thr_x2, thr_y2), channels), filterId=filterId))
  }

  # -- 3. Helper: compute 1D threshold on one channel --------------------------
  threshold_1d <- function(channel, cluster, quantile_val) {
    x  <- exprs(fr)[, channel, drop=FALSE]
    fm <- tryCatch(
      flowMeans::flowMeans(x, varNames=channel, MaxN=2, NumC=2),
      error=function(e) stop(
        "[gate_flowmeans_2d] flowMeans failed on channel '", channel,
        "' for: ", sn, "\n  Error: ", conditionMessage(e))
    )
    df <- as.data.frame(cbind(x, cluster_id=fm@Labels[[1]]))
    colnames(df) <- c("val","cluster_id")
    cs_sum <- df %>% dplyr::group_by(cluster_id) %>%
      dplyr::summarize(med=median(val), .groups="drop") %>% dplyr::arrange(med)
    sel_id    <- if (cluster=="pos") cs_sum$cluster_id[which.max(cs_sum$med)] else
                                     cs_sum$cluster_id[which.min(cs_sum$med)]
    threshold <- stats::quantile(df$val[df$cluster_id==sel_id], quantile_val, names=FALSE)
    # preset_min/max applied AFTER upper_mult by caller
    threshold
  }

  if (is.null(quantile_x)) quantile_x <- if (cluster_x=="pos") 0.004 else 0.999
  if (is.null(quantile_y)) quantile_y <- if (cluster_y=="pos") 0.004 else 0.999

  # -- 4. Compute thresholds for each channel -----------------------------------
  thr_x <- threshold_1d(channels[1], cluster_x, quantile_x) * upper_mult_x
  thr_y <- threshold_1d(channels[2], cluster_y, quantile_y) * upper_mult_y
  if (!is.na(preset_min_x)) thr_x <- max(thr_x, preset_min_x)
  if (!is.na(preset_max_x)) thr_x <- min(thr_x, preset_max_x)
  if (!is.na(preset_min_y)) thr_y <- max(thr_y, preset_min_y)
  if (!is.na(preset_max_y)) thr_y <- min(thr_y, preset_max_y)

  message("  [gate_flowmeans_2d] ", sn,
          " | ", channels[1], "=", round(thr_x, 1),
          " | ", channels[2], "=", round(thr_y, 1))

  # -- 5. Build quad names and return quadGate ---------------------------------
  qx <- if (!is.null(name_x)) name_x else channels[1]
  qy <- if (!is.null(name_y)) name_y else channels[2]
  quad_names <- c(
    paste0(qx, "-", qy, "+"),
    paste0(qx, "+", qy, "+"),
    paste0(qx, "+", qy, "-"),
    paste0(qx, "-", qy, "-")
  )
  # Cache the 4 quad names so apply_gating can use them in gs_pop_add
  assign(paste0(filterId, ":quad_names"), quad_names, envir = .fasta_gate_cache)
  quadGate(.gate=setNames(list(thr_x, thr_y), channels), filterId=filterId)
}



# -----------------------------------------------------------------------------
# SECTION 7: FMO GATING PLUGINS
# -----------------------------------------------------------------------------

# .gate_FMO -------------------------------------------------------------------

#' 1D FMO (Fluorescence Minus One) gate plugin for openCyto
#'
#' Places a threshold at the specified quantile of the FMO control sample
#' distribution. No clustering -- the FMO sample IS the negative population.
#'
#' Use pop="+" to keep events above the threshold (positive cells).
#' Use pop="-" to keep events below the threshold (negative cells).
#'
#' gating_args options:
#'   quantile=<value>        percentile of FMO distribution as threshold (default 0.999)
#'   upper_mult=<value>      multiply threshold by this factor (default 1)
#'   preset_min=<value>      floor for threshold (default NA)
#'   preset_max=<value>      ceiling for threshold (default NA)
#'   method=per_sample|external_control  (default external_control)
#'   control_col, control_val, control_idx  identify the FMO sample
#'   clusterFrom=<pop>       use data from this population (default: parent)
#'
#' @examples
#' #   CD19+  +  DUMP-  CD19  gate_FMO  "control_col='Group',control_val='FMO_CD19'"
.gate_FMO <- function(fr, pp_res, channels, filterId = "gate",
                       method = "external_control",
                       quantile = 0.999, upper_mult = 1,
                       clusterFrom = NULL,
                       preset_min = NA, preset_max = NA,
                       control_col = NULL, control_val = NULL,
                       control_idx = 1L, ...) {

  sn <- identifier(fr)

  if (length(channels) != 1)
    stop("[gate_FMO] Exactly 1 channel required. Got: ", paste(channels, collapse=", "))
  if (!method %in% c("per_sample", "external_control"))
    stop("[gate_FMO] method must be 'per_sample' or 'external_control'. Got: ", method)

  # -- External control: read cached threshold ----------------------------------
  if (method == "external_control") {
    sn_key     <- paste(channels, sn, sep="|")
    cache_key  <- paste(channels, collapse=",")
    active_key <- if (exists(sn_key, envir=.fasta_gate_cache, inherits=FALSE)) sn_key else cache_key
    if (!exists(active_key, envir=.fasta_gate_cache, inherits=FALSE))
      stop("[gate_FMO] [external_control] No cached FMO threshold for channel: ",
           channels, "\n  Ensure apply_gating() pre-computed the FMO threshold.")
    threshold <- get(active_key, envir=.fasta_gate_cache) * upper_mult
    if (!is.na(preset_min)) threshold <- max(threshold, preset_min)
    if (!is.na(preset_max)) threshold <- min(threshold, preset_max)
    message("  [gate_FMO] [external_control] ", sn,
            " | threshold: ", round(threshold, 2),
            if (active_key == sn_key) " (per-sample)" else " (global)")
    gate_list <- setNames(list(c(threshold, Inf)), channels)
    return(rectangleGate(.gate=gate_list, filterId=filterId))
  }

  # -- Per-sample: check clusterFrom cache first, then compute from fr ---------
  sn_cache_key <- paste(channels, sn, sep="|")
  if (exists(sn_cache_key, envir=.fasta_gate_cache, inherits=FALSE)) {
    threshold <- get(sn_cache_key, envir=.fasta_gate_cache) * upper_mult
    if (!is.na(preset_min)) threshold <- max(threshold, preset_min)
    if (!is.na(preset_max)) threshold <- min(threshold, preset_max)
    message("  [gate_FMO] [per_sample] ", sn,
            " | threshold: ", round(threshold, 2), " (from clusterFrom cache)")
    gate_list <- setNames(list(c(threshold, Inf)), channels)
    return(rectangleGate(.gate=gate_list, filterId=filterId))
  }

  x         <- exprs(fr)[, channels, drop=FALSE]
  threshold <- stats::quantile(x[, channels], quantile, names=FALSE) * upper_mult
  if (!is.na(preset_min)) threshold <- max(threshold, preset_min)
  if (!is.na(preset_max)) threshold <- min(threshold, preset_max)
  message("  [gate_FMO] [per_sample] ", sn, " | threshold: ", round(threshold, 2))
  gate_list <- setNames(list(c(threshold, Inf)), channels)
  rectangleGate(.gate=gate_list, filterId=filterId)
}



# .gate_FMO_2d ----------------------------------------------------------------

#' 2D FMO gate plugin for openCyto
#'
#' Places thresholds at specified quantiles of the FMO control sample on two
#' channels independently. Returns a quadGate creating 4 populations.
#' Args follow the _x/_y suffix convention matching gate_flowmeans_2d.
#'
#' gating_args options:
#'   quantile_x, quantile_y      (default 0.999 each)
#'   upper_mult_x, upper_mult_y  (default 1 each)
#'   preset_min_x/y, preset_max_x/y
#'   name_x, name_y              labels for quad populations
#'   method, control_col, control_val, control_idx, clusterFrom
#'
#' @examples
#' #   CD3_G4S  +  DUMP-  CD3,G4S  gate_FMO_2d
#' #     "control_col='Group',control_val='FMO',name_x='CD3',name_y='G4S'"
.gate_FMO_2d <- function(fr, pp_res, channels, filterId = "gate",
                          method       = "external_control",
                          quantile_x   = 0.999, quantile_y   = 0.999,
                          upper_mult_x = 1,    upper_mult_y = 1,
                          clusterFrom  = NULL,
                          preset_min_x = NA,   preset_max_x = NA,
                          preset_min_y = NA,   preset_max_y = NA,
                          name_x = NULL, name_y = NULL,
                          control_col = NULL, control_val = NULL,
                          control_idx = 1L, ...) {

  sn <- identifier(fr)

  if (length(channels) != 2)
    stop("[gate_FMO_2d] Exactly 2 channels required. Got: ", paste(channels, collapse=", "))
  if (!method %in% c("per_sample", "external_control"))
    stop("[gate_FMO_2d] method must be 'per_sample' or 'external_control'. Got: ", method)

  # -- External control: return cached gate -------------------------------------
  if (method == "external_control") {
    sn_key    <- paste(paste(channels, collapse=","), sn, sep="|")
    cache_key <- paste(channels, collapse=",")
    active_key <- if (exists(sn_key, envir=.fasta_gate_cache, inherits=FALSE)) sn_key else cache_key
    if (!exists(active_key, envir=.fasta_gate_cache, inherits=FALSE))
      stop("[gate_FMO_2d] [external_control] No cached FMO gate for channels: ",
           paste(channels, collapse=","))
    cached <- get(active_key, envir=.fasta_gate_cache)
    cached@filterId <- filterId
    message("  [gate_FMO_2d] [external_control] ", sn, " | gate from cache",
            if (active_key == sn_key) " (per-sample)" else " (global)", ".")
    return(cached)
  }

  # -- Per-sample: check clusterFrom cache before computing from fr -------------
  sn_key_x <- paste(channels[1], sn, sep="|")
  sn_key_y <- paste(channels[2], sn, sep="|")
  if (exists(sn_key_x, envir=.fasta_gate_cache, inherits=FALSE) ||
      exists(sn_key_y, envir=.fasta_gate_cache, inherits=FALSE)) {
    thr_x2 <- if (exists(sn_key_x, envir=.fasta_gate_cache, inherits=FALSE))
      get(sn_key_x, envir=.fasta_gate_cache) * upper_mult_x else
      stats::quantile(exprs(fr)[, channels[1]], quantile_x, names=FALSE) * upper_mult_x
    thr_y2 <- if (exists(sn_key_y, envir=.fasta_gate_cache, inherits=FALSE))
      get(sn_key_y, envir=.fasta_gate_cache) * upper_mult_y else
      stats::quantile(exprs(fr)[, channels[2]], quantile_y, names=FALSE) * upper_mult_y
    if (!is.na(preset_min_x)) thr_x2 <- max(thr_x2, preset_min_x)
    if (!is.na(preset_max_x)) thr_x2 <- min(thr_x2, preset_max_x)
    if (!is.na(preset_min_y)) thr_y2 <- max(thr_y2, preset_min_y)
    if (!is.na(preset_max_y)) thr_y2 <- min(thr_y2, preset_max_y)
    message("  [gate_FMO_2d] [per_sample] ", sn, " (clusterFrom cache)",
            " | ", channels[1], "=", round(thr_x2,1), " | ", channels[2], "=", round(thr_y2,1))
    qx2 <- if (!is.null(name_x)) name_x else channels[1]
    qy2 <- if (!is.null(name_y)) name_y else channels[2]
    qn2 <- c(paste0(qx2,"-",qy2,"+"), paste0(qx2,"+",qy2,"+"),
             paste0(qx2,"+",qy2,"-"), paste0(qx2,"-",qy2,"-"))
    assign(paste0(filterId, ":quad_names"), qn2, envir=.fasta_gate_cache)
    return(quadGate(.gate=setNames(list(thr_x2, thr_y2), channels), filterId=filterId))
  }

  # -- Per-sample: compute quantile thresholds from fr --------------------------
  thr_x <- stats::quantile(exprs(fr)[, channels[1]], quantile_x, names=FALSE) * upper_mult_x
  thr_y <- stats::quantile(exprs(fr)[, channels[2]], quantile_y, names=FALSE) * upper_mult_y
  if (!is.na(preset_min_x)) thr_x <- max(thr_x, preset_min_x)
  if (!is.na(preset_max_x)) thr_x <- min(thr_x, preset_max_x)
  if (!is.na(preset_min_y)) thr_y <- max(thr_y, preset_min_y)
  if (!is.na(preset_max_y)) thr_y <- min(thr_y, preset_max_y)

  message("  [gate_FMO_2d] [per_sample] ", sn,
          " | ", channels[1], "=", round(thr_x,1),
          " | ", channels[2], "=", round(thr_y,1))

  qx <- if (!is.null(name_x)) name_x else channels[1]
  qy <- if (!is.null(name_y)) name_y else channels[2]
  quad_names <- c(paste0(qx,"-",qy,"+"), paste0(qx,"+",qy,"+"),
                  paste0(qx,"+",qy,"-"), paste0(qx,"-",qy,"-"))
  assign(paste0(filterId, ":quad_names"), quad_names, envir=.fasta_gate_cache)
  quadGate(.gate=setNames(list(thr_x, thr_y), channels), filterId=filterId)
}



# .gate_tmix_2d ---------------------------------------------------------------
# NOTE: registered as "gate_tmix_2d"

#' 2D quadrant gate plugin wrapping openCyto's gate_quad_tmix (quadGate.tmix)
#'
#' Fits a JOINT bivariate normal mixture (flowClust, K components) on the two
#' channels and returns the FOUR quadrant polygonGates that gate_quad_tmix
#' natively produces (a flowCore `filters` object) -- NOT a single quadGate.
#' Unlike gate_flowmeans_2d / gate_FMO_2d (independent per-axis thresholds), this
#' is a joint fit, so it takes a SINGLE control sample (no control_col_x/y).
#'
#' gating_args options:
#'   method        "per_sample" (default) or "external_control"
#'   K             number of mixture components (default 2)
#'   usePrior      passed to gate_quad_tmix (default "no")
#'   prior         passed to gate_quad_tmix (default list(NA))
#'   quantile1     lower-quadrant boundary quantile (default 0.8)
#'   quantile3     upper-quadrant boundary quantile (default 0.8)
#'   trans         passed to gate_quad_tmix (default 0)
#'   name_x, name_y   quadrant naming (default channel names)
#'   clusterFrom   fit on this population instead of parent (per-sample cache)
#'
#' Returns a flowCore `filters` object of 4 polygonGates named, in gate_quad_tmix
#' order, c("<x>-<y>+","<x>+<y>+","<x>+<y>-","<x>-<y>-") from name_x/name_y. The
#' inline 2D branch in apply_gating adds these 4 polygons per sample.
#'
#' Cache contract:
#'   external_control       -> reads cached `filters` at paste(channels, collapse=",")
#'   clusterFrom/per-sample -> reads cached `filters` at paste(channels, collapse=",")|sn
#'   sets paste0(filterId, ":quad_names") (4-element vector, same order)
.gate_tmix_2d <- function(fr, pp_res, channels, filterId = "gate",
                          method      = "per_sample",
                          K           = 2,
                          usePrior    = "no",
                          prior       = list(NA),
                          quantile1   = 0.8,
                          quantile3   = 0.8,
                          trans       = 0,
                          name_x = NULL, name_y = NULL,
                          clusterFrom = NULL,
                          control_col = NULL, control_val = NULL, ...) {

  sn <- identifier(fr)

  # -- 1. Validate --------------------------------------------------------------
  if (length(channels) != 2)
    stop("[gate_tmix_2d] Exactly 2 channels required. Got: ",
         paste(channels, collapse=", "))

  # Quadrant names (gate_quad_tmix returns the 4 polygons in this fixed order).
  qx <- if (!is.null(name_x)) name_x else channels[1]
  qy <- if (!is.null(name_y)) name_y else channels[2]
  quad_names <- c(paste0(qx, "-", qy, "+"), paste0(qx, "+", qy, "+"),
                  paste0(qx, "+", qy, "-"), paste0(qx, "-", qy, "-"))

  # -- 2. External control: return cached filters (4 polygons) ------------------
  if (method == "external_control") {
    cache_key <- paste(channels, collapse=",")
    if (!exists(cache_key, envir=.fasta_gate_cache, inherits=FALSE))
      stop("[gate_tmix_2d] [external_control] No cached gate for '", filterId,
           "' (channel: ", paste(channels, collapse=','), ").",
           "\n  Check that precompute_external_control ran for the control sample.",
           "\n  Or use clusterFrom=\"<population>\" to compute from a different parent.")
    message("  [gate_tmix_2d] [external_control] ", sn, " | gates from cache.")
    return(get(cache_key, envir=.fasta_gate_cache))
  }

  # -- 3. clusterFrom / per-sample joint cache (filters) ------------------------
  sn_key <- paste(paste(channels, collapse=","), sn, sep="|")
  if (exists(sn_key, envir=.fasta_gate_cache, inherits=FALSE)) {
    message("  [gate_tmix_2d] [per_sample] ", sn, " (clusterFrom/joint cache).")
    return(get(sn_key, envir=.fasta_gate_cache))
  }

  # -- 4. Compute (per_sample): joint tmix fit -> 4 polygon gates ---------------
  if (is.null(K))
    stop("[gate_tmix_2d] K (number of mixture components) is required.")
  # openCyto's quadGate.tmix is deprecated; gate_quad_tmix is the current exported
  # entry point (same signature). It returns a `filters` of 4 polygonGates.
  fobj <- tryCatch(
    openCyto::gate_quad_tmix(fr, channels = channels, K = K,
      usePrior = usePrior, prior = prior,
      quantile1 = quantile1, quantile3 = quantile3,
      trans = trans, plot = FALSE),
    error = function(e) stop(
      "[gate_tmix_2d] gate_quad_tmix failed for: ", sn,
      "\n  Channels: ", paste(channels, collapse=", "),
      "\n  Error: ", conditionMessage(e))
  )
  if (length(fobj) != 4)
    stop("[gate_tmix_2d] expected 4 quadrant gates from gate_quad_tmix, got ",
         length(fobj), " for: ", sn)

  # Rename the 4 polygonGates to the FASTA quad_names, preserving order.
  polys <- lapply(seq_len(4), function(j) {
    g <- fobj[[j]]; g@filterId <- quad_names[j]; g
  })
  out <- flowCore::filters(polys)
  assign(paste0(filterId, ":quad_names"), quad_names, envir = .fasta_gate_cache)
  message("  [gate_tmix_2d] ", sn, " | 4 quad polygons: ",
          paste(quad_names, collapse=", "), " (K=", K, ")")
  out
}



#' .gate_mindensity -- openCyto plugin wrapper with per-sample and external_control support
#'
#' Args (gating_args column):
#'   method          "per_sample" (default) or "external_control"
#'   min             lower bound for density search (default NULL)
#'   max             upper bound for density search (default NULL)
#'   gate_range      c(lo, hi) range to search for minimum density (default NULL)
#'   peaks           number of peaks expected (default NULL = auto)
#'   adjust          bandwidth adjustment for density estimation (default 2)
#'   upper_mult      multiply raw threshold by this factor before preset checks (default 1)
#'   preset_min      floor applied after upper_mult: max(thr*upper_mult, preset_min) (default NA)
#'   preset_max      ceiling applied after upper_mult: min(thr*upper_mult, preset_max) (default NA)
#'   control_col     pData column used to identify control sample (external_control only)
#'   control_val     pData value matching the control sample (external_control only)
#'   groupBy         optional pData column for grouped cache keys
.gate_mindensity <- function(fr, pp_res, channels, filterId = "",
                              method      = "per_sample",
                              min         = NULL,  max       = NULL,
                              gate_range  = NULL,  peaks     = NULL,
                              positive    = TRUE,
                              upper_mult  = 1,
                              preset_min  = NA,    preset_max = NA,
                              control_col = NULL,  control_val = NULL,
                              groupBy     = NULL, ...) {

  sn      <- identifier(fr)
  channel <- channels[1]

  # Helper: run openCyto gate_mindensity, extract threshold, apply multiplier + presets
  compute_threshold <- function(fr_ctrl) {
    g <- tryCatch(
      openCyto::gate_mindensity(fr_ctrl, channel = unname(channels[1]),
                                 filterId = filterId,
                                 min = min, max = max,
                                 gate_range = gate_range, peaks = peaks,
                                 positive = positive),
      error = function(e) stop("[gate_mindensity] Failed: ", conditionMessage(e))
    )
    # Extract lower bound from returned rectangleGate
    bounds <- c(g@min[[unname(channels[1])]], g@max[[unname(channels[1])]])
    thr <- bounds[is.finite(bounds)][1]
    thr <- thr * upper_mult
    if (!is.na(preset_min) && is.finite(preset_min)) thr <- base::max(thr, preset_min)
    if (!is.na(preset_max) && is.finite(preset_max)) thr <- base::min(thr, preset_max)
    thr
  }

  if (method == "external_control") {

    # Build cache key (optionally per-groupBy group)
    sn        <- identifier(fr)
    sn_key    <- paste0(channel, "|", sn)
    cache_key <- channel
    active_key <- if (exists(sn_key, envir=.fasta_gate_cache, inherits=FALSE)) sn_key else cache_key
    if (!exists(active_key, envir=.fasta_gate_cache, inherits=FALSE))
      stop("[gate_mindensity] External control not yet cached for key: ", active_key)
    # Use active_key (per-sample sn_key when present, else the global channel
    # key). The previous duplicate line re-read the global cache_key and
    # discarded the per-sample threshold -- removed.
    thr <- get(active_key, envir=.fasta_gate_cache) * upper_mult
    if (!is.na(preset_min) && is.finite(preset_min)) thr <- base::max(thr, preset_min)
    if (!is.na(preset_max) && is.finite(preset_max)) thr <- base::min(thr, preset_max)
    message("  [gate_mindensity] ext_ctrl | ", channel, " = ", round(thr, 1), " (sn: ", sn, ")")

  } else {
    # per_sample
    thr <- compute_threshold(fr)
    message("  [gate_mindensity] per_sample | ", channel, " = ", round(thr, 1), " (sn: ", sn, ")")
  }

  rectangleGate(.gate = setNames(list(c(thr, Inf)), channel), filterId = filterId)
}


#' .gate_mindensity_2d -- 2D quadGate using independent minimum-density thresholds per axis
.gate_mindensity_2d <- function(fr, pp_res, channels, filterId = "",
                                 method       = "per_sample",
                                 min_x=NULL,   min_y=NULL,
                                 max_x=NULL,   max_y=NULL,
                                 gate_range_x=NULL, gate_range_y=NULL,
                                 peaks_x=NULL, peaks_y=NULL,
                                 positive_x=TRUE, positive_y=TRUE,
                                 upper_mult_x=1,  upper_mult_y=1,
                                 preset_min_x=NA, preset_min_y=NA,
                                 preset_max_x=NA, preset_max_y=NA,
                                 control_col=NULL, control_val=NULL,
                                 name_x=NULL, name_y=NULL,
                                 groupBy=NULL, ...) {
  if (length(channels) != 2)
    stop("[gate_mindensity_2d] Exactly 2 channels required.")
  ch_x <- channels[1]; ch_y <- channels[2]
  sn   <- identifier(fr)
  # Helper: run gate_mindensity on one channel, return numeric threshold
  .get_thr2d <- function(fr_data, ch, mn, mx, gr, pk, pos, umult, pmin, pmax) {
    g <- tryCatch(
      openCyto::gate_mindensity(fr_data, channel=unname(ch), filterId="tmp",
                                 min=mn, max=mx, gate_range=gr, peaks=pk, positive=pos),
      error=function(e) stop("[gate_mindensity_2d] Failed: ", conditionMessage(e)))
    bounds <- c(g@min[[unname(ch)]], g@max[[unname(ch)]])
    thr <- bounds[is.finite(bounds)][1] * umult
    if (!is.na(pmin) && is.finite(pmin)) thr <- base::max(thr, pmin)
    if (!is.na(pmax) && is.finite(pmax)) thr <- base::min(thr, pmax)
    thr
  }
  if (method == "external_control") {
    sn    <- identifier(fr)
    sn_x  <- paste0(ch_x, "|", sn); sn_y  <- paste0(ch_y, "|", sn)
    key_x <- if (exists(sn_x, envir=.fasta_gate_cache, inherits=FALSE)) sn_x else ch_x
    key_y <- if (exists(sn_y, envir=.fasta_gate_cache, inherits=FALSE)) sn_y else ch_y
    if (!exists(key_x, envir=.fasta_gate_cache, inherits=FALSE) ||
        !exists(key_y, envir=.fasta_gate_cache, inherits=FALSE))
      stop("[gate_mindensity_2d] Thresholds not cached for: ", key_x, ", ", key_y)
    thr_x <- get(key_x, envir=.fasta_gate_cache) * upper_mult_x
    thr_y <- get(key_y, envir=.fasta_gate_cache) * upper_mult_y
    if (!is.na(preset_min_x) && is.finite(preset_min_x)) thr_x <- base::max(thr_x, preset_min_x)
    if (!is.na(preset_max_x) && is.finite(preset_max_x)) thr_x <- base::min(thr_x, preset_max_x)
    if (!is.na(preset_min_y) && is.finite(preset_min_y)) thr_y <- base::max(thr_y, preset_min_y)
    if (!is.na(preset_max_y) && is.finite(preset_max_y)) thr_y <- base::min(thr_y, preset_max_y)
    message("  [gate_mindensity_2d] ext_ctrl | ", ch_x, "=", round(thr_x,1),
            " | ", ch_y, "=", round(thr_y,1), " (sn: ", sn, ")")
  } else {
    thr_x <- .get_thr2d(fr, ch_x, min_x, max_x, gate_range_x, peaks_x, positive_x,
                        upper_mult_x, preset_min_x, preset_max_x)
    thr_y <- .get_thr2d(fr, ch_y, min_y, max_y, gate_range_y, peaks_y, positive_y,
                        upper_mult_y, preset_min_y, preset_max_y)
    message("  [gate_mindensity_2d] per_sample | ", ch_x, "=", round(thr_x,1),
            " | ", ch_y, "=", round(thr_y,1), " (sn: ", sn, ")")
  }
  qx <- if (!is.null(name_x)) name_x else channels[1]
  qy <- if (!is.null(name_y)) name_y else channels[2]
  quad_names <- c(paste0(qx,"-",qy,"+"), paste0(qx,"+",qy,"+"),
                  paste0(qx,"+",qy,"-"), paste0(qx,"-",qy,"-"))
  assign(paste0(filterId, ":quad_names"), quad_names, envir=.fasta_gate_cache)
  quadGate(.gate=setNames(list(thr_x, thr_y), channels), filterId=filterId)
}


# 2D quadrant gating methods handled by the inline quad branch (not openCyto's
# gatingTemplate). Single source of truth -- referenced everywhere a gate must
# be recognized as a 2D quad method (rows_for_gt exclusion, inline dispatch,
# precompute control_col exemption, plot annotation, gating QC). Add a new 2D
# quad method name here and to its precompute/plugin handling.
.FASTA_QUAD_2D_METHODS <- c("gate_flowmeans_2d", "gate_FMO_2d", "gate_tmix_2d", "gate_mindensity_2d")

# Parse a gate_logic 'gates' argument into a clean character vector. Accepts both
# the comma-separated-string form  gates='CD3+,CD19+'  and the legacy vector form
# gates=c('CD3+','CD19+'). A single name 'CD3+' returns c('CD3+'). Used by both
# apply_gating and apply_gating2 so either template style works in either engine.
.fasta_parse_pop_list <- function(g) {
  if (is.null(g)) return(character(0))
  if (length(g) > 1) return(trimws(as.character(g)))                  # c('A','B', ...)
  parts <- trimws(strsplit(as.character(g), ",", fixed = TRUE)[[1]])  # 'A,B,...' or 'A'
  parts[nzchar(parts)]
}

apply_plotting <-
function (gs, gating_template) 
{
    .mk_label <- function(ch, mk_map) {
        mk <- if (ch %in% names(mk_map) && !is.na(mk_map[ch])) 
            mk_map[ch]
        else NULL
        if (!is.null(mk) && mk != ch) 
            paste0(mk, "\n(", ch, ")")
        else ch
    }
    .fasta_plot_theme <- theme_bw() + theme(legend.position = "none", 
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1, 
            face = "plain"), strip.text = element_text(size = 6,
            face = "bold", color = "white"), strip.background = element_rect(color = "navy",
            fill = "navy"), plot.subtitle = element_text(size = 6.5, 
            face = "italic", color = "darkgrey"), plot.title = element_text(size = 10, 
            face = "bold"), plot.caption = element_text(size = 6, 
            color = "grey40", hjust = 0))
    gt_df <- as.data.frame(gating_template)
    if (!"PlotType" %in% colnames(gt_df)) 
        stop("[apply_plotting] 'PlotType' column not found in gating_template.")
    .ag_keep <- if ("active_gate" %in% colnames(gt_df)) 
        is.na(gt_df$active_gate) | suppressWarnings(as.logical(gt_df$active_gate)) != 
            FALSE
    else rep(TRUE, nrow(gt_df))
    pt_df <- gt_df[!is.na(gt_df$PlotType) & nzchar(trimws(gt_df$PlotType)) & 
        .ag_keep, , drop = FALSE]
    if (nrow(pt_df) == 0) {
        message("[apply_plotting] No rows with PlotType in gating_template -- nothing to plot.")
        return(list())
    }
    pt_df$gate_alias <- pt_df$alias
    all_pops <- gs_get_pop_paths(gs)
    all_mk <- markernames(gs[[1]])
    all_ch <- colnames(gs[[1]])
    resolve_ch <- function(name) {
        if (name %in% all_ch) 
            return(name)
        n <- names(all_mk)[match(name, all_mk)]
        if (!is.na(n)) 
            return(n)
        stop("[apply_plotting] Cannot resolve channel: '", name, 
            "'")
    }
    resolve_pop <- function(alias, pops = all_pops) {
        if (alias == "root") 
            return("root")
        if (grepl("/", alias, fixed = TRUE)) {
            parts <- Filter(nzchar, trimws(strsplit(alias, "/")[[1]]))
            last_p <- parts[length(parts)]
            cands <- pops[basename(pops) == last_p]
            if (length(parts) > 1 && length(cands) > 0) {
                prior_path <- paste(parts[-length(parts)], collapse = "/")
                cands <- cands[endsWith(dirname(cands), prior_path)]
            }
            if (length(cands) == 0) 
                stop("[apply_plotting] Path '", alias, "' not found in GatingSet.")
            if (length(cands) > 1) {
                warning("[apply_plotting] Multiple matches for path '", 
                  alias, "' -- using first: ", cands[1])
                cands <- cands[1]
            }
            return(cands)
        }
        m <- pops[basename(pops) == alias]
        if (length(m) == 1) 
            return(m)
        if (length(m) > 1) {
            warning("[apply_plotting] Multiple pops match '", 
                alias, "' -- using first.")
            return(m[1])
        }
        stop("[apply_plotting] Population '", alias, "' not in GatingSet.")
    }
    li_col <- names(pt_df)[tolower(names(pt_df)) == "layout_index"]
    if (length(li_col) == 0) 
        pt_df$Layout_index <- 1L
    else if (li_col != "Layout_index") {
        pt_df$Layout_index <- pt_df[[li_col]]
        pt_df[[li_col]] <- NULL
    }
    pt_df$Layout_index <- suppressWarnings(as.integer(pt_df$Layout_index))
    pt_df$Layout_index[is.na(pt_df$Layout_index)] <- 1L
    layout_ids <- sort(unique(pt_df$Layout_index))
    # QC: detect duplicate layout positions (same layout_index + row + col)
    if (all(c("layout_row", "layout_col", "Layout_index") %in% names(pt_df))) {
        .lpos <- pt_df[!is.na(pt_df$layout_row) & !is.na(pt_df$layout_col),
                       c("alias", "layout_row", "layout_col", "Layout_index")]
        .lpos$pos_key <- paste(.lpos$Layout_index, .lpos$layout_row, .lpos$layout_col, sep=":")
        .dups <- .lpos[duplicated(.lpos$pos_key) | duplicated(.lpos$pos_key, fromLast=TRUE), ]
        if (nrow(.dups) > 0) {
            .dup_msg <- paste(apply(.dups, 1, function(r)
                sprintf("  '%s' -> layout %s, row %s, col %s", r["alias"], r["Layout_index"], r["layout_row"], r["layout_col"])),
                collapse = "\n")
            stop("[apply_plotting] Duplicate layout positions detected.\n",
                 "Multiple gates share the same layout_index + layout_row + layout_col.\n",
                 "Only the last entry will be rendered; earlier ones will be dropped.\n",
                 "Fix the layout_row/layout_col values for these gates:\n",
                 .dup_msg)
        }
    }
    build_design <- function(rows_df) {
        if (!all(c("layout_col", "layout_row") %in% tolower(colnames(rows_df)))) 
            return(NULL)
        r_vals <- suppressWarnings(as.integer(rows_df[[names(rows_df)[tolower(names(rows_df)) == 
            "layout_row"]]]))
        c_vals <- suppressWarnings(as.integer(rows_df[[names(rows_df)[tolower(names(rows_df)) == 
            "layout_col"]]]))
        if (all(is.na(r_vals)) || all(is.na(c_vals))) 
            return(NULL)
        max_r <- max(r_vals, na.rm = TRUE)
        max_c <- max(c_vals, na.rm = TRUE)
        mat <- matrix("#", nrow = max_r, ncol = max_c)
        ids <- LETTERS[seq_len(nrow(rows_df))]
        for (k in seq_len(nrow(rows_df))) {
            if (!is.na(r_vals[k]) && !is.na(c_vals[k])) 
                mat[r_vals[k], c_vals[k]] <- ids[k]
        }
        paste(apply(mat, 1, paste, collapse = ""), collapse = "\n")
    }
    plot_list <- setNames(vector("list", length(gs)), sampleNames(gs))
    plot_params_cache <- list()
    for (.pi in seq_len(nrow(pt_df))) {
        .ptr <- pt_df[.pi, ]
        .alias <- trimws(.ptr$alias)
        .gating_method <- trimws(.ptr$gating_method)
        if (.gating_method == "copy_gate" && !is.na(.ptr$gating_args) && 
            nzchar(trimws(.ptr$gating_args))) {
            .cp2 <- tryCatch(eval(parse(text = paste0("list(", 
                gsub("\"", "'", .ptr$gating_args, fixed = TRUE), 
                ")"))), error = function(e) list())
            .src2 <- if (!is.null(.cp2$copy_from)) 
                pt_df[pt_df$alias == .cp2$copy_from, , drop = FALSE]
            else NULL
            if (!is.null(.src2) && nrow(.src2) > 0) {
                .src_meth2 <- trimws(.src2$gating_method[1])
                .ptr <- .src2[1, ]
                .alias <- trimws(pt_df[.pi, "alias"])
            }
        }
        .dims_r <- trimws(strsplit(.ptr$dims[1], ",")[[1]])
        .ch <- tryCatch(sapply(.dims_r, resolve_ch, USE.NAMES = FALSE), 
            error = function(e) NULL)
        .par <- tryCatch(resolve_pop(.ptr$parent), error = function(e) NULL)
        if (!is.null(.ch) && !is.null(.par)) {
            .xch <- .ch[1]
            .ych <- if (length(.ch) >= 2) 
                .ch[2]
            else "FSC-A"
            plot_params_cache[[.alias]] <- tryCatch(computePlotParams(gs, 
                basename(.par), .xch, .ych, adaptiveBins(gs, 
                  basename(.par))), error = function(e) NULL)
        }
    }
    message("[apply_plotting] Axis limits cached for ", length(plot_params_cache), 
        " gate(s).")
    .debug_cache2 <<- plot_params_cache
    .debug_cache <<- plot_params_cache
    for (s in seq_along(gs)) {
        sn <- sampleNames(gs)[s]
        pd_row <- pData(gs)[sn, , drop = FALSE]
        layout_list <- setNames(vector("list", length(layout_ids)), 
            paste0("layout_", layout_ids))
        for (lid in layout_ids) {
            rows_in_layout <- pt_df[pt_df$Layout_index == lid, 
                , drop = FALSE]
            panel_ids <- LETTERS[seq_len(nrow(rows_in_layout))]
            design_str <- build_design(rows_in_layout)
            # Derive caption wrap width from number of layout columns (33 chars per col)
            .lc_col_nm <- names(rows_in_layout)[tolower(names(rows_in_layout)) == "layout_col"]
            .n_layout_cols <- if (length(.lc_col_nm) > 0) {
                .cv <- suppressWarnings(as.integer(rows_in_layout[[.lc_col_nm[1]]]))
                max(1L, max(.cv, na.rm = TRUE))
            } else nrow(rows_in_layout)
            .caption_wrap_width <- max(40L, as.integer(250L / .n_layout_cols))
            panels <- setNames(vector("list", nrow(rows_in_layout)), 
                panel_ids)
            for (i in seq_len(nrow(rows_in_layout))) {
                pt_row <- rows_in_layout[i, ]
                alias <- trimws(pt_row$alias)
                plot_type <- tolower(trimws(pt_row$PlotType))
                lbl_col <- if ("label_col" %in% colnames(pt_row) && 
                  !is.na(pt_row$label_col)) 
                  trimws(pt_row$label_col)
                else NULL
                subtitle <- if (!is.null(lbl_col) && lbl_col %in% 
                  colnames(pd_row)) 
                  as.character(pd_row[[lbl_col]])
                else sn
                if (pt_row$gating_method == "gate_logic" && !is.na(pt_row$gating_args) && 
                  nzchar(trimws(pt_row$gating_args))) {
                  .gl_args <- tryCatch(eval(parse(text = paste0("list(", 
                    trimws(pt_row$gating_args), ")"))), error = function(e) list())
                  .gl_gates <- .gl_args$gates
                  if (!is.null(.gl_gates) && length(.gl_gates) > 
                    0) {
                    .src_dims <- sapply(.gl_gates, function(g) {
                      .sr <- gt_df[gt_df$alias == trimws(g), 
                        , drop = FALSE]
                      if (nrow(.sr) > 0) 
                        trimws(.sr$dims[1])
                      else NA_character_
                    })
                    .src_dims <- .src_dims[!is.na(.src_dims)]
                    if (length(.src_dims) > 0 && length(unique(.src_dims)) == 
                      1) {
                      pt_row$dims <- .src_dims[1]
                    }
                    else if (!is.na(pt_row$dims) && nzchar(trimws(pt_row$dims))) {
                    }
                    else {
                      panels[[panel_ids[i]]] <- patchwork::plot_spacer()
                      next
                    }
                  }
                }
                # layout_dims (comma-sep markernames) overrides the plotting axes
                # for this row; an unresolvable marker is a hard error. Otherwise
                # fall back to the gating dims (warn -> spacer on failure).
                .ld_raw <- if ("layout_dims" %in% names(pt_row) && !is.na(pt_row$layout_dims) &&
                    nzchar(trimws(as.character(pt_row$layout_dims))))
                    trimws(as.character(pt_row$layout_dims)) else NA_character_
                if (!is.na(.ld_raw)) {
                  dims_raw <- trimws(strsplit(.ld_raw, ",")[[1]])
                  channels <- sapply(dims_raw, function(.d) tryCatch(resolve_ch(.d),
                    error = function(e) stop("[apply_plotting] layout_dims: cannot resolve '",
                      .d, "' to a channel (row '", alias, "').", call. = FALSE)),
                    USE.NAMES = FALSE)
                } else {
                  dims_raw <- trimws(strsplit(pt_row$dims[1], ",")[[1]])
                  channels <- tryCatch(sapply(dims_raw, resolve_ch,
                    USE.NAMES = FALSE), error = function(e) {
                    warning(conditionMessage(e))
                    NULL
                  })
                  if (is.null(channels)) {
                    panels[[panel_ids[i]]] <- patchwork::plot_spacer()
                    next
                  }
                }
                parent_path <- tryCatch(resolve_pop(pt_row$parent), 
                  error = function(e) {
                    warning(conditionMessage(e))
                    NULL
                  })
                gate_path <- tryCatch(resolve_pop(alias), error = function(e) NULL)
                if (is.null(parent_path)) {
                  panels[[panel_ids[i]]] <- patchwork::plot_spacer()
                  next
                }
                if (plot_type == "scatter") {
                  x_ch <- channels[1]
                  y_ch <- if (length(channels) >= 2) 
                    channels[2]
                  else "FSC-A"
                  args_str_plot <- if ("gating_args" %in% colnames(pt_row) && 
                    !is.na(pt_row$gating_args) && nzchar(trimws(pt_row$gating_args))) 
                    trimws(pt_row$gating_args)
                  else ""
                  # layout_labels: parse comma-sep values (gate_name, percent, count)
                  .lbl_types <- if ("layout_labels" %in% names(pt_row) &&
                      !is.na(pt_row$layout_labels) && nzchar(trimws(as.character(pt_row$layout_labels)))) {
                    trimws(strsplit(tolower(as.character(pt_row$layout_labels)), "[,|/ ]+")[[1]])
                  } else c("gate_name", "percent", "count")
                  .mk_stat_lbl <- function(.name, .pct, .cnt) {
                    parts <- c(
                      if ("gate_name" %in% .lbl_types) as.character(.name),
                      if ("percent"   %in% .lbl_types) paste0(round(.pct, 1), "%"),
                      if ("count"     %in% .lbl_types) as.character(.cnt)
                    )
                    paste(parts[!is.na(parts) & nzchar(parts)], collapse = "\n")
                  }
                  .params_key <- if (pt_row$gating_method == 
                    "copy_gate" && nzchar(trimws(args_str_plot))) {
                    .cp <- tryCatch(eval(parse(text = paste0("list(", 
                      trimws(args_str_plot), ")"))), error = function(e) list())
                    if (!is.null(.cp$copy_from) && !is.null(plot_params_cache[[.cp$copy_from]])) 
                      .cp$copy_from
                    else alias
                  }
                  else alias
                  # Bypass the alias-keyed cache when layout_dims overrode the
                  # axes -- the cached params were computed for the default dims
                  # and would give wrong limits/binwidth for the new channels.
                  params <- if (is.na(.ld_raw) && !is.null(plot_params_cache[[.params_key]]))
                    plot_params_cache[[.params_key]]
                  else computePlotParams(gs[s], basename(parent_path),
                    x_ch, y_ch, adaptiveBins(gs[s], basename(parent_path)))
                  computePlotParams(gs[s], basename(parent_path), 
                    x_ch, y_ch, adaptiveBins(gs[s], basename(parent_path)))
                  imm_children <- if (!is.null(gate_path)) 
                    all_pops[dirname(all_pops) == gate_path]
                  else character(0)
                  is_quad_children <- length(imm_children) > 
                    0 && all(grepl("[+-]", basename(imm_children)))
                  overlay_pops <- if (!is.null(gate_path)) {
                    if (is_quad_children)
                      imm_children
                    else gate_path
                  }
                  else NULL
                  plot_title <- paste0("Alias: ", alias)
                  # Build grandparent/parent display
                  .par_str <- trimws(as.character(pt_row$parent))
                  .par_parts <- strsplit(.par_str, "/")[[1]]
                  .sub_parent <- if (length(.par_parts) > 1) {
                    paste(tail(.par_parts, 2), collapse = "/")
                  } else {
                    .gp_row <- gt_df[gt_df$alias == .par_str, , drop = FALSE]
                    .gpar <- if (nrow(.gp_row) > 0 && !is.na(.gp_row$parent[1]) &&
                        nzchar(.gp_row$parent[1]) && .gp_row$parent[1] != "root")
                        basename(trimws(.gp_row$parent[1])) else NULL
                    if (!is.null(.gpar)) paste0(.gpar, "/", .par_str) else .par_str
                  }
                  plot_subtitle <- paste0(.sub_parent, " | Method: ", pt_row$gating_method)
                  .grpby_val <- if ("groupBy" %in% colnames(pt_row) && 
                    !is.na(pt_row$groupBy) && nzchar(trimws(pt_row$groupBy))) 
                    trimws(pt_row$groupBy)
                  else NULL
                  raw_cap <- if (nzchar(trimws(args_str_plot))) 
                    paste0("Args: ", args_str_plot)
                  else NULL
                  if (!is.null(.grpby_val))
                    raw_cap <- paste0(if (!is.null(raw_cap))
                      raw_cap
                    else "", if (!is.null(raw_cap))
                      " | "
                    else "", "groupBy: ", .grpby_val)
                  # Indicate collapseDataForGating only when TRUE (one pooled gate
                  # replicated across the group); blank/FALSE adds nothing.
                  if ("collapseDataForGating" %in% colnames(pt_row) &&
                      !is.na(pt_row$collapseDataForGating) &&
                      isTRUE(suppressWarnings(as.logical(pt_row$collapseDataForGating))))
                    raw_cap <- paste0(if (!is.null(raw_cap)) raw_cap else "",
                                      if (!is.null(raw_cap)) " | " else "",
                                      "collapseDataForGating: TRUE")
                  plot_caption <- if (!is.null(raw_cap))
                    paste(strwrap(raw_cap, width = .caption_wrap_width), collapse = "\n")
                  else NULL
                  .src_meth_q <- if (pt_row$gating_method == 
                    "copy_gate" && nzchar(trimws(args_str_plot))) {
                    .cpa_q <- tryCatch(eval(parse(text = paste0("list(", 
                      trimws(args_str_plot), ")"))), error = function(e) list())
                    if (!is.null(.cpa_q$copy_from)) {
                      .sr_q <- gt_df[gt_df$alias == .cpa_q$copy_from, 
                        , drop = FALSE]
                      if (nrow(.sr_q) > 0) 
                        trimws(.sr_q$gating_method[1])
                      else ""
                    }
                    else ""
                  }
                  else ""
                  is_quad_method <- pt_row$gating_method %in%
                    .FASTA_QUAD_2D_METHODS ||
                    (pt_row$gating_method == "copy_gate" && .src_meth_q %in%
                      .FASTA_QUAD_2D_METHODS)
                  if (pt_row$gating_method == "copy_gate" && 
                    nzchar(trimws(args_str_plot))) {
                    .cpargs <- tryCatch(eval(parse(text = paste0("list(", 
                      trimws(args_str_plot), ")"))), error = function(e) list())
                    .src_alias <- .cpargs$copy_from
                    if (!is.null(.src_alias)) {
                      .src_r <- gt_df[gt_df$alias == .src_alias, 
                        , drop = FALSE]
                      if (nrow(.src_r) > 0 && !is.na(.src_r$gating_args[1])) {
                        .src_args <- tryCatch(eval(parse(text = paste0("list(", 
                          .src_r$gating_args[1], ")"))), error = function(e) list())
                        if (!is.null(.src_args$name_x)) 
                          pt_row$name_x_override <- .src_args$name_x
                        if (!is.null(.src_args$name_y)) 
                          pt_row$name_y_override <- .src_args$name_y
                      }
                    }
                  }
                  imm_children <- if (!is.null(gate_path)) 
                    all_pops[dirname(all_pops) == gate_path]
                  else if (is_quad_method) 
                    all_pops[dirname(all_pops) == parent_path]
                  else character(0)
                  if (is_quad_method && length(imm_children) > 
                    0 && "gating_args" %in% colnames(pt_row) && 
                    !is.na(pt_row$gating_args) && nzchar(trimws(pt_row$gating_args))) {
                    .args_q <- tryCatch(eval(parse(text = paste0("list(", 
                      trimws(pt_row$gating_args), ")"))), error = function(e) list())
                    .nx <- if (!is.null(pt_row$name_x_override)) 
                      pt_row$name_x_override
                    else .args_q$name_x
                    .ny <- if (!is.null(pt_row$name_y_override)) 
                      pt_row$name_y_override
                    else .args_q$name_y
                    if (!is.null(.nx) && !is.null(.ny)) {
                      .pat_x <- gsub("[^a-zA-Z0-9]", ".", .nx)
                      .pat_y <- gsub("[^a-zA-Z0-9]", ".", .ny)
                      .matched <- grepl(.pat_x, basename(imm_children)) & 
                        grepl(.pat_y, basename(imm_children))
                      if (any(.matched)) 
                        imm_children <- imm_children[.matched]
                    }
                  }
                  overlay_pops <- if (is_quad_method && length(imm_children) > 
                    0) 
                    imm_children
                  else if (!is.null(gate_path)) 
                    gate_path
                  else NULL
                  if (pt_row$gating_method == "gate_logic") {
                    .gl_args2 <- tryCatch(eval(parse(text = paste0("list(",
                      trimws(pt_row$gating_args), ")"))), error = function(e) list())
                    .gl_gates2 <- .gl_args2$gates
                    if (!is.null(.gl_gates2)) {
                      .gl_paths <- tryCatch(sapply(trimws(.gl_gates2),
                        function(g) resolve_pop(g)), error = function(e) NULL)
                      if (!is.null(.gl_paths))
                        overlay_pops <- .gl_paths
                    }
                  }
                  # layout_gates (comma-sep gate/population names) overrides the
                  # overlay set computed above -- resolved to full paths here
                  # (unresolvable -> hard error). The draw passes render a listed
                  # gate only where its dims match the panel axes; non-matching
                  # ones are silently skipped. When set, it also suppresses the
                  # default primary-gate annotation so ONLY the listed gates show.
                  .lg_override <- FALSE
                  .lg_raw <- if ("layout_gates" %in% names(pt_row) && !is.na(pt_row$layout_gates) &&
                      nzchar(trimws(as.character(pt_row$layout_gates))))
                      trimws(as.character(pt_row$layout_gates)) else NA_character_
                  if (!is.na(.lg_raw)) {
                    .lg_names <- Filter(nzchar, trimws(strsplit(.lg_raw, ",")[[1]]))
                    overlay_pops <- vapply(.lg_names, function(.g) tryCatch(resolve_pop(.g),
                      error = function(e) stop("[apply_plotting] layout_gates: cannot resolve gate '",
                        .g, "' to a population (row '", alias, "').", call. = FALSE)),
                      character(1), USE.NAMES = FALSE)
                    .lg_override <- TRUE
                  } else {
                    # 1D split row (pop='+-'/'-+'): the row creates "<alias>+" and
                    # "<alias>-" -- overlay + label both (shared threshold line).
                    .pop_pm <- if ("pop" %in% names(pt_row) && !is.na(pt_row$pop))
                        trimws(as.character(pt_row$pop)) else "+"
                    if (.pop_pm %in% c("+-", "-+")) {
                      .pm <- tryCatch(c(resolve_pop(paste0(alias, "+")), resolve_pop(paste0(alias, "-"))),
                                      error = function(e) NULL)
                      if (!is.null(.pm)) { overlay_pops <- .pm; .lg_override <- TRUE }
                    }
                  }
                  fr_s <- gh_pop_get_data(gs[[sn]], parent_path)
                  dat_s <- as.data.frame(exprs(fr_s)[, c(x_ch, 
                    y_ch), drop = FALSE])
                  colnames(dat_s) <- c("x_var", "y_var")
                  if (nrow(dat_s) == 0) {
                    message("[apply_plotting] Empty panel: alias='", 
                      alias, "' sample='", sn, "'")
                    dat_empty <- data.frame(x_var = NA_real_, 
                      y_var = NA_real_, label_facet = subtitle)
                    p_empty <- ggplot(dat_empty, aes(x = x_var, 
                      y = y_var)) + facet_wrap(~label_facet) + 
                      coord_cartesian(xlim = params$xlim, ylim = params$ylim) + 
                      .fasta_plot_theme + labs(x = .mk_label(x_ch, 
                      all_mk), y = .mk_label(y_ch, all_mk), title = plot_title, 
                      subtitle = plot_subtitle, caption = plot_caption) + 
                      theme(panel.background = element_rect(fill = "grey97"))
                    if (is_quad_method && length(overlay_pops) > 
                      1) {
                      .qxs <- .qys <- c()
                      for (.gp in overlay_pops) {
                        .gc <- tryCatch(gh_pop_get_gate(gs[[sn]], 
                          .gp), error = function(e) NULL)
                        if (is.null(.gc) || !inherits(.gc, "rectangleGate")) 
                          next
                        if (x_ch %in% names(.gc@min) && is.finite(.gc@min[[x_ch]])) 
                          .qxs <- c(.qxs, .gc@min[[x_ch]])
                        if (x_ch %in% names(.gc@max) && is.finite(.gc@max[[x_ch]])) 
                          .qxs <- c(.qxs, .gc@max[[x_ch]])
                        if (y_ch %in% names(.gc@min) && is.finite(.gc@min[[y_ch]])) 
                          .qys <- c(.qys, .gc@min[[y_ch]])
                        if (y_ch %in% names(.gc@max) && is.finite(.gc@max[[y_ch]])) 
                          .qys <- c(.qys, .gc@max[[y_ch]])
                      }
                      .qx <- if (length(.qxs) > 0) 
                        median(.qxs)
                      else mean(params$xlim)
                      .qy <- if (length(.qys) > 0) 
                        median(.qys)
                      else mean(params$ylim)
                      p_empty <- p_empty + geom_vline(xintercept = .qx, 
                        color = "red", linewidth = 0.7) + geom_hline(yintercept = .qy, 
                        color = "red", linewidth = 0.7)
                      .xlim <- params$xlim
                      .ylim <- params$ylim
                      .nx_ov <- if (!is.null(pt_row$name_x_override)) 
                        pt_row$name_x_override
                      else x_ch
                      .ny_ov <- if (!is.null(pt_row$name_y_override)) 
                        pt_row$name_y_override
                      else y_ch
                      for (.gp in overlay_pops) {
                        .gn <- basename(.gp)
                        .nx_pat <- gsub("[^a-zA-Z0-9]", ".", 
                          .nx_ov)
                        .ny_pat <- gsub("[^a-zA-Z0-9]", ".", 
                          .ny_ov)
                        .xpos <- grepl(paste0(.nx_pat, "\\+"), 
                          .gn)
                        .ypos <- grepl(paste0(.ny_pat, "\\+"), 
                          .gn)
                        .lx <- if (.xpos) 
                          mean(c(.qx, .xlim[2]))
                        else mean(c(.xlim[1], .qx))
                        .ly <- if (.ypos) 
                          mean(c(.qy, .ylim[2]))
                        else mean(c(.ylim[1], .qy))
                        p_empty <- p_empty + annotate("text", 
                          x = .lx, y = .ly, label = paste0(.gn, 
                            "\n0\nNA%"), size = 2.5, color = "grey20", 
                          hjust = 0.5, fontface = "bold")
                      }
                    }
                    else if (!is.null(overlay_pops)) {
                      for (.gp in overlay_pops) {
                        .gc2 <- tryCatch(gh_pop_get_gate(gs[[sn]], 
                          .gp), error = function(e) NULL)
                        if (is.null(.gc2)) 
                          next
                        if (inherits(.gc2, "rectangleGate") && 
                          x_ch %in% names(.gc2@min) && is.finite(.gc2@min[[x_ch]])) 
                          p_empty <- p_empty + geom_vline(xintercept = .gc2@min[[x_ch]], 
                            color = "red", linewidth = 0.7, linetype = "dashed")
                      }
                      p_empty <- p_empty + annotate("text", x = mean(params$xlim), 
                        y = mean(params$ylim), label = "0 events\nNA%", 
                        size = 4, color = "grey50", fontface = "italic", 
                        hjust = 0.5)
                    }
                    else {
                      p_empty <- p_empty + annotate("text", x = mean(params$xlim), 
                        y = mean(params$ylim), label = "0 events\nNA%", 
                        size = 4, color = "grey50", fontface = "italic", 
                        hjust = 0.5)
                    }
                    panels[[panel_ids[i]]] <- p_empty
                    next
                  }
                  dat_s$label_facet <- subtitle
                  .plot_xlim <- params$xlim
                  .plot_ylim <- params$ylim
                  if (!is.null(gate_path)) {
                    .g_tmp <- tryCatch(gh_pop_get_gate(gs[[sn]], 
                      gate_path), error = function(e) NULL)
                    if (!is.null(.g_tmp) && inherits(.g_tmp,
                      "polygonGate")) {
                      .bnd_tmp <- as.data.frame(.g_tmp@boundaries)
                      # Ignore the open-edge sentinels (flowCore stores +/-Inf as
                      # +/-2147483647); only real cutpoints should widen the view,
                      # otherwise the axes blow out to billions (e.g. gate_tmix_2d
                      # polygons). Keep values within a sane fluorescence range.
                      .bx <- .bnd_tmp[is.finite(.bnd_tmp[, 1]) & abs(.bnd_tmp[, 1]) < 1e7, 1]
                      .by <- .bnd_tmp[is.finite(.bnd_tmp[, 2]) & abs(.bnd_tmp[, 2]) < 1e7, 2]
                      if (length(.bx)) {
                        .plot_xlim[1] <- min(.plot_xlim[1], min(.bx))
                        .plot_xlim[2] <- max(.plot_xlim[2], max(.bx))
                      }
                      if (length(.by)) {
                        .plot_ylim[1] <- min(.plot_ylim[1], min(.by))
                        .plot_ylim[2] <- max(.plot_ylim[2], max(.by))
                      }
                    }
                  }
                  p <- ggplot(dat_s, aes(x = x_var, y = y_var)) + 
                    geom_bin2d(binwidth = params$binwidth) + 
                    scale_fill_gradientn(colours = c("navy", 
                      "#1E90FF", "cyan", "#32CD32", "yellow", 
                      "orange", "red"), trans = "log10", na.value = "navy", 
                      name = "Count") + coord_cartesian(xlim = .plot_xlim, 
                    ylim = .plot_ylim, clip = "on") + facet_wrap(~label_facet) + 
                    .fasta_plot_theme + labs(x = .mk_label(x_ch, 
                    all_mk), y = .mk_label(y_ch, all_mk), title = plot_title, 
                    subtitle = plot_subtitle, caption = plot_caption)
                  # gate_logic threshold helper (shared by the primary gate_logic
                  # crosshair AND the layout_gates overlay path): walk a gate_logic's
                  # component pops to the first finite x/y rectangleGate bound,
                  # recursing one level into a referenced gate_logic so NOT(CAR_all)
                  # still finds the CAR thresholds. Returns c(tx, ty).
                  .gl_xh <- function(pops, depth = 0) {
                    tx <- NA_real_; ty <- NA_real_
                    if (depth > 3 || length(pops) == 0) return(c(tx, ty))
                    for (pp in pops) {
                      gg <- tryCatch(gh_pop_get_gate(gs[[sn]], resolve_pop(pp)),
                                     error = function(e) NULL)
                      if (is.null(gg) || !inherits(gg, "rectangleGate")) next
                      if (is.na(tx) && x_ch %in% names(gg@min)) {
                        .vv <- c(gg@min[[x_ch]], gg@max[[x_ch]]); .vv <- .vv[is.finite(.vv)]
                        if (length(.vv)) tx <- .vv[1]
                      }
                      if (is.na(ty) && y_ch %in% names(gg@min)) {
                        .vv <- c(gg@min[[y_ch]], gg@max[[y_ch]]); .vv <- .vv[is.finite(.vv)]
                        if (length(.vv)) ty <- .vv[1]
                      }
                      if ((is.na(tx) || is.na(ty)) && all(!is.finite(c(gg@min, gg@max)))) {
                        .srr <- gt_df[gt_df$alias == basename(resolve_pop(pp)), , drop = FALSE]
                        if (nrow(.srr) > 0) {
                          .saa <- tryCatch(eval(parse(text = paste0("list(",
                            .srr$gating_args[1], ")"))), error = function(e) list())
                          .rr <- .gl_xh(.fasta_parse_pop_list(.saa$gates), depth + 1)
                          if (is.na(tx)) tx <- .rr[1]
                          if (is.na(ty)) ty <- .rr[2]
                        }
                      }
                      if (!is.na(tx) && !is.na(ty)) break
                    }
                    c(tx, ty)
                  }
                  # Auto gate-boundary annotation for the primary population
                  # (suppressed when layout_gates overrides the overlay set).
                  if (!is.null(gate_path) && !.lg_override) {
                    .g_auto <- tryCatch(gh_pop_get_gate(gs[[sn]], gate_path),
                        error = function(e) NULL)
                    if (!is.null(.g_auto)) {
                      # quadGate child: use quadintersection for both crosshairs
                      .qi_auto <- attr(.g_auto, "quadintersection")
                      if (!is.null(.qi_auto)) {
                        .vl_q <- if (x_ch %in% names(.qi_auto) && is.finite(.qi_auto[[x_ch]]))
                            .qi_auto[[x_ch]] else NA_real_
                        .hl_q <- if (y_ch %in% names(.qi_auto) && is.finite(.qi_auto[[y_ch]]))
                            .qi_auto[[y_ch]] else NA_real_
                        if (!is.na(.vl_q)) {
                          p <- p + geom_vline(data = data.frame(xint = .vl_q, label_facet = subtitle),
                              aes(xintercept = xint), color = "red", linewidth = 0.7, inherit.aes = FALSE)
                          p <- p + geom_text(data = data.frame(xint = pmin(.vl_q, .plot_xlim[2] - diff(.plot_xlim)*0.02) - diff(.plot_xlim)*0.015, lbl = round(.vl_q), label_facet = subtitle),
                              aes(x = xint, label = lbl),
                              y = .plot_ylim[2] - diff(.plot_ylim)*0.05,
                              vjust = 0, hjust = 1, angle = 90, size = 2.8,
                              color = "red", inherit.aes = FALSE)
                        }
                        if (!is.na(.hl_q)) {
                          p <- p + geom_hline(data = data.frame(yint = .hl_q, label_facet = subtitle),
                              aes(yintercept = yint), color = "red", linewidth = 0.7, inherit.aes = FALSE)
                          p <- p + geom_text(data = data.frame(yv = pmin(.hl_q + diff(.plot_ylim)*0.020, .plot_ylim[2] - diff(.plot_ylim)*0.04), lbl = round(.hl_q), label_facet = subtitle),
                              aes(y = yv, label = lbl),
                              x = .plot_xlim[2] - diff(.plot_xlim)*0.03,
                              vjust = 0, hjust = 1, size = 2.5,
                              color = "red", inherit.aes = FALSE)
                        }
                      } else if (inherits(.g_auto, "rectangleGate")) {
                        # rectangleGate: a 2D box is drawn as an actual rectangle
                        # (open +/-Inf edges pushed off-view so coord clip removes
                        # them, leaving only the finite boundary edges -- no
                        # full-panel crosshair); a 1D gate stays a single line.
                        .xmn <- if (x_ch %in% names(.g_auto@min)) .g_auto@min[[x_ch]] else -Inf
                        .xmx <- if (x_ch %in% names(.g_auto@max)) .g_auto@max[[x_ch]] else  Inf
                        .ymn <- if (y_ch %in% names(.g_auto@min)) .g_auto@min[[y_ch]] else -Inf
                        .ymx <- if (y_ch %in% names(.g_auto@max)) .g_auto@max[[y_ch]] else  Inf
                        .x_con <- is.finite(.xmn) || is.finite(.xmx)
                        .y_con <- is.finite(.ymn) || is.finite(.ymx)
                        .xlab_r <- if (is.finite(.xmn)) .xmn else if (is.finite(.xmx)) .xmx else NA_real_
                        .ylab_r <- if (is.finite(.ymn)) .ymn else if (is.finite(.ymx)) .ymx else NA_real_
                        if (.x_con && .y_con) {
                          # 2D box: draw ONLY the finite edges as segments (omit open
                          # +/-Inf edges entirely). Open ends extend well past the view
                          # so coord clip trims them with no perpendicular cap -- a
                          # geom_rect would instead draw a spurious edge at the panel
                          # boundary (coord expansion lands it just inside the view).
                          .vx0 <- .plot_xlim[1] - diff(.plot_xlim); .vx1 <- .plot_xlim[2] + diff(.plot_xlim)
                          .vy0 <- .plot_ylim[1] - diff(.plot_ylim); .vy1 <- .plot_ylim[2] + diff(.plot_ylim)
                          .yseg0 <- if (is.finite(.ymn)) .ymn else .vy0
                          .yseg1 <- if (is.finite(.ymx)) .ymx else .vy1
                          .xseg0 <- if (is.finite(.xmn)) .xmn else .vx0
                          .xseg1 <- if (is.finite(.xmx)) .xmx else .vx1
                          .segs <- data.frame(x = numeric(0), xend = numeric(0), y = numeric(0), yend = numeric(0))
                          if (is.finite(.xmn)) .segs <- rbind(.segs, data.frame(x = .xmn, xend = .xmn, y = .yseg0, yend = .yseg1))  # left
                          if (is.finite(.xmx)) .segs <- rbind(.segs, data.frame(x = .xmx, xend = .xmx, y = .yseg0, yend = .yseg1))  # right
                          if (is.finite(.ymn)) .segs <- rbind(.segs, data.frame(x = .xseg0, xend = .xseg1, y = .ymn, yend = .ymn))  # bottom
                          if (is.finite(.ymx)) .segs <- rbind(.segs, data.frame(x = .xseg0, xend = .xseg1, y = .ymx, yend = .ymx))  # top
                          if (nrow(.segs)) {
                            .segs$label_facet <- subtitle
                            p <- p + geom_segment(data = .segs, aes(x = x, xend = xend, y = y, yend = yend),
                                color = "red", linewidth = 0.7, inherit.aes = FALSE)
                          }
                          if (!is.na(.xlab_r)) {
                            p <- p + geom_text(data = data.frame(xint = pmin(.xlab_r, .plot_xlim[2] - diff(.plot_xlim)*0.02) - diff(.plot_xlim)*0.015, lbl = round(.xlab_r), label_facet = subtitle),
                                aes(x = xint, label = lbl), y = .plot_ylim[2] - diff(.plot_ylim)*0.05,
                                vjust = 0, hjust = 1, angle = 90, size = 2.8, color = "red", inherit.aes = FALSE)
                          }
                          if (!is.na(.ylab_r)) {
                            p <- p + geom_text(data = data.frame(yv = pmin(.ylab_r + diff(.plot_ylim)*0.020, .plot_ylim[2] - diff(.plot_ylim)*0.04), lbl = round(.ylab_r), label_facet = subtitle),
                                aes(y = yv, label = lbl), x = .plot_xlim[2] - diff(.plot_xlim)*0.03,
                                vjust = 0, hjust = 1, size = 2.5, color = "red", inherit.aes = FALSE)
                          }
                        } else {
                          # 1D gate: single threshold line at the finite boundary
                          if (.x_con && !is.na(.xlab_r)) {
                            p <- p + geom_vline(data = data.frame(xint = .xlab_r, label_facet = subtitle),
                                aes(xintercept = xint), color = "red", linewidth = 0.7, inherit.aes = FALSE)
                            p <- p + geom_text(data = data.frame(xint = pmin(.xlab_r, .plot_xlim[2] - diff(.plot_xlim)*0.02) - diff(.plot_xlim)*0.015, lbl = round(.xlab_r), label_facet = subtitle),
                                aes(x = xint, label = lbl), y = .plot_ylim[2] - diff(.plot_ylim)*0.05,
                                vjust = 0, hjust = 1, angle = 90, size = 2.8, color = "red", inherit.aes = FALSE)
                          }
                          if (.y_con && !is.na(.ylab_r)) {
                            p <- p + geom_hline(data = data.frame(yint = .ylab_r, label_facet = subtitle),
                                aes(yintercept = yint), color = "red", linewidth = 0.7, inherit.aes = FALSE)
                            p <- p + geom_text(data = data.frame(yv = pmin(.ylab_r + diff(.plot_ylim)*0.020, .plot_ylim[2] - diff(.plot_ylim)*0.04), lbl = round(.ylab_r), label_facet = subtitle),
                                aes(y = yv, label = lbl), x = .plot_xlim[2] - diff(.plot_xlim)*0.03,
                                vjust = 0, hjust = 1, size = 2.5, color = "red", inherit.aes = FALSE)
                          }
                        }
                      }
                    }
                  # gate_logic crosshair: the dummy gate has no geometry, so derive
                  # the dividing thresholds from the component quadrant gates'
                  # rectangleGate bounds. Recurse one level into a referenced
                  # gate_logic so NOT(CAR_all) still finds the CAR thresholds.
                  # Draw vline/hline (span to infinity). Additive only.
                  # (suppressed when layout_gates overrides the overlay set.)
                  if (pt_row$gating_method == "gate_logic" && !.lg_override) {
                    .gla <- tryCatch(eval(parse(text = paste0("list(",
                      trimws(pt_row$gating_args), ")"))), error = function(e) list())
                    .glc <- .gl_xh(.fasta_parse_pop_list(.gla$gates))
                    if (!is.na(.glc[1])) {
                      p <- p + geom_vline(data = data.frame(xint = .glc[1], label_facet = subtitle),
                          aes(xintercept = xint), color = "red", linewidth = 0.7, inherit.aes = FALSE)
                      p <- p + geom_text(data = data.frame(xint = pmin(.glc[1], .plot_xlim[2] - diff(.plot_xlim)*0.02) - diff(.plot_xlim)*0.015, lbl = round(.glc[1]), label_facet = subtitle),
                          aes(x = xint, label = lbl), y = .plot_ylim[2] - diff(.plot_ylim)*0.05,
                          vjust = 0, hjust = 1, angle = 90, size = 2.8, color = "red", inherit.aes = FALSE)
                    }
                    if (!is.na(.glc[2])) {
                      p <- p + geom_hline(data = data.frame(yint = .glc[2], label_facet = subtitle),
                          aes(yintercept = yint), color = "red", linewidth = 0.7, inherit.aes = FALSE)
                      p <- p + geom_text(data = data.frame(yv = pmin(.glc[2] + diff(.plot_ylim)*0.020, .plot_ylim[2] - diff(.plot_ylim)*0.04), lbl = round(.glc[2]), label_facet = subtitle),
                          aes(y = yv, label = lbl), x = .plot_xlim[2] - diff(.plot_xlim)*0.03,
                          vjust = 0, hjust = 1, size = 2.8, color = "red", inherit.aes = FALSE)
                    }
                  }
                  }
                  if (!is.null(overlay_pops)) {
                    if (is_quad_method && length(overlay_pops) >
                      1 && !.lg_override) {
                      .thr_x <- NA_real_
                      .thr_y <- NA_real_
                      for (.gp in overlay_pops) {
                        .g <- tryCatch(gh_pop_get_gate(gs[[sn]], 
                          .gp), error = function(e) NULL)
                        if (is.null(.g) || !inherits(.g, "rectangleGate")) 
                          next
                        if (x_ch %in% names(.g@min) && is.finite(.g@min[[x_ch]])) 
                          .thr_x <- .g@min[[x_ch]]
                        if (y_ch %in% names(.g@min) && is.finite(.g@min[[y_ch]])) 
                          .thr_y <- .g@min[[y_ch]]
                      }
                      if (!is.na(.thr_x)) {
                        p <- p + geom_vline(data = data.frame(xint = .thr_x, 
                          label_facet = subtitle), aes(xintercept = xint), 
                          color = "red", linewidth = 0.7, inherit.aes = FALSE)
                        p <- p + geom_text(data = data.frame(xint = pmin(.thr_x, .plot_xlim[2] - diff(.plot_xlim)*0.02) - diff(.plot_xlim)*0.015, lbl = round(.thr_x), label_facet = subtitle),
                          aes(x = xint, label = lbl),
                          y = .plot_ylim[2] - diff(.plot_ylim)*0.05,
                          vjust = 0, hjust = 1, angle = 90, size = 2.8,
                          color = "red", inherit.aes = FALSE)
                      }
                      if (!is.na(.thr_y)) {
                        p <- p + geom_hline(data = data.frame(yint = .thr_y, 
                          label_facet = subtitle), aes(yintercept = yint), 
                          color = "red", linewidth = 0.7, inherit.aes = FALSE)
                        p <- p + geom_text(data = data.frame(yv = pmin(.thr_y + diff(.plot_ylim)*0.020, .plot_ylim[2] - diff(.plot_ylim)*0.04), lbl = round(.thr_y), label_facet = subtitle),
                          aes(y = yv, label = lbl),
                          x = .plot_xlim[2] - diff(.plot_xlim)*0.03,
                          vjust = 0, hjust = 1, size = 2.8,
                          color = "red", inherit.aes = FALSE)
                      }
                      # polygonGate quad children (e.g. gate_tmix_2d) have no single
                      # crosshair threshold -- draw each polygon outline instead.
                      # Clamp the open-edge sentinels (+/-2147483647) to just beyond
                      # the view: the far edges then sit off-panel (clipped by
                      # coord_cartesian) and, crucially, don't train the x/y scale
                      # to ~4e9 -- which would make stat_bin2d demand billions of
                      # bins and fail ("number of histogram bins must be < 1e6").
                      .xpad <- diff(.plot_xlim) * 0.05
                      .ypad <- diff(.plot_ylim) * 0.05
                      for (.gp in overlay_pops) {
                        .gpoly <- tryCatch(gh_pop_get_gate(gs[[sn]], .gp),
                          error = function(e) NULL)
                        if (is.null(.gpoly) || !inherits(.gpoly, "polygonGate"))
                          next
                        .pbnd <- as.data.frame(.gpoly@boundaries)
                        if (ncol(.pbnd) < 2 || nrow(.pbnd) < 2) next
                        colnames(.pbnd) <- c("px", "py")
                        .pbnd$px <- pmin(pmax(.pbnd$px, .plot_xlim[1] - .xpad), .plot_xlim[2] + .xpad)
                        .pbnd$py <- pmin(pmax(.pbnd$py, .plot_ylim[1] - .ypad), .plot_ylim[2] + .ypad)
                        .pbnd$label_facet <- subtitle
                        .pbnd <- rbind(.pbnd, .pbnd[1, , drop = FALSE])
                        p <- p + geom_path(data = .pbnd, aes(x = px, y = py),
                          color = "red", linewidth = 0.7, inherit.aes = FALSE)
                      }
                      # Label-placement crosshair: rect-quad thresholds if present;
                      # for tmix polygon quads, derive from each polygon's
                      # center-facing finite corner so every label sits in the
                      # middle of its quadrant (not on sparse events).
                      .lbl_thr_x <- .thr_x; .lbl_thr_y <- .thr_y
                      if (is.na(.lbl_thr_x) || is.na(.lbl_thr_y)) {
                        .txs <- numeric(0); .tys <- numeric(0)
                        for (.gp2 in overlay_pops) {
                          .g2 <- tryCatch(gh_pop_get_gate(gs[[sn]], .gp2), error = function(e) NULL)
                          if (is.null(.g2) || !inherits(.g2, "polygonGate")) next
                          .b2 <- as.data.frame(.g2@boundaries)
                          if (ncol(.b2) < 2) next
                          .bx2 <- .b2[is.finite(.b2[, 1]) & abs(.b2[, 1]) < 1e7, 1]
                          .by2 <- .b2[is.finite(.b2[, 2]) & abs(.b2[, 2]) < 1e7, 2]
                          .sg2 <- regmatches(basename(.gp2), gregexpr("[+-]", basename(.gp2)))[[1]]
                          .sx2 <- if (length(.sg2) >= 1) .sg2[1] else "+"
                          .sy2 <- if (length(.sg2) >= 2) .sg2[2] else "+"
                          if (length(.bx2)) .txs <- c(.txs, if (.sx2 == "+") min(.bx2) else max(.bx2))
                          if (length(.by2)) .tys <- c(.tys, if (.sy2 == "+") min(.by2) else max(.by2))
                        }
                        if (is.na(.lbl_thr_x) && length(.txs)) .lbl_thr_x <- median(.txs)
                        if (is.na(.lbl_thr_y) && length(.tys)) .lbl_thr_y <- median(.tys)
                      }
                      for (.gp in overlay_pops) {
                        pop_ff <- tryCatch(gh_pop_get_data(gs[[sn]],
                          .gp), error = function(e) NULL)
                        parent_ff <- tryCatch(gh_pop_get_data(gs[[sn]], 
                          dirname(.gp)), error = function(e) NULL)
                        if (is.null(pop_ff) || is.null(parent_ff) || 
                          nrow(exprs(parent_ff)) == 0) 
                          next
                        pct_val <- round(100 * nrow(exprs(pop_ff))/nrow(exprs(parent_ff)), 
                          1)
                        cnt_val <- nrow(exprs(pop_ff))
                        pop_e <- exprs(pop_ff)
                        if (!x_ch %in% colnames(pop_e) || !y_ch %in% 
                          colnames(pop_e)) 
                          next
                        .gname <- basename(.gp)
                        .x_pos <- if (grepl(paste0(x_ch, "\\+|\\+.*", 
                          y_ch), .gname) || (is.finite(.thr_x) && 
                          !is.na(.thr_x) && any(grepl("\\+", 
                          regmatches(.gname, gregexpr("[^a-zA-Z0-9_][+-]|^[+-]", 
                            .gname)))))) {
                          if (is.finite(.thr_x)) 
                            (.thr_x + .plot_xlim[2])/2
                          else .plot_xlim[2] * 0.75
                        }
                        else {
                          if (is.finite(.thr_x)) 
                            (.plot_xlim[1] + .thr_x)/2
                          else .plot_xlim[1] * 1.25
                        }
                        .signs <- regmatches(.gname, gregexpr("[+-]", 
                          .gname))[[1]]
                        .x_sign <- if (length(.signs) >= 1) 
                          .signs[1]
                        else "+"
                        .y_sign <- if (length(.signs) >= 2) 
                          .signs[2]
                        else "+"
                        # Label position: threshold-based corners for quadGate
                        # quads (finite .thr_x/.thr_y); event centroid for tmix
                        # polygon quads (NA thresholds -- the old corner fallback
                        # overflowed the panel). Always clamped inside the panel
                        # with a margin so the geom_label box is not clipped.
                        .xmargin <- 0.12 * diff(.plot_xlim)
                        .ymargin <- 0.08 * diff(.plot_ylim)
                        lx <- if (is.finite(.lbl_thr_x))
                                (if (.x_sign == "+") (.lbl_thr_x + .plot_xlim[2])/2
                                 else (.plot_xlim[1] + .lbl_thr_x)/2)
                              else stats::median(pop_e[, x_ch], na.rm = TRUE)
                        ly <- if (is.finite(.lbl_thr_y))
                                (if (.y_sign == "+") (.lbl_thr_y + .plot_ylim[2])/2 + 0.09 * diff(.plot_ylim)
                                 else (.plot_ylim[1] + .lbl_thr_y)/2 - 0.09 * diff(.plot_ylim))
                              else stats::median(pop_e[, y_ch], na.rm = TRUE)
                        if (!is.finite(lx)) lx <- mean(.plot_xlim)
                        if (!is.finite(ly)) ly <- mean(.plot_ylim)
                        lx <- max(.plot_xlim[1] + .xmargin, min(.plot_xlim[2] - .xmargin, lx))
                        ly <- max(.plot_ylim[1] + .ymargin, min(.plot_ylim[2] - .ymargin, ly))
                        lbl_df <- data.frame(lx = lx, ly = ly, 
                          stat_lbl = .mk_stat_lbl(basename(.gp), pct_val, cnt_val), label_facet = subtitle)
                        p <- p + geom_label(data = lbl_df, aes(x = lx, 
                          y = ly, label = stat_lbl), fill = "white", 
                          color = "black", alpha = 0.7, size = 2.3, 
                          fontface = "bold", linewidth = 0, inherit.aes = FALSE)
                      }
                    }
                    else {
                      for (.gp in overlay_pops) {
                        g_obj <- tryCatch(gh_pop_get_gate(gs[[sn]],
                          .gp), error = function(e) NULL)
                        if (is.null(g_obj))
                          next
                        # Skip an all-infinite rectangleGate (the gate_logic dummy
                        # gate) -- drawing it renders a misleading full-panel box.
                        # NOTE: the layout_gates+gate_logic crosshair path was
                        # reverted here while diagnosing a Positron session crash --
                        # see Pending Work / History before re-adding it.
                        if (inherits(g_obj, "rectangleGate") &&
                            all(!is.finite(c(g_obj@min, g_obj@max))))
                          next
                        if (inherits(g_obj, "rectangleGate")) {
                          # Draw if the gate constrains EITHER plotted axis (a
                          # missing axis is treated as unconstrained). Skip only
                          # when the gate touches neither -- so a listed gate whose
                          # dims don't match the panel is silently omitted.
                          .has_x <- x_ch %in% names(g_obj@min)
                          .has_y <- y_ch %in% names(g_obj@min)
                          if (!.has_x && !.has_y)
                            next
                          x_min <- if (.has_x) g_obj@min[[x_ch]] else -Inf
                          x_max <- if (.has_x) g_obj@max[[x_ch]] else  Inf
                          y_min <- if (.has_y) g_obj@min[[y_ch]] else -Inf
                          y_max <- if (.has_y) g_obj@max[[y_ch]] else  Inf
                          # 1D gate on x (one finite bound, y unconstrained) ->
                          # vline at the finite bound, whether it is min ([thr,Inf])
                          # or max ([-Inf,thr], e.g. a flipped negative gate). Same
                          # for a 1D gate on y. Only a genuinely 2D box falls
                          # through to geom_rect.
                          .x_one <- is.finite(x_min) != is.finite(x_max)
                          .y_one <- is.finite(y_min) != is.finite(y_max)
                          .x_open <- !is.finite(x_min) && !is.finite(x_max)
                          .y_open <- !is.finite(y_min) && !is.finite(y_max)
                          if (.x_one && .y_open) {
                            xv <- if (is.finite(x_min)) x_min else x_max
                            p <- p + geom_vline(data = data.frame(xint = xv,
                              label_facet = subtitle),
                              aes(xintercept = xint), color = "red",
                              linewidth = 0.7, inherit.aes = FALSE)
                          }
                          else if (.y_one && .x_open) {
                            yv <- if (is.finite(y_min)) y_min else y_max
                            p <- p + geom_hline(data = data.frame(yint = yv,
                              label_facet = subtitle),
                              aes(yintercept = yint), color = "red",
                              linewidth = 0.7, inherit.aes = FALSE)
                          }
                          else {
                            # 2D box: draw ONLY the finite edges as segments (omit
                            # open +/-Inf edges; open ends extended past the view so
                            # coord clip trims them with no perpendicular cap). A
                            # geom_rect would draw spurious top/right edges at the
                            # panel boundary (coord expansion lands them in-view).
                            .vx0 <- .plot_xlim[1] - diff(.plot_xlim); .vx1 <- .plot_xlim[2] + diff(.plot_xlim)
                            .vy0 <- .plot_ylim[1] - diff(.plot_ylim); .vy1 <- .plot_ylim[2] + diff(.plot_ylim)
                            .ys0 <- if (is.finite(y_min)) y_min else .vy0
                            .ys1 <- if (is.finite(y_max)) y_max else .vy1
                            .xs0 <- if (is.finite(x_min)) x_min else .vx0
                            .xs1 <- if (is.finite(x_max)) x_max else .vx1
                            .segs2 <- data.frame(x = numeric(0), xend = numeric(0), y = numeric(0), yend = numeric(0))
                            if (is.finite(x_min)) .segs2 <- rbind(.segs2, data.frame(x = x_min, xend = x_min, y = .ys0, yend = .ys1))
                            if (is.finite(x_max)) .segs2 <- rbind(.segs2, data.frame(x = x_max, xend = x_max, y = .ys0, yend = .ys1))
                            if (is.finite(y_min)) .segs2 <- rbind(.segs2, data.frame(x = .xs0, xend = .xs1, y = y_min, yend = y_min))
                            if (is.finite(y_max)) .segs2 <- rbind(.segs2, data.frame(x = .xs0, xend = .xs1, y = y_max, yend = y_max))
                            if (nrow(.segs2)) {
                              .segs2$label_facet <- subtitle
                              p <- p + geom_segment(data = .segs2, aes(x = x, xend = xend, y = y, yend = yend),
                                  color = "red", linewidth = 0.7, inherit.aes = FALSE)
                            }
                          }
                        }
                        else if (inherits(g_obj, "polygonGate")) {
                          bnd <- as.data.frame(g_obj@boundaries)
                          colnames(bnd) <- c("px", "py")
                          bnd$label_facet <- subtitle
                          bnd <- rbind(bnd, bnd[1, , drop = FALSE])
                          p <- p + geom_path(data = bnd, aes(x = px, 
                            y = py), color = "red", linewidth = 0.7, 
                            inherit.aes = FALSE)
                        }
                        # Stats label per overlaid gate. Normally skipped for a
                        # gate_logic row (its own centroid label is drawn after the
                        # loop), but when layout_gates explicitly lists gates we
                        # label EACH listed gate per its layout_labels.
                        if (.lg_override || pt_row$gating_method != "gate_logic") {
                          pop_ff <- tryCatch(gh_pop_get_data(gs[[sn]],
                            .gp), error = function(e) NULL)
                          parent_ff <- tryCatch(gh_pop_get_data(gs[[sn]],
                            dirname(.gp)), error = function(e) NULL)
                          if (!is.null(pop_ff) && !is.null(parent_ff) && 
                            nrow(exprs(parent_ff)) > 0) {
                            pct_val <- round(100 * nrow(exprs(pop_ff))/nrow(exprs(parent_ff)), 
                              1)
                            cnt_val <- nrow(exprs(pop_ff))
                            pop_e <- exprs(pop_ff)
                            if (x_ch %in% colnames(pop_e) && 
                              y_ch %in% colnames(pop_e)) {
                              lx <- median(pop_e[, x_ch], na.rm = TRUE)
                              ly <- median(pop_e[, y_ch], na.rm = TRUE)
                              lx <- max(.plot_xlim[1], min(.plot_xlim[2], 
                                lx))
                              ly <- max(.plot_ylim[1], min(.plot_ylim[2], 
                                ly))
                              lbl_df <- data.frame(lx = lx, ly = ly, 
                                stat_lbl = .mk_stat_lbl(basename(.gp), pct_val, cnt_val), 
                                label_facet = subtitle)
                              p <- p + geom_label(data = lbl_df, 
                                aes(x = lx, y = ly, label = stat_lbl), 
                                fill = "white", color = "black", 
                                alpha = 0.7, size = 2.3, fontface = "bold", 
                                linewidth = 0, inherit.aes = FALSE)
                            }
                          }
                        }
                      }
                      if (pt_row$gating_method == "gate_logic" &&
                        !is.null(gate_path)) {
                        .gl_pop_ff <- tryCatch(gh_pop_get_data(gs[[sn]],
                          gate_path), error = function(e) NULL)
                        .gl_parent_ff <- tryCatch(gh_pop_get_data(gs[[sn]], 
                          parent_path), error = function(e) NULL)
                        if (!is.null(.gl_pop_ff) && !is.null(.gl_parent_ff) && 
                          nrow(exprs(.gl_parent_ff)) > 0) {
                          .gl_pct <- round(100 * nrow(exprs(.gl_pop_ff))/nrow(exprs(.gl_parent_ff)), 
                            1)
                          .gl_cnt <- nrow(exprs(.gl_pop_ff))
                          # Place the label at the population's event centroid so it
                          # points at the gated region (fallback to panel center when
                          # the population is empty -- nothing to point at).
                          .gl_pe <- exprs(.gl_pop_ff)
                          .gllx <- if (nrow(.gl_pe) > 0 && x_ch %in% colnames(.gl_pe)) median(.gl_pe[, x_ch], na.rm = TRUE) else NA_real_
                          .glly <- if (nrow(.gl_pe) > 0 && y_ch %in% colnames(.gl_pe)) median(.gl_pe[, y_ch], na.rm = TRUE) else NA_real_
                          if (!is.finite(.gllx)) .gllx <- mean(.plot_xlim)
                          if (!is.finite(.glly)) .glly <- mean(.plot_ylim)
                          .gllx <- max(.plot_xlim[1] + diff(.plot_xlim)*0.10, min(.plot_xlim[2] - diff(.plot_xlim)*0.10, .gllx))
                          .glly <- max(.plot_ylim[1] + diff(.plot_ylim)*0.08, min(.plot_ylim[2] - diff(.plot_ylim)*0.08, .glly))
                          .gl_lbl_df <- data.frame(lx = .gllx,
                            ly = .glly, stat_lbl = .mk_stat_lbl(alias, .gl_pct, .gl_cnt),
                            label_facet = subtitle)
                          p <- p + geom_label(data = .gl_lbl_df, 
                            aes(x = lx, y = ly, label = stat_lbl), 
                            fill = "white", color = "black", 
                            alpha = 0.7, size = 2.5, fontface = "bold", 
                            linewidth = 0, inherit.aes = FALSE)
                        }
                      }
                    }
                    panels[[panel_ids[i]]] <- p
                  }
                  else if (plot_type == "histogram") {
                    x_ch <- channels[1]
                    fr <- gh_pop_get_data(gs[[sn]], parent_path)
                    vals <- exprs(fr)[, x_ch]
                    x_lo <- quantile(vals, 0.005, names = FALSE)
                    x_hi <- quantile(vals, 0.995, names = FALSE)
                    vline_x <- NULL
                    if (!is.null(gate_path)) {
                      g_obj <- tryCatch(gh_pop_get_gate(gs[[sn]], 
                        gate_path), error = function(e) NULL)
                      if (!is.null(g_obj) && inherits(g_obj, 
                        "rectangleGate")) 
                        vline_x <- g_obj@min[[x_ch]]
                    }
                    p <- ggplot(data.frame(x = vals), aes(x = x)) + 
                      geom_histogram(aes(y = after_stat(density)), 
                        bins = 100, fill = "steelblue", color = NA, 
                        alpha = 0.5) + geom_density(color = "firebrick", 
                      linewidth = 0.8) + {
                      if (!is.null(vline_x)) 
                        geom_vline(xintercept = vline_x, color = "firebrick", 
                          linewidth = 0.8, linetype = "dashed")
                    } + xlim(x_lo, x_hi) + labs(x = x_ch, y = "Density", 
                      subtitle = subtitle) + .fasta_plot_theme
                    panels[[panel_ids[i]]] <- p
                  }
                  else {
                    warning("[apply_plotting] Unknown PlotType '", 
                      plot_type, "' for '", alias, "'.")
                    panels[[panel_ids[i]]] <- patchwork::plot_spacer()
                  }
                }
            }
            panels <- lapply(panels, function(px) if (is.null(px) || 
                !inherits(px, "gg")) 
                patchwork::plot_spacer()
            else px)
            layout_list[[paste0("layout_", lid)]] <- if (!is.null(design_str)) 
                patchwork::wrap_plots(panels, design = design_str)
            else patchwork::wrap_plots(panels)
        }
        plot_list[[sn]] <- layout_list
    }
    plot_list
}


#' run_gating_qc — comprehensive per-gate, per-sample boundary QC
#'
#' Checks:
#'   1. Identical thresholds across ALL samples (per-sample gate not working)
#'   2. Identical thresholds across groups (groupBy not differentiating)
#'   3. copy_gate boundary mismatch vs source gate
#'   4. Empty populations (0 events)
#'   5. Extreme percentages (>99.5% or <0.1% of parent)
#'   6. Low parent event count (<50 events -> unreliable %)
#'   7. Gate boundary outside data range (threshold outside P1-P99 of channel)
#'   8. gate_logic count vs manual OR/AND of components
#'
#' @param gs         GatingSet
#' @param gt_df      Gating template data.frame
#' @param min_events Threshold for "low parent events" check (default 50)
#'
#' @return data.frame with columns:
#'   alias, parent, gating_method, check, status, group, sample,
#'   threshold_x, threshold_y, pct, n_events, n_parent, detail

# ── Gating QC ─────────────────────────────────────────────────────────────────

#' run_gating_qc — comprehensive gating template compliance and boundary QC
#'
#' Failure types (check column):
#'  Template compliance:
#'   population_missing             alias in template but not in GatingSet
#'   active_gate_false_present      gate marked active_gate=FALSE but still in gs
#'   parent_path_mismatch           gate exists under wrong parent in gs
#'   gate_logic_component_missing   component pop for gate_logic not in gs
#'  Threshold / boundary:
#'   identical_threshold_all_samples   all samples share the same threshold
#'   identical_threshold_across_groups groups share the same threshold (groupBy)
#'   within_group_threshold_variation  samples in same group have different thresholds
#'   threshold_outside_data_range      threshold is outside P1-P99 of channel data
#'  Counts / percentages:
#'   quadrant_children_sum_mismatch   4 quad children do not sum to parent
#'   gate_logic_count_mismatch        gate_logic count differs from manual OR/AND
#'   empty_population                 0 events in population
#'   extreme_pct_high                 >99.5% of parent
#'   extreme_pct_low                  <0.1% of parent (with >100 parent events)
#'   low_parent_events                parent has fewer than min_events
#'
#' @param gs          GatingSet
#' @param gt_df       Gating template data.frame
#' @param min_events  Threshold for low_parent_events check (default 50)
run_gating_qc <- function(gs, gt_df, min_events = 50) {

  all_pops <- gs_get_pop_paths(gs, path = "full")
  all_sns  <- sampleNames(gs)
  pd       <- pData(gs)
  all_ch   <- colnames(gs[[1]])
  all_mk   <- markernames(gs[[1]])

  .rp <- function(alias) {
    h <- all_pops[basename(all_pops) == basename(trimws(alias))]
    if (length(h) == 1) return(h[1])
    if (length(h) > 1) { m <- h[endsWith(h, trimws(alias))]; if (length(m) >= 1) return(m[1]) }
    NA_character_
  }

  .thresh <- function(g) {
    if (is.null(g)) return(c(NA_real_, NA_real_))
    qi <- attr(g, "quadintersection")
    if (!is.null(qi)) return(c(qi[1], if (length(qi) > 1) qi[2] else NA_real_))
    if (inherits(g, "rectangleGate")) {
      mn <- g@min[is.finite(g@min)]; mx <- g@max[is.finite(g@max)]
      v  <- c(mn, mx)
      return(c(if (length(v) > 0) v[1] else NA_real_, if (length(v) > 1) v[2] else NA_real_))
    }
    c(NA_real_, NA_real_)
  }

  .rch <- function(d) {
    d <- trimws(d)
    if (d %in% all_ch) d else { n <- names(all_mk)[match(d, all_mk)]; if (!is.na(n)) n else d }
  }

  .is_ec <- function(args_str) {
    nzchar(args_str) && (grepl("external_control", args_str, ignore.case=TRUE) ||
                          grepl("control_col",     args_str, ignore.case=TRUE))
  }

  rows <- list()
  .row <- function(alias, parent, method, check, status, group=NA, sample=NA,
                   tx=NA, ty=NA, pct=NA, n=NA, npar=NA, detail="")
    data.frame(alias=alias, parent=parent, gating_method=method, check=check,
               status=status, group=as.character(group), sample=as.character(sample),
               threshold_x=round(tx,2), threshold_y=round(ty,2),
               pct=round(pct,2), n_events=as.integer(n), n_parent=as.integer(npar),
               detail=detail, stringsAsFactors=FALSE)

  # ── Check A: active_gate=FALSE populations still in gs ──────────────────────
  if ("active_gate" %in% names(gt_df)) {
    inactive <- gt_df[!is.na(gt_df$active_gate) &
                      tolower(as.character(gt_df$active_gate)) == "false", ]
    for (k in seq_len(nrow(inactive))) {
      al <- trimws(inactive$alias[k])
      if (!is.na(.rp(al)))
        rows[[length(rows)+1]] <- .row(al, inactive$parent[k], inactive$gating_method[k],
          "active_gate_false_present", "FLAG",
          detail = "Gate has active_gate=FALSE but population exists in GatingSet")
    }
  }

  for (i in seq_len(nrow(gt_df))) {
    row    <- gt_df[i, , drop=FALSE]
    alias  <- trimws(row$alias)
    parent <- trimws(row$parent)
    method <- trimws(row$gating_method)
    args_str <- if (!is.na(row$gating_args)) trimws(row$gating_args) else ""
    grpby  <- if ("groupBy" %in% names(row) && !is.na(row$groupBy) &&
                  nzchar(trimws(row$groupBy))) trimws(row$groupBy) else NULL
    active <- if ("active_gate" %in% names(row) && !is.na(row$active_gate))
                tolower(as.character(row$active_gate)) else "true"

    if (active == "false") next   # skip deliberately inactive gates

    pop_path <- .rp(alias)

    # ── Check B: population_missing ─────────────────────────────────────────
    is_quad_method <- method %in% .FASTA_QUAD_2D_METHODS
    is_copy_method <- method == "gate_logic_copy" || method == "copy_gate"

    if (is_quad_method && is.na(pop_path)) {
      # Reconstruct the 4 expected quad population names from name_x / name_y args
      .qa <- tryCatch(eval(parse(text=paste0("list(",args_str,")"))), error=function(e) list())
      .dims_raw <- trimws(strsplit(trimws(row$dims), ",")[[1]])
      .nx <- if (!is.null(.qa$name_x)) .qa$name_x else .rch(.dims_raw[1])
      .ny <- if (!is.null(.qa$name_y)) .qa$name_y else
               if (length(.dims_raw) > 1) .rch(.dims_raw[2]) else ""
      .qn <- c(paste0(.nx,"-",.ny,"+"), paste0(.nx,"+",.ny,"+"),
               paste0(.nx,"+",.ny,"-"), paste0(.nx,"-",.ny,"-"))
      .missing_qn <- .qn[!.qn %in% basename(all_pops)]
      if (length(.missing_qn) > 0)
        rows[[length(rows)+1]] <- .row(alias, parent, method, "population_missing", "FLAG",
          detail=paste0("Expected quad populations not in GatingSet: ",
                        paste(.missing_qn, collapse=", ")))
      # Set pop_path to the first existing quad child (not its parent)
      # so subsequent checks (thresholds, % etc.) read the quad crosshair, not the parent gate
      .found_qp <- all_pops[basename(all_pops) %in% .qn]
      pop_path <- if (length(.found_qp) > 0) .found_qp[1] else NA_character_
      if (is.na(pop_path)) next

    } else if (is_copy_method && is.na(pop_path)) {
      # Resolve children from copy_from source gate
      .ca <- tryCatch(eval(parse(text=paste0("list(",args_str,")"))), error=function(e) list())
      .src_alias <- .ca$copy_from
      if (!is.null(.src_alias)) {
        .src_path <- .rp(.src_alias)
        if (!is.na(.src_path)) {
          .src_children <- basename(all_pops[dirname(all_pops) == .src_path])
          .par_path <- .rp(parent)
          if (!is.na(.par_path) && length(.src_children) > 0) {
            .missing_cc <- .src_children[
              !paste0(.par_path, "/", .src_children) %in% all_pops]
            if (length(.missing_cc) > 0)
              rows[[length(rows)+1]] <- .row(alias, parent, method, "population_missing", "FLAG",
                detail=paste0("copy_gate children not found under parent '", parent, "': ",
                              paste(.missing_cc, collapse=", ")))
          }
        }
      }
      next  # no single pop_path for copy_gate aliases

    } else if (is.na(pop_path)) {
      rows[[length(rows)+1]] <- .row(alias, parent, method, "population_missing", "FLAG",
        detail="Alias in gating template but not found in GatingSet")
      next
    }

    # ── Check C: parent_path_mismatch ───────────────────────────────────────
    # Skip for quad-method gates: pop_path was set to the quad parent, not a child
    # Skip for root/notDebris: "root" resolves to "/" in dirname — false positive
    if (!is_quad_method && parent != "root") {
      actual_parent_path <- dirname(pop_path)
      expected_par_path  <- .rp(parent)
      if (!is.na(expected_par_path) &&
          !identical(actual_parent_path, expected_par_path) &&
          basename(actual_parent_path) != basename(parent)) {
        rows[[length(rows)+1]] <- .row(alias, parent, method,
          "parent_path_mismatch", "WARN",
          detail = paste0("Expected parent: ", expected_par_path,
                          " | Actual parent: ", actual_parent_path))
      }
    }

    # ── Per-sample stats ─────────────────────────────────────────────────────
    sn_stats <- lapply(all_sns, function(sn) {
      g     <- tryCatch(gh_pop_get_gate(gs[[sn]], pop_path), error = function(e) NULL)
      thr   <- .thresh(g)
      n_pop <- tryCatch(gh_pop_get_count(gs[[sn]], pop_path), error = function(e) NA_integer_)
      par_p <- .rp(parent)
      n_par <- if (!is.na(par_p))
                 tryCatch(gh_pop_get_count(gs[[sn]], par_p), error = function(e) NA_integer_)
               else NA_integer_
      pct   <- if (!is.na(n_pop) && !is.na(n_par) && n_par > 0) 100*n_pop/n_par else NA_real_
      grp   <- if (!is.null(grpby) && grpby %in% colnames(pd)) pd[sn, grpby] else NA
      list(sn=sn, tx=thr[1], ty=thr[2], n=n_pop, npar=n_par, pct=pct, grp=grp)
    })

    tx_vec <- sapply(sn_stats, `[[`, "tx")
    ty_vec <- sapply(sn_stats, `[[`, "ty")

    # collapseDataForGating intentionally produces ONE gate replicated across the
    # group (or all samples) -- identical thresholds are BY DESIGN, so skip the
    # identical-threshold checks (1 & 2) for these rows.
    .is_collapse_row <- "collapseDataForGating" %in% names(row) &&
                        !is.na(row$collapseDataForGating) &&
                        isTRUE(suppressWarnings(as.logical(row$collapseDataForGating)))

    # ── Check 1: identical_threshold_all_samples ─────────────────────────────
    if (!.is_collapse_row &&
        method %in% c("gate_FMO","gate_FMO_2d","gate_flowmeans","gate_flowmeans_2d",
                      "gate_tmix_2d",
                      "gate_mindensity","gate_mindensity_2d","gate_singlet")) {
      if (!all(is.na(tx_vec)) && length(unique(round(tx_vec[!is.na(tx_vec)], 1))) == 1)
        rows[[length(rows)+1]] <- .row(alias, parent, method,
          "identical_threshold_all_samples", "WARN",
          detail = paste0("All samples: threshold_x=", round(tx_vec[1],1),
                          if (!is.na(ty_vec[1])) paste0(", threshold_y=",round(ty_vec[1],1)) else ""))
    }

    # ── Check 2: identical_threshold_across_groups ──────────────────────────
    if (!.is_collapse_row && !is.null(grpby) && method != "gate_static") {
      grp_vec     <- sapply(sn_stats, `[[`, "grp")
      grps        <- unique(grp_vec[!is.na(grp_vec)])
      if (length(grps) > 1) {
        grp_meds <- tapply(tx_vec, grp_vec, function(x) median(x, na.rm=TRUE))
        grp_meds <- grp_meds[!is.na(grp_meds)]
        if (length(unique(round(grp_meds, 1))) < length(grp_meds)) {
          dup_grps <- names(grp_meds)[duplicated(round(grp_meds, 1))]
          rows[[length(rows)+1]] <- .row(alias, parent, method,
            "identical_threshold_across_groups", "FLAG",
            detail = paste0("Groups share same threshold: ", paste(dup_grps, collapse=", ")))
        }
      }
    }

    # ── Check 3: within_group_threshold_variation ───────────────────────────
    if (!is.null(grpby) && .is_ec(args_str) && method != "gate_static") {
      grp_vec <- sapply(sn_stats, `[[`, "grp")
      for (grp in unique(grp_vec[!is.na(grp_vec)])) {
        grp_tx <- tx_vec[grp_vec == grp & !is.na(grp_vec)]
        grp_tx <- grp_tx[!is.na(grp_tx)]
        if (length(grp_tx) > 1 && length(unique(round(grp_tx, 1))) > 1)
          rows[[length(rows)+1]] <- .row(alias, parent, method,
            "within_group_threshold_variation", "FLAG", group=grp,
            detail = paste0("Group '", grp, "' thresholds not uniform: [",
                            paste(sort(unique(round(grp_tx,1))), collapse=", "), "]"))
      }
    }

    # ── Check 4: quadrant_children_sum_mismatch ──────────────────────────────
    if (method %in% c(.FASTA_QUAD_2D_METHODS,"copy_gate")) {
      quad_children <- all_pops[dirname(all_pops) == pop_path]
      if (length(quad_children) == 4) {
        for (st in sn_stats) {
          n_parent <- st$npar
          if (is.na(n_parent) || n_parent == 0) next
          # re-get parent count using pop_path as parent
          n_self <- tryCatch(gh_pop_get_count(gs[[st$sn]], pop_path), error=function(e) NA_integer_)
          n_kids <- sum(sapply(quad_children, function(cp)
            tryCatch(gh_pop_get_count(gs[[st$sn]], cp), error=function(e) 0L)))
          if (!is.na(n_self) && n_self > 10 && abs(n_kids - n_self) > max(2, 0.01*n_self))
            rows[[length(rows)+1]] <- .row(alias, parent, method,
              "quadrant_children_sum_mismatch", "FLAG", sample=st$sn,
              n=n_self, detail=paste0("Parent=",n_self," | Sum(4 children)=",n_kids,
                                      " | diff=",abs(n_kids-n_self)))
        }
      }
    }

    # ── Checks 5-7: per-sample event counts / percentages ────────────────────
    for (st in sn_stats) {
      grp_lbl <- if (!is.null(grpby)) st$grp else NA

      if (!is.na(st$n) && st$n == 0)
        rows[[length(rows)+1]] <- .row(alias, parent, method, "empty_population", "FLAG",
          group=grp_lbl, sample=st$sn, n=st$n, npar=st$npar, detail="0 events")

      if (!is.na(st$npar) && st$npar > 0 && st$npar < min_events)
        rows[[length(rows)+1]] <- .row(alias, parent, method, "low_parent_events", "WARN",
          group=grp_lbl, sample=st$sn, n=st$n, npar=st$npar,
          detail=paste0("Parent has only ",st$npar," events"))

      # extreme_pct_high: population is essentially all of its parent (gate not discriminating)
      if (!is.na(st$pct) && !is.na(st$npar) && st$npar > 0 && st$pct > 99.5)
        rows[[length(rows)+1]] <- .row(alias, parent, method, "extreme_pct_high", "WARN",
          group=grp_lbl, sample=st$sn, pct=st$pct, n=st$n, npar=st$npar,
          detail=paste0(round(st$pct,2),"% of parent (>99.5%)"))

      # extreme_pct_low: tiny but non-zero population, with enough parent events to be meaningful
      if (!is.na(st$pct) && !is.na(st$n) && st$n > 0 &&
          !is.na(st$npar) && st$npar > 100 && st$pct < 0.1)
        rows[[length(rows)+1]] <- .row(alias, parent, method, "extreme_pct_low", "WARN",
          group=grp_lbl, sample=st$sn, pct=st$pct, n=st$n, npar=st$npar,
          detail=paste0(round(st$pct,3),"% of parent (<0.1%, ",st$npar," parent events)"))
    }

    # ── Check 8: threshold_outside_data_range ────────────────────────────────
    if (!is.na(tx_vec[1]) && !all(is.na(tx_vec)) && method != "gate_static") {
      dims_raw <- trimws(strsplit(trimws(row$dims), ",")[[1]])
      ch_x <- tryCatch(.rch(dims_raw[1]), error = function(e) NA_character_)
      if (!is.na(ch_x) && ch_x %in% all_ch) {
        # Parse control column/value for FMO source identification
        .ec_args   <- tryCatch(eval(parse(text=paste0("list(",args_str,")"))), error=function(e) list())
        .ctrl_cols <- unique(c(.ec_args$control_col, .ec_args$control_col_x, .ec_args$control_col_y))
        .ctrl_cols <- .ctrl_cols[!sapply(.ctrl_cols, is.null)]
        .ctrl_vals <- unique(c(.ec_args$control_val, .ec_args$control_val_x, .ec_args$control_val_y))
        .ctrl_vals <- .ctrl_vals[!sapply(.ctrl_vals, is.null)]

        # Option A (2026-06-13): per_sample gates only. A control-derived
        # (external/FMO) threshold is placed at the CONTROL's background level, so
        # testing it against each sample's own P1-P99.95 is apples-to-oranges and
        # floods the table with expected, non-actionable flags. Skip EC rows; the
        # check is only meaningful when threshold + data come from the same sample.
        .is_ec_row <- length(.ctrl_cols) > 0 || grepl("external_control", args_str)
        if (!.is_ec_row)
        for (st in sn_stats) {
          if (is.na(st$tx)) next

          # Suppress flag if this sample IS the FMO control used to derive the gate
          .is_fmo_source <- length(.ctrl_cols) > 0 && length(.ctrl_vals) > 0 &&
            any(sapply(seq_along(.ctrl_cols), function(j) {
              col <- .ctrl_cols[[j]]; val <- .ctrl_vals[[min(j, length(.ctrl_vals))]]
              col %in% colnames(pd) && grepl(val, pd[st$sn, col], ignore.case=TRUE)
            }))
          if (.is_fmo_source) next

          par_p <- .rp(parent); if (is.na(par_p)) next
          fr <- tryCatch(gh_pop_get_data(gs[[st$sn]], par_p), error=function(e) NULL)
          if (is.null(fr) || nrow(exprs(fr)) < 10) next
          qs <- quantile(exprs(fr)[, ch_x], c(0.01, 0.9995), na.rm=TRUE)
          if (st$tx < qs[1] || st$tx > qs[2])
            rows[[length(rows)+1]] <- .row(alias, parent, method,
              "threshold_outside_data_range", "WARN", sample=st$sn, tx=st$tx,
              detail=paste0("Threshold ",round(st$tx,1)," outside P1-P99.95 [",
                            round(qs[1],1),"-",round(qs[2],1),"] of ",ch_x))
        }
      }
    }

    # ── Check 9: gate_logic ──────────────────────────────────────────────────
    if (method == "gate_logic") {
      lg_args <- tryCatch(eval(parse(text=paste0("list(",args_str,")"))), error=function(e) list())
      # gates= may be a comma-string ('A,B,C') or c('A','B','C') -- parse the same
      # way apply_gating2 does, else a comma-string is treated as ONE bogus pop
      # name and every multi-component gate_logic falsely flags as "missing".
      lg_gates <- .fasta_parse_pop_list(lg_args$gates)
      if (length(lg_gates) > 0 && !is.null(lg_args$logic)) {
        comp_paths <- sapply(lg_gates, .rp)

        # Check D: missing components
        missing_comp <- lg_gates[is.na(comp_paths)]
        if (length(missing_comp) > 0)
          rows[[length(rows)+1]] <- .row(alias, parent, method,
            "gate_logic_component_missing", "FLAG",
            detail=paste0("Component(s) not in GatingSet: ", paste(missing_comp, collapse=", ")))

        valid_comp <- comp_paths[!is.na(comp_paths)]
        if (length(valid_comp) == length(lg_gates)) {
          for (st in sn_stats) {
            comp_idx <- lapply(valid_comp, function(pp)
              tryCatch(gh_pop_get_indices(gs[[st$sn]], pp), error=function(e) NULL))
            comp_idx <- Filter(Negate(is.null), comp_idx)
            if (length(comp_idx) == length(valid_comp)) {
              .lg <- toupper(lg_args$logic)
              # NOT = parent_idx & !component (matches apply_gating2's index logic);
              # comparing a NOT count against the bare component count flags EVERY
              # NOT gate as a mismatch -- handle it explicitly.
              expected <- if (.lg == "OR") {
                sum(Reduce("|", comp_idx))
              } else if (.lg == "NOT") {
                .par_p <- .rp(parent)
                .par_idx <- if (!is.na(.par_p))
                    tryCatch(gh_pop_get_indices(gs[[st$sn]], .par_p), error=function(e) NULL)
                  else rep(TRUE, length(comp_idx[[1]]))
                if (is.null(.par_idx)) NA_real_ else sum(.par_idx & !comp_idx[[1]])
              } else {
                sum(Reduce("&", comp_idx))
              }
              if (!is.na(expected) && !is.na(st$n) && abs(st$n - expected) > 2)
                rows[[length(rows)+1]] <- .row(alias, parent, method,
                  "gate_logic_count_mismatch", "FLAG", sample=st$sn, n=st$n,
                  detail=paste0("Stored=",st$n," vs manual ",lg_args$logic,"=",expected))
            }
          }
        }
      }
    }
  }

  if (length(rows) == 0) {
    message("[run_gating_qc] All checks passed. No issues found.")
    return(invisible(data.frame()))
  }

  result <- do.call(rbind, rows)
  result$status_rank <- c(FLAG=1L, WARN=2L, INFO=3L)[result$status]
  result <- result[order(result$status_rank, result$check, result$alias), ]
  result$status_rank <- NULL
  rownames(result) <- NULL

  n_flag <- sum(result$status == "FLAG", na.rm=TRUE)
  n_warn <- sum(result$status == "WARN", na.rm=TRUE)
  message("[run_gating_qc] Done: ", n_flag, " FLAG(s), ", n_warn,
          " WARN(s) across ", length(unique(result$alias)), " gate(s).")
  result
}


#' gating_qc_tables — split QC results into one table per check type
#'
#' Returns a named list of data.frames, one per unique check type.
#' Each table is sorted FLAG first, then WARN, then by alias/sample.
#' Use for report sections: one table = one type of issue.
#'
#' @param qc_result  Output of run_gating_qc()
#' @param checks     Optional character vector to subset specific check types
#'
#' @return Named list of data.frames
gating_qc_tables <- function(qc_result, checks = NULL) {

  # All known check types with human-readable labels
  check_labels <- c(
    population_missing                = 'Missing Populations',
    active_gate_false_present         = 'Disabled Gates Still Present',
    parent_path_mismatch              = 'Parent Path Mismatch',
    gate_logic_component_missing      = 'gate_logic: Missing Components',
    identical_threshold_all_samples   = 'Identical Threshold (All Samples)',
    identical_threshold_across_groups = 'Identical Threshold Across Groups',
    within_group_threshold_variation  = 'Within-Group Threshold Variation',
    threshold_outside_data_range      = 'Threshold Outside Data Range',
    quadrant_children_sum_mismatch    = 'Quadrant Children Sum Mismatch',
    gate_logic_count_mismatch         = 'gate_logic Count Mismatch',
    empty_population                  = 'Empty Populations',
    extreme_pct_high                  = 'Extreme % High (>99.5%)',
    extreme_pct_low                   = 'Extreme % Low (<0.1%)',
    low_parent_events                 = 'Low Parent Event Count'
  )

  all_checks <- if (!is.null(checks)) checks else names(check_labels)

  # Empty-row template for 'no issues' entries
  .empty_cols <- c('alias','parent','gating_method','check','status','group',
                   'sample','threshold_x','threshold_y','pct','n_events','n_parent','detail')
  .no_issue_row <- function(ct)
    setNames(as.data.frame(matrix(NA, nrow=1, ncol=length(.empty_cols)),
                           stringsAsFactors=FALSE), .empty_cols) |>
    within({ check = ct; status = 'PASS';
             detail = paste0('No issues found: ', check_labels[ct]) })

  status_rank <- c(FLAG=1L, WARN=2L, PASS=3L, INFO=4L)

  tables <- setNames(lapply(all_checks, function(ct) {
    # Filter to this check type if results exist
    has_data <- !is.null(qc_result) && nrow(qc_result) > 0 && ct %in% qc_result$check
    keep_cols <- c("status","detail","alias","parent","gating_method","check","sample")

    if (!has_data) {
      df <- .no_issue_row(ct)
    } else {
      df <- qc_result[qc_result$check == ct, , drop=FALSE]
      df$status_rank <- status_rank[df$status]
      df <- df[order(df$status_rank, df$alias, df$sample), ]
      df$status_rank <- NULL
      rownames(df) <- NULL
    }
    df[, intersect(keep_cols, names(df)), drop=FALSE]
  }), all_checks)

  # Summary
  message('[gating_qc_tables] ', length(tables), ' check type(s):')
  for (ct in all_checks) {
    lbl   <- check_labels[ct]
    n_f   <- sum(tables[[ct]]$status == 'FLAG',  na.rm=TRUE)
    n_w   <- sum(tables[[ct]]$status == 'WARN',  na.rm=TRUE)
    n_p   <- sum(tables[[ct]]$status == 'PASS',  na.rm=TRUE)
    status_str <- if (n_f > 0 || n_w > 0)
      paste0(n_f, ' FLAG, ', n_w, ' WARN') else 'PASS'
    message('  [', ct, '] ', lbl, ' — ', status_str)
  }
  tables
}

# panel_QC — compare GatingSet markernames/channels against master panel inventory
# Uses Marker_Alias column as fallback when Antigen name differs from GS marker name.
# Checks (1) panel markers present in GatingSet, (2) GatingSet markers present in panel,
# (3) panel fluorochromes matched to a GatingSet channel, (4) GatingSet channels matched
# to a panel fluorochrome. Viability markers/channels are exempt.
# Returns a data.frame of mismatches (empty = all clear).
panel_QC <- function(gs, PanelMasterInventoryLocation, PanelIteration) {
  if (!requireNamespace('readxl', quietly=TRUE))
    stop("[panel_QC] Package 'readxl' required.")

  inv   <- readxl::read_excel(PanelMasterInventoryLocation, sheet = 1)
  panel <- inv[!is.na(inv$Iteration) & inv$Iteration == PanelIteration, , drop=FALSE]
  if (nrow(panel) == 0)
    stop("[panel_QC] No rows for Iteration='", PanelIteration,
         "'. Available:\n  ", paste(sort(unique(na.omit(inv$Iteration))), collapse='\n  '))

  # Strip zero-width / non-breaking invisible characters (the inventory has them)
  .clean <- function(x) trimws(gsub('[\u200B\u200C\u200D\u2060\u00A0\uFEFF]', '', x, perl=TRUE))
  .norm_fluor <- function(x) {
    x <- .clean(x)
    x <- tolower(gsub('[[:space:]]+', '', x))      # drop spaces
    x <- sub('-[ah]$', '', x)                       # drop -A (area) / -H (height) channel suffix FIRST
    x <- gsub('[-_.‐-―]', '', x)          # drop rogue separators: - _ . and dash variants
    x <- sub('^af([0-9])', 'alexafluor\\1', x)      # AF### -> alexafluor###
    x
  }
  # Normalise an antigen NAME for matching: strip separators (_ - . whitespace)
  # so PD-1 == PD1, HLA-DR == HLADR, CD14_CD16 == CD14CD16. Does NOT strip '/'
  # or parentheses -- TCF7/1, CD197(CCR7) etc. still rely on Marker_Alias.
  .norm_ag <- function(x) tolower(gsub('[_[:space:].-]', '', .clean(x)))

  gs_mk       <- markernames(gs); gs_mk <- gs_mk[nzchar(gs_mk)]
  gs_antigens <- unique(.clean(gs_mk))

  pAg   <- .clean(panel$Antigen)
  pAli  <- if ('Marker_Alias' %in% names(panel)) .clean(panel$Marker_Alias) else rep(NA_character_, nrow(panel))
  pFl   <- .clean(panel$Fluorochrome)
  keep  <- !grepl('^unconjugated$', pFl, ignore.case=TRUE) & nzchar(pFl)
  pAg   <- pAg[keep]; pAli <- pAli[keep]; pFl_raw <- pFl[keep]

  # Effective name: Marker_Alias takes precedence when not NA
  pEff  <- ifelse(!is.na(pAli) & nzchar(pAli), pAli, pAg)

  # All names recognised as valid panel markers (both Antigen + Alias)
  pAll_names <- unique(.norm_ag(c(pAg, pAli[!is.na(pAli) & nzchar(pAli)])))

  rows <- list()

  # Viability dyes are named differently in panel vs GatingSet (e.g. "Live Dead
  # eFluor 780" vs marker "Viability") -- exempt them in BOTH directions.
  .is_viability <- function(...) grepl('viabilit|live.?dead|aqua|zombie|fixable|efluor.?780|7.?aad|dapi',
                                       paste(..., collapse=' '), ignore.case=TRUE)

  # 1. Panel markers (by effective name) not found in GatingSet
  for (k in seq_along(pEff)) {
    if (!.norm_ag(pEff[k]) %in% .norm_ag(gs_antigens) && !.is_viability(pEff[k], pFl_raw[k]))
      rows[[length(rows)+1]] <- data.frame(
        source       = 'Panel Inventory',
        flag         = 'Not found in GatingSet markernames',
        antigen      = pAg[k],
        marker_alias = if (!is.na(pAli[k]) && nzchar(pAli[k])) pAli[k] else NA_character_,
        fluorochrome = pFl_raw[k],
        stringsAsFactors = FALSE)
  }

  # 2. GatingSet markers not found in panel (neither Antigen nor Alias); viability exempt
  for (m in gs_antigens) {
    if (!.norm_ag(m) %in% pAll_names &&
        !grepl('viabilit|live.dead|aqua|zombie|fixable', m, ignore.case=TRUE))
      rows[[length(rows)+1]] <- data.frame(
        source       = 'GatingSet markernames',
        flag         = 'Not found in Panel Inventory',
        antigen      = m,
        marker_alias = NA_character_,
        fluorochrome = names(gs_mk)[.clean(gs_mk) == m][1],
        stringsAsFactors = FALSE)
  }

  # Fluorochrome / channel matching (normalised: strips spaces, -A suffix, AF->AlexaFluor)
  gs_ch_raw      <- names(gs_mk)
  gs_fluors_norm <- .norm_fluor(gs_ch_raw)
  pFl_norm       <- .norm_fluor(pFl_raw)

  # 3. Panel fluorochromes not matched to any GatingSet channel (viability exempt)
  for (k in seq_along(pFl_norm)) {
    if (!pFl_norm[k] %in% gs_fluors_norm && !.is_viability(pEff[k], pFl_raw[k]))
      rows[[length(rows)+1]] <- data.frame(
        source       = 'Panel Inventory',
        flag         = 'Fluorochrome not matched to any GatingSet channel',
        antigen      = pAg[k],
        marker_alias = if (!is.na(pAli[k]) && nzchar(pAli[k])) pAli[k] else NA_character_,
        fluorochrome = pFl_raw[k],
        stringsAsFactors = FALSE)
  }

  # 4. GatingSet channels not matched to any panel fluorochrome (viability channels exempt)
  for (j in seq_along(gs_ch_raw)) {
    if (!gs_fluors_norm[j] %in% pFl_norm &&
        !grepl('efluor|viabilit|live.dead|aqua|zombie|fixable|780',
               gs_ch_raw[j], ignore.case=TRUE))
      rows[[length(rows)+1]] <- data.frame(
        source       = 'GatingSet markernames',
        flag         = 'Channel not matched to any panel fluorochrome',
        antigen      = .clean(gs_mk)[j],
        marker_alias = NA_character_,
        fluorochrome = gs_ch_raw[j],
        stringsAsFactors = FALSE)
  }

  if (length(rows) == 0) {
    message('[panel_QC] Panel as expected. Iteration: "', PanelIteration,
            '" | ', length(pEff), ' marker(s) verified.')
    return(invisible(data.frame()))
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  message('[panel_QC] ', nrow(result), ' mismatch(es) found for Iteration: "',
          PanelIteration, '"')
  result
}


# =============================================================================
# SECTION 9 (v0.70 PROTOTYPE): apply_gating2 -- gt_gating-free sequential engine
# -----------------------------------------------------------------------------
# Parallel re-implementation of apply_gating that DROPS openCyto's gt_gating /
# gatingTemplate orchestration. It walks the template top-to-bottom; because
# parents are always created before children, there is NO DAG resolution, NO
# precompute cache, NO deferral, and NO .logic_indices index restoration. Each
# row computes a per-sample (or per-group / control-derived) gate by calling the
# SAME .gate_* plugin bodies (in per_sample mode) used by apply_gating, then adds
# it with flowWorkspace's native gs_pop_add / booleanFilter.
#
# This is a PARALLEL prototype -- apply_gating, the plugins, and everything else
# are untouched. Run both and diff gs_pop_get_stats() to confirm parity before
# retiring apply_gating.
#
# PROTOTYPE SCOPE (parity target = current LYL273 template):
#   methods : gate_mindensity(_2d), gate_flowmeans(_2d), gate_FMO(_2d),
#             gate_tmix_2d, gate_singlet, gate_static, gate_logic, copy_gate,
#             + builtins singletGate / gate_flowclust_2d
#   control : single control_col/control_val (optionally + groupBy)
#   pop="-" : 1D rectangle gates (flipped to [-Inf, thr])
# Anything outside this stops() loudly. Deliberately deferred until basics pass:
#   per-axis control_col_x/y, 2-sided or polygon negatives, quad-children copy.
#
# !! FIRST-RUN VERIFY (things I could not test statically):
#   (a) flowCore::char2booleanFilter() exists / signature -- used for gate_logic
#       and is wrapped in a tryCatch fallback to the booleanFilter() NSE form.
#   (b) gs_pop_add(gs, <named-per-sample list>, parent=) truly applies per sample
#       for rectangleGate/polygonGate (your old workflow says yes; quadGates are
#       sidestepped here by decomposing to 4 rectangleGates, avoiding bug #4).
# =============================================================================

# Pure gate calculator: feed ONE flowFrame to the right plugin in per_sample mode
# and return its gate object. No cache, no control logic -- the driver picks fr.
.fasta_compute_gate <- function(fr, method, channels, filterId, args) {
  drop <- c("method", "control_col", "control_val", "control_idx",
            "control_col_x", "control_col_y", "control_val_x", "control_val_y",
            "control_idx_x", "control_idx_y",
            "clusterFrom", "clusterFrom_x", "clusterFrom_y", "groupBy")
  plugin_args <- args[!names(args) %in% drop]
  plugin_args$method <- "per_sample"
  fn <- tryCatch(get(paste0(".", method), envir = globalenv()), error = function(e) NULL)
  if (!is.null(fn)) {
    # Reject typo'd / unrecognized gating_args LOUDLY -- otherwise they vanish into
    # the plugin's `...` and silently produce a wrong gate (e.g. `x_min` instead of
    # `min_x`, or the earlier `y_min, 5000`). gate_static is exempt: it reads ALL
    # its bounds (xmin/xmax/ymin/ymax/min/max) from `...` by design.
    if (method != "gate_static") {
      .valid   <- names(formals(fn))
      .unknown <- setdiff(names(plugin_args), .valid)
      if (length(.unknown))
        stop(sprintf(
          "[apply_gating2] unrecognized gating_arg(s) for method '%s': %s\n  Valid args: %s",
          method, paste(.unknown, collapse = ", "),
          paste(setdiff(.valid, c("fr", "pp_res", "channels", "filterId", "...")),
                collapse = ", ")), call. = FALSE)
    }
    return(do.call(fn, c(list(fr = fr, pp_res = NULL, channels = channels,
                              filterId = filterId), plugin_args)))
  }
  # Builtins without a FASTA wrapper
  if (method == "singletGate")
    return(flowStats::singletGate(fr, area = channels[1], height = channels[2],
                                  filterId = filterId))
  if (method == "gate_flowclust_2d")
    return(do.call(openCyto::gate_flowclust_2d,
                   c(list(fr = fr, xChannel = channels[1], yChannel = channels[2],
                          filterId = filterId), plugin_args)))
  stop("[apply_gating2] Unsupported gating_method in v0.70 prototype: '", method,
       "'. Add a branch in .fasta_compute_gate().")
}

apply_gating2 <- function(gs, template, mc.cores = 1, verbose = TRUE,
                          collapse_max_events = 5e5L, collapse_seed = 1234L) {

  if (inherits(gs, "cytoset") || inherits(gs, "flowSet")) {
    message("[apply_gating2] cytoset/flowSet provided -- converting to GatingSet.")
    gs <- GatingSet(gs)
  }
  if (!inherits(gs, "GatingSet"))
    stop("[apply_gating2] gs must be a GatingSet, cytoset, or flowSet.")

  # Clear any pre-existing populations so re-runs start clean
  existing <- gs_get_pop_paths(gs)[-1]
  if (length(existing) > 0)
    for (p in rev(existing)) gs_pop_remove(gs, p)

  # Clear .fasta_gate_cache so the plugins' per_sample paths compute fresh from
  # the flowFrame we feed them. A prior apply_gating() run leaves per-sample
  # threshold keys here; without this they'd be read instead of recomputed and
  # silently break parity. apply_gating2 itself never reads this cache.
  rm(list = ls(envir = .fasta_gate_cache), envir = .fasta_gate_cache)

  # -- Normalize template (same prep as apply_gating) -----------------------
  td <- as.data.frame(template)
  td <- td[!is.na(td$alias) & nzchar(trimws(td$alias)), , drop = FALSE]
  if ("gating_args" %in% names(td))
    td$gating_args <- gsub("\"", "'", td$gating_args, fixed = TRUE)
  for (col in c("gating_method", "parent", "alias"))
    if (col %in% names(td)) td[[col]] <- trimws(as.character(td[[col]]))

  all_ch <- colnames(gs[[1]]); all_mk <- markernames(gs[[1]])
  sns <- sampleNames(gs); pd <- pData(gs)

  # gate_logic is index-based (booleanFilter references don't resolve in this
  # flowWorkspace build). Track the manually-set indices so they can be
  # re-applied after recompute() wipes them (recompute re-evaluates the
  # pass-through dummy gate to "all parent events").
  .logic_indices <- list()
  .logic_aliases <- character(0)
  # Re-apply just the gate_logic NODE indices. Called at the top of each loop
  # iteration so a gate parented ON a gate_logic node reads correct parent data.
  restore_logic_nodes <- function() {
    for (al in names(.logic_indices)) {
      li <- .logic_indices[[al]]
      for (sn in names(li$by_sn))
        tryCatch(gh_pop_set_indices(gs[[sn]], li$path, li$by_sn[[sn]]),
                 error = function(e) NULL)
    }
  }

  resolve_ch <- function(name) {
    if (name %in% all_ch) return(name)
    n <- names(all_mk)[match(name, all_mk)]
    if (!is.na(n)) return(n)
    stop("[apply_gating2] Cannot resolve channel/marker: '", name, "'")
  }
  resolve_pop <- function(alias) {
    if (alias == "root") return("root")
    pops <- gs_get_pop_paths(gs)
    if (grepl("/", alias, fixed = TRUE)) {
      parts <- Filter(nzchar, trimws(strsplit(alias, "/")[[1]]))
      cands <- pops[basename(pops) == parts[length(parts)]]
      if (length(parts) > 1 && length(cands) > 0) {
        prior <- paste(parts[-length(parts)], collapse = "/")
        cands <- cands[endsWith(dirname(cands), prior)]
      }
    } else cands <- pops[basename(pops) == alias]
    if (length(cands) == 0) stop("[apply_gating2] Population not found: '", alias, "'")
    if (length(cands) > 1)
      warning("[apply_gating2] Multiple matches for '", alias, "' -- using: ", cands[1])
    cands[1]
  }
  pd_col <- function(col) {  # tolerate merge .x/.y suffixes
    if (col %in% names(pd)) return(col)
    for (sfx in c(".x", ".y")) if (paste0(col, sfx) %in% names(pd)) return(paste0(col, sfx))
    stop("[apply_gating2] Column '", col, "' not in pData.")
  }
  # Parse a gating_args string into a named list. A malformed string (e.g. a
  # missing '=' like `y_min, 5000`, or an unquoted token) must FAIL LOUDLY here
  # -- silently returning list() lets the row fall through to a cryptic crash
  # deep in the external_control/pick_in path ("argument is of length zero").
  parse_args <- function(s, ctx = "?")
    if (!is.na(s) && nzchar(trimws(s)))
      tryCatch(eval(parse(text = paste0("list(", s, ")"))),
               error = function(e)
                 stop("[apply_gating2] could not parse gating_args for '", ctx, "': ",
                      conditionMessage(e), "\n  gating_args = ", s, call. = FALSE))
    else list()

  # Decompose a quadGate into 4 axis-aligned rectangleGates (order matches the
  # plugin quad_names convention: -+, ++, +-, --). Sidesteps the gs_pop_add
  # quadGate per-sample collapse (apply_gating bug #4) -- rectangle lists apply
  # per sample cleanly.
  decompose_quad <- function(qg, channels, qnames) {
    bx <- qg@boundary[[channels[1]]]; by <- qg@boundary[[channels[2]]]
    setNames(list(
      rectangleGate(setNames(list(c(-Inf, bx), c(by, Inf)), channels), filterId = qnames[1]),
      rectangleGate(setNames(list(c(bx, Inf), c(by, Inf)), channels), filterId = qnames[2]),
      rectangleGate(setNames(list(c(bx, Inf), c(-Inf, by)), channels), filterId = qnames[3]),
      rectangleGate(setNames(list(c(-Inf, bx), c(-Inf, by)), channels), filterId = qnames[4])
    ), qnames)
  }
  flip_neg_1d <- function(g, channel) {  # complement of [thr, Inf] within parent
    mn <- g@min[[channel]]; mx <- g@max[[channel]]
    thr <- if (is.finite(mn)) mn else mx
    rectangleGate(setNames(list(c(-Inf, thr)), channel), filterId = g@filterId)
  }

  vmsg <- function(...) if (isTRUE(verbose)) message(...)
  message("[apply_gating2] Gating ", nrow(td), " row(s), sequentially, gt_gating-free.")

  for (i in seq_len(nrow(td))) {
    row    <- td[i, , drop = FALSE]
    alias  <- row$alias
    method <- row$gating_method
    parent <- row$parent
    pop    <- if ("pop" %in% names(row)) trimws(as.character(row$pop)) else "+"
    is_active <- !("active_gate" %in% names(row)) || is.na(row$active_gate) ||
                 isTRUE(suppressWarnings(as.logical(row$active_gate)))
    if (!is_active) { vmsg("[apply_gating2] skip '", alias, "' (active_gate=FALSE)"); next }

    args <- parse_args(row$gating_args, alias)
    parent_path <- if (parent == "root") "root" else resolve_pop(parent)

    # A previous gate's recompute() wipes any gate_logic node back to "all parent
    # events". Restore them now so that if THIS gate is parented on a gate_logic
    # node, it reads/copies against correct parent data. (Children gated under a
    # gate_logic node are also re-evaluated in the end-of-loop restore below.)
    restore_logic_nodes()

    # ---- gate_logic: index-based AND / OR / NOT (no booleanFilter) -----------
    if (method == "gate_logic") {
      if (any(c("clusterFrom", "clusterFrom_x", "clusterFrom_y") %in% names(args)))
        stop("[apply_gating2] [gate_logic] '", alias, "': clusterFrom is not supported for ",
             "gate_logic -- it combines existing populations by index, there is no data fit ",
             "to relocate. Remove clusterFrom from gating_args.", call. = FALSE)
      logic <- toupper(trimws(args$logic %||% "AND"))
      gates <- .fasta_parse_pop_list(args$gates)
      if (length(gates) == 0)
        stop("[apply_gating2] [gate_logic] '", alias,
             "' needs gates= (comma-separated names or c(...)).")
      if (!logic %in% c("AND", "OR", "NOT"))
        stop("[apply_gating2] [gate_logic] logic must be AND, OR, or NOT. Got: '", logic, "'.")
      if (logic == "NOT" && length(gates) != 1)
        stop("[apply_gating2] [gate_logic] logic='NOT' inverts exactly ONE gate; got ",
             length(gates), ".")
      comp_paths <- vapply(gates, function(g) resolve_pop(g), character(1))
      # Pass-through dummy gate so the node exists; membership is set by index.
      dummy_ch <- tryCatch(resolve_ch(trimws(strsplit(row$dims, ",")[[1]])[1]),
                           error = function(e) all_ch[1])
      gs_pop_add(gs, rectangleGate(setNames(list(c(-Inf, Inf)), dummy_ch), filterId = alias),
                 parent = parent_path)
      recompute(gs)
      # recompute() just reset every PRIOR gate_logic node's pass-through dummy to
      # "all parent events", wiping its real membership. Restore them before we
      # read component indices below, so a gate_logic that references ANOTHER
      # gate_logic (e.g. CAR- = NOT(CAR_all)) sees correct membership, not 100%.
      for (.al in names(.logic_indices)) {
        .li <- .logic_indices[[.al]]
        for (.s in names(.li$by_sn))
          tryCatch(gh_pop_set_indices(gs[[.s]], .li$path, .li$by_sn[[.s]]),
                   error = function(e) NULL)
      }
      gl_path <- if (parent_path == "root") paste0("/", alias)
                 else paste0(parent_path, "/", alias)
      .logic_indices[[alias]] <- list(path = gl_path, by_sn = list())
      .logic_aliases <- c(.logic_aliases, alias)
      for (sn in sns) {
        idx_list <- lapply(comp_paths, function(pp) gh_pop_get_indices(gs[[sn]], pp))
        bool_idx <- switch(logic,
          AND = Reduce(`&`, idx_list),
          OR  = Reduce(`|`, idx_list),
          NOT = {
            par_idx <- if (parent_path == "root") rep(TRUE, length(idx_list[[1]]))
                       else gh_pop_get_indices(gs[[sn]], parent_path)
            par_idx & !idx_list[[1]]
          })
        gh_pop_set_indices(gs[[sn]], gl_path, bool_idx)
        .logic_indices[[alias]]$by_sn[[sn]] <- bool_idx
      }
      vmsg("[apply_gating2] [gate_logic] '", alias, "' (", logic, ": ",
           paste(basename(comp_paths), collapse = ", "), ") -> ", parent_path)
      next
    }

    # ---- copy_gate: copy gate(s) from an existing population to a new parent --
    # Handles every gate type:
    #   single-node sources (gate_static/singlet/mindensity/FMO 1D, ...) -> copy
    #     that population's per-sample gate to parent_path as `alias`.
    #   2D quad sources (gate_flowmeans_2d/FMO_2d/mindensity_2d/tmix_2d) have NO
    #     single node -- the 4 children are added directly under their parent --
    #     so copy those 4 children (rectangleGate or polygonGate) by name.
    #   gate_logic sources have no gate geometry -> rejected (recreate via gate_logic).
    if (method == "copy_gate") {
      if (any(c("clusterFrom", "clusterFrom_x", "clusterFrom_y") %in% names(args)))
        stop("[apply_gating2] [copy_gate] '", alias, "': clusterFrom is not supported for ",
             "copy_gate -- it clones an existing gate's geometry, there is no data fit ",
             "to relocate. Remove clusterFrom from gating_args.", call. = FALSE)
      src <- trimws(args$copy_from %||%
        stop("[apply_gating2] [copy_gate] '", alias, "' needs copy_from="))
      src_row <- td[td$alias == src, , drop = FALSE]

      # Source deactivated -> its populations don't exist; skip.
      if (nrow(src_row) > 0 && "active_gate" %in% names(src_row) &&
          !is.na(src_row$active_gate[1]) &&
          isFALSE(suppressWarnings(as.logical(src_row$active_gate[1])))) {
        warning("[apply_gating2] [copy_gate] '", alias, "': source '", src,
                "' has active_gate=FALSE -- skipping.")
        next
      }
      # gate_logic is index-based -> no copyable gate.
      if (src %in% .logic_aliases ||
          (nrow(src_row) > 0 && identical(trimws(src_row$gating_method[1]), "gate_logic")))
        stop("[apply_gating2] [copy_gate] '", alias, "': cannot copy gate_logic source '",
             src, "' (no gate geometry -- recreate it with a gate_logic row).")

      live <- gs_get_pop_paths(gs)
      src_path <- tryCatch(resolve_pop(src), error = function(e) NULL)

      copy_one <- function(from_path, new_name) {
        gl <- setNames(lapply(sns, function(sn) {
          g <- tryCatch(gh_pop_get_gate(gs[[sn]], from_path), error = function(e) NULL)
          if (!is.null(g)) g@filterId <- new_name
          g
        }), sns)
        if (any(vapply(gl, is.null, logical(1)))) return(FALSE)
        gs_pop_add(gs, gl, parent = parent_path)
        TRUE
      }

      if (!is.null(src_path)) {
        # Case A: source is a single existing node.
        src_children <- live[dirname(live) == src_path]
        first_g <- gh_pop_get_gate(gs[[sns[1]]], src_path)
        if (inherits(first_g, "quadGate") && length(src_children) > 0) {
          for (cp in src_children) copy_one(cp, basename(cp))   # quadGate node + children
        } else {
          copy_one(src_path, alias)                              # plain single gate
        }
      } else {
        # Case B: source is a 2D-quad alias with no own node -> copy its children.
        if (nrow(src_row) == 0)
          stop("[apply_gating2] [copy_gate] '", alias, "': source '", src,
               "' is neither a population nor a template alias.")
        src_parent <- resolve_pop(src_row$parent[1])
        src_args   <- parse_args(src_row$gating_args[1], paste0(alias, " (copy source '", src, "')"))
        nx <- src_args$name_x; ny <- src_args$name_y
        src_children <- live[dirname(live) == src_parent]
        if (!is.null(nx) && !is.null(ny))
          src_children <- src_children[
            grepl(gsub("[^A-Za-z0-9]", ".", nx), basename(src_children)) &
            grepl(gsub("[^A-Za-z0-9]", ".", ny), basename(src_children))]
        if (length(src_children) == 0)
          stop("[apply_gating2] [copy_gate] '", alias, "': no quad children found for source '",
               src, "' under '", src_parent, "'.")
        for (cp in src_children) copy_one(cp, basename(cp))
      }
      recompute(gs)
      vmsg("[apply_gating2] [copy_gate] '", alias, "' <- '", src, "' (parent ", parent_path, ")")
      next
    }

    # ---- channel resolution + source population -----------------------------
    channels <- vapply(trimws(strsplit(row$dims, ",")[[1]]), resolve_ch, character(1),
                       USE.NAMES = FALSE)
    src_pop <- if (!is.null(args$clusterFrom)) resolve_pop(args$clusterFrom) else parent_path

    per_axis <- !is.null(args$control_col_x) || !is.null(args$control_col_y)
    is_ec <- !is.null(args$control_col) || per_axis ||
             (!is.na(row$gating_args) && grepl("external_control", row$gating_args))
    .is_collapse <- "collapseDataForGating" %in% names(row) &&
                    !is.na(row$collapseDataForGating) &&
                    isTRUE(suppressWarnings(as.logical(row$collapseDataForGating)))

    # ---- build the named-by-sample gate list --------------------------------
    if (.is_collapse) {
      # collapseDataForGating (mirrors openCyto): pool the parent (or clusterFrom)
      # events ACROSS the group -- or all samples when groupBy is blank -- compute
      # ONE gate on the pooled data with any of our methods, and replicate it to
      # every sample in the group. Mutually exclusive with external_control
      # (collapse pools the union; external_control picks one control sample).
      if (is_ec)
        stop("[apply_gating2] '", alias, "': collapseDataForGating and external_control ",
             "are mutually exclusive -- collapse pools the union of events, external_control ",
             "picks one control sample. Use one.", call. = FALSE)
      grp_spec <- if (!is.null(args$groupBy)) args$groupBy
                  else if ("groupBy" %in% names(row) && !is.na(row$groupBy) &&
                           nzchar(trimws(as.character(row$groupBy)))) trimws(as.character(row$groupBy))
                  else NULL
      # Resolve groups: blank -> one group (all); numeric N -> every N samples;
      # else one or more pData columns (":"/"," separated) -> unique combination.
      groups <- if (is.null(grp_spec)) list(seq_len(nrow(pd))) else {
        .gspec <- trimws(as.character(grp_spec))
        if (grepl("^[0-9]+$", .gspec)) {
          .n <- as.integer(.gspec)
          unname(split(seq_len(nrow(pd)), ceiling(seq_len(nrow(pd)) / .n)))
        } else {
          .cols <- vapply(trimws(strsplit(.gspec, "[:,]")[[1]]), pd_col, character(1))
          .key  <- do.call(paste, c(lapply(.cols, function(cc) as.character(pd[[cc]])), sep = "|"))
          unname(split(seq_len(nrow(pd)), .key))
        }
      }
      # Pool one group's events (only the gating channels), reproducibly
      # downsampled to bound memory/perf. Per-sample cap balances representation;
      # the global RNG is saved/restored so seeding here is local (won't perturb
      # tmix/flowClust randomness elsewhere).
      .pool_group <- function(g_sns) {
        cap_each <- max(1L, floor(collapse_max_events / length(g_sns)))
        .old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
                       get(".Random.seed", envir = .GlobalEnv) else NULL
        set.seed(collapse_seed)
        mats <- lapply(g_sns, function(sn) {
          fr <- gh_pop_get_data(gs[[sn]], src_pop)
          m  <- exprs(fr)[, channels, drop = FALSE]
          if (nrow(m) > cap_each) m <- m[sort(sample.int(nrow(m), cap_each)), , drop = FALSE]
          m
        })
        if (is.null(.old_seed)) suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
        else assign(".Random.seed", .old_seed, envir = .GlobalEnv)
        pooled <- do.call(rbind, mats)
        if (nrow(pooled) == 0)
          stop("[apply_gating2] '", alias, "': 0 pooled events for collapse gating.", call. = FALSE)
        flowFrame(pooled)
      }
      gate_list <- setNames(vector("list", length(sns)), sns)
      for (gi in seq_along(groups)) {
        g_sns <- rownames(pd)[groups[[gi]]]; g_sns <- g_sns[g_sns %in% sns]
        if (!length(g_sns)) next
        pooled_fr <- .pool_group(g_sns)
        g <- .fasta_compute_gate(pooled_fr, method, channels, alias, args)
        vmsg("[apply_gating2] [collapse",
             if (!is.null(grp_spec)) paste0("|grp ", gi, "/", length(groups)) else "",
             "] '", alias, "' from ", length(g_sns), " pooled sample(s), ",
             nrow(exprs(pooled_fr)), " events")
        for (sn in g_sns) gate_list[[sn]] <- g
      }
    } else if (!is_ec) {
      # per_sample (or clusterFrom per_sample): compute from each sample's own data
      gate_list <- setNames(lapply(sns, function(sn) {
        fr <- gh_pop_get_data(gs[[sn]], src_pop)
        if (nrow(exprs(fr)) == 0)
          stop("[apply_gating2] '", alias, "': sample '", sn, "' has 0 events in '", src_pop, "'.")
        .fasta_compute_gate(fr, method, channels, alias, args)
      }), sns)
    } else {
      grp <- if (!is.null(args$groupBy)) pd_col(args$groupBy)
             else if ("groupBy" %in% names(row) && !is.na(row$groupBy) && nzchar(trimws(row$groupBy)))
               pd_col(trimws(row$groupBy)) else NULL
      # Pick one control sample (row index) matching col==val, restricted to
      # idx_pool (the current group, or all samples). control_idx breaks ties.
      pick_in <- function(col, val, idx_pool, idx_arg) {
        if (length(col) == 0 || is.null(col) || (is.character(col) && !nzchar(trimws(col))))
          stop("[apply_gating2] '", alias, "': this row is treated as external_control ",
               "(gating_args mentions it), but no control column was parsed. Check that ",
               "control_col / control_col_x / control_col_y are spelled correctly and the ",
               "gating_args parsed cleanly. gating_args = ", row$gating_args, call. = FALSE)
        hits <- intersect(idx_pool, which(grepl(val, pd[[pd_col(col)]], ignore.case = TRUE)))
        if (length(hits) == 0)
          stop("[apply_gating2] '", alias, "': no control matches ", col, "='", val, "'",
               if (!is.null(grp)) " in group" else "", ".")
        ci <- if (!is.null(idx_arg)) as.integer(idx_arg) else 1L
        hits[min(ci, length(hits))]
      }
      cgate <- function(sn) .fasta_compute_gate(gh_pop_get_data(gs[[sn]], src_pop),
                                                method, channels, alias, args)
      # Build the gate to broadcast over the samples in idx_pool.
      make_gate <- function(idx_pool, tag) {
        if (per_axis) {
          # x and y thresholds come from DIFFERENT control samples. Fit each
          # control with the same plugin and take only that axis's threshold.
          cs_x <- rownames(pd)[pick_in(args$control_col_x, args$control_val_x, idx_pool, args$control_idx_x)]
          cs_y <- rownames(pd)[pick_in(args$control_col_y, args$control_val_y, idx_pool, args$control_idx_y)]
          gx <- cgate(cs_x); gy <- cgate(cs_y)
          if (!inherits(gx, "quadGate") || !inherits(gy, "quadGate"))
            stop("[apply_gating2] per-axis control needs a 2D quad method (alias '", alias, "').")
          thr_x <- gx@boundary[[channels[1]]]; thr_y <- gy@boundary[[channels[2]]]
          vmsg("[apply_gating2] [ext_ctrl", tag, "|per-axis] '", alias,
               "' x<-", cs_x, " (", round(thr_x, 1), ") y<-", cs_y, " (", round(thr_y, 1), ")")
          quadGate(setNames(list(thr_x, thr_y), channels), filterId = alias)
        } else {
          cs <- rownames(pd)[pick_in(args$control_col, args$control_val, idx_pool, args$control_idx)]
          vmsg("[apply_gating2] [ext_ctrl", tag, "] '", alias, "' from control: ", cs)
          cgate(cs)
        }
      }
      gate_list <- setNames(vector("list", length(sns)), sns)
      if (is.null(grp)) {
        g <- make_gate(seq_len(nrow(pd)), "")
        for (sn in sns) gate_list[[sn]] <- g
      } else {
        for (gv in unique(pd[[grp]])) {
          grp_idx <- which(pd[[grp]] == gv)
          g <- make_gate(grp_idx, paste0("|", gv))
          for (sn in rownames(pd)[grp_idx]) gate_list[[sn]] <- g
        }
      }
    }

    # ---- add to the GatingSet, by gate type ---------------------------------
    g1 <- gate_list[[1]]

    if (pop %in% c("+-", "-+")) {
      # 1D split: from one shared per-sample threshold, add BOTH the positive
      # population ([thr, Inf]) and the negative population ([-Inf, thr]) as two
      # rectangleGates, named "<alias>+" and "<alias>-". "+-" and "-+" are
      # equivalent. Strictly 1D (single-channel rectangleGate).
      if (!(inherits(g1, "rectangleGate") && length(channels) == 1))
        stop("[apply_gating2] pop='", pop, "' is only supported for 1D gates ",
             "(a single-channel rectangleGate); alias '", alias, "', method '", method, "'.")
      pos_id <- paste0(alias, "+"); neg_id <- paste0(alias, "-")
      pos <- setNames(lapply(sns, function(sn) { g <- gate_list[[sn]]; g@filterId <- pos_id; g }), sns)
      neg <- setNames(lapply(sns, function(sn) { g <- flip_neg_1d(gate_list[[sn]], channels[1]); g@filterId <- neg_id; g }), sns)
      gs_pop_add(gs, pos, parent = parent_path)
      gs_pop_add(gs, neg, parent = parent_path)
      recompute(gs)

    } else if (inherits(g1, "filters")) {
      # gate_tmix_2d: each sample's gate is a `filters` of 4 polygonGates.
      qn <- vapply(seq_len(4), function(j) g1[[j]]@filterId, character(1))
      for (jq in seq_len(4)) {
        per_sn <- setNames(lapply(sns, function(sn) {
          pg <- gate_list[[sn]][[jq]]; pg@filterId <- qn[jq]; pg
        }), sns)
        gs_pop_add(gs, per_sn, parent = parent_path)
      }
      recompute(gs)

    } else if (inherits(g1, "quadGate")) {
      qx <- args$name_x %||% channels[1]; qy <- args$name_y %||% channels[2]
      qn <- c(paste0(qx, "-", qy, "+"), paste0(qx, "+", qy, "+"),
              paste0(qx, "+", qy, "-"), paste0(qx, "-", qy, "-"))
      for (jq in seq_len(4)) {
        per_sn <- setNames(lapply(sns, function(sn)
          decompose_quad(gate_list[[sn]], channels, qn)[[jq]]), sns)
        gs_pop_add(gs, per_sn, parent = parent_path)
      }
      recompute(gs)

    } else if (inherits(g1, "rectangleGate") && length(channels) == 1 && pop == "-") {
      neg <- setNames(lapply(sns, function(sn) flip_neg_1d(gate_list[[sn]], channels[1])), sns)
      gs_pop_add(gs, neg, parent = parent_path); recompute(gs)

    } else if (inherits(g1, "polygonGate") && pop == "-") {
      stop("[apply_gating2] polygon negative (pop='-') not yet in v0.70 (alias '", alias, "').")

    } else {
      # 1D/2D rectangleGate (pop='+'), singlet/static polygon, etc.
      gs_pop_add(gs, gate_list, parent = parent_path); recompute(gs)
    }
    vmsg("[apply_gating2] done: '", alias, "' (", method, ", parent ", parent_path, ")")
  }

  # Finalize gate_logic: restore each node's manual indices (recompute wiped
  # them), then RE-GATE its descendants against the restored parent -- a child
  # gated under a gate_logic node (e.g. a copy_gate of CD4_CD8 under CAR_all) was
  # computed against the wiped parent during the loop. Re-evaluate each child's
  # gate on root data and intersect with its restored parent indices. Process
  # shallow-to-deep (nchar order) so a parent is restored before its children.
  if (length(.logic_indices) > 0) {
    all_paths <- gs_get_pop_paths(gs)
    for (al in names(.logic_indices)) {
      li <- .logic_indices[[al]]
      gl_children <- all_paths[startsWith(all_paths, paste0(li$path, "/"))]
      for (sn in names(li$by_sn)) {
        gh_pop_set_indices(gs[[sn]], li$path, li$by_sn[[sn]])
        if (length(gl_children) > 0) {
          fr_root <- gh_pop_get_data(gs[[sn]], "root")
          n_root  <- nrow(exprs(fr_root))
          for (cp in gl_children[order(nchar(gl_children))]) {
            g_obj <- tryCatch(gh_pop_get_gate(gs[[sn]], cp), error = function(e) NULL)
            if (is.null(g_obj)) next
            g_res <- tryCatch(as(flowCore::filter(fr_root, g_obj), "logical"),
                              error = function(e) rep(FALSE, n_root))
            par_idx <- tryCatch(gh_pop_get_indices(gs[[sn]], dirname(cp)),
                                error = function(e) li$by_sn[[sn]])
            gh_pop_set_indices(gs[[sn]], cp, par_idx & g_res)
          }
        }
      }
    }
    vmsg("[apply_gating2] Finalized ", length(.logic_indices),
         " gate_logic gate(s) + descendants.")
  }

  pops <- gs_get_pop_paths(gs)
  message("[apply_gating2] Done. ", length(pops) - 1, " population(s) created.")
  invisible(gs)
}
