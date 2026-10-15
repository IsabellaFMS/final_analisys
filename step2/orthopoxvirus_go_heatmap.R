# Builds the GO heatmap figure: top 10 terms per category (BP, MF, CC) for each of the four viruses.

# Keep R diagnostics in English regardless of the machine locale.
Sys.setenv(LANGUAGE = "en")

user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
default_user_lib <- file.path(Sys.getenv("USERPROFILE", unset = path.expand("~")), "R", "win-library", "4.5")

if (!nzchar(user_lib) || !dir.exists(user_lib)) {
  user_lib <- default_user_lib
}
.libPaths(unique(c(user_lib, .libPaths())))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(plotly)
  library(htmlwidgets)
})

VIRUS_ORDER <- c("CPXV", "MPXV", "VACV", "VARV")

# Same column-to-species mapping the Python package keeps in config.py.
SPECIES_COLUMNS <- c(
  cowpox_refseq_Brighton_Red_proteins = "CPXV",
  mpox_refseq_Rivers_State_proteins = "MPXV",
  vaccinia_WR_refseq_proteins = "VACV",
  varv_refseq_proteins = "VARV"
)

GO_NAMESPACES <- c("biological_process", "molecular_function", "cellular_component")

resolve_base_dir <- function() {
  project_root <- Sys.getenv("OPXV_PROJECT_ROOT", unset = "")
  if (nzchar(project_root) && dir.exists(project_root)) {
    return(normalizePath(project_root))
  }


  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1]))))
  }

  normalizePath(getwd())
}

read_bytes_as_text <- function(path) {
  # Reading raw bytes keeps the vendored widget assets intact regardless of the
  # encoding they declare, which readLines would otherwise re-encode.
  rawToChar(readBin(path, what = "raw", n = file.info(path)$size))
}

parse_go_ontology <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")

  markers <- lines %in% c("[Term]", "[Typedef]")
  stanza <- cumsum(markers)
  stanza_kind <- c("", lines[markers])[stanza + 1L]
  in_term <- stanza_kind == "[Term]"
  total <- max(stanza)

  field <- function(prefix, pattern) {
    hits <- which(in_term & startsWith(lines, prefix))
    list(stanza = stanza[hits], value = sub(pattern, "\\1", lines[hits]))
  }

  ids <- field("id: GO:", "^id: (GO:[0-9]+).*$")
  labels <- field("name: ", "^name: (.*)$")
  namespaces <- field("namespace: ", "^namespace: (.*)$")
  alt_ids <- field("alt_id: GO:", "^alt_id: (GO:[0-9]+).*$")
  parents <- field("is_a: GO:", "^is_a: (GO:[0-9]+).*$")

  term_id <- rep(NA_character_, total)
  term_id[ids$stanza] <- ids$value
  term_name <- rep(NA_character_, total)
  term_name[labels$stanza] <- labels$value
  term_namespace <- rep("", total)
  term_namespace[namespaces$stanza] <- namespaces$value
  term_obsolete <- rep(FALSE, total)
  term_obsolete[stanza[which(in_term & lines == "is_obsolete: true")]] <- TRUE

  complete <- which(!is.na(term_id) & !is.na(term_name))

  # An alternate identifier denotes the same term, so it is registered as an
  # extra key pointing at the primary entry.
  alt_owner <- term_id[alt_ids$stanza]
  usable_alt <- !is.na(alt_owner)
  owner_index <- match(alt_owner[usable_alt], term_id[complete])

  keys <- c(term_id[complete], alt_ids$value[usable_alt])
  names_by_key <- c(term_name[complete], term_name[complete][owner_index])
  namespace_by_key <- c(term_namespace[complete], term_namespace[complete][owner_index])
  obsolete_by_key <- c(term_obsolete[complete], term_obsolete[complete][owner_index])

  # A later stanza wins whenever a key repeats.
  last <- !duplicated(keys, fromLast = TRUE)

  list(
    name = setNames(names_by_key[last], keys[last]),
    namespace = setNames(namespace_by_key[last], keys[last]),
    obsolete = setNames(obsolete_by_key[last], keys[last]),
    parents = split(parents$value, factor(term_id[parents$stanza])),
    alias = setNames(alt_owner[usable_alt], alt_ids$value[usable_alt])
  )
}

build_ancestor_lookup <- function(ontology) {
  cache <- new.env(parent = emptyenv())

  resolve <- function(go_id) {
    cached <- cache[[go_id]]
    if (!is.null(cached)) {
      return(cached)
    }
    # Seeding the cache before recursing keeps a malformed cycle from hanging.
    assign(go_id, character(0), envir = cache)

    alias <- unname(ontology$alias[go_id])
    key <- if (is.na(alias)) go_id else alias

    found <- character(0)
    for (parent in ontology$parents[key][[1]]) {
      found <- c(found, parent, resolve(parent))
    }
    found <- unique(found)
    assign(go_id, found, envir = cache)
    found
  }

  resolve
}

load_protein_go_terms <- function(path, ontology) {
  header <- strsplit(readLines(path, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]]
  wanted <- c("Protein_Accession", "GO_BP", "GO_CC", "GO_MF")
  missing_columns <- setdiff(wanted, header)
  if (length(missing_columns) > 0) {
    stop(paste("Columns missing from", path, ":", paste(missing_columns, collapse = ", ")))
  }

  # Only the accession and the three namespace columns are needed, and skipping
  # the rest avoids parsing the long description fields.
  consolidated <- read.delim(
    path,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    colClasses = ifelse(header %in% wanted, "character", "NULL")
  )

  # Terms are stored as GO:0004222 (metalloendopeptidase activity), and proteins
  # without annotations carry a "-", so the identifiers are pulled out directly.
  namespaces <- paste(consolidated$GO_BP, consolidated$GO_CC, consolidated$GO_MF)
  matches <- regmatches(namespaces, gregexpr("GO:[0-9]+", namespaces))

  data.frame(
    Protein_ID = rep(trimws(consolidated$Protein_Accession), lengths(matches)),
    GO_ID = unlist(matches),
    stringsAsFactors = FALSE
  ) %>%
    filter(
      GO_ID %in% names(ontology$name),
      !ontology$obsolete[GO_ID],
      !startsWith(ontology$name[GO_ID], "obsolete ")
    ) %>%
    distinct(Protein_ID, GO_ID)
}

keep_leaf_terms <- function(protein_terms, ancestors_of) {
  protein_terms %>%
    group_by(Protein_ID) %>%
    # A term is dropped when a more specific term on the same protein already
    # implies it, which is what makes the counts reflect leaf annotations.
    filter({
      inherited <- unique(unlist(lapply(GO_ID, ancestors_of)))
      !GO_ID %in% inherited
    }) %>%
    ungroup()
}

read_orthogroup_map <- function(paths) {
  frames <- list()
  for (path in paths) {
    if (!file.exists(path)) {
      next
    }
    table <- read.delim(path, stringsAsFactors = FALSE, quote = "", comment.char = "")
    for (column in intersect(names(SPECIES_COLUMNS), colnames(table))) {
      # A column with no protein at all is read as logical, so it is coerced
      # back before being split.
      cells <- as.character(table[[column]])
      proteins <- trimws(unlist(strsplit(cells[!is.na(cells)], ",", fixed = TRUE)))
      proteins <- proteins[nzchar(proteins)]
      if (length(proteins) > 0) {
        frames[[length(frames) + 1]] <- data.frame(
          Protein_ID = proteins,
          Virus = SPECIES_COLUMNS[[column]],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(frames) == 0) {
    stop(paste(
      "No OrthoFinder table could be read. Expected one of:",
      paste(paths, collapse = ", ")
    ))
  }

  distinct(bind_rows(frames))
}

embed_widget_assets <- function(html_file) {
  # saveWidget can only inline its assets through pandoc, which this project
  # does not require, so the vendored files are folded in directly instead.
  text <- read_bytes_as_text(html_file)
  base <- dirname(html_file)

  inline <- function(text, pattern, open, close) {
    matches <- gregexpr(pattern, text, perl = TRUE)
    tags <- regmatches(text, matches)[[1]]
    if (length(tags) == 0) {
      return(text)
    }
    replacements <- vapply(tags, function(tag) {
      asset <- file.path(base, sub(pattern, "\\1", tag, perl = TRUE))
      if (!file.exists(asset)) {
        return(tag)
      }
      paste0(open, "\n", read_bytes_as_text(asset), "\n", close)
    }, character(1), USE.NAMES = FALSE)
    regmatches(text, matches) <- list(replacements)
    text
  }

  text <- inline(text, '<script src="([^"]+)"[^>]*></script>', "<script>", "</script>")
  text <- inline(text, '<link href="([^"]+)"[^>]*rel="stylesheet"[^>]*/?>', "<style>", "</style>")
  writeBin(charToRaw(text), html_file)
}

save_standalone_widget <- function(widget, output_html) {
  staging <- file.path(tempdir(), "opxv_heatmap_widget")
  unlink(staging, recursive = TRUE)
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)

  staged_html <- file.path(staging, basename(output_html))
  saveWidget(widget, file = staged_html, selfcontained = FALSE, libdir = "lib")
  embed_widget_assets(staged_html)

  dir.create(dirname(output_html), recursive = TRUE, showWarnings = FALSE)
  file.copy(staged_html, output_html, overwrite = TRUE)
  unlink(staging, recursive = TRUE)
}

base_dir <- resolve_base_dir()
output_dir <- file.path(base_dir, "outputs", "visualizations")
tables_dir <- file.path(base_dir, "outputs", "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

annotation_table <- file.path(tables_dir, "step3_final_consolidated_annotations.tsv")
obo_file <- file.path(base_dir, "sources", "go-basic.obo")
# The protein-to-species assignment comes from the step 1 OrthoFinder tables.
orthogroup_files <- file.path(
  dirname(base_dir), "step1", "sources",
  c("orthogroups.tsv", "orthogroups_unassigned_genes.tsv")
)
output_html <- file.path(output_dir, "orthopoxvirus_go_heatmap.html")

cat("Base directory:", base_dir, "\n")

if (!file.exists(annotation_table)) {
  stop(paste(
    "Consolidated annotation table not found:", annotation_table,
    "\nRun 'python run_analysis.py pipeline' before the heatmap step."
  ))
}

if (!file.exists(obo_file)) {
  stop(paste("The file go-basic.obo was not found in", dirname(obo_file)))
}

cat("Reading the Gene Ontology...\n")
ontology <- parse_go_ontology(obo_file)
ancestors_of <- build_ancestor_lookup(ontology)
cat("GO terms loaded:", length(ontology$name), "\n")

cat("Reading the consolidated annotations...\n")
protein_terms <- load_protein_go_terms(annotation_table, ontology)
leaf_terms <- keep_leaf_terms(protein_terms, ancestors_of)
cat("Proteins with GO annotations:", n_distinct(leaf_terms$Protein_ID), "\n")

library_data <- read_orthogroup_map(orthogroup_files)

go_expanded <- leaf_terms %>%
  inner_join(library_data, by = "Protein_ID", relationship = "many-to-many") %>%
  mutate(
    GO_Term = unname(ontology$name[GO_ID]),
    Category = unname(ontology$namespace[GO_ID])
  ) %>%
  filter(
    Category %in% GO_NAMESPACES,
    !tolower(gsub(" ", "_", GO_Term)) %in% GO_NAMESPACES
  )

cat("Proteins with GO annotations and virus labels:", n_distinct(go_expanded$Protein_ID), "\n\n")

go_freq <- count(go_expanded, GO_Term, Category, name = "Frequency")

# Ties are broken by term name so the selection stays reproducible.
top_10_per_cat <- go_freq %>%
  group_by(Category) %>%
  arrange(desc(Frequency), GO_Term, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

cat("Top 10 terms per category selected\n\n")

matrix_data <- go_expanded %>%
  filter(GO_Term %in% top_10_per_cat$GO_Term) %>%
  count(GO_Term, Virus, name = "Frequency") %>%
  mutate(Virus = factor(Virus, levels = VIRUS_ORDER)) %>%
  pivot_wider(
    names_from = Virus,
    values_from = Frequency,
    values_fill = 0,
    names_expand = TRUE
  ) %>%
  arrange(GO_Term)

matrix_numeric <- as.matrix(matrix_data[, VIRUS_ORDER])
rownames(matrix_numeric) <- matrix_data$GO_Term

write.table(
  matrix_data,
  file = file.path(output_dir, "orthopoxvirus_go_heatmap_matrix.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  top_10_per_cat,
  file = file.path(output_dir, "orthopoxvirus_go_heatmap_top10_terms.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Auxiliary matrices written to disk\n\n")

# Each category is introduced by a bold label sitting on an empty row, which is
# what separates the three blocks in the figure.
labels_final <- c()
matrix_final <- matrix(nrow = 0, ncol = ncol(matrix_numeric))
colnames(matrix_final) <- colnames(matrix_numeric)

for (category in GO_NAMESPACES) {
  category_terms <- top_10_per_cat %>%
    filter(Category == category) %>%
    arrange(desc(Frequency), GO_Term) %>%
    pull(GO_Term)

  if (length(category_terms) == 0) {
    next
  }

  labels_final <- c(labels_final, paste0("<b>", toupper(gsub("_", " ", category)), "</b>"))
  spacer <- matrix(NA, nrow = 1, ncol = ncol(matrix_numeric))
  colnames(spacer) <- colnames(matrix_numeric)
  matrix_final <- rbind(matrix_final, spacer)

  labels_final <- c(labels_final, category_terms)
  matrix_final <- rbind(matrix_final, matrix_numeric[category_terms, , drop = FALSE])
}

# Plotly draws the y axis from the bottom up.
labels_final <- rev(labels_final)
matrix_final <- matrix_final[nrow(matrix_final):1, , drop = FALSE]

fig <- plot_ly(
    height = 900,
    width = 1000,
    z = matrix_final,
    x = colnames(matrix_final),
    y = labels_final,
    type = "heatmap",
    colorscale = list(
      c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1),
      c("#eff3ff", "#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c", "#08306b")
    ),
    zmin = 0,
    zmax = max(matrix_final, na.rm = TRUE),
    hovertemplate = "<b>%{y}</b><br><b>%{x}</b><br>Occurrences: %{z}<extra></extra>",
    colorbar = list(
      title = "Term<br>Frequency",
      thickness = 20,
      len = 1.0,
      x = 1.02
    )
  ) %>%
  layout(
    title = list(
      text = "Top 10 GO terms per category",
      x = 0.675,
      xanchor = "center",
      font = list(size = 15, color = "#1a1a1a")
    ),
    font = list(family = "Arial, sans-serif", size = 11),
    xaxis = list(
      title = "<i>Orthopoxvirus</i>",
      tickfont = list(size = 13, color = "#1a1a1a"),
      side = "bottom"
    ),
    yaxis = list(
      title = "Gene Ontology Terms",
      tickfont = list(size = 10),
      titlefont = list(size = 12),
      standoff = 350
    ),
    plot_bgcolor = "white",
    paper_bgcolor = "white",
    margin = list(t = 70, l = 500, r = 150, b = 80)
  )

save_standalone_widget(fig, output_html)

cat("Heatmap written to:", output_html, "\n")
cat("Analysis completed\n\n")
