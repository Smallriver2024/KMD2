#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(parallel)
})

parse_args <- function(args) {
  opts <- list(
    input_dir = "kmd/Survival",
    output_dir = "kmd/Survival_2_0/plot_ggsurvplot",
    pfile = NULL,
    kinase_sql = "kmd/SQL/kinase_basic.sql",
    width = 1000,
    height = 720,
    res = 120,
    break_x_by = "auto",
    short_break_x_by = 1,
    long_break_x_by = 2,
    long_year_threshold = 15,
    days_per_year = 365.25,
    limit = Inf,
    only = NULL,
    workers = 1,
    overwrite = FALSE
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--overwrite")) {
      opts$overwrite <- TRUE
      i <- i + 1
      next
    }
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Invalid argument near: ", key)
    }
    value <- args[[i + 1]]
    name <- sub("^--", "", key)
    name <- gsub("-", "_", name)
    if (!name %in% names(opts)) {
      stop("Unknown option: ", key)
    }
    opts[[name]] <- value
    i <- i + 2
  }

  opts$width <- as.integer(opts$width)
  opts$height <- as.integer(opts$height)
  opts$res <- as.integer(opts$res)
  opts$break_x_by <- as.character(opts$break_x_by)
  if (tolower(opts$break_x_by) != "auto" && is.na(suppressWarnings(as.numeric(opts$break_x_by)))) {
    stop("--break-x-by must be numeric or auto")
  }
  opts$short_break_x_by <- as.numeric(opts$short_break_x_by)
  opts$long_break_x_by <- as.numeric(opts$long_break_x_by)
  opts$long_year_threshold <- as.numeric(opts$long_year_threshold)
  opts$days_per_year <- as.numeric(opts$days_per_year)
  opts$limit <- if (is.infinite(opts$limit)) Inf else as.integer(opts$limit)
  opts$workers <- max(1L, as.integer(opts$workers))
  if (is.null(opts$pfile)) {
    opts$pfile <- file.path(opts$input_dir, "survival_p_clean.txt")
  }
  if (!is.null(opts$only)) {
    opts$only <- unlist(strsplit(opts$only, ",", fixed = TRUE))
    opts$only <- trimws(opts$only)
    opts$only <- opts$only[nzchar(opts$only)]
  }
  opts
}

read_survival_records <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  keys <- vapply(parts, `[[`, character(1), 1)
  pvals <- sub("^Logrank p-value=", "", vapply(parts, `[[`, character(1), 2))
  cancer <- sub("^.*_", "", keys)
  gene <- sub("_[^_]+$", "", keys)
  data.frame(key = keys, gene = gene, cancer = cancer, pval = pvals, stringsAsFactors = FALSE)
}

read_gene_id_map <- function(path) {
  split_sql_tuple <- function(line) {
    line <- sub("^\\(", "", line)
    line <- sub("\\),?$", "", line)
    line <- sub("\\);?$", "", line)
    chars <- strsplit(line, "", fixed = TRUE, useBytes = TRUE)[[1]]
    fields <- character()
    current <- character()
    in_quote <- FALSE
    escape <- FALSE
    for (ch in chars) {
      if (escape) {
        current <- c(current, ch)
        escape <- FALSE
        next
      }
      if (in_quote && ch == "\\") {
        escape <- TRUE
        next
      }
      if (ch == "'") {
        in_quote <- !in_quote
        next
      }
      if (!in_quote && ch == ",") {
        fields <- c(fields, trimws(paste0(current, collapse = "")))
        current <- character()
        next
      }
      current <- c(current, ch)
    }
    fields <- c(fields, trimws(paste0(current, collapse = "")))
    fields[fields == "NULL"] <- NA_character_
    fields
  }

  lines <- readLines(path, warn = FALSE)
  lines <- grep("^\\(", lines, value = TRUE)
  rows <- lapply(lines, split_sql_tuple)
  out <- data.frame(
    gene = vapply(rows, function(x) x[[2]], character(1)),
    entrez_id = vapply(rows, function(x) x[[7]], character(1)),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$gene) & !is.na(out$entrez_id), ]
  out[!duplicated(out$gene), ]
}

choose_break_x_by <- function(max_follow_up_years, break_x_by, short_break_x_by, long_break_x_by, long_year_threshold) {
  if (tolower(break_x_by) != "auto") {
    return(as.numeric(break_x_by))
  }
  if (is.finite(max_follow_up_years) && max_follow_up_years > long_year_threshold) {
    return(long_break_x_by)
  }
  short_break_x_by
}

plot_one <- function(
  record, gene_id_map, input_dir, output_dir, width, height, res,
  break_x_by, short_break_x_by, long_break_x_by, long_year_threshold,
  days_per_year, overwrite
) {
  entrez <- gene_id_map$entrez_id[match(record$gene, gene_id_map$gene)]
  if (is.na(entrez) || !nzchar(entrez)) {
    stop("No Entrez ID for ", record$gene)
  }

  group_file <- file.path(input_dir, "group_data", paste0(record$cancer, "_", entrez, "_50_50.csv"))
  if (!file.exists(group_file)) {
    stop("Missing group_data file: ", group_file)
  }

  out_file <- file.path(output_dir, paste0(record$gene, "_", record$cancer, ".png"))
  if (file.exists(out_file) && !overwrite) {
    return("skipped_exists")
  }

  dat <- read.csv(group_file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("Patient", "Days", "Status", "Expression", "Group")
  if (!all(required %in% names(dat))) {
    stop("Unexpected CSV header in ", group_file)
  }
  dat <- dat[!is.na(dat$Days) & !is.na(dat$Status) & !is.na(dat$Group), ]
  dat$Days <- as.numeric(dat$Days)
  dat$FollowUpYears <- dat$Days / days_per_year
  dat$event <- ifelse(dat$Status == "Dead", 1, 0)
  dat$Group <- factor(dat$Group, levels = c("Low", "High"))
  max_follow_up_years <- max(dat$FollowUpYears, na.rm = TRUE)
  actual_break_x_by <- choose_break_x_by(
    max_follow_up_years,
    break_x_by,
    short_break_x_by,
    long_break_x_by,
    long_year_threshold
  )

  fit <- survfit(Surv(FollowUpYears, event) ~ Group, data = dat)
  p_label <- paste0("Log-rank p = ", record$pval)

  gp <- ggsurvplot(
    fit,
    data = dat,
    conf.int = TRUE,
    pval = p_label,
    surv.median.line = "hv",
    risk.table = TRUE,
    risk.table.height = 0.24,
    risk.table.y.text = FALSE,
    xlab = "Follow up time(years)",
    ylab = "Survival probability",
    legend = c(0.80, 0.75),
    legend.title = "",
    legend.labs = c("Low expression", "High expression"),
    break.x.by = actual_break_x_by,
    palette = c("#2E86AB", "#D1495B"),
    censor.shape = 124,
    censor.size = 3,
    size = 1.05,
    ggtheme = theme_classic(base_size = 14) +
      theme(
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "grey20"),
        legend.background = element_rect(fill = scales::alpha("white", 0.8), color = NA),
        plot.margin = margin(8, 12, 8, 8)
      ),
    tables.theme = theme_cleantable(base_size = 12)
  )
  gp$plot <- gp$plot +
    ggtitle(paste0(record$gene, " in ", record$cancer)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))

  png(out_file, width = width, height = height, res = res, type = "cairo-png")
  on.exit(dev.off(), add = TRUE)
  print(gp)
  paste0("written; break_x_by=", actual_break_x_by, "; max_year=", round(max_follow_up_years, 2))
}

render_record <- function(idx, records, gene_id_map, opts) {
  rec <- records[idx, ]
  status <- "written"
  msg <- ""
  tryCatch(
    {
      status <- plot_one(
        rec, gene_id_map, opts$input_dir, opts$output_dir,
        opts$width, opts$height, opts$res, opts$break_x_by,
        opts$short_break_x_by, opts$long_break_x_by, opts$long_year_threshold,
        opts$days_per_year, opts$overwrite
      )
      if (startsWith(status, "written; ")) {
        msg <- sub("^written; ", "", status)
        status <- "written"
      }
    },
    error = function(e) {
      status <<- "error"
      msg <<- conditionMessage(e)
    }
  )
  data.frame(key = rec$key, status = status, message = msg, stringsAsFactors = FALSE)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

  records <- read_survival_records(opts$pfile)
  if (!is.null(opts$only)) {
    records <- records[records$key %in% opts$only, ]
  }
  if (is.finite(opts$limit)) {
    records <- head(records, opts$limit)
  }
  gene_id_map <- read_gene_id_map(opts$kinase_sql)

  log_path <- file.path(opts$output_dir, "render_log.tsv")
  con <- file(log_path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("key\tstatus\tmessage", con)

  total <- nrow(records)
  message(
    "Rendering ", total, " KM plot(s) to ", opts$output_dir,
    " using ", opts$workers, " worker(s); x-axis unit = years; break_x_by = ", opts$break_x_by
  )

  if (total == 0) {
    results <- data.frame(key = character(), status = character(), message = character())
  } else if (opts$workers == 1L) {
    results <- do.call(
      rbind,
      lapply(seq_len(total), function(idx) {
        out <- render_record(idx, records, gene_id_map, opts)
        if (idx %% 100 == 0 || idx == total) {
          message(idx, "/", total, " complete")
        }
        out
      })
    )
  } else {
    results <- do.call(
      rbind,
      parallel::mclapply(
        seq_len(total),
        render_record,
        records = records,
        gene_id_map = gene_id_map,
        opts = opts,
        mc.cores = opts$workers,
        mc.preschedule = FALSE
      )
    )
  }
  apply(results, 1, function(row) writeLines(paste(row, collapse = "\t"), con))
  message("Done. Log: ", log_path)
}

main()
