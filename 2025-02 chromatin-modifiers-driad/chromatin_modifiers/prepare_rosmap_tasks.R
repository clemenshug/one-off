library(tidyverse)
library(synExtra)
library(data.table)
library(powerjoin)
library(here)
library(qs)

synapser::synLogin()
syn <- synDownloader("~/data", .cache = TRUE)

# fnROSMAP <- DRIAD::wrangleROSMAP("~/data")

rosmap_tdp43_classification <- syn("syn44277761") %>%
  read_csv()

# > rosmap_tdp43_classification %>% count(brain_region, class_low)
# # A tibble: 6 × 3
# brain_region class_low     n
# <chr>        <chr>     <int>
#   1 DLPFC        negative     24
# 2 DLPFC        positive    193
# 3 HCN          negative    157
# 4 HCN          positive    512
# 5 PCC          negative     85
# 6 PCC          positive    461
# > rosmap_tdp43_classification %>% count(brain_region, class_high)
# # A tibble: 6 × 3
# brain_region class_high     n
# <chr>        <chr>      <int>
#   1 DLPFC        negative     102
# 2 DLPFC        positive     115
# 3 HCN          negative     511
# 4 HCN          positive     158
# 5 PCC          negative     342
# 6 PCC          positive     204

rosmap_quant_clinical <- syn("syn44335073") %>%
  read_csv()

rosmap_quant_clinical_split <- rosmap_quant_clinical %>%
  mutate(
    across(-c(ID, brain_region, PMI, AOD, CDR, Braak), ~log10(.x + 1))
  ) %>%
  power_inner_join(
    distinct(rosmap_tdp43_classification, file, class_low, class_high),
    by = c("ID" = "file"),
    check = check_specs(
      duplicate_keys_left = "abort",
      duplicate_keys_right = "abort",
      unmatched_keys_left = "warn",
      unmatched_keys_right = "abort"
    )
  ) %>% {
    bind_rows(
      low_threshold = select(., -class_high) %>%
        rename(class = class_low),
      high_threshold = select(., -class_low) %>%
        rename(class = class_high),
      all = select(., -class_low, -class_high) %>%
        mutate(class = "all"),
      .id = "classification"
    )
  } %>%
  group_nest(brain_region, classification, class)

data_dir <- file.path("driad", "chromatin_modifiers")

prepare_task <- \(x, task) {
  task_map <- list(
    all_pooled = c("0", "0", "0", "1", "1", "2", "2"),
    all = as.character(0:6)
  )
  # browser()
  mutate(
    x,
    Label = task_map[[task]][Braak + 1] %>%
      fct_inseq(ordered = TRUE)
  )
}

prediction_tasks <- rosmap_quant_clinical_split %>%
  crossing(
    comparison = c("all", "all_pooled")
  ) %>%
  # Only use conditions that worked best in cross-validation
  filter(
    brain_region == "PCC",
    classification != "all"
  ) %>%
  mutate(
    task = map2(data, comparison, prepare_task)
  ) %>%
  select(-data) %>%
  mutate(
    task_id = paste("prediction_task", brain_region, classification, class, comparison, sep = "_"),
    task_path = file.path(data_dir, paste0(task_id, ".qs"))
  )

prediction_tasks_pairs <- prediction_tasks %>%
  mutate(
    pairs = map(task, DRIAD::preparePairs)
  )

dir.create(data_dir, showWarnings = FALSE)
qsave(
  prediction_tasks_pairs, file.path(data_dir, "prediction_tasks.qs")
)
# prediction_tasks_pairs <- qread("driad/prediction_tasks.qs")

pwalk(
  prediction_tasks_pairs,
  function(task, pairs, task_path, ...) {
    # browser()
    x <- list(task = task, pairs = pairs)
    qsave(
      x,
      task_path
      # compress = "xz"
    )
  }
)

chrom_mod_lrt_res <- syn("syn63854253") %>%
  read_csv()

chrom_mod_sig <- chrom_mod_lrt_res %>%
  drop_na(hgnc_symbol) %>%
  arrange(pvalue) %>%
  filter(padj < .1) %>%
  group_by(name) %>%
  summarize(
    genes = list(head(hgnc_symbol, n = 300)),
    .groups = "drop"
  ) %>%
  mutate(
    gene_set_size = map_int(genes, length)
  )

write_csv(
  chrom_mod_sig,
  file.path(data_dir, "chromatin_modifier_sig.csv.gz")
)

synStoreMany(
  file.path(data_dir, "chromatin_modifier_sig.csv.gz"),
  "syn63901091", forceVersion = FALSE
)
# DRIAD::evalGeneSets(
#   dge_gmt$genes[1],
#   XY = prediction_tasks_pairs$task[[1]],
#   lP = prediction_tasks_pairs$pairs[[1]],
#   nBK = 0,
#   method = "or"
# ) %>%
#   dplyr::select(-Feats, -BK)

gmt_genes <- chrom_mod_sig$genes %>%
  reduce(union)
setdiff(gmt_genes, colnames(rosmap_quant_clinical))

gene_set_sizes <- map_int(chrom_mod_sig$genes, length) %>%
  unique()

gene_universe <- colnames(rosmap_quant_clinical) %>%
  setdiff(
    c("ID", "brain_region", "PMI", "AOD", "CDR",
      "Braak", "Barcode", "Label")
  )

N_BACKGROUND_GENE_SETS <- 1000

set.seed(42)
background_gene_sets <- crossing(
  gene_set_size = gene_set_sizes,
  seq_id = seq_len(N_BACKGROUND_GENE_SETS)
) %>%
  mutate(
    genes = map(
      gene_set_size,
      \(x) sample(gene_universe, x, replace = FALSE)
    ),
    gs_id = paste0("background_", gene_set_size, "_", seq_id)
  )

qsave(
  background_gene_sets,
  file.path(data_dir, "background_gene_sets.qs")
)
# background_gene_sets <- qread(file.path(data_dir, "background_gene_sets.qs"))

compound_gene_sets <- chrom_mod_sig %>%
  mutate(
    gs_id = paste0("compound_", name)
  )

qsave(
  compound_gene_sets,
  file.path(data_dir, "compound_gene_sets.qs")
)
# compound_gene_sets <- qread(file.path(data_dir, "compound_gene_sets.qs"))

combined_gene_sets <- bind_rows(
  compound = compound_gene_sets,
  background = background_gene_sets,
  .id = "gs_type"
)

# Enforce minimum gene set size of 10
# and use only high_threshold for now
combined_gene_sets_selected <- combined_gene_sets %>%
  filter(
    gene_set_size >= 10
  )

prediction_task_gene_set_pairs <- cross_join(
  prediction_tasks_pairs %>%
    select(-c(task, pairs)),
  combined_gene_sets_selected
)

qsave(
  prediction_task_gene_set_pairs,
  file.path(data_dir, "prediction_task_gene_set_pairs.qs")
)
prediction_task_gene_set_pairs <- qread(file.path(data_dir, "prediction_task_gene_set_pairs.qs"))

# Figuring out running time of each task / gene set size combination by
# running a few samples

set.seed(42)
prediction_task_gene_set_pairs_performance_selected <- prediction_task_gene_set_pairs %>%
  group_by(
    task_id, gene_set_size
  ) %>%
  slice_sample(n = 2) %>%
  ungroup() %>%
  mutate(
    gene_sets = map2(genes, gs_id, \(x, y) set_names(list(dummy = x), y))
  )

run_pred_job <- function(task_path, gene_sets, ...) {
  library(tidyverse)
  library(furrr)
  library(DRIAD)
  library(qs)
  library(here)

  plan(sequential)

  message("Reading task...", task_path)
  task <- qread(task_path)
  message("Evaluating gene sets...", length(gene_sets))
  res <- tryCatch({
    DRIAD::evalGeneSets(
      gene_sets,
      XY = task[["task"]],
      lP = task[["pairs"]],
      nBK = 0,
      method = "or"
    ) %>%
      dplyr::select(-Feats, -BK)
  },
  error = function(e) {
    warning(e)
    NULL
  }
  )
  message("Done...")
  res
}


library(batchtools)

reg <- makeRegistry(
  file.dir = here(paste0("registry_chromatin_modifier_run_time", gsub(" ", "_", Sys.time()))),
  seed = 1
)
# reg <- loadRegistry("reg  ")

batchMap(
  fun = run_pred_job,
  task_path = prediction_task_gene_set_pairs_performance_selected$task_path,
  gene_sets = prediction_task_gene_set_pairs_performance_selected$gene_sets,
  reg = reg
)

# run_bk_job(
#   task_id = prediction_tasks_gene_sets[1:5,][["id"]][[1]],
#   background_sets = prediction_tasks_gene_sets[1:5,][["background_sets"]][[1]],
#   background_task_id = prediction_tasks_gene_sets[1:5,][["background_task_id"]][[1]]
# )

job_table <- findJobs(reg = reg) %>%
  # Chunk jobs into a single array job
  mutate(chunk = 1)

submitJobs(
  job_table,
  resources = list(
    memory = "4gb",
    ncpus = 1L,
    partition = "short",
    walltime = 30*60,
    chunks.as.arrayjobs = TRUE
  )
)

running_times_raw <- getJobTable(findDone(reg = reg), reg = reg) %>%
  bind_cols(
    prediction_task_gene_set_pairs_performance_selected
  )

qsave(
  running_times_raw,
  file.path(data_dir, "running_times_raw.qs")
)

running_times <- running_times_raw %>%
  group_by(
    task_id, gene_set_size
  ) %>%
  summarize(
    est_running_time = mean(time.running),
    .groups = "drop"
  )

set.seed(42)
# Planning jobs so that each one should take ~1h
PLANNED_RUNNING_TIME <- 60*60
prediction_task_job_sets <- prediction_task_gene_set_pairs %>%
  power_inner_join(
    running_times,
    by = c("task_id", "gene_set_size"),
    check = check_specs(
      duplicate_keys_right = "warn",
      unmatched_keys_left = "warn",
      unmatched_keys_right = "warn"
    )
  ) %>%
  slice_sample(prop = 1) %>%
  group_by(task_id) %>%
  mutate(
    batch_id = floor(cumsum(unclass(est_running_time)) / PLANNED_RUNNING_TIME)
  ) %>%
  ungroup() %>%
  select(task_id, task_path, gs_id, batch_id, genes) %>%
  group_nest(task_id, task_path, batch_id) %>%
  mutate(
    gene_sets = map(data, \(x) set_names(x$genes, x$gs_id))
  )

qsave(
  prediction_task_job_sets, file.path(data_dir, "prediction_task_job_sets.qs")
)
# prediction_task_job_sets <- qread(file.path(data_dir, "prediction_task_job_sets.qs"))

reg <- makeRegistry(
  file.dir = here(paste0("registry_chromatin_modifiers_", gsub(" ", "_", Sys.time()))),
  seed = 1
)
# reg <- loadRegistry("reg  ")

batchMap(
  fun = run_pred_job,
  task_path = prediction_task_job_sets$task_path,
  gene_sets = prediction_task_job_sets$gene_sets,
  reg = reg
)

job_table <- findJobs(reg = reg) %>%
  # Chunk jobs into a single array job
  mutate(chunk = 1)

submitJobs(
  job_table[findExpired(reg = reg)],
  resources = list(
    memory = "1gb",
    ncpus = 1L,
    partition = "short",
    walltime = 120*60,
    chunks.as.arrayjobs = TRUE
  ),
  reg = reg
)

# Some jobs error out with the following messages. But only 7, no problem.
# Error in (function (.x, .f, ..., .progress = FALSE)  : ℹ In index: 6.
# ℹ With name: background_212_787_prediction_task_3.
# Caused by error in `dplyr::mutate()`:
# ℹ In argument: `AUC = purrr::map_dbl(Data, lpocv, lP, method)`.
# Caused by error in `purrr::map_dbl()`:
# ℹ In index: 1.
# Caused by error in `purrr::map()`:
# ℹ In index: 47.
# ℹ With name: 47.
# Caused by error in `La.svd()`:
# ! error code 1 from Lapack routine 'dgesdd'


# Check if errored jobs contain compound gene sets
prediction_task_job_sets %>%
  slice(
    # findErrors()$job.id
    -findDone(reg = reg)$job.id
  ) %>%
  select(-gene_sets) %>%
  unnest(data) %>%
  filter(str_detect(gs_id, "compound"))
# Nope!

all_res_raw <- reduceResultsDataTable(
  findDone(reg = reg), reg = reg
)

qsave(
  all_res_raw,
  file.path(data_dir, "all_res_raw.qs")
)
all_res_raw <- qread(file.path(data_dir, "all_res_raw.qs"))

# background_res <- qread("driad/background_results_raw.qs")
synStoreMany(
  file.path(data_dir, "all_res_raw.qs"),
  "syn63901091", forceVersion = FALSE
)

all_res_joined <- all_res_raw %>%
  unnest(result) %>%
  rename(gs_id = Set) %>%
  power_inner_join(
    prediction_task_job_sets %>%
      transmute(
        job.id = seq_len(n()),
        task_id
      ),
    by = c("job.id"),
    check = check_specs(
      duplicate_keys_right = "warn",
      unmatched_keys_left = "warn",
      unmatched_keys_right = "warn"
    )
  ) %>%
  power_inner_join(
    prediction_task_gene_set_pairs %>%
      select(-c(genes, task_path)),
    by = c("task_id", "gs_id"),
    check = check_specs(
      duplicate_keys_left = "warn",
      duplicate_keys_right = "warn",
      unmatched_keys_left = "warn",
      unmatched_keys_right = "warn"
    )
  )

fwrite(
  select(all_res_joined, where(negate(is.list))),
  file.path(data_dir, "compound_and_background_res.csv.gz")
)


synStoreMany(
  c(
    file.path(
      data_dir,
      c(
        "compound_and_background_res.csv.gz",
        "prediction_task_gene_set_pairs.qs",
        "background_gene_sets.qs",
        "prediction_tasks.qs"
      )
    )
  ),
  "syn63901091",
  forceVersion = FALSE
)


set.seed(42)
prediction_tasks_all <- rosmap_quant_clinical_split %>%
  crossing(
    comparison = c("all", "all_pooled")
  ) %>%
  # Only use conditions that worked best in cross-validation
  filter(
    brain_region == "PCC",
    classification == "all"
  ) %>%
  mutate(
    task = map2(data, comparison, prepare_task)
  ) %>%
  select(-data) %>%
  mutate(
    task_id = paste("prediction_task", brain_region, classification, class, comparison, sep = "_"),
    task_path = file.path(data_dir, paste0(task_id, ".qs")),
    pairs = map(task, DRIAD::preparePairs)
  )

dir.create(data_dir, showWarnings = FALSE)
qsave(
  prediction_tasks_all, file.path(data_dir, "prediction_tasks_all_only.qs")
)
# prediction_tasks_all <- qread("driad/prediction_tasks_all_only.qs")

pwalk(
  prediction_tasks_all,
  function(task, pairs, task_path, ...) {
    # browser()
    x <- list(task = task, pairs = pairs)
    qsave(
      x,
      task_path
      # compress = "xz"
    )
  }
)


prediction_task_gene_set_pairs_all <- cross_join(
  prediction_tasks_all %>%
    select(-c(task, pairs)),
  combined_gene_sets_selected
)

qsave(
  prediction_task_gene_set_pairs_all,
  file.path(data_dir, "prediction_task_gene_set_pairs_all_only.qs")
)
# prediction_task_gene_set_pairs_all <- qread(file.path(data_dir, "prediction_task_gene_set_pairs_all_only.qs"))

set.seed(42)
prediction_task_gene_set_pairs_all_performance_selected <- prediction_task_gene_set_pairs_all %>%
  group_by(
    task_id, gene_set_size
  ) %>%
  slice_sample(n = 2) %>%
  ungroup() %>%
  mutate(
    gene_sets = map2(genes, gs_id, \(x, y) set_names(list(dummy = x), y))
  )

library(batchtools)

reg <- makeRegistry(
  file.dir = here(paste0("registry_chromatin_modifier_all_run_time", gsub(" ", "_", Sys.time()))),
  seed = 1
)
# reg <- loadRegistry("reg  ")

batchMap(
  fun = run_pred_job,
  task_path = prediction_task_gene_set_pairs_all_performance_selected$task_path,
  gene_sets = prediction_task_gene_set_pairs_all_performance_selected$gene_sets,
  reg = reg
)

# run_bk_job(
#   task_id = prediction_tasks_gene_sets[1:5,][["id"]][[1]],
#   background_sets = prediction_tasks_gene_sets[1:5,][["background_sets"]][[1]],
#   background_task_id = prediction_tasks_gene_sets[1:5,][["background_task_id"]][[1]]
# )

job_table <- findJobs(reg = reg) %>%
  # Chunk jobs into a single array job
  mutate(chunk = 1)

submitJobs(
  job_table,
  resources = list(
    memory = "4gb",
    ncpus = 1L,
    partition = "short",
    walltime = 30*60,
    chunks.as.arrayjobs = TRUE
  )
)

running_times_raw <- getJobTable(findDone(reg = reg), reg = reg) %>%
  bind_cols(
    prediction_task_gene_set_pairs_all_performance_selected
  )

qsave(
  running_times_raw,
  file.path(data_dir, "running_times_raw.qs")
)

running_times <- running_times_raw %>%
  group_by(
    task_id, gene_set_size
  ) %>%
  summarize(
    est_running_time = mean(time.running),
    .groups = "drop"
  )

set.seed(42)
# Planning jobs so that each one should take ~1h
PLANNED_RUNNING_TIME <- 60*60
prediction_task_job_sets_all <- prediction_task_gene_set_pairs_all %>%
  power_inner_join(
    running_times,
    by = c("task_id", "gene_set_size"),
    check = check_specs(
      duplicate_keys_right = "warn",
      unmatched_keys_left = "warn",
      unmatched_keys_right = "warn"
    )
  ) %>%
  slice_sample(prop = 1) %>%
  group_by(task_id) %>%
  mutate(
    batch_id = floor(cumsum(unclass(est_running_time)) / PLANNED_RUNNING_TIME)
  ) %>%
  ungroup() %>%
  select(task_id, task_path, gs_id, batch_id, genes) %>%
  group_nest(task_id, task_path, batch_id) %>%
  mutate(
    gene_sets = map(data, \(x) set_names(x$genes, x$gs_id))
  )

qsave(
  prediction_task_job_sets_all, file.path(data_dir, "prediction_task_job_sets_all.qs")
)
# prediction_task_job_sets_all <- qread(file.path(data_dir, "prediction_task_job_sets_all.qs"))

reg <- makeRegistry(
  file.dir = here(paste0("registry_chromatin_modifiers_all_", gsub(" ", "_", Sys.time()))),
  seed = 1
)
# reg <- loadRegistry("reg  ")

batchMap(
  fun = run_pred_job,
  task_path = prediction_task_job_sets_all$task_path,
  gene_sets = prediction_task_job_sets_all$gene_sets,
  reg = reg
)

job_table <- findJobs(reg = reg) %>%
  # Chunk jobs into a single array job
  mutate(chunk = 1)

submitJobs(
  job_table[findExpired()],
  resources = list(
    memory = "1gb",
    ncpus = 1L,
    partition = "short",
    walltime = 3*60*60,
    chunks.as.arrayjobs = TRUE
  ),
  reg = reg
)



all_res_raw <- reduceResultsDataTable(
  findDone(reg = reg), reg = reg
)

qsave(
  all_res_raw,
  file.path(data_dir, "all_res_raw.qs")
)
all_res_raw <- qread(file.path(data_dir, "all_res_raw.qs"))

# background_res <- qread("driad/background_results_raw.qs")
synStoreMany(
  file.path(data_dir, "all_res_raw.qs"),
  "syn63901091", forceVersion = FALSE
)

all_res_joined <- all_res_raw %>%
  unnest(result) %>%
  rename(gs_id = Set) %>%
  power_inner_join(
    prediction_task_job_sets_all %>%
      transmute(
        job.id = seq_len(n()),
        task_id
      ),
    by = c("job.id"),
    check = check_specs(
      duplicate_keys_right = "warn",
      unmatched_keys_left = "warn",
      unmatched_keys_right = "warn"
    )
  ) %>%
  power_inner_join(
    prediction_task_gene_set_pairs_all %>%
      select(-c(genes, task_path)),
    by = c("task_id", "gs_id"),
    check = check_specs(
      duplicate_keys_left = "warn",
      duplicate_keys_right = "warn",
      unmatched_keys_left = "warn",
      unmatched_keys_right = "warn"
    )
  )

fwrite(
  select(all_res_joined, where(negate(is.list))),
  file.path(data_dir, "compound_and_background_res_all_only.csv.gz")
)


synStoreMany(
  c(
    file.path(
      data_dir,
      c(
        "compound_and_background_res_all_only.csv.gz",
        "prediction_task_gene_set_pairs_all_only.qs",
        "prediction_tasks_all_only.qs"
      )
    )
  ),
  "syn63901091",
  forceVersion = FALSE
)
