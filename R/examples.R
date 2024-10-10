#' Examples
#'

#' plot-boxplot-two-genes-two-groups
#' d <- query_data(exp_file = "./examples/gene_expression.tsv", sample_info_file = "./examples/sample_info.tsv", which_gene_symbols = c("BDNF", "TP53"), which_groups = c("Female_Control", "Female_MECFS"))
#' boxplotly(d, method = "t.test", log_scale = FALSE, enable_label = FALSE, enable_log2fc = TRUE)


#' plot-boxplot-one-gene-two-groups
#' d <- query_data(exp_file = "./examples/gene_expression.tsv", sample_info_file = "./examples/sample_info.tsv", which_gene_symbols = c("BDNF"), which_groups = c("Female_Control", "Female_MECFS"))
#' boxplotly(d, method = "t.test", log_scale = FALSE, enable_label = FALSE, enable_log2fc = TRUE)


#' plot-boxplot-multiple-genes-two-groups
#' d <- query_data(exp_file = "./examples/gene_expression.tsv", sample_info_file = "./examples/sample_info.tsv", which_gene_symbols = c("BDNF", "TP53", "EGFR"), which_groups = c("Female_Control", "Female_MECFS"))
#' boxplotly(d, method = "t.test", log_scale = FALSE, enable_label = FALSE, enable_log2fc = TRUE)


#' @name plot-boxplot-multiple-genes-multiple-groups
NULL
