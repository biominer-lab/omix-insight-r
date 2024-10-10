#' It's a wrapper of limma and DESeq2 for differential expression analysis.
#' Jingcheng Yang <yjcyxky@163.com>
#' 2024-10-08

#' Run differential expression analysis using limma
#'
#' @param exp_data_file Path to the expression data file: the columns are samples, the rows are genes; the first column is the gene/protein id (id), the second column is the gene symbol (gene_symbol)
#' @param sample_info_file Path to the sample information file: the file should have at least two columns, sample_id and group and the sample_id should be the same as the column names of the expression data file
#' @param test_group_name Name of the test group: the name of the test group should be one of the group in the sample information file
#' @param control_group_name Name of the control group: the name of the control group should be one of the group in the sample information file
#' @param output_file Path to the output file: the output file will be a tab-separated file with the following columns: id, gene_symbol, and the rest of the columns are the log-fold changes, t-statistics, p-values, and adjusted p-values
#' @param enableLog Whether to log-transform the expression data: if TRUE, the expression data will be log-transformed
#' 
#' @return None
#' 
#' @importFrom dplyr select
#' 
#' @examples
#' run_deg_with_limma("exp_data.txt", "sample_info.txt", "test_group", "control_group", "output.txt", TRUE)
#' 
#' @export
run_deg_with_limma <- function(exp_data_file, sample_info_file, test_group_name, control_group_name, output_file, enableLog = FALSE) {
    # check if the input files exist
    if (!file.exists(exp_data_file)) {
        stop("Expression data file does not exist")
    }
    if (!file.exists(sample_info_file)) {
        stop("Sample information file does not exist")
    }

    # Check if the input columns exist
    if (!all(c("sample_id", "group") %in% colnames(sample_info_file))) {
        stop("Sample information file does not have the required columns: sample_id and group")
    }

    # Check if the test group and control group exist in the sample information file
    if (!test_group_name %in% sample_info_file$group) {
        stop("Test group name does not exist in the sample information file")
    }
    if (!control_group_name %in% sample_info_file$group) {
        stop("Control group name does not exist in the sample information file")
    }

    # Check if the expression data file has the required columns
    if (!all(c("id", "gene_symbol") %in% colnames(exp_data_file))) {
        stop("Expression data file does not have the required columns: id and gene_symbol")
    }

    # Check if the expression data file has unique ids and these sample ids are consistent with the sample information file
    if (any(duplicated(exp_data_file$id))) {
        stop("Expression data file has duplicate ids")
    }

    if (!all(exp_data_file$id %in% sample_info_file$sample_id)) {
        stop("Expression data file has sample ids that do not exist in the sample information file")
    }

    # Read expression data and sample information
    expression_data <- read.table(exp_data_file, header = TRUE, sep = "\t")
    row.names(expression_data) <- expression_data$ID
    anno_data <- expression_data %>% select(c("gene_symbol", "id"))
    expression_data <- expression_data %>% select(-c("gene_symbol", "id"))
    if (enableLog) {
        expression_data <- log2(expression_data)
    }
    sample_info <- read.table(sample_info_file, header = TRUE, sep = "\t")

    # Ensure the sample order is consistent
    sample_info <- sample_info[match(colnames(expression_data), sample_info$sample_id), ]

    # Create design matrix
    group <- factor(sample_info$group)
    print("group:")
    print(group)
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)

    # Create contrast matrix (based on the input test group and control group names)
    contrast_formula <- paste(test_group_name, "-", control_group_name, sep = "")
    print("contrast_formula:")
    print(contrast_formula)
    contrast.matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

    # Perform differential expression analysis using limma
    fit <- limma::lmFit(expression_data, design)
    fit2 <- limma::contrasts.fit(fit, contrast.matrix)
    fit2 <- limma::eBayes(fit2)

    # Get differential expression results
    results <- limma::topTable(fit2, adjust = "fdr", number = Inf)
    results <- merge(results, anno_data, by.x = "row.names", by.y = "id", all.x = TRUE)

    results$id <- results$Row.names
    results <- results %>% select(-Row.names)

    # Reorder the columns
    results <- results %>% select(id, gene_symbol, everything())

    # Save results
    write.table(results, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}


#' [Not Ready] Run differential expression analysis using DESeq2
#'
#' @param counts_file Path to the counts data file: the columns are samples, the rows are genes; the first column is the gene/protein id (id), the second column is the gene symbol (gene_symbol), the remaining columns are the counts
#' @param sample_info_file Path to the sample information file: the file should have at least two columns, sample_id and group and the sample_id should be the same as the column names of the counts data file
#' @param test_group_name Name of the test group: the name of the test group should be one of the group in the sample information file
#' @param control_group_name Name of the control group: the name of the control group should be one of the group in the sample information file
#' @param output_file Path to the output file: the output file will be a tab-separated file with the following columns: id, gene_symbol, and the rest of the columns are the log-fold changes, t-statistics, p-values, and adjusted p-values
#' 
#' @return None
#' 
#' @importFrom dplyr select
#' 
#' @examples
#' run_deg_with_deseq2("counts.txt", "sample_info.txt", "test_group", "control_group", "output.txt")
#' 
#' @export
run_deg_with_deseq2 <- function(counts_file, sample_info_file, test_group_name, control_group_name, output_file) {
    # Check if the input files exist
    if (!file.exists(counts_file)) {
        stop("Counts data file does not exist")
    }
    if (!file.exists(sample_info_file)) {
        stop("Sample information file does not exist")
    }

    # Check if the input columns exist
    if (!all(c("sample_id", "group") %in% colnames(sample_info_file))) {
        stop("Sample information file does not have the required columns: sample_id and group")
    }

    # Check if the counts data file has the required columns
    if (!all(c("id", "gene_symbol") %in% colnames(counts_file))) {
        stop("Counts data file does not have the required columns: id and gene_symbol")
    }

    # Check if the counts data file has unique ids and these sample ids are consistent with the sample information file
    if (any(duplicated(counts_file$id))) {
        stop("Counts data file has duplicate ids")
    }

    if (!all(counts_file$id %in% sample_info_file$sample_id)) {
        stop("Counts data file has sample ids that do not exist in the sample information file")
    }

    # Read counts data and sample information
    counts_data <- read.table(counts_file, header = TRUE, sep = "\t")
    row.names(counts_data) <- counts_data$id
    anno_data <- counts_data %>% select(c("gene_symbol", "id"))
    counts_data <- counts_data %>% select(-c("gene_symbol", "id"))

    # Ensure the sample order is consistent
    sample_info <- sample_info[match(colnames(counts_data), sample_info$sample_id), ]

    # Create DESeq2 object
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_data, colData = sample_info, design = ~ group)

    # Run DESeq2
    dds <- DESeq2::DESeq(dds)

    # Get differential expression results
    results <- DESeq2::results(dds, contrast = c("group", test_group_name, control_group_name))
    results <- merge(results, anno_data, by.x = "row.names", by.y = "id", all.x = TRUE)

    results$id <- results$Row.names
    results <- results %>% select(-Row.names)

    # Reorder the columns
    results <- results %>% select(id, gene_symbol, everything())

    # Save results
    write.table(results, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}
