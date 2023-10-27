library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  counts <- read_tsv(counts_csv)
  colData <- read_csv(metafile_csv)
  rowData <- counts["gene"]
  counts <- as.data.frame(
    dplyr::select(counts, -gene)
  )
  rownames(counts) <- rowData$gene
  se <- SummarizedExperiment(assays=list(counts=as.matrix(counts)),colData=colData,rowData =rowData)
  se <- se[,se$timepoint %in% c("vP0","vAd")]
  se$timepoint <- relevel(factor(se$timepoint),ref="vP0")
  return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  deseq <- DESeqDataSet(se, design = design)
  deseq <- DESeq (deseq)
  result <- results(deseq) 
  return(list(deseq,data.frame(result)))
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`.
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  deseq2_res$volc_plot_status[deseq2_res$padj<padj_threshold & deseq2_res$log2FoldChange>0] <- 'UP'
  deseq2_res$volc_plot_status[deseq2_res$padj<padj_threshold & deseq2_res$log2FoldChange<0] <- 'DOWN'
  deseq2_res$volc_plot_status[is.na(deseq2_res$volc_plot_status)] <- "NS"
  return(tibble(deseq2_res))
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  pval_plot <- ggplot(data = labeled_results, aes(x = pvalue)) + geom_histogram()
  return(pval_plot)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  sig = na.omit(res[res['padj']>padj_threshold,])
  log2FC_plot <- ggplot(data = sig, aes(x = log2FoldChange)) + geom_histogram()
  return(log2FC_plot)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  labeled_results <- arrange(labeled_results,padj)
  genes <- labeled_results[1:num_genes,]$genes
  rows_to_keep <- rowData(dds_obj)$gene %in% genes
  dds_subset <- dds[rows_to_keep, ]
  final = tibble(data.frame(assays(dds_subset)$counts))
  final1 = pivot_longer(final,cols=colnames(final))
  sp <- ggplot(data = final1, aes(x = name, y = value)) +geom_point() +scale_x_discrete(name = "Names")
  return(sp)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  vp <- ggplot(data=labeled_results, aes(x=log2FoldChange, y=-log10(padj),col=volc_plot_status)) + geom_point() + theme_minimal()
  return(vp)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  meta <- read_table(id2gene_path,col_names=FALSE)
  colnames(meta) <- c("GeneId","GeneSymbol")
  labeled_results$GeneId = labeled_results$genes
  final <- merge(labeled_results,meta,by='GeneId')
  final <- drop_na(final)
  #print(final)
  res <- c(final$log2FoldChange)
  names(res) <- t(final$GeneSymbol)[1:length(final$GeneSymbol)]
  return(sort(res,decreasing=TRUE))
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  pathways <- gmtPathways(gmt_file_path)
  fgseaRes <- fgsea(pathways, rnk_list, minSize=15, maxSize=500)
  return(tibble(fgseaRes))
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  fgsea_results <- fgsea_results[order(-fgsea_results$NES),]
  top_10_positive <- head(fgsea_results[fgsea_results$NES > 0,], num_paths)
  top_10_negative <- head(fgsea_results[fgsea_results$NES < 0,], num_paths)
  top_20 <- rbind(top_10_positive, top_10_negative)
  top_20$NES_sign <- ifelse(top_20$NES > 0, "Positive", "Negative")
  ggplot(top_20, aes(x = reorder(pathway, NES), y = NES, fill = NES_sign)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    labs(title = "Top 10 Pathways by NES Score",
         x = "Pathway",
         y = "NES Score") +
    scale_fill_manual(values = c("Positive" = "blue", "Negative" = "red"))
}

