
#' Function for mapping query cells to a HemaScape reference using Symphony algorithm
#'
#' This function maps query cells to a pre-built HemaScape reference of the DensityPath trajectory using the Symphony algorithm.
#' It performs EE projection, and predicts branch assignments, density clusters, refined density clusters along each trajectory branch, and pseudotime.
#'
#' @param input A Seurat obejct with count data.
#' @param vars Query batch variable(s) to integrate over (column names in metadata of the Seurat object).
#' @param sigma Fuzziness parameter for soft clustering (sigma = 1 is hard clustering; default: 0.1).
#'
#' @importFrom raster raster
#' @importFrom FNN knn.reg
#' @importFrom gdistance costDistance
#' @importFrom symphony mapQuery knnPredict
#'
#' @return The Seurat object containing mapped results, including predicted branch assignments, density clusters, density cluster divided by branch, and pseudotime.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(HemaScribe)
#' library(GEOquery)
#' library(Seurat)
#' library(fs)
#'
#' ## Step 1: Download and load the LK dataset from GEO accession GSM7712015
#' # Download and rename files associated with the GEO accession
#' geo_accession <- "GSM7712015"
#' supp_files <- getGEOSuppFiles(geo_accession)
#' files <- dir_ls(geo_accession, regexp = "GSM7712015_AA003_")
#' file_move(files, path(geo_accession, path_file(files) %>% str_replace("GSM7712015_AA003_", "")))
#'
#' # Read count data and create a Seurat object
#' count_mat <- Seurat::Read10X(geo_accession)
#' Seurat_obj <- CreateSeuratObject(counts = count_mat)
#'
#' ## Step 2: Specify meta column indicating the sequencing batch (default using "orig.ident" column)
#' Seurat_obj$Seq_batch_ID <- rep(1, ncol(Seurat_obj))
#'
#' ## Step 3: Perform HemaScape mapping
#' # Run the HemaScape function on the Seurat object
#' Seurat_obj <- HemaScape(
#'   input = Seurat_obj,
#'   vars = 'Seq_batch_ID'
#' )
#'
#' ## Step 4: Inspect the results
#' # View predicted branch assignments on predicted EE embedding
#' DimPlot(Seurat_obj,reduction = 'EE', group.by='branch_pred')
#'
#' # View predicted density clusters on predicted EE embedding
#' DimPlot(Seurat_obj,reduction = 'EE', group.by='density_cluster_pred')
#'
#' # View predicted refined density clusters along each trajectory branch on predicted EE embedding
#' DimPlot(Seurat_obj,reduction = 'EE', group.by='branch_segment_clusters_pred')
#'
#' # View predicted pseudotime on predicted EE embedding
#' FeaturePlot(Seurat_obj, features = "pseudotime_pred", reduction = "EE")
#' }
#'
#' @seealso \code{\link{HemaScape}} for more details on the function.
HemaScape <- function(input, vars='orig.ident', sigma = 0.1){
  # Access the internal reference dataset
  # reference <- get("reference", envir = asNamespace("HemaScribe"))
  reference <- HemaScribeData$mapping_reference

  if (inherits(input, "Seurat")) {
    exp_query <- input[["RNA"]]$counts
    metadata_query <- input@meta.data
  }else{
    stop("Only Seurat formats are supported.")
  }

  if (all(startsWith(rownames(exp_query), "ENSMUSG"))) {
    message("Row names are Ensembl IDs. Map to Gene Symbols.")
    gene_symbols <- AnnotationDbi::mapIds(
      org.Mm.eg.db::org.Mm.eg.db,
      keys = rownames(exp_query),  # Row names of the matrix
      column = "SYMBOL",           # Column to map to (gene symbols)
      keytype = "ENSEMBL"          # Key type (Ensembl IDs)
    )

    # Remove rows with NA or duplicated gene symbols
    gene_symbols <- gene_symbols[!is.na(gene_symbols)]
    gene_symbols <- gene_symbols[!duplicated(gene_symbols)]

    # Subset the matrix to include only rows with valid gene symbols
    exp_query <- exp_query[names(gene_symbols), ]

    # Update row names with gene symbols
    rownames(exp_query) <- as.vector(gene_symbols)
  }

  # Remove cell cycle and histone family genes
  Remin_Gene <- setdiff(rownames(exp_query), rownames(exp_query)[stringr::str_detect(rownames(exp_query), 'Rpl|Rps|Mrpl|Mrps|mt-')])
  Remin_Gene <- setdiff(Remin_Gene, reference$cellcycle_genes_removed)
  all_gene <- rownames(exp_query)
  Histone_genes <- all_gene[grep('Hist', all_gene)]
  Remin_Gene <- setdiff(Remin_Gene, Histone_genes)
  exp_query <- exp_query[Remin_Gene,]

  Z_ref = reference$Z_corr

  query_donor = symphony::mapQuery(exp_query, metadata_query, reference, vars, verbose=FALSE,
                         do_normalize=TRUE, do_umap=FALSE, sigma)

  rlang::inform("Predict EE coordinates")
  EE1 = FNN::knn.reg(t(Z_ref), t(query_donor$Z), reference$meta_data$EE1, k = 10)$pred
  EE2 = FNN::knn.reg(t(Z_ref), t(query_donor$Z), reference$meta_data$EE2, k = 10)$pred

  EE <- cbind(EE1, EE2)
  colnames(EE) <- c('EE_1', 'EE_2')
  rownames(EE) <- rownames(query_donor$meta_data)

  input[["EE"]] <- Seurat::CreateDimReducObject(
    embeddings = EE,
    key = "EE_",
    assay = Seurat::DefaultAssay(input)
  )

  rlang::inform("Predict branch assignment")
  set.seed(0)
  query_donor = symphony::knnPredict(query_donor, reference,
                           train_labels = reference$meta_data$branch.labels, k = 30, confidence = TRUE)
  query_donor$meta_data$branch_pred <- query_donor$meta_data$cell_type_pred_knn
  query_donor$meta_data$branch_pred_prob <- query_donor$meta_data$cell_type_pred_knn_prob

  rlang::inform("Predict density cluster")
  set.seed(0)
  query_donor = symphony::knnPredict(query_donor, reference,
                           train_labels = reference$meta_data$density_cluster, k = 30, confidence = TRUE)
  query_donor$meta_data$density_cluster_pred <- query_donor$meta_data$cell_type_pred_knn
  query_donor$meta_data$density_cluster_pred_prob <- query_donor$meta_data$cell_type_pred_knn_prob

  rlang::inform("Predict refined density clusters along each trajectory branch")
  set.seed(0)
  query_donor = symphony::knnPredict(query_donor, reference,
                           train_labels = reference$meta_data$branch_segment_cluster, k = 30, confidence = TRUE)
  query_donor$meta_data$branch_segment_clusters_pred <- query_donor$meta_data$cell_type_pred_knn
  query_donor$meta_data$branch_segment_clusters_pred_prob <- query_donor$meta_data$cell_type_pred_knn_prob

  query_donor$meta_data$cell_type_pred_knn <- NULL
  query_donor$meta_data$cell_type_pred_knn_prob <- NULL

  rlang::inform("Predict pseudotime")
  pseudotime_pred <- double()
  for(j in 1:nrow(EE)){
    point <- EE[j,]
    dis <- gdistance::costDistance(reference$T, point, reference$traj_globalpoints[,1:2])
    pseudotime_pred[j] <- reference$traj_globalpoints_dis[min(which(dis==min(dis)))]
  }
  pseudotime_pred <- pseudotime_pred/max(reference$meta_data$pseudotime_orig)
  query_donor$meta_data$pseudotime_pred <- pseudotime_pred

  rlang::inform("Returning final mapping results")
  input@meta.data <- query_donor$meta_data
  input
}
