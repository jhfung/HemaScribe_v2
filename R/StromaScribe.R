#' Stroma classifier
#'
#' This performs annotation of mouse bone marrow stromal cell types using the
#' single cell references from Swann et al. (2024), and also optionally the
#' broad HemaScribe annotation on hematopoietic cells.
#'
#' @param input A SingleCellExperiment or Seurat object
#' @param broad.classify Whether to run "broad annotation" on hematopoietic cells (default: TRUE)
#' @param return.full Whether to return full predictions directly (default: FALSE)
#' @returns An object of the same type as input with added annotations, or a dataframe with classifier results
#' @export
stroma_classify <- function(input, broad.classify = TRUE, return.full = FALSE) {
  if (inherits(input, "Seurat")) {
    exp_query <- input[["RNA"]]$counts
    metadata_query <- input@meta.data
  } else if (inherits(input, "SingleCellExperiment")) {
    exp_query <- SingleCellExperiment::counts(input)
    metadata_query <- data.frame(SingleCellExperiment::colData(input))
  } else {
    rlang::abort("Only SingleCellExperiment and Seurat formats are supported.")
  }

  if (!startsWith(rownames(exp_query)[1], "ENSMUSG")) {
    suppressMessages({
      gene_ids <- AnnotationDbi::mapIds(
        org.Mm.eg.db::org.Mm.eg.db,
        keys = rownames(exp_query),
        column = "ENSEMBL",
        keytype = "SYMBOL"
      )
    })
    gene_ids <- gene_ids[!is.na(gene_ids)]
    gene_ids <- gene_ids[!duplicated(unlist(gene_ids))]
    exp_query <- exp_query[names(gene_ids),]
    rownames(exp_query) <- gene_ids
  }

  ref.stroma <- HemaScribeData$ref.stroma
  query <- symphony::mapQuery(
    exp_query,
    metadata_query,
    ref.stroma,
    do_normalize=TRUE,
    do_umap=FALSE
  )

  rlang::inform("Predicting cell type annotations")
  mapped <- symphony::knnPredict(query, ref.stroma, ref.stroma$meta_data$Type, save_as="stroma.annot")
  pred <- mapped$meta_data[,c("stroma.annot", "stroma.annot_prob")]

  if (broad.classify) {
    hem.cats <- c("EryP1", "EryP2", "EryP3", "EryP4", "EryP5", "HSPC", "Lymph", "MkP", "Mono", "Neut1", "Neut2")
    input.hem <- input[,which(pred$stroma.annot %in% hem.cats)]

    pred.broad <- broad_classify(input.hem, return.full=TRUE)
    pred.broad <- data.frame(pred.broad)
    pred.broad$max.score <- apply(pred.broad[,grep("^orig\\.results.*scores", colnames(pred.broad), value=TRUE)], 1, max)
    pred.broad <- pred.broad[,c("pruned.labels", "max.score")]
    colnames(pred.broad) <- c("broad.annot", "broad.annot_score")
    pred.broad <- pred.broad[pred.broad$broad.annot != "",]

    pred <- merge(pred, pred.broad, by="row.names", all.x=TRUE, all.y=FALSE)
    rownames(pred) <- pred$Row.names
    pred$Row.names <- NULL

    if (return.full) {
      return(pred)
    }

    if (inherits(input, "Seurat")) {
      output <- Seurat::AddMetaData(input, metadata = pred[c("stroma.annot", "broad.annot")])
      return(output)
    } else if (inherits(input, "SingleCellExperiment")) {
      output <- input
      pred <- pred[match(colnames(output), rownames(pred)),]
      colData(output) <- cbind(SingleCellExperiment::colData(output), pred)
      return(output)
    }
  }
}

