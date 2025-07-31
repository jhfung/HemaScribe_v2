#' Stroma classifier
#'
#' This performs annotation of mouse bone marrow stromal cell types using the
#' single cell references from Swann et al. (2024), and also optionally the
#' broad HemaScribe annotation on hematopoietic cells. Currently only supports
#' Seurat object as inputs.
#'
#' @param input A Seurat object
#' @param broad.classify Whether to run "broad annotation" on hematopoietic cells (default: TRUE)
#' @param return.seurat Whether to return Seurat object, otherwise table (default: TRUE)
#' @returns A Seurat object with added annotations, or a table of classifier results
#' @export
stroma_classify <- function(input, broad.classify = TRUE, return.seurat = TRUE) {
  if (inherits(input, "Seurat")) {
    exp_query <- input[["RNA"]]$counts
    metadata_query <- input@meta.data
  } else {
    rlang::abort("Input must be a Seurat object.")
  }

  if (startsWith(rownames(exp_query)[1], "ENSMUSG")) {
    gene_symbols <- AnnotationDbi::mapIds(
      org.Mm.eg.db::org.Mm.eg.db,
      keys = rownames(exp_query),
      column = "SYMBOL",
      keytype = "ENSEMBL"
    )
    gene_symbols <- unique(gene_symbols[!is.na(gene_symbols)])

    exp_query <- exp_query[names(gene_symbols),]
    rownames(exp_query) <- gene_symbols
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
  input <- Seurat::AddMetaData(input, metadata=pred)

  if (broad.classify) {
    hem.cats <- c("EryP1", "EryP2", "EryP3", "EryP4", "EryP5", "HSPC", "Lymph", "MkP", "Mono", "Neut1", "Neut2")
    input.hem <- subset(input, subset=(input$stroma.annot %in% hem.cats))
    if (!("data" %in% SeuratObject::Layers(input.hem))) {
      input.hem <- Seurat::NormalizeData(input.hem, scale.factor=10000, verbose=FALSE)
    }

    pred.broad <- broad_classify(input.hem, return.full=TRUE)
    pred.broad <- data.frame(pred.broad)
    pred.broad$max.score <- apply(pred.broad[,grep("^orig\\.results.*scores", colnames(pred.broad), value=TRUE)], 1, max)
    pred.broad <- pred.broad[,c("pruned.labels", "max.score")]
    colnames(pred.broad) <- c("broad.annot", "broad.annot_score")
    pred.broad <- pred.broad[pred.broad$broad.annot != "",]

    pred <- merge(pred, pred.broad, by="row.names", all.x=TRUE, all.y=FALSE)
    rownames(pred) <- pred$Row.names
    pred$Row.names <- NULL

    input <- Seurat::AddMetaData(input, metadata=pred.broad)
  }

  if (return.seurat) {
    return(input)
  } else {
    return(pred)
  }
}

