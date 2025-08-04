test_that("hematopoietic_score works (Seurat, gene symbols)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_1.rds"))
  hem.scores <- hematopoietic_score(aa003)
  expect_equal(hem.scores, aa003@meta.data["hematopoietic.score"])
})

test_that("hematopoietic_score works (Seurat, gene IDs)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_2.rds"))
  hem.scores <- hematopoietic_score(aa003)
  expect_equal(hem.scores, aa003@meta.data["hematopoietic.score"])
})

test_that("hematopoietic_score works (SingleCellExperiment, gene IDs)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_3.rds"))
  hem.scores <- hematopoietic_score(aa003)
  expect_equal(hem.scores, data.frame(SingleCellExperiment::colData(aa003)["hematopoietic.score"]))
})

test_that("hematopoietic_score works (SingleCellExperiment, gene symbols)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_4.rds"))
  hem.scores <- hematopoietic_score(aa003)
  expect_equal(hem.scores, aa003$hematopoietic.score)
})

test_that("broad_classify works (Seurat, gene symbols)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_1.rds"))
  broad.annot <- broad_classify(aa003, return.full = TRUE)
  broad.annot <- data.frame(broad.annot["pruned.labels"])
  colnames(broad.annot) <- "broad.annot"
  expect_equal(broad.annot, aa003@meta.data["broad.annot"])
})

test_that("broad_classify works (Seurat, gene IDs)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_2.rds"))
  broad.annot <- broad_classify(aa003, return.full = TRUE)
  broad.annot <- data.frame(broad.annot["pruned.labels"])
  colnames(broad.annot) <- "broad.annot"
  expect_equal(broad.annot, aa003@meta.data["broad.annot"])
})

test_that("broad_classify works (SingleCellExperiment, gene IDs)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_3.rds"))
  broad.annot <- broad_classify(aa003, return.full = TRUE)
  broad.annot <- data.frame(broad.annot["pruned.labels"])
  colnames(broad.annot) <- "broad.annot"
  expect_equal(broad.annot, data.frame(SingleCellExperiment::colData(aa003)["broad.annot"]))
})

test_that("broad_classify works (SingleCellExperiment, gene symbols)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_4.rds"))
  broad.annot <- broad_classify(aa003, return.full = TRUE)
  broad.annot <- data.frame(broad.annot["pruned.labels"])
  colnames(broad.annot) <- "broad.annot"
  expect_equal(broad.annot, data.frame(SingleCellExperiment::colData(aa003)["broad.annot"]))
})

test_that("fine_classify works (Seurat, gene symbols)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_1.rds"))
  fine.annot <- fine_classify(aa003, return.full = TRUE)
  fine.annot <- data.frame(fine.annot["predicted.id"])
  colnames(fine.annot) <- "fine.annot"
  expect_equal(fine.annot, aa003@meta.data["fine.annot"])
})

test_that("fine_classify works (Seurat, gene IDs)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_2.rds"))
  fine.annot <- fine_classify(aa003, return.full = TRUE)
  fine.annot <- data.frame(fine.annot["predicted.id"])
  colnames(fine.annot) <- "fine.annot"
  expect_equal(fine.annot, aa003@meta.data["fine.annot"])
})

test_that("fine_classify works (SingleCellExperiment, gene IDs)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_3.rds"))
  fine.annot <- fine_classify(aa003, return.full = TRUE)
  fine.annot <- data.frame(fine.annot["predicted.id"])
  colnames(fine.annot) <- "fine.annot"
  expect_equal(fine.annot, data.frame(SingleCellExperiment::colData(aa003)["fine.annot"]))
})

test_that("fine_classify works (SingleCellExperiment, gene symbols)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_4.rds"))
  fine.annot <- fine_classify(aa003, return.full = TRUE)
  fine.annot <- data.frame(fine.annot["predicted.id"])
  colnames(fine.annot) <- "fine.annot"
  expect_equal(fine.annot, data.frame(SingleCellExperiment::colData(aa003)["fine.annot"]))
})

test_that("HemaScribe works (Seurat, gene symbols)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_1.rds"))
  suppressMessages({
    results <- HemaScribe(aa003, return.full = TRUE)
  })
  expect_equal(results$combined["combined.annot"], aa003@meta.data["combined.annot"])
  expect_equal(results$combined["HSPC.annot"], aa003@meta.data["HSPC.annot"])
})

test_that("HemaScribe works (Seurat, gene IDs)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_2.rds"))
  suppressMessages({
    results <- HemaScribe(aa003, return.full = TRUE)
  })
  expect_equal(results$combined["combined.annot"], aa003@meta.data["combined.annot"])
  expect_equal(results$combined["HSPC.annot"], aa003@meta.data["HSPC.annot"])
})

test_that("HemaScribe works (SingleCellExperiment, gene IDs)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_3.rds"))
  suppressMessages({
    results <- HemaScribe(aa003, return.full = TRUE)
  })
  expected <- data.frame(SingleCellExperiment::colData(aa003))
  expected <- expected[order(rownames(expected)),]
  expect_equal(results$combined["combined.annot"], expected["combined.annot"])
  expect_equal(results$combined["HSPC.annot"], expected["HSPC.annot"])
})

test_that("HemaScribe works (SingleCellExperiment, gene symbols)", {
  aa003 <- readRDS(test_path("fixtures", "AA003_4.rds"))
  suppressMessages({
    results <- HemaScribe(aa003, return.full = TRUE)
  })
  expected <- data.frame(SingleCellExperiment::colData(aa003))
  expected <- expected[order(rownames(expected)),]
  expect_equal(results$combined["combined.annot"], expected["combined.annot"])
  expect_equal(results$combined["HSPC.annot"], expected["HSPC.annot"])
})

test_that("fine_classify throws warning with too few HSPCs", {
  # expect_error(...)
})


