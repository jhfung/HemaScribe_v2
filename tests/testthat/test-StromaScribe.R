test_that("stroma_classify works (Seurat, gene symbols)", {
  pbsem <- readRDS(test_path("fixtures", "PBSEM_1.rds"))
  suppressMessages({
    results <- stroma_classify(pbsem, return.full=TRUE)
  })
  stroma.sim <- mean(results["stroma.annot"] == pbsem@meta.data["stroma.annot"])
  # message(stroma.sim)
  expect_true(stroma.sim > 0.95)
  broad.sim <- mean(results["broad.annot"] == pbsem@meta.data["broad.annot"], na.rm=TRUE)
  # message(broad.sim)
  expect_true(broad.sim > 0.95)
  print(sessionInfo())
})

test_that("stroma_classify works (Seurat, gene IDs)", {
  pbsem <- readRDS(test_path("fixtures", "PBSEM_2.rds"))
  suppressMessages({
    results <- stroma_classify(pbsem, return.full=TRUE)
  })
  stroma.sim <- mean(results["stroma.annot"] == pbsem@meta.data["stroma.annot"])
  # message(stroma.sim)
  expect_true(stroma.sim > 0.95)
  broad.sim <- mean(results["broad.annot"] == pbsem@meta.data["broad.annot"], na.rm=TRUE)
  # message(broad.sim)
  expect_true(broad.sim > 0.95)
})

test_that("stroma_classify works (SingleCellExperiment, gene IDs)", {
  pbsem <- readRDS(test_path("fixtures", "PBSEM_3.rds"))
  suppressMessages({
    results <- stroma_classify(pbsem, return.full=TRUE)
  })
  expected <- data.frame(SingleCellExperiment::colData(pbsem))
  expected <- expected[order(rownames(expected)),]
  stroma.sim <- mean(results["stroma.annot"] == expected["stroma.annot"])
  # message(stroma.sim)
  expect_true(stroma.sim > 0.95)
  broad.sim <- mean(results["broad.annot"] == expected["broad.annot"], na.rm=TRUE)
  # message(broad.sim)
  expect_true(broad.sim > 0.95)
})

test_that("stroma_classify works (SingleCellExperiment, gene symbols)", {
  pbsem <- readRDS(test_path("fixtures", "PBSEM_4.rds"))
  suppressMessages({
    results <- stroma_classify(pbsem, return.full=TRUE)
  })
  expected <- data.frame(SingleCellExperiment::colData(pbsem))
  expected <- expected[order(rownames(expected)),]
  stroma.sim <- mean(results["stroma.annot"] == expected["stroma.annot"])
  # message(stroma.sim)
  expect_true(stroma.sim > 0.95)
  broad.sim <- mean(results["broad.annot"] == expected["broad.annot"], na.rm=TRUE)
  # message(broad.sim)
  expect_true(broad.sim > 0.95)
})
