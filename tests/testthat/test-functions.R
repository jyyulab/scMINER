context("Unit testing non-graphing functions in scMINER's 'functions.R' file\n")
# Testing set from: https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_protein_v3
library("scMINER")
library(tools)
library('rhdf5')

# Reading input and preprocessing data
d.curr <- readscRNAseqData(file="../../testing_set/",is.10x = T, CreateSparseEset = F, add.meta=F)
e.curr <- CreateSparseEset(data=d.curr$raw.data,feature.data = d.curr$feature.data, add.meta = T)

cutoffs <- list(
  nCell_cutoff = 3,
  umi_cf_lo = 2340,
  umi_cf_hi = 57773,
  nGene_cf = 378,
  ERCC_cf = 0,
  mito_cf = 0.146
)

e.sel.curr <- preMICA.filtering(SparseEset = e.curr, cutoffs = cutoffs)
norm = 1e6

# Normalization and transformation
exp.norm.curr <- sweep(exprs(e.sel.curr), 2, norm/unname(Matrix::colSums(exprs(e.sel.curr))), '*')
exp.log2.curr <- log(exp.norm.curr+1,base=2)
eset.log2.curr <- CreateSparseEset(data=exp.log2.curr, meta.data = pData(e.sel.curr), feature.data = fData(e.sel.curr), add.meta = F)

# Reading output from MICA
eset.curr <- readMICAoutput(eset = eset.log2, load_ClusterRes = TRUE, output_file = "../../testing_set/answerkey/mica_output/test_k4_tsne_ClusterMem.txt")

# Getting activity from SJARACNe run
acs.curr <- GetActivityFromSJARACNe(SJARACNe_output_path = "../../testing_set/answerkey/sjaracne_output",
                                    SJARACNe_input_eset = eset.curr,
                                    group_name = "ClusterRes",
                                    save_network_file = FALSE,
                                    save_path = NULL)

# Load the answer-key from a pre-made R environment
load('../../testing_set/answerkey/testEnv.RData')
load('../../testing_set/answerkey/testEnv2.RData')

# Tests the contents of the matrix read from the input data
test_that("readscRNAseqData works", {
  expect_equal(d, d.curr)
})


# Tests the contents of sparse matrix
test_that("CreateSparseEset works", {
  expect_equal(e, e.curr)
})


test_that("preMICA.filtering works", {
  expect_equal(eset.sel, e.sel.curr)
})


test_that("generateMICAinput works on txt files", {
  generateMICAinput(d = exp.log2.curr, filename="test_input.txt")
  match <- 0
  hash1 <- md5sum("test_input.txt")
  hash2 <- md5sum("../../testing_set/answerkey/mica_1k_input.txt")
  if (hash1 == hash2) {
    match <- 1
  }
  expect_equal(match, 1)
  file.remove("test_input.txt")
})


# Currently test below does not work because h5 files have different checksum values,
# must find a reasonable way to compare them

test_that("generateMICAinput works on h5 files", {
  generateMICAinput(d = exp.log2.curr, filename="test_input.h5")
  match <- 0
  ID.curr <- h5read("test_input.h5", "ID")
  FN.curr <- h5read("test_input.h5", "fn")
  input.curr <- h5read("test_input.h5", "input")
  expect_equal(ID, ID.curr)
  expect_equal(FN, FN.curr)
  expect_equal(input, input.curr)
  file.remove("test_input.h5")
})


test_that("generateMICAinput errors when inputting invalid file type", {
  expect_error(generateMICAinput(d = exp.log2.curr, filename="improper_file_name.foo"), "Your filename *")
})


test_that("readMICAoutput works", {
  expect_equal(eset, eset.curr)
})


test_that("generateSJARACNeInput outputs correct .txt files", {
  generateSJARACNeInput(input_eset = eset.curr, funcType = "TF",
                        ref = "hg",
                        wd.src = "./sj_test/",
                        group_name = "ClusterRes")
  match <- 0
  sj_files <- list.files(path = "./sj_test/", recursive = TRUE)
  all_exp <- grep("*.exp", sj_files)
  g1 <- sj_files[all_exp[1]]
  hash1 <- md5sum(paste("./sj_test/", g1, sep=""))
  hash2 <- md5sum("../../testing_set/answerkey/sjaracne_input/Group1/G1.exp")
  if (hash1 == hash2) {
    match <- 1
  }
  expect_equal(match, 1)
  unlink("./sj_test/", recursive=TRUE)
})

test_that("GenerateActivityFromSJARACNe works", {
  expect_equal(acs, acs.curr)
})

test_that("get_activity works", {
  skip("TODO")
})


