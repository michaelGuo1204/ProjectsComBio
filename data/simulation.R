library(tidyverse)
library(dyngen)
library(anndata)
#%%

#backbone <- backbone_linear_simple()
#config <-
#  initialise_model(
#    backbone = backbone,
#    num_cells = 500,
#    num_tfs = nrow(backbone$module_info),
#    num_targets = 50,
#    num_hks = 50,
#    verbose = FALSE,
#    download_cache_dir = tools::R_user_dir("dyngen", "data"),
#)
#%%
#out <- generate_dataset(
#  config,
#  format = "anndata",
#  make_plots = TRUE
#)
##%%
#linear_ad <- out$dataset
#linear_ad$write_h5ad("./linear.h5ad")
#ggsave("./linear.png", out$plot, width = 20, height = 20)
#%%

bifurbackbone <- backbone_bifurcating()
bifurconfig <-
  initialise_model(
    backbone = bifurbackbone,
    num_cells = 1000,
    num_tfs = nrow(bifurbackbone$module_info),
    num_targets = 50,
    num_hks = 50,
    verbose = FALSE,
    download_cache_dir = tools::R_user_dir("dyngen", "data"),
)
bifurout <- generate_dataset(
  bifurconfig,
  format = "anndata",
  make_plots = TRUE
)
#%%
bifur_ad <- bifurout$dataset
#bifur_ad$write_h5ad("./bifurcation.h5ad")
ggsave("./bifurcation.png", bifurout$plot, width = 20, height = 20)