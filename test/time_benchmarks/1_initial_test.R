
# inital run 04-23-26 


output <- SEraster::rasterizeGeneExpression(
  spe_list, 
  resolution = 200, 
  BPPARAM = BiocParallel::SerialParam()
)

plt1 <- SEraster::plotRaster(output[[1]], name = "total counts") + ggplot2::theme_void()
plt2 <- SEraster::plotRaster(output[[2]], name = "total counts") + ggplot2::theme_void()

library(patchwork)
plt1 + plt2


set.seed(0)
start_time <- Sys.time()
deltaList <- replicate(length(rownames(output$target)), c(0.01, 0.05, seq(0.1, 0.9, .1)))
merfishCorrelation <- spatialCorrelationGeneExpIterPermutations(
  output, 
  deltaX = deltaList, 
  deltaY = deltaList, 
  nThreads = 5
)
end_time <- Sys.time()
print(end_time - start_time) 

save(merfishCorrelation, file = "test/time_benchmarks/outputs/0423initial.Rdata")

# NOTE: this run returned a lot of correct p-values with 0.000000 

# Finding resolutions for test ====================================================


resolutions <- c(50, 100, 200, 250, 300, 350, 400, 500, 600)

dir.create("test/time_benchmarks/raster_plots", recursive = TRUE, showWarnings = FALSE)

library(ggplot2)

for (res in resolutions) {
  
  output <- SEraster::rasterizeGeneExpression(
    spe_list,
    resolution = res,
    BPPARAM = BiocParallel::MulticoreParam(workers = 22)
  )
  
  plots <- lapply(seq_along(output), function(i) {
    SEraster::plotRaster(output[[i]]) +
      ggplot2::labs(title = paste0(names(spe_list)[i], " | resolution = ", res))
  })
  
  combined_plot <- wrap_plots(plots)
  
  ggsave(
    filename = file.path(
      "test/time_benchmarks/raster_plots",
      paste0("raster_res_", res, ".png")
    ),
    plot = combined_plot,
    width = 10,
    height = 5,
    dpi = 300
  )
}


