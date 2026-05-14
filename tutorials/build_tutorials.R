
dir.create("docs/tutorials", recursive = TRUE, showWarnings = FALSE)

rmarkdown::render(
  "tutorials/acute-kidney-injury-10x-visium-rasterized.Rmd", 
  output_file = "acute-kidney-injury-10x-visium-rasterized.html", 
  output_dir = "docs/tutorials"
)

rmarkdown::render(
  "tutorials/tutorial-2.Rmd",
  output_file = "tutorial-2.html",
  output_dir = "docs/tutorials"
)