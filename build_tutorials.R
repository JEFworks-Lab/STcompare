
dir.create("docs/tutorials", recursive = TRUE, showWarnings = FALSE)

rmarkdown::render(
  "tutorials/tutorial-1.Rmd",
  output_file = "tutorial-1.html",
  output_dir = "docs/tutorials"
)

rmarkdown::render(
  "tutorials/tutorial-2.Rmd",
  output_file = "tutorial-2.html",
  output_dir = "docs/tutorials"
)