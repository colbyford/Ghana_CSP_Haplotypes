
# Update the libraries
library(devtools)
library(git2r)

path <- file.path(tempfile(pattern="Rswarm-"), "Rswarm")
dir.create(path, recursive=TRUE)
repo <- clone("https://github.com/lerch-a/Rswarm.git", path)
clone("https://github.com/torognes/swarm.git", file.path(path, "src", "swarm"))
install(path)

path <- file.path(tempfile(pattern="Rvsearch-"), "Rvsearch")
dir.create(path, recursive=TRUE)
repo <- clone("https://github.com/lerch-a/Rvsearch.git", path)
clone("https://github.com/torognes/vsearch.git", file.path(path, "src", "vsearch"), branch=) #v2.15.0
install(path)

path <- file.path(tempfile(pattern="NGmergeR-"), "NGmergeR")
dir.create(path, recursive=TRUE)
repo <- clone("https://github.com/lerch-a/NGmergeR.git", path)
clone("https://github.com/harvardinformatics/NGmerge.git", file.path(path, "src", "NGmerge"))
install(path)

# remove.packages("HaplotypR")
detach("package:HaplotypR", unload=TRUE)
devtools::install_github("lerch-a/HaplotypR")

install.packages(c("ape", "phangorn"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msa")
