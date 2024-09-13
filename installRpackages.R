install.packages("/home/ragg_1.3.2.9000.tar.gz", repos = NULL, type = "source")

packages <- list("devtools", "ggplot2", "magrittr", "data.table", "bayestestR")
for (package in packages) {
	install.packages(package)
}
devtools::install_github("Illumina/happyR")

packages <- list("ggplot2", "magrittr", "data.table", "bayestestR", "devtools", "happyR")
sink("InstalledRpackages.txt")
for (package in packages) {
	cat(package, ": ")
	tryCatch({
		version <- paste(unlist(packageVersion(package)), collapse='.')
		cat(version, '\n')
	}, error=function(e){
		cat("Not installed\n")
	})
}
sink()
