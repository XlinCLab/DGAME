#!/usr/bin/env Rscript

# Function to silently check if a package is installed
is_package_installed <- function(package_name) {
	return(require(package_name, quietly = TRUE, character.only = TRUE))
}

# Function to load a package without printing its startup messages
load_package <- function(package_name) {
	suppressPackageStartupMessages(library(package_name, quietly = TRUE, character.only = TRUE))
}

# Function to check if a package is installed and install it if necessary
install_if_needed <- function(package_name) {
    if (!is_package_installed(package_name)) {
        install.packages(package_name, repos = "https://cran.r-project.org")
    }
	load_package(package_name)
}

# Function to install a list of packages if they are not already installed
install_packages_if_needed <- function(packages) {
	for (package in packages) {
		install_if_needed(package)
	}
}
