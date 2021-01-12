#!/usr/bin/env Rscript

packs <- c("optparse","data.table","R.utils","dplyr","ggplot2")

for (p in packs) {
  if( !require(p, character.only = T)) {
    print(p)
    install.packages( p,  repos = "https://cloud.r-project.org" )
  }
}