setwd("~/Documents/GitHub/Ancestral-Condition-Test/cleaned/Analyses/RPackage")


Rd2roxygen::Rd2roxygen("AncCond")

roxygen2::roxygenise("AncCond")



formatR::tidy_dir("AncCond/R")

styler::style_dir("AncCond/R")





devtools::check("AncCond")

devtools::build_manual("AncCond")





setwd("~/Desktop/projects/reporting/AncCond")

usethis::use_vignette("figures")





pkgdown::build_site()





formatR::tidy_dir("survivalFigs/R")

styler::style_dir("survivalFigs/R")

roxygen2::roxygenise("survivalFigs")

devtools::build_readme("survivalFigs")

devtools::check("survivalFigs")

devtools::build_manual("survivalFigs")

pkgdown::build_site("survivalFigs")
