# Precompiled vignettes that depend on API key
# Must manually move image files from pkg/ to pkg/vignettes/ after knit

library(knitr)
knit("vignettes/phylotaR.Rmd.orig", "vignettes/phylotaR.Rmd")
