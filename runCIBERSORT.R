# first, download CIBERSORT.R from CIBERSORT website
# setup run directory and files
source("CIBERSORT.R")

tdir <- tempfile('CIBERSORT')
dir.create(tdir)
owd <- setwd( tdir )
on.exit({
  unlink(tdir, recursive = TRUE)
  setwd(owd)
}) 
message("* Writing input files ... ", appendLF = FALSE)
write.table(reference, file = xf <- 'reference.tsv', sep = "\t", row.names = TRUE, col.names = NA)
write.table(y, file = yf <- 'mixture.tsv', sep = "\t", row.names = TRUE, col.names = NA)
message('OK')

# run
message("* Running CIBERSORT ... ", appendLF = FALSE)
res <- CIBERSORT(xf, yf, ...)
message('OK')

