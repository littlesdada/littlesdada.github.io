tpm <- function(data,species) {
  geneName <- rownames(data)
  infoFlag <- 1
  if (species == 'hsa') {
    geneLength <- read.delim('hg19_genes_length.tsv', header = TRUE, stringsAsFactors = TRUE)
  } else if (species == 'mmu') {
    geneLength <- read.delim('mm9_genes_length.tsv', header = TRUE, stringsAsFactors = TRUE)
  } else {
    print('no gene length information')
    infoFlag <- 0
  }
  # Calculate TPM
  if(infoFlag == 1) {
    rpk <- data
    gLength <- vector(mode="numeric", length=dim(data)[1])
    for (i in 1:dim(data)[1]) {
      dx <- which(geneLength$geneName == geneName[i])
      if(length(dx) > 0) {
        gLength[i] <- max(geneLength$geneLength[dx])
      } else {
        gLength[i] <- -1
      }
    }
    dx <- which(gLength == -1)
    rpk[dx, ] <- 0
    rpk <- rpk * 1000/kronecker(matrix(1, 1, dim(data)[2]), gLength)
    sscale <- apply(rpk, 2, sum)/10^6
    tpm <- rpk / kronecker(matrix(1, dim(data)[1], 1), sscale)
    write.table(tpm,"tpm.txt",sep="\t",quote = FALSE)
  }
}