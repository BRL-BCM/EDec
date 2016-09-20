# EDec
EDec (Epigenomic Deconvolution) is a technique that, starting from methylation and gene expression profiles of bulk tissue samples, infers cell type composition of each input sample as well as DNA methylation and gene transcription profiles of constituent cell types. 

#Installing EDec in R
You can install `EDec` straght from github. To install it:

**1**. If you haven't already done so, instal the `devtools` library.

    install.packages("devtools")
    
**2**. Install the `EDec` package.

    devtools::install.github("BRL-BCM/EDec")
    
**3**. Load the `EDec` package.

    library(EDec)
