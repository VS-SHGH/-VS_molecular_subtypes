# Exact composition proxies used by the historical 71.8% variance model.
# These are ESTIMATE-like signature means, not the official ESTIMATE output.

ESTIMATE_LIKE_STROMAL_GENES <- c(
  "ABI3BP", "ADAM12", "ADAM33", "ADAMTS2", "ADAMTS5", "ANGPT1", "ANGPTL1",
  "ANTXR1", "ASPN", "BGN", "C11orf96", "CCDC80", "CD248", "CDH11", "CLEC11A",
  "COL10A1", "COL11A1", "COL12A1", "COL1A1", "COL1A2", "COL3A1", "COL4A1",
  "COL5A1", "COL5A2", "COL5A3", "COL6A1", "COL6A2", "COL6A3", "COL8A1",
  "COMP", "CRISPLD2", "CSF1", "CYR61", "DCN", "DPT", "EFEMP2", "ELN",
  "FAP", "FBN1", "FBLN1", "FBLN2", "FBLN5", "FNDC1", "FN1", "GLT8D2",
  "HTRA1", "HTRA3", "ISLR", "LAMA4", "LAMB1", "LAMC1", "LOX", "LOXL1",
  "LOXL2", "LRRC15", "LUM", "MATN2", "MFAP5", "MGP", "MMP14", "MXRA5",
  "MXRA8", "NBL1", "NTM", "OLFML2B", "OMD", "PCOLCE", "PDGFRB", "PDLIM3",
  "PLAT", "PRRX1", "PTN", "RCN3", "RIN2", "SFRP4", "SPARC", "SPOCK1",
  "SRPX2", "SULF1", "TAGLN", "TGFB1", "THBS2", "THY1", "TIMP3", "TNFSF4",
  "VCAN", "WISP1", "ZFHX4"
)

ESTIMATE_LIKE_IMMUNE_GENES <- c(
  "APC", "ARHGAP25", "BTK", "C1QA", "C1QB", "C2", "C3AR", "C5AR1", "CCL13",
  "CCL19", "CCL22", "CCL5", "CCR2", "CCR5", "CD14", "CD163", "CD1C", "CD2",
  "CD209", "CD28", "CD3D", "CD3E", "CD4", "CD40", "CD40LG", "CD48", "CD5",
  "CD52", "CD53", "CD6", "CD68", "CD7", "CD74", "CD79A", "CD79B", "CD8A",
  "CD8B", "CD86", "CD96", "CIITA", "CLEC4A", "CORO1A", "CRTAM", "CSF2RB",
  "CTLA4", "CTSS", "CX3CR1", "CXCL10", "CXCL13", "CXCL9", "CYBB", "DOCK2",
  "EVI2B", "FASLG", "FLT3", "FPR3", "FYN", "GIMAP5", "GPR18", "GZMA",
  "HAVCR2", "HCLS1", "HLA-DMB", "HLA-E", "HLA-F", "HMMR", "ICOS", "IKZF1",
  "IL10RA", "IL2RA", "IL2RG", "IL7R", "IRF8", "ITGAL", "ITK", "KLRB1",
  "KLRK1", "LAG3", "LAIR1", "LCK", "LCP2", "LY86", "LY96", "MARCO", "MS4A1",
  "NCKAP1L", "NKG7", "P2RY13", "PDCD1", "PDCD1LG2", "PRF1", "PSMB9", "PTPRC",
  "PYHIN1", "RAC2", "SELPLG", "SLAMF6", "SLAMF8", "SPN", "ST8SIA4", "TARP",
  "TCL1A", "TLR10", "TLR7", "TLR8", "TNFRSF13B", "TNFRSF17", "TRAT1",
  "TRBC1", "TRBC2", "TRBV28", "TREM2", "UBA7", "WDFY4"
)

signature_mean <- function(expr_mat, signature) {
  available <- intersect(signature, rownames(expr_mat))
  if (length(available) < 5L) {
    stop("Too few genes available for signature: ", length(available))
  }
  colMeans(expr_mat[available, , drop = FALSE])
}
