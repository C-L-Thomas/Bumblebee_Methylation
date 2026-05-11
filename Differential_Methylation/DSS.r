library(DSS)
library(readr)
library(data.table)
library(rtracklayer)
library(topGO)
library(GO.db)
library(bsseq)
library(ggplot2)

################################ Preamble ################################

in_dir <- "." 		#Setwd
out_prefix <- "DSS"	#Set prefix
fdr_thresh <- 0.05	#Set FDR threshold

gene_anno_file <- "/path/genes_with_start_and_end.txt"	#File stating gene start and end
gff_file       <- "/path/GCF_910591885.1_iyBomTerr1.2_genomic.gff"	#gff
go_map_file       <- "/path/BomTerr1.2_Eggnog_go_map.tsv"	#Go terms

nodeSize <- 5	#Minimum node size for topgo
topNodes <- 200	#Number of go terms to print (if sig)

contrast_names <- c("Control_vs_5aza", "Control_vs_6aza", "5aza_vs_6aza", "repro_vs_sterile")

sample_names <- c("Col_1_6aza","Col_1_6aza","Col_1_control","Col_1_control","Col_1_5aza","Col_1_5aza",
  "Col_3_6aza","Col_3_6aza","Col_3_control","Col_3_control","Col_3_5aza","Col_3_5aza",
  "Col_4_6aza","Col_4_6aza","Col_4_control","Col_4_control","Col_4_5aza","Col_4_5aza",
  "Col_6_6aza","Col_6_6aza","Col_6_control","Col_6_control","Col_6_5aza","Col_6_5aza",
  "Col_7_6aza","Col_7_6aza","Col_7_control","Col_7_control","Col_7_5aza","Col_7_5aza",
  "Col_8_6aza","Col_8_6aza","Col_8_control","Col_8_control","Col_8_5aza","Col_8_5aza")

Drug   <- rep(c("5aza","5aza","Control","Control","6aza","6aza"), 6)
#Colony <- rep(c("1","3","4","6","7","8"), each = 6)
Colony <- rep(c("1","2","3","4","5","6"), each = 6)


Status <- c("sterile","repro","sterile","repro","repro","sterile",
  "sterile","repro","sterile","repro","sterile","repro",
  "repro","sterile","sterile","repro","repro","sterile",
  "sterile","repro","repro","sterile","sterile","repro",
  "sterile","repro","sterile","repro","repro","sterile",
  "repro","sterile","sterile","repro","sterile","repro")

file.list <- sort(list.files(in_dir, pattern="*binomial_results.txt", full.names=TRUE))	#Load files

samples <- lapply(file.list, function(x) {
  read_delim(x, "\t", col_names = TRUE,trim_ws = TRUE)})	#Load files

BSobj <- makeBSseqData(samples, sample_names)	#Convert samples to BSseq object

design <- data.frame(Drug   = factor(Drug),
  Colony = factor(Colony),
  Status = factor(Status))

################################ PCA ################################

meth <- getMeth(BSobj, type = "raw")
meth_filt <- meth[complete.cases(meth), , drop = FALSE] # remove rows with any NA first

row_sd <- apply(meth_filt, 1, sd)	#Calculate SD
meth_filt <- meth_filt[row_sd > 0, , drop = FALSE]	#Remove sites with no SD from PCA
pca <- prcomp(t(meth_filt), scale. = TRUE) #Perform PCA
var_exp <- (pca$sdev^2 / sum(pca$sdev^2)) * 100	#Calculate PCA variance 

df <- data.frame(pca$x, sample = colnames(meth_filt),
  Drug = design$Drug, Colony = design$Colony, Status = design$Status) #Combine metadata with PCA

df$Drug <- factor(df$Drug, levels = c("Control", "5aza", "6aza")) #Set levels
df$Status <- factor(df$Status, levels = c("sterile", "repro"))	#Set levels
df$Colony <- factor(df$Colony)	#Set levels

# plot
p <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(color = Colony, shape = Drug, fill = Status),
	size = 4, stroke = 1.2) +
  scale_shape_manual(values = c("Control" = 24, "5aza" = 21, "6aza" = 22)) +
  scale_fill_manual(values = c("sterile" = "white", "repro" = "black")) +
  guides(fill = guide_legend(override.aes = list(shape = 21, color = "black"))) +
  labs(x = paste0("PC1 (", round(var_exp[1],1), "%)"),y = paste0("PC2 (", round(var_exp[2],1), "%)"),
    color = "Colony",shape = "Drug",fill = "Status") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),legend.position = "right"
)	#Plot PCA

ggsave("pca_plot.png", p, width = 7, height = 5, dpi = 300) #Save

################################ Models ################################

design$Status <- relevel(design$Status, ref="sterile")

#Model A
designA <- copy(design)
designA$Drug <- relevel(designA$Drug, ref="Control") #Set reference
fitA <- DMLfit.multiFactor(BSobj, design=designA, formula=~Drug + Status + Colony) #Run model

#Reproductive status
dms_repro <- DMLtest.multiFactor(fitA, coef="Statusrepro") #Perform Wald tests
dms_repro_dt <- as.data.table(dms_repro)

setnames(dms_repro_dt, old = intersect("fdrs", names(dms_repro_dt)), new = "fdr")
dms_repro_sig <- dms_repro_dt[!is.na(fdr) & fdr < fdr_thresh] #Select significant sites 
dms_repro_sig[, contrast := "repro_vs_sterile"]

fwrite(dms_repro_sig, paste0(out_prefix, "_repro_vs_sterile_sigCpGs.tsv"), sep="\t") #Write results
cat("[repro_vs_sterile] sig CpGs:", nrow(dms_repro_sig), "\n")

#6 aza vs control
dms_6aza_ctl <- DMLtest.multiFactor(fitA, coef="Drug6aza") #Perform Wald tests
dms_6aza_ctl_dt <- as.data.table(dms_6aza_ctl)

setnames(dms_6aza_ctl_dt, old = intersect("fdrs", names(dms_6aza_ctl_dt)), new = "fdr")
dms_ctl_6aza_sig <- dms_6aza_ctl_dt[!is.na(fdr) & fdr < fdr_thresh]  #Select significant sites    
dms_ctl_6aza_sig[, contrast := "Control_vs_6aza"]

fwrite(dms_ctl_6aza_sig, paste0(out_prefix, "_Control_vs_6aza_sigCpGs.tsv"), sep="\t") #Write results
cat("[Control_vs_6aza] sig CpGs:", nrow(dms_ctl_6aza_sig), "\n")

#Colony
dms_colony <- DMLtest.multiFactor(fitA, term="Colony") #multi-degree-of-freedom Wald test across all colony levels
dms_colony_dt <- as.data.table(dms_colony)

setnames(dms_colony_dt, old = intersect("pvals", names(dms_colony_dt)), new = "pval")
dms_colony_dt[, fdr := p.adjust(pval, method = "BH")] #Need to perform fdr because of multi-degree freedom wald test
dms_colony_sig <- dms_colony_dt[!is.na(fdr) & fdr < fdr_thresh] #Select significant sites
dms_colony_sig[, contrast := "Colony_effect"]

fwrite(dms_colony_sig, paste0(out_prefix, "_Colony_effect_sigCpGs.tsv"), sep = "\t") #Write results
cat("[Colony_effect] sig CpGs:", nrow(dms_colony_sig), "\n")

#Model B
designB <- copy(design)
designB$Drug <- relevel(designB$Drug, ref="5aza")
fitB <- DMLfit.multiFactor(BSobj, design=designB, formula=~Drug + Status + Colony)

# Control vs 5aza
dms_ctl_5aza <- DMLtest.multiFactor(fitB, coef="DrugControl") #Perform Wald tests
dms_ctl_5aza_dt <- as.data.table(dms_ctl_5aza)

setnames(dms_ctl_5aza_dt, old = intersect("fdrs", names(dms_ctl_5aza_dt)), new = "fdr")
dms_ctl_5aza_sig <- dms_ctl_5aza_dt[!is.na(fdr) & fdr < fdr_thresh] #Select significant sites
dms_ctl_5aza_sig[, contrast := "Control_vs_5aza"]

fwrite(dms_ctl_5aza_sig, paste0(out_prefix, "_Control_vs_5aza_sigCpGs.tsv"), sep="\t") #Write results
cat("[Control_vs_5aza] sig CpGs:", nrow(dms_ctl_5aza_sig), "\n")

#  6aza vs 5aza
dms_6aza_5aza <- DMLtest.multiFactor(fitB, coef="Drug6aza") #Perform Wald tests
dms_6aza_5aza_dt <- as.data.table(dms_6aza_5aza)

setnames(dms_6aza_5aza_dt, old = intersect("fdrs", names(dms_6aza_5aza_dt)), new = "fdr")
dms_5aza_6aza_sig <- dms_6aza_5aza_dt[!is.na(fdr) & fdr < fdr_thresh] #Select significant sites
dms_5aza_6aza_sig[, contrast := "5aza_vs_6aza"]

fwrite(dms_5aza_6aza_sig, paste0(out_prefix, "_5aza_vs_6aza_sigCpGs.tsv"), sep="\t") #Write results
cat("[5aza_vs_6aza] sig CpGs:", nrow(dms_5aza_6aza_sig), "\n")

################################ Gene Annotation ################################

genes <- fread(gene_anno_file) #Read gene regions file

genes[, `:=`(
  chr   = as.character(chr),
  start = as.integer(start),
  end   = as.integer(end),
  geneID = sub("^gene-", "", as.character(geneID))
)]

setkey(genes, chr, start, end)

# Control_vs_5aza
dt_ctl_5aza <- fread(paste0(out_prefix, "_Control_vs_5aza_sigCpGs.tsv")) #Read DMS file
dt_ctl_5aza[, `:=`(chr = as.character(chr), pos = as.integer(pos),
                   cpg_start = as.integer(pos), cpg_end = as.integer(pos))]
setkey(dt_ctl_5aza, chr, cpg_start, cpg_end)

dt_ctl_5aza_withGenes <- foverlaps(
  dt_ctl_5aza, genes,
  by.x = c("chr","cpg_start","cpg_end"),
  by.y = c("chr","start","end"),
  type = "within",
  nomatch = NA)

fwrite(dt_ctl_5aza_withGenes,
       paste0(out_prefix, "_Control_vs_5aza_sigCpGs_withGenes.tsv"),
       sep="\t") #Write Results

cat("[Control_vs_5aza] CpGs within genes:", sum(!is.na(dt_ctl_5aza_withGenes$geneID)), "\n")

# Control_vs_6aza
dt_ctl_6aza <- fread(paste0(out_prefix, "_Control_vs_6aza_sigCpGs.tsv")) #Read DMS file
dt_ctl_6aza[, `:=`(chr = as.character(chr), pos = as.integer(pos),
                   cpg_start = as.integer(pos), cpg_end = as.integer(pos))]
setkey(dt_ctl_6aza, chr, cpg_start, cpg_end)

dt_ctl_6aza_withGenes <- foverlaps(
  dt_ctl_6aza, genes,
  by.x = c("chr","cpg_start","cpg_end"),
  by.y = c("chr","start","end"),
  type = "within",
  nomatch = NA)

fwrite(dt_ctl_6aza_withGenes,
       paste0(out_prefix, "_Control_vs_6aza_sigCpGs_withGenes.tsv"),
       sep="\t") #Write Results

cat("[Control_vs_6aza] CpGs within genes:", sum(!is.na(dt_ctl_6aza_withGenes$geneID)), "\n")

#5aza_vs_6aza
dt_5aza_6aza <- fread(paste0(out_prefix, "_5aza_vs_6aza_sigCpGs.tsv"))
dt_5aza_6aza[, `:=`(chr = as.character(chr), pos = as.integer(pos),
                    cpg_start = as.integer(pos), cpg_end = as.integer(pos))]
setkey(dt_5aza_6aza, chr, cpg_start, cpg_end)

dt_5aza_6aza_withGenes <- foverlaps(
  dt_5aza_6aza, genes,
  by.x = c("chr","cpg_start","cpg_end"),
  by.y = c("chr","start","end"),
  type = "within",
  nomatch = NA
)

fwrite(dt_5aza_6aza_withGenes,
       paste0(out_prefix, "_5aza_vs_6aza_sigCpGs_withGenes.tsv"),
       sep="\t") #Write Results
cat("[5aza_vs_6aza] CpGs within genes:", sum(!is.na(dt_5aza_6aza_withGenes$geneID)), "\n")

#repro_vs_sterile
dt_repro <- fread(paste0(out_prefix, "_repro_vs_sterile_sigCpGs.tsv"))
dt_repro[, `:=`(chr = as.character(chr), pos = as.integer(pos),
                cpg_start = as.integer(pos), cpg_end = as.integer(pos))]
setkey(dt_repro, chr, cpg_start, cpg_end)

dt_repro_withGenes <- foverlaps(
  dt_repro, genes,
  by.x = c("chr","cpg_start","cpg_end"),
  by.y = c("chr","start","end"),
  type = "within",
  nomatch = NA
)

fwrite(dt_repro_withGenes,
       paste0(out_prefix, "_repro_vs_sterile_sigCpGs_withGenes.tsv"),
       sep="\t") #Write Results
cat("[repro_vs_sterile] CpGs within genes:", sum(!is.na(dt_repro_withGenes$geneID)), "\n")

# Colony_effect

dt_colony <- fread(paste0(out_prefix, "_Colony_effect_sigCpGs.tsv")) #Read DMS file
dt_colony[, `:=`(chr = as.character(chr), pos = as.integer(pos),
                 cpg_start = as.integer(pos), cpg_end = as.integer(pos))]
setkey(dt_colony, chr, cpg_start, cpg_end)

dt_colony_withGenes <- foverlaps(
  dt_colony, genes,
  by.x = c("chr","cpg_start","cpg_end"),
  by.y = c("chr","start","end"),
  type = "within",
  nomatch = NA)

fwrite(dt_colony_withGenes,
       paste0(out_prefix, "_Colony_effect_sigCpGs_withGenes.tsv"),
       sep = "\t") #Write Results

cat("[Colony_effect] CpGs within genes:", sum(!is.na(dt_colony_withGenes$geneID)), "\n")

################################ Gene Annotation ################################
go_map <- fread(go_map_file, sep="\t", header=FALSE) #Read file
setnames(go_map, c("geneID","GO","Aspect"))
go_map <- unique(go_map)

gene2go_BP <- lapply(split(go_map[Aspect=="P"]$GO, go_map[Aspect=="P"]$geneID), unique)
gene2go_MF <- lapply(split(go_map[Aspect=="F"]$GO, go_map[Aspect=="F"]$geneID), unique)
gene2go_CC <- lapply(split(go_map[Aspect=="C"]$GO, go_map[Aspect=="C"]$geneID), unique)

universe_all <- unique(go_map$geneID)
universe_BP <- intersect(universe_all, names(gene2go_BP))
universe_MF <- intersect(universe_all, names(gene2go_MF))
universe_CC <- intersect(universe_all, names(gene2go_CC))

cat("GO universe sizes: all=", length(universe_all),
    " BP=", length(universe_BP),
    " MF=", length(universe_MF),
    " CC=", length(universe_CC), "\n")

# Control_vs_5aza
dt <- fread(paste0(out_prefix, "_Control_vs_5aza_sigCpGs_withGenes.tsv")) #Read File
sig_genes <- unique(dt[!is.na(geneID) & geneID != "", geneID])

sig_genes_use <- intersect(sig_genes, universe_BP)
  geneList <- factor(as.integer(universe_BP %in% sig_genes_use))
  names(geneList) <- universe_BP
  GOdata <- new("topGOdata", ontology="BP", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2go_BP, nodeSize=nodeSize)
  res_classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
  res_elim    <- runTest(GOdata, algorithm="elim",    statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=res_classic, elimFisher=res_elim,
                  orderBy="elimFisher", ranksOf="classicFisher", topNodes=topNodes)
  fwrite(as.data.table(tab), paste0(out_prefix, "_topGO_Control_vs_5aza_BP.tsv"), sep="\t")

sig_genes_use <- intersect(sig_genes, universe_MF)
geneList <- factor(as.integer(universe_MF %in% sig_genes_use))
  names(geneList) <- universe_MF
  GOdata <- new("topGOdata", ontology="MF", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2go_MF, nodeSize=nodeSize)
  res_classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
  res_elim    <- runTest(GOdata, algorithm="elim",    statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=res_classic, elimFisher=res_elim,
                  orderBy="elimFisher", ranksOf="classicFisher", topNodes=topNodes)
  fwrite(as.data.table(tab), paste0(out_prefix, "_topGO_Control_vs_5aza_MF.tsv"), sep="\t")

sig_genes_use <- intersect(sig_genes, universe_CC)
  geneList <- factor(as.integer(universe_CC %in% sig_genes_use))
  names(geneList) <- universe_CC
  GOdata <- new("topGOdata", ontology="CC", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2go_CC, nodeSize=nodeSize)
  res_classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
  res_elim    <- runTest(GOdata, algorithm="elim",    statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=res_classic, elimFisher=res_elim,
                  orderBy="elimFisher", ranksOf="classicFisher", topNodes=topNodes)
  fwrite(as.data.table(tab), paste0(out_prefix, "_topGO_Control_vs_5aza_CC.tsv"), sep="\t")

# Control_vs_6aza
dt <- fread(paste0(out_prefix, "_Control_vs_6aza_sigCpGs_withGenes.tsv"))
sig_genes <- unique(dt[!is.na(geneID) & geneID != "", geneID])

sig_genes_use <- intersect(sig_genes, universe_BP)
  geneList <- factor(as.integer(universe_BP %in% sig_genes_use))
  names(geneList) <- universe_BP
  GOdata <- new("topGOdata", ontology="BP", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2go_BP, nodeSize=nodeSize)
  res_classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
  res_elim    <- runTest(GOdata, algorithm="elim",    statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=res_classic, elimFisher=res_elim,
                  orderBy="elimFisher", ranksOf="classicFisher", topNodes=topNodes)
  fwrite(as.data.table(tab), paste0(out_prefix, "_topGO_Control_vs_6aza_BP.tsv"), sep="\t")

sig_genes_use <- intersect(sig_genes, universe_MF)
  geneList <- factor(as.integer(universe_MF %in% sig_genes_use))
  names(geneList) <- universe_MF
  GOdata <- new("topGOdata", ontology="MF", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2go_MF, nodeSize=nodeSize)
  res_classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
  res_elim    <- runTest(GOdata, algorithm="elim",    statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=res_classic, elimFisher=res_elim,
                  orderBy="elimFisher", ranksOf="classicFisher", topNodes=topNodes)
  fwrite(as.data.table(tab), paste0(out_prefix, "_topGO_Control_vs_6aza_MF.tsv"), sep="\t")

sig_genes_use <- intersect(sig_genes, universe_CC)
  geneList <- factor(as.integer(universe_CC %in% sig_genes_use))
  names(geneList) <- universe_CC
  GOdata <- new("topGOdata", ontology="CC", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2go_CC, nodeSize=nodeSize)
  res_classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
  res_elim    <- runTest(GOdata, algorithm="elim",    statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=res_classic, elimFisher=res_elim,
                  orderBy="elimFisher", ranksOf="classicFisher", topNodes=topNodes)
  fwrite(as.data.table(tab), paste0(out_prefix, "_topGO_Control_vs_6aza_CC.tsv"), sep="\t")

# 5aza_vs_6aza
dt <- fread(paste0(out_prefix, "_5aza_vs_6aza_sigCpGs_withGenes.tsv"))
sig_genes <- unique(dt[!is.na(geneID) & geneID != "", geneID])

sig_genes_use <- intersect(sig_genes, universe_BP)
  geneList <- factor(as.integer(universe_BP %in% sig_genes_use))
  names(geneList) <- universe_BP
  GOdata <- new("topGOdata", ontology="BP", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2go_BP, nodeSize=nodeSize)
  res_classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
  res_elim    <- runTest(GOdata, algorithm="elim",    statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=res_classic, elimFisher=res_elim,
                  orderBy="elimFisher", ranksOf="classicFisher", topNodes=topNodes)
  fwrite(as.data.table(tab), paste0(out_prefix, "_topGO_5aza_vs_6aza_BP.tsv"), sep="\t")

sig_genes_use <- intersect(sig_genes, universe_MF)
  geneList <- factor(as.integer(universe_MF %in% sig_genes_use))
  names(geneList) <- universe_MF
  GOdata <- new("topGOdata", ontology="MF", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2go_MF, nodeSize=nodeSize)
  res_classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
  res_elim    <- runTest(GOdata, algorithm="elim",    statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=res_classic, elimFisher=res_elim,
                  orderBy="elimFisher", ranksOf="classicFisher", topNodes=topNodes)
  fwrite(as.data.table(tab), paste0(out_prefix, "_topGO_5aza_vs_6aza_MF.tsv"), sep="\t")

sig_genes_use <- intersect(sig_genes, universe_CC)
  geneList <- factor(as.integer(universe_CC %in% sig_genes_use))
  names(geneList) <- universe_CC
  GOdata <- new("topGOdata", ontology="CC", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2go_CC, nodeSize=nodeSize)
  res_classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
  res_elim    <- runTest(GOdata, algorithm="elim",    statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=res_classic, elimFisher=res_elim,
                  orderBy="elimFisher", ranksOf="classicFisher", topNodes=topNodes)
  fwrite(as.data.table(tab), paste0(out_prefix, "_topGO_5aza_vs_6aza_CC.tsv"), sep="\t")

# repro_vs_sterile
dt <- fread(paste0(out_prefix, "_repro_vs_sterile_sigCpGs_withGenes.tsv"))
sig_genes <- unique(dt[!is.na(geneID) & geneID != "", geneID])

sig_genes_use <- intersect(sig_genes, universe_BP)
  geneList <- factor(as.integer(universe_BP %in% sig_genes_use))
  names(geneList) <- universe_BP
  GOdata <- new("topGOdata", ontology="BP", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2go_BP, nodeSize=nodeSize)
  res_classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
  res_elim    <- runTest(GOdata, algorithm="elim",    statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=res_classic, elimFisher=res_elim,
                  orderBy="elimFisher", ranksOf="classicFisher", topNodes=topNodes)
  fwrite(as.data.table(tab), paste0(out_prefix, "_topGO_repro_vs_sterile_BP.tsv"), sep="\t")

sig_genes_use <- intersect(sig_genes, universe_MF)
  geneList <- factor(as.integer(universe_MF %in% sig_genes_use))
  names(geneList) <- universe_MF
  GOdata <- new("topGOdata", ontology="MF", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2go_MF, nodeSize=nodeSize)
  res_classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
  res_elim    <- runTest(GOdata, algorithm="elim",    statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=res_classic, elimFisher=res_elim,
                  orderBy="elimFisher", ranksOf="classicFisher", topNodes=topNodes)
  fwrite(as.data.table(tab), paste0(out_prefix, "_topGO_repro_vs_sterile_MF.tsv"), sep="\t")

sig_genes_use <- intersect(sig_genes, universe_CC)
  geneList <- factor(as.integer(universe_CC %in% sig_genes_use))
  names(geneList) <- universe_CC
  GOdata <- new("topGOdata", ontology="CC", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2go_CC, nodeSize=nodeSize)
  res_classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
  res_elim    <- runTest(GOdata, algorithm="elim",    statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=res_classic, elimFisher=res_elim,
                  orderBy="elimFisher", ranksOf="classicFisher", topNodes=topNodes)
  fwrite(as.data.table(tab), paste0(out_prefix, "_topGO_repro_vs_sterile_CC.tsv"), sep="\t")

# Colony_effect
dt <- fread(paste0(out_prefix, "_Colony_effect_sigCpGs_withGenes.tsv"))
sig_genes <- unique(dt[!is.na(geneID) & geneID != "", geneID])

sig_genes_use <- intersect(sig_genes, universe_BP)
  geneList <- factor(as.integer(universe_BP %in% sig_genes_use))
  names(geneList) <- universe_BP
  GOdata <- new("topGOdata", ontology="BP", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2go_BP, nodeSize=nodeSize)
  res_classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
  res_elim    <- runTest(GOdata, algorithm="elim",    statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=res_classic, elimFisher=res_elim,
                  orderBy="elimFisher", ranksOf="classicFisher", topNodes=topNodes)
  fwrite(as.data.table(tab), paste0(out_prefix, "_topGO_Colony_effect_BP.tsv"), sep="\t")

sig_genes_use <- intersect(sig_genes, universe_MF)
  geneList <- factor(as.integer(universe_MF %in% sig_genes_use))
  names(geneList) <- universe_MF
  GOdata <- new("topGOdata", ontology="MF", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2go_MF, nodeSize=nodeSize)
  res_classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
  res_elim    <- runTest(GOdata, algorithm="elim",    statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=res_classic, elimFisher=res_elim,
                  orderBy="elimFisher", ranksOf="classicFisher", topNodes=topNodes)
  fwrite(as.data.table(tab), paste0(out_prefix, "_topGO_Colony_effect_MF.tsv"), sep="\t")

sig_genes_use <- intersect(sig_genes, universe_CC)
  geneList <- factor(as.integer(universe_CC %in% sig_genes_use))
  names(geneList) <- universe_CC
  GOdata <- new("topGOdata", ontology="CC", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2go_CC, nodeSize=nodeSize)
  res_classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
  res_elim    <- runTest(GOdata, algorithm="elim",    statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=res_classic, elimFisher=res_elim,
                  orderBy="elimFisher", ranksOf="classicFisher", topNodes=topNodes)
  fwrite(as.data.table(tab), paste0(out_prefix, "_topGO_Colony_effect_CC.tsv"), sep="\t")

## ==========================================
## 4) Check overlaps
## =========================================

ctl_5aza <- fread("DSS_Control_vs_5aza_sigCpGs.tsv")[, paste0(chr,":",pos)]
ctl_6aza <- fread("DSS_Control_vs_6aza_sigCpGs.tsv")[, paste0(chr,":",pos)]
a5_a6    <- fread("DSS_5aza_vs_6aza_sigCpGs.tsv")[, paste0(chr,":",pos)]
repro    <- fread("DSS_repro_vs_sterile_sigCpGs.tsv")[, paste0(chr,":",pos)]

length(intersect(ctl_5aza, ctl_6aza))
length(intersect(ctl_5aza, a5_a6))
length(intersect(ctl_5aza, repro))
length(intersect(repro, a5_a6))
length(intersect(ctl_6aza, a5_a6))
length(intersect(repro, ctl_6aza))
