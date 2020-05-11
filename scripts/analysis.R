library(Seurat, quietly = T)
library(dplyr, quietly = T)
library(Matrix, quietly = T)
library(readr, quietly = T)
library(cowplot)
library(org.Mm.eg.db)
library(gdata)
library(ggpubr)
library(ggrepel)
library(biomaRt)
require(org.Mm.eg.db)
library(patchwork)
library(ggplot2)
source("functions.R")
library('purrr')
source("tmod_prep.R")

feature_annot<-readRDS("feature_annot_gencode.vM6.rds")

matrix_dir = paste("data/", "Brain_4C", sep="")
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = make.unique(feature.names$V2)

Brain_4C <- CreateSeuratObject(counts = mat, project = "Brain_4C", min.cells=5)
Brain_4C <- PercentageFeatureSet(Brain_4C, pattern = "^mt-", col.name = "percent.mt")
Brain_4C <- PercentageFeatureSet(Brain_4C, pattern = "^Rp[sl]", col.name = "percent.ribo")

Brain_4C <- PercentageFeatureSet(Brain_4C, pattern = "^mt-", col.name = "percent.mt")
Brain_4C$stim <- "Brain_4C"

write_rds(Brain_4C, "results/Brain_4C.rds")
Brain_4C<-read_rds("results/Brain_4C.rds")
Brain_4Cs <- subset(Brain_4C, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25)
Brain_4Cs <- NormalizeData(Brain_4Cs, verbose = FALSE)
Brain_4Cs <- FindVariableFeatures(Brain_4Cs, selection.method = "vst", nfeatures = 2000)


matrix_dir = paste("data/", "Brain_37C", sep="")
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = make.unique(feature.names$V2)

Brain_37C <- CreateSeuratObject(counts = mat, project = "Brain_37C", min.cells=5)
Brain_37C <- PercentageFeatureSet(Brain_37C, pattern = "^mt-", col.name = "percent.mt")
Brain_37C <- PercentageFeatureSet(Brain_37C, pattern = "^Rp[sl]", col.name = "percent.ribo")
vln37c<-VlnPlot(Brain_37C, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol=1, pt.size = 0.5)
Brain_37C <- PercentageFeatureSet(Brain_37C, pattern = "^mt-", col.name = "percent.mt")
Brain_37C$stim <- "Brain_37C"

write_rds(Brain_37C, "results/Brain_37C.rds")
Brain_37C<-read_rds("results/Brain_37C.rds")
Brain_37Cs <- subset(Brain_37C, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25)
Brain_37Cs <- NormalizeData(Brain_37Cs, verbose = FALSE)
Brain_37Cs <- FindVariableFeatures(Brain_37Cs, selection.method = "vst", nfeatures = 2000)



rds_file <-"Seurat3_clusters.rds"

if (!file.exists(rds_file)) {
  
  
  brain.anchors <- FindIntegrationAnchors(object.list = list(Brain_4Cs, Brain_37Cs), dims = 1:15)
  brain.combined <- IntegrateData(anchorset = brain.anchors, dims = 1:15)
  DefaultAssay(brain.combined) <- "integrated"

  brain.combined <- ScaleData(brain.combined, verbose = FALSE)
  brain.combined <- RunPCA(brain.combined,  verbose = FALSE)
  #VizDimLoadings(brain.combined, dims = 1:2, reduction = "pca")
  #DimHeatmap(brain.combined, dims = 1:15, cells = 500, balanced = TRUE)
  #ElbowPlot(brain.combined)
  brain.combined <- RunUMAP(brain.combined, reduction = "pca", dims = 1:15)
  brain.combined <- FindNeighbors(brain.combined, reduction = "pca", dims = 1:15)
  brain.combined_bf.Clust<-brain.combined
  brain.combined <- FindClusters(brain.combined, resolution = 0.41)
  
  
  write_rds(brain.combined, rds_file)
}else{
  brain.combined = readRDS(rds_file)
}

DefaultAssay(brain.combined) <- "RNA"

brain.combined_clean <- clean_up(brain.combined, nsd=3, nnn=10)



brain.combined_clean.rn<-RenameIdents(brain.combined_clean, '0'='Pyramidal', '1'='DG', '2'='Microglia',
                                      '3'='Microglia', '4'='Astrocytes', '5'='Oligodendr.', '6'='NPC', '7'='Fibroblast-like',
                                      '8'='OPC','9'='Microglia','10'='Oligodendr.', 
                                      '11'='Endothelial', '12'='BAM', '13'='Mural', 
                                      '14'='Astrocytes', '15'='Microglia', '16'='Polydendr.', 
                                      '17'='Astrocytes', '18'='Oligodendr.','19'='Endothelial')


brain.combined2<-brain.combined_clean.rn
brain.combined2$celltype.stim <- paste(Idents(brain.combined2), brain.combined2$stim, sep = "_")
brain.combined2$celltype <- Idents(brain.combined2)
Idents(brain.combined2) <- "celltype.stim"


cellcount<-table(brain.combined2$seurat_clusters, brain.combined2@meta.data$orig.ident) %>% as.data.frame() %>% tbl_df() %>% tidyr::spread(Var2, Freq) %>% mutate(cluster=as.character(Var1)) %>% dplyr::select(-Var1) 

cells<-c('Pyramidal','DG','Microglia', 'Microglia', 
         'Astrocytes', 'Oligodendr.', 'NPC', 'Fibroblast-like',
         'OPC', 'Microglia','10'='Oligodendr.', 
         'Endothelial', 'BAM', 'Mural', 
         'Astrocytes', 'Microglia','Polydendr.', 
         'Astrocytes', 'Oligodendr.','Endothelial')
cellcount$cells<-cells


colors=c("royalblue1", "royalblue3", "seagreen3", "goldenrod1", "firebrick2", "plum4", "darkorchid4", "violet", "cadetblue4", "yellow3", "peru","black")

cellcount_sum<-cellcount %>% group_by(cells) %>% summarise(Brain_4C=sum(Brain_4C), Brain_37C=sum(Brain_37C))

color_table<-cbind(cells %>% unique(), colors )
colnames(color_table)<-c("cells", "colors")
color_table<-tbl_df(color_table)
cellcount_bar<- cellcount_sum %>% mutate('37C'=Brain_37C/sum(Brain_37C),'4C'=Brain_4C/sum(Brain_4C) ) %>% 
  dplyr::select(-Brain_4C, -Brain_37C) %>% 
  tidyr::gather(Protocol, count, `37C`:`4C`) %>% left_join(color_table)



rds_file <- "Seurat3_markers.rds"

if (!file.exists(rds_file)) {
  allmarkers <- FindAllMarkers(brain.combined_clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write_rds(allmarkers, rds_file)
}else{
  allmarkers = readRDS(rds_file)
}




pdf("cell_count_fraction_Fig1.pdf", height=5, width=2)
cellcount_bar%>% ggplot(aes(x = Protocol, y = count, fill = cells))+
  geom_bar(stat="identity")+ylab("Cell fraction")+theme(legend.position = "none")+scale_fill_manual(values=cellcount_bar$colors)
dev.off()


pdf("umap_Fig1.pdf", height = 5, width = 10)
DimPlot(brain.combined_clean.rn,  reduction = "umap", pt.size=1.2, cols=colors, split.by = "stim")+theme(legend.position = "none")
dev.off()


pdf("markers_Fig1.pdf")

colors_vln<-c("royalblue1", "royalblue1", "royalblue2", "royalblue2", "seagreen3", "seagreen3", "goldenrod1", 
              "goldenrod1", "firebrick2", "firebrick2","plum4", "plum4", "darkorchid4","darkorchid4",
              "violet","violet", "cadetblue4", "cadetblue4","peru", "peru")

vln_features=c("Igsf9b","Sez6", "Prox1","Synpr", "P2ry12", "Cx3cr1", "Slc1a2", "Slc1a3",
               "Mag", "Mbp", "Hopx", "Sox9", "Dcn","Col1a2", "Pdgfra", "Vcan",
               "Flt1","Ly6c1", "Fcgr1","Axl", "Abcc9", "Ube2c")


StackedVlnPlot(obj = brain.combined_clean.rn, features = vln_features)
dev.off()




pdf("vln_Plots_suppl.pdf", height = 24, width = 18)
plots<-VlnPlot(brain.combined2, features = c("Jun", "Egr1", "Meg3", "Hspa1a", "Hspa1b", "Jund", "Junb", "Hspa8", "Plp1", "Ly6h", "Malat1", "Snrpg"), split.by = "stim", group.by = "celltype", pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 2)
dev.off()

pdf("featureplot_suppl.pdf", height = 5, width = 20)
plots<-FeaturePlot(brain.combined_clean.rn, features=c("Fos", "Rplp0", "Rplp1", "Rps28" ,"Rps29"), combine=F, split.by = "stim")
CombinePlots(plots = plots, ncol=5)
dev.off()



pdf("microglia_zoom_Fig2.pdf", height = 5, width = 10)
DimPlot(brain.combined_clean.rn,  reduction = "umap", pt.size=1.2, cols=colors, split.by = "stim")+theme(legend.position = "none")+ylim(-3,5)+xlim(2,15)
dev.off()


pdf("microglia_Jun_Fos_Fig2.pdf", height = 10, width = 10)

plots<-FeaturePlot(brain.combined_clean.rn, features=c("Jun", "Fos"), split.by = "stim", col=c("grey", "seagreen"), combine = F)
plots <- lapply(X = plots, FUN = function(p) p + xlim(c(2, 15)) + ylim(c(-3, 5)))
CombinePlots(plots = plots)

dev.off()


rds_file <-"DE_list.rds"

if (!file.exists(rds_file)) {

  
  de_list<-list()
  
  for (i in levels(brain.combined_clean.rn)){
    print (i)
    if (i=='OPC') { # no OPCs found in 37C protocol
      next
    }
    
    de_list[[i]]<-FindMarkers(brain.combined2, ident.1 = paste(i, "Brain_4C", sep="_"), ident.2 = paste(i, "Brain_37C", sep="_"), verbose = FALSE, logfc.threshold=0)
    
  }
  
  de_list[["Neurons"]] <- FindMarkers(brain.combined2, ident.1 = c("DG_Brain_4C", "Pyramidal_Brain_4C"), ident.2 = c("DG_Brain_37C", "Pyramidal_Brain_37C"), verbose = FALSE, logfc.threshold=0)
  
  WriteXLS::WriteXLS(de_list,
                     ExcelFileName = "DE_list.xls", SheetNames = names(de_list), col.names = TRUE, row.names = TRUE)
  
  
  
  write_rds(de_list,  rds_file)
}else{
  de_list<-readRDS(rds_file)
}


############# DE in Neurons ###################

neurons_s <- subset(brain.combined_clean.rn, idents = c("Pyramidal", "DG"))
Idents(neurons_s) <- "stim"
avg.neurons <- log1p(AverageExpression(neurons_s, verbose = FALSE)$RNA)
avg.neurons$gene <- rownames(avg.neurons)

neurons.genes.to.label <- (avg.neurons %>% filter(abs(Brain_4C-Brain_37C)>1.5))$gene



neurons_plot_scatter <- ggplot(avg.neurons, aes(Brain_4C, Brain_37C)) + geom_point() + ggtitle("Neurons")
neurons_plot_scatter<-neurons_plot_scatter + geom_text_repel(data=subset(avg.neurons, gene %in% neurons.genes.to.label),aes(label=gene), size=4, colour="red3")


pdf("neurons_scatter_Fig2.pdf")
neurons_plot_scatter
dev.off()


############# DE in Astrocytes ##################

astrocytes_s <- subset(brain.combined_clean.rn, idents = c("Astrocytes"))
Idents(astrocytes_s) <- "stim"
avg.astrocytes <- log1p(AverageExpression(astrocytes_s, verbose = FALSE)$RNA)
avg.astrocytes$gene <- rownames(avg.astrocytes)

astrocytes.genes.to.label <- (avg.astrocytes %>% filter(abs(Brain_4C-Brain_37C)>1.5))$gene

astrocytes_plot_scatter <- ggplot(avg.astrocytes, aes(Brain_4C, Brain_37C)) + geom_point() + ggtitle("Astrocytes")
#astrocytes_plot_scatter <- LabelPoints(plot = astrocytes_plot_scatter, points = astrocytes.genes.to.label, repel = TRUE)
astrocytes_plot_scatter<-astrocytes_plot_scatter + geom_text_repel(data=subset(avg.astrocytes, gene %in% astrocytes.genes.to.label),aes(label=gene), size=4, colour="red3")



pdf("astrocytes_scatter_Fig2.pdf")
astrocytes_plot_scatter
dev.off()


Astrocytes_all_diff<- de_list$Astrocytes %>% tibble::rownames_to_column("gene") 

avg.astrocytes_dt <- avg.astrocytes  %>% tbl_df()
colnames(avg.astrocytes_dt) <- c("Brain_4C.logAvgExpr", "Brain_37C.logAvgExpr", "gene")


astro <- Astrocytes_all_diff %>% left_join(avg.astrocytes_dt, by=c("gene"))
DT::datatable(astro) 



################### Astrocytes protein and RNA expression #############


pa <-read.csv("Astrocyte_proteins.csv", header=T) %>% tbl_df()%>% mutate(Protein=as.character(Protein))

paa <- pa 

astrocyte_plot_p <- paa %>% ggplot(aes(log2_intensity_37-log2_intensity_4C, -log10(adj.pvalue)))+geom_point(size=0.7)+xlab("Log2 FC(37C/4C)")+ylab("-log10(adj.P-value)")+ggtitle("Astrocytes")

pa <- pa %>% mutate(Protein = strsplit(as.character(Protein), ";")) %>% tidyr::unnest(Protein) %>% dplyr::select(Protein, everything())


UNIPROT <- mapIds(org.Mm.eg.db, pa$Protein , "SYMBOL", 'UNIPROT') 
UNIPROT<-cbind(UNIPROT, names(UNIPROT)) 

colnames(UNIPROT)<-c("gene", "Protein")
UNIPROT <- UNIPROT %>% tbl_df()

protein_a <- pa %>% left_join(UNIPROT) 

protein_a <- protein_a %>% dplyr::select(gene, adj.pvalue, log2FC, log2_intensity_37, log2_intensity_4C) %>% group_by(gene) %>% summarise(adj.pvalue=max(adj.pvalue), log2FC=median(log2_intensity_37 - log2_intensity_4C)) %>% arrange(adj.pvalue)

pt_astro <- protein_a %>% full_join(astro)

astrocyte_top100_tr <- pt_astro  %>% arrange(p_val_adj) %>% filter(!is.na(log2FC)) %>% head(50)

pt_a <- pt_astro %>% mutate(Detected = ifelse(gene %in% astrocyte_top100_tr$gene, "Top 50 sign. RNA",  "All"))

#  pt_a %>% filter(gene %in% astrocyte_top100_tr) %>% as.data.frame()

plot.pt_a<-pt_a %>% tidyr::drop_na()%>% ggplot(aes(x=log2FC, y=-avg_logFC, color=Detected)) +geom_point() + stat_cor(method = "pearson", aes(color = Detected), size=5) + geom_hline(yintercept = 0, linetype="dashed")+ geom_vline(xintercept = 0, linetype="dashed")+ggtitle("Astrocytes")+xlab("Protein LFC: 37C vs 4C")+ylab("RNA LFC: 37C vs 4C")+ geom_text_repel(data=subset(pt_a, gene %in% astrocyte_top100_tr$gene),aes(label=gene), size=4)+xlim(-4,4)+ylim(-4,4)+ scale_color_manual(values=c("grey", "black"))

pdf("Astrocyte_prot_rna.pdf", height=8, width = 11)
plot.pt_a
dev.off()

cor(pt_a$log2FC, -pt_a$avg_logFC, use="pairwise.complete.obs")


###################### DE in Microglia #####################

microglia_s <- subset(brain.combined_clean.rn, idents = c("Microglia"))
Idents(microglia_s) <- "stim"
avg.microglia <- log1p(AverageExpression(microglia_s, verbose = FALSE)$RNA)
avg.microglia$gene <- rownames(avg.microglia)

microglia.genes.to.label <- (avg.microglia %>% filter(abs(Brain_4C-Brain_37C)>1.5))$gene



microglia_plot_scatter <- ggplot(avg.microglia, aes(Brain_4C, Brain_37C)) + geom_point() + ggtitle("Microglia")
#microglia_plot_scatter <- LabelPoints(plot = microglia_plot_scatter, points = microglia.genes.to.label, repel = TRUE)
microglia_plot_scatter<-microglia_plot_scatter + geom_text_repel(data=subset(avg.microglia, gene %in% microglia.genes.to.label),aes(label=gene), size=4, colour="red3")


pdf("microglia_scatter_Fig2.pdf")
microglia_plot_scatter
dev.off()

microglia_all_diff<- de_list$Microglia %>% tibble::rownames_to_column("gene")

avg.microglia_dt <- avg.microglia  %>% tbl_df()
colnames(avg.microglia_dt) <- c("Brain_4C.logAvgExpr", "Brain_37C.logAvgExpr", "gene")


micro <- microglia_all_diff %>% left_join(avg.microglia_dt, by=c("gene"))



pm<-read.csv("Microglia_proteins.csv", header=T) %>% tbl_df()%>% mutate(Protein=as.character(Protein))

pmm <- pm

microglia_plot_p <- pmm %>% ggplot(aes(log2_intensity_37 - log2_intensity_4C, -log10(adj.pvalue)))+geom_point(size=0.7)+xlab("Log2 FC(37C/4C)")+ylab("-log10(adj.P-value)")+ggtitle("Microglia")

pm <- pm %>% mutate(Protein = strsplit(as.character(Protein), ";")) %>% tidyr::unnest(Protein) %>% dplyr::select(Protein, everything())
UNIPROT <- mapIds(org.Mm.eg.db, pm$Protein , "SYMBOL", 'UNIPROT') 
UNIPROT<-cbind(UNIPROT, names(UNIPROT)) 
colnames(UNIPROT)<-c("gene", "Protein")
UNIPROT <- UNIPROT %>% tbl_df()
protein_m <- pm %>% left_join(UNIPROT) 
protein_m <- protein_m %>% dplyr::select(gene, adj.pvalue, log2FC, log2_intensity_37, log2_intensity_4C) %>% group_by(gene) %>% summarise(adj.pvalue=max(adj.pvalue), log2FC=median(log2_intensity_37-log2_intensity_4C)) %>% arrange(adj.pvalue)


pt_micro <- protein_m %>% full_join(micro)


microglia_top100_tr <- pt_micro  %>% arrange(p_val_adj) %>% filter(!is.na(log2FC)) %>% head(50)

pt_m <- pt_micro %>% mutate(Detected = ifelse(gene %in% microglia_top100_tr$gene, "Top 50 sign. RNA",  "All"))

pt_m %>% filter(gene %in% microglia_top100_tr) %>% as.data.frame()


plot.pt_m <- pt_m %>% tidyr::drop_na() %>% ggplot(aes(x=log2FC, y=-avg_logFC, color=Detected)) + geom_point() + stat_cor(method = "pearson", aes(color = Detected), size=5) + geom_hline(yintercept = 0, linetype="dashed")+ geom_vline(xintercept = 0, linetype="dashed")+ggtitle("Microglia")+xlab("Protein LFC: 37C vs 4C")+ylab("RNA LFC: 37C vs 4C")+ geom_text_repel(data=subset(pt_m, gene %in% microglia_top100_tr$gene),aes(label=gene), size=4)+xlim(-4,4)+ylim(-4,4) + scale_color_manual(values=c("grey", "black"))

pdf("Microglia_prot_rna.pdf", height=8, width = 11)
plot.pt_m

dev.off()

cor(pt_m$log2FC, -pt_m$avg_logFC, use="pairwise.complete.obs")

################# GO proteins ###############

#The enrichment tests are carried out on Gorilla GO web application 

Go_micro_prot<-gdata::read.xls("GO_Proteomics.xlsx", sheet=1)
go_micro_p_plot<-Go_micro_prot %>% ggplot(aes(x=b/B, y=Description, color=-log10(FDR.q.value), size=b, shape=in.37C))+geom_point()+ ggplot2::scale_color_gradient(high = "purple", low = "darkred") + ggplot2::theme_bw() + ggplot2::labs(size = "# annotated genes", y = "GO term", color = "-log10(FDR q-value)", x = "gene-ratio", title = "Microglia", shape="in 37C")
go_micro_p_plot

Go_astro_prot<-gdata::read.xls("GO_Proteomics.xlsx", sheet=2)
go_astro_p_plot<-Go_astro_prot %>% ggplot(aes(x=b/B, y=Description, color=-log10(FDR.q.value), size=b, shape=in.37C))+geom_point()+ ggplot2::scale_color_gradient(high = "purple", low = "darkred") + ggplot2::theme_bw() + ggplot2::labs(size = "# annotated genes", y = "GO term", color = "-log10(FDR q-value)", x = "gene-ratio", title = "Astrocytes", shape="in 37C")
go_astro_p_plot

pdf("results/go_proteomics.pdf", height=6, width=20)
plot_grid(plotlist = list(go_micro_p_plot, go_astro_p_plot))
dev.off()

################# RNA GO terms with TMOD ################
#similar GO terms are enriched when data analyzed with Gorilla


require ('tmod')

tmod_dbs <- list(MSigDB=readRDS("data/msigdb_HS_v6.2.rds"))
tmod_mapping <- list(
  default=feature_annot %>% dplyr::select(from_gene_id=gene_id, to_gene_id=gene_name)
)
tmod_mapping$default$to_gene_id = toupper(tmod_mapping$default$to_gene_id)

msigdb <- tmod_dbs$MSigDB


go <- msigdb[ msigdb$MODULES$Category == "C5" ]


neurons<- de_list$Neurons %>% tibble::rownames_to_column("Gene") %>% mutate(Gene=toupper(Gene))
pyramidal <- de_list$Pyramidal %>% tibble::rownames_to_column("Gene") %>% mutate(Gene=toupper(Gene))
microglia <- de_list$Microglia %>% tibble::rownames_to_column("Gene") %>% mutate(Gene=toupper(Gene))#%>% left_join(mouse_human, by=c("Gene"="Mm_symbol"))
astrocytes<- de_list$Astrocytes %>% tibble::rownames_to_column("Gene")  %>% mutate(Gene=toupper(Gene))#%>% left_join(mouse_human, by=c("Gene"="Mm_symbol"))
DG <- de_list$DG %>% tibble::rownames_to_column("Gene")  %>% mutate(Gene=toupper(Gene))# %>% left_join(mouse_human, by=c("Gene"="Mm_symbol"))
oligo <- de_list$Oligodendr. %>% tibble::rownames_to_column("Gene") %>% mutate(Gene=toupper(Gene)) #%>% left_join(mouse_human, by=c("Gene"="Mm_symbol"))
endothelial <- de_list$Endothelial %>% tibble::rownames_to_column("Gene") %>% mutate(Gene=toupper(Gene))
npc <- de_list$NPC %>% tibble::rownames_to_column("Gene") %>% mutate(Gene=toupper(Gene))
bam <- de_list$BAM %>% tibble::rownames_to_column("Gene") %>% mutate(Gene=toupper(Gene))
mural <- de_list$Mural %>% tibble::rownames_to_column("Gene") %>% mutate(Gene=toupper(Gene))
polydendrocytes <- de_list$Polydendr. %>% tibble::rownames_to_column("Gene") %>% mutate(Gene=toupper(Gene))
fibroblast <- de_list$`Fibroblast-like` %>% tibble::rownames_to_column("Gene") %>% mutate(Gene=toupper(Gene))


pie_list<-list()

pie_list$oligo   <- tmodDecideTests(oligo$Gene, pval=oligo$p_val_adj, lfc=oligo$avg_logFC, lfc.thr=0.5, pval.thr=0.01, mset=go, labels="oligo")
pie_list$microglia   <- tmodDecideTests(microglia$Gene, pval=microglia$p_val_adj, lfc=microglia$avg_logFC, lfc.thr=0.5, pval.thr=0.01, mset=go, labels="microglia")
pie_list$astrocytes   <- tmodDecideTests(astrocytes$Gene, pval=astrocytes$p_val_adj, lfc=astrocytes$avg_logFC, lfc.thr=0.5, pval.thr=0.01, mset=go, labels="astrocytes")
pie_list$bam   <- tmodDecideTests(bam$Gene, pval=bam$p_val_adj, lfc=bam$avg_logFC, lfc.thr=0.5, pval.thr=0.01, mset=go, labels="bam")
pie_list$endothelial   <- tmodDecideTests(endothelial$Gene, pval=endothelial$p_val_adj, lfc=endothelial$avg_logFC, lfc.thr=0.5, pval.thr=0.01, mset=go, labels="endothelial")
pie_list$npc   <- tmodDecideTests(npc$Gene, pval=npc$p_val_adj, lfc=npc$avg_logFC, lfc.thr=0.5, pval.thr=0.01, mset=go, labels="npc")
pie_list$mural   <- tmodDecideTests(mural$Gene, pval=mural$p_val_adj, lfc=mural$avg_logFC, lfc.thr=0.5, pval.thr=0.01, mset=go, labels="mural")
pie_list$fibroblast   <- tmodDecideTests(fibroblast$Gene, pval=fibroblast$p_val_adj, lfc=fibroblast$avg_logFC, lfc.thr=0.5, pval.thr=0.01, mset=go, labels="fibroblast")
pie_list$neurons   <- tmodDecideTests(neurons$Gene, pval=neurons$p_val_adj, lfc=neurons$avg_logFC, lfc.thr=0.5, pval.thr=0.01, mset=go, labels="neurons")

tmod_list<-list()

tmod_list$oligo <-tmodUtest(l=toupper(oligo$Gene), mset=go, qval=0.1) %>% filter(AUC>0.75)
tmod_list$microglia <- tmodUtest(l=toupper(microglia$Gene), mset=go, qval=0.1)%>% filter(AUC>0.75)
tmod_list$astrocytes<-tmodUtest(l=toupper(astrocytes$Gene), mset=go, qval=0.1)%>% filter(AUC>0.75)
tmod_list$npc<-tmodUtest(l=toupper(npc$Gene), mset=go, qval=0.1)%>% filter(AUC>0.75)
tmod_list$bam <-tmodUtest(l=toupper(bam$Gene), mset=go, qval=0.1)%>% filter(AUC>0.75)
tmod_list$endothelial <- tmodUtest(l=toupper(endothelial$Gene), mset=go, qval=0.1)%>% filter(AUC>0.75)
tmod_list$mural <- tmodUtest(l=toupper(mural$Gene), mset=go, qval=0.1)%>% filter(AUC>0.75)
tmod_list$fibroblast <- tmodUtest(l=toupper(fibroblast$Gene), mset=go, qval=0.1)%>% filter(AUC>0.75)
tmod_list$neurons <- tmodUtest(l=toupper(neurons$Gene), mset=go, qval=0.1)%>% filter(AUC>0.75)



n <- names(pie_list)
pie_list <- unlist(pie_list, recursive=FALSE)
names(pie_list) <- n



main_pie<-pie_list
main_pie$oligo <- NULL
main_pie$endothelial <- NULL
main_pie$fibroblast <- NULL
main_pie$bam <- NULL
main_pie$mural <- NULL
main_pie$npc <- NULL
sup_pie<-pie_list
sup_pie$microglia<-NULL
sup_pie$neurons<-NULL
sup_pie$astrocytes<-NULL


main_tmod<-tmod_list
main_tmod$oligo <- NULL
main_tmod$endothelial <- NULL
main_tmod$fibroblast <- NULL
main_tmod$bam <- NULL
main_tmod$mural <- NULL
main_tmod$npc <- NULL
sup_tmod<-tmod_list
sup_tmod$microglia<-NULL
sup_tmod$neurons<-NULL
sup_tmod$astrocytes<-NULL


tmod_modules_main <-c("M11500","M17781","M17827", "M18911","M15503", "M18505","M11272", "M12750", "M12919","M11077", "M18023","M18135") 
tmod_modules_sup<- unique(c("M13134", "M17218","M15151","M18363","M16669","M18911","M12388", "M17218","M17360","M12919","M12259","M11055","M12120","M12089", "M15478", "M11036", "M16429"))

tmod_modules<-unique(tmod_modules)



pdf("tmod.pdf", height=14, width = 24)
tmodPanelPlot(x=tmod_list, pie=pie_list,  select=tmod_modules, pie.colors=c("royalblue3", "grey", "red3"), filter.rows.pval = 0.1, pval.thr = 10^-1)
dev.off()

pdf("tmod_Fig2.pdf", height=14, width = 24)
tmodPanelPlot(x=main_tmod, pie=main_pie,  select=tmod_modules_main, pie.colors=c("royalblue3", "grey", "red3"), filter.rows.pval = 0.1, pval.thr = 10^-1, text.cex = 2)
dev.off()


pdf("tmod_suppl.pdf", height=14, width = 24)
tmodPanelPlot(x=sup_tmod, pie=sup_pie,  select=tmod_modules_sup, pie.colors=c("royalblue3", "grey", "red3"), filter.rows.pval = 0.1, pval.thr = 10^-1, text.cex = 2)
dev.off()



cex=1.9

pdf("evidence_plots_Fig2.pdf", height=15, width=10)
par(mfrow=c(3, 2))

tmod::evidencePlot(toupper(neurons$Gene), mset=go, m="M18023", gene.labels=TRUE, col="royalblue3", col.main="royalblue3", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex, main=" cytoskeletal adaptor activity")
tmod::evidencePlot(toupper(neurons$Gene), mset=go, m="M18135", gene.labels=TRUE, col="royalblue3", col.main="royalblue3", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex, main="atpase activity coupled to transmembr.\n movement of ions phosphorylative mech.")

tmod::evidencePlot(toupper(astrocytes$Gene), mset=go, m="M12388", gene.labels=TRUE, col="goldenrod1", col.main="goldenrod1", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex, main="protein targeting to membrane")
tmod::evidencePlot(toupper(astrocytes$Gene), mset=go, m="M11077", gene.labels=TRUE, col.main="goldenrod1", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex, col="goldenrod1", main="translational initiation")

tmod::evidencePlot(toupper(microglia$Gene), mset=go, m="M11209", gene.labels=TRUE, col="seagreen3", col.main="seagreen3", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex, main="establishment of protein localization\n to endoplasmic reticulum")
tmod::evidencePlot(toupper(microglia$Gene), mset=go, m="M17827", gene.labels=TRUE, col="seagreen3", col.main="seagreen3", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex, main="anchored component of membrane")

dev.off()



pdf("evidence_plots_supp.pdf", height=15, width=20)
par(mfrow=c(3, 4))


tmod::evidencePlot(toupper(oligo$Gene), mset=go, m="M12919", gene.labels=TRUE,col="firebrick2", col.main="firebrick2", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex, main="oxidative phosphorylation")
tmod::evidencePlot(toupper(oligo$Gene), mset=go, m="M17114", gene.labels=TRUE, col="firebrick2", col.main="firebrick2", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex, main="cytosolic part")

tmod::evidencePlot(toupper(endothelial$Gene), mset=go, m="M12919", gene.labels=TRUE, col="cadetblue4", col.main="cadetblue4", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex, main="oxidative phosphorylation")
tmod::evidencePlot(toupper(endothelial$Gene), mset=go, m="M17339", gene.labels=TRUE, col="cadetblue4", col.main="cadetblue4", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex, main="cytosolic large ribosomal subunit")

tmod::evidencePlot(toupper(mural$Gene), mset=go, m="M12040", gene.labels=TRUE, col="peru", col.main="peru", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex, main="protein localization to endoplasmic reticulum")
tmod::evidencePlot(toupper(mural$Gene), mset=go, m="M12331", gene.labels=TRUE, col="peru", col.main="peru", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex, main="glucan metabolic process")

tmod::evidencePlot(toupper(fibroblast$Gene), mset=go, m="M18422", gene.labels=TRUE, col="darkorchid4", col.main="darkorchid4", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex, main="structural constituent of ribosome")
tmod::evidencePlot(toupper(fibroblast$Gene), mset=go, m="M17302", gene.labels=TRUE, col="darkorchid4", col.main="darkorchid4", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex,main="basement membrane")


tmod::evidencePlot(toupper(bam$Gene), mset=go, m="M14922", gene.labels=TRUE, col="yellow3", col.main="yellow3", lwd=cex, gl.cex=cex, cex.axis=cex, cex.main=cex, cex.lab=cex, main="multi organism metabolic process")
tmod::evidencePlot(toupper(bam$Gene), mset=go, m="M17781", gene.labels=TRUE, col="yellow3", col.main="yellow3", lwd=cex, gl.cex=cex, cex.axis=cex, cex.main=cex, cex.lab=cex, main="cytosolic ribosome")

tmod::evidencePlot(toupper(npc$Gene), mset=go, m="M18911", gene.labels=TRUE, col="plum4", col.main="plum4", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex, main="extracellular ligand gated\n ion channel activity")
tmod::evidencePlot(toupper(npc$Gene), mset=go, m="M18276", gene.labels=TRUE, col="plum4", col.main="plum4", lwd=cex, gl.cex=cex, cex.axis=cex, cex.lab=cex, cex.main=cex, main="glutamate receptor activity")


dev.off()


