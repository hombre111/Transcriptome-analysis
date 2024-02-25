
library(tidyverse)
library(BiocManager)
library(edgeR)
library(devtools)
library(annotables)
library(ggrepel)
library(pheatmap)
library(ggalt)
library(viridis)

#### Data loading and cleanup #######

dir <- "_____SET_DIR____"
setwd(dir)

cleaned_raw <- read.csv("_____FILE_NAME____", header = TRUE, row.names = 1)

col_names<- c("C1","C2","C3","C4", "DHA1","DHA2","DHA3","DHA4","X1", "X2","X3","X4")

colnames(cleaned_raw)[1:12] <- col_names

cleaned_raw <- cleaned_raw |> 
  relocate(starts_with("DHA"), .after = "X4")



#Normalisation

group <- factor(c(rep("A-C", 4), rep("B-X", 4), rep("C-DHA",4)))
lecone  <- DGEList(counts = cleaned_raw, group = group)

norma <- cleaned_raw[filterByExpr(lecone), ]
sumy_sampli <- colSums(norma)
norma_matryca <- apply(norma, 1, FUN=function(X) log((X+1)/(sumy_sampli+1)))
scaled_norma_matryca <- scale(norma_matryca)
dim(scaled_norma_matryca)

#PCA
astro_pca <- prcomp(scaled_norma_matryca)


## linear model to explore factors associated with PCs

pca.df <- as.data.frame(astro_pca$x)
pca.df$Sample <- rownames(pca.df)
meta.df <- data.frame("Sample"=rownames(pca.df),
                      "Condition"=gsub(rownames(pca.df), 
                                       pattern="([A-Za-z]+)([0-9]+)",
                                       replacement="\\1"),
                      "Culture"= gsub(rownames(pca.df), 
                                      pattern="([A-Za-z]+)([0-9]+)",
                                      replacement="\\2"),
                      "Time"= c(51,51,21,31,51,51,21,31,51,51,21,31))

pca.meta.df <- merge(pca.df, meta.df, by="Sample")

pvalues_time <- numeric()
pvalues_con <- numeric()

for (i in as.vector(colnames(pca.meta.df[,2:13]))) {
  formula_str <- paste(i, "~ Time")
  formula_obj <- as.formula(formula_str)
  lmp <- lm(formula_obj, data = pca.meta.df)
  coef_summary <- summary(lmp)
  pvalues_time <- c(pvalues_time, coef_summary$coefficients[2, 4])
}

plot(log10(pvalues_time))



#PC variance scree plot


PC_components<- colnames(data.frame(astro_pca$x))
variance_explained <- data.frame(Variance_pc = round(astro_pca$sdev^2 / sum(astro_pca$sdev^2) * 100, 2),
                                 PC = PC_components)
variance_explained$PC <- factor(variance_explained$PC, levels = PC_components)
table(variance_explained)

ggplot(variance_explained, aes(x = PC, y = Variance_pc, fill = PC))+
  geom_col(linetype = 1, linewidth = 20)+
  theme_bw()+
  scale_fill_viridis(discrete = TRUE, option = "E", direction = 1)+
  labs(title = "Percentage of variance explained by each PC", y = "Variance explained %")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#PCA-Score plot


astro_pca_data <- as.data.frame(astro_pca$x)
head(astro_pca_data)
astro_pca_data$Sample <- rownames(astro_pca)

custom_palette <- c("Non-stimulated" = "#21908CFF", "Stimulated" = "#440154FF", "Stimulated+DHA" = "#FFA500")

astro_score_data <- data.frame(Sample = rownames(astro_pca$x), 
                               PC1 = astro_pca$x[, 1], 
                               PC2 = astro_pca$x[, 2])

color_control <- c("C1","C2","C3","C4")
color_X <- c("X1","X2","X3","X4")
color_DHA <- c("DHA1","DHA2","DHA3","DHA4")

astro_score_data$Condition <- ifelse(
  astro_score_data$Sample %in% color_control, 
  "Non-stimulated", 
  ifelse(
    astro_score_data$Sample %in% color_X,
    "Stimulated",
    "Stimulated+DHA"
  )
)


group_summary <- data.frame(
  group_x = c(-70, 70, 5),
  group_y = c(0, 10, 20),
  Condition = c("Non-stimulated", "Stimulated", "Stimulated+DHA") 
)

score_plot<- ggplot(astro_score_data, 
                    aes(x = PC1, y = PC2, fill = Condition)) +
  geom_point(pch = 21, 
             size = 4, 
             alpha = 0.8,
             color = "black", 
             aes(colour = Condition)) +
  scale_color_manual(values = custom_palette)+
  scale_fill_manual(values = custom_palette)+
  geom_encircle(inherit.aes = TRUE, 
                alpha = 0.1, 
                show.legend = FALSE, 
                expand = 0.15, 
                spread = 0.2)+
  labs(title = "Score Plot", x = "PC1: 42.85% variance", y = "PC2: 21.93% variance") +
  theme_test()+
  theme(plot.title = element_text(face = "bold"))+
  geom_text(data = group_summary, aes(x = group_x, y = group_y, label = Condition), vjust = -0.5)

score_plot


#Correlation matrix of samples based on PCA results


corelation_matrix<- cor(t(astro_pca_data))
corelation_matrix

cor_color_palette <- magma(100)
gradient_palette <- colorRampPalette(c("#21908CFF", "#240154FF", "#FFA500"))(n = 100)

pheatmap(corelation_matrix, color = gradient_palette, border_color = FALSE)



# Differential gene expression Analysis with edgeR
#Creating design matrix

jazda <- cleaned_raw

group <- factor(c(rep("A-C", 4), rep("B-X", 4), rep("C-DHA",4)))
culture_time <- factor(rep(c(51,51,21,31),3))

design_matrix <- data.frame("Condition" = group,
                            "Time" = culture_time)
row.names(design_matrix) <- colnames(jazda)

lecone  <- DGEList(counts = jazda, group = group)

matrix <- cbind(lecone$samples, design_matrix)
matrix <- matrix[, c(4,5)]

condition <- factor(design_matrix$Condition)
time <- factor(design_matrix$Time)

model_matrix <- model.matrix(~time+condition)


# Filtering and normalisation
zatrzymaj <- filterByExpr(lecone)
filtered_counts <- lecone[zatrzymaj, ]

normalised <- normLibSizes(filtered_counts) 
normalised$design <- model_matrix


###Estimate dispersion
dge <- estimateDisp(normalised, design = normalised$design)

count_mean_rows <- rowMeans(dge$counts)
disp <- dge$tagwise.dispersion
trended_disp <- dge$trended.dispersion
data_frame_eh <- data.frame("Mean_counts" = count_mean_rows, "Disp" = disp, "Trended" = trended_disp)

### Plot dispersion vs counts
ggplot(data_frame_eh, aes(x = log10(Mean_counts), y = log10(Disp)))+
  geom_point(aes(alpha = 0.01))+
  geom_line(aes(x = log10(Mean_counts), y = log10(Trended), color = "red"), linewidth = 1.2)+
  labs(x = "log10 mean gene counts", y = "log10 tagwise dispersion")


### Run DGE 
fit <- glmQLFit(dge, dge$design)

xvsc <- data.frame(topTags(glmQLFTest(fit, coef = 4), n = 50000))
dhavsx <- data.frame(topTags(glmQLFTest(fit, contrast = c(0,0,0,-1,1)), n = 50000))


#### Annotate genes
#Find symbols for top genes
listaaa_genow_all<- annotables::rnor6 |> 
  filter(ensgene %in% row.names(xvsc))

listaaa_genow_all <- listaaa_genow_all[!duplicated(listaaa_genow_all$ensgene), ]
xvsc <- xvsc[row.names(xvsc) %in% listaaa_genow_all$ensgene, ]
xvsc$ensgene <- row.names(xvsc)
xvsc <- merge(xvsc, listaaa_genow_all[,c("ensgene", "symbol")], by.x = "ensgene", all.y = TRUE)



lista_genow_dha <-  annotables::rnor6 |> 
  filter(ensgene %in% row.names(dhavsx))

lista_genow_dha <- lista_genow_dha[!duplicated(lista_genow_dha$ensgene), ]
dhavsx <- dhavsx[row.names(dhavsx) %in% lista_genow_dha$ensgene, ]
dhavsx$ensgene <- row.names(dhavsx)
dhavsx <- merge(dhavsx, lista_genow_dha[,c("ensgene", "symbol")], by.x = "ensgene", all.y = TRUE)


#Volcano plot - TIC vs C

xvsc$diffexpressed <- "Not regulated"
xvsc$diffexpressed[xvsc$logFC > 0.6 & xvsc$FDR < 0.05] <- "Upregulated"
xvsc$diffexpressed[xvsc$logFC < -0.6 & xvsc$FDR < 0.05] <- "Downregulated"

labels <- xvsc %>% 
  arrange(desc(logFC)) %>% 
  select("symbol")

labels <- rbind(
  head(labels, 10),
  tail(labels, 10)
)

xvsc$delabel <- ifelse(xvsc$symbol %in% labels$symbol, xvsc$symbol, NA)


ggplot(xvsc, aes(x = logFC, y = abs(log10(FDR)), col = diffexpressed))+
  geom_point(alpha = 0.5, size = 2)+
  labs(y ="-log10 adjusted pvalue", legend = FALSE)+
  geom_vline(xintercept=c(-0.6, 0.6), col="grey", linetype = "longdash") +
  geom_hline(yintercept=-log10(0.05), col="grey", linetype = "longdash") + 
  scale_color_manual(values=c("#21908CFF", "grey", "#FFA500"))+
  geom_text_repel(aes(label = delabel), color = "black")+
  theme_bw() +
  theme(axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))+
  guides(color = guide_legend(title = NULL))



#Volcano plot - DHA vs TIC

dhavsx$diffexpressed <- "Not regulated"
dhavsx$diffexpressed[dhavsx$logFC > 0.6 & dhavsx$FDR < 0.05] <- "Upregulated"
dhavsx$diffexpressed[dhavsx$logFC < -0.6 & dhavsx$FDR < 0.05] <- "Downregulated"

labels_dha <- dhavsx %>% 
  arrange(desc(logFC)) %>% 
  select("symbol")

labels_dha <- rbind(
  head(labels_dha, 10),
  tail(labels_dha, 10)
)

dhavsx$delabel <- ifelse(dhavsx$symbol %in% labels_dha$symbol, dhavsx$symbol, NA)


ggplot(dhavsx, aes(x = logFC, y = abs(log10(FDR)), col = diffexpressed))+
  geom_point(alpha = 0.5, size = 2)+
  labs(y ="-log10 adjusted pvalue", legend = FALSE)+
  geom_vline(xintercept=c(-0.6, 0.6), col="grey", linetype = "longdash") +
  geom_hline(yintercept=-log10(0.05), col="grey", linetype = "longdash") + 
  scale_color_manual(values=c("#21908CFF", "grey", "#FFA500"))+
  geom_text_repel(aes(label = delabel), color = "black")+
  theme_bw() +
  theme(axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))+
  guides(color = guide_legend(title = NULL))




#### Heatmap of top genes 

scaled_norma_matryca <- data.frame(t(scaled_norma_matryca))
scaled_norma_matryca$ensgene <- row.names(scaled_norma_matryca)


scaled_norma_matryca <- scaled_norma_matryca[row.names(scaled_norma_matryca) %in% listaaa_genow_all$ensgene, ]
scaled_genes_symbols <- merge(scaled_norma_matryca, listaaa_genow_all[,c("ensgene", "symbol")], .by= "ensgene" )


##### TIC vs C heatmap 

symbols_xvsc <- rbind(
  head(xvsc %>% 
         arrange( desc(logFC)), 50),
  tail(xvsc %>% 
         arrange( desc(logFC)), 50)
)

xvsc_heatmap <- scaled_genes_symbols %>% 
  filter(ensgene %in% symbols_xvsc$ensgene)

xvsc_heatmap <- xvsc_heatmap[xvsc_heatmap$symbol != "", ]
rownames(xvsc_heatmap) <- xvsc_heatmap$symbol

gradient_palette <- colorRampPalette(c("#21908CFF", "#240154FF", "#FFA500"))(n = 100)

annotation_colorz <- list(Condition = c(Control = "#85908CFF", Xtail = "blue", DHA = "#FFA500"))
sample_group <- data.frame(Condition = rep(c("Control", "Xtail", "DHA"), c(4,4,4)))

head(xvsc_heatmap[,2:13])

pheatmap(xvsc_heatmap[,2:13], 
         color = gradient_palette,
         show_rownames = TRUE,
         border_color = FALSE,
         cluster_cols = FALSE, 
         annotation_col = sample_group,
         annotation_colors = annotation_colorz,
         cutree_rows = 1,
         scale = "none")



##### DHA vs TIC heatmap 
dha_heatmap <-  scaled_genes_symbols %>% 
  filter(ensgene %in% dhavsx[dhavsx$FDR < 0.05, ]$ensgene) 

dha_heatmap <- dha_heatmap[dha_heatmap$symbol != "", ]
rownames(dha_heatmap) <- dha_heatmap$symbol


pheatmap(dha_heatmap[,2:13], 
         color = gradient_palette,
         show_rownames = TRUE,
         border_color = FALSE,
         cluster_cols = FALSE, 
         annotation_col = sample_group,
         annotation_colors = annotation_colorz,
         cutree_rows = 1,
         scale = "none")




