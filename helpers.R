##########################################################################
# Copyright © 2024 Daniel Zimmermann. <daniel.zimmermann.dz@outlook.com> #
#                                                                        #
# You may use, distribute and modify this code under the terms of the    #
# MIT license.                                                           #
#                                                                        #
# You should have received a copy of the MIT license with this file.     #
# If not, please visit: https://opensource.org/license/mit               #
##########################################################################

libs <- c("DESeq2", "dplyr", "magrittr", "tidyverse", "RColorBrewer",
          "stringr", "ggplot2", "ggrepel", "ggpubr", "org.Hs.eg.db", 
          "GenomicTools.fileHandler", "ComplexHeatmap", "circlize", "plotly")

lapply(libs, library, character.only = TRUE)

read_counts <- function(filenames, data_path, sample_IDs) {
  countdata <- filenames %>%
    file.path(data_path, .) %>%
    lapply(importFeatureCounts) %>% 
    lapply(function(x) x[[1]]) %>%
    reduce(function(x, y) full_join(x, y, by="Geneid")) %>%
    set_colnames(c("Geneid", sample_IDs))
  
  countdata
}

read_genes <- function(filepath, genes) {
  genedata <- importGTF(filepath, 
                        level = "transcript", 
                        features = c("gene_id", 
                                     "gene_biotype", 
                                     "gene_name",
                                     "gene_source"),
                        merge.all = T)
  
  genedata %>%
    as_tibble() %>%
    dplyr::select(starts_with("gene")) %>%
    group_by(gene_id) %>%
    summarise(across(everything(), first)) %>%
    arrange(match(gene_id, genes)) %>%
    set_colnames(c("ENSEMBL", "Biotype", "Name", "Source"))
}

filter_counts <- function(countdata) {
  countdata %>% 
    # Remove genes with all zero counts
    filter(rowSums(across(where(is.numeric))) > 0) %>% 
    # Merge duplicate entries
    group_by(Geneid) %>%
    summarise(across(!where(is.numeric), ~ dplyr::first(.)), across(where(is.numeric), ~ sum(.))) %>%
    ungroup() %>%
    as.data.frame() %>%
    column_to_rownames("Geneid")
}


plot_pca <- function(data, n = 500, maxPC = 2, center = TRUE, scale = TRUE, plot_all = F) {
  samples <- colData(data) %>% 
    as.data.frame()
  
  
  data.pca <- assay(data) %>%
    t() %>%
    prcomp(center = center, scale. = scale) 
  
  data.pca.plot <- data.pca$x %>%
    as.data.frame() %>%
    rownames_to_column("Sample_ID") %>%
    left_join(samples, by = "Sample_ID")
  
  data.pca.var <- as.data.frame(t(summary(data.pca)$importance)) %>%
    set_colnames(c("sd", "var_prop", "var_cumprop")) %>%
    mutate(across(c(var_prop, var_cumprop), ~ .x * 100), Component = seq.int(1, nrow(.), 1))
  
  p1 <- ggplot(data.pca.var, aes(x = Component)) +
    geom_bar(aes(y = var_prop), stat = "identity") +
    geom_line(aes(y = var_cumprop)) +
    geom_point(aes(y = var_cumprop)) +
    scale_x_continuous(limits = c(0.5, nrow(data.pca.var) + 0.5), 
                       breaks = seq(1, nrow(data.pca.var), by = 1),
                       expand = expansion(add=0)) +
    scale_y_continuous(expand = expansion(mult=c(0, 0.05))) +
    ylab("Variance (%)") +
    theme_pubr() +
    theme(panel.grid.major.y = element_line(linewidth = 0.5))
  
  print(p1)
  
  
  p2 <- ggplot(data.pca.plot, aes(x = PC1, y = PC2, color = Class)) +
    geom_point(size = 2) +
    labs(color = "") +
    theme_pubr() +
    theme(panel.grid.major = element_line(linewidth = 0.5),
          legend.position = "right")
  
  print(p2)
  
  if (plot_all) {
    
    p3 <- ggplot(data.pca.plot, aes(color = Class)) +
      geom_autopoint(size = 2) +
      geom_autodensity(position = "identity", fill = NA) +
      facet_matrix(vars(paste0("PC", seq.int(1, maxPC, 1))), layer.diag = 2) +
      labs(color = "") +
      theme_bw() +
      theme(panel.grid.major = element_line(linewidth = 0.5))
    
    print(p3)
  }
}


plot_heatmap <- function(data, color_scale, column_colors, genedata, samples, 
                         pvalues, title = "Rel. Expression", cluster_rows = TRUE, 
                         cluster_cols = TRUE, max_rows = 20,scale_by_row = TRUE) {
  
  annotation.row <- genedata %>%
    filter(ENSEMBL %in% rownames(data)) 
  
  annotation.col <- samples %>%
    filter(Sample_ID %in% colnames(data)) %>%
    dplyr::select(c(Class))
  
  annotation_labels <- unique(annotation.col$Class)
  annotation_colors <- column_colors[annotation_labels]
  
  if (scale_by_row) {
    data %<>%
      t() %>%
      scale() %>%
      t()
  }
  
  if (nrow(data) <= max_rows) {
    hm <- Heatmap(as.matrix(data), col = color_scale,  
                  name = title, cluster_rows = cluster_rows, column_split = annotation.col$class,
                  show_row_dend = TRUE, show_column_names = FALSE, row_labels = annotation.row$GeneSymbol,
                  show_row_names = TRUE, row_title = NULL, column_title = NULL, cluster_columns = cluster_cols, 
                  top_annotation = HeatmapAnnotation(class = anno_block(labels = annotation_labels,
                                                                        gp = gpar(fill = annotation_colors))),
                  heatmap_legend_param = list(legend_height = unit(10, "cm"),
                                              legend_width = unit(1.5, "cm"),
                                              title_position = "lefttop-rot"))
  
    draw(hm)
    
  } else {
    hm <- Heatmap(as.matrix(data), col = color_scale,  
                  name = title, cluster_rows = cluster_rows, column_split = annotation.col$Class, 
                  show_row_dend = FALSE, show_column_names = FALSE, show_row_names = FALSE, 
                  row_title = NULL, column_title = NULL, cluster_columns = cluster_cols, 
                  top_annotation = HeatmapAnnotation(Class = anno_block(labels = annotation_labels,
                                                                        gp = gpar(fill = annotation_colors))),
                  heatmap_legend_param = list(legend_height = unit(10, "cm"),
                                              legend_width = unit(1.5, "cm"),
                                              title_position = "lefttop-rot"))
    
    draw(hm)
    
    data.top <- data %>% 
      as.data.frame() %>%
      rownames_to_column("ENSEMBL") %>%
      left_join(data.frame(pvalues, ENSEMBL = names(pvalues)), by = "ENSEMBL") %>%
      column_to_rownames("ENSEMBL") %>%
      arrange(pvalues) %>%
      dplyr::select(!pvalues) %>%
      head(max_rows) 
    
    annotation.row.top <- annotation.row %>%
      filter(ENSEMBL %in% rownames(data.top))
    
    hm <- Heatmap(as.matrix(data.top), col = color_scale,  
                  name = title, cluster_rows = cluster_rows, column_split = annotation.col$Class, 
                  show_row_dend = TRUE, show_column_names = FALSE, row_labels = annotation.row.top$Name,
                  show_row_names = TRUE, row_title = NULL, column_title = NULL, cluster_columns = cluster_cols, 
                  top_annotation = HeatmapAnnotation(Class = anno_block(labels = annotation_labels,
                                                                        gp = gpar(fill = annotation_colors))),
                  heatmap_legend_param = list(legend_height = unit(10, "cm"),
                                              legend_width = unit(1.5, "cm"),
                                              title_position = "lefttop-rot"))
    
    draw(hm)
    
  }
}


volcanoplot <- function(data, title = NA, max_labels = NA, lfc_limit = NA) {
  data.plot <- data %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    mutate(genelabels = "", 
           out_of_bounds = (ifelse(is.na(vp_lfc_limit), 0, abs(log2FoldChange) > vp_lfc_limit) * sign(log2FoldChange)),
           log2FoldChange_capped = ifelse(out_of_bounds != 0, vp_lfc_limit * sign(log2FoldChange), log2FoldChange),
           padj = ifelse(padj == 0, .Machine$double.xmin, padj))
  
  if (is.na(max_labels) || sum(data.plot$threshold) < max_labels) {
    data.plot$genelabels[data.plot$threshold] <- data.plot$ENSEMBL[data.plot$threshold]
  } else if (max_labels > 0) {
    data.plot$genelabels[data.plot$threshold][1:max_labels] <- data.plot$ENSEMBL[data.plot$threshold][1:max_labels]
  }
  
  maxFC <- max(abs(data.plot$log2FoldChange))
  if (!is.na(lfc_limit) && lfc_limit < maxFC) {
    xlim <- lfc_limit
  } else {
    xlim <- maxFC * 1.04
  }
  
  breaks_y <- 10^-seq.int(0, -floor(log10(min(data.plot$padj))), by=20)
  
  
  ggplot(data.plot, aes(x = log2FoldChange_capped, y = padj)) +
    geom_point(data = subset(data.plot, out_of_bounds == 0), aes(colour=threshold), alpha = 0.5) +
    geom_point(data = subset(data.plot, out_of_bounds == -1), aes(colour=threshold), shape = "\u25c4", size=2) +
    geom_point(data = subset(data.plot, out_of_bounds == 1), aes(colour=threshold), shape = "\u25BA", size=2) +
    geom_hline(yintercept = 0.05, linetype = 2) +
    geom_vline(xintercept = c(-1, 1), linetype = 2) +
    geom_text_repel(aes(label = genelabels)) +
    scale_x_continuous(breaks = scales::pretty_breaks(), limits = c(-xlim, xlim), expand = expansion(0.01)) +
    scale_y_continuous(trans = c("log10", "reverse"), breaks = breaks_y) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "darkred")) +
    ggtitle(title) +
    xlab("log2 Fold Change") +
    ylab("adjusted p-value") +
    theme_pubr() +
    theme(legend.position = "none")
}


volcanoplot_interactive <- function(data, file = NA, alpha = 0.05, minFC = 1, title = NA) {
  data.plot <- data %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj),
           threshold = padj < alpha & abs(log2FoldChange) >= minFC)
  
  
  p <- plot_ly(data = data.plot, x = ~log2FoldChange, y = ~padj, color = ~threshold, text = ~Name, 
               opacity = 0.6, type = "scatter", mode = "markers", colors = c("black", "darkred"),
               hoverinfo = "none",
               hovertemplate = paste("<b>Gene:</b> %{text}",
                                     "<br><b>Log2 Fold Change:</b> %{x:.3r}",
                                     "<br><b>Adj. p-Value:</b> %{y:.2e}<extra></extra>")) %>%
    layout(title = list(text = title),
           xaxis = list(title = "Log2 Fold Change", zeroline = F),
           yaxis = list(title = "Adj. p-Value", type = "log", autorange="reversed",
                        exponentformat = "power", showexponent = "all"),
           shapes = list(list(type = "line",
                              x0 = 0,
                              x1 = 1,
                              xref = "paper",
                              y0 = alpha,
                              y1 = alpha,
                              line = list(color = "black", width = 1, dash = "dash")),
                         list(type = "line",
                              x0 = -minFC,
                              x1 = -minFC,
                              y0 = 0,
                              y1 = 1,
                              yref = "paper",
                              line = list(color = "black", width = 1, dash = "dash")),
                         list(type = "line",
                              x0 = minFC,
                              x1 = minFC,
                              y0 = 0,
                              y1 = 1,
                              yref = "paper",
                              line = list(color = "black", width = 1, dash = "dash"))),
           hoverlabel=list(bgcolor="white"),
           showlegend = F,
           margin = list(l = 50, r = 50, b = 50, t = 50, pad = 20)) %>%
    toWebGL()
  if (is.na(file)) {
    p
  } else {
    p %>% htmlwidgets::saveWidget(file, selfcontained = T)
  }
}
