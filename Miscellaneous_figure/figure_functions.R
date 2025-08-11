################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
################################################################################
# LIST OF REQUIRED PACKAGES
pkgs <- c("tidyverse", "data.table", "BiocManager", "here",
          "ggplot2", "ggrepel", "svglite", "eulerr", "ggpubr", "cowplot",
          "ggsci", "circlize", "RColorBrewer", "colorspace")
bio_pkgs <- c("biomaRt", "DESeq2", "ComplexHeatmap", "BSgenome.Mmusculus.UCSC.mm10",
              "plotgardener", "Gviz")

# INSTALL PACKAGES IF NECESSARY
#install.packages(pkgs)
#BiocManager::install(bio_pkgs)
#remotes::install_github("EvaYiwenWang/PLSDAbatch")

# LOAD PACKAGES
lapply(pkgs, require, character.only = TRUE)
lapply(bio_pkgs, require, character.only = TRUE)


# CLEANING
rm(pkgs, bio_pkgs)

# ensembl.v102 <- useMart(host = "https://nov2020.archive.ensembl.org",
#                         biomart = "ENSEMBL_MART_ENSEMBL",
#                         dataset = "mmusculus_gene_ensembl")
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
################################################################################
label_kb_mb <- function(x) {
  ifelse(x >= 1000000, paste0(x / 1000000, "Mb"), paste0(x / 1000, "kb"))
}
create_volcanoMA <- function(rdsDir, figDir,
                             dds, output_prefix, note,
                             alpha, foldchange, log2FC_cutoff, geneList,
                             width, height,
                             maxOverlap_volcano = 10,
                             maxOverlap_MA = 10){
  
  # PRE-PROCESSING
  resLFCSort.tb <- readRDS(here(rdsDir, paste0("resLFCSort.tb_", output_prefix, ".rds")))
  resLFCSort.tb <- resLFCSort.tb %>% dplyr::mutate(diffExpressed = case_when(padj < alpha & shrinked_log2FC >= foldchange ~ "UP_strong",
                                                                             padj < alpha & shrinked_log2FC < foldchange & shrinked_log2FC >= 0 ~ "UP_weak",
                                                                             padj < alpha & shrinked_log2FC < 0 & shrinked_log2FC >= -foldchange ~ "DOWN_weak",
                                                                             padj < alpha & shrinked_log2FC < -foldchange ~ "DOWN_strong",
                                                                             TRUE ~ "NO"))
  resLFCSort.tb$diffExpressed <- factor(resLFCSort.tb$diffExpressed, levels = c("NO", "DOWN_weak", "UP_weak", "DOWN_strong", "UP_strong"))
  resLFCSort.tb <- resLFCSort.tb %>% arrange(diffExpressed)
  
  # Apply limit to log2FC
  resLFCSort.tb <- resLFCSort.tb %>% dplyr::mutate(log2FoldChange_MAX = ifelse(abs(log2FoldChange) >= log2FC_cutoff,
                                                                               sign(log2FoldChange)*log2FC_cutoff,
                                                                               log2FoldChange),
                                                   shrinked_log2FC_MAX = ifelse(abs(shrinked_log2FC) >= log2FC_cutoff,
                                                                                sign(shrinked_log2FC)*log2FC_cutoff,
                                                                                shrinked_log2FC),
                                                   rank_metric_up = -log10(padj) * log2FoldChange,
                                                   rank_metric_down = -log10(padj) * (-log2FoldChange))
  
  # Selecting genes to label
  n_top_genes <- 5
  top_genes <- (resLFCSort.tb %>%
    filter(padj < alpha) %>%  # Consider only significant genes
    arrange(desc(rank_metric_up)) %>%
    slice_head(n = n_top_genes))$external_gene_name
  bottom_genes <- (resLFCSort.tb %>%
    filter(padj < alpha) %>%  # Consider only significant genes
    arrange(desc(rank_metric_down)) %>%
    slice_head(n = n_top_genes))$external_gene_name
  
  genesToLabel <- unique(c(top_genes, bottom_genes, geneList))
  genesToLabel.tb <- resLFCSort.tb %>% dplyr::filter(external_gene_name %in% genesToLabel)
  
  
  # Parameters for plotting
  mycolors <- c(strong_red, weak_red, weak_blue, strong_blue, no_grey)
  names(mycolors) <- c("UP_strong", "UP_weak", "DOWN_weak", "DOWN_strong", "NO")
  
  downStrong.num <- table(resLFCSort.tb$diffExpressed)["DOWN_strong"]
  downWeak.num <- table(resLFCSort.tb$diffExpressed)["DOWN_weak"]
  upStrong.num <- table(resLFCSort.tb$diffExpressed)["UP_strong"]
  upWeak.num <- table(resLFCSort.tb$diffExpressed)["UP_weak"]
  
  
  # Volcano
  plot <- ggplot(
    data = resLFCSort.tb,
    aes(
      x = log2FoldChange_MAX,
      y = -log10(padj),
      color = diffExpressed
    )
  ) +
    geom_point(
      size = 1,
      alpha = 1,
      stroke = 0,
      shape = ifelse(abs(resLFCSort.tb$log2FoldChange) >= log2FC_cutoff, 17, 16)) +
    geom_vline(
      xintercept = c(-foldchange, foldchange),
      linetype = "dashed",
      color = "black",
      size = lineThick*mmToLineUnit,
      lineend = "square"
    ) +
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      color = "black",
      size = lineThick*mmToLineUnit,
      lineend = "square"
    ) +
    geom_text_repel(
      data = genesToLabel.tb,
      aes(label = external_gene_name),
      size = fontSizeS*ptToMM,
      hjust = 0,
      color = "black",
      force = 0.1,
      force_pull = 5,
      max.overlaps = maxOverlap_volcano,
      min.segment.length = 0,
      segment.color = 'black',
      segment.size = lineMedium*mmToLineUnit,
      family = "Helvetica"
    ) +
    scale_color_manual(values = mycolors) +
    ggtitle(
      paste0(
        output_prefix, "\n",
        "alpha: ", alpha, ", abs(log2fc_cutoff): ", foldchange, "\n",
        "strong_down: ", downStrong.num, ", ",
        "strong_up: ", upStrong.num, "\n",
        "weak_down: ", downWeak.num, ", ",
        "weak_up: ", upWeak.num
      )
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = fontSizeS,
        family = fontType
      ),
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent"),
    ) +
    xlab("log2(fold change)") +
    ylab("-log10(adjusted p-value)") +
    coord_cartesian(clip = "off") + xlim(c(-log2FC_cutoff, log2FC_cutoff))
  
  fileName <- here(figDir, paste0("volcano_", output_prefix, "_", alpha, "_", foldchange, "_log2FC"))
  
  svglite(paste0(fileName, ".svg"),
          width = width,
          height = height)
  print(plot)
  dev.off()
  
  png(
    paste0(fileName, ".png"),
    width = width,
    height = height,
    res = 600,
    units = "in"
  )
  print(plot)
  dev.off()
  
  # MA
  plot <- ggplot(
    data = resLFCSort.tb,
    aes(
      x = baseMean,
      y = log2FoldChange_MAX,
      color = diffExpressed
    )
  ) +
    geom_point(
      size = 1,
      alpha = 1,
      stroke = 0,
      shape = ifelse(abs(resLFCSort.tb$log2FoldChange) >= log2FC_cutoff, 17, 16)
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "solid",
      color = "black",
      size = lineThick*mmToLineUnit,
      lineend = "square"
    ) +
    # Top labels
    geom_text_repel(
      data = genesToLabel.tb,
      aes(label = external_gene_name),
      size = fontSizeS*ptToMM,
      hjust = 0,
      color = "black",
      force = 0.1,
      force_pull = 5,
      max.overlaps = maxOverlap_MA,
      min.segment.length = 0,
      segment.color = 'black',
      segment.size = lineMedium*mmToLineUnit,
      family = "Helvetica"
    ) +
    scale_color_manual(values = mycolors) +
    ggtitle(
      paste0(
        output_prefix, "\n",
        "alpha: ", alpha, ", abs(log2fc_cutoff): ", foldchange, "\n",
        "strong_down: ", downStrong.num, ", ",
        "strong_up: ", upStrong.num, "\n",
        "weak_down: ", downWeak.num, ", ",
        "weak_up: ", upWeak.num
      )
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = fontSizeS,
        family = fontType
      ),
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent"),
    ) +
    xlab("baseMean") +
    ylab("log2(fold change)") +
    coord_cartesian(clip = "off") + scale_x_log10() + ylim(c(-log2FC_cutoff, log2FC_cutoff))

  fileName <- here(figDir, paste0("ma_", output_prefix, "_", alpha, "_", foldchange, "_log2FC"))

  svglite(paste0(fileName, ".svg"),
          width = width,
          height = height)
  print(plot)
  dev.off()

  png(
    paste0(fileName, ".png"),
    width = width,
    height = height,
    res = 600,
    units = "in"
  )
  print(plot)
  dev.off()
  
  
}

slideBinNormReads <- function(reads.all, binSize, shiftSize, refGenome = BSgenome.Mmusculus.UCSC.mm10){
  # Get genome info and chr size
  chr = as.character(unique(seqnames(reads.all)))
  x = seqinfo(refGenome)
  seqlevels(x) = chr
  
  # Make bin & sliding bin
  gr.bin <- tileGenome(x, tilewidth = binSize, cut.last.tile.in.chrom = TRUE)
  numBins <- floor((max(end(gr.bin))-binSize)/shiftSize)+1
  gr.bin.sliding <- GRanges(
    seqnames = rep(seqlevels(x), numBins),
    ranges = IRanges(
      start = seq(1, by = shiftSize, length.out = numBins),
      end = seq(binSize, by = shiftSize, length.out = numBins)
    ),
    strand = "*"
  )
  
  # Using only non-blind fragment reads
  capture <- reads.all[reads.all$type == "non_blind"]
  capture.rds <- capture[,1]
  # Using normalized non-smoothened readcounts
  capture.rds$nreads <- capture$normReads
  
  
  # Find overlap to bin
  hits <- findOverlaps(gr.bin.sliding, capture.rds)
  hits_dt <- data.table(queryHits = queryHits(hits), subjectHits = subjectHits(hits))
  
  # Convert capture.rds to data.table
  capture_dt <- data.table(nreads = capture.rds$nreads)
  
  # Set keys for efficient join
  setkey(hits_dt, queryHits)
  
  # Perform a join and calculate the mean using data.table's efficient grouping and summarization
  result = hits_dt[, .(average_score = if (.N > 0) mean(capture_dt$nreads[subjectHits]) else 0), by = queryHits]
  
  
  # Find overlaps to bin (assuming findOverlaps() returns hits as a data.frame-like structure)
  hits <- findOverlaps(gr.bin.sliding, capture.rds)
  hits.tb <- tibble(queryHits = queryHits(hits), subjectHits = subjectHits(hits))
  
  temp <- tibble(nreads = capture.rds$nreads, 
                 subjectHits = 1:nrow(as_tibble(capture.rds)))
  
  
  result <- hits.tb %>% left_join(temp, by = "subjectHits") %>%
    group_by(queryHits) %>%
    summarize(avgNormReads = mean(nreads))
  
  gr.bin.sliding.tb <- as_tibble(gr.bin.sliding)
  gr.bin.sliding.tb$queryHits = 1:nrow(gr.bin.sliding.tb)
  
  gr.bin.sliding.tb <- gr.bin.sliding.tb %>% left_join(result, by = "queryHits") %>%
    dplyr::mutate(reads = if_else(is.na(avgNormReads), 0, avgNormReads)) %>%
    dplyr::select(seqnames, start, end, width, strand, reads)
  
  gr.bin.sliding <- makeGRangesFromDataFrame(gr.bin.sliding.tb, keep.extra.columns = TRUE)
  
  return(gr.bin.sliding)
}

createDataTrackNoYlim <- function(file_name, color, track_name) {
  track <- DataTrack(
    import(here(refDir, file_name)),
    strand = "*",
    genome = genome,
    chromosome = chromosome,
    col.histogram = color,
    fill.histogram = color,
    name = track_name,
    col.axis = "black",
    cex = 1, col = "black", 
    fontcolor.title = "black", 
    col.baseline = "black", lwd.baseline = lineThick,
  )
  # Set font size and line width for each data track
  displayPars(track) <- list(fontsize = fontSizeS, lwd = 0.75)
  return(track)
}

createDataTrack <- function(file_name, color, track_name, ymaxEpi) {
  ylim = c(0, ymaxEpi)
  track <- DataTrack(
    import(here(refDir, file_name)),
    strand = "*",
    genome = genome,
    chromosome = chromosome,
    col.histogram = color,
    fill.histogram = color,
    name = track_name,
    col.axis = "black",
    ylim = ylim,
    cex = 1, col = "black", 
    fontcolor.title = "black", 
    col.baseline = "black", lwd.baseline = lineThick,
  )
  # Set font size and line width for each data track
  displayPars(track) <- list(fontsize = fontSizeS, lwd = 0.75)
  return(track)
}

make4Cfigure_type1_4c_withRep <- function(vp.pos, VP, binSize, shiftSize, wSize,
                               range.plus, range.minus, yrange, yrange.min, yrange.max,
                               genome){
  # Import RDS
  reads.ESC.D1 <- readRDS(here(rdsDir, paste0(VP, "_ESC-G1DMSO_Exp1-Rep1.rds")))$reads
  reads.ESC.D2 <- readRDS(here(rdsDir, paste0(VP, "_ESC-G1DMSO_Exp2-Rep1.rds")))$reads
  reads.EPI.D1 <- readRDS(here(rdsDir, paste0(VP, "_Epi-G1DMSO_Exp1-Rep1.rds")))$reads
  reads.EPI.D2 <- readRDS(here(rdsDir, paste0(VP, "_Epi-G1DMSO_Exp2-Rep1.rds")))$reads
  
  vpChr <- as.character(unique(seqnames(reads.ESC.D1)))
  
  # Make sliding bin
  bin.ESC.D1 <- as_tibble(slideBinNormReads(reads.ESC.D1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.ESC.D1 = reads) %>% dplyr::select(start, reads.ESC.D1)
  bin.ESC.D2 <- as_tibble(slideBinNormReads(reads.ESC.D2, binSize, shiftSize)) %>%
    dplyr::mutate(reads.ESC.D2 = reads) %>% dplyr::select(start, reads.ESC.D2)
  bin.EPI.D1 <- as_tibble(slideBinNormReads(reads.EPI.D1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.EPI.D1 = reads) %>% dplyr::select(start, reads.EPI.D1)
  bin.EPI.D2 <- as_tibble(slideBinNormReads(reads.EPI.D2, binSize, shiftSize)) %>%
    dplyr::mutate(reads.EPI.D2 = reads) %>% dplyr::select(start, reads.EPI.D2)
  
  bin.data <- full_join(bin.ESC.D1, bin.ESC.D2, by = "start") %>% 
    full_join(bin.EPI.D1, by = "start") %>%
    full_join(bin.EPI.D2, by = "start")
  
  bin.data <- bin.data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(reads.avg.ESC.D = (reads.ESC.D1 + reads.ESC.D2)/2,
                  reads.avg.EPI.D = (reads.EPI.D1 + reads.EPI.D2)/2,
                  reads.ESC.D.lower = min(reads.ESC.D1, reads.ESC.D2),
                  reads.ESC.D.upper = max(reads.ESC.D1, reads.ESC.D2),
                  reads.EPI.D.lower = min(reads.EPI.D1, reads.EPI.D2),
                  reads.EPI.D.upper = max(reads.EPI.D1, reads.EPI.D2),
                  pos = start + binSize/2)
  
  #### Visualization
  accentColor <- "darkolivegreen1"
  pvalueCutoff <- 0.05
  
  
  
  # Filtering data
  bin.data.ranged <- bin.data %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus)
  
  bin.data.ranged <- bin.data.ranged %>% dplyr::mutate_at(vars(-pos), 
                                                          ~ifelse(pos > vp.pos - binSize &
                                                                    pos < vp.pos + binSize,
                                                                  NA, .))
  # Adding significant box 1, EPI vs ESC
  
  # temp <- bin.data.ranged %>%
  #   dplyr::filter(complete.cases(reads.ESC.D1, reads.ESC.D2, reads.EPI.D1, reads.EPI.D2)) %>%
  #   dplyr::rowwise() %>%
  #   dplyr::mutate(pvalue = (t.test(c(reads.ESC.D1, reads.ESC.D2), c(reads.EPI.D1, reads.EPI.D2),
  #                                  alternative = "two.sided", paired = FALSE))$p.value)
  # 
  # temp2 = temp %>% dplyr::filter(pvalue < pvalueCutoff)
  # rect.df = data.frame(
  #   xmin = temp2$start,
  #   xmax = temp2$start + binSize
  # )
  
  # Panel 1: ESC vs EPI DMSO, binned sliding
  gg1 <- ggplot(bin.data.ranged) +
    # Significant box
    # geom_rect(data = rect.df, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = yrange),
    #           fill = accentColor, alpha = 1) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    # Line plot
    geom_line(aes(x = pos, y = reads.avg.ESC.D), 
              linetype = "solid", color = "black",
              linewidth = lineThick*mmToLineUnit, lineend = "square") +
    geom_line(aes(x = pos, y = reads.avg.EPI.D), 
              linetype = "solid", color = strong_purple,
              linewidth = lineThick*mmToLineUnit, lineend = "square") +
    geom_ribbon(aes(x = pos, ymin = reads.ESC.D.lower, ymax = reads.ESC.D.upper), 
                color = NA, fill = "black", alpha = 0.3) +
    geom_ribbon(aes(x = pos, ymin = reads.EPI.D.lower, ymax = reads.EPI.D.upper), 
                color = NA, fill = strong_purple, alpha = 0.3) +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(0, yrange, by = 500),
                       labels = scales::number_format(accuracy = 1), limits = c(0, yrange)) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = fontSizeS,
        family = fontType
      ),
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    ) +
    labs(x = NULL, y = "normalized 4C signal")  +
    annotate("text", x = vp.pos - range.minus, y = yrange, label = " G1.ESC.DMSO (black), G1.Epi.DMSO (purple)", 
             size = fontSizeM/3, family = fontType, hjust = 0, vjust = 1)
  
  # Panel 2: EPI DMSO - ESC DMSO
  yrange.max2 <- ceiling(max(bin.data.ranged$reads.avg.EPI.D - bin.data.ranged$reads.avg.ESC.D, na.rm = TRUE)/100)*100
  yrange.min2 <- floor(min(bin.data.ranged$reads.avg.EPI.D - bin.data.ranged$reads.avg.ESC.D, na.rm = TRUE)/100)*100
  
  gg2 <- ggplot(bin.data.ranged) +
    # geom_rect(data = rect.df, aes(xmin = xmin, xmax = xmax, ymin = yrange.min2, ymax = yrange.max2),
    #           fill = accentColor, alpha = 1) +
    geom_ribbon(aes(x = pos, ymin = 0, ymax = pmax(0, reads.avg.EPI.D - reads.avg.ESC.D)), 
                fill = strong_red, alpha = 1) +
    geom_ribbon(aes(x = pos, ymin = 0, ymax = pmin(0, reads.avg.EPI.D - reads.avg.ESC.D)), 
                fill = strong_blue,alpha = 1) +
    geom_hline(yintercept = 0,
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(yrange.min2, yrange.max2, by = 200),
                       labels = scales::number_format(accuracy = 1), limits = c(yrange.min2, yrange.max2)) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = fontSizeS,
        family = fontType
      ),
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    ) +
    labs(x = NULL, y =  "normalized 4C signal") +
    annotate("text", x = vp.pos - range.minus, y = yrange.max2, label = " G1.Epi.DMSO - G1.ESC.DMSO", 
             size = fontSizeM/3, family = fontType, hjust = 0, vjust = 1)
  
  
  ######### Plot
  width <- panelSize(3)*mmToInch
  height <- panelSize(1.5)*mmToInch
  
  
  pdf(here(figDir, paste0("visualization_", VP, "_4C_EpivsESC_bin", floor(binSize/1000),
                          "kb_wSize", wSize, 
                          "_range_L", floor(range.minus/1000),
                          "_R", floor(range.plus/1000), ".pdf")), onefile = TRUE, width = width, height = height)
  print(plot_grid(gg1, gg2, align = 'v', axis = 'b', ncol = 1))
  dev.off()
}
make4Cfigure_type1_track <- function(vp.pos, VP, binSize, shiftSize, wSize,
                                  range.plus, range.minus, yrange, yrange.min, yrange.max,
                                  genome, ymaxEpi, genomePanelSize = 2){
  reads.ESC.D1 <- readRDS(here(rdsDir, paste0(VP, "_ESC-G1DMSO_Exp1-Rep1.rds")))$reads
  vpChr <- as.character(unique(seqnames(reads.ESC.D1)))
  
  chromosome <- vpChr
  start_pos <- vp.pos - range.minus
  end_pos <- vp.pos + range.plus
  
  purple1 <- lighten(strong_purple, amount = 0.08*7)
  purple2 <- lighten(strong_purple, amount = 0.08*6)
  purple3 <- lighten(strong_purple, amount = 0.08*5)
  purple4 <- lighten(strong_purple, amount = 0.08*4)
  purple5 <- lighten(strong_purple, amount = 0.08*3)
  purple6<- lighten(strong_purple, amount = 0.08*2)
  purple7 <- lighten(strong_purple, amount = 0.08*1)
  
  # Read data
  axisTrack <- GenomeAxisTrack(labelPos = "below", 
                               col = "black",
                               fontsize = fontSizeS*5/4, fontcolor = "black",
                               lwd = lineThick)
  
  fm <- Gviz:::.getBMFeatureMap()
  fm["symbol"] <- "external_gene_name"
  bm <-BiomartGeneRegionTrack(chromosome=chromosome, genome=genome, 
                              start=start_pos, end = end_pos, 
                              biomart=ensembl.v102, filters=list("biotype" = c("protein_coding", "lncRNA", "miRNA", "lincRNA")),
                              collapseTranscripts="longest",
                              featureMap=fm,
                              name="Gene",
                              cex=1, col = "black", 
                              fontcolor.title = "black",
                              fontcolor = "black", fontfamily = fontType, fontsize = fontSizeS*5/3,
                              fontcolor.group = "black", fontfamily.group = fontType, fontsize.group = fontSizeS*5/3,
                              lwd = lineThick, showId = TRUE)
  
  
  track0 <- createDataTrackNoYlim(
    "GSM2438476_EC-DG-3458-H3K27AC_ASYN_1_bin50bp.bw",
    "black",
    "H3K27ac ESC"
  )
  track1 <- createDataTrack(
    "GSM3314701_H3K27ac-0h_mm10Lifsted.black_bin50bp.bw",
    purple1,
    "H3K27ac EpiLC 0h", ymaxEpi
  )
  track2 <- createDataTrack(
    "GSM3314702_H3K27ac-1h_mm10Lifsted.black_bin50bp.bw",
    purple2,
    "H3K27ac EpiLC 1h", ymaxEpi
  )
  track3 <- createDataTrack(
    "GSM3314703_H3K27ac-6h_mm10Lifsted.black_bin50bp.bw",
    purple3,
    "H3K27ac EpiLC 6h", ymaxEpi
  )
  track4 <- createDataTrack(
    "GSM3314704_H3K27ac-12h_mm10Lifsted.black_bin50bp.bw",
    purple4,
    "H3K27ac EpiLC 12h", ymaxEpi
  )
  track5 <- createDataTrack(
    "GSM3314705_H3K27ac-24h_mm10Lifsted.black_bin50bp.bw",
    purple5,
    "H3K27ac EpiLC 24h", ymaxEpi
  )
  track6 <- createDataTrack(
    "GSM3314706_H3K27ac-36h_mm10Lifsted.black_bin50bp.bw",
    purple6,
    "H3K27ac EpiLC 36h", ymaxEpi
  )
  track7 <- createDataTrack(
    "GSM3314707_H3K27ac-48h_mm10Lifsted.black_bin50bp.bw",
    purple7,
    "H3K27ac EpiLC 48h", ymaxEpi
  )
  track8 <- createDataTrack(
    "GSM3314708_H3K27ac-72h_mm10Lifsted.black_bin50bp.bw",
    strong_purple,
    "H3K27ac EpiLC 72h", ymaxEpi
  )
  
  width <- panelSize(3.28)*mmToInch
  height <- panelSize(2.5)*mmToInch
  
  pdf(here(figDir, paste0("visualization_", VP, "_ChIPseq_bin", floor(binSize/1000),
                          "kb_wSize", wSize, 
                          "_range_L", floor(range.minus/1000),
                          "_R", floor(range.plus/1000), ".pdf")), onefile = TRUE, width = width, height = height)
  plotTracks(c(axisTrack, bm, track0, track1, track2, track3, track4, track5, track6, track7, track8 ),
             from = start_pos, to = end_pos, type = "histogram",
             sizes = c(2, genomePanelSize, 2, 2, 2, 2, 2, 2, 2, 2, 2))
  dev.off()
}

make4Cfigure_type1_4c_noRep <- function(vp.pos, VP, binSize, shiftSize, wSize,
                                  range.plus, range.minus, yrange, yrange.min, yrange.max,
                                  genome){
  # Import RDS
  reads.ESC.D1 <- readRDS(here(rdsDir, paste0(VP, "_ESC-G1DMSO_Exp1-Rep1.rds")))$reads
  reads.EPI.D1 <- readRDS(here(rdsDir, paste0(VP, "_Epi-G1DMSO_Exp1-Rep1.rds")))$reads
  
  vpChr <- as.character(unique(seqnames(reads.ESC.D1)))
  
  # Make sliding bin
  bin.ESC.D1 <- as_tibble(slideBinNormReads(reads.ESC.D1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.ESC.D1 = reads) %>% dplyr::select(start, reads.ESC.D1)
  bin.EPI.D1 <- as_tibble(slideBinNormReads(reads.EPI.D1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.EPI.D1 = reads) %>% dplyr::select(start, reads.EPI.D1)
  
  bin.data <- full_join(bin.ESC.D1, bin.EPI.D1, by = "start")
  
  bin.data <- bin.data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(reads.avg.ESC.D = (reads.ESC.D1),
                  reads.avg.EPI.D = (reads.EPI.D1),
                  pos = start + binSize/2)
  
  #### Visualization
  accentColor <- "darkolivegreen1"
  pvalueCutoff <- 0.05
  
  
  
  # Filtering data
  bin.data.ranged <- bin.data %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus)
  
  bin.data.ranged <- bin.data.ranged %>% dplyr::mutate_at(vars(-pos), 
                                                          ~ifelse(pos > vp.pos - binSize &
                                                                    pos < vp.pos + binSize,
                                                                  NA, .))

  # Panel 1: ESC vs EPI DMSO, binned sliding
  gg1 <- ggplot(bin.data.ranged) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    # Line plot
    geom_line(aes(x = pos, y = reads.avg.ESC.D), 
              linetype = "solid", color = "black",
              linewidth = lineThick*mmToLineUnit, lineend = "square") +
    geom_line(aes(x = pos, y = reads.avg.EPI.D), 
              linetype = "solid", color = strong_purple,
              linewidth = lineThick*mmToLineUnit, lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(0, yrange, by = 500),
                       labels = scales::number_format(accuracy = 1), limits = c(0, yrange)) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = fontSizeS,
        family = fontType
      ),
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    ) +
    labs(x = NULL, y = "normalized 4C signal")  +
    annotate("text", x = vp.pos - range.minus, y = yrange, label = " G1.ESC.DMSO (black), G1.Epi.DMSO (purple)", 
             size = fontSizeM/3, family = fontType, hjust = 0, vjust = 1)
  
  # Panel 2: EPI DMSO - ESC DMSO
  yrange.max2 <- ceiling(max(bin.data.ranged$reads.avg.EPI.D - bin.data.ranged$reads.avg.ESC.D, na.rm = TRUE)/100)*100
  yrange.min2 <- floor(min(bin.data.ranged$reads.avg.EPI.D - bin.data.ranged$reads.avg.ESC.D, na.rm = TRUE)/100)*100
  
  gg2 <- ggplot(bin.data.ranged) +
    geom_ribbon(aes(x = pos, ymin = 0, ymax = pmax(0, reads.avg.EPI.D - reads.avg.ESC.D)), 
                fill = strong_red, alpha = 1) +
    geom_ribbon(aes(x = pos, ymin = 0, ymax = pmin(0, reads.avg.EPI.D - reads.avg.ESC.D)), 
                fill = strong_blue,alpha = 1) +
    geom_hline(yintercept = 0,
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(yrange.min2, yrange.max2, by = 200),
                       labels = scales::number_format(accuracy = 1), limits = c(yrange.min2, yrange.max2)) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = fontSizeS,
        family = fontType
      ),
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    ) +
    labs(x = NULL, y =  "normalized 4C signal") +
    annotate("text", x = vp.pos - range.minus, y = yrange.max2, label = " G1.Epi.DMSO - G1.ESC.DMSO", 
             size = fontSizeM/3, family = fontType, hjust = 0, vjust = 1)
  
  
  ######### Plot
  width <- panelSize(3)*mmToInch
  height <- panelSize(1.5)*mmToInch
  
  
  pdf(here(figDir, paste0("visualization_", VP, "_4C_EpivsESC_bin", floor(binSize/1000),
                          "kb_wSize", wSize, 
                          "_range_L", floor(range.minus/1000),
                          "_R", floor(range.plus/1000), "_noRep.pdf")), onefile = TRUE, width = width, height = height)
  print(plot_grid(gg1, gg2, align = 'v', axis = 'b', ncol = 1))
  dev.off()
}

make4Cfigure_type2_2vs1 <- function(vp.pos, VP, binSize, shiftSize, wSize,
                                        range.plus, range.minus, yrange, yrange.min, yrange.max,
                                        genome, name1, name2, note, sampleColor, width = 3, height = 1.5){
  # Import RDS
  reads.D1 <- readRDS(here(rdsDir, paste0(VP, "_", name1, "_Exp1-Rep1.rds")))$reads
  reads.D2 <- readRDS(here(rdsDir, paste0(VP, "_", name1, "_Exp2-Rep1.rds")))$reads
  reads.T1 <- readRDS(here(rdsDir, paste0(VP, "_", name2, "_Exp1-Rep1.rds")))$reads
  
  vpChr <- as.character(unique(seqnames(reads.D1)))
  
  # Make sliding bin
  bin.D1 <- as_tibble(slideBinNormReads(reads.D1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D1 = reads) %>% dplyr::select(start, reads.D1)
  bin.D2 <- as_tibble(slideBinNormReads(reads.D2, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D2 = reads) %>% dplyr::select(start, reads.D2)
  bin.T2 <- as_tibble(slideBinNormReads(reads.T1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.T1 = reads) %>% dplyr::select(start, reads.T1)
  
  bin.data <- full_join(full_join(bin.D1, bin.D2, by = "start"),
                        bin.T2, by = "start")
  bin.data <- bin.data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(reads.avg.D = (reads.D1 + reads.D2)/2,
                  reads.avg.T = reads.T1,
                  reads.D.lower = min(reads.D1, reads.D2),
                  reads.D.upper = max(reads.D1, reads.D2),
                  pos = start + binSize/2)
  
  #### Visualization
  accentColor <- "darkolivegreen1"
  pvalueCutoff <- 0.05
  
  
  
  # Filtering data
  bin.data.ranged <- bin.data %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus)
  
  bin.data.ranged <- bin.data.ranged %>% dplyr::mutate_at(vars(-pos), 
                                                          ~ifelse(pos > vp.pos - binSize &
                                                                    pos < vp.pos + binSize,
                                                                  NA, .))
  
  # Panel 1: 
  gg1 <- ggplot(bin.data.ranged) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    # Line plot
    geom_line(aes(x = pos, y = reads.avg.D), 
              linetype = "solid", color = "black",
              linewidth = lineThick*mmToLineUnit, lineend = "square") +
    geom_line(aes(x = pos, y = reads.avg.T), 
              linetype = "solid", color = sampleColor,
              linewidth = lineThick*mmToLineUnit, lineend = "square") +
    geom_ribbon(aes(x = pos, ymin = reads.D.lower, ymax = reads.D.upper), 
                color = NA, fill = "black", alpha = 0.3) +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(0, yrange, by = 500),
                       labels = scales::number_format(accuracy = 1), limits = c(0, yrange)) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = fontSizeS,
        family = fontType
      ),
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    ) +
    labs(x = NULL, y = "normalized 4C signal")  +
    annotate("text", x = vp.pos - range.minus, y = yrange, label = note, 
             size = fontSizeM/3, family = fontType, hjust = 0, vjust = 1)
  
  # Panel 2:
  yrange.max2 <- ceiling(max(bin.data.ranged$reads.avg.T - bin.data.ranged$reads.avg.D, na.rm = TRUE)/100)*100
  yrange.min2 <- floor(min(bin.data.ranged$reads.avg.T - bin.data.ranged$reads.avg.D, na.rm = TRUE)/100)*100
  
  gg2 <- ggplot(bin.data.ranged) +
    # geom_rect(data = rect.df, aes(xmin = xmin, xmax = xmax, ymin = yrange.min2, ymax = yrange.max2),
    #           fill = accentColor, alpha = 1) +
    geom_ribbon(aes(x = pos, ymin = 0, ymax = pmax(0, reads.avg.T - reads.avg.D)), 
                fill = strong_red, alpha = 1) +
    geom_ribbon(aes(x = pos, ymin = 0, ymax = pmin(0, reads.avg.T - reads.avg.D)), 
                fill = strong_blue,alpha = 1) +
    geom_hline(yintercept = 0,
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(yrange.min2, yrange.max2, by = 200),
                       labels = scales::number_format(accuracy = 1), limits = c(yrange.min2, yrange.max2)) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = fontSizeS,
        family = fontType
      ),
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    ) +
    labs(x = NULL, y =  "normalized 4C signal") +
    annotate("text", x = vp.pos - range.minus, y = yrange.max2, label = "delta", 
             size = fontSizeM/3, family = fontType, hjust = 0, vjust = 1)
 
  ######### Plot
  width <- panelSize(width)*mmToInch
  height <- panelSize(height)*mmToInch
  
  
  pdf(here(figDir, paste0("visualization_", VP, "_4C_", note,  "_bin", floor(binSize/1000),
                          "kb_wSize", wSize, 
                          "_range_L", floor(range.minus/1000),
                          "_R", floor(range.plus/1000), ".pdf")), onefile = TRUE, width = width, height = height)
  print(plot_grid(gg1, gg2, align = 'v', axis = 'b', ncol = 1))
  dev.off()
}

make4Cfigure_type2_1vs1 <- function(vp.pos, VP, binSize, shiftSize, wSize,
                                    range.plus, range.minus, yrange, yrange.min, yrange.max,
                                    genome, name1, name2, note, sampleColor){
  # Import RDS
  reads.D1 <- readRDS(here(rdsDir, paste0(VP, "_", name1, "_Exp1-Rep1.rds")))$reads
  reads.T1 <- readRDS(here(rdsDir, paste0(VP, "_", name2, "_Exp1-Rep1.rds")))$reads
  
  vpChr <- as.character(unique(seqnames(reads.D1)))
  
  # Make sliding bin
  bin.D1 <- as_tibble(slideBinNormReads(reads.D1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D1 = reads) %>% dplyr::select(start, reads.D1)
  bin.T2 <- as_tibble(slideBinNormReads(reads.T1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.T1 = reads) %>% dplyr::select(start, reads.T1)
  
  bin.data <- full_join(bin.D1, bin.T2, by = "start")
  bin.data <- bin.data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(reads.avg.D = reads.D1,
                  reads.avg.T = reads.T1,
                  pos = start + binSize/2)
  
  #### Visualization
  accentColor <- "darkolivegreen1"
  pvalueCutoff <- 0.05
  
  
  
  # Filtering data
  bin.data.ranged <- bin.data %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus)
  
  bin.data.ranged <- bin.data.ranged %>% dplyr::mutate_at(vars(-pos), 
                                                          ~ifelse(pos > vp.pos - binSize &
                                                                    pos < vp.pos + binSize,
                                                                  NA, .))
  
  # Panel 1: 
  gg1 <- ggplot(bin.data.ranged) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    # Line plot
    geom_line(aes(x = pos, y = reads.avg.D), 
              linetype = "solid", color = "black",
              linewidth = lineThick*mmToLineUnit, lineend = "square") +
    geom_line(aes(x = pos, y = reads.avg.T), 
              linetype = "solid", color = sampleColor,
              linewidth = lineThick*mmToLineUnit, lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(0, yrange, by = 500),
                       labels = scales::number_format(accuracy = 1), limits = c(0, yrange)) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = fontSizeS,
        family = fontType
      ),
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    ) +
    labs(x = NULL, y = "normalized 4C signal")  +
    annotate("text", x = vp.pos - range.minus, y = yrange, label = note, 
             size = fontSizeM/3, family = fontType, hjust = 0, vjust = 1)
  
  # Panel 2:
  yrange.max2 <- ceiling(max(bin.data.ranged$reads.avg.T - bin.data.ranged$reads.avg.D, na.rm = TRUE)/100)*100
  yrange.min2 <- floor(min(bin.data.ranged$reads.avg.T - bin.data.ranged$reads.avg.D, na.rm = TRUE)/100)*100
  
  gg2 <- ggplot(bin.data.ranged) +
    geom_rect(data = rect.df, aes(xmin = xmin, xmax = xmax, ymin = yrange.min2, ymax = yrange.max2),
              fill = accentColor, alpha = 1) +
    geom_ribbon(aes(x = pos, ymin = 0, ymax = pmax(0, reads.avg.T - reads.avg.D)), 
                fill = strong_red, alpha = 1) +
    geom_ribbon(aes(x = pos, ymin = 0, ymax = pmin(0, reads.avg.T - reads.avg.D)), 
                fill = strong_blue,alpha = 1) +
    geom_hline(yintercept = 0,
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(yrange.min2, yrange.max2, by = 200),
                       labels = scales::number_format(accuracy = 1), limits = c(yrange.min2, yrange.max2)) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = fontSizeS,
        family = fontType
      ),
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    ) +
    labs(x = NULL, y =  "normalized 4C signal") +
    annotate("text", x = vp.pos - range.minus, y = yrange.max2, label = "delta", 
             size = fontSizeM/3, family = fontType, hjust = 0, vjust = 1)
  
  ######### Plot
  width <- panelSize(3)*mmToInch
  height <- panelSize(1.5)*mmToInch
  
  
  pdf(here(figDir, paste0("visualization_", VP, "_4C_", note,  "_bin", floor(binSize/1000),
                          "kb_wSize", wSize, 
                          "_range_L", floor(range.minus/1000),
                          "_R", floor(range.plus/1000), "_noRep.pdf")), onefile = TRUE, width = width, height = height)
  print(plot_grid(gg1, gg2, align = 'v', axis = 'b', ncol = 1))
  dev.off()
}

make4Cfigure_pub_tracks_esc <- function(vp.pos, VP, binSize, shiftSize, wSize,
                                        range.plus, range.minus, yrange, yrange.min, yrange.max,
                                        genome, genomePanelSize = 2){
  # Setting up the chr
  reads.ESC.D1 <- readRDS(here(rdsDir, paste0(VP, "_ESC-G1DMSO_Exp1-Rep1.rds")))$reads
  vpChr <- as.character(unique(seqnames(reads.ESC.D1)))
  
  chromosome <- vpChr
  start_pos <- vp.pos - range.minus
  end_pos <- vp.pos + range.plus
  
  # Read data
  axisTrack <- GenomeAxisTrack(labelPos = "below", 
                               col = "black",
                               fontsize = fontSizeS*5/4, fontcolor = "black",
                               lwd = lineThick)
  
  fm <- Gviz:::.getBMFeatureMap()
  fm["symbol"] <- "external_gene_name"
  bm <-BiomartGeneRegionTrack(chromosome=chromosome, genome=genome, 
                              start=start_pos, end = end_pos, 
                              biomart=ensembl.v102, filters=list("biotype" = c("protein_coding", "lncRNA", "miRNA", "lincRNA")),
                              collapseTranscripts="longest",
                              featureMap=fm,
                              name="Gene",
                              cex=1, col = "black", 
                              fontcolor.title = "black",
                              fontcolor = "black", fontfamily = fontType, fontsize = fontSizeS*5/3,
                              fontcolor.group = "black", fontfamily.group = fontType, fontsize.group = fontSizeS*5/3,
                              lwd = lineThick, showId = TRUE)
  
  track0 <- createDataTrackNoYlim(
    "33248_CTCF_07-729_Bruce-4_trim_q20_dedup_black_depthNorm_bin50bp.bw",
    "#8DA0CB",
    "CTCF"
  )
  
  track1 <- createDataTrackNoYlim(
    "33250_RAD21_ab992_Bruce-4_trim_q20_dedup_black_depthNorm_bin50bp.bw",
    "#1F78B4",
    "RAD21"
  )
  
  track2 <- createDataTrackNoYlim(
    "33255_H3K4me3_04-745_Bruce-4_trim_q20_dedup_black_depthNorm_bin50bp.bw",
    "#66C2A5",
    "H3K4me3"
  )
  
  track3 <- createDataTrackNoYlim(
    "GSM2438476_EC-DG-3458-H3K27AC_ASYN_1_bin50bp.bw",
    "#FC8D62",
    "H3K27ac"
  )
  
  width <- panelSize(3)*mmToInch
  height <- panelSize(1.5)*mmToInch
  
  pdf(here(figDir, paste0("pub_visualization_", VP, "_ChIPseq_bin", floor(binSize/1000),
                          "kb_wSize", wSize, 
                          "_range_L", floor(range.minus/1000),
                          "_R", floor(range.plus/1000), ".pdf")), onefile = TRUE, width = width, height = height)
  plotTracks(c(axisTrack, bm, track3, track2, track1, track0),
             from = start_pos, to = end_pos, type = "histogram",
             sizes = c(2, genomePanelSize, 1, 1, 1, 1))
  dev.off()
}

make4Cfigure_pub_tracks_epi <- function(vp.pos, VP, binSize, shiftSize, wSize,
                                        range.plus, range.minus, yrange, yrange.min, yrange.max,
                                        genome, genomePanelSize = 2){
  # Setting up the chr
  reads.ESC.D1 <- readRDS(here(rdsDir, paste0(VP, "_ESC-G1DMSO_Exp1-Rep1.rds")))$reads
  vpChr <- as.character(unique(seqnames(reads.ESC.D1)))
  
  chromosome <- vpChr
  start_pos <- vp.pos - range.minus
  end_pos <- vp.pos + range.plus
  
  # Read data
  axisTrack <- GenomeAxisTrack(labelPos = "below", 
                               col = "black",
                               fontsize = fontSizeS*5/4, fontcolor = "black",
                               lwd = lineThick)
  
  fm <- Gviz:::.getBMFeatureMap()
  fm["symbol"] <- "external_gene_name"
  bm <-BiomartGeneRegionTrack(chromosome=chromosome, genome=genome, 
                              start=start_pos, end = end_pos, 
                              biomart=ensembl.v102, filters=list("biotype" = c("protein_coding", "lncRNA", "miRNA", "lincRNA")),
                              collapseTranscripts="longest",
                              featureMap=fm,
                              name="Gene",
                              cex=1, col = "black", 
                              fontcolor.title = "black",
                              fontcolor = "black", fontfamily = fontType, fontsize = fontSizeS*5/3,
                              fontcolor.group = "black", fontfamily.group = fontType, fontsize.group = fontSizeS*5/3,
                              lwd = lineThick, showId = TRUE)
  


  
  track0 <- createDataTrackNoYlim(
    "GSM3314701_H3K27ac-0h_mm10Lifsted.black_bin50bp.bw",
    "#FC8D62",
    "H3K27ac"
  )
  
  track1 <- createDataTrackNoYlim(
    "GSM3314703_H3K27ac-6h_mm10Lifsted.black_bin50bp.bw",
    darken("#FC8D62", amount = 0.5),
    "H3K4me3"
  )
  
  width <- panelSize(3)*mmToInch
  height <- panelSize(1.1)*mmToInch/6*(4+genomePanelSize)
  
  pdf(here(figDir, paste0("pub_visualization_", VP, "_ChIPseq_bin", floor(binSize/1000),
                          "kb_wSize", wSize, 
                          "_range_L", floor(range.minus/1000),
                          "_R", floor(range.plus/1000), "epi.pdf")), onefile = TRUE, width = width, height = height)
  plotTracks(c(axisTrack, bm, track0, track1),
             from = start_pos, to = end_pos, type = "histogram",
             sizes = c(2, genomePanelSize, 1, 1))
  dev.off()
}

make4Cfigure_pub_type2_2vs1_v4C <- function(vp.chr, vp.pos, VP, binSize, shiftSize, wSize,
                                    range.plus, range.minus, 
                                    yrange, yrange2,
                                    name1, name2, note,
                                    sampleColor, 
                                    v4c1, v4c2,
                                    width = 3, height = 1.5){
  ##############################################################################
  #### 4C
  ##############################################################################
  
  # Import RDS
  reads.D1 <- readRDS(here(rdsDir, paste0(VP, "_", name1, "_Exp1-Rep1.rds")))$reads
  reads.D2 <- readRDS(here(rdsDir, paste0(VP, "_", name1, "_Exp2-Rep1.rds")))$reads
  reads.T1 <- readRDS(here(rdsDir, paste0(VP, "_", name2, "_Exp1-Rep1.rds")))$reads
  
  # Make sliding bin
  bin.D1 <- as_tibble(slideBinNormReads(reads.D1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D1 = reads) %>% dplyr::select(start, reads.D1)
  bin.D2 <- as_tibble(slideBinNormReads(reads.D2, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D2 = reads) %>% dplyr::select(start, reads.D2)
  bin.T2 <- as_tibble(slideBinNormReads(reads.T1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.T1 = reads) %>% dplyr::select(start, reads.T1)
  
  bin.data <- full_join(full_join(bin.D1, bin.D2, by = "start"),
                        bin.T2, by = "start")
  bin.data <- bin.data %>%
    dplyr::mutate(
      pos = start + binSize/2-1,
      reads.avg.D = (reads.D1 + reads.D2)/2,
      reads.avg.T = reads.T1,
      reads.min = pmin(reads.avg.D, reads.avg.T),
      reads.max = pmax(reads.avg.D, reads.avg.T),
      condition = case_when(
        reads.avg.D > reads.avg.T ~ "D_greater",    # A is greater than B
        reads.avg.T > reads.avg.D ~ "T_greater",    # B is greater than A
        TRUE ~ "overlap"                    # Both are equal
      )
      
    )
  
  # Filtering data
  bin.data.ranged <- bin.data %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus) %>%
    dplyr::mutate(reads.D_greater = if_else(condition == "D_greater", reads.avg.D, reads.min),
                  reads.T_greater = if_else(condition == "T_greater", reads.avg.T, reads.min))
  
  # Plotting
  gg1 <- ggplot() +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.D_greater), fill = "black", alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.T_greater), fill = sampleColor, alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.min), fill = lighten(no_grey, amount = 0.5), alpha = 1) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    geom_vline(xintercept = vp.pos,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(0, yrange, by = 500),
                       labels = scales::number_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, yrange)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    )
  
  #############################################################################
  #### v4C
  ##############################################################################
  # v4C
  DMSO.df <- as_tibble(import.bw(here(refDir, v4c1), as = "GRanges")) %>% dplyr::mutate(pos = (end + start-2)/2) %>%
    dplyr::select(seqnames, pos, score)
  dTAG.df <- as_tibble(import.bw(here(refDir, v4c2), as = "GRanges")) %>% dplyr::mutate(pos = (end + start-2)/2) %>%
    dplyr::select(seqnames, pos, score)
  
  combined.df <- full_join(DMSO.df, dTAG.df, by = c("seqnames", "pos"), suffix = c(".DMSO", ".dTAG")) %>%
    mutate(
      reads.min = pmin(score.DMSO, score.dTAG),
      reads.max = pmax(score.DMSO, score.dTAG),
      condition = case_when(
        score.DMSO > score.dTAG ~ "D_greater",    # A is greater than B
        score.dTAG > score.DMSO ~ "T_greater",    # B is greater than A
        TRUE ~ "overlap"                    # Both are equal
      )
    )
  
  bin.data.ranged <- combined.df %>%
    dplyr::filter(seqnames == vp.chr) %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus) %>%
    dplyr::mutate(reads.D_greater = if_else(condition == "D_greater", score.DMSO, reads.min),
                  reads.T_greater = if_else(condition == "T_greater", score.dTAG, reads.min))
  
  # Plotting
  gg2 <- ggplot() +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.D_greater), fill = "black", alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.T_greater), fill = sampleColor, alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.min), fill = lighten(no_grey, amount = 0.5), alpha = 1) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    # scale_y_continuous(breaks = seq(0, yrange, by = 500),
    #                    labels = scales::number_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, yrange2)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    )
  
  
    ######### Plot
  width <- panelSize(width)*mmToInch
  height <- panelSize(height)*mmToInch
  
  pdf(here(figDir, paste0("visualization_", VP, "_4C_", note,  "_bin", floor(binSize/1000),
                          "kb_wSize", wSize, 
                          "_range_L", floor(range.minus/1000),
                          "_R", floor(range.plus/1000), ".pdf")), onefile = TRUE, width = width, height = height)
  print(plot_grid(gg1, gg2, align = 'v', axis = 'b', ncol = 1))
  dev.off()
}

make4Cfigure_pub_type2_epi_2vs1_v4C <- function(vp.chr, vp.pos, VP, binSize, shiftSize, wSize,
                                            range.plus, range.minus, 
                                            yrange, yrange2,
                                            name1, name2, name3, name4, note,
                                            sampleColor, sampleColor2,
                                            v4c1, v4c2,
                                            width = 3, height = 1.5){
  ##############################################################################
  #### 4C
  ##############################################################################
  
  # Import RDS
  reads.D1 <- readRDS(here(rdsDir, paste0(VP, "_", name1, "_Exp1-Rep1.rds")))$reads
  reads.D2 <- readRDS(here(rdsDir, paste0(VP, "_", name1, "_Exp2-Rep1.rds")))$reads
  reads.T1 <- readRDS(here(rdsDir, paste0(VP, "_", name2, "_Exp1-Rep1.rds")))$reads
  
  # Make sliding bin
  bin.D1 <- as_tibble(slideBinNormReads(reads.D1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D1 = reads) %>% dplyr::select(start, reads.D1)
  bin.D2 <- as_tibble(slideBinNormReads(reads.D2, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D2 = reads) %>% dplyr::select(start, reads.D2)
  bin.T2 <- as_tibble(slideBinNormReads(reads.T1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.T1 = reads) %>% dplyr::select(start, reads.T1)
  
  bin.data <- full_join(full_join(bin.D1, bin.D2, by = "start"),
                        bin.T2, by = "start")
  bin.data <- bin.data %>%
    dplyr::mutate(
      pos = start + binSize/2-1,
      reads.avg.D = (reads.D1 + reads.D2)/2,
      reads.avg.T = reads.T1,
      reads.min = pmin(reads.avg.D, reads.avg.T),
      reads.max = pmax(reads.avg.D, reads.avg.T),
      condition = case_when(
        reads.avg.D > reads.avg.T ~ "D_greater",    # A is greater than B
        reads.avg.T > reads.avg.D ~ "T_greater",    # B is greater than A
        TRUE ~ "overlap"                    # Both are equal
      )
      
    )
  
  # Filtering data
  bin.data.ranged <- bin.data %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus) %>%
    dplyr::mutate(reads.D_greater = if_else(condition == "D_greater", reads.avg.D, reads.min),
                  reads.T_greater = if_else(condition == "T_greater", reads.avg.T, reads.min))
  
  # Plotting
  gg1 <- ggplot() +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.D_greater), fill = "black", alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.T_greater), fill = sampleColor, alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.min), fill = lighten(no_grey, amount = 0.5), alpha = 1) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    geom_vline(xintercept = vp.pos,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(0, yrange, by = 500),
                       labels = scales::number_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, yrange)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    )
  
  ##############################################################################
  #### 4C
  ##############################################################################
  
  # Import RDS
  reads.D1 <- readRDS(here(rdsDir, paste0(VP, "_", name3, "_Exp1-Rep1.rds")))$reads
  reads.D2 <- readRDS(here(rdsDir, paste0(VP, "_", name3, "_Exp2-Rep1.rds")))$reads
  reads.T1 <- readRDS(here(rdsDir, paste0(VP, "_", name4, "_Exp1-Rep1.rds")))$reads
  
  # Make sliding bin
  bin.D1 <- as_tibble(slideBinNormReads(reads.D1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D1 = reads) %>% dplyr::select(start, reads.D1)
  bin.D2 <- as_tibble(slideBinNormReads(reads.D2, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D2 = reads) %>% dplyr::select(start, reads.D2)
  bin.T2 <- as_tibble(slideBinNormReads(reads.T1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.T1 = reads) %>% dplyr::select(start, reads.T1)
  
  bin.data <- full_join(full_join(bin.D1, bin.D2, by = "start"),
                        bin.T2, by = "start")
  bin.data <- bin.data %>%
    dplyr::mutate(
      pos = start + binSize/2-1,
      reads.avg.D = (reads.D1 + reads.D2)/2,
      reads.avg.T = reads.T1,
      reads.min = pmin(reads.avg.D, reads.avg.T),
      reads.max = pmax(reads.avg.D, reads.avg.T),
      condition = case_when(
        reads.avg.D > reads.avg.T ~ "D_greater",    # A is greater than B
        reads.avg.T > reads.avg.D ~ "T_greater",    # B is greater than A
        TRUE ~ "overlap"                    # Both are equal
      )
      
    )
  
  # Filtering data
  bin.data.ranged <- bin.data %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus) %>%
    dplyr::mutate(reads.D_greater = if_else(condition == "D_greater", reads.avg.D, reads.min),
                  reads.T_greater = if_else(condition == "T_greater", reads.avg.T, reads.min))
  
  # Plotting
  gg3 <- ggplot() +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.D_greater), fill = "black", alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.T_greater), fill = sampleColor, alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.min), fill = lighten(no_grey, amount = 0.5), alpha = 1) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    geom_vline(xintercept = vp.pos,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(0, yrange, by = 500),
                       labels = scales::number_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, yrange)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    )
  
  ##############################################################################
  #### 4C
  ##############################################################################
  
  # Import RDS
  reads.D1 <- readRDS(here(rdsDir, paste0(VP, "_", name1, "_Exp1-Rep1.rds")))$reads
  reads.D2 <- readRDS(here(rdsDir, paste0(VP, "_", name1, "_Exp2-Rep1.rds")))$reads
  reads.T1 <- readRDS(here(rdsDir, paste0(VP, "_", name3, "_Exp1-Rep1.rds")))$reads
  
  # Make sliding bin
  bin.D1 <- as_tibble(slideBinNormReads(reads.D1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D1 = reads) %>% dplyr::select(start, reads.D1)
  bin.D2 <- as_tibble(slideBinNormReads(reads.D2, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D2 = reads) %>% dplyr::select(start, reads.D2)
  bin.T2 <- as_tibble(slideBinNormReads(reads.T1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.T1 = reads) %>% dplyr::select(start, reads.T1)
  
  bin.data <- full_join(full_join(bin.D1, bin.D2, by = "start"),
                        bin.T2, by = "start")
  bin.data <- bin.data %>%
    dplyr::mutate(
      pos = start + binSize/2-1,
      reads.avg.D = (reads.D1 + reads.D2)/2,
      reads.avg.T = reads.T1,
      reads.min = pmin(reads.avg.D, reads.avg.T),
      reads.max = pmax(reads.avg.D, reads.avg.T),
      condition = case_when(
        reads.avg.D > reads.avg.T ~ "D_greater",    # A is greater than B
        reads.avg.T > reads.avg.D ~ "T_greater",    # B is greater than A
        TRUE ~ "overlap"                    # Both are equal
      )
      
    )
  
  # Filtering data
  bin.data.ranged <- bin.data %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus) %>%
    dplyr::mutate(reads.D_greater = if_else(condition == "D_greater", reads.avg.D, reads.min),
                  reads.T_greater = if_else(condition == "T_greater", reads.avg.T, reads.min))
  
  # Plotting
  gg4 <- ggplot() +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.D_greater), fill = "black", alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.T_greater), fill = sampleColor2, alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.min), fill = lighten(no_grey, amount = 0.5), alpha = 1) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    geom_vline(xintercept = vp.pos,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(0, yrange, by = 500),
                       labels = scales::number_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, yrange)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    )
  
  #############################################################################
  #### v4C
  ##############################################################################
  # v4C
  DMSO.df <- as_tibble(import.bw(here(refDir, v4c1), as = "GRanges")) %>% dplyr::mutate(pos = (end + start-2)/2) %>%
    dplyr::select(seqnames, pos, score)
  dTAG.df <- as_tibble(import.bw(here(refDir, v4c2), as = "GRanges")) %>% dplyr::mutate(pos = (end + start-2)/2) %>%
    dplyr::select(seqnames, pos, score)
  
  combined.df <- full_join(DMSO.df, dTAG.df, by = c("seqnames", "pos"), suffix = c(".DMSO", ".dTAG")) %>%
    mutate(
      reads.min = pmin(score.DMSO, score.dTAG),
      reads.max = pmax(score.DMSO, score.dTAG),
      condition = case_when(
        score.DMSO > score.dTAG ~ "D_greater",    # A is greater than B
        score.dTAG > score.DMSO ~ "T_greater",    # B is greater than A
        TRUE ~ "overlap"                    # Both are equal
      )
    )
  
  bin.data.ranged <- combined.df %>%
    dplyr::filter(seqnames == vp.chr) %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus) %>%
    dplyr::mutate(reads.D_greater = if_else(condition == "D_greater", score.DMSO, reads.min),
                  reads.T_greater = if_else(condition == "T_greater", score.dTAG, reads.min))
  
  # Plotting
  gg2 <- ggplot() +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.D_greater), fill = "black", alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.T_greater), fill = sampleColor, alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.min), fill = lighten(no_grey, amount = 0.5), alpha = 1) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    # scale_y_continuous(breaks = seq(0, yrange, by = 500),
    #                    labels = scales::number_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, yrange2)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    )
  
  
  ######### Plot
  width <- panelSize(width)*mmToInch
  height <- panelSize(height)*mmToInch
  
  pdf(here(figDir, paste0("visualization_", VP, "_4C_", note,  "_bin", floor(binSize/1000),
                          "kb_wSize", wSize, 
                          "_range_L", floor(range.minus/1000),
                          "_R", floor(range.plus/1000), "epi.pdf")), onefile = TRUE, width = width, height = height)
  print(plot_grid(gg4, gg1, gg3, gg2, align = 'v', axis = 'b', ncol = 1))
  dev.off()
}

make4Cfigure_pub_type2_1vs1_v4C <- function(vp.chr, vp.pos, VP, binSize, shiftSize, wSize,
                                            range.plus, range.minus, 
                                            yrange, yrange2,
                                            name1, name2, note,
                                            sampleColor, 
                                            v4c1, v4c2,
                                            width = 3, height = 1.5){
  ##############################################################################
  #### 4C
  ##############################################################################
  
  # Import RDS
  reads.D1 <- readRDS(here(rdsDir, paste0(VP, "_", name1, "_Exp1-Rep1.rds")))$reads
  reads.T1 <- readRDS(here(rdsDir, paste0(VP, "_", name2, "_Exp1-Rep1.rds")))$reads
  
  # Make sliding bin
  bin.D1 <- as_tibble(slideBinNormReads(reads.D1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D1 = reads) %>% dplyr::select(start, reads.D1)
  bin.T1 <- as_tibble(slideBinNormReads(reads.T1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.T1 = reads) %>% dplyr::select(start, reads.T1)
  
  bin.data <- full_join(bin.D1, bin.T1, by = "start")
  bin.data <- bin.data %>%
    dplyr::mutate(
      pos = start + binSize/2-1,
      reads.avg.D = reads.D1,
      reads.avg.T = reads.T1,
      reads.min = pmin(reads.avg.D, reads.avg.T),
      reads.max = pmax(reads.avg.D, reads.avg.T),
      condition = case_when(
        reads.avg.D > reads.avg.T ~ "D_greater",    # A is greater than B
        reads.avg.T > reads.avg.D ~ "T_greater",    # B is greater than A
        TRUE ~ "overlap"                    # Both are equal
      )
      
    )
  
  # Filtering data
  bin.data.ranged <- bin.data %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus) %>%
    dplyr::mutate(reads.D_greater = if_else(condition == "D_greater", reads.avg.D, reads.min),
                  reads.T_greater = if_else(condition == "T_greater", reads.avg.T, reads.min))
  
  # Plotting
  gg1 <- ggplot() +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.D_greater), fill = "black", alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.T_greater), fill = sampleColor, alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.min), fill = lighten(no_grey, amount = 0.5), alpha = 1) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    geom_vline(xintercept = vp.pos,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(0, yrange, by = 500),
                       labels = scales::number_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, yrange)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    )
  
  #############################################################################
  #### v4C
  ##############################################################################
  # v4C
  DMSO.df <- as_tibble(import.bw(here(refDir, v4c1), as = "GRanges")) %>% dplyr::mutate(pos = (end + start-2)/2) %>%
    dplyr::select(seqnames, pos, score)
  dTAG.df <- as_tibble(import.bw(here(refDir, v4c2), as = "GRanges")) %>% dplyr::mutate(pos = (end + start-2)/2) %>%
    dplyr::select(seqnames, pos, score)
  
  combined.df <- full_join(DMSO.df, dTAG.df, by = c("seqnames", "pos"), suffix = c(".DMSO", ".dTAG")) %>%
    mutate(
      reads.min = pmin(score.DMSO, score.dTAG),
      reads.max = pmax(score.DMSO, score.dTAG),
      condition = case_when(
        score.DMSO > score.dTAG ~ "D_greater",    # A is greater than B
        score.dTAG > score.DMSO ~ "T_greater",    # B is greater than A
        TRUE ~ "overlap"                    # Both are equal
      )
    )
  
  bin.data.ranged <- combined.df %>%
    dplyr::filter(seqnames == vp.chr) %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus) %>%
    dplyr::mutate(reads.D_greater = if_else(condition == "D_greater", score.DMSO, reads.min),
                  reads.T_greater = if_else(condition == "T_greater", score.dTAG, reads.min))
  
  # Plotting
  gg2 <- ggplot() +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.D_greater), fill = "black", alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.T_greater), fill = sampleColor, alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.min), fill = lighten(no_grey, amount = 0.5), alpha = 1) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    # scale_y_continuous(breaks = seq(0, yrange, by = 500),
    #                    labels = scales::number_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, yrange2)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    )
  
  
  ######### Plot
  width <- panelSize(width)*mmToInch
  height <- panelSize(height)*mmToInch
  
  pdf(here(figDir, paste0("visualization_", VP, "_4C_", note,  "_bin", floor(binSize/1000),
                          "kb_wSize", wSize, 
                          "_range_L", floor(range.minus/1000),
                          "_R", floor(range.plus/1000), ".pdf")), onefile = TRUE, width = width, height = height)
  print(plot_grid(gg1, gg2, align = 'v', axis = 'b', ncol = 1))
  dev.off()
}

make4Cfigure_pub_type2_epi_1vs1_v4C <- function(vp.chr, vp.pos, VP, binSize, shiftSize, wSize,
                                                range.plus, range.minus, 
                                                yrange, yrange2,
                                                name1, name2, name3, name4, note,
                                                sampleColor, sampleColor2,
                                                v4c1, v4c2,
                                                width = 3, height = 1.5){
  ##############################################################################
  #### 4C
  ##############################################################################
  
  # Import RDS
  reads.D1 <- readRDS(here(rdsDir, paste0(VP, "_", name1, "_Exp1-Rep1.rds")))$reads
  reads.T1 <- readRDS(here(rdsDir, paste0(VP, "_", name2, "_Exp1-Rep1.rds")))$reads
  
  # Make sliding bin
  bin.D1 <- as_tibble(slideBinNormReads(reads.D1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D1 = reads) %>% dplyr::select(start, reads.D1)
  bin.T2 <- as_tibble(slideBinNormReads(reads.T1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.T1 = reads) %>% dplyr::select(start, reads.T1)
  
  bin.data <- full_join(bin.D1, bin.T2, by = "start")
  bin.data <- bin.data %>%
    dplyr::mutate(
      pos = start + binSize/2-1,
      reads.avg.D = (reads.D1),
      reads.avg.T = reads.T1,
      reads.min = pmin(reads.avg.D, reads.avg.T),
      reads.max = pmax(reads.avg.D, reads.avg.T),
      condition = case_when(
        reads.avg.D > reads.avg.T ~ "D_greater",    # A is greater than B
        reads.avg.T > reads.avg.D ~ "T_greater",    # B is greater than A
        TRUE ~ "overlap"                    # Both are equal
      )
      
    )
  
  # Filtering data
  bin.data.ranged <- bin.data %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus) %>%
    dplyr::mutate(reads.D_greater = if_else(condition == "D_greater", reads.avg.D, reads.min),
                  reads.T_greater = if_else(condition == "T_greater", reads.avg.T, reads.min))
  
  # Plotting
  gg1 <- ggplot() +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.D_greater), fill = "black", alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.T_greater), fill = sampleColor, alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.min), fill = lighten(no_grey, amount = 0.5), alpha = 1) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    geom_vline(xintercept = vp.pos,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(0, yrange, by = 500),
                       labels = scales::number_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, yrange)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    )
  
  ##############################################################################
  #### 4C
  ##############################################################################
  
  # Import RDS
  reads.D1 <- readRDS(here(rdsDir, paste0(VP, "_", name3, "_Exp1-Rep1.rds")))$reads
  reads.T1 <- readRDS(here(rdsDir, paste0(VP, "_", name4, "_Exp1-Rep1.rds")))$reads
  
  # Make sliding bin
  bin.D1 <- as_tibble(slideBinNormReads(reads.D1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D1 = reads) %>% dplyr::select(start, reads.D1)
  bin.T2 <- as_tibble(slideBinNormReads(reads.T1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.T1 = reads) %>% dplyr::select(start, reads.T1)
  
  bin.data <- full_join(bin.D1, bin.T2, by = "start")
  bin.data <- bin.data %>%
    dplyr::mutate(
      pos = start + binSize/2-1,
      reads.avg.D = (reads.D1),
      reads.avg.T = reads.T1,
      reads.min = pmin(reads.avg.D, reads.avg.T),
      reads.max = pmax(reads.avg.D, reads.avg.T),
      condition = case_when(
        reads.avg.D > reads.avg.T ~ "D_greater",    # A is greater than B
        reads.avg.T > reads.avg.D ~ "T_greater",    # B is greater than A
        TRUE ~ "overlap"                    # Both are equal
      )
      
    )
  
  # Filtering data
  bin.data.ranged <- bin.data %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus) %>%
    dplyr::mutate(reads.D_greater = if_else(condition == "D_greater", reads.avg.D, reads.min),
                  reads.T_greater = if_else(condition == "T_greater", reads.avg.T, reads.min))
  
  # Plotting
  gg3 <- ggplot() +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.D_greater), fill = "black", alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.T_greater), fill = sampleColor, alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.min), fill = lighten(no_grey, amount = 0.5), alpha = 1) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    geom_vline(xintercept = vp.pos,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(0, yrange, by = 500),
                       labels = scales::number_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, yrange)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    )
  
  ##############################################################################
  #### 4C
  ##############################################################################
  
  # Import RDS
  reads.D1 <- readRDS(here(rdsDir, paste0(VP, "_", name1, "_Exp1-Rep1.rds")))$reads
  reads.T1 <- readRDS(here(rdsDir, paste0(VP, "_", name3, "_Exp1-Rep1.rds")))$reads
  
  # Make sliding bin
  bin.D1 <- as_tibble(slideBinNormReads(reads.D1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D1 = reads) %>% dplyr::select(start, reads.D1)
  bin.T2 <- as_tibble(slideBinNormReads(reads.T1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.T1 = reads) %>% dplyr::select(start, reads.T1)
  
  bin.data <- full_join(bin.D1, bin.T2, by = "start")
  bin.data <- bin.data %>%
    dplyr::mutate(
      pos = start + binSize/2-1,
      reads.avg.D = (reads.D1),
      reads.avg.T = reads.T1,
      reads.min = pmin(reads.avg.D, reads.avg.T),
      reads.max = pmax(reads.avg.D, reads.avg.T),
      condition = case_when(
        reads.avg.D > reads.avg.T ~ "D_greater",    # A is greater than B
        reads.avg.T > reads.avg.D ~ "T_greater",    # B is greater than A
        TRUE ~ "overlap"                    # Both are equal
      )
      
    )
  
  # Filtering data
  bin.data.ranged <- bin.data %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus) %>%
    dplyr::mutate(reads.D_greater = if_else(condition == "D_greater", reads.avg.D, reads.min),
                  reads.T_greater = if_else(condition == "T_greater", reads.avg.T, reads.min))
  
  # Plotting
  gg4 <- ggplot() +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.D_greater), fill = "black", alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.T_greater), fill = sampleColor2, alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.min), fill = lighten(no_grey, amount = 0.5), alpha = 1) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    geom_vline(xintercept = vp.pos,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(0, yrange, by = 500),
                       labels = scales::number_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, yrange)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    )
  
  #############################################################################
  #### v4C
  ##############################################################################
  # v4C
  DMSO.df <- as_tibble(import.bw(here(refDir, v4c1), as = "GRanges")) %>% dplyr::mutate(pos = (end + start-2)/2) %>%
    dplyr::select(seqnames, pos, score)
  dTAG.df <- as_tibble(import.bw(here(refDir, v4c2), as = "GRanges")) %>% dplyr::mutate(pos = (end + start-2)/2) %>%
    dplyr::select(seqnames, pos, score)
  
  combined.df <- full_join(DMSO.df, dTAG.df, by = c("seqnames", "pos"), suffix = c(".DMSO", ".dTAG")) %>%
    mutate(
      reads.min = pmin(score.DMSO, score.dTAG),
      reads.max = pmax(score.DMSO, score.dTAG),
      condition = case_when(
        score.DMSO > score.dTAG ~ "D_greater",    # A is greater than B
        score.dTAG > score.DMSO ~ "T_greater",    # B is greater than A
        TRUE ~ "overlap"                    # Both are equal
      )
    )
  
  bin.data.ranged <- combined.df %>%
    dplyr::filter(seqnames == vp.chr) %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus) %>%
    dplyr::mutate(reads.D_greater = if_else(condition == "D_greater", score.DMSO, reads.min),
                  reads.T_greater = if_else(condition == "T_greater", score.dTAG, reads.min))
  
  # Plotting
  gg2 <- ggplot() +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.D_greater), fill = "black", alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.T_greater), fill = sampleColor, alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.min), fill = lighten(no_grey, amount = 0.5), alpha = 1) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    # scale_y_continuous(breaks = seq(0, yrange, by = 500),
    #                    labels = scales::number_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, yrange2)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    )
  
  
  ######### Plot
  width <- panelSize(width)*mmToInch
  height <- panelSize(height)*mmToInch
  
  pdf(here(figDir, paste0("visualization_", VP, "_4C_", note,  "_bin", floor(binSize/1000),
                          "kb_wSize", wSize, 
                          "_range_L", floor(range.minus/1000),
                          "_R", floor(range.plus/1000), "epi.pdf")), onefile = TRUE, width = width, height = height)
  print(plot_grid(gg4, gg1, gg3, gg2, align = 'v', axis = 'b', ncol = 1))
  dev.off()
}


make4Cfigure_pub_epiVSesc_2vs1_v4C <- function(vp.chr, vp.pos, VP, binSize, shiftSize, wSize,
                                            range.plus, range.minus, 
                                            yrange, yrange2,
                                            name1, name2, note,
                                            sampleColor, 
                                            v4c1, v4c2,
                                            width = 3, height = 1.5){
  ##############################################################################
  #### 4C
  ##############################################################################
  
  # Import RDS
  reads.D1 <- readRDS(here(rdsDir, paste0(VP, "_", name1, "_Exp1-Rep1.rds")))$reads
  reads.D2 <- readRDS(here(rdsDir, paste0(VP, "_", name1, "_Exp2-Rep1.rds")))$reads
  reads.T1 <- readRDS(here(rdsDir, paste0(VP, "_", name2, "_Exp1-Rep1.rds")))$reads
  
  # Make sliding bin
  bin.D1 <- as_tibble(slideBinNormReads(reads.D1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D1 = reads) %>% dplyr::select(start, reads.D1)
  bin.D2 <- as_tibble(slideBinNormReads(reads.D2, binSize, shiftSize)) %>%
    dplyr::mutate(reads.D2 = reads) %>% dplyr::select(start, reads.D2)
  bin.T2 <- as_tibble(slideBinNormReads(reads.T1, binSize, shiftSize)) %>%
    dplyr::mutate(reads.T1 = reads) %>% dplyr::select(start, reads.T1)
  
  bin.data <- full_join(full_join(bin.D1, bin.D2, by = "start"),
                        bin.T2, by = "start")
  bin.data <- bin.data %>%
    dplyr::mutate(
      pos = start + binSize/2-1,
      reads.avg.D = (reads.D1 + reads.D2)/2,
      reads.avg.T = reads.T1,
      reads.min = pmin(reads.avg.D, reads.avg.T),
      reads.max = pmax(reads.avg.D, reads.avg.T),
      condition = case_when(
        reads.avg.D > reads.avg.T ~ "D_greater",    # A is greater than B
        reads.avg.T > reads.avg.D ~ "T_greater",    # B is greater than A
        TRUE ~ "overlap"                    # Both are equal
      )
      
    )
  
  # Filtering data
  bin.data.ranged <- bin.data %>%
    dplyr::filter(pos > vp.pos - range.minus) %>%
    dplyr::filter(pos < vp.pos + range.plus) %>%
    dplyr::mutate(reads.D_greater = if_else(condition == "D_greater", reads.avg.D, reads.min),
                  reads.T_greater = if_else(condition == "T_greater", reads.avg.T, reads.min))
  
  # Plotting
  gg1 <- ggplot() +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.D_greater), fill = "black", alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.T_greater), fill = sampleColor, alpha = 1) +
    geom_ribbon(data = bin.data.ranged,
                aes(x = pos, ymin = 0, ymax = reads.min), fill = lighten(no_grey, amount = 0.5), alpha = 1) +
    # Zero line
    geom_hline(yintercept = 0,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    geom_vline(xintercept = vp.pos,
               color = "black",
               size = lineThick*mmToLineUnit,
               lineend = "square") +
    scale_x_continuous(expand = c(0, 0), labels = label_kb_mb,
                       limits = c(vp.pos - range.minus, vp.pos + range.plus)) +
    scale_y_continuous(breaks = seq(0, yrange, by = 500),
                       labels = scales::number_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, yrange)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.text = element_text(
        size = fontSizeS,
        family = fontType,
        color = "#000000"
      ),
      axis.line = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "#000000",
        size = lineThick*mmToLineUnit,
        lineend = "square"
      ),
      panel.background = element_rect(fill = "transparent")
    )
  
  ######### Plot
  width <- panelSize(width)*mmToInch
  height <- panelSize(height)*mmToInch
  
  pdf(here(figDir, paste0("visualization_", VP, "_4C_", note,  "_bin", floor(binSize/1000),
                          "kb_wSize", wSize, 
                          "_range_L", floor(range.minus/1000),
                          "_R", floor(range.plus/1000), "esc_vs_epi.pdf")), onefile = TRUE, width = width, height = height)
  plot(gg1)
  dev.off()
}
