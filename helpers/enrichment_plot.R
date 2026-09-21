# ===========================================================================
#  helpers/enrichment_plot.R
#
#  GO / KEGG enrichment of a gene-symbol list, drawn as a horizontal
#  rounded-bar plot (top 5 terms per ontology, plus the top 5 KEGG pathways).
#
#  This helper is used by steps 19 and 20.  It is factored out because the
#  enrichment step is run twice - once on the proteins associated with each
#  CpG, and once on the proteins associated with breast cancer - with
#  identical code.
#
#      source("config/config.R")
#      source("helpers/enrichment_plot.R")
# ===========================================================================

#' Run GO (BP/CC/MF) and KEGG enrichment on a vector of gene symbols.
#'
#' @param genes  Character vector of HGNC gene symbols.
#' @return A data frame with the selected terms, ready for plotting, or NULL
#'         when no gene maps to an Entrez ID.
enrich_go_kegg <- function(genes) {

  genes <- unique(as.vector(genes))

  entrez <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
  entrez <- as.character(entrez)
  geneid <- entrez[entrez != "NA"]

  if (length(geneid) == 0) return(NULL)

  # GO over all three ontologies; no filtering here, the top terms are
  # selected afterwards so that the same code path handles every input size.
  go <- enrichGO(gene = geneid, OrgDb = org.Hs.eg.db,
                 pvalueCutoff = 1, qvalueCutoff = 1, ont = "all", readable = TRUE)
  GO <- as.data.frame(go)

  kk <- enrichKEGG(gene = geneid, organism = "hsa",
                   pvalueCutoff = 1, qvalueCutoff = 1)
  KEGG <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  KEGG <- as.data.frame(KEGG)

  use_pathway <- group_by(GO, ONTOLOGY) %>%
    top_n(5, wt = -p.adjust) %>%
    group_by(p.adjust) %>%
    top_n(1, wt = Count) %>%
    rbind(
      top_n(KEGG, 5, -p.adjust) %>%
        group_by(p.adjust) %>%
        top_n(1, wt = Count) %>%
        mutate(ONTOLOGY = "KEGG")
    ) %>%
    ungroup() %>%
    mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")))) %>%
    dplyr::arrange(ONTOLOGY, p.adjust) %>%
    mutate(Description = factor(Description, levels = Description)) %>%
    tibble::rowid_to_column("index")

  use_pathway
}


#' Draw the enrichment plot and write it to a PDF.
#'
#' @param use_pathway  Output of enrich_go_kegg().
#' @param out_file     Destination PDF path.
#' @param pal          Four colours for BP / CC / MF / KEGG.
#' @param width,height Figure dimensions in inches.
plot_enrichment <- function(use_pathway, out_file,
                            pal = c("#7bc4e2", "#acd372", "#fbb05b", "#ed6ca4"),
                            width = 12, height = 11) {

  if (is.null(use_pathway) || nrow(use_pathway) == 0) {
    message("  no enrichment terms to plot")
    return(invisible(NULL))
  }

  width_bar  <- 0.5                              # width of the left-hand label block
  xaxis_max  <- max(-log10(use_pathway$p.adjust)) + 1

  rect_data <- group_by(use_pathway, ONTOLOGY) %>%
    reframe(n = n()) %>%
    ungroup() %>%
    mutate(xmin = -3 * width_bar, xmax = -2 * width_bar,
           ymax = cumsum(n),
           ymin = lag(ymax, default = 0) + 0.6,
           ymax = ymax + 0.4)

  p <- use_pathway %>%
    ggplot(aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
    geom_round_col(aes(y = Description), width = 0.6, alpha = 0.8) +
    geom_text(aes(x = 0.05, label = Description), hjust = 0, size = 5) +
    geom_text(aes(x = 0.1, label = geneID, colour = ONTOLOGY),
              hjust = 0, vjust = 2.6, size = 3.5,
              fontface = "italic", show.legend = FALSE) +
    geom_point(aes(x = -width_bar, size = Count), shape = 21) +
    geom_text(aes(x = -width_bar, label = Count)) +
    scale_size_continuous(name = "Count", range = c(5, 16)) +
    geom_round_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                        fill = ONTOLOGY),
                    data = rect_data, radius = unit(2, "mm"),
                    inherit.aes = FALSE) +
    geom_text(aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
                  label = ONTOLOGY),
              data = rect_data, inherit.aes = FALSE) +
    geom_segment(aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
                 linewidth = 1.5, inherit.aes = FALSE) +
    labs(y = NULL) +
    scale_fill_manual(name = "Category", values = pal) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(breaks = seq(0, xaxis_max, 2), expand = expansion(c(0, 0))) +
    theme_prism() +
    theme(axis.text.y = element_blank(),
          axis.line   = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_text())

  ggsave(out_file, p, width = width, height = height)
  invisible(p)
}
