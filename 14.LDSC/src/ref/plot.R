options(repr.plot.width = 80, repr.plot.height = 40)
plot <- ggplot(ldsc_table_plot, aes(x = index, y = phenotype, fill = enrichment_score, label = stars)) +
    geom_tile(color = "black", linewidth = 2) +
    geom_text(size = 20) +
    scale_fill_distiller(palette = "RdBu") + 
    coord_equal() +
    theme_minimal() +
    theme(
        text = element_text(family = "mono", size = 50),
        axis.text.x = element_text(hjust = 1, angle = 70),
        # axis.text.y = element_text(size = 40),
        # axis.title.x = element_text(size = 50),
        # axis.title.y = element_text(size = 50),
        legend.key.height = unit(2, "cm"),
        # legend.text = element_text(size = 40),
        # legend.title = element_text(size = 50)
    ) +
    scale_y_discrete(label = sapply(strsplit(levels(ldsc_table_plot$phenotype), split = "[.]"), tail, 1))
plot
ggsave(
    filename = sprintf("%s/combined-ldsc.pdf", figure_dir), 
    plot = plot, 
    dpi = 300, 
    width = 80, 
    height = 40, 
    device = "pdf",
    limitsize = FALSE
)
