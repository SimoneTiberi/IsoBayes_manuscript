memory_plot = function(df, title){
  ggplot(df, aes(x = Model, y = RAM, label = RAM, fill = Model)) +
    geom_bar(stat = "identity") +
    scale_fill_manual("Model", values = PALETTE_MODELS) +
    geom_text(size = 4.25, colour = "black",  vjust=-0.25) +
    labs(x = NULL,
         y = "Memory (MB)",
         title = title) +
    theme_classic() +
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 12,
                                     angle = 45,
                                     margin = margin(t = 40, r = 0,
                                                     b = -40, l = 0,
                                                     unit = "pt")),
          axis.text.y = element_text(size = 10),
          legend.position = "none")
}
