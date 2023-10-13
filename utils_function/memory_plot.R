memory_plot = function(df, title){
  ggplot(df, aes(x = Model, y = RAM, label = RAM, fill = Model)) +
    geom_bar(stat = "identity") +
    scale_fill_manual("Model", values = PALETTE_MODELS) +
    geom_text(size = 4.25, colour = "black",  vjust=-0.25) +
    labs(x = "",
         y = "Memory (MB)",
         title = title) +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title.y = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 12, face = "bold",
                                     angle = 45, margin = margin(t = 40, r = 0,
                                                                 b = 0, l = 0,
                                                                 unit = "pt")),
          axis.text.y = element_text(size = 10, face = "bold"),
          legend.position = "none") #+
    #geom_errorbar(aes(ymin=RAM-RAM_sd, ymax=RAM+RAM_sd), width=.1, size = 0.5,
     #             position=position_dodge(.9))
}
