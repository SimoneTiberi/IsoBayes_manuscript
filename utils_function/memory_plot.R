memory_plot = function(df, title){
  ggplot(df, aes(x = Model, y = RAM, label = RAM, fill = Model)) +
    geom_bar(stat = "identity") +
    scale_fill_manual("Model", values = PALETTE_MODELS) +
    geom_text(size = 4.25, colour = "black",  vjust=-0.25) +
    labs(y = "RAM (MB)",
         x = "Tool",
         title = title) +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold", angle = 10, margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
          axis.text.y = element_text(size = 10, face = "bold"),
          legend.position = "none")
}


  labs(y = "RAM (MB)",
       x = "Tool",
       title = glue("RAM - {DATA}")) +