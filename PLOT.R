

library("heatmaply")
plot_folder <- "../PLOT/"
heatmaply(df, 
          main = "title",
          dendrogram = "none",
          cellnote = df,
          cellnote_size = 12,
          cellnote_textposition="middle center",
          fontsize_row = 10,
          fontsize_col = 10,
          height=520,width=900,
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "springgreen4",  high = "deeppink3"),
          file = paste0(plot_folder,"plot_name.png")
)
