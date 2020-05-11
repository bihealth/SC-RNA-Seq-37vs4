modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, cols=colors,... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}
clean_up <- function(sobj, nnn=10, nsd=3) {
  
  df <- FetchData(sobj,c('UMAP_1','UMAP_2','ident'))
  df$ident <- gsub(' [0-9]$','',df$ident)
  
  good.cells <- df %>%
    tibble::rownames_to_column('cell') %>%
    group_by(ident) %>%
    dplyr::mutate(x=mean(UMAP_1),
                  y=mean(UMAP_2)) %>%
    dplyr::mutate(spread=sqrt(mean((UMAP_1-x)**2+
                                     (UMAP_2-y)**2))) %>%
    dplyr::mutate(center.dist=sqrt((UMAP_1-x)**2+(UMAP_2-y)**2)/spread) %>%
    dplyr::filter(center.dist < nsd) %>%
    dplyr::pull(cell)
  
  dist.mat <- as.matrix(dist(df[c('UMAP_1','UMAP_2')]))
  
  keep <- c()
  for (cell in good.cells) {
    neighbors <- names(sort(dist.mat[cell,])[2:nnn])
    if (names(which.max(table(df[neighbors,'ident'])))==df[cell,'ident']) {
      keep <- c(keep,cell)
    }
  }
  #  print(paste('keeping',length(keep),'of',length(sobj@cell.names),'cells'))
  subset(sobj,cells=keep)
}



