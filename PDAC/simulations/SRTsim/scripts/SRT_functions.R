library('patchwork')
library('ggplot2')

base_plot <- function(data, color){
  plot <- ggplot(data=data, aes(x=x, y=y, color=color))+
          geom_point(size=0.3)+
          scale_y_reverse() + 
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            aspect.ratio = 1)

  return(plot)
}

compare_ref_sim <- function(gene, simSRT_obj){
  # Gene Expression
  if(!grepl('ucell', tolower(gene))){
    ref_count_matrix <- cbind.data.frame(simSRT_obj@refcolData[,c("x","y")], sel_gene = simSRT_obj@refCounts[gene,])
    sim_count_matrix <- cbind.data.frame(simSRT_obj@refcolData[,c("x","y")], sel_gene = simSRT_obj@simCounts[gene,])
    limits_color <- c(NA, NA)
  } 
  
  # Uscore
  if(grepl('ucell', tolower(gene))){

    library('UCell')
    gene.sets <- list() 
    gene.sets$classical <- c("TFF1","TFF2","TFF3","CEACAM6", "LGALS4", "ST6GALNAC1", "PLA2G10","TSPAN8","LYZ","MYO1A", "VSIG2", "CLRN3", "CDH17", "AGR3", "AGR2", "BTNL8", "ANXA10", "FAM3D", "CTSE", "REG4")
    gene.sets$basal <- c("SERPINB3", "SPRR3","SERPINB4", "VGLL1","DHRS9", "SPRR1B", "KRT17", "KRT15", "TNS4", "SCEL", "KRT6A", "KRT7", "CST6", "LY6D", "FAM83A", "AREG", "FGFBP1", "GPR87", "LEMD1","S100A2","SLC2A1")

    ref_Ucell <- ScoreSignatures_UCell(simSRT_obj@refCounts, features = gene.sets)
    sim_Ucell <- ScoreSignatures_UCell(simSRT_obj@simCounts, features = gene.sets)
    limits_color <- c(0,1)
    
    if(grepl('classical', tolower(gene))){
      ref_count_matrix <- cbind.data.frame(simSRT_obj@refcolData[,c("x","y")], ref_Ucell[,'classical_UCell'])
      sim_count_matrix <- cbind.data.frame(simSRT_obj@simcolData[,c("x","y")], sim_Ucell[,'classical_UCell'])
    }
    if(grepl('basal', tolower(gene))){
      ref_count_matrix <- cbind.data.frame(simSRT_obj@refcolData[,c("x","y")], sel_gene = ref_Ucell[,'basal_UCell'])
      sim_count_matrix <- cbind.data.frame(simSRT_obj@simcolData[,c("x","y")], sel_gene = sim_Ucell[,'basal_UCell'])
    }
  }
  
  colnames(ref_count_matrix) <- colnames(sim_count_matrix) <- c("x","y","sel_gene")
  
  # Vizualisation
  ref <- base_plot(ref_count_matrix, ref_count_matrix[['sel_gene']])+ 
    ggtitle(paste0("Reference - ", gene))+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_color_gradientn(name=gene, limits = limits_color, colors = c('blue',  'yellow', 'red'))+
    coord_fixed()
  
  sim <- base_plot(sim_count_matrix, sim_count_matrix[['sel_gene']])+ 
    ggtitle(paste0("Simulation - ", gene))+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_color_gradientn(name=gene, limits = limits_color, colors = c('blue',  'yellow', 'red'))+
    coord_fixed()

  wrap_plots(ref, sim)
}

redesign_visualisation <- function(data){
  color_palette <- c('basal' = 'red', 'classical' = 'yellow', 'classical and basal' = 'orange', 'other' = 'blue')

  ref <- base_plot(data, data$label)+scale_color_manual(values = color_palette)
  target <- base_plot(data, data$target_label)+scale_color_manual(values = color_palette)
  wrap_plots(ref, target)
}
