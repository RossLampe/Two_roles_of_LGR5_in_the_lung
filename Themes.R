color_Mesc = "#198BFD"   #neutral blue
color_Epi = "#EC6464"    #neutral red

colors_DotPlot <- c("blue","red")
colors_FeaturePlot <- c("yellow", "blue")

theme_DotPlot <- function(p_DotPlot){
  #rt <- theme_get() #retrieves default theme
  rt <- p_DotPlot + 
    	theme(plot.title = element_text(face = "bold", size = 12, vjust = 3)) +
    	theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=10)) +
    	theme(legend.key.height = unit(0.3, 'cm')) + 
    	theme(legend.key.width = unit(0.3, 'cm')) + 
    	theme(legend.title=element_text(size=10, face = "bold", color="black")) +  
    	theme(legend.text=element_text(size=8, color="black")) + 
    	guides(color = guide_colorbar(title = 'AE')) +  #Average Expression
	guides(size = guide_legend(title = 'PE')) +	#Percent Expression
    	xlab('') + ylab('')   #axis labels
  return(rt)
}


theme_GO_dotplot <- function(p_GO_dotplot){
  rt <- (p_GO_dotplot) + 
    	theme(plot.title = element_text(face = "bold", size = 12, vjust = 3)) +
 	 theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=10)) +
 	 theme(legend.key.height = unit(0.3, 'cm')) + 
 	 theme(legend.key.width = unit(0.3, 'cm')) + 
 	 theme(legend.title=element_text(size=10, face = "bold", color="black")) +  
 	 theme(legend.text=element_text(size=8, color="black")) 

  return(rt)
}


theme_DoHeatmap <- function(p_DoHeatmap){
  rt <- p_DoHeatmap + 
      	theme(plot.title = element_text(size = 12, vjust = 3)) +
  	#theme(axis.text.x=element_text(size=5)) +
  	theme(axis.text.y=element_text(size=7)) +
  	theme(legend.key.height = unit(0.4, 'cm')) + 
  	theme(legend.key.width = unit(0.4, 'cm')) + 
  	theme(legend.title=element_text(size=12, face = "bold", color="black")) +
  	theme(legend.text=element_text(size=10, color="black")) +
  	labs(fill='Expression') +   
  	xlab('') + ylab('')   #axis labels
  return(rt)
}

