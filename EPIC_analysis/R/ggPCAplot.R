#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# R function to create ggplot of PCAs
#
# author. Ayden Saffari <ayden.saffari@lshtm.ac.uk>
# affiliations. MRC ING, LSHTM, BBK
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
pal <- c("#F87E86","#0E7597","#2E7B39","#FDA80A","#B31918","#9CD893","#E6BD9B","#B7BAAA","#683A2E","#211318")
ggPCAplot <- function(df,pcs=c(1,2),batch="Sex",xlim = c(-500,500),ylim = c(-500,500),pal=NA) {

	#x and y labels
	xlab <- paste("PC",pcs[1],sep="")
	ylab <- paste("PC",pcs[2],sep="")
	
	df <- data.frame(PC1=pmax(xlim[1], pmin(xlim[2], 
        df[,which(colnames(df) == xlab)])),PC2=pmax(ylim[1], pmin(ylim[2], 
        df[,which(colnames(df) == ylab)])),batch=paste(batch,df[,which(colnames(df) == batch)],sep="_"))
        
	#get factor levels
	levs <- levels(df[,which(colnames(df) == batch)])
	
    #set shapes
    df$pch <- ifelse(
    	   df$PC1 == xlim[1] | df$PC2 == ylim[1], "down",
    	   ifelse(df$PC1 == xlim[2] | df$PC2 == ylim[2], "up", "norm"))
    
    #plot
    ggplot(data=df,aes(PC1, PC2)) + 
    	   geom_point(aes(fill=batch, color=batch, shape=pch)) +
           xlim(xlim[1],xlim[2]) +
           ylim(ylim[1],ylim[2]) +
           scale_fill_manual(values=pal) +
           scale_color_manual(values=pal) +
    	   scale_shape_manual(values=c("down"=25,"up"=24,"norm"=21)) +           
    	   xlab(xlab) +
    	   ylab(ylab) +
    	   guides(shape=FALSE)

}
