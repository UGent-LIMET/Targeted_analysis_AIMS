# all plotting
library(moments)
library(ggplot2)
library(car)
library(RColorBrewer)

#https://ggplot2.tidyverse.org/reference/scale_shape.html
shapes <- c(21, 24, 22, 23, 16, 17, 15, 18, 1, 2, 0, 5, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 19, 20) #25 available
shapes <- rep(shapes,100)

#default setting here, before module loaded so no error
SIZE_POINTS <- 3

#color hexa rgb codes used color palette
#library(RColorBrewer)
#my.cols <- brewer.pal(9, "Blues") #9 max present
#my.cols

plot_TICchromatogram <- function(TIC_chromatogram) {
  ggplot(TIC_chromatogram, aes(x=TIC_chromatogram$scan, y=TIC_chromatogram$scan_TIC)) + geom_line() +
    geom_line(colour = "black") +
    geom_text(aes(label = TIC_chromatogram$scan), size = 3, colour = "red", angle = 0, vjust = -.05) +
    labs(title=paste("TIC chromatogram ",filename , sep=""), x="Scan", y="Intensity") +
    theme_customgridbox()
}

plot_Counts_MZ_values <- function(TIC_chromatogram) {
  ggplot(TIC_chromatogram, aes(x=TIC_chromatogram$scan, y=TIC_chromatogram$scan_TIC)) + geom_line() +
    geom_line(colour = "black") +
    geom_text(aes(label = TIC_chromatogram$scan), size = 3, colour = "red", angle = 0, vjust = -.05) +
    labs(title=paste("MZ value counts ",filename , sep=""), x="Scan", y="Counts") +
    theme_customgridbox()
}

plot_surface <- function(roi){
  library(RColorBrewer)
  #possible colorpallette:
  #display.brewer.all()
  #http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  
  ggplot(above_fwhm, aes(x=round_time, y=approx_mz, fill=intensity)) + theme_bw() +
    geom_tile() +
    geom_text(aes(label=round(intensity, 2)), alpha = 0.8, size = 2) +
    scale_fill_gradientn(colours = rev(brewer.pal(9, "Blues")))+
    #stat_density_2d(colour="black", alpha = 0.8) +
    labs(title=paste0("Contour plot of STD-", STD_ID, " with intensities above 50% peaktop"),x="Time (min)", y="MZ (Da)", fill="Intensity")+
    theme_customnogridbox()
  
  #theme(panel.background = element_rect(fill = "blue", size = 0.5, linetype = "solid"), #brewer.pal(1, "Blues")
  #panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"), 
  #panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white")
  #)
}

plot_MZplane <- function(data_melt, ppm){
  #plot mz-chromatogram of ROI (reims data) for max 15 samples
  ggplot(data_melt, aes(x = MZ, y = value)) + 
    geom_line(aes(color = variable)) +
    geom_point(aes(color = variable, shape=variable)) +
    scale_shape_manual(values=shapes) +
    geom_vline(xintercept = standard, color = 'black', size = .1)+ #line location theoretical MZ
    labs(title=paste0("MZ window ", name_standard, " (", ppm, " ppm)"),x="MZ (Da)", y="Intensity")+    
    theme_customnogridbox()+
    theme(legend.position="right", 
          legend.title = element_text(size=0), #"variable" == samplename
          legend.text = element_text(size=8))
}

plot_MZplaneWOlegend <- function(data_melt, ppm){
  #plot mz-chromatogram of ROI (reims data)
  ggplot(data_melt, aes(x = MZ, y = value)) + 
    geom_line(aes(color = variable)) +
    geom_point(aes(color = variable, shape=variable)) +
    scale_shape_manual(values=shapes) +
    geom_vline(xintercept = standard, color = 'black', size = .1)+ #line location theoretical MZ
    labs(title=paste0("MZ window ", name_standard, " (", ppm, " ppm)"),x="MZ (Da)", y="Intensity")+    
    theme_customnogridbox()+
    theme(legend.position='none')
}

plot_linear_regression <- function(std_metadata){
  #plot trendline point + linear regression line with formula/r2
  ggplot(std_metadata, aes(x = area_ratio, y = Concentration)) + 
    geom_point(aes(area_ratio)) +
    stat_smooth(method=lm, se=TRUE)+ #with standard error show
    geom_text(x = -Inf, y = Inf, label = equation, parse = TRUE, vjust = 1.5, hjust=-0.1) + #info with form/r2 on plot, not easy on top since diff intensity; so put above x-axis
    labs(title=paste0("Trendline ", name_standard), x="Area ratio", y="Concentration (ng/ul)")+  #microsymbol nok in rbox  
    theme_customnogridbox()+
    theme(legend.position="right", 
          legend.title = element_text(size=0), 
          legend.text = element_text(size=8))
}

plot_skewkurt <- function(sampleMatrix) {
  skew <- as.data.frame(sapply(sampleMatrix,skewness))
  kurt <- as.data.frame(sapply(sampleMatrix,kurtosis))   #Pearson's measure of kurtosis, actually measures outliers (>1 stdev from mean)
  skewkurt <- cbind(skew,kurt)

  skewkurtplot <- ggplot(skewkurt)+
	geom_line(aes(x=c(1:dim(skewkurt)[1]),y=skewkurt[,1],color='red'))+
	geom_line(aes(x=c(1:dim(skewkurt)[1]),y=skewkurt[,2],color='black'))+
	geom_hline(yintercept=2,color='blue',size=1)+
	geom_hline(yintercept = -2,color='blue',size=1)+
	theme(legend.position = "right",
		legend.title = element_text(face="bold",size=11),
		legend.key = element_rect(fill=NA),
		axis.title = element_text(face="bold",size=12),
		axis.text = element_text(size=12),
		#legend.key.size=2,
		plot.title = element_text(face="bold",size=12,hjust=0.5),
		panel.background = element_rect(fill=NA),
		panel.grid.major = element_line(colour="grey80"),
		panel.border = element_rect(fill=NA))+
	labs(title="Skewness and kurtosis",
	   x="Variable",
	   y="Skewness and kurtosis value")+
	scale_color_manual(name="Measures",
					 labels=c("Kurtosis", "Skewness"),
					 values=c("black","red"))+
	guides(alpha=FALSE,color=guide_legend(override.aes = list(size=2)))
}

plot_pca <- function(scores, df) {
  ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(df)))+ #substring(rownames(df),12,16))) + #todo: adjust to funn if shorter: substring(rownames(df),1,6)))
    geom_hline(yintercept = 0, colour = "gray65") +
    geom_vline(xintercept = 0, colour = "gray65") +
    geom_text(colour = "black", alpha = 0.5, size = 2) +
	theme_customgridbox()+
    ggtitle("PCA plot")
}

plot_oplsda <- function(oplsda_scores, comp, opls_comp, HotellingellipseOPLSDA) {

  ggplot(oplsda_scores,aes(x = oplsda_scores[,1],y=oplsda_scores[,2]))+
    stat_ellipse(aes(x=oplsda_scores[,1],y=oplsda_scores[,2],color=comp),size=1,
                 type="t",level=0.95,show.legend = FALSE,segments=200)+
    geom_path(data=HotellingellipseOPLSDA,aes(x=HotellingellipseOPLSDA$x,y=HotellingellipseOPLSDA$y),size=1,color='black',linejoin = 'round')+
    geom_point(aes(fill=comp,shape=comp),color='black',size=SIZE_POINTS)+
    theme_customnogridbox()+
    labs(title=paste0("OPLS-DA score plot ", pairwise_comparison), #title="", 
         x=paste("Predictive component 1 (",format(round(opls_comp@modelDF$R2X[1]*100,digits=1),nsmall = 1),"% explained var.)",sep=""),
         y=paste("Orthogonal component 1 (",format(round(opls_comp@modelDF$R2X[2]*100, digits=1),nsmall = 1),"% explained var.)",sep=""),
         fill="Group",color="Group",shape="Group")+
    #geom_label(aes(x=Inf,y=Inf),size=5, hjust=1,vjust=1,
    #           label=paste("Q2Y(cum) = ",format(round(opls_comp@summaryDF$`Q2(cum)`, digits=3),nsmall = 3),
    #                   "\nR2X(cum) = ",format(round(opls_comp@summaryDF$`R2X(cum)`, digits=3),nsmall = 3),
    #                   "\nR2Y(cum) = ",format(round(opls_comp@summaryDF$`R2Y(cum)`, digits=3),nsmall = 3),sep=""))+
    scale_shape_manual(values=shapes)+
    scale_fill_brewer(type='seq',palette='Blues')+
    scale_color_brewer(type='seq',palette='Blues')+
    guides(alpha=FALSE)
    #labs(x =bquote(Q^2~Y(cum) == .(format(opls_comp@summaryDF$`Q2(cum)`),digits=3)))#+ #works only as lab, not in label      #expression(paste("Q"^2,"Y(cum) = ",sep="")),
}

plot_vip <- function(vipdf) {
  ggplot(vipdf,aes(x=vipdf[,"pvip"],y=vipdf[,"vip_comp"]))+
    geom_point(aes(alpha=0.9,color='black'))+
    geom_hline(yintercept = 1, color = 'red', size = 1)+
    geom_vline(xintercept = 0.05, color = 'red', size = 1)+
    theme_customgridbox()+
    labs(title="VIP-plot",x="p-value", y="VIP")+
    scale_color_manual(values="black")+
    theme(legend.position="none") #no legend
}

plot_Splot <- function(Splotframe_comp) {
  ggplot(Splotframe_comp,aes(x=Splotframe_comp[,1],y=Splotframe_comp[,2]))+
    geom_point(aes(alpha=0.9,color=Corr_cutoff))+
    geom_hline(yintercept = Cutoffvalue_corr,color='red',size=1) +
    geom_hline(yintercept = -Cutoffvalue_corr,color='red',size=1)+
    theme_customgridbox()+
    labs(title="S-plot",x="Covariance", y="Correlation")+
    scale_color_manual(values=c("black","red"))+
    theme(legend.position="none") #no legend
}

plot_loading <- function(myCIsframe_comp) {
  ggplot(myCIsframe_comp,aes(x=c(1:(length(samples_matrix_comp_no0)-0)),y=myCIsframe_comp[,1]))+
    geom_pointrange(mapping=aes(x=c(1:(length(samples_matrix_comp_no0)-0)),y=myCIsframe_comp[,1],ymin=myCIsframe_comp[,2],ymax=myCIsframe_comp[,3]))+
    theme_customgridbox()+
    labs(title="Loading plot",x="Variable", y="Loading",color="Groups",shape="Groups")+
    guides(alpha=FALSE,color=guide_legend(override.aes = list(size=3))) 
}

plot_ConfusionMatrix_comp <- function(confusion_matrix_long){
  #makes confusion matrix with actual vs predicted values
  ggplot(data =  confusion_matrix_long, mapping = aes(x = ActualValid, y = predValid)) +
    geom_tile(aes(fill = value), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f",value)), vjust = 1) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title=paste0("Confusion matrix comparison ", pairwise_comparison),x="Actual", y="Predicted")+    
    theme_customnogridbox()+
    theme(legend.position='none')
}

plot_plsda <- function(plsda_scores, comp, pls_comp, HotellingellipsePLSDA) {
  ggplot(plsda_scores,aes(x = plsda_scores[,1],y=plsda_scores[,2]))+
    stat_ellipse(aes(x=plsda_scores[,1],y=plsda_scores[,2],color=comp),size=1,
		             type="t",level=0.95,show.legend = FALSE,segments=200)+
    geom_path(data=HotellingellipsePLSDA,aes(x=HotellingellipsePLSDA$x,y=HotellingellipsePLSDA$y),size=1,color='black',linejoin = 'round')+
    geom_point(aes(fill=comp,shape=comp),color='black',size=SIZE_POINTS)+
    theme_customnogridbox()+
    labs(title=paste0("PLS-DA score plot ", multiple_comparison), #title="", #
    	   x=paste("Predictive component 1 (",format(round(pls_comp@modelDF$R2X[1]*100,digits=1),nsmall = 1),"% explained var.)",sep=""),
    	   y=paste("Predictive component 2 (",format(round(pls_comp@modelDF$R2X[2]*100,digits=1),nsmall = 1),"% explained var.)",sep=""),
    	   fill="Group",color="Group",shape="Group")+
    #geom_label(aes(x=Inf,y=Inf),size=5,hjust=1,vjust=1,
		#	 label=paste("Q2Y(cum) = ",format(round(pls_comp@summaryDF$`Q2(cum)`, digits=3),nsmall = 3),
		#				 "\nR2X(cum) = ",format(round(pls_comp@summaryDF$`R2X(cum)`, digits=3),nsmall = 3),
		#				 "\nR2Y(cum) = ",format(round(pls_comp@summaryDF$`R2Y(cum)`, digits=3),nsmall = 3),sep=""))+
    scale_shape_manual(values=shapes)+
    scale_fill_brewer(type='seq',palette='Blues')+
    scale_color_brewer(type='seq',palette='Blues')+
    guides(alpha=FALSE)
}

plot_ConfusionMatrix_mcomp <- function(confusion_matrix_long){
  #makes confusion matrix with actual vs predicted values
  ggplot(data =  confusion_matrix_long, mapping = aes(x = ActualValid, y = predValid)) +
    geom_tile(aes(fill = value), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f",value)), vjust = 1) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title=paste0("Confusion matrix multiple comparison ", multiple_comparison),x="Actual", y="Predicted")+    
    theme_customnogridbox()+
    theme(legend.position='none')
}

plot_pcascores <- function(PCAscores_comp, comp) {
  ggplot(PCAscores_comp,aes(x = PCAscores_comp[,1],y=PCAscores_comp[,2]))+
    stat_ellipse(aes(x=PCAscores_comp[,1],y=PCAscores_comp[,2],color=comp),size=1,
                 type="t",level=0.95,show.legend = FALSE,segments=200)+
    geom_path(data=HotellingellipsePCA_comp,aes(x=HotellingellipsePCA_comp$x,y=HotellingellipsePCA_comp$y),size=1,color='black',linejoin = 'round')+
    geom_point(aes(fill=comp,shape=comp),color='black',size=SIZE_POINTS)+
    theme_customnogridbox()+
    labs(title=paste("PCA score plot ",projection , sep=""), #title="", 
         x=paste("PC 1 (",format(round(PCAoptnumber@R2[1]*100,digits=1),nsmall = 1),"% explained var.)",sep = ""),
         y=paste("PC 2 (",format(round(PCAoptnumber@R2[2]*100,digits=1),nsmall = 1),"% explained var.)",sep = ""),
         fill="Group",color="Group",shape="Group")+
    scale_shape_manual(values=shapes)+ 
    scale_fill_brewer(type='seq',palette='Blues')+
    scale_color_brewer(type='seq',palette='Blues')+
    guides(alpha=FALSE) 
}

Hotellingellipse <- function(samples,scorematrix){
  library(car)
  Ellipseaxisscore1<-var(scorematrix[,1])*
    qf(0.95,2,samples-2)*
    ((2*(samples^2-1))/sqrt(samples*(samples-2)))*0.05
  Shape1<-c(Ellipseaxisscore1,0)  #Divide by number of samples to get on the right scale?
  Ellipseaxisscore2<-var(scorematrix[,2])*
    qf(0.95,2,samples-2)*
    ((2*(samples^2-1))/sqrt(samples*(samples-2)))*0.05
  Shape2<-c(0,Ellipseaxisscore2)
  Shapematrix=cbind(Shape1,Shape2)
  Ellipse1<-as.data.frame(ellipse(c(mean(scorematrix[,1]), mean(scorematrix[,2])), shape=Shapematrix, radius=1,segments=200,add=FALSE))
  return(Ellipse1)
}

plot_heatmap_hclust <- function(data_melt){
  #make heatmap, colorblind friendly, order samplesaccording to comp + gropu sec axis
  library(RColorBrewer)
  colors<-brewer.pal(11,name="RdYlBu")
  pal<-colorRampPalette(colors)
  
  ggplot(data_melt, aes(x=SampleName, y=variable)) +    
    geom_tile(aes(fill = value)) +    
    scale_x_discrete (limits = comp_names2$SampleName[ord]) +  #according to clustering!
    scale_fill_gradientn(colours=brewer.pal(n=9, name="PuBuGn")) +  #https://www.r-bloggers.com/2013/10/creating-colorblind-friendly-figures/
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=SIZE_HEATMAP_YAXIS), 
          #text = element_text(size = 1),
          axis.text.y = element_text(size = SIZE_HEATMAP_YAXIS))+    #https://statisticsglobe.com/change-font-size-of-ggplot2-plot-in-r-axis-text-main-title-legend 
    labs(x="SampleName", y="Variable", fill="Value") #1 funct for both comp/Mcomp
}

plot_heatmap_w_group <- function(data_melt){
  #make heatmap, colorblind friendly, order samplesaccording to comp + gropu sec axis
  library(RColorBrewer)
  colors<-brewer.pal(11,name="RdYlBu")
  pal<-colorRampPalette(colors)
  
  ggplot(data_melt, aes(x=order_plot, y=variable)) +    
    geom_tile(aes(fill = value)) +    
    scale_x_continuous(breaks = 1:length(SampleName_reorder),
                       labels = SampleName_reorder,
                       sec.axis = dup_axis(name = "Group",
                                           breaks = 1:length(SampleName_reorder),
                                           labels = comp_names$comp))+    
    scale_fill_gradientn(colours=brewer.pal(n=9, name="PuBuGn")) +  #https://www.r-bloggers.com/2013/10/creating-colorblind-friendly-figures/
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=SIZE_HEATMAP_YAXIS), 
          axis.text.y = element_text(size = SIZE_HEATMAP_YAXIS))+    #https://statisticsglobe.com/change-font-size-of-ggplot2-plot-in-r-axis-text-main-title-legend 
    labs(x="SampleName", y="Variable", fill="Value") #1 funct for both comp/Mcomp
}

plot_duoboxplot <- function(samples_matrix_comp_no0) {
  #https://statisticsglobe.com/boxplot-in-r
  
  ggplot(samples_matrix_comp_no0, aes(x = comp, y = samples_matrix_comp_no0[,standard_nr], fill = comp)) +    # Create boxplot chart in ggplot2; 
    geom_boxplot() +
    labs(title="", #title=paste0("Boxplot comparison ", pairwise_comparison), 
         x="Pairwise comparison", 
         y=name_standard)+
    scale_fill_brewer(type='seq',palette='Blues')+
    theme_customgridbox()+
    theme(legend.position="none") #no legend
}

plot_multiboxplot <- function(samples_matrix_comp_no0) {
  #https://statisticsglobe.com/boxplot-in-r
  
  ggplot(samples_matrix_comp_no0, aes(x = comp, y = samples_matrix_comp_no0[,standard_nr], fill = comp)) +    # Create boxplot chart in ggplot2; 
    geom_boxplot() +
    labs(title="", #title=paste0("Boxplot multiple comparison ", multiple_comparison), 
         x="Multiple comparison", 
         y=name_standard)+
    scale_fill_brewer(type='seq',palette='Blues')+
    theme_customgridbox()+
    theme(legend.position="none") #no legend
}

plot_volcano <- function(volcano_df){
  #https://erikaduan.github.io/posts/2021-01-02-volcano-plots-with-ggplot2/
  
  ggplot(volcano_df, aes(x = log2(as.numeric(Fold.change)+0.000001), y = -log10(as.numeric(p.value)+0.000001) )) + 
    geom_point() +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed") + 
    geom_vline(xintercept = c(log2(0.5), log2(2)),
               linetype = "dashed") +
    labs(title= "",  #"p-value versus fold change",  
         x="log2 fold change", 
         y="-log10 p-value") +
    theme_customnogridbox()+
    theme(legend.position="none") #no legend
}

plot_boxplot <- function(corr_info) {
  #https://statisticsglobe.com/boxplot-in-r
  #https://www.tutorialgateway.org/r-ggplot2-jitter/
  
  ggplot(corr_info, aes(x = "Varset1 vs Varset2", y = corr_info$ABSrho)) +    # Create boxplot chart in ggplot2; 
    #geom_jitter(position=position_jitter(0.2), alpha = 0.8, color="#9ECAE1")+
    geom_boxplot(fill="#9ECAE1", outlier.color="black") +
    labs(title="", #title=paste0("Boxplot comparison ", pairwise_comparison), 
         x="", 
         y="Rho") +
    scale_fill_brewer(type='seq',palette='Blues')+
    theme_customgridbox()+
    theme(legend.position="none") #no legend
}

plot_dotplot <- function(corr_info_top) {
  #https://www.tutorialgateway.org/r-ggplot2-jitter/
  
  ggplot(corr_info_top, aes(x = "Varset1 vs Varset2", y = corr_info_top$ABSrho)) +    # Create boxplot chart in ggplot2; 
    geom_jitter(alpha = 0.8, color="#9ECAE1", height = 0, width = 0.1)+
    #geom_jitter(position=position_jitter(0.1), alpha = 0.8, color="#9ECAE1")+ #no! get other value on y-axis!
    geom_point(alpha = 0.8, color="#9ECAE1")+ #overlay ALL points if same value
    #geom_boxplot(fill="#9ECAE1", outlier.color="black") +
    labs(title="", #title=paste0("Boxplot comparison ", pairwise_comparison), 
         x="", 
         y="Rho") +
    scale_fill_brewer(type='seq',palette='Blues')+
    theme_customgridbox()+
    theme(legend.position="none") #no legend
}

theme_customgridbox<- function () {
  theme(legend.position = "right",
  legend.title = element_text(face="bold",size=18),
  legend.text = element_text(size=18),
  legend.key = element_rect(fill=NA),
  axis.title = element_text(size=16),
  axis.text = element_text(size=16,colour="black"),
  axis.ticks = element_line(colour = "black"),
  plot.title = element_text(face="bold",size=18,hjust=0.5),
  panel.background = element_rect(fill=NA),
  panel.grid.major = element_line(colour="grey80"),
  panel.border = element_rect(fill=NA,size=1))
}

theme_customnogridbox<- function () {
  theme(legend.position = "right",
  legend.title = element_text(face="bold",size=18),
  legend.text = element_text(size=18),
  legend.key = element_rect(fill=NA),
  axis.title = element_text(size=16),
  axis.text = element_text(size=16,colour="black"),
  axis.ticks = element_line(colour = "black"),
  plot.title = element_text(face="bold",size=18,hjust=0.5),
  panel.background = element_rect(fill=NA),
  panel.border = element_rect(fill=NA,size=1))
}


  
