library(dplyr,warn.conflicts=F)
library(pheatmap)

args<-commandArgs(T)
hm_file <- args[1]
anim_file <- args[2]
disease_csv <- args[3]
prefix <- args[4]


print(paste(prefix,'Start'))

hm_csv <- read.csv(hm_file,header=T,stringsAsFactors=F)
anim_csv <- read.csv(anim_file, header=T,stringsAsFactors=F)
disease_csv <- read.csv(disease_csv,header=T,stringsAsFactors=F)
colnames(disease_csv)[1] <- 'gene'

hm_csv$rowMeans.df2.<-log2(hm_csv$rowMeans.df2.)
anim_csv$rowMeans.df2.<-log2(anim_csv$rowMeans.df2.)



disease <- sort(unique(disease_csv$disease))
ctype <- sort(unique(intersect(hm_csv$type,anim_csv$type)))

### Function ###
selectVal<-function(mydata,disease_csv,cell,d){
         d_gene <- filter(disease_csv,disease==d)$gene
         g_val <- filter(mydata,type==cell) %>% filter(Human.gene.name %in% d_gene) %>% 
                         select(Human.gene=Human.gene.name,value=rowMeans.df2.)
         return(g_val)
}
#####
rho_tab <- list()
for(cell in ctype){
	rho_tab[[cell]]<-c()
	for(d in disease){
		hm_g_val <- selectVal(hm_csv, disease_csv, cell, d)
		anim_g_val <- selectVal(anim_csv, disease_csv, cell, d)
		hm_anim_val<-merge(hm_g_val,anim_g_val,by=1)
		rho<-tryCatch({cor.test(
		  hm_anim_val$value.x,
		  hm_anim_val$value.y,
		  method='spearman')$estimate},
				error=function(e){print(paste(cell,disease));print(hm_anim_val)}
			      )
		rho_tab[[cell]]<-c(rho_tab[[cell]],rho)
	}
}
rho_tab <- data.frame(rho_tab)
rownames(rho_tab) <- disease
write.csv(rho_tab,paste0(prefix,'_spearman.csv'))

pheatmap(rho_tab,cluster_cols=F,cluster_rows=F,border_color=NA,
         fontsize=8,cellheight=15,cellwidth=15,
         width=5,height=5,filename = paste0(prefix,'_spearman_heatmap.pdf'))

print('Done')

#cell<-'AC'
#d<-"Diabetic retinopathy"
library(ggplot2)
test<-merge(hm_csv,anim_csv,by=c('Human.gene.name','type'))
test2<-merge(test,disease_csv,by=1)[,c(1,2,3,6,7)]
colnames(test2)<-c('gene','cell','hm_val','anim_val','disease')

text_loc <- data.frame(
  label=paste0('r=',round(as.vector(as.matrix(rho_tab)),2)),
  #label=paste0('r=',1),
  x=rep(-6,length(ctype)*length(disease)),
  y=rep(6,length(ctype)*length(disease)),
  cell=rep(ctype,each=length(disease)),
  disease=rep(disease,length(ctype)),
  stringsAsFactors = F
)

pic<-ggplot(test2,aes(x=anim_val,y=hm_val))+geom_point(size=0.5)+
  stat_smooth(method = lm,formula = y ~ poly(x,1),colour='red',se=F)+
  facet_grid(cell~disease)+
  xlim(-10,8)+ylim(-10,8)+
  theme(
    title = element_blank(),
    panel.spacing = unit(0,'cm'),
    strip.text.x=element_text(size=3),
    axis.text.x = element_text(angle=-90, hjust=0, vjust=0.25),
    axis.text = element_text(size=7),
    panel.background=element_blank(),
    panel.border = element_rect(colour = "black",fill='transparent')
  )+geom_text(data=text_loc,mapping = aes(x=x,y=y,label=label),size = 3,check_overlap=T)
ggsave(paste0('Fig5D_',prefix,'_rho.pdf'),pic,width = 15,height = 8)
