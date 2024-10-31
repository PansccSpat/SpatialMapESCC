source('~/xenium/header.R')
source('~/xenium/xenium_function.R')



##### input #####
job({
  xe1_1 = read_rds('data/xenium_seurat_processed/xenium_t1_s1_seurat.rds')
  xe1_1@meta.data = read_rds('data/xenium_seurat_processed/xenium_t1_s1_meta.rds')
  xe1_2 = read_rds('data/xenium_seurat_processed/xenium_t1_s2_seurat.rds')
  xe1_2@meta.data = read_rds('data/xenium_seurat_processed/xenium_t1_s2_meta.rds')
  xe1_3 = read_rds('data/xenium_seurat_processed/xenium_t1_s3_seurat.rds')
  xe1_3@meta.data = read_rds('data/xenium_seurat_processed/xenium_t1_s3_meta.rds')
  
  xe2_1 = read_rds('data/xenium_seurat_processed/xenium_t2_s1_seurat.rds')
  xe2_1@meta.data = read_rds('data/xenium_seurat_processed/xenium_t2_s1_meta.rds')
  xe2_2 = read_rds('data/xenium_seurat_processed/xenium_t2_s2_seurat.rds')
  xe2_2@meta.data = read_rds('data/xenium_seurat_processed/xenium_t2_s2_meta.rds')
  xe2_3 = read_rds('data/xenium_seurat_processed/xenium_t2_s3_seurat.rds')
  xe2_3@meta.data = read_rds('data/xenium_seurat_processed/xenium_t2_s3_meta.rds')
  
  xe3_1 = read_rds('data/xenium_seurat_processed/xenium_t3_s1_seurat.rds')
  xe3_1@meta.data = read_rds('data/xenium_seurat_processed/xenium_t3_s1_meta.rds')
  xe3_2 = read_rds('data/xenium_seurat_processed/xenium_t3_s2_seurat.rds')
  xe3_2@meta.data = read_rds('data/xenium_seurat_processed/xenium_t3_s2_meta.rds')
  
  xe4_1 = read_rds('data/xenium_seurat_processed/xenium_t4_s1_seurat.rds')
  xe4_1@meta.data = read_rds('data/xenium_seurat_processed/xenium_t4_s1_meta.rds')
  xe4_2 = read_rds('data/xenium_seurat_processed/xenium_t4_s2_seurat.rds')
  xe4_2@meta.data = read_rds('data/xenium_seurat_processed/xenium_t4_s2_meta.rds')
  
  xe5_1 = read_rds('data/xenium_seurat_processed/xenium_t5_s1_seurat.rds')
  xe5_1@meta.data = read_rds('data/xenium_seurat_processed/xenium_t5_s1_meta.rds')
  xe5_2 = read_rds('data/xenium_seurat_processed/xenium_t5_s2_seurat.rds')
  xe5_2@meta.data = read_rds('data/xenium_seurat_processed/xenium_t5_s2_meta.rds')
  xe5_3 = read_rds('data/xenium_seurat_processed/xenium_t5_s3_seurat.rds')
  xe5_3@meta.data = read_rds('data/xenium_seurat_processed/xenium_t5_s3_meta.rds')
  xe5_2_split = read_rds('data/xenium_seurat_processed/xenium_t5_s2_split_seurat.rds')
  xe5_2_split@meta.data = read_rds('data/xenium_seurat_processed/xenium_t5_s2_split_meta.rds')
  xe5_3_split = read_rds('data/xenium_seurat_processed/xenium_t5_s3_split_seurat.rds')
  xe5_3_split@meta.data = read_rds('data/xenium_seurat_processed/xenium_t5_s3_split_meta.rds')
  
  xe6_1 = read_rds('data/xenium_seurat_processed/xenium_t6_s1_seurat.rds')
  xe6_1@meta.data = read_rds('data/xenium_seurat_processed/xenium_t6_s1_meta.rds')
  xe6_2 = read_rds('data/xenium_seurat_processed/xenium_t6_s2_seurat.rds')
  xe6_2@meta.data = read_rds('data/xenium_seurat_processed/xenium_t6_s2_meta.rds')
  xe6_3 = read_rds('data/xenium_seurat_processed/xenium_t6_s3_seurat.rds')
  xe6_3@meta.data = read_rds('data/xenium_seurat_processed/xenium_t6_s3_meta.rds')
},import = 'auto')



cd4marker = c('CD4','CD8A','CD8B',
              'FOXP3','IL2RA','CTLA4', # reg
              'TCF7','CCR7','SELL', # tn
              'IL7R','CCL5', # mem
              'CXCL13','CD200' # tfh
)
cd4marker1 = c('CXCL13','CD200', # tfh
               'FOXP3','IL2RA','CTLA4' # reg
)

cd8marker = c('CD4','CD8A','CD8B',
              'TCF7','CCR7','SELL', # tn
              'IL7R','CXCR4',#mem
              'GNLY', 'NKG7', 'GZMB', 'GZMA', 'GZMH', 'GZMK', #Teff
              'PDCD1','CXCL13','CTLA4','HAVCR2' # tex 
)
cd8marker1 = c('GNLY', 'NKG7', 'GZMB', 'GZMA', 'GZMH', #Teff
              'PDCD1','CXCL13','CTLA4','HAVCR2' # tex 
)

tmarker = c('CD200', # tfh
            'GNLY', 'NKG7', 'GZMB', 'GZMA', 'GZMH', #Teff
            'FOXP3','IL2RA','CTLA4', # reg
            'PDCD1','CXCL13','HAVCR2' # tex 
            )


cytogene = c('GNLY', 'NKG7', 'GZMB', 'GZMA', 'GZMH', 'GZMK', 'KLRB1','TYROBP')  
exhgene = c('PDCD1','HAVCR2','CTLA4','TIGIT','LAG3','CXCL13') 
reggene = c('FOXP3','TIGIT','CD274','IL2RA')

fibmarker = c('PLA2G2A','GPX3','IGFBP6','PI16','ADH1B','TNXB','DCN','CXCL12', # NF
              'IGF1','RGS5','HIF1A','CCND1','POSTN','STAT1','FAP','SERPINE1','DPT', # iCAF
              'MMP11','MMP1','IL7R','STMN1','PDPN','RUNX1' # myCAF
)



primcol = brewer.pal(10, 'Paired')
primcol[7] = '#999999'
names(primcol) = c('B','Epi','Endo','Fib','Gland','Mye','Mus','Mast','Pla','T')

primcol2 = c('#1f78b4','#ffc089','#B5EFB5','#793c1b',
             '#6a3d9a','#333333','#ffff33','#df928e')
names(primcol2) = c('B','Endo','Fib','Gland','Mye','Mus','T','Epi')

primcol3 = c('#793c1b','#333333','#ffc089','#df928e',
             '#1f78b4','#6a3d9a','#df65b0','#1cbe4f','#e31a1c','#ffdcbd',
             '#BDCDFF','#ffff68','#325A9B','#AAF400','#ff9344')
names(primcol3) = c('Gland','Mus','Endo','Epi',
                    'B','Mye','NF-PI16','NF','CAF','T-others',
                    'T-naive','CD8T-eff','CD8T-ex','CD4T-fh','CD4T-reg')
areacol = c('Infiltrated'='#fccde5','Out'='#80b1d3','Others'='#888888')

tmpcol = c('#df928e','#888888','#1cbe4f','#238b45','#B5EFB5',
           '#df65b0','#ef3b2c','#a50f15')
names(tmpcol) = c('Epi','Others','CD8T-n/m','CD8T-eff','CD8T-ex',
                  'CD4T-fh','CD4T-reg','CD4T-n/m')

allcellcol = c('Basal'='#0339f8','Proliferation'='#ff9408','Differentiation'='#75bbfd',
               'Terminal'='#4682B4','HT_Prolif'='#ff9408','HT_Diff/Term'='#4682B4',
               'Invasive'='#de0c62','T cells'='#96f97b','B cells'='#028f1e',
               'Myeloid cells'='#7e1e9c','Fibroblast'='#696969','CAF'='#FFFFFF',
               'Endothelial cells'='#7b002c','Gland cells'='#c69f59','Myocytes'='#1b2431')

tccol = c('CD8T-n/m'='#BDCDFF','CD8T-eff'='#ffff68','CD8T-ex'='#325A9B',
          'CD4T-n/m'='#beaed4','CD4T-fh'='#AAF400','CD4T-reg'='#ff9344')

getp = function(srt, features, group.by){
  p1 = (DimPlot(srt, group.by = group.by, label = T)+NoLegend())+DimPlot(srt, group.by = group.by, label = T)
  p2 = DotPlot(srt, features = features, group.by = group.by)+RotatedAxis()
  return(p1/p2)
}

srtpreproc = function(srt){
  srt = SCTransform(srt)
  srt = RunPCA(srt)
  srt = RunUMAP(srt, dims = 1:10)
  srt = FindNeighbors(srt,dims = 1:10)
  srt = FindClusters(srt, resolution = seq(0.1,1,0.1))
  return(srt)
}


allmeta = read_rds('data/xenium_seurat_processed/allcell_meta.rds')
epipoly = read_rds('stp3_tme/data/epiregion/allsample_poly.rds')

srttc = read_rds('stp3_tme/data/tc/all_seurat.rds')
srttc@meta.data = read_rds('stp3_tme/data/tc/all_meta.rds')

newmt = c('cellID','cellType_merge')
newmt = rbind(srttc@meta.data[,newmt],srtfib@meta.data[,newmt])
newmt$cellType_merge = newmt$cellType_merge %>% 
  mapvalues(c('CD8T-n/m','CD4T-n/m'),c('T-naive','T-naive'))

srttc$celltype = srttc$celltype %>%
  factor(levels = c('CD4T-n/m','CD4T-fh','CD4T-reg','CD8T-n/m','CD8T-eff','CD8T-ex','T-others'))
srttc$cellType_merge = srttc$cellType_merge %>%
  factor(levels = c('T-others','T-naive','CD4T-fh','CD4T-reg','CD8T-eff','CD8T-ex'))

# job({
#   cd4 = subset(srttc, cells=srttc$cellID[grep('CD4',srttc$celltype)])
#   cd8 = subset(srttc, cells=srttc$cellID[grep('CD8',srttc$celltype)])
#   cd4 = cd4 %>% SCTransform(ncells=3000) %>% RunPCA %>% RunUMAP(dims = 1:10)
#   cd8 = cd8 %>% SCTransform(ncells=3000) %>% RunPCA %>% RunUMAP(dims = 1:10)
#   saveRDS(cd4@meta.data, 'stp3_tme/data/tc/cd4_meta.rds')
#   saveRDS(cd4, 'stp3_tme/data/tc/cd4_seurat.rds')
#   saveRDS(cd8@meta.data, 'stp3_tme/data/tc/cd8_meta.rds')
#   saveRDS(cd8, 'stp3_tme/data/tc/cd8_seurat.rds')
# },import = 'auto')
cd4 = read_rds('stp3_tme/data/tc/cd4_seurat.rds')
cd4@meta.data = read_rds('stp3_tme/data/tc/cd4_meta.rds')
cd8 = read_rds('stp3_tme/data/tc/cd8_seurat.rds')
cd8@meta.data = read_rds('stp3_tme/data/tc/cd8_meta.rds')




##### umap #####
# DimPlot(srttc,group.by = 'cellType_merge',cols = primcol3)
cd4$cellType_merge = cd4$cellType_merge %>% factor(levels = c('CD4T-reg','CD4T-fh','T-naive'))
cd8$cellType_merge = cd8$cellType_merge %>% factor(levels = c('CD8T-ex','CD8T-eff','T-naive'))

mydimp(cd4,groupby = 'cellType_merge',cols = primcol3)+NoLegend()
ggsave('stp3_tme/plot/tc_analysis/umap_cd4.pdf',width = 4,height = 4)
mydimp(cd8,groupby = 'cellType_merge',cols = primcol3)+NoLegend()
ggsave('stp3_tme/plot/tc_analysis/umap_cd8.pdf',width = 4,height = 4)


names(areacol) = c('Out','Infiltrated')
mydimp(cd4,groupby = 'epiarea',cols = areacol)+NoLegend()
ggsave('stp3_tme/plot/tc_analysis/umap_area_cd4.pdf',width = 4,height = 4)
mydimp(cd8,groupby = 'epiarea',cols = areacol)+NoLegend()
ggsave('stp3_tme/plot/tc_analysis/umap_area_cd8.pdf',width = 4,height = 4)

mydimp(cd4,groupby = 'cellType_merge',cols = primcol3,splitby = 'epiarea')+NoLegend()
ggsave('stp3_tme/plot/tc_analysis/umap_splitarea_cd4.pdf',width = 8,height = 4)
mydimp(cd8,groupby = 'cellType_merge',cols = primcol3,splitby = 'epiarea')+NoLegend()
ggsave('stp3_tme/plot/tc_analysis/umap_splitarea_cd8.pdf',width = 8,height = 4)

mydimp(cd4,groupby = 'cellType_merge',cols = primcol3,splitby = 'stage1')+NoLegend()
ggsave('stp3_tme/plot/tc_analysis/umap_splitstage_cd4.pdf',width = 6,height = 6)
mydimp(cd8,groupby = 'cellType_merge',cols = primcol3,splitby = 'stage1')+NoLegend()
ggsave('stp3_tme/plot/tc_analysis/umap_splitstage_cd8.pdf',width = 6,height = 6)




##### t marker dotp #####
{# v1
  cd4$cellType_merge = cd4$cellType_merge %>% factor(levels = rev(c('CD4T-reg','CD4T-fh','T-naive')))
  cd8$cellType_merge = cd8$cellType_merge %>% factor(levels = rev(c('CD8T-ex','CD8T-eff','T-naive')))
  mydotp(srttc,features = c('CD4','CD8A','CD8B'),groupby = 'celltype',cols = rev(brewer.pal(9,'RdBu')))+
    RotatedAxis()+theme(axis.line = element_blank(),panel.border = element_rect(fill=NA,color='black',linewidth=0.5))
  mydotp(cd4,features = cd4marker1,groupby = 'cellType_merge',cols = rev(brewer.pal(9,'RdBu')))+
    RotatedAxis()+NoLegend()+xlab('')+ylab('')+
    theme(axis.line = element_blank(),panel.border = element_rect(fill=NA,color='black',linewidth=0.5))->p1
  mydotp(cd8,features = cd8marker1,groupby = 'cellType_merge',cols = rev(brewer.pal(9,'RdBu')))+
    RotatedAxis()+NoLegend()+xlab('')+ylab('')+scale_y_discrete(position = 'right')+
    theme(axis.line = element_blank(),panel.border = element_rect(fill=NA,color='black',linewidth=0.5))->p2
  p1+p2
  ggsave('stp3_tme/plot/tc_analysis/tc_marker_dotp.pdf',width = 6,height = 4)
}

{# v2
  cd4$celltype = cd4$celltype %>% factor(levels = rev(c('CD4T-reg','CD4T-fh','CD4T-n/m')))
  cd8$celltype = cd8$celltype %>% factor(levels = rev(c('CD8T-ex','CD8T-eff','CD8T-n/m')))
  mydotp(srttc,features = c('CD4','CD8A','CD8B'),groupby = 'celltype',cols = rev(brewer.pal(9,'RdBu')))+
    RotatedAxis()+theme(axis.line = element_blank(),panel.border = element_rect(fill=NA,color='black',linewidth=0.5))
  mydotp(cd4,features = cd4marker1,groupby = 'celltype',cols = rev(brewer.pal(9,'RdBu')))+
    RotatedAxis()+NoLegend()+xlab('')+ylab('')+
    theme(axis.line = element_blank(),panel.border = element_rect(fill=NA,color='black',linewidth=0.5))->p1
  mydotp(cd8,features = cd8marker1,groupby = 'celltype',cols = rev(brewer.pal(9,'RdBu')))+
    RotatedAxis()+NoLegend()+xlab('')+ylab('')+scale_y_discrete(position = 'right')+
    theme(axis.line = element_blank(),panel.border = element_rect(fill=NA,color='black',linewidth=0.5))->p2
  p1+p2
  ggsave('stp3_tme/plot/tc_analysis/tc_marker_dotp.pdf',width = 6,height = 4)
}

{# v3
  cd4$cellType_merge = cd4$cellType_merge %>% factor(levels = rev(c('CD4T-reg','CD4T-fh','CD4T-n/m')))
  cd8$cellType_merge = cd8$cellType_merge %>% factor(levels = rev(c('CD8T-ex','CD8T-eff','CD8T-n/m')))
  DotPlot(cd4,features = cd4marker1,group.by = 'cellType_merge',cols = 'RdBu')+
    RotatedAxis()+xlab('')+ylab('')
  ggsave('stp3_tme/plot/tc_analysis/marker_dotp_cd4.pdf',width = 4.5,height = 2.5)
  DotPlot(cd8,features = cd8marker1,group.by = 'cellType_merge',cols = 'RdBu')+
    RotatedAxis()+xlab('')+ylab('')
  ggsave('stp3_tme/plot/tc_analysis/marker_dotp_cd8.pdf',width = 5,height = 2.5)
}



##### t marker heat area #####
tmp = unique(srttc@meta.data[,c('epiarea','stage1')])
tmp$tmp = paste0(tmp$epiarea,'-',tmp$stage1)
tmp = arrange(tmp,epiarea,stage1)

anncol = data.frame(area=c(rep('Infiltrated',4),rep('Out',4)))
anncolor = list(area=c('Out'='#fccde5','Infiltrated'='#80b1d3'))

tmp1 = cd4
tmp1$tmp = paste0(tmp1$epiarea,'-',tmp1$stage1)
tmp1$tmp = tmp1$tmp %>% factor(levels = unique(tmp$tmp))
tmp1 = aggrExpr(tmp1,features = cd4marker1,aggrby = tmp1$tmp,returnmatrix = T)

pdf('stp3_tme/plot/tc_analysis/cd4marker_heat_stage-area.pdf',width = 5,height = 5)
tmp1 %>% t %>% scale() %>% t %>% .[c('CXCL13','CD200','FOXP3','CTLA4','IL2RA'),] %>% 
  ComplexHeatmap::pheatmap(color=rev(brewer.pal(9,'RdBu')),cluster_cols = F,cluster_rows = F,
                           gaps_col = 4,breaks = seq(-2,2,length=100),
                           annotation_col = anncol,annotation_colors = anncolor,border_color = NA)
dev.off()

tmp1 = cd8
tmp1$tmp = paste0(tmp1$epiarea,'-',tmp1$stage1)
tmp1$tmp = tmp1$tmp %>% factor(levels = unique(tmp$tmp))
tmp1 = aggrExpr(tmp1,features = cd8marker1,aggrby = tmp1$tmp,returnmatrix = T)

pdf('stp3_tme/plot/tc_analysis/cd8marker_heat_stage-area.pdf',width = 5,height = 5)
tmp1 %>% t %>% scale() %>% t %>% .[c('GZMB','GNLY','GZMA','NKG7','GZMH','CTLA4','CXCL13','PDCD1','HAVCR2'),] %>% 
  ComplexHeatmap::pheatmap(color=rev(brewer.pal(9,'RdBu')),cluster_cols = F,cluster_rows = F,
                           gaps_col = 4,breaks = seq(-2,2,length=100),
                           annotation_col = anncol,annotation_colors = anncolor,border_color = NA)
dev.off()




##### spatial plot #####
{
  objlist = setNames(list('xe1_1','xe1_2','xe1_3','xe2_1','xe2_2','xe2_3','xe3_1','xe3_2',
                          'xe4_1','xe4_2','xe5_1','xe5_2_split','xe5_3_split','xe6_1','xe6_2','xe6_3'),
                     c('T1S1','T1S2','T1S3','T2S1','T2S2','T2S3','T3S1','T3S2',
                       'T4S1','T4S2','T5S1','T5S2','T5S3','T6S1','T6S2','T6S3'))
  lapply(names(objlist), function(x){
    coltmp = c(tccol,'Epi'='#dddddd','Others'='#888888')
    dt = get(objlist[[x]])
    
    dt$celltype = 'Others'
    dt$celltype[dt$ct3=='Epi'] = 'Epi'
    dt = renewMeta(dt, srttc@meta.data[,c('cellID','celltype')])
    
    flip=F;reversex=F;reversey=F
    g1 = c('T1S1','T1S2','T1S3','T3S1','T3S2','T4S1')
    g2 = c('T2S1','T2S2','T2S3')
    g3 = c('T4S2')
    g4 = c('T5S1','T5S2','T5S3')
    g5 = c('T6S1')
    g6 = c('T6S2','T6S3')
    if(x%in%g1){flip=T;reversex=T}
    if(x%in%g2){flip=T}
    if(x%in%g3){flip=T;reversey=T}
    if(x%in%g5){reversey=T}
    if(x%in%g6){reversex=T}
    
    fname = paste0('stp3_tme/plot/tc_spatial/',x,'.png')
    p = xedimplot(dt, groupby = 'celltype',flip = flip,ptlevels = c(names(tccol),'Epi','Others'),
                  reversex = reversex,reversey = reversey,color = coltmp,bgcol = 'black')
    ggsave(fname,p,width=40,height=20)
  })
}




##### proportion #####
srttc@meta.data %>% filter(epiarea=='Infiltrated') %>% filter(cellType_merge!='T-others') %>% propplot3()

cd4$cellType_merge = cd4$cellType_merge %>% factor(levels = rev(c('CD4T-reg','CD4T-fh','T-naive')))
cd8$cellType_merge = cd8$cellType_merge %>% factor(levels = rev(c('CD8T-ex','CD8T-eff','T-naive')))
propplot2(cd4@meta.data,color = primcol3,groupby = 'stage1',breaks = seq(0,1,0.25))+NoLegend()->p1
propplot2(cd8@meta.data,color = primcol3,groupby = 'stage1')+NoLegend()+
  theme(axis.text.y = element_blank())+ylab('')->p2
p1+p2
ggsave('stp3_tme/plot/tc_analysis/prop_tc.pdf',width = 7,height = 4)

p1=cd4@meta.data %>% filter(epiarea=='Infiltrated') %>% 
  propplot2(color = primcol3,groupby = 'stage1')+NoLegend()+ylab('')
p2=cd4@meta.data %>% filter(epiarea=='Out') %>% 
  propplot2(color = primcol3,groupby = 'stage1')+NoLegend()+ylab('')
ggsave('stp3_tme/plot/tc_analysis/prop_cd4-area.pdf',p1/p2,width = 4,height = 8)

p1=cd8@meta.data %>% filter(epiarea=='Infiltrated') %>% 
  propplot2(color = primcol3,groupby = 'stage1')+NoLegend()+ylab('')
p2=cd8@meta.data %>% filter(epiarea=='Out') %>%  
  propplot2(color = primcol3,groupby = 'stage1')+NoLegend()+ylab('')
ggsave('stp3_tme/plot/tc_analysis/prop_cd8-area.pdf',p1/p2,width = 4,height = 8)


propcd4 = getprop(cd4$cellType_merge,cd4@meta.data[,c('stage1','sp_stg')]) %>% 
  gather('celltype','prop',all_of(unique(cd4$cellType_merge)))
propcd8 = getprop(cd8$cellType_merge,cd8@meta.data[,c('stage1','sp_stg')]) %>% 
  gather('celltype','prop',all_of(unique(cd8$cellType_merge)))
mediancd4 = aggregate(list(dt=propcd4$prop),by=as.list(propcd4[,c('group_1','celltype'),drop=F]),median)
mediancd8 = aggregate(list(dt=propcd8$prop),by=as.list(propcd8[,c('group_1','celltype'),drop=F]),median)

propcd4 %>%
  ggplot(aes(group_1,prop))+
  theme1+
  theme(axis.ticks.x = element_blank())+rotate_x_text(angle = 30)+
  NoLegend()+
  ylab('Proportion')+xlab('')+
  scale_y_continuous(labels = c('0','0.5','1.0'),breaks = seq(0,1,0.5),limits = c(0,1.4))+
  scale_fill_manual(values = col_stage)+scale_color_manual(values = col_stage)+
  stat_compare_means(comparisons = makecomparison(4))+
  facet_wrap(~celltype)+
  geom_boxplot(aes(fill=group_1),alpha=0.5,outlier.color = NA)+
  geom_quasirandom(aes(color=group_1))
ggsave('stp3_tme/plot/tc_analysis/prop_cd4.pdf',width = 6,height = 3.5)
aggregate(list(dt=propcd4$prop),by=as.list(propcd4[,c('group_1','celltype'),drop=F]),median)

propcd8 %>%
  ggplot(aes(group_1,prop))+
  theme1+
  theme(axis.ticks.x = element_blank())+rotate_x_text(angle = 30)+
  NoLegend()+
  ylab('Proportion')+xlab('')+
  scale_y_continuous(labels = c('0','0.5','1.0'),breaks = seq(0,1,0.5),limits = c(0,1.4))+
  scale_fill_manual(values = col_stage)+scale_color_manual(values = col_stage)+
  stat_compare_means(comparisons = makecomparison(4))+
  facet_wrap(~celltype)+
  geom_boxplot(aes(fill=group_1),alpha=0.5,outlier.color = NA)+
  geom_quasirandom(aes(color=group_1))
ggsave('stp3_tme/plot/tc_analysis/prop_cd8.pdf',width = 6,height = 3.5)
aggregate(list(dt=propcd8$prop),by=as.list(propcd8[,c('group_1','celltype'),drop=F]),median)





##### proportion split by area #####
propbar = function(meta){
  meta$celltype = meta$celltype %>% droplevels()
  propnum = table(meta[,c('celltype','stage1','epiarea','sp_stg')]) %>% data.frame()
  
  totalnum = aggregate(list(total=propnum$Freq),by=as.list(propnum[,c('stage1','epiarea','sp_stg')]),sum) %>% 
    filter(total>0)
  propnum = propnum %>% inner_join(totalnum)
  propnum$prop = propnum$Freq/propnum$total
  propnum$prop[propnum$epiarea=='Out'] = -propnum$prop[propnum$epiarea=='Out']
  
  propnum_mean = aggregate(list(meanvalue=propnum$prop),by=as.list(propnum[,c('stage1','epiarea','celltype')]),median)
  
  propnum %>%
    ggplot(aes(stage1,prop))+
    theme1+theme(axis.line.x = element_blank(),axis.ticks.x = element_blank())+
    NoLegend()+
    facet_wrap(~celltype)+
    scale_fill_manual(values = tccol)+scale_color_manual(values = tccol)+
    # scale_y_continuous(expand = c(0,0.3))+
    # stat_mean(geom = 'bar',aes(fill=celltype),alpha=0.5)+
    xlab('')+ylab('Cell proportion')+
    geom_bar(data=propnum_mean,aes(stage1,meanvalue,fill=celltype),color='black',stat='identity',position = position_dodge(),alpha=0.5)+
    # geom_point(aes(fill=celltype),color='black',size=1.5,shape=21,stroke=0.8,position=position_jitter())+
    geom_point(aes(color=celltype),size=1.5,position=position_jitter())+
    stat_summary(data=propnum[propnum$epiarea=='Infiltrated',],geom = 'errorbar',
                 fun.max = function(x){quantile(x,0.75)},
                 fun.min = function(x){quantile(x,0.25)},
                 width=0.5,linewidth=0.6,
                 aes(group=celltype),position = position_dodge(width = 0.9))+
    stat_compare_means(data=propnum[propnum$epiarea=='Infiltrated',],
                       comparisons = makecomparison(4),tip.length = 0,label.y = c(1,1.2,1.4))+
    stat_summary(data=propnum[propnum$epiarea=='Out',],geom = 'errorbar',
                 fun.max = function(x){quantile(x,0.75)},
                 fun.min = function(x){quantile(x,0.25)},
                 width=0.5,linewidth=0.6,
                 aes(group=celltype),position = position_dodge(width = 0.9))+
    stat_compare_means(data=propnum[propnum$epiarea=='Out',],
                       comparisons = makecomparison(4),tip.length = 0,label.y = c(-1.2,-1.4,-1.6),vjust = 1.6)+
    geom_hline(yintercept = 0)
}

propbar(cd4@meta.data)
ggsave('stp3_tme/plot/tc_analysis/prop_cd4_splitbyarea.pdf',width = 8,height = 6)
propbar(cd8@meta.data)
ggsave('stp3_tme/plot/tc_analysis/prop_cd8_splitbyarea.pdf',width = 8,height = 6)

prop1 = getprop(cd4$cellType_merge,cd4@meta.data[,c('stage1','sp_stg','epiarea')]) %>% 
  gather('celltype','prop',all_of(unique(cd4$cellType_merge)))
prop1$group_1 = prop1$group_1 %>% factor(levels = c('NOR','LGIN','HGIN','ESCC'))
prop1 %>% 
  ggplot(aes(group_1,prop))+
  facet_wrap(~group_3+celltype)+
  stat_compare_means(comparisons = makecomparison(4))+
  geom_boxplot()
ggsave('stp3_tme/plot/tc_analysis/prop_cd4_splitbyarea_pval.pdf',width = 8,height = 8)
prop2 = getprop(cd8$cellType_merge,cd8@meta.data[,c('stage1','sp_stg','epiarea')]) %>% 
  gather('celltype','prop',all_of(unique(cd8$cellType_merge)))
prop2$group_1 = prop2$group_1 %>% factor(levels = c('NOR','LGIN','HGIN','ESCC'))
prop2 %>% 
  ggplot(aes(group_1,prop))+
  facet_wrap(~group_3+celltype)+
  stat_compare_means(comparisons = makecomparison(4))+
  geom_boxplot()
ggsave('stp3_tme/plot/tc_analysis/prop_cd8_splitbyarea_pval.pdf',width = 8,height = 8)





##### epi area #####
{## t1s2
  poly12 = read_rds('stp3_tme/data/epiregion/T1S2_poly.rds')
  
  polyall = lapply(1:length(poly12), function(x){
    dt = data.frame(poly12[[x]])
    dt$grp = x
    return(dt)
  }) %>% Reduce(rbind,.)
  colnames(polyall) = c('x','y','grp')
  
  c12 = xe1_2@images$fov@boundaries$centroids@coords
  rownames(c12) = colnames(xe1_2)
  c12 = c12[rownames(c12)%in%srttc$cellID,]
  c12 = c12[,c(2,1)]
  colnames(c12) = c('x','y')
  
  system.time({
    tmp = lapply(poly12, function(x){
      rst = isinner2(c12, x)
    }) %>% Reduce('|',.)
  })
  
  ggplot(polyall, aes(x,y))+
    geom_polygon(aes(group=grp),alpha=0.25)+
    geom_point(data = data.frame(c12[tmp,]))
  xedimplot(xe1_2, groupby = 'ct3',flip = T,reversex = T,color = primcol2,bgcol = 'black')+
    geom_polygon(data = polyall,aes(group=grp),alpha=0.5,fill='lightblue',color='white',linewidth=0.5)
  ggsave('stp3_tme/xe12_epiregion.png',width = 20,height = 10)
  
  tc12 = subset(xe1_2, cells = xe1_2$cellID[xe1_2$ct3=='T'])
  tc12inepi = rownames(c12)[tmp]
  tc12$epiarea = 'out'
  tc12$epiarea[tc12$cellID%in%tc12inepi] = 'infiltrated'
  xedimplot(tc12,groupby = 'epiarea',flip = T,reversex = T,bgcol = 'black')+
    geom_polygon(data = polyall,aes(group=grp),alpha=0.5,fill='lightblue',color='white',linewidth=0.5)
  ggsave('stp3_tme/xe12_tc_epiregion.png',width = 20,height = 10)
  
  tmp1 = subset(tc12, cells = tc12$cellID[tc12$cellType_merge!='T-others'])
  tprop = getprop(tmp1$epiarea,tmp1@meta.data[,c('stage')])
  propplot(tmp1$epiarea,tmp1$stage)
  subprop_inf = getprop(tmp1$cellType_merge,tmp1@meta.data[,c('stage','epiarea')])
}

{## get infiltrated t
  objlist = setNames(list('xe1_1','xe1_2','xe1_3','xe2_1','xe2_2','xe2_3','xe3_1','xe3_2',
                          'xe4_1','xe4_2','xe5_1','xe5_2_split','xe5_3_split','xe6_1','xe6_2','xe6_3'),
                     c('T1S1','T1S2','T1S3','T2S1','T2S2','T2S3','T3S1','T3S2',
                       'T4S1','T4S2','T5S1','T5S2','T5S3','T6S1','T6S2','T6S3'))
  allepitc=lapply(names(objlist), function(x){
    dt = get(objlist[[x]])
    dtt = subset(dt,cells=intersect(dt$cellID,srttc$cellID))
    
    epipoly = paste0('stp3_tme/data/epiregion/',x,'_poly.rds')
    epipoly = read_rds(epipoly)
    epipolyall = lapply(1:length(epipoly), function(x){
      dt = data.frame(epipoly[[x]])
      dt$grp = x
      return(dt)
    }) %>% Reduce(rbind,.)
    colnames(epipolyall) = c('x','y','grp')
    
    flip=F;reversex=F;reversey=F
    g1 = c('T1S1','T1S2','T1S3','T3S1','T3S2','T4S1')
    g2 = c('T2S1','T2S2','T2S3')
    g3 = c('T4S2')
    g4 = c('T5S1','T5S2','T5S3')
    g5 = c('T6S1')
    g6 = c('T6S2','T6S3')
    if(x%in%g1){flip=T;reversex=T}
    if(x%in%g2){flip=T}
    if(x%in%g3){flip=T;reversey=T}
    if(x%in%g5){reversey=T}
    if(x%in%g6){reversex=T}
    
    crd = dtt@images$fov@boundaries$centroids@coords
    rownames(crd) = colnames(dtt)
    if(flip){
      crd = crd[,c(2,1)]
      colnames(crd) = c('x','y')
    }
    inepi = lapply(epipoly, function(x){
      rst = isinner2(crd, x)
    }) %>% Reduce('|',.)
    
    epit = colnames(dtt)[inepi]
    
    dt$tmp = 'Others'
    dt$tmp[dt$cellID%in%dtt$cellID] = 'Out'
    dt$tmp[dt$cellID%in%epit] = 'Infiltrated'
    
    fname = paste0('stp3_tme/plot/tc_infiltrated/',x,'.png')
    p = xedimplot(dt, groupby = 'tmp',flip = flip,reversex = reversex,reversey = reversey,bgcol = 'black',ptlevels = c('Infiltrated','Out','Others'),
                  color = c('#888888',c(hue_pal()(2))))+geom_polygon(data=epipolyall,aes(group=grp),fill='lightblue',alpha=0.3)
    ggsave(fname,p,width=30,height=12)
    
    return(epit)
  }) %>% Reduce(append,.)
  
  srttc$epiarea = 'Out'
  srttc$epiarea[srttc$cellID%in%allepitc] = 'Infiltrated'
  srttc$stage = srttc$stage1
  srttc$sampleID = srttc$sp_stg
  srttc$celltype = srttc$cellType_merge
  srttc$cellType_merge = srttc$cellType_merge %>% mapvalues(c('CD4T-n/m','CD8T-n/m'),c('T-naive','T-naive'))
  
  saveRDS(srttc@meta.data, 'stp3_tme/data/tc/all_meta.rds')
}

{## get infiltrated plot 
  objlist = setNames(list('xe1_1','xe1_2','xe1_3','xe2_1','xe2_2','xe2_3','xe3_1','xe3_2',
                          'xe4_1','xe4_2','xe5_1','xe5_2_split','xe5_3_split','xe6_1','xe6_2','xe6_3'),
                     c('T1S1','T1S2','T1S3','T2S1','T2S2','T2S3','T3S1','T3S2',
                       'T4S1','T4S2','T5S1','T5S2','T5S3','T6S1','T6S2','T6S3'))
  lapply(names(objlist), function(x){
    dt = get(objlist[[x]])
    
    polye = epipoly[[x]]
    epipolyall = lapply(1:length(polye), function(y){
      dt = data.frame(polye[[y]])
      dt$grp = y
      return(dt)
    }) %>% Reduce(rbind,.)
    colnames(epipolyall) = c('x','y','grp')
    
    flip=F;reversex=F;reversey=F
    g1 = c('T1S1','T1S2','T1S3','T3S1','T3S2','T4S1')
    g2 = c('T2S1','T2S2','T2S3')
    g3 = c('T4S2')
    g4 = c('T5S1','T5S2','T5S3')
    g5 = c('T6S1')
    g6 = c('T6S2','T6S3')
    if(x%in%g1){flip=T;reversex=T}
    if(x%in%g2){flip=T}
    if(x%in%g3){flip=T;reversey=T}
    if(x%in%g5){reversey=T}
    if(x%in%g6){reversex=T}
    
    
    dt$epiarea = 'Others'
    dt = renewMeta(dt,srttc@meta.data[,c('cellID','epiarea')])
    
    fname = paste0('stp3_tme/plot/tc_infiltrated/',x,'.png')
    p = xedimplot(dt, groupby = 'epiarea',flip = flip,reversex = reversex,reversey = reversey,bgcol = 'black',ptlevels = c('Infiltrated','Out','Others'),
                  color = areacol)+geom_polygon(data=epipolyall,aes(group=grp),color='white',linetype='dashed',fill=NA)
    ggsave(fname,p,width=30,height=12)
  }) %>% Reduce(append,.)
}






##### epi area subtype #####
{
  objlist = setNames(list('xe1_1','xe1_2','xe1_3','xe2_1','xe2_2','xe2_3','xe3_1','xe3_2',
                          'xe4_1','xe4_2','xe5_1','xe5_2_split','xe5_3_split','xe6_1','xe6_2','xe6_3'),
                     c('T1S1','T1S2','T1S3','T2S1','T2S2','T2S3','T3S1','T3S2',
                       'T4S1','T4S2','T5S1','T5S2','T5S3','T6S1','T6S2','T6S3'))
  allepitc=lapply(names(objlist), function(x){
    dt = get(objlist[[x]])
    dt$cellType_merge = 'Others'
    dt = renewMeta(dt,srttc@meta.data)
    dt$cellType_merge = dt$cellType_merge %>% mapvalues('T-others','Others')
    
    epipoly = paste0('stp3_tme/data/epiregion/',x,'_poly.rds')
    epipoly = read_rds(epipoly)
    epipolyall = lapply(1:length(epipoly), function(x){
      dt = data.frame(epipoly[[x]])
      dt$grp = x
      return(dt)
    }) %>% Reduce(rbind,.)
    colnames(epipolyall) = c('x','y','grp')
    
    flip=F;reversex=F;reversey=F
    g1 = c('T1S1','T1S2','T1S3','T3S1','T3S2','T4S1')
    g2 = c('T2S1','T2S2','T2S3')
    g3 = c('T4S2')
    g4 = c('T5S1','T5S2','T5S3')
    g5 = c('T6S1')
    g6 = c('T6S2','T6S3')
    if(x%in%g1){flip=T;reversex=T}
    if(x%in%g2){flip=T}
    if(x%in%g3){flip=T;reversey=T}
    if(x%in%g5){reversey=T}
    if(x%in%g6){reversex=T}
    
    tmpcol = append(c('Others'='#888888'),primcol3)
    
    fname = paste0('stp3_tme/plot/tc_subtype/',x,'.png')
    p = xedimplot(dt, groupby = 'cellType_merge',flip = flip,reversex = reversex,reversey = reversey,
                  bgcol = 'black',ptlevels = c(rev(levels(srttc$cellType_merge)),'Others'),color = tmpcol)+
      geom_polygon(data=epipolyall,aes(group=grp),color='white',fill=NA,alpha=0.2,linetype='dashed')
    ggsave(fname,p,width=20,height=10)
  })
}




##### infiltrated density #####
{
  objlist = setNames(list('xe1_1','xe1_2','xe1_3','xe2_1','xe2_2','xe2_3','xe3_1','xe3_2',
                          'xe4_1','xe4_2','xe5_1','xe5_2_split','xe5_3_split','xe6_1','xe6_2','xe6_3'),
                     c('T1S1','T1S2','T1S3','T2S1','T2S2','T2S3','T3S1','T3S2',
                       'T4S1','T4S2','T5S1','T5S2','T5S3','T6S1','T6S2','T6S3'))
  tdens=lapply(names(objlist), function(x){
    dt = get(objlist[[x]])
    
    inft = srttc$cellID[srttc$epiarea=='Infiltrated']
    dt$tmp = 'Others'
    dt$tmp[dt$cellID%in%srttc$cellID] = 'Out'
    dt$tmp[dt$cellID%in%inft] = 'Infiltrated'
    
    crd_all = dt@images$fov@boundaries$centroids@coords
    rownames(crd_all) = colnames(dt)
    
    tdens_sp = lapply(unique(dt$sp_stg),function(y){
      dt1 = dt@meta.data %>% filter(sp_stg==y)
      crd = crd_all[dt1$cellID,]
      
      areaall = crd %>% data.frame %>% dplyr::slice(chull(crd))
      areaall_ = areaall %>% as.matrix %>% splancs::areapl()
      areaall_ = areaall_/1e6
      
      areaepi = crd[dt1$cellID[dt1$ct3=='Epi'],] 
      areaepi = areaepi %>% data.frame %>% dplyr::slice(chull(areaepi))
      areaepi_ = areaepi %>% as.matrix %>% splancs::areapl()
      areaepi_ = areaepi_/1e6
      
      areames_ = areaall_-areaepi_
      
      infnum = sum(dt1$tmp=='Infiltrated')
      outnum = sum(dt1$tmp=='Out')
      infdens = infnum/areaepi_
      outdens = outnum/areames_
      rst = data.frame(sp_stg=y,inf=infdens,out=outdens)
      return(rst)
    }) %>% Reduce(rbind,.)
    
    return(tdens_sp)
  }) %>% Reduce(rbind,.)
  
  stg = unique(allmeta[,c('sp_stg','stage1')])
  tdens = tdens %>% left_join(stg)
  tdens$out[is.infinite(tdens$out)] = 0
  tdens = tdens %>% gather('epiarea','normnum',inf,out)
  saveRDS(tdens,'stp3_tme/data/tc/tc_density.rds')
}


tdens = read_rds('stp3_tme/data/tc/tc_density.rds')
tdens %>% 
  filter(normnum<1000) %>% 
  ggplot(aes(stage1,normnum))+
  theme1+
  scale_fill_manual(values = col_stage)+scale_color_manual(values = col_stage)+
  NoLegend()+
  geom_jitter(aes(color=stage1))+
  geom_boxplot(aes(fill=stage1),alpha=0.5,outlier.color = NA)+
  stat_compare_means(comparisons = append(makecomparison(4),list(c(2,3),c(2,4))),label = 'p.signif')+
  # ylim(c(0,1000))+
  facet_wrap(~epiarea)
ggsave('stp3_tme/plot/tc_analysis/Tdensity_box.pdf',width = 6,height = 5)
tdens %>% 
  filter(normnum<1000) %>% 
  ggplot(aes(epiarea,normnum))+
  theme1+
  scale_fill_manual(values = c("#80b1d3","#fccde5"))+scale_color_manual(values = c("#80b1d3","#fccde5"))+
  NoLegend()+
  geom_jitter(aes(color=epiarea))+
  geom_boxplot(aes(fill=epiarea),alpha=0.5,outlier.color = NA)+
  stat_compare_means(comparisons = list(c(1,2)),label = 'p.signif')+
  # ylim(c(0,1000))+
  facet_wrap(~stage1,nrow = 1)
ggsave('stp3_tme/plot/tc_analysis/Tdensity_areacompare_box.pdf',width = 6,height = 4)




##### infiltrated num normalize to epi num #####
{
  objlist = setNames(list('xe1_1','xe1_2','xe1_3','xe2_1','xe2_2','xe2_3','xe3_1','xe3_2',
                          'xe4_1','xe4_2','xe5_1','xe5_2_split','xe5_3_split','xe6_1','xe6_2','xe6_3'),
                     c('T1S1','T1S2','T1S3','T2S1','T2S2','T2S3','T3S1','T3S2',
                       'T4S1','T4S2','T5S1','T5S2','T5S3','T6S1','T6S2','T6S3'))
  tnorm=lapply(names(objlist), function(x){
    dt = get(objlist[[x]])
    dt$epiarea = 'Out'
    epimeta = srttc@meta.data[,c('cellID','epiarea')]
    dt = renewMeta(dt,epimeta)
    
    tnorm_sp = lapply(unique(dt$sp_stg),function(y){
      dt1 = dt@meta.data %>% filter(sp_stg==y)
      
      tnum = sum(dt1$epiarea=='Infiltrated')
      epinum = sum(dt1$ct3=='Epi')
      
      tnorm = tnum/epinum
      rst = data.frame(sp_stg=y,tnorm=tnorm)
      return(rst)
    }) %>% Reduce(rbind,.)
    
    return(tnorm_sp)
  }) %>% Reduce(rbind,.)
  
  stg = unique(allmeta[,c('sp_stg','stage1')])
  tnorm = tnorm %>% left_join(stg)
  saveRDS(tnorm,'stp3_tme/data/tc/tc_normtoEpi.rds')
}

tnorm %>% 
  ggplot(aes(stage1, tnorm))+
  theme1+
  scale_fill_manual(values = col_stage)+scale_color_manual(values = col_stage)+
  NoLegend()+
  geom_boxplot(aes(fill=stage1),alpha=0.5,outlier.color = NA)+
  geom_jitter(aes(color=stage1))+
  stat_compare_means(comparisons = makecomparison(4))

infprop = srttc@meta.data %>% filter(epiarea=='Infiltrated')
infprop = getprop(infprop$cellType_merge, infprop[,c('sp_stg','stage1')])
colnames(infprop)[1:2] = c('sp_stg','stage1')
infprop = left_join(infprop, tnorm)
infprop = infprop %>% gather('celltype','cellprop',all_of(unique(srttc$cellType_merge)))

infprop %>% 
  ggplot(aes(stage1, tnorm*cellprop))+
  geom_boxplot()+
  stat_compare_means(comparisons = makecomparison(4))+
  facet_wrap(~celltype,scales = 'free')




##### inf analysis #####
srttc@meta.data %>% filter(cellType_merge!='T-others') %>% propplot3(cellcolumn='epiarea')
ggsave('stp3_tme/plot/tc_analysis/prop_inf_box.pdf',width = 4,height = 8)

propplot2(srttc@meta.data,col_cells,'epiarea','stage1')
srttc@meta.data %>% filter(epiarea=='Infiltrated') %>% filter(cellType_merge!='T-others') %>% propplot3(plotrow = 1)
ggsave('stp3_tme/plot/tc_analysis/prop_inf-subtype_box.pdf',width = 8,height = 8)

cd4@meta.data %>% filter(epiarea=='Infiltrated') %>% propplot3(plotrow = 1)
cd8@meta.data %>% filter(epiarea=='Infiltrated') %>% propplot3(plotrow = 1)

tmp = cd8@meta.data %>% filter(epiarea=='Infiltrated')
tmp1 = getprop(tmp$cellType_merge,tmp[,c('stage1','sp_stg')])





##### t cell density #####
{## calculate by grid
  gridid = read_rds('stp1_summary/data/allsample_gridid.rds')
  tcgrid = gridid %>% dplyr::select(-cellType_merge) %>% inner_join(srttc@meta.data[,c('cellID','cellType_merge','epiarea')])
  
  tcnum = table(tcgrid$epiarea,tcgrid$cellType_merge,tcgrid$gridid,tcgrid$stage1,tcgrid$sp_stg) %>% 
    data.frame() %>% `colnames<-`(c('epiarea','celltype','gridid','stage','sp_stg','cellnum')) %>% 
    filter(cellnum>0)
  tcnum = table(tcgrid$cellType_merge,tcgrid$gridid,tcgrid$stage1,tcgrid$sp_stg) %>% 
    data.frame() %>% `colnames<-`(c('celltype','gridid','stage','sp_stg','cellnum')) %>% 
    filter(cellnum>0)

  tcmeannum = aggregate(list(meannum=tcnum$cellnum),
                          by=as.list(tcnum[,c('celltype','stage','sp_stg')]),
                          # by=as.list(tcnum[,c('celltype','stage','sp_stg','epiarea')]),
                          mean)
  tmp = aggregate(list(meannum=tcmeannum$meannum),
                  by=as.list(tcmeannum[,c('celltype','stage')]),
                  mean)
  
  tcmeannum %>% 
    # filter(epiarea=='Infiltrated') %>% 
    ggplot(aes(stage,meannum))+
    theme1+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
    xlab('')+
    NoLegend()+
    stat_compare_means(comparisons = makecomparison(4))+
    facet_wrap(~celltype,scales = 'free')+
    geom_jitter(size=2,aes(color=celltype))+
    geom_boxplot(aes(fill=celltype),outlier.color = NA,alpha=0.5)
  # ggsave('stp3_tme/plot/tc_analysis/density-grid_tcsubtype.pdf',width = 16,height = 3.2)
}

{## density per cell area
  dens_all = read_rds('stp1_summary/data/allsample_density_primtype.rds')
  dens_all = allmeta %>% left_join(dens_all[,c('cellID','dens')])
  dens_t = dens_all %>% filter(ct3=='T') %>% left_join(srttc@meta.data[,c('cellID','epiarea')])
  meandens_t = aggregate(list(dens=dens_t$dens),by=as.list(dens_t[,c('sp_stg','stage1','epiarea')]),mean)
  
  meandens_t %>% 
    ggplot(aes(stage1,dens))+
    theme1+theme(axis.ticks.x = element_blank())+
    xlab('')+ylab('Cell density')+
    NoLegend()+rotate_x_text(angle = 30)+
    stat_compare_means(comparisons = makecomparison(4))+
    scale_fill_manual(values = areacol)+scale_color_manual(values = areacol)+
    scale_y_continuous(expand = c(0.1,0.3))+
    facet_wrap(~epiarea,scales = 'free',nrow = 2)+
    geom_jitter(size=2,aes(color=epiarea))+
    geom_boxplot(aes(fill=epiarea),outlier.color = NA,alpha=0.5)
  ggsave('stp3_tme/plot/tc_analysis/tc_density_area.pdf',width = 3,height = 7)
}

{## density per cell subtype
  objlist = setNames(list('xe1_1','xe1_2','xe1_3','xe2_1','xe2_2','xe2_3','xe3_1','xe3_2',
                          'xe4_1','xe4_2','xe5_1','xe5_2_split','xe5_3_split','xe6_1','xe6_2','xe6_3'),
                     c('T1S1','T1S2','T1S3','T2S1','T2S2','T2S3','T3S1','T3S2',
                       'T4S1','T4S2','T5S1','T5S2','T5S3','T6S1','T6S2','T6S3'))
  job({
    dens_tc = lapply(names(objlist), function(x){
      message(paste0('running sample ',x))
      crd = crds[[x]]
      crd = crd[,c('x','y')]
      meta = srttc@meta.data[srttc$cellID%in%rownames(crd),]
      
      dens_ct = lapply(unique(meta$celltype),function(y){
        message(paste0('running ct ',y))
        meta1 = meta[meta$celltype==y,]
        crd1 = crd[rownames(crd)%in%meta1$cellID,]
        
        cdens = calcdens(crd1, 100)
        rst = data.frame(cellID=rownames(crd1),dens = cdens, celltype = y)
        return(rst)
      }) %>% Reduce(rbind,.)
      
      return(dens_ct)
    }) %>% Reduce(rbind,.)
    saveRDS(dens_tc,'stp3_tme/data/tc/tc_subtypeDensity.rds')
  },import = 'auto')
  
  dens_tc = srttc@meta.data %>% left_join(dens_tc[,c('cellID','dens')])
  
  ## subtype
  meandens_tc = aggregate(list(dens=dens_tc$dens),by=as.list(dens_tc[,c('sp_stg','stage1','celltype')]),mean)
  meandens_tc %>% 
    filter(celltype%in%names(tccol)) %>% 
    ggplot(aes(stage1,dens))+
    theme1+theme(axis.ticks.x = element_blank())+
    xlab('')+ylab('Cell density')+
    NoLegend()+rotate_x_text(angle = 30)+
    stat_compare_means(comparisons = makecomparison(4))+
    scale_fill_manual(values = tccol)+scale_color_manual(values = tccol)+
    scale_y_continuous(expand = c(0.1,0.3))+
    facet_wrap(~celltype,scales = 'free',nrow = 1)+
    geom_jitter(size=2,aes(color=celltype))+
    geom_boxplot(aes(fill=celltype),outlier.color = NA,alpha=0.5)
  ggsave('stp3_tme/plot/tc_analysis/tc_density_subtype.pdf',width = 13,height = 3.5)
  
  ## subtype & area
  meandens_tc_area = aggregate(list(dens=dens_tc$dens),by=as.list(dens_tc[,c('sp_stg','stage1','celltype','epiarea')]),mean)
  meandens_tc_area %>% 
    filter(celltype%in%names(tccol)) %>% 
    ggplot(aes(stage1,dens))+
    theme1+theme(axis.ticks.x = element_blank())+
    xlab('')+ylab('Cell density')+
    NoLegend()+rotate_x_text(angle = 30)+
    stat_compare_means(comparisons = makecomparison(4))+
    scale_fill_manual(values = tccol)+scale_color_manual(values = tccol)+
    scale_y_continuous(expand = c(0.1,0.3))+
    facet_wrap(~epiarea+celltype,scales = 'free',nrow = 2)+
    geom_jitter(size=2,aes(color=celltype))+
    geom_boxplot(aes(fill=celltype),outlier.color = NA,alpha=0.5)
  ggsave('stp3_tme/plot/tc_analysis/tc_density_subtype-area.pdf',width = 13,height = 7)
}




##### subtype function score #####
gnlist = list(cytogene = c('GNLY', 'NKG7', 'GZMB', 'GZMA', 'GZMH', 'KLRB1', 'GZMK','TYROBP'),  
              exhgene = c('PDCD1','HAVCR2','CTLA4','TIGIT','LAG3','CXCL13'), 
              reggene = c('FOXP3','TIGIT','CD274','IL2RA'),
              naivegene = c("SELL","CCR7","LEF1","TCF7"),
              actgene = c('FCGR3A','STAT5A'),
              prolifegene = c('MKI67','TOP2A')
)

tscore = getscore(srttc, gnlist)
tscore = tscore %>% gather('func','score',all_of(names(gnlist)))
tscore_aggr = aggregate(list(score=tscore$score),by = as.list(tscore[,c('sampleID','stage1','celltype','epiarea','func')]),mean)

## groupby celltype 
tscore %>% 
  ggplot(aes(celltype,score))+
  theme1+RotatedAxis()+
  facet_wrap(~func)+
  geom_boxplot(aes(fill=celltype))



## groupby celltype + stage
tscore %>% 
  ggplot(aes(stage1,score))+
  theme1+RotatedAxis()+
  facet_grid(func~celltype)+
  geom_boxplot(aes(fill=celltype))



## groupby infiltrated 
tscore_aggr %>% 
  filter(func%in%c('naivegene','actgene','exhgene','cytogene')) %>% 
  filter(celltype%in%c('CD8T-n/m','CD8T-eff','CD8T-ex')) %>% 
  ggplot(aes(epiarea, score))+
  theme1+RotatedAxis()+scale_y_continuous(expand = expansion(c(0.05,0.1)))+
  scale_fill_manual(values = tccol)+scale_color_manual(values = tccol)+
  stat_compare_means(comparisons = makecomparison(2))+
  facet_grid(celltype~func)+
  # geom_violin(aes(fill=celltype),alpha=0.5,scale = 'width')+
  ggbeeswarm::geom_quasirandom(aes(color=celltype),shape=16)+
  stat_summary(geom='errorbar',fun.min = function(x) quantile(x,0.25),fun.max = function(x) quantile(x,0.75),
               linewidth=0.5,width=0.3,lineend='butt')+
  stat_summary(geom='errorbar',fun.min = median,fun.max = median,
               linewidth=1,width=0.5,lineend='butt')
ggsave('stp3_tme/plot/tc_analysis/funcscore_cd8.pdf',width = 12,height = 8)

tscore_aggr %>% 
  filter(func%in%c('naivegene','actgene','reggene')) %>% 
  filter(celltype%in%c('CD4T-n/m','CD4T-fh','CD4T-reg')) %>% 
  ggplot(aes(epiarea, score))+
  theme1+RotatedAxis()+scale_y_continuous(expand = expansion(c(0.05,0.1)))+
  scale_fill_manual(values = tccol)+scale_color_manual(values = tccol)+
  stat_compare_means(comparisons = makecomparison(2))+
  facet_grid(celltype~func)+
  # geom_violin(aes(fill=celltype),alpha=0.5,scale = 'width')+
  ggbeeswarm::geom_quasirandom(aes(color=celltype),shape=16)+
  stat_summary(geom='errorbar',fun.min = function(x) quantile(x,0.25),fun.max = function(x) quantile(x,0.75),
               linewidth=0.5,width=0.3,lineend='butt')+
  stat_summary(geom='errorbar',fun.min = median,fun.max = median,
               linewidth=1,width=0.5,lineend='butt')
ggsave('stp3_tme/plot/tc_analysis/funcscore_cd4.pdf',width = 9,height = 8)



## score heatmap
tscore_aggr1 = aggregate(list(score=tscore$score),by = as.list(tscore[,c('celltype','func')]),mean)
tscore_aggr1 = tscore_aggr1[grep('CD[48]',tscore_aggr1$celltype),]
tscore_aggr1$tctype = tscore_aggr1$celltype %>% as.character %>%  strsplit('-',fixed = T) %>% sapply('[',1)

tscore_aggr1 = tscore_aggr1 %>% spread(func,score)
rownames(tscore_aggr1) = tscore_aggr1$celltype

pheatmap(tscore_aggr1[,-1:-2],scale = 'column',annotation_row = tscore_aggr1[,2,drop=F],cluster_cols = F,cluster_rows = F)






