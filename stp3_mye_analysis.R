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

allmeta = read_rds('data/xenium_seurat_processed/allcell_meta.rds')

srtmye = read_rds('data/xenium_seurat_processed/ct_mye_srt.rds')
srtmye@meta.data = read_rds('data/xenium_seurat_processed/ct_mye_meta.rds')


myemarker = c('SFN','DCN','PLN','VWF','CD19','CD68','CD2', # pan marker
              'MS4A2','KIT','CPA3', # mast
              'GZMB','JCHAIN', # pDC
              'CD1C','CD1E','FCER1A', # APC
              'CCR7', # tDC
              'FCN1','S100A9','S100A8','FCGR3A', # Mono
              'IL1B', # Macro
              'LYVE1','C1QC','C1QB','C1QA','SPP1','FN1' # Macro
)

primcol2 = c('#1f78b4','#ffc089','#B5EFB5','#793c1b',
             '#6a3d9a','#333333','#ffff33','#df928e')
names(primcol2) = c('B','Endo','Fib','Gland','Mye','Mus','T','Epi')

primcol3 = c('#793c1b','#333333','#ffc089','#df928e',
             '#1f78b4','#6a3d9a','#df65b0','#1cbe4f','#e31a1c','#ffdcbd',
             '#BDCDFF','#ffff68','#325A9B','#AAF400','#ff9344')
names(primcol3) = c('Gland','Mus','Endo','Epi',
                    'B','Mye','NF-PI16','NF','CAF','T-others',
                    'T-naive','CD8T-eff','CD8T-ex','CD4T-fh','CD4T-reg')
areacol = c('#fccde5','#80b1d3')

tmpcol = c('#df928e','#888888','#1cbe4f','#238b45','#B5EFB5',
           '#df65b0','#ef3b2c','#a50f15')
names(tmpcol) = c('Epi','Others','CD8T-n/m','CD8T-eff','CD8T-ex',
                  'CD4T-fh','CD4T-reg','CD4T-n/m')

allcellcol = c('Basal'='#0339f8','Proliferation'='#ff9408','Differentiation'='#75bbfd',
               'Terminal'='#4682B4','HT_Prolif'='#ff9408','HT_Diff/Term'='#4682B4',
               'Invasive'='#de0c62','T cells'='#96f97b','B cells'='#028f1e',
               'Myeloid cells'='#7e1e9c','Fibroblast'='#696969','CAF'='#FFFFFF',
               'Endothelial cells'='#7b002c','Gland cells'='#c69f59','Myocytes'='#1b2431')

myecol = c(Mono='#7fc97f',Mac='#1f78b4','Mac-2'='#e31a1c','Mac-3'='#ffff33',DC='#beaed4',Mast='#b15928')


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




##### umap #####
{## subset 
  myesbst = sample(srtmye$cellID, dim(srtmye)[2]/6)
  myesbst = subset(srtmye, cells = myesbst)
  
  myesbst = SCTransform(myesbst, ncells=3000)
  myesbst = RunPCA(myesbst)
  myesbst = RunUMAP(myesbst, dims=1:10)
  
  mydimp(myesbst,groupby = 'cellType_merge',cols = brewer.pal(9,'Set1'))+NoLegend()
}


mydimp(srtmye,groupby = 'cellType_merge',cols = brewer.pal(9,'Set1'))+NoLegend()




##### dotp #####
myemarker2 = c('MS4A2','KIT','CPA3', # mast
               'JCHAIN','CD1C','FCER1A','CCR7',
               'FCN1','S100A9','S100A8','FCGR3A', # Mono
               'LYVE1','C1QC','C1QB','C1QA','IL1B','SPP1')
mydotp(srtmye,myemarker2,'cellType_merge',cols = rev(brewer.pal(9,'RdBu')))+coord_flip()+RotatedAxis()
ggsave('stp3_tme/plot/mye_analysis/mye_marker_dot.pdf',width = 6,height = 4)
DotPlot(srtmye,features = myemarker2,group.by = 'cellType_merge',cols = 'RdBu')+RotatedAxis()+xlab('')+ylab('')
ggsave('stp3_tme/plot/mye_analysis/mye_marker_dot.pdf',width = 7,height = 3.2)



##### proportion #####
propplot2(srtmye@meta.data,myecol,groupby = 'stage1')


propmye = getprop(srtmye$cellType_merge,srtmye@meta.data[,c('stage1','sp_stg')]) %>% 
  gather('celltype','prop',all_of(unique(srtmye$cellType_merge)))

propmye %>%
  ggplot(aes(group_1,prop))+
  theme1+
  theme(axis.ticks.x = element_blank())+rotate_x_text(angle = 30)+
  NoLegend()+
  ylab('Proportion')+xlab('')+
  scale_y_continuous(labels = c('0','0.5','1.0'),breaks = seq(0,1,0.5),limits = c(0,1.1))+
  scale_fill_manual(values = myecol)+scale_color_manual(values = myecol)+
  stat_compare_means(comparisons = makecomparison(4))+
  facet_wrap(~celltype)+
  geom_boxplot(aes(fill=celltype),alpha=0.5,outlier.color = NA)+
  geom_quasirandom(aes(color=celltype))






##### proportion split by area #####
meta=srtmye@meta.data
meta$cellType_merge = meta$cellType_merge %>% factor %>% droplevels()
propnum = table(meta[,c('cellType_merge','stage1','inf','sp_stg')]) %>% data.frame()

totalnum = aggregate(list(total=propnum$Freq),by=as.list(propnum[,c('stage1','inf','sp_stg')]),sum) %>% 
  filter(total>0)
propnum = propnum %>% inner_join(totalnum)
propnum$prop = propnum$Freq/propnum$total
propnum$prop[propnum$inf=='Out'] = -propnum$prop[propnum$inf=='Out']

propnum_mean = aggregate(list(meanvalue=propnum$prop),by=as.list(propnum[,c('stage1','inf','cellType_merge')]),median)

propnum %>%
  ggplot(aes(stage1,prop))+
  theme1+theme(axis.line.x = element_blank(),axis.ticks.x = element_blank())+
  NoLegend()+
  facet_wrap(~cellType_merge,nrow = 1)+
  scale_fill_manual(values = myecol)+scale_color_manual(values = myecol)+
  # scale_y_continuous(expand = c(0,0.3))+
  # stat_mean(geom = 'bar',aes(fill=cellType_merge),alpha=0.5)+
  xlab('')+ylab('Cell proportion')+
  geom_bar(data=propnum_mean,aes(stage1,meanvalue,fill=cellType_merge),color='black',stat='identity',position = position_dodge(),alpha=0.5)+
  geom_point(aes(color=cellType_merge),size=1.5,position=position_jitter())+
  stat_summary(data=propnum[propnum$inf=='Infiltrated',],geom = 'errorbar',
               fun.max = function(x){quantile(x,0.75)},
               fun.min = function(x){quantile(x,0.25)},
               width=0.5,linewidth=0.6,
               aes(group=cellType_merge),position = position_dodge(width = 0.9))+
  stat_compare_means(data=propnum[propnum$inf=='Infiltrated',],
                     comparisons = makecomparison(4),tip.length = 0,label.y = c(1,1.2,1.4))+
  stat_summary(data=propnum[propnum$inf=='Out',],geom = 'errorbar',
               fun.max = function(x){quantile(x,0.75)},
               fun.min = function(x){quantile(x,0.25)},
               width=0.5,linewidth=0.6,
               aes(group=cellType_merge),position = position_dodge(width = 0.9))+
  stat_compare_means(data=propnum[propnum$inf=='Out',],
                     comparisons = makecomparison(4),tip.length = 0,label.y = c(-1.2,-1.4,-1.6),vjust = 1.6)+
  geom_hline(yintercept = 0)

ggsave('stp3_tme/plot/mye_analysis/prop_mye_splitbyarea.pdf',width = 15,height = 6)





##### spatial #####
{
  objlist = setNames(list('xe1_1','xe1_2','xe1_3','xe2_1','xe2_2','xe2_3','xe3_1','xe3_2',
                          'xe4_1','xe4_2','xe5_1','xe5_2_split','xe5_3_split','xe6_1','xe6_2','xe6_3'),
                     c('T1S1','T1S2','T1S3','T2S1','T2S2','T2S3','T3S1','T3S2',
                       'T4S1','T4S2','T5S1','T5S2','T5S3','T6S1','T6S2','T6S3'))
  lapply(names(objlist), function(x){
    coltmp = c(myecol,'Epi'='#dddddd','Others'='#888888')
    dt = get(objlist[[x]])
    
    dt$cellType_merge = 'Others'
    dt$cellType_merge[dt$ct3=='Epi'] = 'Epi'
    dt = renewMeta(dt, srtmye@meta.data[,c('cellID','cellType_merge')])
    
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
    
    fname = paste0('stp3_tme/plot/mye_spatial/',x,'.png')
    p = xedimplot(dt, groupby = 'cellType_merge',flip = flip,ptlevels = c(names(myecol),'Epi','Others'),
                  reversex = reversex,reversey = reversey,color = coltmp,bgcol = 'black')
    ggsave(fname,p,width=40,height=20)
  })
}





##### density by area #####
{## grid dens
  gridid = read_rds('stp1_summary/data/allsample_gridid.rds')
  myedt = inner_join(gridid[,c('cellID','gridid')],srtmye@meta.data)
  
  gridmyenum = table(myedt$inf,myedt$gridid,myedt$stage1,myedt$sp_stg, myedt$cellType_merge) %>% 
    data.frame() %>% `colnames<-`(c('inf','gridid','stage','sp_stg','celltype','cellnum')) %>% 
    filter(cellnum>0)
  
  meanmyenum = aggregate(list(meannum=gridmyenum$cellnum),
                         by=as.list(gridmyenum[,c('inf','stage','sp_stg')]),
                         mean)
  meanmyenum_ct = aggregate(list(meannum=gridmyenum$cellnum),
                            by=as.list(gridmyenum[,c('inf','stage','sp_stg','celltype')]),
                            mean)
  meanmyenum_ct %>% 
    ggplot(aes(stage,meannum))+
    theme1+theme(axis.ticks.x = element_blank())+
    xlab('')+ylab('')+
    NoLegend()+RotatedAxis()+
    # scale_fill_manual(values = epicol)+scale_color_manual(values = epicol)+
    scale_y_continuous(limits = function(x){
      mindt = 0
      if(x[2]<=30){
        maxdt = 30
      }else{
        maxdt = 60
      }
      return(c(mindt,maxdt))
    },breaks = function(x){
      if(x[2]<=30){
        dt = seq(0,30,10)
      }else{
        dt = seq(0,60,20)
      }
      return(dt)
    },expand = c(0,0))+
    # facet_wrap(~inf,scales = 'free_y',nrow = 1)+
    facet_wrap(~inf+celltype,scales = 'free_y',nrow=2)+
    ylim(c(0,NA))+stat_compare_means(comparisons = makecomparison(4))+
    geom_jitter(size=2,aes(color=inf))+
    geom_boxplot(aes(fill=inf),outlier.color = NA,alpha=0.5)
}


{## density per cell area
  dens_all = read_rds('stp1_summary/data/allsample_density_primtype.rds')
  dens_all = allmeta %>% left_join(dens_all[,c('cellID','dens')])
  dens_mye = dens_all %>% filter(ct3=='Mye') %>% left_join(srtmye@meta.data[,c('cellID','inf')])
  meandens_mye = aggregate(list(dens=dens_mye$dens),by=as.list(dens_mye[,c('sp_stg','stage1','inf')]),mean)
  
  meandens_mye %>% 
    ggplot(aes(stage1,dens))+
    theme1+theme(axis.ticks.x = element_blank())+
    xlab('')+ylab('Cell density')+
    NoLegend()+rotate_x_text(angle = 30)+
    stat_compare_means(comparisons = makecomparison(4))+
    scale_fill_manual(values = areacol)+scale_color_manual(values = areacol)+
    scale_y_continuous(expand = c(0.1,0.3))+
    facet_wrap(~inf,scales = 'free',nrow = 2)+
    geom_jitter(size=2,aes(color=inf))+
    geom_boxplot(aes(fill=inf),outlier.color = NA,alpha=0.5)
  ggsave('stp3_tme/plot/mye_analysis/mye_density_area.pdf',width = 3,height = 7)
}


{## density per cell subtype
  objlist = setNames(list('xe1_1','xe1_2','xe1_3','xe2_1','xe2_2','xe2_3','xe3_1','xe3_2',
                          'xe4_1','xe4_2','xe5_1','xe5_2_split','xe5_3_split','xe6_1','xe6_2','xe6_3'),
                     c('T1S1','T1S2','T1S3','T2S1','T2S2','T2S3','T3S1','T3S2',
                       'T4S1','T4S2','T5S1','T5S2','T5S3','T6S1','T6S2','T6S3'))
  job({
    dens_mye = lapply(names(objlist), function(x){
      message(paste0('running sample ',x))
      crd = crds[[x]]
      crd = crd[,c('x','y')]
      meta = srtmye@meta.data[srtmye$cellID%in%rownames(crd),]
      
      dens_ct = lapply(unique(meta$cellType_merge),function(y){
        message(paste0('running ct ',y))
        meta1 = meta[meta$cellType_merge==y,]
        crd1 = crd[rownames(crd)%in%meta1$cellID,]
        
        cdens = calcdens(crd1, 100)
        rst = data.frame(cellID=rownames(crd1),dens = cdens, cellType_merge = y)
        return(rst)
      }) %>% Reduce(rbind,.)
      
      return(dens_ct)
    }) %>% Reduce(rbind,.)
    saveRDS(dens_mye,'stp3_tme/data/mye/mye_subtypeDensity.rds')
  },import = 'auto')
  
  dens_mye = srtmye@meta.data %>% left_join(dens_mye[,c('cellID','dens')])
  
  ## subtype
  meandens_mye_sub = aggregate(list(dens=dens_mye$dens),by=as.list(dens_mye[,c('sp_stg','stage1','cellType_merge')]),mean)
  meandens_mye_sub %>% 
    ggplot(aes(stage1,dens))+
    theme1+theme(axis.ticks.x = element_blank())+
    xlab('')+ylab('Cell density')+
    NoLegend()+rotate_x_text(angle = 30)+
    stat_compare_means(comparisons = makecomparison(4))+
    scale_fill_manual(values = myecol)+scale_color_manual(values = myecol)+
    scale_y_continuous(expand = c(0.1,0.3))+
    facet_wrap(~cellType_merge,scales = 'free',nrow = 1)+
    geom_jitter(size=2,aes(color=cellType_merge))+
    geom_boxplot(aes(fill=cellType_merge),outlier.color = NA,alpha=0.5)
  ggsave('stp3_tme/plot/mye_analysis/mye_density_subtype.pdf',width = 13,height = 3.5)
  
  ## subtype area
  meandens_mye_sub_area = aggregate(list(dens=dens_mye$dens),by=as.list(dens_mye[,c('sp_stg','stage1','cellType_merge','inf')]),mean)
  meandens_mye_sub_area %>% 
    ggplot(aes(stage1,dens))+
    theme1+theme(axis.ticks.x = element_blank())+
    xlab('')+ylab('Cell density')+
    NoLegend()+rotate_x_text(angle = 30)+
    stat_compare_means(comparisons = makecomparison(4))+
    scale_fill_manual(values = myecol)+scale_color_manual(values = myecol)+
    scale_y_continuous(expand = c(0.1,0.3))+
    facet_wrap(~inf+cellType_merge,scales = 'free',nrow = 2)+
    geom_jitter(size=2,aes(color=cellType_merge))+
    geom_boxplot(aes(fill=cellType_merge),outlier.color = NA,alpha=0.5)
}




##### density compare with t cell #####
{## density per cell area
  dens_all = read_rds('stp1_summary/data/allsample_density_primtype.rds')
  dens_all = allmeta %>% left_join(dens_all[,c('cellID','dens')])
  inftype = srtmye@meta.data[,c('cellID','inf')] %>% `colnames<-`(c('cellID','epiarea'))
  inftype = srttc@meta.data[,c('cellID','epiarea')] %>% rbind(inftype)
  dens_myet = dens_all %>% filter(ct3%in%c('Mye','T')) %>% left_join(inftype)
  meandens_myet = aggregate(list(dens=dens_myet$dens),by=as.list(dens_myet[,c('sp_stg','stage1','epiarea','ct3')]),mean)
  
  meandens_myet %>% 
    ggplot(aes(ct3,dens))+
    theme1+theme(axis.ticks.x = element_blank())+
    xlab('')+ylab('Cell density')+
    NoLegend()+
    stat_compare_means(comparisons = makecomparison(2))+
    scale_fill_manual(values = primcol2)+scale_color_manual(values = primcol2)+
    scale_y_continuous(expand = c(0.1,0.3))+
    facet_grid(epiarea~stage1,scales = 'free')+
    geom_jitter(size=2,aes(color=ct3))+
    geom_boxplot(aes(fill=ct3),outlier.color = NA,alpha=0.5)
  ggsave('stp3_tme/plot/mye_analysis/mye_density_compareWithT.pdf',width = 8,height = 5)
}









