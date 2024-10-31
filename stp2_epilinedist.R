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



srtepi = read_rds('data/xenium_seurat_processed/ct_epi_srt.rds')
srtepi$stage1 = srtepi$stage1 %>% factor(levels = c('NOR','LGIN','HGIN','ESCC'))

epimeta = read_rds('data/xenium_seurat_processed/ct_epi_meta.rds')
fibmeta = srtfib@meta.data
crds = read_rds('stp1_summary/data/allsample_coords.rds')
propriac = fibmeta$cellID[fibmeta$fibarea=='Propria']

epitype = fread('stp3_tme/data/all_epi_anno_meta_20240120.csv',data.table = F)
epitype = epitype %>% dplyr::select('cellID','meta_anno')





##### color #####
primcol = brewer.pal(10, 'Paired')
primcol[7] = '#999999'
names(primcol) = c('B','Epi','Endo','Fib','Gland','Mye','Mus','Mast','Pla','T')

primcol2 = c('#1f78b4','#ffc089','#B5EFB5','#793c1b',
             '#6a3d9a','#333333','#ffff33','#df928e')
names(primcol2) = c('B','Endo','Fib','Gland','Mye','Mus','T','Epi')



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




##### epi border #####
objlist = setNames(list('xe1_1','xe1_2','xe1_3','xe2_1','xe2_2','xe2_3','xe3_1','xe3_2',
                        'xe4_1','xe4_2','xe5_1','xe5_2_split','xe5_3_split','xe6_1','xe6_2','xe6_3'),
                   c('T1S1','T1S2','T1S3','T2S1','T2S2','T2S3','T3S1','T3S2',
                     'T4S1','T4S2','T5S1','T5S2','T5S3','T6S1','T6S2','T6S3'))
# bd = lapply(setNames(names(objlist),names(objlist)),function(x){
#   nms = tolower(x)
#   bd = list(d=fread(paste0('stp3_tme/data/epiborder/',nms,'_D.csv'),data.table = F), u=fread(paste0('stp3_tme/data/epiborder/',nms,'_U.csv'),data.table = F))
#   return(bd)
# })
# saveRDS(bd,'stp3_tme/data/epiborder/allborder.rds')
bd = read_rds('stp3_tme/data/epiborder/allborder.rds')



{
  x = 'T6S2'
  nms = tolower(x)
  d=fread(paste0('stp3_tme/data/epiborder/',nms,'_D.csv'),data.table = F)
  u=fread(paste0('stp3_tme/data/epiborder/',nms,'_U.csv'),data.table = F)
  borderseg = list(d = d,u = u)
  
  tmp = crds[[x]]
  (tmp %>% 
    filter(x<7000,x>6000,y<4000,y>3000,ct3=='Epi') %>% 
    ggplot(aes(x,y))+
    geom_point(size=0.5,color='lightblue')+coord_fixed()+
    geom_path(data=borderseg[['d']])) %>% plotly::ggplotly()
  
  tmpdt = c(6390.497,3557.635,6556.084,3453.649,6661.388,3354.170) %>% matrix(byrow = T,ncol = 2)
  
  poly = read_rds(paste0('stp3_tme/data/epiregion/',x,'_poly.rds'))
  tmp1 = poly[2]
  plotpoly(tmp1)
  
  f1(tmp1) %>% plotly::ggplotly()
  tmp2 = tmp1[[1]][c(511:nrow(tmp1[[1]]),1:97),] %>% 
    rbind(.,tmpdt) %>% rbind(tmp1[[1]][134:231,])
  f2(tmp2)
  d = tmp2 %>% data.frame() %>% `colnames<-`(c('x','y'))
  fwrite(d,paste0('stp3_tme/data/epiborder/',nms,'_D.csv'))
}




##### calculate epi distance #####
objlist = setNames(list('xe1_1','xe1_2','xe1_3','xe2_1','xe2_2','xe2_3','xe3_1','xe3_2',
                        'xe4_1','xe4_2','xe5_1','xe5_2_split','xe5_3_split','xe6_1','xe6_2','xe6_3'),
                   c('T1S1','T1S2','T1S3','T2S1','T2S2','T2S3','T3S1','T3S2',
                     'T4S1','T4S2','T5S1','T5S2','T5S3','T6S1','T6S2','T6S3'))
lapply(names(objlist),function(x){
  message(x)
  
  crd = crds[[x]]
  epic = crd[rownames(crd)%in%srtepi$cellID,]
  epic = epic[,c('x','y')]
  
  {# label outlier
    poly = read_rds(paste0('stp3_tme/data/epiregion/',x,'_poly.rds'))
    
    polyall = lapply(1:length(poly), function(x){
      dt = data.frame(poly[[x]])
      dt$grp = as.character(x)
      return(dt)
    }) %>% Reduce(rbind,.)
    colnames(polyall) = c('x','y','grp')
    
    tmp = lapply(poly, function(x){
      rst = isinner2(epic, x)
    }) %>% Reduce('|',.)
    
    dt = epic %>% mutate(grp=tmp)
    ggplot(polyall, aes(x,y))+
      geom_point(data = dt,aes(color=grp),size=0.1)+
      geom_polygon(aes(group=grp),alpha=0.25,color='black')+coord_fixed()
    ggsave(paste0('stp3_tme/plottmp/',x,'_outlier.png'),width = 20,height = 10)
  }
  
  
  borderseg = bd[[x]]
  segs = lapply(names(borderseg),function(x){
    seg = borderseg[[x]][,c('x','y')]
    seg$grp=x
    return(seg)
  }) %>% Reduce(rbind,.)
  
  {# get invasive epi
    if(x%in%c("T5S1","T5S2","T5S3")){
      dt$inv = FALSE
      dt %>% 
        ggplot(aes(x,y))+
        geom_point(aes(color=inv),size=0.1)+
        geom_path(data=borderseg[['u']])+
        geom_path(data=borderseg[['d']])+
        coord_fixed()
      ggsave(paste0('stp3_tme/plottmp/',x,'_invasive.png'),width = 20,height = 10)
    }else{
      segs %>% 
        ggplot(aes(x,y))+
        geom_point(data = dt,aes(color=grp),size=0.1)+
        geom_path(aes(group=grp))+coord_fixed()
      
      segd = borderseg[['d']]
      pt1 = c(segd[1,'x'],0)
      pt2 = c(segd[nrow(segd),'x'],0)
      if(x%in%c('T4S2','T6S1')){
        pt1 = c(segd[1,'x'],max(epic$y))
        pt2 = c(segd[nrow(segd),'x'],max(epic$y))
      }
      polybelow=rbind(pt1,segd) %>% rbind(pt2)
      
      system.time({
        tmp1 = isinner2(epic, polybelow)
      })
      
      dt$inv = tmp1
      dt %>% 
        ggplot(aes(x,y))+
        geom_point(aes(color=inv),size=0.1)+
        geom_path(data=borderseg[['u']])+
        geom_polygon(data=polybelow,alpha=0.3,color='black')+
        coord_fixed()
      ggsave(paste0('stp3_tme/plottmp/',x,'_invasive.png'),width = 20,height = 10)
    }
  }
  
  {# calculate dist to border
    ddist = apply(epic,1,function(x){
      return(disttoseg(x,borderseg[['d']]))
    })
    udist = apply(epic,1,function(x){
      return(disttoseg(x,borderseg[['u']]))
    })
    dt$ddist = ddist
    dt$udist = udist
    dt %>% 
      filter(grp==TRUE) %>% 
      filter(inv==FALSE) %>% 
      ggplot(aes(x,y))+
      geom_point(aes(color=ddist),size=0.1)+
      scale_color_gradientn(colors = brewer.pal(9,'Reds'))+
      geom_path(data=segs,aes(group=grp))+coord_fixed()
    ggsave(paste0('stp3_tme/plottmp/',x,'_downDist.png'),width = 20,height = 10)
    dt %>% 
      filter(grp==TRUE) %>% 
      filter(inv==FALSE) %>% 
      ggplot(aes(x,y))+
      geom_point(aes(color=udist),size=0.1)+
      scale_color_gradientn(colors = brewer.pal(9,'Reds'))+
      geom_path(data=segs,aes(group=grp))+coord_fixed()
    ggsave(paste0('stp3_tme/plottmp/',x,'_upDist.png'),width = 20,height = 10)
    
    reldist = ifelse(dt$inv,ddist/(udist-ddist),ddist/(ddist+udist))
    
    dt$reldist = reldist
    dt %>% 
      filter(grp==TRUE) %>% 
      ggplot(aes(x,y))+
      geom_point(aes(color=reldist),size=0.1)+
      scale_color_gradientn(colors = brewer.pal(9,'Reds'))+
      geom_path(data=segs,aes(group=grp))+coord_fixed()
    ggsave(paste0('stp3_tme/plottmp/',x,'_relativeDist.png'),width = 20,height = 10)
  }
  saveRDS(dt,paste0('stp3_tme/data/epi_borderdist/',x,'.rds'))
  return(NULL)
})

epidist = lapply(names(objlist),function(x){
  dt = read_rds(paste0('stp3_tme/data/epi_borderdist/',x,'.rds'))
  return(dt)
}) %>% Reduce(rbind,.)
colnames(epidist) = c('x','y','inEpi','invasive','downDist','upDist','relativeDist')
epidist$cellID = rownames(epidist)
saveRDS(epidist,'stp3_tme/data/epi_borderdist/allepi_borderdist.rds')




##### analysis #####
epidist = read_rds('stp3_tme/data/epi_borderdist/allepi_borderdist.rds')
epidist = epidist %>% left_join(srtepi@meta.data) %>% inner_join(epitype)
epidist$epithick = epidist$upDist+epidist$downDist
epidist = epidist[epidist$inEpi,]
epidist = epidist[!epidist$invasive,]
epidist = epidist %>% filter(stage1!='ESCC')
epidist$sample1 = epidist$cellID %>% str_split('_',simplify= T) %>% .[,1]

epidist_aggr = aggregate(epidist$epithick,by=as.list(epidist[,c('stage1','sp_stg')]),mean)
# epidist_aggr = aggregate(epidist$epithick,by=as.list(epidist[,c('stage1','sample1')]),mean)
epidist_aggr$stage1 = epidist_aggr$stage1 %>% factor(levels = c('NOR','LGIN','HGIN','ESCC'))
epidist_aggr$reldist = epidist_aggr$x*0.2125
epidist_aggr %>%
  .[grep('T3|T4',epidist_aggr$sp_stg,invert = T),] %>%
  # .[grep('T3',epidist_aggr$sample1,invert = T),] %>%
  ggplot(aes(stage1,reldist,color=stage1))+
  theme1+
  stat_compare_means(comparisons = makecomparison(3))+
  scale_color_manual(values = col_stage)+
  # scale_y_continuous(limits = c(45,220))+
  scale_y_continuous(limits = c(50,220))+
  NoLegend()+
  geom_jitter(width = 0.3)+
  geom_boxplot(fill=NA,outlier.color = NA)
ggsave('stp3_tme/plottmp/epithick_sp.pdf',width = 3,height = 3)




epidist_aggr1 = aggregate(as.list(epidist[,c('relativeDist','epithick')]),by=as.list(epidist[,c('stage1','sp_stg','meta_anno')]),mean)
epidist_aggr1$meta_anno = epidist_aggr1$meta_anno %>% factor(levels = c('Basal','Proliferation','Differentiation','Terminal','Non-invasive','Invasive'))
# epidist %>% 
epidist_aggr1 %>% 
  filter(stage1!='HGIN',!meta_anno%in%c('Invasive','Non-invasive')) %>% 
  ggplot(aes(stage1,relativeDist))+
  # geom_jitter()+
  geom_jitter(aes(fill=meta_anno),shape=21)+
  geom_boxplot(aes(fill=meta_anno),outlier.color = NA,alpha=0.5)+
  theme1+
  RotatedAxis()+
  NoLegend()+
  scale_color_manual(values = brewer.pal(9,'Set3'))+
  scale_fill_manual(values = brewer.pal(9,'Set3'))+
  scale_y_continuous(limits = c(0,1.1),breaks = seq(0,1,0.25))+
  stat_compare_means(comparisons = makecomparison(2))+
  facet_wrap(~meta_anno,scales = 'free_x',nrow = 1)
ggsave('stp3_tme/plottmp/reldist_program.pdf',width = 8,height = 4)



tmp = aggregate(as.list(epidist[,c('relativeDist'),drop=F]),by=as.list(epidist[,c('stage1','sp_stg','meta_anno')]),mean)
epidist_aggr2 = tmp %>% `colnames<-`(c('stage1','sp_stg','meta_anno','relativeDist-mean'))
tmp = aggregate(as.list(epidist[,c('relativeDist'),drop=F]),by=as.list(epidist[,c('stage1','sp_stg','meta_anno')]),median)
epidist_aggr2 = tmp %>% `colnames<-`(c('stage1','sp_stg','meta_anno','relativeDist-median')) %>% left_join(epidist_aggr2,.)
tmp = aggregate(as.list(epidist[,c('relativeDist'),drop=F]),by=as.list(epidist[,c('stage1','sp_stg','meta_anno')]),min)
epidist_aggr2 = tmp %>% `colnames<-`(c('stage1','sp_stg','meta_anno','relativeDist-min')) %>% left_join(epidist_aggr2,.)
tmp = aggregate(as.list(epidist[,c('relativeDist'),drop=F]),by=as.list(epidist[,c('stage1','sp_stg','meta_anno')]),max)
epidist_aggr2 = tmp %>% `colnames<-`(c('stage1','sp_stg','meta_anno','relativeDist-max')) %>% left_join(epidist_aggr2,.)
tmp = aggregate(as.list(epidist[,c('relativeDist'),drop=F]),by=as.list(epidist[,c('stage1','sp_stg','meta_anno')]),function(x){return(quantile(x,0.25))})
epidist_aggr2 = tmp %>% `colnames<-`(c('stage1','sp_stg','meta_anno','relativeDist-q25')) %>% left_join(epidist_aggr2,.)
tmp = aggregate(as.list(epidist[,c('relativeDist'),drop=F]),by=as.list(epidist[,c('stage1','sp_stg','meta_anno')]),function(x){return(quantile(x,0.75))})
epidist_aggr2 = tmp %>% `colnames<-`(c('stage1','sp_stg','meta_anno','relativeDist-q75')) %>% left_join(epidist_aggr2,.)
epidist_aggr2$meta_anno = epidist_aggr2$meta_anno %>% factor(levels = c('Basal','Proliferation','Differentiation','Terminal','Non-invasive','Invasive'))
fwrite(epidist_aggr2,'stp3_tme/allsample_layerPercent.csv',sep = ',')





##### t1s1 example #####
crds = read_rds('stp1_summary/data/allsample_coords.rds')

{## t1s1
  poly11 = read_rds('stp3_tme/data/epiregion/T1S1_poly.rds')
  
  polyall = lapply(1:length(poly11), function(x){
    dt = data.frame(poly11[[x]])
    dt$grp = as.character(x)
    return(dt)
  }) %>% Reduce(rbind,.)
  colnames(polyall) = c('x','y','grp')
  
  c11 = crds[['T1S1']]
  epic11 = c11[rownames(c11)%in%srtepi$cellID,]
  epic11 = epic11[,c('x','y')]
  
  system.time({
    tmp = lapply(poly11, function(x){
      rst = isinner2(epic11, x)
    }) %>% Reduce('|',.)
  })
  
  dt = epic11 %>% mutate(grp=tmp)
  ggplot(polyall, aes(x,y))+
    geom_polygon(aes(group=grp),alpha=0.25)+
    geom_point(data = dt,aes(color=grp),size=0.1)+coord_fixed()
  ggsave('stp3_tme/xe11_epi-outlier.png',width = 20,height = 10)
  
  
  
  borderseg = list(d=fread('stp3_tme/data/epiborder/t1s1_D.csv'),u=fread('stp3_tme/data/epiborder/t1s1_U.csv'))
  segs = lapply(names(borderseg),function(x){
    seg = borderseg[[x]][,c('x','y')]
    seg$grp=x
    return(seg)
  }) %>% Reduce(rbind,.)
  
  ## get invasive epi
  segs %>% 
    ggplot(aes(x,y))+
    geom_point(data = dt,aes(color=grp),size=0.1)+
    geom_path(aes(group=grp))+coord_fixed()
  ptbelow = c(max(epic11$x),0,0,0) %>% matrix(nrow = 2,byrow = T) %>% `colnames<-`(c('x','y'))
  ptup = c(max(epic11$x),max(epic11$y),0,max(epic11$y)) %>% matrix(nrow = 2,byrow = T) %>% `colnames<-`(c('x','y'))
  polybelow=rbind(borderseg[['d']][,c('x','y')],ptbelow)
  polybelow %>% 
    ggplot(aes(x,y))+
    geom_point(data = dt,aes(color=grp),size=0.1)+
    geom_polygon(alpha=0.5)+coord_fixed()
  
  system.time({
    tmp1 = isinner2(epic11, polybelow)
  })
  
  dt$inv = tmp1
  polybelow %>% 
    ggplot(aes(x,y))+
    geom_point(data = dt,aes(color=inv),size=0.1)+
    geom_polygon(alpha=0.5)+coord_fixed()
  
  ## calculate dist to border
  ddist = apply(epic11,1,function(x){
    return(disttoseg(x,borderseg[['d']][,c('x','y')]))
  })
  udist = apply(epic11,1,function(x){
    return(disttoseg(x,borderseg[['u']][,c('x','y')]))
  })
  dt %>% 
    filter(grp==TRUE) %>% 
    filter(inv==FALSE) %>% 
    filter(x>8000) %>% 
    ggplot(aes(x,y))+
    geom_point(aes(color=ddist),size=0.1)+
    scale_color_gradientn(colors = brewer.pal(9,'Reds'))+
    geom_path(data=segs,aes(group=grp))+coord_fixed()
  ggsave('stp3_tme/xe11_epi-downDist.png',width = 40,height = 15)
  dt %>% 
    filter(grp==TRUE) %>% 
    filter(inv==FALSE) %>% 
    filter(x>8000) %>% 
    ggplot(aes(x,y))+
    geom_point(aes(color=udist),size=0.1)+
    scale_color_gradientn(colors = brewer.pal(9,'Reds'))+
    geom_path(data=segs,aes(group=grp))+coord_fixed()
  ggsave('stp3_tme/xe11_epi-upDist.png',width = 40,height = 15)
  
  reldist = ifelse(dt$inv,ddist/(udist-ddist),ddist/(ddist+udist))
  
  dt = dt %>% mutate(udist=udist,ddist=ddist,reldist=reldist)
  dt %>% 
    filter(grp==TRUE) %>% 
    ggplot(aes(x,y))+
    geom_point(aes(color=reldist),size=0.1)+
    scale_color_gradientn(colors = brewer.pal(9,'Reds'))+
    geom_path(data=segs,aes(group=grp))+coord_fixed()
  ggsave('stp3_tme/xe11_epi-relativeDist.png',width = 40,height = 15)
  
  saveRDS(dt,'stp3_tme/xe11_epidist.rds')
  
  ### get quantile line 
  tmp = dt %>% filter(grp==TRUE) %>% filter(inv==FALSE)
  tmp = tmp %>% filter(reldist<=0.25)
  as_tmp = alphahull::ashape(tmp[,c('x','y')],alpha = 150)
  poly25 = extractpoly(as_tmp)
  tmp = dt %>% filter(grp==TRUE) %>% filter(inv==FALSE)
  tmp = tmp %>% filter(reldist<=0.50)
  as_tmp = alphahull::ashape(tmp[,c('x','y')],alpha = 150)
  poly50 = extractpoly(as_tmp)
  tmp = dt %>% filter(grp==TRUE) %>% filter(inv==FALSE)
  tmp = tmp %>% filter(reldist<=0.75)
  as_tmp = alphahull::ashape(tmp[,c('x','y')],alpha = 150)
  poly75 = extractpoly(as_tmp)
  
  polytmp = append(poly25,poly50) %>% append(poly75)
  plotpoly(polytmp)
  polytmp = lapply(1:length(polytmp),function(x){
    dt = polytmp[[x]] %>% data.frame
    colnames(dt) = c('x','y')
    dt$grp = x
    return(dt)
  }) %>% Reduce(rbind,.)
  
  dt %>% 
    filter(grp==TRUE) %>% 
    filter(inv==FALSE) %>% 
    filter(x>8000) %>% 
    ggplot(aes(x,y))+
    geom_point(aes(color=reldist),size=0.1)+
    scale_color_gradientn(colors = brewer.pal(9,'Reds'))+
    coord_fixed()+
    geom_polygon(data=polytmp,aes(group=grp),fill=NA,color='black')+
    geom_path(data=borderseg[['u']])
  ggsave('stp3_tme/xe11_epi-relativeDist_withquantline.png',width = 40,height = 15)
  
  
  ### calculate dist range 
  epimeta = read.csv('stp3_tme/data/all_epi_anno_meta_20240120.csv')
  epimeta$stage1 = epimeta$stage1 %>% factor(levels = c('NOR','LGIN','HGIN','ESCC'))
  
  epimeta11 = dt %>% mutate(cellID=rownames(.)) %>% inner_join(epimeta)
  unique(epimeta11$meta_anno)
  epimeta11 %>% 
    filter(stage1%in%c('NOR','LGIN')) %>% 
    ggplot(aes(stage1,reldist))+
    theme1+
    RotatedAxis()+
    scale_fill_manual(values = brewer.pal(9,'Set3'))+
    stat_compare_means(comparisons = makecomparison(2))+
    facet_wrap(~meta_anno)+
    geom_boxplot(aes(fill=meta_anno))
}




##### get dist line ##### 
epidist = read_rds('stp3_tme/data/epi_borderdist/allepi_borderdist.rds')
getDistline = function(dt, distances){
  crd = dt[,c('x','y')]
  
  loopundone=T
  alp=10
  while(loopundone){
    for(i in distances){
      dt1 = dt %>% filter(downDist<i)
      crd1 = dt1[,c('x','y')]
      as = alphahull::ashape(crd1,alpha = alp)
      R.utils::withTimeout({
        polytmp = extractpoly(as)
      },)
      plotpoly(polytmp)+geom_point(data=crd1,aes(x,y),size=0.2)
      
    }
    alp=alp+10
  }
}


epitype = fread('stp3_tme/data/all_epi_anno_meta_20240120.csv',data.table = F)
tmp = epitype %>% filter(sp_stg%in%c('T1-NOR-2','T1-LGIN-2'))

propplot(tmp$meta_anno,tmp$stage1)




