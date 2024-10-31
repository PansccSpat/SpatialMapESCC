source('~/xenium/header.R')
source('~/xenium/xenium_function.R')
library(scatterpie)




# input #####
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


srtepi = read_rds('data/xenium_seurat_processed/ct_epi_srt_2.rds')
epidt = read_rds('stp2_epi/data/epiid_correctstage.rds')
crds = read_rds('stp1_summary/data/allsample_coords.rds')
epipoly = read_rds('stp3_tme/data/epiregion/allsample_poly.rds')

borderdist = read_rds('stp3_tme/data/epi_borderdist/allepi_borderdist.rds')
allborder = read_rds('stp3_tme/data/epiborder/allborder.rds')
t1dist = read_rds('stp3_tme/data/epi_borderdist/T1epi_borderdist_v2.rds')
t1border = read_rds('stp3_tme/data/epiborder/T1border_v2.rds')
t3dist = read_rds('stp3_tme/data/epi_borderdist/T3ESCC_borderdist_v2.rds')
t3border = read_rds('stp3_tme/data/epiborder/T3ESCCborder_v2.rds')
t4dist = read_rds('stp3_tme/data/epi_borderdist/T4partE_borderdist_v2.rds')
t4dtlist = read_rds('stp3_tme/data/epiborder/T4partEborder_v2.rds')

epianno = fread('stp1_summary/data/xenium_epi.csv')

epigeneannot = fread('stp3_tme/data/hm_gene_anno.txt')

probelist = fread('data/probe_list.txt')

gridid = read_rds('stp1_summary/data/allsample_gridid.rds')

epiblacklist = fread('stp2_epi/data/gene_otherTissue.txt')
epiblacklist = probelist[probelist$gene%in%epiblacklist$Gene_otherTissue,]


allmk = c('EPCAM','SFN','KRT4','KRT5', # epi
          'FN1','DCN','COL1A1','COL1A2','COL3A1','COL6A1', # fibro
          'PLN','MYH11', # muscle
          'VWF','PECAM1','ENG','CDH5', # endo
          'PTPRC', # immune
          'CD2','CD3D','CD3E', # t cell
          'CD1E','CD68','LYZ','CD14','FCGR3A', # myeloid
          'MS4A2','KIT','CPA3', # mast
          'CD19','CD79A','MS4A1', # B
          'JCHAIN','IGHG2' # plasma
)

epimarker = c('KRT15','MYC', # QP
              'TOP2A','MKI67', # CY
              'ANXA1','KRT13','KRT4','CSTA','S100A8','TRIM29', # MD
              'SPRR3','SPRR2E','ECM1','EMP1', # TD
              'HIF1A','SERPINE1','VEGFA', # HY
              'KRT17', # RO
              'SPP1','CES1', 'AKR1C3','AKR1C2','CBR1', # DO
              'HLA-DRB1','HLA-DQB1'# AP
)

epicol = c('#0339f8','#ff9408','#75bbfd','#de0c62')
names(epicol) = c('Basal','Proliferation','Differentiation','Invasive')

progcol = c('Cell adhesion'="#0339f8",
            'Cycling'="#ff9408",
            'Mucosal defense'="#75bbfd",
            'Inflammation'="#B3DE69",
            'DNA repair'="#a6eee6",
            'Oncogenic TFs'="#9467bd",
            'Angiogenesis'="#41ab5d",
            'EMT'="#ff0000")

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



supp2 = fread('data/supptable/supp2_0422_reshaped.txt',data.table = F)
supp2_epig = supp2 %>% filter(inxe=='Yes')
supp2_epig = supp2_epig[grep('^Marker|Other',supp2_epig$ann_new,invert = T),]

prognew = read_rds('stp2_epi/data/epiProgram/epigene_annotation_v5.rds')


allborder = read_rds('stp3_tme/data/epiborder/allborder.rds')
t1dist = read_rds('stp3_tme/data/epi_borderdist/T1epi_borderdist_v2.rds')
t1border = read_rds('stp3_tme/data/epiborder/T1border_v2.rds')
t3dist = read_rds('stp3_tme/data/epi_borderdist/T3ESCC_borderdist_v2.rds')
t3border = read_rds('stp3_tme/data/epiborder/T3ESCCborder_v2.rds')
t4dist = read_rds('stp3_tme/data/epi_borderdist/T4partE_borderdist_v2.rds')
t4dtlist = read_rds('stp3_tme/data/epiborder/T4partEborder_v2.rds')





# gene expr spatial #####
{# balck bg
  objlist = setNames(list('xe1_1','xe1_2','xe1_3','xe2_1','xe2_2','xe2_3','xe3_1','xe3_2',
                          'xe4_1','xe4_2','xe5_1','xe5_2_split','xe5_3_split','xe6_1','xe6_2','xe6_3'),
                     c('T1S1','T1S2','T1S3','T2S1','T2S2','T2S3','T3S1','T3S2',
                       'T4S1','T4S2','T5S1','T5S2','T5S3','T6S1','T6S2','T6S3'))
  lapply(names(objlist), function(x){
    dt = get(objlist[[x]])
    
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
    
    fname = paste0('stp1_summary/plot/',x,'_primarytype.png')
    p = xedimplot(dt, groupby = 'ct3',flip = flip,reversex = reversex,reversey = reversey,color = primcol2,bgcol = 'black',pt.size = 0.1)+NoLegend()

    fname = paste0('stp1_summary/stage-prim_plot/',x,'_primarytype.png')
    dt$tmp = paste0(dt$ct3,'-',dt$stage1)
    p = xedimplot(dt, groupby = 'tmp',flip = flip,reversex = reversex,reversey = reversey,color = allcol,bgcol = 'black',pt.size = 0.1)+NoLegend()
  })
}



{
  exprcol_gene1 = c('#dfe0fc','white','#E41A1C')
  exprcol_gene1 = c(colorRampPalette(c('white','#6baed6','#08519c'))(67),colorRampPalette(c('#c994c7','#ae017e'))(33))
  exprcol_gene1 = c('white','#E41A1C','#312271')
  exprcol_gene1 = c('white','#E41A1C')
  
  gns = c('KRT5','COL17A1','NOTCH1','JAG1','MKI67','TOP2A','MYC','SOX2','TP63','ANXA1','ECM1','EMP1','HIF1A','VEGFA','MMP11','MMP2')
  
  epnl = xe3_1 %>% subset(cells=intersect(.$cellID[.$ct3=='Epi'],epidt$cellID))
  eplh = xe3_2 %>% subset(cells=intersect(.$cellID[.$ct3=='Epi'],epidt$cellID))
  epht = xe1_1 %>% subset(cells=intersect(.$cellID[.$ct3=='Epi'],epidt$cellID))
  
  roinl = c(10000,12000,2100,3500)
  roilh = c(6100,8200,2800,4000)
  roiht = c(6800,8800,3500,5300)
  
  p1 = lapply(gns,function(x){
    xefeatureplot(epnl,feature = x,bgcol = 'white',zoomin = roinl,flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = exprcol_gene1,collim = c(0,4))
  })
  p1stg = xedimplot(epnl,groupby = 'stage1',bgcol = 'white',zoomin = roinl,flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = col_stage)
  ggsave('stp2_epi/plot/epiProgram_spatial/pergene/nl.png',cowplot::plot_grid(plotlist = append(p1,list(p1stg)),nrow=5),width = 20,height = 20,limitsize = F)
  
  p2 = lapply(gns,function(x){
    xefeatureplot(eplh,feature = x,bgcol = 'white',zoomin = roilh,flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = exprcol_gene1,collim = c(0,4))
  })
  p2stg = xedimplot(eplh,groupby = 'stage1',bgcol = 'white',zoomin = roilh,flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = col_stage)
  ggsave('stp2_epi/plot/epiProgram_spatial/pergene/lh.png',cowplot::plot_grid(plotlist = append(p2,list(p2stg)),nrow=5),width = 20,height = 20,limitsize = F)
  
  p3 = lapply(gns,function(x){
    xefeatureplot(epht,feature = x,bgcol = 'white',zoomin = roiht,flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = exprcol_gene1,collim = c(0,4))
  })
  p3stg = xedimplot(epht,groupby = 'stage1',bgcol = 'white',zoomin = roiht,flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = col_stage)
  ggsave('stp2_epi/plot/epiProgram_spatial/pergene/ht.png',cowplot::plot_grid(plotlist = append(p3,list(p3stg)),nrow=5),width = 20,height = 20,limitsize = F)
}





# prog score spatial new #####
prognew = read_rds('stp2_epi/data/epiProgram/epigene_annotation_v5.rds')

{## subset data ####
  t4dist = read_rds('stp3_tme/data/epi_borderdist/T4partE_borderdist_v2.rds')
  t4dtlist = read_rds('stp3_tme/data/epiborder/T4partEborder_v2.rds')
  t1dist = read_rds('stp3_tme/data/epi_borderdist/T1epi_borderdist_v2.rds')
  t1border = read_rds('stp3_tme/data/epiborder/T1border_v2.rds')
  
  dt1 = xe4_1 %>% subset(cells=intersect(.$cellID[.$ct3=='Epi'],epidt$cellID))
  dt2 = xe5_2_split %>% subset(cells=intersect(.$cellID[.$ct3=='Epi'],epidt$cellID))
  dt3 = xe1_2 %>% subset(cells=intersect(.$cellID[.$ct3=='Epi'],epidt$cellID))
  dt4 = xe4_2 %>% subset(cells=rownames(t4dtlist$T4S2$E2$crd))
  dt1@meta.data = dt1@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% `rownames<-`(.$cellID)
  dt2@meta.data = dt2@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% `rownames<-`(.$cellID)
  dt3@meta.data = dt3@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% `rownames<-`(.$cellID)
  dt4@meta.data = dt4@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% `rownames<-`(.$cellID)
  
  roin = c(13750,14750,4100,4900)
  roil = c(2500,3500,2500,4000)
  roih = c(15500,17000,2500,4500)
  roit = getcoordsrange(dt4,flip = T)
  xedimplot(dt1,'celltype',flip = T,reversex = T,axis = T,zoomin = roin)
  xedimplot(dt2,'celltype',axis = T,zoomin = roil)
  xedimplot(dt3,'celltype',flip = T,reversex = T,axis = T,zoomin = roih)
  xedimplot(dt4,'celltype',flip = T,reversey = T,axis = T,zoomin = roit)
  
  dt1 = subset(dt1, cells = getcellinrange(dt1,roin,T))
  dt2 = subset(dt2, cells = getcellinrange(dt2,roil))
  dt3 = subset(dt3, cells = getcellinrange(dt3,roih,T))
  dt4 = subset(dt4, cells = getcellinrange(dt4,roit,T))
}

{## high expr program ####
  epiprogscore = read_rds('stp2_epi/data/epiProgram/epiprogscore_v5.rds')
  epiprogscore$stage1 = epiprogscore$stage1 %>% factor(levels = c('NOR','LGIN','HGIN','ESCC'))
  
  epiprogscore_scale = epiprogscore[,unique(prognew$annotation)]
  epiprogscore_scale = scale(epiprogscore_scale)
  progname = colnames(epiprogscore_scale)
  epiprog_max = apply(epiprogscore_scale,1,function(x){return(progname[which.max(x)])})
  
  epiprogscore$maxprog = epiprog_max
  
  epiprogscore %>% 
    filter(cellID%in%epidt$cellID) %>% 
    ggplot(aes(maxprog))+
    theme1+RotatedAxis()+scale_fill_manual(values = progcol)+
    geom_histogram(stat = 'count',aes(fill=maxprog))+
    facet_wrap(~stage1+celltype,scale='free',nrow=2)
  
  {# plot max program
    dt1@meta.data = dt1@meta.data %>% left_join(epiprogscore[,c('cellID','maxprog')]) %>% `rownames<-`(.$cellID)
    dt2@meta.data = dt2@meta.data %>% left_join(epiprogscore[,c('cellID','maxprog')]) %>% `rownames<-`(.$cellID)
    dt3@meta.data = dt3@meta.data %>% left_join(epiprogscore[,c('cellID','maxprog')]) %>% `rownames<-`(.$cellID)
    dt4@meta.data = dt4@meta.data %>% left_join(epiprogscore[,c('cellID','maxprog')]) %>% `rownames<-`(.$cellID)
    
    xedimplot(dt1,'celltype',flip = T,reversex = T,axis = T,gridline = T,zoomin = roin)
    roi = c(14100,14200,4250,4400)
    xedimplot(dt1,'celltype',flip = T,reversex = T,axis = T,zoomin = roi)
    tmp = subset(dt1, cells = getcellinrange(dt1,roi,T))
    
    xedimplot(dt1,'maxprog',color=progcol,flip = T,reversex = T,axis = F,pt.size = 2)
    xedimplot(dt2,'maxprog',color=progcol,axis = F,pt.size = 2)
    xedimplot(dt3,'maxprog',color=progcol,flip = T,reversex = T,axis = F,pt.size = 2)
    xedimplot(dt4,'maxprog',color=progcol,flip = T,reversey = T,axis = F,pt.size = 1.5)
  }
}

{##  prop of high expr program ####
  library(scatterpie)
  epiprogscore = read_rds('stp2_epi/data/epiProgram/epiprogscore_v5.rds')
  
  tmp = epiprogscore %>% filter(cellID%in%epidt$cellID) %>% 
    select(unique(prognew$annotation),cellID, stage1, sp_stg, celltype) %>% 
    gather('prog','score',all_of(unique(prognew$annotation)))
  
  epiprog_prop = aggregate(list(thresh=tmp$score),by=as.list(tmp[,c('prog'),drop=F]),function(x) mean(x)+0.5*sd(x))
  epiprog_prop = left_join(tmp,epiprog_prop) %>% mutate(highexpr=ifelse(score>thresh,1,0))
  epiprog_prop = epiprog_prop %>% select(cellID,prog,highexpr) %>% spread(prog,highexpr)
  
  {# plot program prop
    plotdt = xe1_2 %>% subset(cells=epidt$cellID)
    plotdt@meta.data = plotdt@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% `rownames<-`(.$cellID)
    xedimplot(plotdt,'celltype',flip = T,reversex = T,axis = T,gridline = T)
    roi = c(22150,22250,3400,3500)
    xedimplot(plotdt,'celltype',color=epicol,flip = T,bgcol = 'black',gridline = T,gridcol = 'white',reversex = T,axis = T,zoomin = c(22000,22500,3400,4000))
    xedimplot(plotdt,'celltype',flip = T,reversex = T,axis = T,zoomin = roi)
    plotdt = subset(plotdt, cells = getcellinrange(plotdt,roi,T))
    
    plotdt = crds$T1S2[plotdt$cellID,c('x','y','cellID')] %>% left_join(epiprog_prop)
    plotdt$radius = 3
    ggplot()+
      theme1+
      coord_fixed()+scale_x_reverse()+
      scale_fill_manual(values = progcol)+
      geom_scatterpie(data=plotdt,aes(x=x,y=y,group=cellID,r=radius),cols=colnames(epiprog_prop)[-1],color=NA)
    
    plotdt = dt4
    xedimplot(plotdt,'celltype',flip = T,reversey = T,axis = T,gridline = T)
    roi = c(18650,18800,2450,2600)
    xedimplot(plotdt,'celltype',flip = T,reversey = T,axis = T,xintercept = roi[1:2],yintercept = roi[3:4])
    xedimplot(plotdt,'celltype',flip = T,reversey = T,axis = T,zoomin = roi)
    plotdt = subset(plotdt, cells = getcellinrange(plotdt,roi,T))
    
    plotdt = crds$T4S2[plotdt$cellID,c('x','y','cellID')] %>% left_join(epiprog_prop)
    plotdt$radius = 3
    ggplot()+
      theme1+
      coord_fixed()+scale_y_reverse()+
      scale_fill_manual(values = progcol)+
      geom_scatterpie(data=plotdt,aes(x=x,y=y,group=cellID,r=radius),cols=colnames(epiprog_prop)[-1],color=NA)
  }
}

{## spatial score 2 ####
  {# preprocessing
    t4tmp = aggregate(list(mx=t4dist$borderDist),by=as.list(t4dist[,c('sid','areaID'),drop=F]),max)
    t4tmp = left_join(t4dist,t4tmp)
    t4tmp$relativeDist = t4tmp$borderDist / t4tmp$mx
    t4dist = t4dist %>% left_join(t4tmp[,c('cellID','relativeDist')])
    
    distdt = read_rds('stp3_tme/data/epi_borderdist/allepi_borderdist.rds') %>% filter(inEpi,cellID%in%t1dist$cellID=='FALSE')
    distdt = t1dist %>% filter(inEpi) %>% rbind(distdt)
    distdt = distdt %>% inner_join(epidt[,c('cellID','celltype','stage1','sp_stg')])
    distdt$dst = distdt$relativeDist 
    distdt = distdt %>% filter(stage1!='ESCC',dst<=1)
    distdt = distdt %>% dplyr::select(all_of(c('cellID','celltype','stage1','sp_stg','dst')))
    
    distdt1 = t4dist[,c('cellID','relativeDist')] %>% `colnames<-`(c('cellID','dst')) %>% 
      inner_join(epidt[,c('cellID','celltype','stage1','sp_stg')]) %>% dplyr::select(all_of(c('cellID','celltype','stage1','sp_stg','dst')))
    distdt = rbind(distdt,distdt1)
    
    grouplength = 0.01
    rownames(distdt) = distdt$cellID
    distval = distdt$dst
    groupdt = cut(distval,breaks = seq(floor(min(distval)/grouplength)*grouplength, ceiling(max(distval)/grouplength)*grouplength,grouplength),include.lowest = T)
    groupdt = data.frame(group=groupdt)
    groupdt$dn = groupdt$group %>% gsub('\\(|\\]|\\[','',.) %>% strsplit(',') %>% sapply('[',1) %>% as.numeric()
    groupdt$up = groupdt$group %>% gsub('\\(|\\]|\\[','',.) %>% strsplit(',') %>% sapply('[',2) %>% as.numeric()
    groupdt$grouploc = rowMeans(groupdt[,c('up','dn')])
    distdt$grouploc = groupdt$grouploc
  }
  
  prognew = read_rds('stp2_epi/data/epiProgram/epigene_annotation_v5.rds')
  epiprogscore = read_rds('stp2_epi/data/epiProgram/epiprogscore_v5.rds')
  progorder = c('Cell adhesion','Mucosal defense','Cycling','Inflammation','DNA repair','EMT','Oncogenic TFs','Angiogenesis')
  
  collist = lapply(names(progcol),function(p){
    col = progcol[[p]]
    col1 = colorRampPalette(c(col,'black'))(3)[2]
    cols = c('#dfe0fc','#ececfd','white',col,col1)
    return(cols)
  }) %>% `names<-`(names(progcol))
  pprog = function(dt, pname, flip=T,reversex=F,reversey=F,collim=NULL){
    xefeatureplot(dt,pname,color=collist[[pname]],flip = flip,reversey = reversey,reversex = reversex,bgcol = 'white',collim=collim,
                  gridline = T,pt.size = 0.5,pt.shape = 16)
  }
  
  
  {# hgin t1s3
    roih = c(10300,12600,2100,3500)
    dth = subset(xe1_3, cells=intersect(getcellinrange(xe1_3,roih,T),t1dist$cellID[t1dist$inEpi]))
    dth = subset(dth,cells=epidt$cellID)
    dth@meta.data = dth@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% 
      left_join(epiprogscore[,c('cellID',progorder)]) %>% `rownames<-`(.$cellID)
    xedimplot(dth,'celltype',flip = T,reversex = T,axis=T,gridline = T,xintercept = roih[1:2],yintercept = roih[3:4],zoomin = roih)
    
    ## generate contour line 
    pts = data.frame(x=runif(300*300,0,2300)+roih[1],y=runif(300*300,0,1400)+roih[3])
    bd1=data.frame(t1border$T1S3$u) %>% `colnames<-`(c('x','y')); bd2=data.frame(t1border$T1S3$d) %>% `colnames<-`(c('x','y'))
    poly = rbind(bd2,bd1[nrow(t1border$T1S3$u):1,])
    plotpoly(list(poly))
    job({pts = pts[isinner2(pts,poly),]},import = 'auto')
    
    xedimplot(dth,'celltype',flip = T,reversex = T,axis=T,xintercept = roih[1:2],yintercept = roih[3:4],zoomin = roih)+
      geom_point(data=pts,size=0.1)
    bd2_ = bd2 %>% filter(x<=roih[2]+500,x>=roih[1]-500)
    job({ptsdist = apply(pts,1,function(x){disttoseg(x,bd2_)})},import = 'auto')
    pts %>% ggplot(aes(x,y))+geom_point(aes(color=ptsdist))+
      scale_color_gradient(low = 'white',high = 'red')+
      geom_path(data=bd2_)
    
    distgroup = max(ptsdist)/10
    bdlist = list()
    for(i in 1:10){
      print(i)
      tmp = pts[ptsdist>distgroup*(i-1),]
      as = alphahull::ashape(tmp,alpha = 20)
      plot(as)
      bdlist[[i]] = extractpoly(as)[[1]] %>% data.frame()
      colnames(bdlist[[i]]) = c('x','y')
      bdlist[[i]]$group = paste0('g',i)
    }
    bds_h = Reduce(rbind,bdlist)
    
    ## spatial plot 
    plist_h = list()
    for(prg in progorder){
      collim=NULL
      collim=switch(prg,
             # 'Cell adhesion'=c(0,1.5),
             # 'Mucosal defense'=c(0,1.5),
             'Inflammation'=c(0,0.8),
             # 'Cycling'=c(0,1.5),
             # 'DNA repair'=c(0,1),
             'Angiogenesis'=c(0,1.2),
             # 'Oncogenic TFs'=c(0,2),
             'EMT'=c(0,0.5))

      plist_h[[prg]] = pprog(dth,prg,T,T,collim = collim)
    }
    cowplot::plot_grid(plotlist = plist_h,ncol = 4)
    ggsave('stp2_epi/plot/epiProgram_spatial/progscore/hgin_progSpatial.pdf',width = 20,height = 10)
    
    ## trend plot 
    tmp_h = dth
    tmp_h@meta.data = tmp_h@meta.data %>% left_join(distdt[dth$cellID,c('cellID','dst','grouploc')])
    tmp_h = tmp_h@meta.data %>% gather('prog','score',unique(prognew$annotation))
    tmp_h = aggregate(list(score=tmp_h$score),by=as.list(tmp_h[,c('prog','grouploc','celltype')]),mean)
    
    trendtest = lapply(unique(tmp_h$prog),function(prg){
      print(prg)
      dt2 = tmp_h %>% filter(prog==prg)
      testrst = mk.test(dt2$score)
      rst2 = data.frame(sval=testrst$estimates[1],zval=testrst$statistic,pval=testrst$p.value,prog=prg,stage1=stg)
    }) %>% Reduce(rbind,.)
    trendtest$pval_adj = p.adjust(trendtest$pval,method = 'fdr')
    
    labdt_h = tmp_h
    labdt_h$grouploc = labdt_h$grouploc %>% as.character %>% as.numeric()
    labdt_h = aggregate(list(yloc=labdt_h$score),by=as.list(labdt_h[,c('prog'),drop=F]),max) %>% left_join(trendtest)
    
    tmp_h %>% 
      filter(celltype=='Invasive') %>% 
      ggplot(aes(grouploc,score))+
      geom_point()+
      stat_cor(method = 'kendall')+
      # geom_text(data=labdt_h,x=0.1,aes(y=yloc*1.1,label=paste0('p = ',signif(pval_adj,2),'\nz = ',signif(zval,2))),hjust=0,vjust=1)+
      facet_wrap(~prog,scales = 'free')
  }
  
  {# escc t4s1
    dte = xe4_2 %>% subset(cells=rownames(t4dtlist$T4S2$E2$crd))
    dte@meta.data = dte@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% 
      left_join(epiprogscore[,c('cellID',unique(prognew$annotation))]) %>% `rownames<-`(.$cellID)
    # roit = c(19650,19950,2650,2950)
    roit = c(18400,20400,1500,2300)
    xedimplot(dte,'celltype',flip = T,axis=T,gridline = T,xintercept = roit[1:2],yintercept = roit[3:4])
    dte_ = subset(dte, cells = getcellinrange(dte,roit,T))
    
    ## plot complete fov
    plist_t = list()
    for(prg in progorder){
      collim=NULL
      collim=switch(prg,
                    # 'Cell adhesion'=c(0,1.5),
                    # 'Mucosal defense'=c(0,1.5),
                    'Inflammation'=c(0,0.7),
                    # 'Cycling'=c(0,1.5),
                    # 'DNA repair'=c(0,1),
                    'Angiogenesis'=c(0,0.7),
                    # 'Oncogenic TFs'=c(0,2),
                    'EMT'=c(0,0.4))
      plist_t[[prg]] = pprog(dte,prg,T,collim = collim,reversey = T)+theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(color=NA,fill='transparent'))
    }
    cowplot::plot_grid(plotlist = plist_t,ncol = 4)
    ggsave('stp2_epi/plot/epiProgram_spatial/progscore/escc_progSpatial_full.pdf',width = 20,height = 10)
    
    
    ## generate contour line 
    pts = data.frame(x=runif(300*300,0,2000)+roit[1],y=runif(300*300,0,800)+roit[3])
    poly = t4dtlist$T4S2$E2$poly %>% data.frame() %>% `colnames<-`(c('x','y'))
    plotpoly(list(poly))
    job({pts = pts[isinner2(pts,poly),]},import = 'auto')
    
    xedimplot(dte,'celltype',flip = T,reversey=T,axis=T,xintercept = roit[1:2],yintercept = roit[3:4],zoomin = roit)+
      geom_point(data=pts,size=0.1)
    bd_ = poly
    job({ptsdist = apply(pts,1,function(x){disttoseg(x,bd_)})},import = 'auto')
    pts %>% ggplot(aes(x,y))+geom_point(aes(color=ptsdist))+
      theme1+theme(panel.background = element_rect(fill='black'))+
      scale_color_gradient(low = 'white',high = 'red')+
      geom_path(data=bd_,color='grey')
    
    distgroup = max(ptsdist)/10
    bdlist = list()
    for(i in 1:10){
      print(i)
      tmp = pts[ptsdist>distgroup*(i-1),]
      as = alphahull::ashape(tmp,alpha = 20)
      plot(as)
      bdlist[[i]] = extractpoly(as)[[1]] %>% data.frame()
      colnames(bdlist[[i]]) = c('x','y')
      bdlist[[i]]$group = paste0('g',i)
    }
    bds_t = Reduce(rbind,bdlist)
    
    plist_t = list()
    for(prg in progorder){
      if(prg=='EMT'){
        collim=c(0,0.4)
      }else{collim=NULL}
      plist_t[[prg]] = pprog(dte_,prg,T,collim = collim)+geom_path(data=bds_t,aes(group=group),linetype='dashed',linewidth=0.5)
    }
    cowplot::plot_grid(plotlist = plist_t,ncol = 4)
    ggsave('stp2_epi/plot/epiProgram_spatial/progscore/escc_progSpatial.pdf',width = 20,height = 10)
    
    ## trend plot 
    tmp_e = dte_
    tmp_e@meta.data = tmp_e@meta.data %>% left_join(distdt[dte_$cellID,c('cellID','dst','grouploc')])
    tmp_e = tmp_e@meta.data %>% gather('prog','score',unique(prognew$annotation))
    tmp_e = aggregate(list(score=tmp_e$score),by=as.list(tmp_e[,c('prog','grouploc','celltype')]),mean)
    
    trendtest = lapply(unique(tmp_e$prog),function(prg){
      print(prg)
      dt2 = tmp_e %>% filter(prog==prg)
      dt2 = dt2 %>% filter(celltype=='Invasive')
      testrst = mk.test(dt2$score)
      rst2 = data.frame(sval=testrst$estimates[1],zval=testrst$statistic,pval=testrst$p.value,prog=prg,stage1=stg)
    }) %>% Reduce(rbind,.)
    trendtest$pval_adj = p.adjust(trendtest$pval,method = 'fdr')
    
    labdt_h = tmp_e
    labdt_h$grouploc = labdt_h$grouploc %>% as.character %>% as.numeric()
    labdt_h = aggregate(list(yloc=labdt_h$score),by=as.list(labdt_h[,c('prog'),drop=F]),max) %>% left_join(trendtest)
    
    tmp_e %>% 
      filter(celltype=='Invasive') %>% 
      ggplot(aes(grouploc,score))+
      geom_point()+
      stat_cor(method = 'kendall')+
      # geom_text(data=labdt_h,x=0.1,aes(y=yloc*1.1,label=paste0('p = ',signif(pval_adj,2),'\nz = ',signif(zval,2))),hjust=0,vjust=1)+
      facet_wrap(~prog,scales = 'free')
  }
}

{## spatial score ####
  epiprogscore = read_rds('stp2_epi/data/epiProgram/epiprogscore_v5.rds')
  
  dt1@meta.data = dt1@meta.data %>% dplyr::select(-contains(unique(prognew$annotation))) %>% left_join(epiprogscore[,c('cellID',unique(prognew$annotation))]) %>% `rownames<-`(.$cellID)
  dt2@meta.data = dt2@meta.data %>% dplyr::select(-contains(unique(prognew$annotation))) %>% left_join(epiprogscore[,c('cellID',unique(prognew$annotation))]) %>% `rownames<-`(.$cellID)
  dt3@meta.data = dt3@meta.data %>% dplyr::select(-contains(unique(prognew$annotation))) %>% left_join(epiprogscore[,c('cellID',unique(prognew$annotation))]) %>% `rownames<-`(.$cellID)
  dt4@meta.data = dt4@meta.data %>% dplyr::select(-contains(unique(prognew$annotation))) %>% left_join(epiprogscore[,c('cellID',unique(prognew$annotation))]) %>% `rownames<-`(.$cellID)
  
  
  exprcol_gene1 = c('#dfe0fc','white',colorRampPalette(c('white','#E41A1C'))(5)[-1])
  exprcol_gene1 = c('#dfe0fc','white','#E41A1C')
  exprcol_gene1 = c('#ccebc5','white','#08519c')
  exprcol_gene1 = c('#dfe0fc','#ececfd','white','#E41A1C','#a50f15')
  
  proi = function(dt,pathname,color,flip,reversex=F,reversey=F,roi,collim=NULL,bds=NULL){
    # color=colorRampPalette(c(color,'black'))(100)[1:60]
    p1=xefeatureplot(dt,pathname,color=color,flip = flip,reversey = reversey,reversex = reversex,bgcol = 'white',collim=collim,
                     gridline = T,axis = T,pt.size = 0.5,xintercept = roi[1:2],yintercept = roi[3:4])
    if(!is.null(bds)){
      p2=xefeatureplot(dt,pathname,color=color,flip = flip,reversey = reversey,reversex = reversex,bgcol = 'white',collim=collim,
                       pt.size = 1,zoomin = roi,cellborder = T,celllinecol = 'black',linewidth = 0.2,areaalpha = 0.7)+
        geom_path(data=bds,aes(x,y,group=group),linetype='dashed',linewidth=0.5)
    }else{
      p2=xefeatureplot(dt,pathname,color=color,flip = flip,reversey = reversey,reversex = reversex,bgcol = 'white',collim=collim,
                       pt.size = 1,zoomin = roi,cellborder = T,celllinecol = 'black',linewidth = 0.2,areaalpha = 0.7)
    }
    p1+p2
  }
  
  roih = c(16150,16300,2980,3130)
  roih = c(16150,16450,2700,3000)
  proi(dt3,'Oncogenic TFs',exprcol_gene1,T,T,F,roih)
  proi(dt3,'EMT',exprcol_gene1,T,T,F,roih)
  proi(dt3,'Cycling',exprcol_gene1,T,T,F,roih)
  
  
  {# hgin t1s3
    dth = xe1_3 %>% subset(cells=epidt$cellID)
    dth@meta.data = dth@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% 
      left_join(epiprogscore[,c('cellID',unique(prognew$annotation))]) %>% `rownames<-`(.$cellID)
    roih = c(10300,12600,2100,3500)
    xedimplot(dth,'celltype',flip = T,reversex = T,axis=T,gridline = T,xintercept = roih[1:2],yintercept = roih[3:4],zoomin = roih)
    dth = subset(dth, cells=intersect(getcellinrange(dth,roih,T),t1dist$cellID[t1dist$inEpi]))
    
    roih = c(11300,11600,2200,2500)
    
    ## generate contour line 
    pts = data.frame(x=runif(150*150,0,300)+roih[1],y=runif(150*150,0,300)+roih[3])
    bd1=data.frame(t1border$T1S3$u) %>% `colnames<-`(c('x','y')); bd2=data.frame(t1border$T1S3$d) %>% `colnames<-`(c('x','y'))
    poly = rbind(bd2,bd1[nrow(t1border$T1S3$u):1,])
    plotpoly(list(poly))
    job({pts = pts[isinner2(pts,poly),]},import = 'auto')
    
    xedimplot(dth,'celltype',flip = T,reversex = T,axis=T,xintercept = roih[1:2],yintercept = roih[3:4],zoomin = roih)+
      geom_point(data=pts,size=0.1)
    bd2_ = bd2 %>% filter(x<=roih[2]+500,x>=roih[1]-500)
    job({ptsdist = apply(pts,1,function(x){disttoseg(x,bd2_)})},import = 'auto')
    pts %>% ggplot(aes(x,y))+geom_point(aes(color=ptsdist))+
      scale_color_gradient(low = 'white',high = 'red')+
      geom_path(data=bd2_)
    
    distgroup = max(ptsdist)/10
    bdlist = list()
    for(i in 1:10){
      print(i)
      tmp = pts[ptsdist>distgroup*(i-1),]
      as = alphahull::ashape(tmp,alpha = 20)
      plot(as)
      bdlist[[i]] = extractpoly(as)[[1]] %>% data.frame()
      colnames(bdlist[[i]]) = c('x','y')
      bdlist[[i]]$group = paste0('g',i)
    }
    bds_h = Reduce(rbind,bdlist)
    
    
    proi(dth,'Cell adhesion',exprcol_gene1,T,T,F,roih,bds=bds_h)
    proi(dth,'Mucosal defense',exprcol_gene1,T,T,F,roi = c(11600,11900,3100,3400))
    proi(dth,'Inflammation',exprcol_gene1,T,T,F,roih,bds=bds_h)
    proi(dth,'Cycling',exprcol_gene1,T,T,F,roih,bds=bds_h)
    proi(dth,'DNA repair',exprcol_gene1,T,T,F,roih,bds=bds_h)
    proi(dth,'Angiogenesis',exprcol_gene1,T,T,F,roi=c(11200,11500,2200,2500))
    proi(dth,'Oncogenic TFs',exprcol_gene1,T,T,F,roih,collim = c(0,2),bds=bds_h)
    proi(dth,'EMT',exprcol_gene1,T,T,F,roih,collim = c(0,0.6),bds=bds_h)
  }
  
  {# escc t4s1
    dte = xe4_2 %>% subset(cells=rownames(t4dtlist$T4S2$E2$crd))
    dte@meta.data = dte@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% 
      left_join(epiprogscore[,c('cellID',unique(prognew$annotation))]) %>% `rownames<-`(.$cellID)
    
    # roit = c(19600,19750,2750,2900)
    # roit = c(20100,20400,2550,2850)
    roit = c(19650,19950,2650,2950)
    
    ## generate contour line 
    pts = data.frame(x=runif(150*150,0,300)+roit[1],y=runif(150*150,0,300)+roit[3])
    poly = t4dtlist$T4S2$E2$poly %>% data.frame() %>% `colnames<-`(c('x','y'))
    plotpoly(list(poly))
    job({pts = pts[isinner2(pts,poly),]},import = 'auto')
    
    xedimplot(dte,'celltype',flip = T,reversey=T,axis=T,xintercept = roit[1:2],yintercept = roit[3:4],zoomin = roit)+
      geom_point(data=pts,size=0.1)
    bd_ = poly %>% filter(x<=roit[2]+500,x>=roit[1]-500,y>=roit[3]-500,y<=roi[4]+500)
    job({ptsdist = apply(pts,1,function(x){disttoseg(x,bd_)})},import = 'auto')
    pts %>% ggplot(aes(x,y))+geom_point(aes(color=ptsdist))+
      theme1+theme(panel.background = element_rect(fill='black'))+
      scale_color_gradient(low = 'white',high = 'red')+
      geom_path(data=bd_,color='grey')
    
    distgroup = max(ptsdist)/10
    bdlist = list()
    for(i in 1:10){
      print(i)
      tmp = pts[ptsdist>distgroup*(i-1),]
      as = alphahull::ashape(tmp,alpha = 20)
      plot(as)
      bdlist[[i]] = extractpoly(as)[[1]] %>% data.frame()
      colnames(bdlist[[i]]) = c('x','y')
      bdlist[[i]]$group = paste0('g',i)
    }
    bds_t = Reduce(rbind,bdlist)
    
    proi(dte,'Oncogenic TFs',exprcol_gene1,T,F,T,roit,collim = c(0,1.5),bds=bds_t)
    proi(dte,'EMT',exprcol_gene1,T,F,T,roit,collim = c(0,0.5),bds=bds_t)
    proi(dte,'Cycling',exprcol_gene1,T,F,T,roit,bds=bds_t)
  }
}





# prog score spatial #####
# {# subset data
#   epnl = xe3_1 %>% subset(cells=intersect(.$cellID[.$ct3=='Epi'],epidt$cellID))
#   eplh = xe3_2 %>% subset(cells=intersect(.$cellID[.$ct3=='Epi'],epidt$cellID))
#   epht = xe1_1 %>% subset(cells=intersect(.$cellID[.$ct3=='Epi'],epidt$cellID))
#   
#   roinl = c(10000,12000,2100,3500)
#   roilh = c(6100,8200,2800,4000)
#   roiht = c(6800,8800,3500,5300)
#   
#   dt1 = subset(epnl, cells = getcellinrange(epnl,roinl,T))
#   dt1 = subset(dt1, cells = dt1$cellID[dt1$stage1!='ESCC'])
#   dt2 = subset(eplh, cells = getcellinrange(eplh,roilh,T))
#   dt3 = subset(epht, cells = getcellinrange(epht,roiht,T))
# }
{# subset data
  t4dist = read_rds('stp3_tme/data/epi_borderdist/T4partE_borderdist_v2.rds')
  t4dtlist = read_rds('stp3_tme/data/epiborder/T4partEborder_v2.rds')
  
  dt1 = xe3_1 %>% subset(cells=intersect(.$cellID[.$ct3=='Epi'],epidt$cellID))
  dt1 = subset(dt1, cells = dt1$cellID[dt1$stage1!='ESCC'])
  # dt1 = xe1_2 %>% subset(cells=intersect(.$cellID[.$ct3=='Epi'],epidt$cellID))
  dt2 = xe1_2 %>% subset(cells=intersect(.$cellID[.$ct3=='Epi'],epidt$cellID))
  # dt3 = dt2
  dt3 = xe4_2 %>% subset(cells=rownames(t4dtlist$T4S2$E2$crd))
  # dt3 = xe4_1 %>% subset(cells=rownames(t4dtlist$T4S1$E1$crd))
  dt1@meta.data = dt1@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% `rownames<-`(.$cellID)
  dt2@meta.data = dt2@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% `rownames<-`(.$cellID)
  dt3@meta.data = dt3@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% `rownames<-`(.$cellID)
  
  roinl = c(10000,12000,2100,3500)
  # roinl = c(19900,20900,3550,4050)
  roih = c(15700,16900,2600,4400)
  # roit = c(2800,4800,1800,3800)
  roit = getcoordsrange(dt3,flip = T)
  xedimplot(dt1,'celltype',flip = T,reversex = T,axis = T,zoomin = roinl)
  xedimplot(dt2,'celltype',flip = T,reversex = T,axis = T,zoomin = roih)
  xedimplot(dt3,'celltype',flip = T,reversex = T,axis = T,zoomin = roit)
  
  dt1 = subset(dt1, cells = getcellinrange(dt1,roinl,T))
  dt2 = subset(dt2, cells = getcellinrange(dt2,roih,T))
  dt3 = subset(dt3, cells = getcellinrange(dt3,roit,T))
  
  bd1 = allborder$T3S1$d %>% cutpoly(roinl)
  bd2 = allborder$T1S2$d %>% cutpoly(roih)
  bd3 = t4dtlist$T4S2$E2$poly
}

{# get score 
  gnlist = prognew %>% 
    filter(gene%in%c('KRT5','KRT15')==FALSE) %>% 
    filter(annotation!='Cellular metabolism')
  gnlist = split(gnlist$gene,gnlist$annotation)
  
  dt1@meta.data = getscore(dt1,gnlist) %>% `rownames<-`(.$cellID)
  dt2@meta.data = getscore(dt2,gnlist) %>% `rownames<-`(.$cellID)
  dt3@meta.data = getscore(dt3,gnlist) %>% `rownames<-`(.$cellID)
}

{# plot v3
  exprcol_gene1 = c('#dfe0fc','white',colorRampPalette(c('white','#E41A1C'))(5)[-1])
  exprcol_gene1 = c('#dfe0fc','white','#E41A1C')
  exprcol_gene1 = c('#dfe0fc','#ececfd','white','#E41A1C','#a50f15')
  # exprcol_gene1 = c('#ccebc5','white','#08519c')
  
  p = lapply(names(gnlist),function(x){
    lim1 = floor(range(dt1@meta.data[dt1$cellID%in%epidt$cellID,x])/0.5)*0.5
    lim2 = floor(range(dt2@meta.data[dt2$cellID%in%epidt$cellID,x])/0.5)*0.5
    lim3 = floor(range(dt3@meta.data[dt3$cellID%in%epidt$cellID,x])/0.5)*0.5
    lims = range(c(lim1,lim2,lim3))
    if(x=='Adhesion'){
      lims = c(0,1.5)
    }
    
    p1=xefeatureplot(dt1,feature = x,bgcol = 'white',flip = T,reversex = T,gridline = F,axis = F,pt.size = 1,color = exprcol_gene1,collim = lims)+theme(plot.background = element_blank(),panel.background = element_blank())
    p2=xefeatureplot(dt2,feature = x,bgcol = 'white',flip = T,reversex = T,gridline = F,axis = F,pt.size = 1,color = exprcol_gene1,collim = lims)+theme(plot.background = element_blank(),panel.background = element_blank())
    p3=xefeatureplot(dt3,feature = x,bgcol = 'white',flip = T,reversey = T,gridline = F,axis = F,pt.size = 1,color = exprcol_gene1,collim = lims)+theme(plot.background = element_blank(),panel.background = element_blank())
    # p = patchwork::wrap_plots(p1,p2,p3,ncol = 3,guides = 'collect')
    p = list(p1,p2,p3)
    return(p)
  })
  p = lapply(1:3,function(x){
    plist = list()
    for(i in 1:8){
      plist[[i]] = p[[i]][[x]]
    }

    prst = cowplot::plot_grid(plotlist = plist,nrow = 2)
  })

  ggsave('stp2_epi/plot/epiProgram_spatial/progscore/exprspatial_nl.pdf',p[[1]],width = 60,height = 10,limitsize = F)
  ggsave('stp2_epi/plot/epiProgram_spatial/progscore/exprspatial_h.pdf',p[[2]],width = 60,height = 10,limitsize = F)
  ggsave('stp2_epi/plot/epiProgram_spatial/progscore/exprspatial_t.pdf',p[[3]],width = 60,height = 10,limitsize = F)
  ggsave('stp2_epi/plot/epiProgram_spatial/progscore/exprspatial_nl.png',p[[1]],width = 60,height = 10,limitsize = F)
  ggsave('stp2_epi/plot/epiProgram_spatial/progscore/exprspatial_h.png',p[[2]],width = 60,height = 10,limitsize = F)
  ggsave('stp2_epi/plot/epiProgram_spatial/progscore/exprspatial_t.png',p[[3]],width = 60,height = 10,limitsize = F)
  ggsave('stp2_epi/plot/epiProgram_spatial/progscore/exprspatial_legend.pdf',p[[3]],width = 60,height = 10,limitsize = F)
  
  
  # p1stg = xedimplot(dt1,groupby = 'stage1',bgcol = 'white',flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = col_stage)
  # p2stg = xedimplot(dt2,groupby = 'stage1',bgcol = 'white',flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = col_stage)
  # p3stg = xedimplot(dt3,groupby = 'stage1',bgcol = 'white',flip = T,reversey = T,gridline = T,axis = F,pt.size = 1,color = col_stage)
  # pstg = patchwork::wrap_plots(p1stg,p2stg,p3stg,ncol = 3,guides = 'collect')
  # 
  # ggsave('stp2_epi/plot/epiProgram_spatial/progscore/allprog_v3.png',cowplot::plot_grid(plotlist = append(p,list(pstg)),ncol = 3),width = 90,height = 20,limitsize = F)
  
  {# plot celltype
    xedimplot(dt1,'celltype',bgcol = 'white',flip = T,reversex = T,pt.size = 1,color = epicol,gridline = T)
    ggsave('stp2_epi/plot/epiProgram_spatial/progscore/epitype_nl.pdf',width = 20,height = 8)
    xedimplot(dt2,'celltype',bgcol = 'white',flip = T,reversex = T,pt.size = 1,color = epicol,gridline = T)
    ggsave('stp2_epi/plot/epiProgram_spatial/progscore/epitype_h.pdf',width = 20,height = 8)
    xedimplot(dt3,'celltype',bgcol = 'white',flip = T,reversey = T,pt.size = 1,color = epicol,gridline = T)
    ggsave('stp2_epi/plot/epiProgram_spatial/progscore/epitype_t.pdf',width = 20,height = 8)
  }
}

{# plot v2
  exprcol_gene1 = c('#dfe0fc','white',colorRampPalette(c('white','#E41A1C'))(5)[-1])
  exprcol_gene1 = c('#dfe0fc','white','#E41A1C')
  
  p = lapply(names(gnlist),function(x){
    lim1 = floor(range(dt1@meta.data[,x])/0.5)*0.5
    lim2 = floor(range(dt2@meta.data[,x])/0.5)*0.5
    lim3 = floor(range(dt3@meta.data[,x])/0.5)*0.5
    lims = range(c(lim1,lim2,lim3))
    if(x=='Adhesion'){
      lims = c(0,1)
    }
    
    p1=xefeatureplot(dt1,feature = x,bgcol = 'white',flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = exprcol_gene1,collim = lims)
    p2=xefeatureplot(dt2,feature = x,bgcol = 'white',flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = exprcol_gene1,collim = lims)
    p3=xefeatureplot(dt3,feature = x,bgcol = 'white',flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = exprcol_gene1,collim = lims)
    p = patchwork::wrap_plots(p1,p2,p3,ncol = 3,guides = 'collect')
    return(p)
  })
  
  p1stg = xedimplot(dt1,groupby = 'stage1',bgcol = 'white',flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = col_stage)
  p2stg = xedimplot(dt2,groupby = 'stage1',bgcol = 'white',flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = col_stage)
  p3stg = xedimplot(dt3,groupby = 'stage1',bgcol = 'white',flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = col_stage)
  pstg = patchwork::wrap_plots(p1stg,p2stg,p3stg,ncol = 3,guides = 'collect')
  
  ggsave('stp2_epi/plot/epiProgram_spatial/progscore/allprog.png',cowplot::plot_grid(plotlist = append(p,list(pstg)),ncol = 3),width = 90,height = 20,limitsize = F)
}

{# plot by stage
  exprcol_gene1 = c(colorRampPalette(c('white','#6baed6','#08519c'))(67),colorRampPalette(c('#c994c7','#ae017e'))(33))
  exprcol_gene1 = c('white','#E41A1C','#312271')
  exprcol_gene1 = c('white','#E41A1C')
  exprcol_gene1 = c('#dfe0fc','white','#E41A1C')
  exprcol_gene1 = c('#dfe0fc','white',colorRampPalette(c('white','#E41A1C'))(5)[-1])
  
  roinl = c(10000,12000,2100,3500)
  roilh = c(6100,8200,2800,4000)
  roiht = c(6800,8800,3500,5300)
  
  p1 = lapply(names(gnlist),function(x){
    exprlim = floor(range(epnl@meta.data[,x])/0.5)*0.5
    xefeatureplot(epnl,feature = x,bgcol = 'white',zoomin = roinl,flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = exprcol_gene1,collim = exprlim)
  })
  p1stg = xedimplot(epnl,groupby = 'stage1',bgcol = 'white',zoomin = roinl,flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = col_stage)
  ggsave('stp2_epi/plot/epiProgram_spatial/progscore/nl.png',cowplot::plot_grid(plotlist = append(p1,list(p1stg)),nrow=3),width = 30,height = 20,limitsize = F)
  
  p2 = lapply(names(gnlist),function(x){
    exprlim = floor(range(eplh@meta.data[,x])/0.5)*0.5
    xefeatureplot(eplh,feature = x,bgcol = 'white',zoomin = roilh,flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = exprcol_gene1,collim = exprlim)
  })
  p2stg = xedimplot(eplh,groupby = 'stage1',bgcol = 'white',zoomin = roilh,flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = col_stage)
  ggsave('stp2_epi/plot/epiProgram_spatial/progscore/lh.png',cowplot::plot_grid(plotlist = append(p2,list(p2stg)),nrow=3),width = 30,height = 20,limitsize = F)
  
  p3 = lapply(names(gnlist),function(x){
    exprlim = floor(range(epht@meta.data[,x])/0.5)*0.5
    xefeatureplot(epht,feature = x,bgcol = 'white',zoomin = roiht,flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = exprcol_gene1,collim = exprlim)
  })
  p3stg = xedimplot(epht,groupby = 'stage1',bgcol = 'white',zoomin = roiht,flip = T,reversex = T,gridline = T,axis = F,pt.size = 1,color = col_stage)
  ggsave('stp2_epi/plot/epiProgram_spatial/progscore/ht.png',cowplot::plot_grid(plotlist = append(p3,list(p3stg)),nrow=3),width = 30,height = 20,limitsize = F)
}




# fov spatial plot #####
## t2 nor fov #####
{# t2 nor fov
  xe2_1 = read_rds('data/xenium_seurat_processed/xenium_t2_s1_seurat.rds')
  xe2_1@meta.data = read_rds('data/xenium_seurat_processed/xenium_t2_s1_meta.rds')
  borderdist = read_rds('stp3_tme/data/epi_borderdist/allepi_borderdist.rds')
  roi = c(5400,6400,1700,2300)
  xedimplot(xe2_1,groupby = 'ct3',color=primcol2,flip = T,axis = T,zoomin = roi,gridline = T)
  
  cid = getcellinrange(xe2_1,roi,flip = T)
  cid = intersect(cid,borderdist$cellID[borderdist$inEpi]) %>% intersect(epidt$cellID)
  
  distdt = borderdist %>% filter(cellID%in%cid)
  distdt$dst = distdt$downDist
  
  # plot spatial  
  srtsp = subset(xe2_1,cells=cid)
  srtsp@meta.data = srtsp@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% left_join(distdt) %>% 
    `rownames<-`(.$cellID)
  
  ## generate contour line 
  pts = data.frame(x=runif(300*180,0,1000)+roi[1],y=runif(300*180,0,600)+roi[3])
  bd1=data.frame(allborder$T2S1$u) %>% `colnames<-`(c('x','y')); bd2=data.frame(allborder$T2S1$d) %>% `colnames<-`(c('x','y'))
  poly = rbind(bd2,bd1[nrow(bd1):1,])
  plotpoly(list(poly))
  job({pts = pts[isinner2(pts,poly),]},import = 'auto')
  
  xedimplot(srtsp,'celltype',flip = T,reversex = T,axis=T,xintercept = roi[1:2],yintercept = roi[3:4],zoomin = roi)+
    geom_point(data=pts,size=0.1)
  bd1_ = bd1 %>% filter(x<=roi[2]+500,x>=roi[1]-500)
  bd2_ = bd2 %>% filter(x<=roi[2]+500,x>=roi[1]-500)
  job({ptsdist1 = apply(pts,1,function(x){disttoseg(x,bd1_)})},import = 'auto')
  job({ptsdist2 = apply(pts,1,function(x){disttoseg(x,bd2_)})},import = 'auto')
  reldist = ptsdist2/(ptsdist1+ptsdist2)
  pts %>% ggplot(aes(x,y))+geom_point(aes(color=reldist))+
    scale_color_gradient(low = 'white',high = 'red')+
    geom_path(data=bd2_)
  
  distgroup = 0.1
  bdlist = list()
  for(i in 1:10){
    tmp = pts[reldist>distgroup*(i-1),]
    as = alphahull::ashape(tmp,alpha = 20)
    print(i)
    # plot(as)
    bdlist[[i]] = extractpoly(as)[[1]] %>% data.frame()
    colnames(bdlist[[i]]) = c('x','y')
    bdlist[[i]]$group = paste0('g',i)
  }
  bds_n = Reduce(rbind,bdlist)
  
  
  xedimplot(srtsp,groupby = 'celltype',flip = T,reversex = T,color=epicol,bgcol = 'white',pt.size = 0.8,gridline = T,gridcol = 'black',zoomin = roi)+
    geom_path(data=bds_n,aes(group=group),color='black',linetype='dashed')
  ggsave('stp2_epi/plot/epiProgram_spatial/fov_spatial/Epispa_NOR.pdf',width = 6,height = 6)
}

## t5 lgin fov #####
{# t5 lgin fov
  xe5_2_split = read_rds('data/xenium_seurat_processed/xenium_t5_s2_split_seurat.rds')
  xe5_2_split@meta.data = read_rds('data/xenium_seurat_processed/xenium_t5_s2_split_meta.rds')
  roi = c(2500,3500,2930,3530)
  xedimplot(xe5_2_split,groupby = 'ct3',color=primcol2,axis = T,zoomin = roi)
  
  cid = getcellinrange(xe5_2_split,roi)
  cid = intersect(cid,borderdist$cellID) %>% intersect(epidt$cellID)
  
  distdt = borderdist %>% filter(cellID%in%cid)
  distdt$dst = distdt$downDist
  
  # plot spatial  
  srtsp = subset(xe5_2_split,cells=cid)
  srtsp@meta.data = srtsp@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% left_join(distdt) %>% 
    `rownames<-`(.$cellID)
  xedimplot(srtsp,groupby = 'ct3',color=primcol2,axis = T,zoomin = roi)
  
  ## generate contour line 
  pts = data.frame(x=runif(300*180,0,1000)+roi[1],y=runif(300*180,0,600)+roi[3])
  bd1=data.frame(allborder$T5S2$u) %>% `colnames<-`(c('x','y')); bd2=data.frame(allborder$T5S2$d) %>% `colnames<-`(c('x','y'))
  poly = rbind(bd2,bd1)
  plotpoly(list(poly))
  job({pts = pts[isinner2(pts,poly),]},import = 'auto')
  
  xedimplot(srtsp,'celltype',axis=T,xintercept = roi[1:2],yintercept = roi[3:4],zoomin = roi)+
    geom_point(data=pts,size=0.1)
  bd1_ = bd1 %>% filter(x<=roi[2]+500,x>=roi[1]-500)
  bd2_ = bd2 %>% filter(x<=roi[2]+500,x>=roi[1]-500)
  job({ptsdist1 = apply(pts,1,function(x){disttoseg(x,bd1_)})},import = 'auto')
  job({ptsdist2 = apply(pts,1,function(x){disttoseg(x,bd2_)})},import = 'auto')
  reldist = ptsdist2/(ptsdist1+ptsdist2)
  pts %>% ggplot(aes(x,y))+geom_point(aes(color=reldist))+
    scale_color_gradient(low = 'white',high = 'red')+
    geom_path(data=bd2_)
  
  distgroup = 0.1
  bdlist = list()
  for(i in 1:10){
    tmp = pts[reldist>distgroup*(i-1),]
    as = alphahull::ashape(tmp,alpha = 25)
    print(i)
    # plot(as)
    bdlist[[i]] = extractpoly(as)[[1]] %>% data.frame()
    colnames(bdlist[[i]]) = c('x','y')
    bdlist[[i]]$group = paste0('g',i)
  }
  bds = Reduce(rbind,bdlist)
  
  xedimplot(srtsp,groupby = 'celltype',color=epicol,bgcol = 'white',pt.size = 0.8,gridline = T,gridcol = 'black')+
    geom_path(data=bds,aes(group=group),color='black',linetype='dashed')
  ggsave('stp2_epi/plot/epiProgram_spatial/fov_spatial/Epispa_LGIN.pdf',width = 6,height = 6)
}

## t1 hgin fov #####
{# hgin t1s3
  roih = c(10300,12600,2100,3500)
  dth = subset(xe1_3, cells=intersect(getcellinrange(xe1_3,roih,T),t1dist$cellID[t1dist$inEpi]))
  dth = subset(dth,cells=epidt$cellID)
  dth@meta.data = dth@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% 
    left_join(epiprogscore[,c('cellID',progorder)]) %>% `rownames<-`(.$cellID)
  xedimplot(dth,'celltype',flip = T,reversex = T,axis=T,gridline = T,xintercept = roih[1:2],yintercept = roih[3:4],zoomin = roih)
  
  ## generate contour line 
  pts = data.frame(x=runif(300*300,0,2300)+roih[1],y=runif(300*300,0,1400)+roih[3])
  bd1=data.frame(t1border$T1S3$u) %>% `colnames<-`(c('x','y')); bd2=data.frame(t1border$T1S3$d) %>% `colnames<-`(c('x','y'))
  poly = rbind(bd2,bd1[nrow(t1border$T1S3$u):1,])
  plotpoly(list(poly))
  job({pts = pts[isinner2(pts,poly),]},import = 'auto')
  
  xedimplot(dth,'celltype',flip = T,reversex = T,axis=T,xintercept = roih[1:2],yintercept = roih[3:4],zoomin = roih)+
    geom_point(data=pts,size=0.1)
  bd1_ = bd1 %>% filter(x<=roih[2]+500,x>=roih[1]-500)
  bd2_ = bd2 %>% filter(x<=roih[2]+500,x>=roih[1]-500)
  job({ptsdist1 = apply(pts,1,function(x){disttoseg(x,bd1_)})},import = 'auto')
  job({ptsdist2 = apply(pts,1,function(x){disttoseg(x,bd2_)})},import = 'auto')
  reldist = ptsdist2/(ptsdist1+ptsdist2)
  pts %>% ggplot(aes(x,y))+geom_point(aes(color=reldist))+
    scale_color_gradient(low = 'white',high = 'red')+
    geom_path(data=bd2_)
  
  distgroup = 0.1
  bdlist = list()
  for(i in 1:10){
    tmp = pts[reldist>distgroup*(i-1),]
    as = alphahull::ashape(tmp,alpha = 20)
    print(i)
    # plot(as)
    bdlist[[i]] = extractpoly(as)[[1]] %>% data.frame()
    colnames(bdlist[[i]]) = c('x','y')
    bdlist[[i]]$group = paste0('g',i)
  }
  bds_h = Reduce(rbind,bdlist)
  
  xedimplot(dth,groupby = 'celltype',flip = T,reversex = T,color=epicol,bgcol = 'white',pt.size = 0.8,gridline = T,gridcol = 'black')+
    geom_path(data=bds_h,aes(group=group),color='black',linetype='dashed')
  ggsave('stp2_epi/plot/epiProgram_spatial/fov_spatial/Epispa_HGIN.pdf',width = 12,height = 6)
}

## t4 escc fov #####
{# escc t4s1
  t4dist = read_rds('stp3_tme/data/epi_borderdist/T4partE_borderdist_v2.rds')
  t4dtlist = read_rds('stp3_tme/data/epiborder/T4partEborder_v2.rds')
  
  dte = xe4_2 %>% subset(cells=rownames(t4dtlist$T4S2$E2$crd))
  dte@meta.data = dte@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% 
    left_join(epiprogscore[,c('cellID',unique(prognew$annotation))]) %>% `rownames<-`(.$cellID)
  xedimplot(dte,'celltype',flip = T,reversey = T,axis=T,gridline = T)
  roit = getcoordsrange(dte,flip = T)

  ## generate contour line 
  pts = data.frame(x=runif(300*300,0,2000)+roit[1],y=runif(300*300,0,1500)+roit[3])
  poly = t4dtlist$T4S2$E2$poly %>% data.frame() %>% `colnames<-`(c('x','y'))
  plotpoly(list(poly))
  job({pts = pts[isinner2(pts,poly),]},import = 'auto')
  
  xedimplot(dte,'celltype',flip = T,reversey=T,axis=T,xintercept = roit[1:2],yintercept = roit[3:4],zoomin = roit)+
    geom_point(data=pts,size=0.1)
  bd_ = poly
  job({ptsdist = apply(pts,1,function(x){disttoseg(x,bd_)})},import = 'auto')
  pts %>% ggplot(aes(x,y))+geom_point(aes(color=ptsdist))+
    theme1+theme(panel.background = element_rect(fill='black'))+
    scale_color_gradient(low = 'white',high = 'red')+
    geom_path(data=bd_,color='grey')
  
  distgroup = max(ptsdist)/10
  bdlist = list()
  for(i in 1:10){
    tmp = pts[ptsdist>distgroup*(i-1),]
    as = alphahull::ashape(tmp,alpha = 20)
    print(i)
    # plot(as)
    bdlist[[i]] = extractpoly(as)[[1]] %>% data.frame()
    colnames(bdlist[[i]]) = c('x','y')
    bdlist[[i]]$group = paste0('g',i)
  }
  bds_t = Reduce(rbind,bdlist)
  
  xedimplot(dte,groupby = 'celltype',flip = T,reversey = T,color=epicol,bgcol = 'white',pt.size = 0.8,gridline = T,gridcol = 'black')+
    geom_path(data=bds_t,aes(group=group),color='black',linetype='dashed')
  ggsave('stp2_epi/plot/epiProgram_spatial/fov_spatial/Epispa_ESCC.pdf',width = 16,height = 8)
}









