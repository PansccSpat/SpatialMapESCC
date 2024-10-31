source('~/xenium/header.R')
source('~/xenium/xenium_function.R')
library(ComplexHeatmap)

##### input #####
srtepi = read_rds('data/xenium_seurat_processed/ct_epi_SCT_2.rds')
epimeta = read_rds('data/xenium_seurat_processed/ct_epi_meta_2.rds')
epidt = read_rds('stp2_epi/data/epiid_correctstage.rds')

# lamexpr = fread('stp2_epi/data/lamfish_epi12_jagnotchexpr.txt')
lam = read_rds('stp2_epi/data/lam_epi12_seurat.rds')
lam@meta.data = read_rds('stp2_epi/data/lam_epi12_meta.rds')

epianno = fread('stp1_summary/data/xenium_epi.csv',data.table = F)

epigeneannot = fread('stp3_tme/data/hm_gene_anno.txt',data.table = F)

probelist = fread('data/probe_list.txt',data.table = F)

gridid = read_rds('stp1_summary/data/allsample_gridid.rds')

epiblacklist = fread('stp2_epi/data/gene_otherTissue.txt',data.table = F)
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
epicol1 = epicol 
names(epicol1) = names(epicol1) %>% mapvalues(c('Proliferation','Differentiation'),c('Proliferative','Differentiated'))

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


mp_malignant = fread('~/metaprogram_data/program_malignant.csv',data.table = F)
mp_epi = fread('~/metaprogram_data/program_epi.csv',data.table = F)
mpall = mp_epi %>% gather('prog','gene') %>% mutate(ct='epi')
mpall = mp_malignant %>% gather('prog','gene') %>% mutate(ct='malignant') %>% rbind(mpall)

supp2 = fread('data/supptable/supp2_0422_reshaped.txt',data.table = F)
supp2_epig = supp2 %>% filter(inxe=='Yes')
supp2_epig = supp2_epig[grep('^Marker|Other',supp2_epig$ann_new,invert = T),]

prognew = read_rds('stp2_epi/data/epiProgram/epigene_annotation.rds')





##### progscore statplot v3 #####
prognew = read_rds('stp2_epi/data/epiProgram/epigene_annotation_v5.rds')
# epiprogscore = getscore(srtepi,split(prognew$gene,prognew$annotation))
# saveRDS(epiprogscore,'stp2_epi/data/epiProgram/epiprogscore_v5.rds')
epiprogscore = read_rds('stp2_epi/data/epiProgram/epiprogscore_v5.rds')
progorder = c('Cell adhesion','Mucosal defense','Inflammation','Cycling','DNA repair','Angiogenesis','Oncogenic TFs','Metastasis')
progorder = c('Cell adhesion','Mucosal defense','Inflammation','Cycling','DNA repair','Angiogenesis','Oncogenic TFs','EMT')

epiprogscore$stage1 = epiprogscore$stage1 %>% factor(levels = c('NOR','LGIN','HGIN','ESCC'))
epiprogscore_aggr = epiprogscore %>% filter(cellID%in%epidt$cellID) %>% select(all_of(unique(prognew$annotation)), stage1, sp_stg, celltype) %>% 
  gather('prog','score',all_of(unique(prognew$annotation)))
epiprogscore_aggr = aggregate(list(score=epiprogscore_aggr$score),by=as.list(epiprogscore_aggr[,c('stage1','sp_stg','prog')]), mean)
epiprogscore_aggr$prog = epiprogscore_aggr$prog %>% factor(levels = progorder)

# progcol = c("#8DD3C7","#41ab5d","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5")
# names(progcol) = c('Cell adhesion','Cycling','Mucosal defense','DNA repair','Oncogenic TFs','Inflammation','Angiogenesis','Metastasis')

epiprogscore_aggr %>% 
  ggplot(aes(stage1,score))+
  theme1+scale_color_manual(values = progcol)+scale_fill_manual(values = progcol)+
  RotatedAxis()+xlab('')+ylab('')+NoLegend()+
  facet_wrap(~prog,scales = 'free',nrow = 2)+
  scale_y_continuous(expand = expansion(c(0.05,0.1)))+
  stat_compare_means(comparisons = makecomparison(4))+
  ggbeeswarm::geom_quasirandom(aes(color=prog),shape=16)+
  geom_boxplot(aes(fill=prog),alpha=0.3,outliers = F)
ggsave('stp2_epi/plot/epiProgram/progscore_stage_all.pdf',width = 10,height = 7)





##### progscore bubble #####
prognew = read_rds('stp2_epi/data/epiProgram/epigene_annotation_v5.rds')
progorder = c('Cell adhesion','Mucosal defense','Inflammation','Cycling','DNA repair','Angiogenesis','Oncogenic TFs','Metastasis')
progorder = c('Cell adhesion','Mucosal defense','Inflammation','Cycling','DNA repair','Angiogenesis','Oncogenic TFs','EMT')

epiprogscore = read_rds('stp2_epi/data/epiProgram/epiprogscore_v5.rds')
epiprogscore$stage1 = epiprogscore$stage1 %>% factor(levels = c('NOR','LGIN','HGIN','ESCC'))

## aggr score
tmp = epiprogscore %>% filter(cellID%in%epidt$cellID) %>% select(all_of(unique(prognew$annotation)), stage1, sp_stg, celltype) %>% 
  gather('prog','score',all_of(unique(prognew$annotation)))
# tmp1 = aggregate(list(score=tmp$score),by=as.list(tmp[,c('stage1','celltype','prog')]), mean)
tmp1 = aggregate(list(score=tmp$score),by=as.list(tmp[,c('celltype','prog')]), mean)

## calc dot size
{
  # tmp2 = aggregate(list(thresh=tmp$score),by=as.list(tmp[,c('prog','stage1'),drop=F]),function(x) mean(x)+0.5*sd(x))
  # tmp2 = left_join(tmp,tmp2) %>% mutate(highexpr=ifelse(score>thresh,1,0))
  # tmp2 = aggregate(list(prop=tmp2$highexpr),by=as.list(tmp2[,c('stage1','celltype','prog')]),mean)
  
  tmp2 = aggregate(list(thresh=tmp$score),by=as.list(tmp[,c('prog'),drop=F]),function(x) mean(x)+0.5*sd(x))
  tmp2 = left_join(tmp,tmp2) %>% mutate(highexpr=ifelse(score>thresh,1,0))
  tmp2 = aggregate(list(prop=tmp2$highexpr),by=as.list(tmp2[,c('celltype','prog')]),mean)
  
  # tmp2 = getprop(epidt$celltype,epidt[,c('stage1'),drop=F]) %>% gather('celltype','prop',unique(epidt$celltype)) %>% `colnames<-`(c('stage1','celltype','prop'))
}

## plot bubble
{
  tmp1_scale = tmp1 %>% spread(prog,score)
  tmp1_scale[,unique(prognew$annotation)] = tmp1_scale[,unique(prognew$annotation)] %>% scale
  tmp1_scale = tmp1_scale %>% gather('prog','score_scale',unique(prognew$annotation))
  tmp1_scale = tmp1_scale %>% left_join(tmp2)
  tmp1_scale$prog = tmp1_scale$prog %>% factor(levels = rev(progorder)) %>% 
    mapvalues(c('Adhesion','Cell cycle','Hypoxia & Inflammation'),c('Cell adhesion','Cycling','Inflammation'))
  tmp1_scale$celltype = tmp1_scale$celltype %>% factor(levels = c('Basal','Proliferation','Differentiation','Invasive')) %>% 
    mapvalues(c('Proliferation','Differentiation'),c('Proliferative','Differentiated'))
  
  
  bubcol = c(colorRampPalette(c('#440154','#482576'))(40),colorRampPalette(c('#482576','#21908c','#bbdf27'))(20),colorRampPalette(c('#bbdf27','#fde725'))(40))
  bubcol = c(rep('#440154',20),colorRampPalette(c('#440154','#21908c','#fde725'))(60),rep('#fde725',20))
  
  # tmp1_scale %>% 
  #   mutate(prop = ifelse(prop>0.5,0.5,prop)) %>%
  #   ggplot(aes(celltype,prog))+
  #   theme1+RotatedAxis()+xlab('')+ylab('')+
  #   # theme(panel.spacing.x = unit(0.2,'cm'))+
  #   scale_fill_gradientn(name = 'Z-score',colors = bubcol)+
  #   scale_size_continuous(name = 'Proportion',range = c(1,8),labels = function(x){return(paste0(x*100,'%'))})+
  #   geom_point(shape=21,aes(fill=score_scale,size=prop))+
  #   facet_wrap(~stage1,scales = 'free_x',nrow = 1)
  # ggsave('stp2_epi/plot/epiProgram/progscore_bubble.pdf',width = 8,height = 5)
  # 
  # tmp1_scale %>% 
  #   mutate(prop = ifelse(prop>0.5,0.5,prop)) %>%
  #   ggplot(aes(stage1,prog))+
  #   theme1+RotatedAxis()+xlab('')+ylab('')+
  #   # theme(panel.spacing.x = unit(0.2,'cm'))+
  #   scale_fill_gradientn(name = 'Z-score',colors = bubcol)+
  #   scale_size_continuous(name = 'Proportion',range = c(1,8),labels = function(x){return(paste0(x*100,'%'))})+
  #   geom_point(shape=21,aes(fill=score_scale,size=prop))+
  #   facet_wrap(~celltype,scales = 'free_x',nrow = 1)
  # ggsave('stp2_epi/plot/epiProgram/progscore_bubble.pdf',width = 8,height = 5)
  
  tmp1_scale %>% 
    mutate(prop = ifelse(prop>0.5,0.5,prop)) %>%
    ggplot(aes(celltype,prog))+
    theme1+RotatedAxis()+xlab('')+ylab('')+
    # theme(panel.spacing.x = unit(0.2,'cm'))+
    scale_fill_gradientn(name = 'Z-score',colors = bubcol)+
    scale_size_continuous(name = 'Proportion',range = c(1,8),labels = function(x){return(paste0(x*100,'%'))})+
    geom_point(shape=21,aes(fill=score_scale,size=prop))
  ggsave('stp2_epi/plot/epiProgram/progscore_bubble.pdf',width = 5,height = 5)
}






##### gene heatmap new #####
prognew = read_rds('stp2_epi/data/epiProgram/epigene_annotation_v5.rds')
epiprogscore = read_rds('stp2_epi/data/epiProgram/epiprogscore_v5.rds')
rownames(epiprogscore) = epiprogscore$cellID

progorder = unique(prognew$annotation)
progorder = c('Mucosal defense','Cell death','Metastasis','Hypoxia','Inflammation','DNA damage repair',
              'Angiogenesis','Transcription regulation','Cell cycle','Stemness')
progorder = c('Basal','Cell cycle','Mucosal defense','P53 related','Angiogenesis & Metastasis','Hypoxia & Inflammation',
              'Transcription regulation','Cell death')
progorder = c('Adhesion','Cell cycle','Mucosal defense','P53 pathway','Oncogenic TFs','Hypoxia & Inflammation','Angiogenesis','Metastasis')
progorder = c('Cell adhesion','Mucosal defense','Inflammation','Cycling','DNA repair','Angiogenesis','Oncogenic TFs','Metastasis')
progorder = c('Cell adhesion','Mucosal defense','Inflammation','Cycling','DNA repair','Angiogenesis','Oncogenic TFs','EMT')
progorder = c('Cell adhesion','Mucosal defense','Cycling','Inflammation','DNA repair','EMT','Oncogenic TFs','Angiogenesis')

geneord2 = prognew %>% filter(annotation%in%progorder)
geneord2$annotation = geneord2$annotation %>% factor(levels = progorder)


{# cell order
  ordcell2 = function(cid,progname){
    score1 = epiprogscore[cid,progname,drop=F]
    cellord = apply(score1,1,function(x) return(progname[which.max(x)]))
    cellord = data.frame(cellID=cid,prog=cellord,score=matrixStats::rowMaxs(as.matrix(score1)))
    
    cellord$prog = cellord$prog %>% factor(levels = progorder)
    cellord = cellord %>% arrange(prog,score)
    
    return(cellord)
  }
  
  cellsubset = lapply(levels(epidt$stage1),function(x){
    dt1 = epidt[epidt$stage1==x,]
    dt1$celltype = dt1$celltype %>% factor %>% droplevels()
    cellid = lapply(levels(dt1$celltype),function(y){
      cellid1 = dt1$cellID[dt1$celltype==y]
      return(sample(cellid1,1e4))
    }) %>% unlist()
    return(cellid)
  }) %>% unlist
  
  metadt = srtepi@meta.data[cellsubset,]
  metadt$celltype = metadt$celltype %>% factor(levels = c('Basal','Proliferation','Invasive','Differentiation'))
  
  # sort by stage, celltype 
  # cellord = lapply(levels(metadt$stage1),function(x){
  #   dt1 = metadt %>% filter(stage1==x)
  #   dt1$celltype = dt1$celltype %>% droplevels()
  #   cellord1 = lapply(levels(dt1$celltype),function(y){
  #     dt2 = dt1 %>% filter(celltype==y)
  #     pname = switch(y,
  #                    'Differentiation'='Mucosal defense',
  #                    'Basal'='Adhesion',
  #                    'Proliferation'='Cell cycle',
  #                    'Invasive'=setdiff(progorder,c('Mucosal defense','Adhesion','Cell cycle')))
  #     cellord2 = ordcell2(dt2$cellID,pname)
  #     
  #     cellord2 = data.frame(cellID=dt2$cellID,prog=pname[[1]])
  #     return(cellord2)
  #   }) %>% Reduce(rbind,.)
  # }) %>% Reduce(rbind,.)
  
  # sort by celltype, stage 
  cellord = lapply(levels(metadt$celltype),function(x){
    dt1 = metadt %>% filter(celltype==x)
    dt1$stage1 = dt1$stage1 %>% droplevels()
    pname = switch(x,
                   'Differentiation'='Mucosal defense',
                   'Basal'='Adhesion',
                   'Proliferation'='Cell cycle',
                   'Invasive'=setdiff(progorder,c('Mucosal defense','Adhesion','Cell cycle')))
    cellord1 = lapply(levels(dt1$stage1),function(y){
      dt2 = dt1 %>% filter(stage1==y)
      # cellord2 = ordcell2(dt2$cellID,pname)
      cellord2 = data.frame(cellID=dt2$cellID,prog=pname[[1]])
      return(cellord2)
    }) %>% Reduce(rbind,.)
  }) %>% Reduce(rbind,.)
  
  metadt = metadt %>% inner_join(cellord) %>% `rownames<-`(.$cellID) %>% .[cellord$cellID,]
  exprdt = srtepi@assays$SCT@data[geneord2$gene,metadt$cellID] %>% as.matrix()
}

{# gene order
  exprscore = exprdt %>% t %>% scale %>% t
  exprscore[exprscore<0] = 0
  exprscore = rowSums(exprscore) %>% data.frame(gene=names(.),score=.)
  geneord3 = geneord2 %>% left_join(exprscore)
  geneord3 = geneord3 %>% arrange(annotation,desc(score))
  saveRDS(geneord3,'stp2_epi/data/epiProgram/epigene_heatorder.rds')
}

{# plot 
  exprdt = exprdt[geneord3$gene,metadt$cellID]
  
  ## plot together
  cellann = columnAnnotation(df=metadt[,c('stage1','celltype'),drop=F],
                             col=list(stage1=col_stage,celltype=epicol))
  geneann = rowAnnotation(df=geneord3[,c('annotation'),drop=F],
                          col=list(annotation=progcol))
  
  exprdt = exprdt %>% t %>% scale %>% t
  tmp = ComplexHeatmap::Heatmap(exprdt,cluster_rows = F,cluster_columns = F,show_column_names = F,
                                # row_split = geneord3$annotation,
                                # column_split = metadt$celltype,
                                row_gap = unit(1,'mm'),
                                column_gap = unit(1,'mm'),
                                top_annotation = cellann,
                                left_annotation = geneann,
                                heatmap_legend_param = list(title='Z-score',direction='vertical',title_position='topleft',legend_height=unit(3,'cm')),
                                col = circlize::colorRamp2(c(-1.5,0,1.5),c('purple','black','yellow')))
  
  pdf('stp2_epi/plot/epiProgram/progheatmap/progene_heatmap_subset2_ctstage.pdf',width = 15,height = 20)
  ComplexHeatmap::draw(tmp,heatmap_legend_side='right')
  dev.off()
}





##### gene heatmap order by program #####
prognew = read_rds('stp2_epi/data/epiProgram/epigene_annotation_v3.rds')
epiprogscore = read_rds('stp2_epi/data/epiProgram/epiprogscore_v3.rds')
rownames(epiprogscore) = epiprogscore$cellID

progorder = unique(prognew$annotation)
progorder = c('Mucosal defense','Cell death','Metastasis','Hypoxia','Inflammation','DNA damage repair',
              'Angiogenesis','Transcription regulation','Cell cycle','Stemness')
progorder = c('Basal','Cell cycle','Mucosal defense','P53 related','Angiogenesis & Metastasis','Hypoxia & Inflammation',
              'Transcription regulation','Cell death')

geneord2 = prognew %>% filter(annotation%in%progorder)
geneord2$annotation = geneord2$annotation %>% factor(levels = progorder)


{# cell order
  ordcell = function(cid){
    cellord = lapply(cid,function(x){
      dt1 = epiprogscore[x,progorder] %>% unlist()
      rst = data.frame(cellID = x, prog=progorder[which.max(dt1)],score=dt1[which.max(dt1)])
      return(rst)
    }) %>% Reduce(rbind,.)
    
    cellord$prog = cellord$prog %>% factor(levels = rev(progorder))
    cellord = cellord %>% arrange(prog,score)
    
    return(cellord)
  }
  ordcell2 = function(cid,progname){
    score1 = epiprogscore[cid,progname,drop=F]
    cellord = apply(score1,1,function(x) return(progname[which.max(x)]))
    cellord = data.frame(cellID=cid,prog=cellord,score=matrixStats::rowMaxs(as.matrix(score1)))
    
    cellord$prog = cellord$prog %>% factor(levels = rev(progorder))
    cellord = cellord %>% arrange(prog,score)
    
    return(cellord)
  }
  
  cellsubset = lapply(unique(epidt$stage1),function(x){
    cellid = epidt$cellID[epidt$stage1==x]
    return(sample(cellid,1e5))
  }) %>% unlist
  # cellsubset = epidt$cellID
  
  metadt = srtepi@meta.data[cellsubset,]
  metadt$celltype = metadt$celltype %>% factor(levels = c('Basal','Proliferation','Invasive','Differentiation'))
  
  # sort by stage
  cellord = lapply(levels(metadt$stage1),function(x){
    dt1 = metadt %>% filter(stage1==x)

    cellord1 = ordcell(dt1$cellID)
    return(cellord1)
  }) %>% Reduce(rbind,.)
  
  # sort by stage, celltype 
  cellord = lapply(levels(metadt$stage1),function(x){
    dt1 = metadt %>% filter(stage1==x)
    dt1$celltype = dt1$celltype %>% droplevels()
    cellord1 = lapply(levels(dt1$celltype),function(y){
      dt2 = dt1 %>% filter(celltype==y)
      pname = switch(y,
                     'Differentiation'='Mucosal defense',
                     'Basal'='Stemness',
                     'Proliferation'='Cell cycle',
                     'Invasive'=setdiff(progorder,c('Mucosal defense','Stemness','Cell cycle')))
      cellord2 = ordcell2(dt2$cellID,pname)
      return(cellord2)
    }) %>% Reduce(rbind,.)
  }) %>% Reduce(rbind,.)
  
  # sort by celltype, stage 
  cellord = lapply(levels(metadt$celltype),function(x){
    dt1 = metadt %>% filter(celltype==x)
    dt1$stage1 = dt1$stage1 %>% droplevels()
    pname = switch(x,
                   'Differentiation'='Mucosal defense',
                   'Basal'='Stemness',
                   'Proliferation'='Cell cycle',
                   'Invasive'=setdiff(progorder,c('Mucosal defense','Stemness','Cell cycle')))
    cellord1 = lapply(levels(dt1$stage1),function(y){
      dt2 = dt1 %>% filter(stage1==y)
      cellord2 = ordcell2(dt2$cellID,pname)
      return(cellord2)
    }) %>% Reduce(rbind,.)
  }) %>% Reduce(rbind,.)
  
  metadt = metadt %>% inner_join(cellord) %>% `rownames<-`(.$cellID) %>% .[cellord$cellID,]
  exprdt = srtepi@assays$SCT@data[geneord2$gene,metadt$cellID] %>% as.matrix()
  
  
  ## plot together
  cellann = columnAnnotation(df=metadt[,c('prog','celltype','stage1'),drop=F],
                             col=list(prog=setNames(brewer.pal(10,'Set3'),levels(geneord2$annotation)),celltype=epicol,stage1=col_stage))
  # cellann = columnAnnotation(df=metadt[,c('stage1'),drop=F],
  #                            col=list(stage1=col_stage))
  exprdt = exprdt %>% t %>% scale %>% t
  # exprdt = exprdt %>% apply(1,function(x){rescale(x,c(-2,2))}) %>% t
  tmp = ComplexHeatmap::Heatmap(exprdt,cluster_rows = F,cluster_columns = F,show_column_names = F,
                                row_split = geneord2$annotation,
                                column_split = metadt$celltype,
                                row_gap = unit(2,'mm'),
                                column_gap = unit(2,'mm'),
                                top_annotation = cellann,
                                heatmap_legend_param = list(title='Z-score',direction='vertical',title_position='topleft',legend_height=unit(3,'cm')),
                                col = circlize::colorRamp2(c(-2,0,2),c('#377eb8','#f0f0f0','#e41a1c')))
  
  pdf('stp2_epi/plot/epiProgram/progheatmap/progene_heatmap_subset_ctstage.pdf',width = 20,height = 20)
  ComplexHeatmap::draw(tmp,heatmap_legend_side='right')
  dev.off()
  
  
  ## plot per stage 
  lapply(unique(metadt$stage1),function(x){
    metadt = metadt %>% filter(stage1==x)
    exprdt = exprdt[,metadt$cellID]
    
    cellann = columnAnnotation(df=metadt[,c('prog','celltype','stage1'),drop=F],
                               col=list(prog=setNames(brewer.pal(10,'Set3'),levels(geneord2$annotation)),celltype=epicol,stage1=col_stage))
    exprdt = exprdt %>% t %>% scale %>% t
    tmp = ComplexHeatmap::Heatmap(exprdt,cluster_rows = F,cluster_columns = F,show_column_names = F,
                                  row_split = geneord2$annotation,
                                  column_split = metadt$stage1,
                                  top_annotation = cellann,
                                  heatmap_legend_param = list(title='Z-score',direction='vertical',title_position='topleft',legend_height=unit(3,'cm')),
                                  col = circlize::colorRamp2(c(-2,0,2),c('#377eb8','#f0f0f0','#e41a1c')))
    
    fname = paste0('stp2_epi/plot/epiProgram/progene_heatmap_',x,'.pdf')
    pdf(fname,width = 20,height = 20)
    ComplexHeatmap::draw(tmp,heatmap_legend_side='right')
    dev.off()
  })
}




