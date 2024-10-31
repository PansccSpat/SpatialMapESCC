source('~/xenium/header.R')
source('~/xenium/xenium_function.R')
library(ggstar)

##### input #####
srtepi = read_rds('data/xenium_seurat_processed/ct_epi_srt_2.rds')

# lamexpr = fread('stp2_epi/data/lamfish_epi12_jagnotchexpr.txt')
lam = read_rds('stp2_epi/data/lam_epi12_seurat.rds')
lam@meta.data = read_rds('stp2_epi/data/lam_epi12_meta.rds')

epidt = read_rds('stp2_epi/data/epiid_correctstage.rds')

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




##### celltype spatial #####
{# balck bg
  objlist = setNames(list('xe1_1','xe1_2','xe1_3','xe2_1','xe2_2','xe2_3','xe3_1','xe3_2',
                          'xe4_1','xe4_2','xe5_1','xe5_2_split','xe5_3_split','xe6_1','xe6_2','xe6_3'),
                     c('T1S1','T1S2','T1S3','T2S1','T2S2','T2S3','T3S1','T3S2',
                       'T4S1','T4S2','T5S1','T5S2','T5S3','T6S1','T6S2','T6S3'))
  lapply(names(objlist), function(x){
    dt = get(objlist[[x]])
    dt = subset(dt,cells=intersect(dt$cellID,epidt$cellID))
    dt@meta.data = dt@meta.data %>% left_join(epidt[,c('cellID','celltype')]) %>% `rownames<-`(.$cellID)
    
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
    
    figsize = getcoordsrange(dt,dt$cellID,flip)
    wd = (figsize[2]-figsize[1])/500
    ht = (figsize[4]-figsize[3])/300
    fname = paste0('stp2_epi/plot_epitype/',x,'_celltype.png')
    p = xedimplot(dt, groupby = 'celltype',flip = flip,reversex = reversex,reversey = reversey,color = epicol,bgcol = 'black',pt.size = 0.1,gridline = T,gridcol = 'white')+NoLegend()
    ggsave(fname,p,width = 40, height = 15, limitsize = F)
  })
}






##### gene expr spatial #####
{# molecular index
  subsetmolecular = function(dt, genes, cells=NULL)
  {
    if(is.null(cells)){
      cells=colnames(dt)
    }
    cells=intersect(colnames(dt),cells)
    
    cellpoly = dt@images$fov@boundaries$segmentation@polygons
    cellpoly = lapply(cellpoly[cells],function(x){
      pl = x@Polygons[[1]]@coords %>% data.frame()
      dup = paste0(pl[,1],'-',pl[,2])
      pl = pl[!duplicated(dup),]
      return(pl)
    })
    cellpolyall = cellpoly %>% rbindlist() %>% data.frame()
    molrange = c(min(cellpolyall$x),max(cellpolyall$x),min(cellpolyall$y),max(cellpolyall$y))
    
    mols = dt@images$fov@molecules$molecules[genes]
    mols = lapply(mols,function(x){
      crd = x@coords %>% data.frame()
      crd = crd %>% filter(x>=molrange[1],x<=molrange[2],y>=molrange[3],y<=molrange[4])
      inrange1 = lapply(cellpoly,function(y){
        rst = isinner2(crd,y)
      }) %>% Reduce('|',.)
      rst = crd[inrange1,]
    })
    
    return(mols)
  }
  
  dtlist = list(nl=list(dt=epnl,roi=roinl),lh=list(dt=eplh,roi=roilh),ht=list(dt=epht,roi=roiht))
  job({
    isinner2
    mols = lapply(dtlist,function(x){
      rst = subsetmolecular(x$dt,cycgene,getcellinrange(x$dt,x$roi,T))
      return(rst)
    })
    saveRDS(mols,'stp2_epi/data/3roi_13genes_epimolecularcoords_v2.rds')
  },import='auto')
  
  mols = read_rds('stp2_epi/data/3roi_13genes_epimolecularcoords_v2.rds')
  epi_cycmol = lapply(cycgene,function(x){
    dt = lapply(names(mols),function(y){
      dt1 = mols[[y]][[x]]
      dt1$gene=x
      dt1$stage=y
      return(dt1)
    }) %>% Reduce(rbind,.)
    return(dt)
  }) %>% Reduce(rbind,.)
  epi_cycmol = epi_cycmol %>% dplyr::rename(c('x'='y','y'='x'))
  epi_cycmol$gene = epi_cycmol$gene %>% factor(levels = rev(c('JAG1','NOTCH1','EGFR','ERBB2','MDM2','MYC','STMN1','AURKA','CCND1','SPP1','VEGFA','SOX2','TP63')))
  epi_cycmol$stage = epi_cycmol$stage %>% factor(levels = c('nl','lh','ht'))
  p=epi_cycmol %>% 
    arrange(gene) %>% 
    # filter(gene%in%c('JAG1','NOTCH1','EGFR','ERBB2','MDM2','MYC','STMN1','AURKA','CCND1')) %>% 
    # filter(gene%in%c('SPP1','VEGFA','SOX2','TP63')) %>% 
    filter(x>=11900) %>% 
    ggplot(aes(x,y))+
    geom_point(aes(color=gene),size=0.1,alpha=0.5)+
    theme1+scale_x_reverse()+
    # facet_wrap(~stage,scales = 'free')+
    theme(panel.background = element_blank(),plot.background = element_blank())+
    scale_color_manual(values = exprcol_gene)+guides(color=guide_legend(override.aes = list(size=3,alpha=1)))
  ggsave('stp2_epi/plot/epiexpr_spatial2/cycexpr_spatial_legend.pdf',p,width = 5,height = 5,limitsize = F)
  ggsave('stp2_epi/plot/tmp_cycscore/cycexpr.png',p,width = 40,height = 10,limitsize = F)
  
  for(stg in unique(epi_cycmol$stage)){
    p=epi_cycmol %>%
      filter(stage==stg) %>% 
      arrange(gene) %>% 
      ggplot(aes(x,y))+
      geom_point(aes(color=gene),size=0.1,alpha=1)+
      theme1+scale_x_reverse()+coord_fixed()+
      theme(panel.background = element_blank(),plot.background = element_blank())+
      scale_color_manual(values = exprcol_gene)+guides(color=guide_legend(override.aes = list(size=3,alpha=1)))
    ggsave(paste0('stp2_epi/plot/epiexpr_spatial2/cycexpr_spatial_',stg,'.png'),p,width = 14,height = 10,limitsize = F)
  }
  
  {# pergene
    rg = list(c(10600,11100,2600,3000),
              c(10700,11300,2600,3000),
              c(10700,11300,2600,3000))
    for(sp in 1:3)
    {
      dtrange = rg[[sp]]
      stg = unique(epi_cycmol$stage)[[sp]]
      dt = epi_cycmol %>% filter(stage==stg)
      # dt = dt %>% filter(x>=dtrange[1],x<=dtrange[2],y>=dtrange[3],y<=dtrange[4])
      p = lapply(as.character(unique(dt$gene)),function(g){
        dt1 = dt
        dt1$gene = dt1$gene %>% as.character()
        dt1$gene = ifelse(dt1$gene==g,g,'others') %>% factor(levels = c('others',g))
        dt1 %>% 
          arrange(gene) %>% 
          ggplot(aes(x,y))+
          geom_point(aes(color=gene),size=0.1,alpha=1)+
          theme1+scale_x_reverse()+NoLegend()+
          scale_color_manual(values = c(exprcol_gene,'others'='#cccccc'))+
          guides(color=guide_legend(override.aes = list(size=3,alpha=1)))
      }) %>% cowplot::plot_grid(plotlist = .,nrow = 4,ncol = 4)
      
      ggsave(paste0('stp2_epi/plot/tmp_cycscore/',stg,'_cycexpr.png'),p,width = 50,height = 40,limitsize = F)
    }
  }
}





##### gene expr, score box #####
epicol = c('#0339f8','#ff9408','#75bbfd','#de0c62')
names(epicol) = c('Basal','Proliferation','Differentiation','Invasive')

epidt = read_rds('stp2_epi/data/epiid_correctstage.rds')
cycscore = getscore(srtepi,list(cycling=gns)) %>% filter(cellID%in%epidt$cellID)
cycscore$celltype = cycscore$celltype %>% factor(c('Proliferation','Invasive','Basal','Differentiation'))

cycscore1 = aggregate(list(cycling=cycscore$cycling),by=as.list(cycscore[,c('sp_stg','stage1')]),mean)
cycscore1 %>% 
  ggplot(aes(stage1, cycling))+
  theme1+NoLegend()+
  # theme(panel.border = element_rect(fill=NA,color='black',linewidth=0.5),
        # axis.line = element_blank())+
  xlab('')+ylab('Proliferation score')+
  scale_fill_manual(values = col_stage)+
  stat_boxplot(geom = 'errorbar',size=0.5,width=0.5)+
  geom_boxplot(aes(fill=stage1),outlier.shape = 21,outlier.size = 3,outlier.stroke = 0.1)+
  stat_compare_means(comparisons = makecomparison(4))
ggsave('stp2_epi/plot/prolifeinv_score_stage_box.pdf',width = 4,height = 4)
pv1 = wilcox.test(cycscore$cycling[cycscore$stage1=='NOR'],cycscore$cycling[cycscore$stage1=='LGIN'])$p.value
pv2 = wilcox.test(cycscore$cycling[cycscore$stage1=='NOR'],cycscore$cycling[cycscore$stage1=='HGIN'])$p.value
pv3 = wilcox.test(cycscore$cycling[cycscore$stage1=='NOR'],cycscore$cycling[cycscore$stage1=='ESCC'])$p.value

aggregate(list(cycling=cycscore1$cycling),by=as.list(cycscore1[,c('stage1'),drop=F]),median)


cycscore2 = aggregate(list(cycling=cycscore$cycling),by=as.list(cycscore[,c('sp_stg','celltype')]),mean)
cycscore2$celltype = cycscore2$celltype %>% factor(c('Basal','Proliferation','Differentiation','Invasive'))
cycscore %>% 
  ggplot(aes(celltype, cycling))+
  theme1+NoLegend()+
  theme(panel.border = element_rect(fill=NA,color='black',size=0.5),
        axis.line = element_blank())+
  scale_fill_manual(values = epicol)+
  stat_compare_means(comparisons = makecomparison(4))+
  stat_boxplot(geom = 'errorbar',size=0.5,width=0.5)+
  geom_boxplot(aes(fill=celltype),outlier.shape = 21,outlier.size = 3,outlier.stroke = 0.1)
ggsave('stp2_epi/plot/prolifeinv_score_celltype_box.pdf',width = 4,height = 4)
pv1 = wilcox.test(cycscore$cycling[cycscore$celltype=='Proliferation'],cycscore$cycling[cycscore$celltype=='Invasive'])$p.value
pv2 = wilcox.test(cycscore$cycling[cycscore$celltype=='Proliferation'],cycscore$cycling[cycscore$celltype=='Basal'])$p.value
pv3 = wilcox.test(cycscore$cycling[cycscore$celltype=='Proliferation'],cycscore$cycling[cycscore$celltype=='Differentiation'])$p.value



aggregate(cycscore1$cycling,by=as.list(cycscore1[,c('stage1'),drop=F]),median)


# diffgene = c('ECM1','EMP1','SPRR2E','ANXA1','S100A8','CSTA')
# difscore = getscore(srtepi,list(diff=diffgene)) %>% filter(cellID%in%epidt$cellID)
# 
# difscore1 = aggregate(list(diff=difscore$diff),by=as.list(difscore[,c('sp_stg','celltype')]),mean)
# difscore1$celltype = difscore1$celltype %>% factor(c('Differentiation','Basal','Proliferation','Invasive'))
# difscore1 %>% 
#   ggplot(aes(celltype, diff))+
#   theme1+NoLegend()+
#   theme(panel.border = element_rect(fill=NA,color='black',size=0.5),
#         axis.line = element_blank())+
#   scale_fill_manual(values = epicol)+
#   geom_boxplot(aes(fill=celltype))+
#   stat_compare_means(comparisons = makecomparison(4))
# ggsave('stp2_epi/plot/diff_score_box.pdf',width = 3.5,height = 4)


load('stp2_epi/data/predDiffRFxenium.Rd')
predscore = data.frame(cellID=gsub('.','-',names(predDiffRFxenium.v),fixed = T),score=predDiffRFxenium.v)
predscore$dediffscore = 1-(predscore$score-min(predscore$score))/(max(predscore$score)-min(predscore$score))
predscore1 = predscore %>% inner_join(srtepi@meta.data)
predscore1$celltype = predscore1$celltype %>% factor(c('Differentiation','Proliferation','Basal','Invasive'))
predscore2 = aggregate(list(score=predscore1$dediffscore),by=as.list(predscore1[,c('sp_stg','stage1')]),mean)

predscore2 %>% 
  ggplot(aes(stage1, score))+
  theme1+NoLegend()+
  # theme(panel.border = element_rect(fill=NA,color='black',linewidth=0.5),
  #       axis.line = element_blank())+
  xlab('')+ylab('Differentiation score')+
  scale_fill_manual(values = col_stage)+
  stat_boxplot(geom = 'errorbar',linewidth=0.5,width=0.5)+
  geom_boxplot(aes(fill=stage1),outlier.shape = 21,outlier.size = 3,outlier.stroke = 0.1)+
  stat_compare_means(comparisons = makecomparison(4))
ggsave('stp2_epi/plot/dediff_score_stage_box.pdf',width = 4,height = 4)
pv1 = wilcox.test(predscore1$score[cycscore$stage1=='NOR'],cycscore$score[cycscore$celltype=='LGIN'])$p.value
pv2 = wilcox.test(predscore1$score[cycscore$stage1=='NOR'],cycscore$score[cycscore$celltype=='HGIN'])$p.value
pv3 = wilcox.test(predscore1$score[cycscore$stage1=='NOR'],cycscore$score[cycscore$celltype=='ESCC'])$p.value

aggregate(list(score=predscore2$score),by=as.list(predscore2[,c('stage1'),drop=F]),median)


predscore1 %>% 
  ggplot(aes(celltype, score))+
  theme1+NoLegend()+
  theme(panel.border = element_rect(fill=NA,color='black',size=0.5),
        axis.line = element_blank())+
  scale_fill_manual(values = epicol)+
  stat_boxplot(geom = 'errorbar',size=0.5,width=0.5)+
  geom_boxplot(aes(fill=celltype),outlier.shape = 21,outlier.size = 3,outlier.stroke = 0.1)+
  stat_compare_means(comparisons = makecomparison(4))
ggsave('stp2_epi/plot/diff_score_celltype_box.pdf',width = 4,height = 4)

pv1 = wilcox.test(predscore1$score[cycscore$celltype=='Differentiation'],cycscore$score[cycscore$celltype=='Proliferation'])$p.value
pv2 = wilcox.test(predscore1$score[cycscore$celltype=='Differentiation'],cycscore$score[cycscore$celltype=='Basal'])$p.value
pv3 = wilcox.test(predscore1$score[cycscore$celltype=='Differentiation'],cycscore$score[cycscore$celltype=='Invasive'])$p.value


predscore3 = aggregate(list(score=predscore1$dediffscore),by=as.list(predscore1[,c('sp_stg','stage1','celltype','pid')]),mean)
predscore3 %>% 
  filter(celltype%in%c('Proliferation')) %>% 
  ggplot(aes(stage1, score))+
  theme1+NoLegend()+
  # theme(panel.border = element_rect(fill=NA,color='black',size=0.5),
  #       axis.line = element_blank())+
  xlab('')+ylab('Stemness score')+
  scale_fill_manual(values = col_stage)+
  # scale_y_continuous(limits = c(0,1))+
  facet_wrap(~celltype,scales = 'free')+
  stat_boxplot(geom = 'errorbar',size=0.5,width=0.5)+
  geom_boxplot(aes(fill=stage1),outlier.shape = 21,outlier.size = 3,outlier.stroke = 0.1)+
  stat_compare_means(comparisons = makecomparison(4))
ggsave('stp2_epi/plot/dediff_score_stagecelltype_box.pdf',width = 4,height = 4)

pv1 = wilcox.test(predscore1$score[cycscore$celltype=='Differentiation'],cycscore$score[cycscore$celltype=='Proliferation'])$p.value
pv2 = wilcox.test(predscore1$score[cycscore$celltype=='Differentiation'],cycscore$score[cycscore$celltype=='Basal'])$p.value
pv3 = wilcox.test(predscore1$score[cycscore$celltype=='Differentiation'],cycscore$score[cycscore$celltype=='Invasive'])$p.value

aggregate(predscore3$score,by=as.list(predscore3[,c('stage1','celltype')]),median)



{# cell level stage violin
  p = lapply(c('Basal','Proliferation','Invasive'),function(x){
    dt = predscore1 %>% filter(celltype%in%x)
    compnum = unique(dt$stage1) %>% length
    
    dt %>%
      ggplot(aes(stage1, dediffscore))+
      theme1+NoLegend()+
      xlab('')+ylab('Stemness score')+ggtitle(x)+
      stat_compare_means(comparisons = makecomparison(compnum))+
      scale_fill_manual(values = col_stage)+
      scale_y_continuous(limits = c(0,1.3),breaks = seq(0,1,0.5),labels = c('0','0.5','1.0'))+
      geom_violin(aes(fill=stage1),color=NA)+
      stat_summary(geom='errorbar',fun.min = median,fun.max = median)
    
    # lapply(unique(dt$stage1)[-1],function(y){
    #   dt1 = dt %>% filter(stage1==unique(dt$stage1)[1])
    #   dt2 = dt %>% filter(stage1==y)
    #   pv=wilcox.test(dt1$dediffscore,dt2$dediffscore)$p.value
    # })
  })
  ggsave('stp2_epi/plot/dediff_score_stagecelltype_celllevel_violin.pdf',cowplot::plot_grid(plotlist = p,nrow = 1,rel_widths = c(2,3,2)),width = 8,height = 4)
}

{# cell level stage violin, patient t1 t3
  tmp = predscore1 %>% filter(celltype%in%c('Proliferation','Basal')) %>% filter(pid%in%c('T1','T3'))
  tmp$grp = tmp$stage1 %>% factor(levels = c('Basal','NOR','LGIN','HGIN','ESCC'))
  tmp$grp[tmp$celltype=='Basal'] = 'Basal'
  
  tmp %>% 
    ggplot(aes(grp, dediffscore))+
    theme1+NoLegend()+
    xlab('')+ylab('Stemness score')+
    facet_wrap(~pid)+
    stat_compare_means(comparisons = makecomparison(4,startnum = 2))+
    scale_fill_manual(values = c(epicol,col_stage))+
    scale_y_continuous(limits = c(0,1.3),breaks = seq(0,1,0.5),labels = c('0','0.5','1.0'))+
    geom_violin(aes(fill=grp),color=NA)+
    stat_summary(geom='errorbar',fun.min = median,fun.max = median)
  ggsave('stp2_epi/plot/dediff_score_stageprolife_celllevel_violin.pdf',width = 6,height = 4)
  
  p = lapply(c('T2','T4','T5','T6'),function(x){
    dt = predscore1 %>% filter(pid%in%x) %>% filter(celltype=='Proliferation')
    compnum = unique(dt$stage1) %>% length
    
    dt %>%
      ggplot(aes(stage1, dediffscore))+
      theme1+NoLegend()+
      xlab('')+ylab('Stemness score')+
      facet_wrap(~pid)+
      stat_compare_means(comparisons = makecomparison(compnum))+
      scale_fill_manual(values = col_stage)+
      scale_y_continuous(limits = c(0,1.3),breaks = seq(0,1,0.5),labels = c('0','0.5','1.0'))+
      geom_violin(aes(fill=stage1),color=NA)+
      stat_summary(geom='errorbar',fun.min = median,fun.max = median)
    
    # lapply(unique(dt$stage1)[-1],function(y){
    #   dt1 = dt %>% filter(stage1==unique(dt$stage1)[1])
    #   dt2 = dt %>% filter(stage1==y)
    #   pv=wilcox.test(dt1$dediffscore,dt2$dediffscore)$p.value
    # })
  })
  ggsave('stp2_epi/plot/dediff_score_stageprolife_celllevel_violin_supp.pdf',cowplot::plot_grid(plotlist = p,nrow = 1,rel_widths = c(2,3,2.5,2.5)),width = 10,height = 4)
}






##### prolife diff 2d plot #####
library(ggstar)
cyc_diffscore = inner_join(cycscore,predscore)
cyc_diffscore$celltype = cyc_diffscore$celltype %>% factor(c('Basal','Proliferation','Differentiation','Invasive'))
cyc_diffscore = cyc_diffscore %>% dplyr::rename(c('prolifescore'='cycling','diffscore'='dediffscore'))
cyc_diffscore1 = aggregate(list(prolifescore=cyc_diffscore$prolifescore,diffscore=cyc_diffscore$diffscore),
                           by=as.list(cyc_diffscore[,c('celltype','sp_stg')]),mean)

cyc_diff_median = aggregate(list(prolifescore=cyc_diffscore$prolifescore,diffscore=cyc_diffscore$diffscore),
                            by=as.list(cyc_diffscore[,c('celltype'),drop=F]),mean)
rownames(cyc_diff_median) = cyc_diff_median$celltype

cyc_diffscore %>% 
  ggplot(aes(prolifescore,diffscore))+
  scale_color_manual(values = epicol)+scale_fill_manual(values = epicol)+
  geom_point(size=0.1,aes(color=celltype))+
  geom_star(data=cyc_diff_median,aes(fill=celltype),size=5)+
  geom_hline(yintercept = median(cyc_diffscore$diffscore))+
  geom_vline(xintercept = median(cyc_diffscore$prolifescore))+
  theme1+
  xlim(range(cyc_diffscore$prolifescore))+ylim(range(cyc_diffscore$diffscore))
ggsave('stp2_epi/plot/prolife_diff_scatter.png',width = 6,height=4)
# ggsave('stp2_epi/plot/prolife_diff_scatter_foreg.pdf',width = 6,height=4)

cyc_diffscore %>% 
  ggplot(aes(prolifescore,diffscore))+
  scale_color_manual(values = epicol)+scale_fill_manual(values = epicol)+
  facet_wrap(~celltype)+
  # geom_point(size=0.1,aes(color=celltype))+
  # geom_density_2d_filled()+
  geom_density_2d(binwidth=0.2,aes(color=celltype))+
  geom_star(data=cyc_diff_median,aes(fill=celltype),size=5)+
  geom_hline(yintercept = median(cyc_diffscore$diffscore))+
  geom_vline(xintercept = median(cyc_diffscore$prolifescore))+
  theme1+
  xlim(range(cyc_diffscore$prolifescore))+ylim(range(cyc_diffscore$diffscore))
ggsave('stp2_epi/plot/prolife_diff_scatter_facet.png',width = 10,height=8)

cyc_diffscore %>%
  ggplot(aes(prolifescore))+
  theme1+NoLegend()+
  xlim(c(0,1.5))+
  scale_color_manual(values = epicol)+scale_fill_manual(values = epicol)+
  # geom_point(size=0.1,aes(color=celltype))+
  geom_density(aes(fill=celltype),alpha=0.5)
ggsave('stp2_epi/plot/prolifeinv_score_celltype_density.pdf',width = 9,height = 3)
cyc_diffscore %>%
  ggplot(aes(diffscore))+
  theme1+NoLegend()+
  # xlim(c(0,1.5))+
  scale_color_manual(values = epicol)+scale_fill_manual(values = epicol)+
  # geom_point(size=0.1,aes(color=celltype))+
  geom_density(aes(fill=celltype),alpha=0.5)
ggsave('stp2_epi/plot/dediff_score_celltype_density.pdf',width = 9,height = 3)





heatdt = cyc_diff_median[,2:3] %>% t
heatdt = heatdt[,c('Differentiation','Basal','Proliferation','Invasive')]

col1 = c(brewer.pal(9,'Reds'))
col2 = c(brewer.pal(9,'Blues'))

# ComplexHeatmap::pheatmap(heatdt[c(1,1),],scale = 'row',color = col1,cluster_cols = F,cluster_rows = F)
# ComplexHeatmap::pheatmap(heatdt[c(2,2),],scale = 'row',color = col2,cluster_cols = F,cluster_rows = F)

heatdt1 = heatdt %>% t %>% scale
heatdt1 = heatdt1 %>% data.frame %>% 
  mutate(celltype=rownames(.)) %>% gather('prog','score',-celltype)
heatdt1$celltype = heatdt1$celltype %>% factor(levels = c('Differentiation','Basal','Proliferation','Invasive'))
heatdt1 %>% 
  ggplot(aes(celltype,prog))+
  theme1+RotatedAxis()+
  scale_size_continuous(range = c(2,8))+
  scale_fill_gradientn(colors = col2)+
  geom_point(data=heatdt1[heatdt1$prog=='diffscore',],shape=21,aes(fill=score,size=score))+
  new_scale('fill')+
  scale_fill_gradientn(colors = col1)+
  geom_point(data=heatdt1[heatdt1$prog=='prolifescore',],shape=21,aes(fill=score,size=score))
ggsave('stp2_epi/plot/prolife_diff_dotplot.pdf',width = 5,height = 5)


heatdt2 = aggregate(list(prolifescore=cyc_diffscore$prolifescore,diffscore=cyc_diffscore$diffscore),
                    by=as.list(cyc_diffscore[,c('celltype','stage1'),drop=F]),mean)
heatdt2$celltype = heatdt2$celltype %>% factor(levels = c('Basal','Proliferation','Differentiation','Invasive'))
heatdt2 = heatdt2 %>% gather('prog','score',prolifescore,diffscore)

heatdt2 %>% 
  ggplot(aes(celltype,prog))+
  theme1+RotatedAxis()+
  facet_wrap(~stage1,nrow=1,scales='free_x')+
  scale_size_continuous(range = c(2,8))+
  scale_fill_gradientn(colors = col2)+
  geom_point(data=heatdt2[heatdt2$prog=='diffscore',],shape=21,aes(fill=score,size=score))+
  new_scale('fill')+
  scale_fill_gradientn(colors = col1)+
  geom_point(data=heatdt2[heatdt2$prog=='prolifescore',],shape=21,aes(fill=score,size=score))
ggsave('stp2_epi/plot/prolife_diff_dotplot_stagect.pdf',width = 10,height = 5)

