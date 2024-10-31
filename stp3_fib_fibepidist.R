source('~/xenium/header.R')
source('~/xenium/xenium_function.R')
library(scatterpie)




# input #####
xe1_3 = read_rds('data/xenium_seurat_processed/xenium_t1_s3_seurat.rds')
xe1_3@meta.data = read_rds('data/xenium_seurat_processed/xenium_t1_s3_meta.rds')

xe3_1 = read_rds('data/xenium_seurat_processed/xenium_t3_s1_seurat.rds')
xe3_1@meta.data = read_rds('data/xenium_seurat_processed/xenium_t3_s1_meta.rds')

xe4_2 = read_rds('data/xenium_seurat_processed/xenium_t4_s2_seurat.rds')
xe4_2@meta.data = read_rds('data/xenium_seurat_processed/xenium_t4_s2_meta.rds')

xe6_1 = read_rds('data/xenium_seurat_processed/xenium_t6_s1_seurat.rds')
xe6_1@meta.data = read_rds('data/xenium_seurat_processed/xenium_t6_s1_meta.rds')

srtfib = read_rds('stp3_tme/data/fib/fibro_all_seurat.rds')
srtfib@meta.data = read_rds('stp3_tme/data/fib/fibro_all_meta.rds')
srtfib$stage1 = srtfib$stage1 %>% factor(levels = c('NOR','LGIN','HGIN','ESCC'))
srtepi = read_rds('data/xenium_seurat_processed/ct_epi_srt_2.rds')
fibtoepi = read_rds('stp3_tme/data/fib-epi_dist/fib-epi_dist.rds')

epidt = read_rds('stp2_epi/data/epiid_correctstage.rds')
crds = read_rds('stp1_summary/data/allsample_coords.rds')
epipoly = read_rds('stp3_tme/data/epiregion/allsample_poly.rds')
allborder = read_rds('stp3_tme/data/epiborder/allborder.rds')

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
            'Metastasis'="#ff0000")

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

prognew = read_rds('stp2_epi/data/epiProgram/epigene_annotation_v4.rds')


allborder = read_rds('stp3_tme/data/epiborder/allborder.rds')
t1dist = read_rds('stp3_tme/data/epi_borderdist/T1epi_borderdist_v2.rds')
t1border = read_rds('stp3_tme/data/epiborder/T1border_v2.rds')
t3dist = read_rds('stp3_tme/data/epi_borderdist/T3ESCC_borderdist_v2.rds')
t3border = read_rds('stp3_tme/data/epiborder/T3ESCCborder_v2.rds')
t4dist = read_rds('stp3_tme/data/epi_borderdist/T4partE_borderdist_v2.rds')
t4dtlist = read_rds('stp3_tme/data/epiborder/T4partEborder_v2.rds')




# fib dist expr #####
fibtoinv = read_rds('stp3_tme/data/fib-epi_dist/fib-inv_dist.rds')
grp = fibtoepi$fibdist1
groupdist = 1
fibtoepi$group1 = (grp%/%groupdist+1)*groupdist
fibtoepi$group = paste0(fibtoepi$group1,'um') %>% 
  factor(levels = paste0(sort(unique(fibtoepi$group1)),'um'))


srtfib = read_rds('stp3_tme/data/fib/fibro_all_seurat.rds')
srtfib@meta.data = read_rds('stp3_tme/data/fib/fibro_all_meta.rds')
srtfib$stage1 = srtfib$stage1 %>% factor(levels = c('NOR','LGIN','HGIN','ESCC'))
srtfib$celltype = srtfib$cellType_merge %>% mapvalues('NF-PI16','NF')

cafgn = c('FAP','MMP1','MMP11','POSTN')
expr = getexpr(srtfib,cafgn) %>% left_join(srtfib@meta.data) %>% 
  filter(stage1%in%c('HGIN','ESCC')) %>% inner_join(fibtoinv)

expr_aggr = aggregate(list(expression=expr$expression),by=as.list(expr[,c('group1','stage1','celltype','gene')]),mean)
expr_aggr = expr_aggr %>% 
  filter(group1<150) %>% 
  # filter(gene=='POSTN') %>% 
  filter(celltype=='NF')

# trend test
library(trend)

trendtest = lapply(unique(expr_aggr$stage1),function(stg){
  print(stg)
  dt1 = expr_aggr %>% filter(stage1==stg)
  rst1 = lapply(unique(dt1$gene),function(g){
    print(g)
    dt2 = dt1 %>% filter(gene==g)
    dt2 = dt2 %>% arrange(group1)
    testrst = mk.test(dt2$expr)
    rst2 = data.frame(sval=testrst$estimates[1],zval=testrst$statistic,pval=testrst$p.value,gene=g,stage1=stg)
  }) %>% Reduce(rbind,.)
}) %>% Reduce(rbind,.)
trendtest$pval_adj = p.adjust(trendtest$pval,method = 'fdr')
labdt = expr_aggr
labdt = aggregate(list(yloc=labdt$expr),by=as.list(labdt[,c('stage1','gene')]),min) %>% left_join(trendtest)

expr_aggr %>% 
  # filter(group1>20) %>% 
  ggplot(aes(group1,expression))+
  theme1+xlab('Distance to Invasive cell')+ylab('Gene expression')+
  scale_color_manual(values = c('POSTN'='#e4282a','FAP'='#ffdd2e','MMP11'='#f3a058','MMP1'='#206ed0'))+
  geom_point(aes(group=gene,color=gene),shape=16)+
  geom_smooth(method = 'gam')+
  geom_text(data=labdt,x=5,aes(y=yloc,label=paste0('Trend test\np = ',signif(pval_adj,2),'\nz = ',signif(zval,2))),hjust=0,vjust=0)+
  # stat_cor(method = 'spearman')+
  facet_wrap(~stage1+gene,nrow = 2,scales = 'free')
ggsave('stp3_tme/plot/fibepi_dist/Nf_cafexpr_epidist_trend.pdf',width = 12,height = 6)





# progscore #####
{## cor with dist
  epicaf = read_rds('stp3_tme/data/fib-epi_dist/epi-CAF_dist.rds')
  
  prognew = read_rds('stp2_epi/data/epiProgram/epigene_annotation_v5.rds')
  gnlist = prognew
  gnlist = split(gnlist$gene,gnlist$annotation)
  
  exprdt = read_rds('stp2_epi/data/epiProgram/epiprogscore_v5.rds')
  rownames(exprdt) = exprdt$cellID
  exprdt = exprdt[,c('cellID','pid','stage1',names(gnlist)),drop=F] %>% inner_join(epicaf)
  exprdt = exprdt %>% gather('prog','score',all_of(names(gnlist)))
  
  grp = exprdt$epidist1
  groupdist = 1
  exprdt$group1 = (grp%/%groupdist+1)*groupdist
  exprdt$group1 = ifelse(exprdt$group1>150,min(exprdt$group1[exprdt$group1>150]),exprdt$group1)
  exprdt$group = paste0(exprdt$group1,'um') %>% 
    factor(levels = paste0(sort(unique(exprdt$group1)),'um'))
  
  exprdt_aggr = aggregate(list(score=exprdt$score),by=as.list(exprdt[,c('stage1','prog','group1','group')]),mean)
}


## plot #####
{
  exprdt_aggr %>% 
    ggplot(aes(group1,score))+
    theme1+
    geom_point()+
    facet_wrap(~stage1+prog,nrow = 4,scales = 'free')
  
  ## trend test, mann kendall
  library(trend)
  
  trendtest = lapply(unique(exprdt_aggr$stage1),function(stg){
    print(stg)
    dt1 = exprdt_aggr %>% filter(stage1==stg)
    rst1 = lapply(unique(dt1$prog),function(prg){
      print(prg)
      dt2 = dt1 %>% filter(prog==prg)
      testrst = mk.test(dt2$score)
      rst2 = data.frame(rval=testrst$estimates[3],zval=testrst$statistic,pval=testrst$p.value,prog=prg,stage1=stg)
    }) %>% Reduce(rbind,.)
  }) %>% Reduce(rbind,.) 
  trendtest$pval_adj = p.adjust(trendtest$pval,method = 'fdr')
  
  labdt = exprdt_aggr
  labdt$group1 = labdt$group1 %>% as.character %>% as.numeric()
  labdt = aggregate(list(yloc=labdt$score),by=as.list(labdt[,c('prog'),drop=F]),max) %>% left_join(trendtest)
  
  exprdt_aggr %>% 
    ggplot(aes(group1,score))+
    theme1+theme(panel.spacing.y = unit(0.1,'cm'),panel.border = element_rect(fill=NA,color='black',linewidth = 0.5),strip.background.y = element_rect(fill='grey'))+
    xlab('Distance')+ylab('Program score')+
    scale_y_continuous(expand = expansion(c(0.05,0.1)))+
    scale_color_manual(values = progcol)+scale_fill_manual(values = progcol)+
    geom_point(aes(color=prog),shape=16,size=0.5)+
    geom_smooth(linewidth=0.8)+
    # stat_cor(method = 'kendal')+
    geom_text(data=labdt,x=0.1,aes(y=yloc*1.1,label=paste0('r = ',signif(rval,2),'\nP = ',signif(pval_adj,2))),hjust=0,vjust=1)+
    facet_grid(prog~stage1,scales = 'free')
  ggsave('stp3_tme/plot/fibepi_dist/epiprog_corwithcafdist.pdf',width = 10,height = 9)
}




# JAG1 NOTCH1 expr #####
{## cor with dist
  epicaf = read_rds('stp3_tme/data/fib-epi_dist/epi-CAF_dist.rds')
  
  tmp = srtepi
  tmp@meta.data = tmp@meta.data %>% left_join(epicaf) %>% `rownames<-`(.$cellID)
  
  grp = tmp$epidist1
  groupdist = 1
  tmp$group1 = (grp%/%groupdist+1)*groupdist
  tmp$group1 = ifelse(tmp$group1>150,min(tmp$group1[tmp$group1>150]),tmp$group1)
  tmp$group = paste0(tmp$group1,'um') %>% 
    factor(levels = paste0(sort(unique(tmp$group1)),'um'))
  
  isout = function(x){
    # return(x < quantile(x,0.25) - 1.5 * IQR(x) | x > quantile(x,0.75) + 1.5 * IQR(x))
    return(x < mean(x)-1*sd(x) | x > mean(x)+1*sd(x))
  }
  
  # gnexpr = aggrExpr(tmp,c('NOTCH1','JAG1'),aggrby = tmp@meta.data[,c('sp_stg','stage1','group1','group','pid')]) %>% gather('gene','expr',NOTCH1,JAG1)
  # gnexpr = aggrExpr(tmp,c('NOTCH1','JAG1'),aggrby = tmp@meta.data[,c('stage1','group1','group')]) %>% gather('gene','expr',NOTCH1,JAG1)
  gnexpr = getexpr(tmp,c('NOTCH1','JAG1')) %>% `colnames<-`(c('cellID','gene','expr')) %>% left_join(tmp@meta.data)
  gnexpr = gnexpr %>% filter(pid!='T2')
  gnexpr = aggregate(list(expr=gnexpr$expr),by=as.list(gnexpr[,c('stage1','group1','group','gene')]),mean)
  gnexpr = gnexpr %>% 
    # group_by(interaction(gene,stage1,group)) %>% 
    # mutate(outlier=ifelse(isout(expr),pid,NA)) %>% 
    # filter(pid%in%c('T2')==F) %>% 
    # filter(is.na(outlier)) %>%
    filter(stage1%in%c('HGIN','ESCC')) %>% 
    filter(group1<=200)
  # gnexpr = getexpr(tmp,c('NOTCH1','JAG1')) %>% `colnames<-`(c('cellID','gene','expr')) %>% left_join(tmp@meta.data)
  
  # trend test
  library(trend)
  
  trendtest = lapply(unique(gnexpr$stage1),function(stg){
    print(stg)
    dt1 = gnexpr %>% filter(stage1==stg)
    rst1 = lapply(unique(dt1$gene),function(g){
      print(g)
      dt2 = dt1 %>% filter(gene==g)
      dt2 = dt2 %>% arrange(group1)
      testrst = mk.test(dt2$expr)
      rst2 = data.frame(sval=testrst$estimates[1],zval=testrst$statistic,pval=testrst$p.value,gene=g,stage1=stg)
    }) %>% Reduce(rbind,.)
  }) %>% Reduce(rbind,.)
  trendtest$pval_adj = p.adjust(trendtest$pval,method = 'fdr')
  labdt = gnexpr
  labdt = aggregate(list(yloc=labdt$expr),by=as.list(labdt[,c('stage1','gene')]),min) %>% left_join(trendtest)
  
  gnexpr %>% 
    ggplot(aes(group1,expr))+
    theme1+xlab('Distance to CAF')+ylab('Gene expression')+
    geom_point(aes(color=gene))+
    # scale_color_gradientn(colours = distcol,name='Distance to CAF')+
    scale_color_manual(values = c('JAG1'='#7fbc41','NOTCH1'='#cab2d6'))+
    geom_smooth(method='gam')+
    geom_text(data=labdt,x=5,aes(y=yloc,label=paste0('Trend test\np = ',signif(pval_adj,2),'\nz = ',signif(zval,2))),hjust=0,vjust=0)+
    # stat_cor(method='spearman')+
    facet_wrap(~stage1+gene,scales = 'free')
  ggsave('stp3_tme/plot/fibepi_dist/jagnotch_corwithcafdist.pdf',width = 7,height = 6)
  
  # pval
  # lapply(c('JAG1','NOTCH1'),function(gn){
  #   dt1 = gnexpr %>% filter(gene==gn)
  #   lapply(c('HGIN','ESCC'),function(stg){
  #     dt2 = dt1 %>% filter(stage1==stg)
  #     dt2$group1 = dt2$group1 %>% as.character %>% as.numeric()
  #     tst = cor.test(dt2$expr,dt2$group1,method='spearman')
  #     return(data.frame(r=tst$estimate,pv=tst$p.value,stage=stg,gene=gn))
  #   }) %>% Reduce(rbind,.)
  # }) %>% Reduce(rbind,.)
}



