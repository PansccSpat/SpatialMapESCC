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
allmk2 = c('SFN','KRT4','KRT5','EPCAM', # epi
           'DCN','COL1A1','COL1A2','COL3A1', # fibro
           'PTPRC', # immune
           'CD2','CD3D','CD3E', # t cell
           'CD19','CD79A','MS4A1', # B
           'JCHAIN','IGHG2', # plasma
           'VWF','PECAM1','ENG','CDH5', # endo
           'CD68','LYZ','CD14','FCGR3A', # myeloid
           'CPA3', # mast
           'PLN','MYH11','ACTG2','CNN1','MYL9','TAGLN', # muscle
           'AGR2','KRT23','WFDC2' # gland
)
allmk3 = c('SFN','KRT4','KRT5','EPCAM', # epi
           'VWF','PECAM1','ENG','CDH5', # endo
           'COL1A1','COL1A2','COL3A1','PI16','FAP', # fibro
           'PTPRC', # immune
           'CD2','CD3D','CD3E','CD4','CD8A','CD8B','GZMH','PDCD1','HAVCR2','CD200','FOXP3', # t cell
           'CD19','CD79A','MS4A1', # B
           'JCHAIN','IGHG2', # plasma
           'CD68','LYZ','CD14','FCGR3A', # myeloid
           'CPA3', # mast
           'AGR2','KRT23','WFDC2', # gland
           'PLN','MYH11','ACTG2','CNN1','MYL9','TAGLN' # muscle
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

myemarker = c('SFN','DCN','VWF','CD19','CD68','CD2', # pan marker
              'GZMB','JCHAIN', # pDC
              'CD1C','CD1E','FCER1A', # APC
              'CCR7', # tDC
              'FCN1','S100A9','S100A8','FCGR3A', # Mono
              'IL1B', # Macro
              'LYVE1','C1QC','C1QB','C1QA','SPP1','FN1' # Macro
)

cd4marker = c('CD4','CD8A','CD8B',
              'FOXP3','IL2RA','CTLA4', # reg
              'TCF7','CCR7','SELL', # tn
              'IL7R','CCL5', # mem
              'CXCL13','CD200' # tfh
)

cd8marker = c('CD4','CD8A','CD8B',
              'TCF7','CCR7','SELL', # tn
              'IL7R','CXCR4',#mem
              'GNLY', 'NKG7', 'GZMB', 'GZMA', 'GZMH', 'GZMK', #Teff
              'PDCD1','CXCL13','CTLA4','HAVCR2' # tex 
)
cytogene = c('GNLY', 'NKG7', 'GZMB', 'GZMA', 'GZMH', 'GZMK', 'KLRB1','TYROBP')  
exhgene = c('PDCD1','HAVCR2','CTLA4','TIGIT','LAG3','CXCL13') 
reggene = c('FOXP3','TIGIT','CD274','IL2RA')
tmarker = c('CD4','CD8A','CD8B',
            'GNLY', 'NKG7', 'GZMB', 'GZMA', 'GZMH', 'GZMK', 'KLRB1','TYROBP',
            'PDCD1','HAVCR2','CTLA4','TIGIT','LAG3','CXCL13',
            'FOXP3','CD274','IL2RA')

fibmarker = c('PLA2G2A','GPX3','IGFBP6','PI16','ADH1B','TNXB','DCN','CXCL12', # NF
              'IGF1','RGS5','HIF1A','CCND1','POSTN','STAT1','FAP','SERPINE1','DPT', # iCAF
              'MMP11','MMP1','IL7R','STMN1','PDPN','RUNX1' # myCAF
)



##### colors #####
primcol = brewer.pal(10, 'Paired')
primcol[7] = '#333333'
names(primcol) = c('B','Epi','Endo','Fib','Gland','Mye','Mus','Mast','Pla','T')
primcol1 = c('#1f78b4','#ffc089','#B5EFB5','#793c1b',
             '#6a3d9a','#333333','#ffff33',
             '#fbb4ae','#4a8cff','#ff8820','#a94e4e','#aa0000','#d7000588')
names(primcol1) = c('B','Endo','Fib','Gland','Mye','Mus','T','NOR','LGIN','HGIN','ESCC')

primcol2 = c('#1f78b4','#ffc089','#B5EFB5','#793c1b',
             '#6a3d9a','#333333','#ffff33',
             '#df928e')
names(primcol2) = c('B','Endo','Fib','Gland','Mye','Mus','T','Epi')

# primcol2 = c('#006837','#ffc089','white','#793c1b',
#              '#6a3d9a','#333333','#08519c',
#              '#df928e')
# names(primcol2) = c('B','Endo','Fib','Gland','Mye','Mus','T','Epi')

primcol3_ = c('#793c1b','#333333','#ffc089','#df928e',
              '#1f78b4','#6a3d9a','#df65b0','#1cbe4f','#e31a1c','#ffdcbd',
              '#ffffd4','#ffff68','#bbbb25','#ffe4cc','#ff9344','#a15008')
names(primcol3_) = c('Gland','Mus','Endo','Epi',
                     'B','Mye','NF-PI16','NF','CAF','T-others',
                     'CD8T-n/m','CD8T-eff','CD8T-ex','CD4T-n/m','CD4T-fh','CD4T-reg')
primcol3 = c('#793c1b','#333333','#ffc089','#df928e',
             '#1f78b4','#6a3d9a','#1cbe4f','#df65b0','#e31a1c','#ffdcbd',
             '#BDCDFF','#ffff68','#325A9B','#AAF400','#ff9344')
names(primcol3) = c('Gland','Mus','Endo','Epi',
                    'B','Mye','NF','NF-PI16','CAF','T-others',
                    'T-naive','CD8T-eff','CD8T-ex','CD4T-fh','CD4T-reg')

allcellcol = c('Basal'='#0339f8','Proliferation'='#ff9408','Differentiation'='#75bbfd',
               'Terminal'='#4682B4','HT_Prolif'='#ff9408','HT_Diff/Term'='#4682B4',
               'Invasive'='#de0c62','T cells'='#96f97b','B cells'='#028f1e',
               'Myeloid cells'='#7e1e9c','Fibroblast'='#696969','CAF'='#FFFFFF',
               'Endothelial cells'='#7b002c','Gland cells'='#c69f59','Myocytes'='#1b2431')



redcol = c('#412626','#d62728','#fbb1b2') %>% colorRampPalette()
yelcol = c('#77732d','#ede659','#faf9d5') %>% colorRampPalette()
bluecol = c('#266065','#17becf','#7aecf8') %>% colorRampPalette()
greencol = c('#202f20','#2ca02c','#b9fbb9') %>% colorRampPalette()
orgcol = c('#7c4d24','#ff7f0e','#ffdcbd') %>% colorRampPalette()
pinkcol = c('#a34085','#e377c2','#f29ed8') %>% colorRampPalette()
purpcol = c('#7f52a9','#c5b0d5') %>% colorRampPalette()
c('#b03533','#ff9896','#fed2d2')
polychrome = c("#5A5156", 
               "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", 
               "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", 
               "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", 
               "#1C8356", "#85660D", "#B10DA1", "#FBE426", "#1CBE4F", 
               "#FA0087", "#FC1CBF", "#F7E1A0", "#C075A6", "#782AB6", 
               "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5", "#7ED7D1", 
               "#1C7F93", "#D85FF7", "#683B79", "#66B0FF", "#3B00FB")

ctcol = c(bluecol(7),
          yelcol(3)[-1],
          redcol(9),
          pinkcol(3),
          purpcol(3),'#8a4071',
          orgcol(3)[-1],
          greencol(8),'#999999')
names(ctcol) = c(paste0('NF-',1:3),'iCAF-1',paste0('myCAF-',1:3),
                 'Pla','B',
                 paste0('Mac-',c('SPP1','IL1B','C1QC','LYVE1')),'cDC','pDC','Mo','Mo-Endo','Mast',
                 paste0('CD4T-',c('n/m','fh','reg')),
                 paste0('CD8T-',c('n/m','eff','ex')),'T-others',
                 'Endo','Gland',
                 'HY','DO-spp1','DO','AP','CY','QP','MD','TD',
                 'Mus')
ctcol1 = c(bluecol(7),
          yelcol(3)[-3],
          redcol(9),
          pinkcol(3),
          purpcol(3),'#8a4071',
          orgcol(3)[-1],
          '#2ca02c','#999999')
names(ctcol1) = c(paste0('myCAF-',3:1),'iCAF-1',paste0('NF-',3:1),
                 'Pla','B',
                 paste0('Mac-',c('SPP1','IL1B','C1QC','LYVE1')),'cDC','pDC','Mo','Mo-Endo','Mast',
                 paste0('CD4T-',c('n/m','fh','reg')),
                 paste0('CD8T-',c('n/m','eff','ex')),'T-others',
                 'Endo','Gland',
                 'Epi',
                 'Mus')

## stage col
colepi = c('#fdddda','#d2817d','#984646','#4c2323') %>% `names<-`(paste0('Epi-',c('NOR','LGIN','HGIN','ESCC')))
colfib = c('#e8fae8','#b5efb5','#7ea77e','#364736') %>% `names<-`(paste0('Fib-',c('NOR','LGIN','HGIN','ESCC')))
colendo = c('#ffecdb','#ffc089','#987352','#4c3929') %>% `names<-`(paste0('Endo-',c('NOR','LGIN','HGIN','ESCC')))
colmus = c('#d6d6d6','#adadad','#5b5b5b','#333333') %>% `names<-`(paste0('Mus-',c('NOR','LGIN','HGIN','ESCC')))
colgland = c('#e4d8d1','#ae8a76','#793c1b','#3c1e0d') %>% `names<-`(paste0('Gland-',c('NOR','LGIN','HGIN','ESCC')))
colt = c('#ffffd6','#ffff33','#b2b223','#32320a') %>% `names<-`(paste0('T-',c('NOR','LGIN','HGIN','ESCC')))
colb = c('#d2e3f0','#8fbbd9','#185f8f','#0f3c5a') %>% `names<-`(paste0('B-',c('NOR','LGIN','HGIN','ESCC')))
colmye = c('#e1d8ea','#9677b8','#54307b','#351e4d') %>% `names<-`(paste0('Mye-',c('NOR','LGIN','HGIN','ESCC')))
allcol = c(colepi,colfib,colendo,colmus,colgland,colt,colb,colmye)




##### spatial plot #####
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
    ggsave(fname,p,width = 40, height = 15, limitsize = F)
    
    fname = paste0('stp1_summary/stage-prim_plot/',x,'_primarytype.png')
    dt$tmp = paste0(dt$ct3,'-',dt$stage1)
    p = xedimplot(dt, groupby = 'tmp',flip = flip,reversex = reversex,reversey = reversey,color = allcol,bgcol = 'black',pt.size = 0.1)+NoLegend()
    ggsave(fname,p,width = 40, height = 15, limitsize = F)
  })
}
{# white bg 
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
    
    fname = paste0('stp1_summary/plot_whitebg/',x,'_primarytype.png')
    p = xedimplot(dt, groupby = 'ct3',flip = flip,reversex = reversex,reversey = reversey,color = primcol2,bgcol = 'white',pt.size = 0.1)+NoLegend()
    ggsave(fname,p,width = 30, height = 15, limitsize = F)
  })
}

xedimplot(xe1_2, groupby = 'stage',reversex = T,flip = T,color = primcol1)
xedimplot(xe2_1, groupby = 'stage',flip = T,color = primcol1)
xedimplot(xe3_1, groupby = 'stage1',reversex = T,flip = T,color = primcol1)


xe1_2@meta.data = read_rds('data/xenium_seurat_processed/xenium_t1_s2_meta.rds')
newmt = c('cellID','cellType_merge')
newmt = rbind(srttc@meta.data[,newmt],srtfib@meta.data[,newmt])
newmt$cellType_merge = newmt$cellType_merge %>% 
  mapvalues(c('CD8T-n/m','CD4T-n/m'),c('T-naive','T-naive'))
xe1_2$cellType_merge = xe1_2$ct3
xe1_2 = renewMeta(xe1_2, newmt)
p = xedimplot(xe1_2,'cellType_merge',flip = T,reversex = T,color = primcol3,bgcol = 'black',
              ptlevels = rev(names(primcol3)))
ggsave('xe12_grd.png',p,width = 30, height = 15, limitsize = F)
p = xedimplot(xe1_2,'cellType_merge',flip = T,reversex = T,color = primcol3_,bgcol = 'black',
              ptlevels = rev(names(primcol3_)))
ggsave('xe12_dis.png',p,width = 30, height = 15, limitsize = F)




##### spatial stage ####
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
  
  fname = paste0('stp1_summary/stage-spatial_plot/',x,'.png')
  p = xedimplot(dt, groupby = 'stage1',flip = flip,reversex = reversex,reversey = reversey,
                color = col_stage,bgcol = 'black',pt.size = 0.1)+NoLegend()
  ggsave(fname,p,width = 20, height = 10, limitsize = F)
})




##### region of intrested #####
roi = data.frame(xmin=c(900,7500,15000,18500,19000,19800,21500),
                 xmax=c(1900,8500,16000,19500,20000,20800,22500),
                 ymin=c(3000,4000,3500,3300,3300,3300,3200),
                 ymax=c(4000,500,4500,4300,4300,4300,4200))
rownames(roi) = c('T','HT','H','LH','L','NL','N')

xedimplot(xenium1_2, 'primaryType', flip = T, reversex = T, color = primcol,pt.size=0.1)+NoLegend()+
  geom_rect(aes(xmin=900,xmax=2000,ymin=3000,ymax=4000), fill=NA,color='yellow')+
  geom_rect(aes(xmin=15000,xmax=16000,ymin=3800,ymax=4500), fill=NA,color='yellow')+
  geom_rect(aes(xmin=19100,xmax=20000,ymin=3300,ymax=4200), fill=NA,color='yellow')+
  geom_rect(aes(xmin=21000,xmax=22500,ymin=3000,ymax=4200), fill=NA,color='yellow')

ggplot()+
  geom_rect(aes(xmin=c(900,7000,15000,18500,19100,19500,21000),
                xmax=c(2000,9000,16000,19400,20000,21000,22500),
                ymin=c(3000,4000,3800,3300,3300,3300,3000),
                ymax=c(4000,5100,4500,4200,4200,4200,4200)), fill=NA,color='yellow')


{## stage roi
  ## t1
  roi = data.frame(xmin=c(900,7500,15000,18500,19000,19800,21500),
                   xmax=c(1900,8500,16000,19500,20000,20800,22500),
                   ymin=c(3000,4000,3500,3300,3300,3300,3200),
                   ymax=c(4000,500,4500,4300,4300,4300,4200))
  
  p1 = xedimplot(xenium1_2, 'tmp', flip = T, reversex = T, color = ctcol1,pt.size=0.5,zoomin = roi[7,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'NOR Region')
  p2 = xedimplot(xenium1_2, 'tmp', flip = T, reversex = T, color = ctcol1,pt.size=0.5,zoomin = roi[5,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'LGIN Region')
  p3 = xedimplot(xenium1_2, 'tmp', flip = T, reversex = T, color = ctcol1,pt.size=0.5,zoomin = roi[3,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'HGIN Region')
  p4 = xedimplot(xenium1_2, 'tmp', flip = T, reversex = T, color = ctcol1,pt.size=0.5,zoomin = roi[1,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'ESCC Region')
  ggsave('stp1_summary/stage_plot/stage_NOR_1.png',p1,width = 10,height = 8)
  ggsave('stp1_summary/stage_plot/stage_LGIN_1.png',p2,width = 10,height = 8)
  ggsave('stp1_summary/stage_plot/stage_HGIN_1.png',p3,width = 10,height = 8)
  ggsave('stp1_summary/stage_plot/stage_ESCC_1.png',p4,width = 10,height = 8)
  
  ## t2
  roi = matrix(c(7500,8500,1200,2200,
                 14000,15000,2600,3600),ncol = 4,byrow = T)
  p1 = xedimplot(xenium2_1, 'tmp', flip = T, color = ctcol1,pt.size=0.5,zoomin = roi[1,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'NOR Region')
  p2 = xedimplot(xenium2_1, 'tmp', flip = T, color = ctcol1,pt.size=0.5,zoomin = roi[2,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'HGIN Region')
  ggsave('stp1_summary/stage_plot/stage_NOR_2.png',p1,width = 10,height = 8)
  ggsave('stp1_summary/stage_plot/stage_HGIN_2.png',p2,width = 10,height = 8)
  
  ## t3
  roi = matrix(c(16700,17700,2100,3100,
                 13400,14400,1700,2700,
                 4900,5900,3300,4300,
                 5900,6900,600,1600),ncol = 4,byrow = T)
  p1 = xedimplot(xenium3_1, 'tmp', flip = T,reversex = T, color = ctcol1,pt.size=0.5,zoomin = roi[1,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'NOR Region')
  p2 = xedimplot(xenium3_1, 'tmp', flip = T,reversex = T, color = ctcol1,pt.size=0.5,zoomin = roi[2,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'LGIN Region')
  p3 = xedimplot(xenium3_1, 'tmp', flip = T,reversex = T, color = ctcol1,pt.size=0.5,zoomin = roi[3,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'HGIN Region')
  p4 = xedimplot(xenium3_1, 'tmp', flip = T,reversex = T, color = ctcol1,pt.size=0.5,zoomin = roi[4,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'ESCC Region')
  ggsave('stp1_summary/stage_plot/stage_NOR_3.png',p1,width = 10,height = 8)
  ggsave('stp1_summary/stage_plot/stage_LGIN_3.png',p2,width = 10,height = 8)
  ggsave('stp1_summary/stage_plot/stage_HGIN_3.png',p3,width = 10,height = 8)
  ggsave('stp1_summary/stage_plot/stage_ESCC_3.png',p4,width = 10,height = 8)
  
  ## t4
  roi = matrix(c(13500,14500,3800,4800,
                 10800,11800,4000,5000,
                 6800,7800,0,1000),ncol = 4,byrow = T)
  p1 = xedimplot(xenium4_1, 'tmp', flip = T,reversex = T, color = ctcol1,pt.size=0.5,zoomin = roi[1,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'NOR Region')
  p2 = xedimplot(xenium4_1, 'tmp', flip = T,reversex = T, color = ctcol1,pt.size=0.5,zoomin = roi[2,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'HGIN Region')
  p3 = xedimplot(xenium4_1, 'tmp', flip = T,reversex = T, color = ctcol1,pt.size=0.5,zoomin = roi[3,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'ESCC Region')
  ggsave('stp1_summary/stage_plot/stage_NOR_4.png',p1,width = 10,height = 8)
  ggsave('stp1_summary/stage_plot/stage_HGIN_4.png',p2,width = 10,height = 8)
  ggsave('stp1_summary/stage_plot/stage_ESCC_4.png',p3,width = 10,height = 8)
  
  ## t5
  roi = matrix(c(1700,2700,1000,2000,
                 3500,4500,3200,4200),ncol = 4,byrow = T)
  p1 = xedimplot(xenium5_1, 'tmp', reversey = T, color = ctcol1,pt.size=0.5,zoomin = roi[1,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'LGIN Region')
  p2 = xedimplot(xenium5_1, 'tmp', reversey = F, color = ctcol1,pt.size=0.5,zoomin = roi[2,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'HGIN Region')
  ggsave('stp1_summary/stage_plot/stage_LGIN_5.png',p1,width = 10,height = 8)
  ggsave('stp1_summary/stage_plot/stage_HGIN_5.png',p2,width = 10,height = 8)
  
  ## t6
  roi = matrix(c(7500,8500,1700,2700,
                 2800,3800,1600,2600,
                 3800,4800,1400,2400),ncol = 4,byrow = T)
  p1 = xedimplot(xenium6_1, 'tmp', reversey = T, color = ctcol1,pt.size=0.5,zoomin = roi[1,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'NOR Region')
  p2 = xedimplot(xenium6_1, 'tmp', reversey = T, color = ctcol1,pt.size=0.5,zoomin = roi[2,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'HGIN Region')
  p3 = xedimplot(xenium6_1, 'tmp', reversey = T, color = ctcol1,pt.size=0.5,zoomin = roi[3,],
                 cellborder = T,areaalpha = 0.3,linewidth = 0.2)+
    theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+xlab('')+ylab('')+
    NoLegend()+ggtitle(label = 'ESCC Region')
  ggsave('stp1_summary/stage_plot/stage_NOR_6.png',p1,width = 10,height = 8)
  ggsave('stp1_summary/stage_plot/stage_HGIN_6.png',p2,width = 10,height = 8)
  ggsave('stp1_summary/stage_plot/stage_ESCC_6.png',p3,width = 10,height = 8)
}




### zoomin epi 
xenium1_2$tmp = xenium1_2$cellType_merge
xenium1_2$tmp[xenium1_2$primaryType!='Epi'] = xenium1_2$primaryType[xenium1_2$primaryType!='Epi'] %>% as.character()

epicol = c(brewer.pal(7,'Set3'),brewer.pal(9,'Set1'))
names(epicol) = c(unique(xenium1_2$tmp[xenium1_2$primaryType=='Epi']),unique(xenium1_2$tmp[xenium1_2$primaryType!='Epi']))

xenium1_2$tmp = xenium1_2$tmp %>% 
  factor(levels = names(epicol))



xedimplot(xenium1_2, 'tmp', flip = T, reversex = T, color = epicol,pt.size=0.5,zoomin = roi[7,])+ggtitle(label = 'NOR Region')
xedimplot(xenium1_2, 'tmp', flip = T, reversex = T, color = epicol,pt.size=0.5,zoomin = roi[6,])+ggtitle(label = 'NOR-LGIN')
xedimplot(xenium1_2, 'tmp', flip = T, reversex = T, color = epicol,pt.size=0.5,zoomin = roi[5,])+ggtitle(label = 'LGIN Region')
xedimplot(xenium1_2, 'tmp', flip = T, reversex = T, color = epicol,pt.size=0.5,zoomin = roi[3,])+ggtitle(label = 'HGIN Region')
xedimplot(xenium1_2, 'tmp', flip = T, reversex = T, color = epicol,pt.size=0.5,zoomin = roi[1,])+ggtitle(label = 'ESCC Region')




##### umi violin #####
VlnPlot(srt, features = c('nCount_RNA','nFeature_RNA'),pt.size = 0,
        group.by = 'sp_stg',ncol = 1)




##### marker dotplot #####
DotPlot(srt, features = allmk2, group.by = 'ct3')+RotatedAxis()
mydotp(srt,allmk2,'ct3',rev(brewer.pal(9,'RdBu')),lims = c(-2,2))+RotatedAxis()
ggsave('stp1_summary/summary_plot/marker_dotp.pdf',width = 5,height = 8)

srt$ct3 = srt$ct3 %>% factor(levels = c('Epi','Fib','T','B','Endo','Mye','Mus','Gland'))
mydotp(srt,allmk2,'ct3',rev(brewer.pal(9,'RdBu')),lims = c(-2,2))+rotate_x_text()
ggsave('stp1_summary/summary_plot/marker_dotp.pdf',width = 5,height = 8)

srt$ct3 = srt$ct3 %>% factor(levels = rev(c('Epi','Fib','T','B','Endo','Mye','Mus','Gland')))
DotPlot(srt, features = allmk2, group.by = 'ct3',cols = 'RdBu')+rotate_x_text()
ggsave('stp1_summary/summary_plot/marker_dotp_2.pdf',width = 10,height = 3.5)


newmt = c('cellID','cellType_merge')
newmt = rbind(srttc@meta.data[,newmt],srtfib@meta.data[,newmt])
newmt$cellType_merge = newmt$cellType_merge %>% 
  mapvalues(c('CD8T-n/m','CD4T-n/m'),c('T-naive','T-naive'))

srt$cellType_merge = srt$ct3 %>% as.character()
srt@meta.data[newmt$cellID,'cellType_merge'] = newmt$cellType_merge
srt$cellType_merge = srt$cellType_merge %>% 
  factor(levels = c('Epi','Endo','NF','NF-PI16','CAF','T-others',
                    'T-naive','CD8T-eff','CD8T-ex','CD4T-fh','CD4T-reg',
                    'B','Mye','Gland','Mus'))
srt$cellType_merge = srt$cellType_merge %>% 
  mapvalues(c('Epi','Endo','B','Mye','Gland','Mus'),
            c('Epithelial cells','Endothelial cells ','B cells ','Myeloid cells','Gland cells','Myocytes'))
p = mydotp(srt,allmk3,'cellType_merge',rev(brewer.pal(9,'RdBu')),lims = c(-2,2))+RotatedAxis()
ggsave('stp1_summary/summary_plot/marker_dotp_subtype.pdf',width = 7,height = 14)






