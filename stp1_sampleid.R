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

probelist = fread('data/235gene_label.txt',data.table = F)
# probelist = fread('data/probe_list.txt',data.table = F)


##### color #####
primcol = brewer.pal(10, 'Paired')
primcol[7] = '#999999'
names(primcol) = c('B','Epi','Endo','Fib','Gland','Mye','Mus','Mast','Pla','T')
primcol2 = c('#1f78b4','#ffc089','#B5EFB5','#793c1b',
             '#6a3d9a','#333333','#ffff33',
             '#df928e')
names(primcol2) = c('B','Endo','Fib','Gland','Mye','Mus','T','Epi')


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




##### sample-stage id #####
getSamplestage = function(sampleID){
  sp1 = paste0('xe',sampleID,'_1')
  sp2 = paste0('xe',sampleID,'_2')
  sp3 = paste0('xe',sampleID,'_3')
  
  sp1 = get(sp1)@meta.data[,c('cellID','stage','primaryType')] %>% mutate(ID=paste0('T',sampleID,'S1'))
  sp2 = get(sp2)@meta.data[,c('cellID','stage','primaryType')] %>% mutate(ID=paste0('T',sampleID,'S2'))
  if(exists(sp3)){
    sp3 = get(sp3)@meta.data[,c('cellID','stage','primaryType')] %>% mutate(ID=paste0('T',sampleID,'S3'))
  }else{sp3 = NULL}
  allmeta = Reduce(rbind, list(sp1,sp2,sp3))
  allmeta$stage1 = allmeta$stage %>% str_split(fixed('-'),simplify = T) %>% .[,1]
  allmeta$ct3 = allmeta$primaryType %>% mapvalues(c('Mast','Pla'),c('Mye','B'))
  allmeta = allmeta %>% arrange(ID,stage)

  nmeta = allmeta %>% filter(stage1=='NOR') %>% select(ID,stage,stage1) %>% unique %>% mutate(num=row_number(.))
  lmeta = allmeta %>% filter(stage1=='LGIN') %>% select(ID,stage,stage1) %>% unique %>% mutate(num=row_number(.))
  hmeta = allmeta %>% filter(stage1=='HGIN') %>% select(ID,stage,stage1) %>% unique %>% mutate(num=row_number(.))
  emeta = allmeta %>% filter(stage1=='ESCC') %>% select(ID,stage,stage1) %>% unique %>% mutate(num=row_number(.))
  
  spstg = Reduce(rbind,list(nmeta,lmeta,hmeta,emeta))
  spstg$sp_stg = paste0('T',sampleID,'-',spstg$stage1,'-',spstg$num)
  allmeta = spstg %>% dplyr::select(-num) %>% left_join(allmeta,.)
  rownames(allmeta) = allmeta$cellID
  return(allmeta[,c('cellID','stage1','sp_stg','ct3')])
}

{
  spstg1 = getSamplestage(1)
  spstg2 = getSamplestage(2)
  spstg3 = getSamplestage(3)
  spstg4 = getSamplestage(4)
  spstg5 = getSamplestage(5)
  spstg6 = getSamplestage(6)
  
  xe1_1@meta.data = xe1_1@meta.data %>% dplyr::select(-contains(colnames(spstg1))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg1)
  rownames(xe1_1@meta.data) = xe1_1$cellID
  xe1_2@meta.data = xe1_2@meta.data %>% dplyr::select(-contains(colnames(spstg1))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg1)
  rownames(xe1_2@meta.data) = xe1_2$cellID
  xe1_3@meta.data = xe1_3@meta.data %>% dplyr::select(-contains(colnames(spstg1))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg1)
  rownames(xe1_3@meta.data) = xe1_3$cellID
  
  xe2_1@meta.data = xe2_1@meta.data %>% dplyr::select(-contains(colnames(spstg2))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg2)
  rownames(xe2_1@meta.data) = xe2_1$cellID
  xe2_2@meta.data = xe2_2@meta.data %>% dplyr::select(-contains(colnames(spstg2))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg2)
  rownames(xe2_2@meta.data) = xe2_2$cellID
  xe2_3@meta.data = xe2_3@meta.data %>% dplyr::select(-contains(colnames(spstg2))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg2)
  rownames(xe2_3@meta.data) = xe2_3$cellID
  
  xe3_1@meta.data = xe3_1@meta.data %>% dplyr::select(-contains(colnames(spstg3))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg3)
  rownames(xe3_1@meta.data) = xe3_1$cellID
  xe3_2@meta.data = xe3_2@meta.data %>% dplyr::select(-contains(colnames(spstg3))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg3)
  rownames(xe3_2@meta.data) = xe3_2$cellID

  xe4_1@meta.data = xe4_1@meta.data %>% dplyr::select(-contains(colnames(spstg4))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg4)
  rownames(xe4_1@meta.data) = xe4_1$cellID
  xe4_2@meta.data = xe4_2@meta.data %>% dplyr::select(-contains(colnames(spstg4))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg4)
  rownames(xe4_2@meta.data) = xe4_2$cellID

  xe5_1@meta.data = xe5_1@meta.data %>% dplyr::select(-contains(colnames(spstg5))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg5)
  rownames(xe5_1@meta.data) = xe5_1$cellID
  xe5_2@meta.data = xe5_2@meta.data %>% dplyr::select(-contains(colnames(spstg5))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg5)
  rownames(xe5_2@meta.data) = xe5_2$cellID
  xe5_3@meta.data = xe5_3@meta.data %>% dplyr::select(-contains(colnames(spstg5))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg5)
  rownames(xe5_3@meta.data) = xe5_3$cellID
  xe5_2_split@meta.data = xe5_2_split@meta.data %>% dplyr::select(-contains(colnames(spstg5))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg5)
  rownames(xe5_2_split@meta.data) = xe5_2_split$cellID
  xe5_3_split@meta.data = xe5_3_split@meta.data %>% dplyr::select(-contains(colnames(spstg5))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg5)
  rownames(xe5_3_split@meta.data) = xe5_3_split$cellID

  xe6_1@meta.data = xe6_1@meta.data %>% dplyr::select(-contains(colnames(spstg6))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg6)
  rownames(xe6_1@meta.data) = xe6_1$cellID
  xe6_2@meta.data = xe6_2@meta.data %>% dplyr::select(-contains(colnames(spstg6))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg6)
  rownames(xe6_2@meta.data) = xe6_2$cellID
  xe6_3@meta.data = xe6_3@meta.data %>% dplyr::select(-contains(colnames(spstg6))) %>% 
    mutate(cellID=rownames(.)) %>% left_join(spstg6)
  rownames(xe6_3@meta.data) = xe6_3$cellID
  
  saveRDS(xe1_1@meta.data, 'data/xenium_seurat_processed/xenium_t1_s1_meta.rds')
  saveRDS(xe1_2@meta.data, 'data/xenium_seurat_processed/xenium_t1_s2_meta.rds')
  saveRDS(xe1_3@meta.data, 'data/xenium_seurat_processed/xenium_t1_s3_meta.rds')
  saveRDS(xe2_1@meta.data, 'data/xenium_seurat_processed/xenium_t2_s1_meta.rds')
  saveRDS(xe2_2@meta.data, 'data/xenium_seurat_processed/xenium_t2_s2_meta.rds')
  saveRDS(xe2_3@meta.data, 'data/xenium_seurat_processed/xenium_t2_s3_meta.rds')
  saveRDS(xe3_1@meta.data, 'data/xenium_seurat_processed/xenium_t3_s1_meta.rds')
  saveRDS(xe3_2@meta.data, 'data/xenium_seurat_processed/xenium_t3_s2_meta.rds')
  saveRDS(xe4_1@meta.data, 'data/xenium_seurat_processed/xenium_t4_s1_meta.rds')
  saveRDS(xe4_2@meta.data, 'data/xenium_seurat_processed/xenium_t4_s2_meta.rds')
  saveRDS(xe5_1@meta.data, 'data/xenium_seurat_processed/xenium_t5_s1_meta.rds')
  saveRDS(xe5_2@meta.data, 'data/xenium_seurat_processed/xenium_t5_s2_meta.rds')
  saveRDS(xe5_3@meta.data, 'data/xenium_seurat_processed/xenium_t5_s3_meta.rds')
  saveRDS(xe5_2_split@meta.data, 'data/xenium_seurat_processed/xenium_t5_s2_split_meta.rds')
  saveRDS(xe5_3_split@meta.data, 'data/xenium_seurat_processed/xenium_t5_s3_split_meta.rds')
  saveRDS(xe6_1@meta.data, 'data/xenium_seurat_processed/xenium_t6_s1_meta.rds')
  saveRDS(xe6_2@meta.data, 'data/xenium_seurat_processed/xenium_t6_s2_meta.rds')
  saveRDS(xe6_3@meta.data, 'data/xenium_seurat_processed/xenium_t6_s3_meta.rds')
}

{## allmeta
  metalist = list(xe1_1@meta.data,xe1_2@meta.data,xe1_3@meta.data,
                  xe2_1@meta.data,xe2_2@meta.data,xe2_3@meta.data,
                  xe3_1@meta.data,xe3_2@meta.data,
                  xe4_1@meta.data,xe4_2@meta.data,
                  xe5_1@meta.data,xe5_2@meta.data,xe5_3@meta.data,
                  xe6_1@meta.data,xe6_2@meta.data,xe6_3@meta.data)
  clm = lapply(metalist,colnames) %>% Reduce(intersect,.)
  allmeta = lapply(metalist,function(x){return(x[,clm])}) %>% Reduce(rbind,.)
  allmeta$stage1 = allmeta$stage1 %>% factor(levels = c('NOR','LGIN','HGIN','ESCC'))
  saveRDS(allmeta,'data/xenium_seurat_processed/allmeta.rds')
  
  ord = allmeta[,c('sp_stg','cellID','stage','stage1')]
  ord$sample = ord$cellID %>% str_split('_',simplify = T) %>% .[,1]
  ord = ord %>% select(-cellID) %>% unique()
  ord = ord %>% arrange(stage1,sample,stage)
  
  allmeta$sp_stg = allmeta$sp_stg %>% factor(levels = ord$sp_stg)
  allmeta$pid = allmeta$sp_stg %>% str_split(fixed('-'),simplify = T) %>% .[,1]
  allmeta$sid = allmeta$cellID %>% str_split('_',simplify = T) %>% .[,1]
  saveRDS(allmeta,'data/xenium_seurat_processed/allmeta.rds')
}




