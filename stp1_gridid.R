source('~/xenium/header.R')
source('~/xenium/xenium_function.R')

##### input #####
allmeta = read_rds('data/xenium_seurat_processed/allcell_meta.rds')
crds = read_rds('stp1_summary/data/allsample_coords.rds')

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





##### grid id #####
gridlength = 100
gridid=lapply(names(crds), function(x){
  crd = crds[[x]]
  
  idbyspstg = lapply(unique(crd$sp_stg), function(y){
    dt = crd[crd$sp_stg==y,]
    gridx = dt$x%/%gridlength + 1
    gridy = dt$y%/%gridlength + 1
    gridid = paste0(x,'-',gridx,'-',gridy)
    
    dt$gridid = gridid %>% factor(levels = sort(unique(gridid)))
    return(dt)
  }) %>% Reduce(rbind,.)
  
  idbyspstg$slideid = x
  return(idbyspstg)
}) %>% Reduce(rbind,.)
saveRDS(gridid,'stp1_summary/data/allsample_gridid.rds')  






