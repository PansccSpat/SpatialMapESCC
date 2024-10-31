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

srttc = read_rds('data/xenium_seurat_processed/ct_t_srt.rds')



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



##### epi area #####
areap = function(coords,ap){
  ap = ap[rownames(coords),]
  coords = cbind(coords,ap)
  
  coords %>% 
    ggplot(aes(x,y))+
    theme1+
    geom_point(aes(color=type))
}

{# t1s1
  tmp = subset(xe1_1, cells = xe1_1$cellID[xe1_1$ct3=='Epi'])
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)
  colnames(crd) = c('y','x')

  crd_tmp = crd %>% filter(x<=3000)
  plotly::plot_ly(crd_tmp,x=~x,y=~y)
  cp = c('T1S1_clgbikjg-1','T1S1_mkmhefia-1','T1S1_kcaoffke-1','T1S1_mlocnknh-1')
  job({ap_tmp = areaPoint(dist(crd_tmp),cpoint = cp,r = 50)},import='auto')
  areap(crd_tmp,ap_tmp)
  id1 = ap_tmp$id[ap_tmp$type!='isolated']
  crd_tmp = crd_tmp[id1,c('x','y')]
  as_tmp = alphahull::ashape(crd_tmp,alpha = 100)
  plot(as_tmp)
  poly1 = extractpoly(as_tmp)
  plotpoly(poly1)
  
  crd_tmp = crd %>% filter(x<=6000,x>2900)
  plotly::plot_ly(crd_tmp,x=~x,y=~y)
  cp = c('T1S1_jajnkhib-1','T1S1_ofjhckbd-1')
  job({ap_tmp = areaPoint(dist(crd_tmp),cpoint = cp,r = 30)},import='auto')
  areap(crd_tmp,ap_tmp)
  id1 = ap_tmp$id[ap_tmp$type!='isolated']
  crd_tmp = crd_tmp[id1,c('x','y')]
  as_tmp = alphahull::ashape(crd_tmp,alpha = 100)
  plot(as_tmp)
  poly2 = extractpoly(as_tmp)
  plotpoly(poly2)
  
  crd_tmp = crd %>% filter(x<=9000,x>5900)
  plotly::plot_ly(crd_tmp,x=~x,y=~y)
  cp = c('T1S1_gaenjdkg-1','T1S1_ofjienlh-1')
  job({ap_tmp = areaPoint(dist(crd_tmp),cpoint = cp,r = 30)},import='auto')
  areap(crd_tmp,ap_tmp)
  id1 = ap_tmp$id[ap_tmp$type!='isolated']
  crd_tmp = crd_tmp[id1,c('x','y')]
  as_tmp = alphahull::ashape(crd_tmp,alpha = 100)
  plot(as_tmp)
  poly3 = extractpoly(as_tmp)
  plotpoly(poly3)
  
  crd_tmp = crd %>% filter(x<=12000,x>8900)
  plotly::plot_ly(crd_tmp,x=~x,y=~y)
  cp = c('T1S1_jgnbndfl-1','T1S1_abgljkme-1')
  job({ap_tmp = areaPoint(dist(crd_tmp),cpoint = cp,r = 30)},import='auto')
  areap(crd_tmp,ap_tmp)
  id1 = ap_tmp$id[ap_tmp$type!='isolated']
  crd_tmp = crd_tmp[id1,c('x','y')]
  as_tmp = alphahull::ashape(crd_tmp,alpha = 100)
  plot(as_tmp)
  poly4 = extractpoly(as_tmp)
  plotpoly(poly4)
  
  crd_tmp = crd %>% filter(x<=15000,x>11900)
  plotly::plot_ly(crd_tmp,x=~x,y=~y)
  cp = c('T1S1_lecadonb-1')
  job({ap_tmp = areaPoint(dist(crd_tmp),cpoint = cp,r = 30)},import='auto')
  areap(crd_tmp,ap_tmp)
  id1 = ap_tmp$id[ap_tmp$type!='isolated']
  crd_tmp = crd_tmp[id1,c('x','y')]
  as_tmp = alphahull::ashape(crd_tmp,alpha = 100)
  plot(as_tmp)
  poly5 = extractpoly(as_tmp)
  plotpoly(poly5)
  
  crd_tmp = crd %>% filter(x<=18000,x>14900)
  plotly::plot_ly(crd_tmp,x=~x,y=~y)
  cp = c('T1S1_ofnedcgj-1','T1S1_oikjdmkd-1')
  job({ap_tmp = areaPoint(dist(crd_tmp),cpoint = cp,r = 30)},import='auto')
  areap(crd_tmp,ap_tmp)
  id1 = ap_tmp$id[ap_tmp$type!='isolated']
  crd_tmp = crd_tmp[id1,c('x','y')]
  as_tmp = alphahull::ashape(crd_tmp,alpha = 100)
  plot(as_tmp)
  poly6 = extractpoly(as_tmp)
  plotpoly(poly6)
  
  crd_tmp = crd %>% filter(x>17900)
  plotly::plot_ly(crd_tmp,x=~x,y=~y)
  cp = c('T1S1_kgfolccj-1')
  job({ap_tmp = areaPoint(dist(crd_tmp),cpoint = cp,r = 30)},import='auto')
  areap(crd_tmp,ap_tmp)
  id1 = ap_tmp$id[ap_tmp$type!='isolated']
  crd_tmp = crd_tmp[id1,c('x','y')]
  as_tmp = alphahull::ashape(crd_tmp,alpha = 100)
  plot(as_tmp)
  poly7 = extractpoly(as_tmp)
  plotpoly(poly7)
  
  poly_11 = Reduce(append,list(poly1,poly2,poly3,poly4,poly5,poly6,poly7))
  poly_11 = poly_11[-c(2,4)]
  plotpoly(poly_11)
  saveRDS(poly_11,'stp3_tme/data/epiregion/T1S1_poly.rds')
}
{# t1s2
  tmp = subset(xe1_2, cells = xe1_2$cellID[xe1_2$ct3=='Epi'])
  xedimplot(tmp,groupby = 'ct3',flip = T,reversex = T)
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)
  colnames(crd) = c('y','x')
  
  getcoordsrange(tmp,tmp$cellID[tmp$stage1=='ESCC'],flip = T)
  crd_escc = crd %>% filter(x<=7800)
  plotly::plot_ly(crd_escc,x=~x,y=~y)
  cp = c('T1S2_fiaehfid-1','T1S2_bndohmei-1')
  job({ap_escc = areaPoint(dist(crd_escc),cpoint = cp,r = 50)},import='auto')
  areap(crd_escc,ap_escc)
  id1 = ap_escc$id[ap_escc$type!='isolated']
  crd_escc = crd_escc[id1,c('x','y')]
  as_escc = alphahull::ashape(crd_escc,alpha = 100)
  plot(as_escc)
  polyescc = extractpoly(as_escc)
  plotpoly(polyescc)
  
  getcoordsrange(tmp,tmp$cellID[tmp$stage1=='HGIN'],flip = T)
  crd_hgin = crd %>% filter(x<=19000,x>7700)
  plotly::plot_ly(crd_hgin,x=~x,y=~y)
  cp = c('T1S2_kpfhfmmm-1','T1S2_eknacomk-1','T1S2_hcbfnfgk-1')
  job({ap_hgin = areaPoint(dist(crd_hgin),cpoint = cp,r = 50)},import='auto')
  areap(crd_hgin,ap_hgin)
  id1 = ap_hgin$id[ap_hgin$type!='isolated']
  crd_hgin = crd_hgin[id1,c('x','y')]
  as_hgin = alphahull::ashape(crd_hgin,alpha = 80)
  plot(as_hgin)
  polyhgin = extractpoly(as_hgin)
  plotpoly(polyhgin)
  
  getcoordsrange(tmp,tmp$cellID[tmp$stage1%in%c('LGIN','NOR')],flip = T)
  crd_ln = crd %>% filter(x>18900)
  plotly::plot_ly(crd_ln,x=~x,y=~y)
  cp = c('T1S2_lccialec-1')
  job({ap_ln = areaPoint(dist(crd_ln),cpoint = cp,r = 50)},import='auto')
  areap(crd_ln,ap_ln)
  id1 = ap_ln$id[ap_ln$type!='isolated']
  crd_ln = crd_ln[id1,c('x','y')]
  as_ln = alphahull::ashape(crd_ln,alpha = 100)
  plot(as_ln)
  polyln = extractpoly(as_ln)
  plotpoly(polyln)

  poly_12 = Reduce(append,list(polyescc,polyhgin,polyln))
  poly_12 = poly_12[-2]
  plotpoly(poly_12)
  saveRDS(poly_12,'stp3_tme/data/epiregion/T1S2_poly.rds')
}
{# t1s3
  tmp = subset(xe1_3, cells = xe1_3$cellID[xe1_3$ct3=='Epi'])
  xedimplot(tmp,groupby = 'stage',flip = T,reversex = T)
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)
  colnames(crd) = c('y','x')
  
  getcoordsrange(tmp,tmp$cellID[tmp$stage1=='ESCC'],flip = T)
  crd_escc = crd %>% filter(x<=8050)
  plotly::plot_ly(crd_escc,x=~x,y=~y)
  cp = c('T1S3_mejdelan-1','T1S3_hjepljbp-1','T1S3_chmccjfm-1','T1S3_cgmfdafb-1')
  job({ap_escc = areaPoint(dist(crd_escc),cpoint = cp,r = 50)},import='auto')
  areap(crd_escc,ap_escc)
  id1 = ap_escc$id[ap_escc$type!='isolated']
  crd_escc = crd_escc[id1,c('x','y')]
  as_escc = alphahull::ashape(crd_escc,alpha = 100)
  plot(as_escc)
  polyescc = extractpoly(as_escc)
  plotpoly(polyescc)
  
  getcoordsrange(tmp,tmp$cellID[tmp$stage1=='HGIN'],flip = T)
  crd_hgin = crd %>% filter(x<=19327,x>7950)
  plotly::plot_ly(crd_hgin,x=~x,y=~y)
  cp = c('T1S3_kgfeokol-1','T1S3_bielbglm-1')
  job({ap_hgin = areaPoint(dist(crd_hgin),cpoint = cp,r = 50)},import='auto')
  areap(crd_hgin,ap_hgin)
  id1 = ap_hgin$id[ap_hgin$type!='isolated']
  crd_hgin = crd_hgin[id1,c('x','y')]
  as_hgin = alphahull::ashape(crd_hgin,alpha = 100)
  plot(as_hgin)
  polyhgin = extractpoly(as_hgin)
  plotpoly(polyhgin)
  
  getcoordsrange(tmp,tmp$cellID[tmp$stage1%in%c('LGIN','NOR')],flip = T)
  crd_ln = crd %>% filter(x>19227)
  plotly::plot_ly(crd_ln,x=~x,y=~y)
  cp = c('T1S3_gkjognbf-1','T1S3_iekjbggf-1')
  job({ap_ln = areaPoint(dist(crd_ln),cpoint = cp,r = 50)},import='auto')
  areap(crd_ln,ap_ln)
  id1 = ap_ln$id[ap_ln$type!='isolated']
  crd_ln = crd_ln[id1,c('x','y')]
  as_ln = alphahull::ashape(crd_ln,alpha = 80)
  plot(as_ln)
  polyln = extractpoly(as_ln)
  plotpoly(polyln)

  poly_13 = Reduce(append,list(polyescc,polyhgin,polyln))
  poly_13 = poly_13[-c(2:3)]
  plotpoly(poly_13)
  saveRDS(poly_13,'stp3_tme/data/epiregion/T1S3_poly.rds')
}
{# t2s1
  tmp = subset(xe2_1, cells = xe2_1$cellID[xe2_1$ct3=='Epi'])
  xedimplot(tmp,groupby = 'stage',flip = T)
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)
  colnames(crd) = c('y','x')
  
  plotly::plot_ly(crd,x=~x,y=~y)
  cp = c('T2S1_jkfchill-1')
  job({ap = areaPoint(dist(crd),cpoint = cp,r = 50)},import='auto')
  areap(crd,ap)
  id1 = ap$id[ap$type!='isolated']
  crd = crd[id1,c('x','y')]
  as = alphahull::ashape(crd,alpha = 100)
  plot(as)
  polytmp = extractpoly(as)
  plotpoly(polytmp)
  
  poly_21 = polytmp
  plotpoly(poly_21)
  saveRDS(poly_21,'stp3_tme/data/epiregion/T2S1_poly.rds')
}
{# t2s2
  tmp = subset(xe2_2, cells = xe2_2$cellID[xe2_2$ct3=='Epi'])
  xedimplot(tmp,groupby = 'stage',flip = T)
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)
  colnames(crd) = c('y','x')
  
  plotly::plot_ly(crd,x=~x,y=~y)
  cp = c('T2S2_hcjginok-1')
  job({ap = areaPoint(dist(crd),cpoint = cp,r = 50)},import='auto')
  areap(crd,ap)
  id1 = ap$id[ap$type!='isolated']
  crd = crd[id1,c('x','y')]
  as = alphahull::ashape(crd,alpha = 100)
  plot(as)
  polytmp = extractpoly(as)
  plotpoly(polytmp)
  
  poly_22 = polytmp
  plotpoly(poly_22)
  saveRDS(poly_22,'stp3_tme/data/epiregion/T2S2_poly.rds')
}
{# t2s3
  tmp = subset(xe2_3, cells = xe2_3$cellID[xe2_3$ct3=='Epi'])
  xedimplot(tmp,groupby = 'stage',flip = T)
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)
  colnames(crd) = c('y','x')
  
  plotly::plot_ly(crd,x=~x,y=~y)
  cp = c('T2S3_ochlpbdm-1')
  job({ap = areaPoint(dist(crd),cpoint = cp,r = 50)},import='auto')
  areap(crd,ap)
  id1 = ap$id[ap$type!='isolated']
  crd = crd[id1,c('x','y')]
  as = alphahull::ashape(crd,alpha = 100)
  plot(as)
  polytmp = extractpoly(as)
  plotpoly(polytmp)
  
  poly_23 = polytmp
  plotpoly(poly_23)
  saveRDS(poly_23,'stp3_tme/data/epiregion/T2S3_poly.rds')
}
{# t3s1
  tmp = subset(xe3_1, cells = xe3_1$cellID[xe3_1$ct3=='Epi'])
  xedimplot(tmp,groupby = 'stage',flip = T,reversex = T)
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)
  colnames(crd) = c('y','x')
  
  crd_escc = crd[tmp$cellID[tmp$stage1=='ESCC'],]
  plotly::plot_ly(crd_escc,x=~x,y=~y)
  crd_escc = crd_escc[,c('x','y')]
  as_escc = alphahull::ashape(crd_escc,alpha = 200)
  plot(as_escc)
  polyescc = extractpoly(as_escc)
  plotpoly(polyescc)
  
  crd_hln = crd[tmp$cellID[tmp$stage1!='ESCC'],]
  plotly::plot_ly(crd_hln,x=~x,y=~y)
  cp = c('T3S1_clmiefdi-1')
  job({ap_hln = areaPoint(dist(crd_hln),cpoint = cp,r = 50)},import='auto')
  areap(crd_hln,ap_hln)
  id1 = ap_hln$id[ap_hln$type!='isolated']
  crd_hln = crd_hln[id1,c('x','y')]
  as_hln = alphahull::ashape(crd_hln,alpha = 100)
  plot(as_hln)
  polyhln = extractpoly(as_hln)
  plotpoly(polyhln)
  polyhln = polyhln[-2]
  
  poly_31 = Reduce(append,list(polyescc,polyhln))
  plotpoly(poly_31)
  saveRDS(poly_31,'stp3_tme/data/epiregion/T3S1_poly.rds')
}
{# t3s2
  tmp = subset(xe3_2, cells = xe3_2$cellID[xe3_2$ct3=='Epi'])
  xedimplot(tmp,groupby = 'stage',flip = T,reversex = T)
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)
  colnames(crd) = c('y','x')
  
  crd_escc = crd[tmp$cellID[tmp$stage1=='ESCC'],]
  plotly::plot_ly(crd_escc,x=~x,y=~y)
  crd_escc = crd_escc[,c('x','y')]
  as_escc = alphahull::ashape(crd_escc,alpha = 200)
  plot(as_escc)
  polyescc = extractpoly(as_escc)
  plotpoly(polyescc)
  
  crd_hln = crd[tmp$cellID[tmp$stage1!='ESCC'],]
  plotly::plot_ly(crd_hln,x=~x,y=~y)
  cp = c('T3S2_dkbnfmna-1')
  job({ap_hln = areaPoint(dist(crd_hln),cpoint = cp,r = 50)},import='auto')
  areap(crd_hln,ap_hln)
  id1 = ap_hln$id[ap_hln$type!='isolated']
  crd_hln = crd_hln[id1,c('x','y')]
  as_hln = alphahull::ashape(crd_hln,alpha = 150)
  plot(as_hln)
  polyhln = extractpoly(as_hln)
  plotpoly(polyhln)
  
  poly_32 = Reduce(append,list(polyescc,polyhln))
  plotpoly(poly_32)
  saveRDS(poly_32,'stp3_tme/data/epiregion/T3S2_poly.rds')
}
{# t4s1
  tmp = subset(xe4_1, cells = xe4_1$cellID[xe4_1$ct3=='Epi'])
  xedimplot(tmp,groupby = 'stage',flip = T,reversex = T)
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)
  colnames(crd) = c('y','x')
  
  crd_escc1 = crd[tmp$cellID[tmp$stage=='ESCC-1'],]
  plotly::plot_ly(crd_escc1,x=~x,y=~y)
  cp = c('T4S1_fkmgnkeb-1','T4S1_didhijpk-1','T4S1_nhfnhcfm-1')
  job({ap_escc1 = areaPoint(dist(crd_escc1),cpoint = cp,r = 50)},import='auto')
  areap(crd_escc1,ap_escc1)
  id1 = ap_escc1$id[ap_escc1$type!='isolated']
  crd_escc1 = crd_escc1[id1,c('x','y')]
  as_escc1 = alphahull::ashape(crd_escc1,alpha = 100)
  plot(as_escc1)
  polyescc1 = extractpoly(as_escc1)
  plotpoly(polyescc1)
  polyescc1 = polyescc1[-2]
  
  crd_escc2 = crd[tmp$cellID[tmp$stage=='ESCC-2'],]
  plotly::plot_ly(crd_escc2,x=~x,y=~y)
  cp = c('T4S1_jajcljcm-1','T4S1_hbhjbjjk-1')
  job({ap_escc2 = areaPoint(dist(crd_escc2),cpoint = cp,r = 50)},import='auto')
  areap(crd_escc2,ap_escc2)
  id1 = ap_escc2$id[ap_escc2$type!='isolated']
  crd_escc2 = crd_escc2[id1,c('x','y')]
  as_escc2 = alphahull::ashape(crd_escc2,alpha = 100)
  plot(as_escc2)
  polyescc2 = extractpoly(as_escc2)
  plotpoly(polyescc2)
  polyescc2 = polyescc2[-2]
  
  crd_escc3 = crd[tmp$cellID[tmp$stage=='ESCC-3'],]
  crd_escc3 = crd_escc3[,c('x','y')]
  as_escc3 = alphahull::ashape(crd_escc3,alpha = 300)
  plot(as_escc3)
  polyescc3 = extractpoly(as_escc3)
  plotpoly(polyescc3)
  
  crd_hln = crd[tmp$cellID[tmp$stage1!='ESCC'],]
  plotly::plot_ly(crd_hln,x=~x,y=~y)
  cp = c('T4S1_gpkmnace-1')
  job({ap_hln = areaPoint(dist(crd_hln),cpoint = cp,r = 50)},import='auto')
  areap(crd_hln,ap_hln)
  id1 = ap_hln$id[ap_hln$type!='isolated']
  crd_hln = crd_hln[id1,c('x','y')]
  as_hln = alphahull::ashape(crd_hln,alpha = 100)
  plot(as_hln)
  polyhln = extractpoly(as_hln)
  plotpoly(polyhln)
  
  poly_41 = Reduce(append,list(polyescc1,polyescc2,polyescc3,polyhln))
  plotpoly(poly_41)
  saveRDS(poly_41,'stp3_tme/data/epiregion/T4S1_poly.rds')
}
{# t4s2
  tmp = subset(xe4_2, cells = xe4_2$cellID[xe4_2$ct3=='Epi'])
  xedimplot(tmp,groupby = 'stage',flip = T,reversey = T)
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)
  colnames(crd) = c('y','x')
  
  crd_escc1 = crd[tmp$cellID[tmp$stage=='ESCC-1'],]
  (crd_escc1 %>% mutate(grp=rownames(.)) %>% ggplot(aes(x,y,group=grp))+geom_point(size=0.5)) %>% plotly::ggplotly()
  cp = c('T4S2_iaoeikje-1','T4S2_lmamdlih-1','T4S2_cdoanplb-1','T4S2_hcoldbkn-1')
  job({ap_escc1 = areaPoint(dist(crd_escc1),cpoint = cp,r = 50)},import='auto')
  areap(crd_escc1,ap_escc1)
  id1 = ap_escc1$id[ap_escc1$type!='isolated']
  crd_escc1 = crd_escc1[id1,c('x','y')]
  as_escc1 = alphahull::ashape(crd_escc1,alpha = 100)
  plot(as_escc1)
  polyescc1 = extractpoly(as_escc1)
  plotpoly(polyescc1)
  polyescc1 = polyescc1[-c(4,5,6)]
  
  crd_escc2 = crd[tmp$cellID[tmp$stage=='ESCC-2'],]
  crd_escc2 = crd_escc2[,c('x','y')]
  as_escc2 = alphahull::ashape(crd_escc2,alpha = 300)
  plot(as_escc2)
  polyescc2 = extractpoly(as_escc2)
  plotpoly(polyescc2)

  crd_escc3 = crd[tmp$cellID[tmp$stage=='ESCC-3'],]
  crd_escc3 = crd_escc3[,c('x','y')]
  as_escc3 = alphahull::ashape(crd_escc3,alpha = 300)
  plot(as_escc3)
  polyescc3 = extractpoly(as_escc3)
  plotpoly(polyescc3)
  
  crd_hln = crd[tmp$cellID[tmp$stage1!='ESCC'],]
  plotly::plot_ly(crd_hln,x=~x,y=~y)
  cp = c('T4S2_bcapfhbj-1')
  job({ap_hln = areaPoint(dist(crd_hln),cpoint = cp,r = 50)},import='auto')
  areap(crd_hln,ap_hln)
  id1 = ap_hln$id[ap_hln$type!='isolated']
  crd_hln = crd_hln[id1,c('x','y')]
  as_hln = alphahull::ashape(crd_hln,alpha = 100)
  plot(as_hln)
  polyhln = extractpoly(as_hln)
  plotpoly(polyhln)
  
  poly_42 = Reduce(append,list(polyescc1,polyescc2,polyescc3,polyhln))
  plotpoly(poly_42)
  saveRDS(poly_42,'stp3_tme/data/epiregion/T4S2_poly.rds')
}
{# t5s1
  tmp = subset(xe5_1, cells = xe5_1$cellID[xe5_1$ct3=='Epi'])
  xedimplot(tmp,groupby = 'stage')
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)

  plotly::plot_ly(crd,x=~x,y=~y)
  cp = c('T5S1_aecemjfk-1')
  job({ap = areaPoint(dist(crd),cpoint = cp,r = 50)},import='auto')
  areap(crd,ap)
  id1 = ap$id[ap$type!='isolated']
  crd = crd[id1,c('x','y')]
  as = alphahull::ashape(crd,alpha = 100)
  plot(as)
  polytmp = extractpoly(as)
  plotpoly(polytmp)
  
  poly_51 = polytmp
  plotpoly(poly_51)
  saveRDS(poly_51,'stp3_tme/data/epiregion/T5S1_poly.rds')
}
{# t5s2
  tmp = subset(xe5_2_split, cells = xe5_2_split$cellID[xe5_2_split$ct3=='Epi'])
  xedimplot(tmp,groupby = 'stage')
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)

  plotly::plot_ly(crd,x=~x,y=~y)
  cp = c('T5S2_bgpiikon-1')
  job({ap = areaPoint(dist(crd),cpoint = cp,r = 50)},import='auto')
  areap(crd,ap)
  id1 = ap$id[ap$type!='isolated']
  crd = crd[id1,c('x','y')]
  as = alphahull::ashape(crd,alpha = 100)
  plot(as)
  polytmp = extractpoly(as)
  plotpoly(polytmp)
  
  poly_52 = polytmp
  plotpoly(poly_52)
  saveRDS(poly_52,'stp3_tme/data/epiregion/T5S2_poly.rds')
}
{# t5s3
  tmp = subset(xe5_3_split, cells = xe5_3_split$cellID[xe5_3_split$ct3=='Epi'])
  xedimplot(tmp,groupby = 'stage')
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)

  plotly::plot_ly(crd,x=~x,y=~y)
  cp = c('T5S3_belmiekd-1')
  job({ap = areaPoint(dist(crd),cpoint = cp,r = 50)},import='auto')
  areap(crd,ap)
  id1 = ap$id[ap$type!='isolated']
  crd = crd[id1,c('x','y')]
  as = alphahull::ashape(crd,alpha = 100)
  plot(as)
  polytmp = extractpoly(as)
  plotpoly(polytmp)
  
  poly_53 = polytmp
  plotpoly(poly_53)
  saveRDS(poly_53,'stp3_tme/data/epiregion/T5S3_poly.rds')
}
{# t6s1
  tmp = subset(xe6_1, cells = xe6_1$cellID[xe6_1$ct3=='Epi'])
  xedimplot(tmp,groupby = 'stage',reversey = T)
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)

  crd_escc = crd[tmp$cellID[tmp$stage1=='ESCC'],]
  plotly::plot_ly(crd_escc,x=~x,y=~y)
  cp = c('T6S1_hdolbmdo-1')
  job({ap_escc = areaPoint(dist(crd_escc),cpoint = cp,r = 50)},import='auto')
  areap(crd_escc,ap_escc)
  id1 = ap_escc$id[ap_escc$type!='isolated']
  crd_escc = crd_escc[id1,c('x','y')]
  as_escc = alphahull::ashape(crd_escc,alpha = 100)
  plot(as_escc)
  polyescc = extractpoly(as_escc)
  plotpoly(polyescc)

  crd_hn = crd[tmp$cellID[tmp$stage1!='ESCC'],]
  plotly::plot_ly(crd_hn,x=~x,y=~y)
  cp = c('T6S1_npbjhlcc-1')
  job({ap_hn = areaPoint(dist(crd_hn),cpoint = cp,r = 50)},import='auto')
  areap(crd_hn,ap_hn)
  id1 = ap_hn$id[ap_hn$type!='isolated']
  crd_hn = crd_hn[id1,c('x','y')]
  as_hn = alphahull::ashape(crd_hn,alpha = 100)
  plot(as_hn)
  polyhn = extractpoly(as_hn)
  plotpoly(polyhn)
  
  poly_61 = Reduce(append,list(polyescc,polyhn))
  plotpoly(poly_61)
  saveRDS(poly_61,'stp3_tme/data/epiregion/T6S1_poly.rds')
}
{# t6s2
  tmp = subset(xe6_2, cells = xe6_2$cellID[xe6_2$ct3=='Epi'])
  xedimplot(tmp,groupby = 'stage')
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)
  
  crd_escc = crd[tmp$cellID[tmp$stage%in%c('ESCC-1','ESCC-2')],]
  plotly::plot_ly(crd_escc,x=~x,y=~y)
  cp = c('T6S2_ckoclijc-1')
  job({ap_escc = areaPoint(dist(crd_escc),cpoint = cp,r = 50)},import='auto')
  areap(crd_escc,ap_escc)
  id1 = ap_escc$id[ap_escc$type!='isolated']
  crd_escc = crd_escc[id1,c('x','y')]
  as_escc = alphahull::ashape(crd_escc,alpha = 100)
  plot(as_escc)
  polyescc = extractpoly(as_escc)
  plotpoly(polyescc)
  polyescc = polyescc[-2]

  crd_hn = crd[tmp$cellID[!tmp$stage%in%c('ESCC-1','ESCC-2')],]
  plotly::plot_ly(crd_hn,x=~x,y=~y)
  cp = c('T6S2_ednopeek-1')
  job({ap_hn = areaPoint(dist(crd_hn),cpoint = cp,r = 50)},import='auto')
  areap(crd_hn,ap_hn)
  id1 = ap_hn$id[ap_hn$type!='isolated']
  crd_hn = crd_hn[id1,c('x','y')]
  as_hn = alphahull::ashape(crd_hn,alpha = 100)
  plot(as_hn)
  polyhn = extractpoly(as_hn)
  plotpoly(polyhn)
  
  poly_62 = Reduce(append,list(polyescc,polyhn))
  plotpoly(poly_62)
  saveRDS(poly_62,'stp3_tme/data/epiregion/T6S2_poly.rds')
}
{# t6s3
  tmp = subset(xe6_3, cells = xe6_3$cellID[xe6_3$ct3=='Epi'])
  xedimplot(tmp,groupby = 'stage')
  crd = tmp@images$fov@boundaries$centroids@coords %>% data.frame
  rownames(crd) = colnames(tmp)
  
  crd_escc = crd[tmp$cellID[tmp$stage%in%c('ESCC-1','ESCC-2')],]
  plotly::plot_ly(crd_escc,x=~x,y=~y)
  cp = c('T6S3_klmannch-1','T6S3_ckgnepco-1')
  job({ap_escc = areaPoint(dist(crd_escc),cpoint = cp,r = 50)},import='auto')
  areap(crd_escc,ap_escc)
  id1 = ap_escc$id[ap_escc$type!='isolated']
  crd_escc = crd_escc[id1,c('x','y')]
  as_escc = alphahull::ashape(crd_escc,alpha = 60)
  plot(as_escc)
  polyescc = extractpoly(as_escc)
  plotpoly(polyescc)
  polyescc = polyescc[-3]

  crd_hn = crd[tmp$cellID[!tmp$stage%in%c('ESCC-1','ESCC-2')],]
  plotly::plot_ly(crd_hn,x=~x,y=~y)
  cp = c('T6S3_licblgcj-1')
  job({ap_hn = areaPoint(dist(crd_hn),cpoint = cp,r = 50)},import='auto')
  areap(crd_hn,ap_hn)
  id1 = ap_hn$id[ap_hn$type!='isolated']
  crd_hn = crd_hn[id1,c('x','y')]
  as_hn = alphahull::ashape(crd_hn,alpha = 100)
  plot(as_hn)
  polyhn = extractpoly(as_hn)
  plotpoly(polyhn)
  
  poly_63 = Reduce(append,list(polyescc,polyhn))
  plotpoly(poly_63)
  saveRDS(poly_63,'stp3_tme/data/epiregion/T6S3_poly.rds')
}











