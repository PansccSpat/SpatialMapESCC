source('~/xenium/header.R')
source('~/xenium/xenium_function.R')




# input #####
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

allborder = read_rds('stp3_tme/data/epiborder/allborder.rds')


primcol2 = c('#1f78b4','#ffc089','#B5EFB5','#793c1b',
             '#6a3d9a','#333333','#ffff33','#df928e')
names(primcol2) = c('B','Endo','Fib','Gland','Mye','Mus','T','Epi')





# gene expr spatial full #####
objlist = setNames(list('xe1_1','xe1_2','xe1_3','xe2_1','xe2_2','xe2_3','xe3_1','xe3_2',
                        'xe4_1','xe4_2','xe5_1','xe5_2_split','xe5_3_split','xe6_1','xe6_2','xe6_3'),
                   c('T1S1','T1S2','T1S3','T2S1','T2S2','T2S3','T3S1','T3S2',
                     'T4S1','T4S2','T5S1','T5S2','T5S3','T6S1','T6S2','T6S3'))
job({
  lapply(names(objlist), function(x){
    print(x)
    dt = get(objlist[[x]])
    dt = subset(dt, cells = dt$cellID[dt$ct3!='Mus'])
    
    bd = allborder[[x]][['d']]
    
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
    
    exprcol = rev(brewer.pal(9,'RdBu'))
    exprcol = brewer.pal(9,'Reds')
    fname = paste0('stp3_tme/plottmp_expr/',x,'_expr.png')
    p1 = xefeatureplot(dt, feature = 'JAG1',flip = flip,reversex = reversex,reversey = reversey,color = exprcol,bgcol = 'white',pt.size = 0.1)+geom_path(data=bd,linetype='dashed')
    p2 = xefeatureplot(dt, feature = 'FAP',flip = flip,reversex = reversex,reversey = reversey,color = exprcol,bgcol = 'white',pt.size = 0.1)+geom_path(data=bd,linetype='dashed')
    p3 = xefeatureplot(dt, feature = 'ACTA2',flip = flip,reversex = reversex,reversey = reversey,color = exprcol,bgcol = 'white',pt.size = 0.1)+geom_path(data=bd,linetype='dashed')
    
    p = p1/p2/p3
    
    print('Saving')
    ggsave(fname,p,width = 30,height = 50,limitsize = F)
  })
},import = 'auto')





# gene expr select #####
exprcol = brewer.pal(9,'Reds')
themetrans=theme(panel.background=element_rect(fill='transparent'),plot.background=element_rect(fill='transparent',color = NA))

exprcol1 = c('white','white',brewer.pal(9,'Blues'))
exprcol2 = c('white','white',brewer.pal(9,'Reds'))
exprcol3 = c(rep('white',4),brewer.pal(9,'Greens'))


roin = c(5500,7100,1500,2300)
xedimplot(xe2_1,groupby = 'ct3',color = primcol2,flip = T,reversex = T,axis = T,gridline = T,bgcol = 'white',zoomin = roin,pt.size = 1)
dtn = setdiff(getcellinrange(xe2_1,roin,flip = T),xe2_1$cellID[xe2_1$ct3=='Mus'])
dtn = subset(xe2_1,cells = dtn)
p1 = xedimplot(dtn,zoomin = roin,groupby = 'ct3',color = primcol2,flip = T,reversex = T,gridline = T,bgcol = 'white',pt.size = 1)+NoLegend()
p2 = xefeatureplot(dtn,zoomin = roin, feature = 'JAG1',flip = T,reversex = T,gridline = T,color = exprcol,bgcol = 'white',pt.size = 1)
p3 = xefeatureplot(dtn,zoomin = roin, feature = 'FAP',flip = T,reversex = T,gridline = T,color = exprcol,bgcol = 'white',pt.size = 1)
p4 = xefeatureplot(dtn,zoomin = roin, feature = 'ACTA2',flip = T,reversex = T,gridline = T,color = exprcol,bgcol = 'white',pt.size = 1)
(p1+p2)/(p3+p4)


# roil = c(9000,10600,2700,3500)
# xedimplot(xe3_1,groupby = 'ct3',color = primcol2,flip = T,reversex = T,axis = T,gridline = T,bgcol = 'white',zoomin = roil,pt.size = 1)
# dtl = setdiff(getcellinrange(xe3_1,roil,flip = T),xe3_1$cellID[xe3_1$ct3=='Mus'])
# dtl = subset(xe3_1,cells = dtl)
roil = c(19000,20600,3200,4000)
xedimplot(xe1_2,groupby = 'ct3',color = primcol2,flip = T,reversex = T,axis = T,gridline = T,bgcol = 'white',zoomin = roil,pt.size = 1)
dtl = setdiff(getcellinrange(xe1_2,roil,flip = T),xe1_2$cellID[xe1_2$ct3=='Mus'])
dtl = subset(xe1_2,cells = dtl)
p1 = xedimplot(dtl,zoomin = roil,groupby = 'ct3',color = primcol2,flip = T,reversex = T,gridline = T,bgcol = 'white',pt.size = 1)+NoLegend()
p2 = xefeatureplot(dtl,zoomin = roil, feature = 'JAG1',flip = T,reversex = T,gridline = T,color = exprcol,bgcol = 'white',pt.size = 1)
p3 = xefeatureplot(dtl,zoomin = roil, feature = 'FAP',flip = T,reversex = T,gridline = T,color = exprcol,bgcol = 'white',pt.size = 1)
p4 = xefeatureplot(dtl,zoomin = roil, feature = 'ACTA2',flip = T,reversex = T,gridline = T,color = exprcol,bgcol = 'white',pt.size = 1)
(p1+p2)/(p3+p4)


roih = c(8800,11800,2100,3600)
xedimplot(xe1_3,groupby = 'ct3',color = primcol2,flip = T,reversex = T,axis = T,gridline = T,bgcol = 'white',zoomin = roih,pt.size = 1)
dth = setdiff(getcellinrange(xe1_3,roih,flip = T),xe1_3$cellID[xe1_3$ct3=='Mus'])
dth = subset(xe1_3,cells = dth)
p1 = xedimplot(dth,zoomin = roih,groupby = 'ct3',color = primcol2,flip = T,reversex = T,gridline = T,bgcol = 'white',pt.size = 1)+NoLegend()
p2 = xefeatureplot(dth,zoomin = roih, feature = 'JAG1',flip = T,reversex = T,gridline = T,color = exprcol,bgcol = 'white',pt.size = 1)
p3 = xefeatureplot(dth,zoomin = roih, feature = 'FAP',flip = T,reversex = T,gridline = T,color = exprcol,bgcol = 'white',pt.size = 1)
p4 = xefeatureplot(dth,zoomin = roih, feature = 'ACTA2',flip = T,reversex = T,gridline = T,color = exprcol,bgcol = 'white',pt.size = 1)
(p1+p2)/(p3+p4)


roit = c(8600,11600,700,2200)
xedimplot(xe3_1,groupby = 'ct3',color = primcol2,flip = T,reversex = T,axis = T,gridline = T,bgcol = 'white',zoomin = roit,pt.size = 1)
dtt = setdiff(getcellinrange(xe3_1,roit,flip = T),xe3_1$cellID[xe3_1$ct3=='Mus'])
dtt = subset(xe3_1,cells = dtt)
p1 = xedimplot(dtt,zoomin = roit,groupby = 'ct3',color = primcol2,flip = T,reversex = T,gridline = T,bgcol = 'white',pt.size = 1)+NoLegend()
p2 = xefeatureplot(dtt,zoomin = roit, feature = 'JAG1',flip = T,reversex = T,gridline = T,color = exprcol,bgcol = 'white',pt.size = 1)
p3 = xefeatureplot(dtt,zoomin = roit, feature = 'FAP',flip = T,reversex = T,gridline = T,color = exprcol,bgcol = 'white',pt.size = 1)
p4 = xefeatureplot(dtt,zoomin = roit, feature = 'ACTA2',flip = T,reversex = T,gridline = T,color = exprcol,bgcol = 'white',pt.size = 1)
(p1+p2)/(p3+p4)



p1n = xedimplot(dtn,zoomin = roin,groupby = 'ct3',color = primcol2,flip = T,reversex = T,gridline = F,bgcol = 'transparent',pt.size = 1)+themetrans+NoLegend()
p1l = xedimplot(dtl,zoomin = roil,groupby = 'ct3',color = primcol2,flip = T,reversex = T,gridline = F,bgcol = 'transparent',pt.size = 1)+themetrans+NoLegend()
p1h = xedimplot(dth,zoomin = roih,groupby = 'ct3',color = primcol2,flip = T,reversex = T,gridline = F,bgcol = 'transparent',pt.size = 1)+themetrans+NoLegend()
p1t = xedimplot(dtt,zoomin = roit,groupby = 'ct3',color = primcol2,flip = T,reversex = T,gridline = F,bgcol = 'transparent',pt.size = 1)+themetrans+NoLegend()
p1 = p1n/p1l/p1h/p1t
ggsave('stp3_tme/plot/expr_spatial/expr_fovct.png',p1,width = 15,height = 20,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_fovct_n.png',p1n,width = 15,height = 5,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_fovct_l.png',p1l,width = 15,height = 5,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_fovct_h.png',p1h,width = 15,height = 5,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_fovct_t.png',p1t,width = 15,height = 5,bg = 'transparent')

p2n = xefeatureplot(dtn,zoomin = roin, feature = 'JAG1',flip = T,reversex = T,gridline = F,color = exprcol1,bgcol = 'white',pt.size = 1,collim = c(0,3))+themetrans
p2l = xefeatureplot(dtl,zoomin = roil, feature = 'JAG1',flip = T,reversex = T,gridline = F,color = exprcol1,bgcol = 'white',pt.size = 1,collim = c(0,3))+themetrans
p2h = xefeatureplot(dth,zoomin = roih, feature = 'JAG1',flip = T,reversex = T,gridline = F,color = exprcol1,bgcol = 'white',pt.size = 1,collim = c(0,3))+themetrans
p2t = xefeatureplot(dtt,zoomin = roit, feature = 'JAG1',flip = T,reversex = T,gridline = F,color = exprcol1,bgcol = 'white',pt.size = 1,collim = c(0,3))+themetrans
p2 = p2n/p2l/p2h/p2t
# ggsave('stp3_tme/plot/expr_spatial/expr_jag1.pdf',p2,width = 15,height = 20)
ggsave('stp3_tme/plot/expr_spatial/expr_jag1.pdf',p2n,width = 15,height = 5,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_jag1_n.png',p2n,width = 15,height = 5,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_jag1_l.png',p2l,width = 15,height = 5,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_jag1_h.png',p2h,width = 15,height = 5,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_jag1_t.png',p2t,width = 15,height = 5,bg = 'transparent')

p3n = xefeatureplot(dtn,zoomin = roin, feature = 'FAP',flip = T,reversex = T,gridline = F,color = exprcol2,bgcol = 'white',pt.size = 1,collim = c(0,2.5))+themetrans
p3l = xefeatureplot(dtl,zoomin = roil, feature = 'FAP',flip = T,reversex = T,gridline = F,color = exprcol2,bgcol = 'white',pt.size = 1,collim = c(0,2.5))+themetrans
p3h = xefeatureplot(dth,zoomin = roih, feature = 'FAP',flip = T,reversex = T,gridline = F,color = exprcol2,bgcol = 'white',pt.size = 1,collim = c(0,2.5))+themetrans
p3t = xefeatureplot(dtt,zoomin = roit, feature = 'FAP',flip = T,reversex = T,gridline = F,color = exprcol2,bgcol = 'white',pt.size = 1,collim = c(0,2.5))+themetrans
p3 = p3n/p3l/p3h/p3t
# ggsave('stp3_tme/plot/expr_spatial/expr_fap.pdf',p3,width = 15,height = 20)
ggsave('stp3_tme/plot/expr_spatial/expr_fap.pdf',p3n,width = 15,height = 5,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_fap_n.png',p3n,width = 15,height = 5,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_fap_l.png',p3l,width = 15,height = 5,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_fap_h.png',p3h,width = 15,height = 5,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_fap_t.png',p3t,width = 15,height = 5,bg = 'transparent')

p4n = xefeatureplot(dtn,zoomin = roin, feature = 'ACTA2',flip = T,reversex = T,gridline = F,color = exprcol3,bgcol = 'white',pt.size = 1,collim = c(0,3.5))+themetrans
p4l = xefeatureplot(dtl,zoomin = roil, feature = 'ACTA2',flip = T,reversex = T,gridline = F,color = exprcol3,bgcol = 'white',pt.size = 1,collim = c(0,3.5))+themetrans
p4h = xefeatureplot(dth,zoomin = roih, feature = 'ACTA2',flip = T,reversex = T,gridline = F,color = exprcol3,bgcol = 'white',pt.size = 1,collim = c(0,3.5))+themetrans
p4t = xefeatureplot(dtt,zoomin = roit, feature = 'ACTA2',flip = T,reversex = T,gridline = F,color = exprcol3,bgcol = 'white',pt.size = 1,collim = c(0,3.5))+themetrans
p4 = p4n/p4l/p4h/p4t
# ggsave('stp3_tme/plot/expr_spatial/expr_acta2.pdf',p4,width = 15,height = 20)
ggsave('stp3_tme/plot/expr_spatial/expr_acta2.pdf',p4n,width = 15,height = 5,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_acta2_n.png',p4n,width = 15,height = 5,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_acta2_l.png',p4l,width = 15,height = 5,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_acta2_h.png',p4h,width = 15,height = 5,bg = 'transparent')
ggsave('stp3_tme/plot/expr_spatial/expr_acta2_t.png',p4t,width = 15,height = 5,bg = 'transparent')









