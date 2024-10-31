source('~/xenium/header.R')
source('~/xenium/xenium_function.R')



# input #####
km_test = read_rds('stp3_tme/data/Cellular_Neighborhood/cellfrac_50/kmeans_15.rds')

allmeta_m = read_rds('stp5_batch2process/data/allcell/allmerge_allmeta.rds')

epidt = read_rds('stp2_epi/data/epiid_correctstage.rds')
crds = read_rds('stp1_summary/data/allsample_coords.rds')
epipoly = read_rds('stp3_tme/data/epiregion/allsample_poly.rds')
allborder = read_rds('stp3_tme/data/epiborder/allborder.rds')

gridid = read_rds('stp1_summary/data/allsample_gridid.rds')

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

coltmp = c('#80b1d3','#8dd3c7','#ccebc5','#fccde5','red','#fdb462','#bc80bd','#b3de69','#ffffb3','#bebada','#d9d9d9')
coltmp = c('#fdb462','#8dd3c7','#ccebc5','#fccde5','red','#80b1d3','#bc80bd','#ffff33','#bebada','#a6761d','#d9d9d9')
coltmp = c('#ccebc5','#c6dbef','#238b45','#de0c62','#ff9408','#666666','#793c1b','#BC9DCC','#8dd3c7','#d9d9d9','#bbbbbb')
names(coltmp) = c('Normal epi','Normal epi & Invasive epi','Normal epi & NF','Invasive epi','CAF-epi niche','CAF','TLS','NF','Endo & Fibro','Gland','Muscle')
coltmp = coltmp[c('CAF','TLS','NF','Endo & Fibro','Gland','Muscle','Normal epi','Normal epi & Invasive epi','Normal epi & NF','Invasive epi','CAF-epi niche')]


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


allborder = read_rds('stp3_tme/data/epiborder/allborder.rds')
t1dist = read_rds('stp3_tme/data/epi_borderdist/T1epi_borderdist_v2.rds')
t1border = read_rds('stp3_tme/data/epiborder/T1border_v2.rds')
t3dist = read_rds('stp3_tme/data/epi_borderdist/T3ESCC_borderdist_v2.rds')
t3border = read_rds('stp3_tme/data/epiborder/T3ESCCborder_v2.rds')
t4dist = read_rds('stp3_tme/data/epi_borderdist/T4partE_borderdist_v2.rds')
t4dtlist = read_rds('stp3_tme/data/epiborder/T4partEborder_v2.rds')




# func #####
plotcluster = function(meta, groupby, 
                       fixed=T, color=NULL, pt.size=0.2,pt.shape=16,bgcol = 'white',ptlevels=NULL,
                       reversex=F, reversey=F,flip=F,
                       highlight = NULL,xintercept=NULL,yintercept=NULL,
                       zoomin = NULL, cellborder = F, areaalpha = 0.5, cellgap = 0, cellradius = 0, linewidth = 0.5, celllinecol=NULL,
                       gridline=F,gridcol='black',gridgap=100,gridlinewidth=0.03,axis=F,plotborder=T){
  ### coords
  {
    coords = meta[,c('x','y',groupby)]
    rownames(coords) = rownames(meta)
    coords[,groupby] = coords[,groupby] %>% factor()
    
    if(flip){
      coords = coords %>% dplyr::rename(c('x'='y','y'='x'))
    }
    
    xmin = min(coords$x)
    ymin = min(coords$y)
    xmax = max(coords$x)
    ymax = max(coords$y)
    
    if(length(zoomin)==4){
      zoomin = zoomin %>% as.numeric()
      coords = coords %>% filter(x>=zoomin[1],x<=zoomin[2],y>=zoomin[3],y<=zoomin[4])
      xmin = zoomin[1]
      xmax = zoomin[2]
      ymin = zoomin[3]
      ymax = zoomin[4]
    }
    
    if(!is.null(highlight)){
      # coords[rownames(coords)%in%highlight,groupby] = 'group_1'
      levels(coords[,groupby]) = levels(coords[,groupby]) %>% append('Others')
      coords[!rownames(coords)%in%highlight,groupby] = 'Others'
      coords[,groupby] = coords[,groupby] %>% droplevels()
      
      getPalette = colorRampPalette(brewer.pal(8, "Set1"))
      
      if(is.null(color)){
        cname = unique(coords[,groupby]) %>% setdiff('Others')
        num = cname %>% length()
        color = c(getPalette(num),'#d3d3d3')
        names(color) = c(cname,'Others')
      }else{
        nms = c(names(color),'Others')
        color=c(color,'#d3d3d3')
        names(color) = nms
      }
      
      if(is.null(ptlevels)){
        ptlevels = levels(coords[,groupby])
      }else{
        ptlevels = unique(c(ptlevels,'Others'))
      }
    }
  }
  
  
  ### make plot
  {
    p_bg = 
      coords %>% 
      arrange(groupby) %>% 
      ggplot(aes(x,y))+
      geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill=NA,color=NA)+
      theme1+
      theme(panel.background = element_rect(fill=bgcol))+
      labs(color=groupby)
    
    if(gridline){
      glx = (xmax - xmin)%/%gridgap
      if(glx>0){
        glx1 = c(0,glx+1)*gridgap+xmin
        glx = (0:(glx+1))*gridgap+xmin
      }else{
        glx1 = c(0,glx+1)*gridgap+xmin
        glx=NULL
      }
      gly = (ymax - ymin)%/%gridgap
      if(gly>0){
        gly1 = c(0,gly+1)*gridgap+ymin
        gly = (0:(gly+1))*gridgap+ymin
      }else{
        gly1 = c(0,gly+1)*gridgap+ymin
        gly=NULL
      }
      
      p_bg = p_bg+
        xlab('')+ylab('')+
        scale_x_continuous(expand = c(0,0))+
        scale_y_continuous(expand = c(0,0))+
        geom_hline(yintercept = gly,linewidth=gridlinewidth,color=gridcol)+
        geom_vline(xintercept = glx,linewidth=gridlinewidth,color=gridcol)
    }
    
    
    if(!is.null(ptlevels)){
      p_layer = lapply(rev(ptlevels), function(x){
        dt = coords[coords[,groupby]==x,]
        rst = geom_point(data=dt, aes(x,y,color=.data[[groupby]]),shape=pt.shape,size=pt.size)
        return(rst)
      })
      p_layer = append(list(p_bg),p_layer)
      p = Reduce('+',p_layer)
    }else{p = p_bg+geom_point(aes(color=.data[[groupby]]),shape=pt.shape,size = pt.size)}
    
    
    if(cellborder){
      bd = xenium@images$fov@boundaries$segmentation@polygons
      bd = bd[rownames(coords)] %>% 
        lapply(function(x){
          dt = x@Polygons[[1]]@coords %>% data.frame()
          dt$ID = x@ID
          return(dt)
        }) %>% rbindlist() %>% data.frame()
      
      celltype = data.frame(rownames(coords), coords[,3])
      colnames(celltype) = c('ID',groupby)
      bd = bd %>% left_join(celltype)
      
      if(flip){
        bd = bd %>% dplyr::rename(c('x'='y','y'='x'))
      }
      
      if(!is.null(ptlevels)){
        if(is.null(celllinecol)){
          p_layer1 = lapply(rev(ptlevels), function(x){
            dt = bd[bd[,groupby]==x,]
            rst = ggforce::geom_shape(data=dt,
                                      aes(x,y,fill = .data[[groupby]],group=ID,color=.data[[groupby]]),
                                      linewidth = linewidth,
                                      alpha = areaalpha,
                                      expand = cellgap,
                                      radius = cellradius)
            return(rst)
          })
        }else{
          p_layer1 = lapply(rev(ptlevels), function(x){
            dt = bd[bd[,groupby]==x,]
            rst = ggforce::geom_shape(data=dt,
                                      aes(x,y,fill = .data[[groupby]],group=ID),
                                      color=celllinecol,
                                      linewidth = linewidth,
                                      alpha = areaalpha,
                                      expand = cellgap,
                                      radius = cellradius)
            return(rst)
          })
        }
        
        p_layer2 = list()
        p_layer2[1] = p_layer[1]
        p_layer2[1:length(ptlevels)*2] = p_layer[1:length(ptlevels)+1]
        p_layer2[1:length(ptlevels)*2+1] = p_layer1
        p = Reduce('+',p_layer2)
      }else{
        if(is.null(celllinecol)){
          p = p+
            ggforce::geom_shape(data=bd,
                                aes(x,y,fill = .data[[groupby]],group=ID,color=.data[[groupby]]),
                                linewidth = linewidth,
                                alpha = areaalpha,
                                expand = cellgap,
                                radius = cellradius)
        }else{
          p = p+
            ggforce::geom_shape(data=bd,
                                aes(x,y,fill = .data[[groupby]],group=ID),
                                color=celllinecol,
                                linewidth = linewidth,
                                alpha = areaalpha,
                                expand = cellgap,
                                radius = cellradius)
        }
      }
    }
    if(!is.null(color)){
      p = p+scale_color_manual(values = color)+scale_fill_manual(values = color)
    }
    if(reversex){
      p=p+scale_x_reverse()
    }
    if(reversey){
      p=p+scale_y_reverse()
    }
    if(fixed){
      p = p+coord_fixed()
    }
    if(!is.null(xintercept)){
      p = p+geom_vline(xintercept = xintercept, color='black', linetype='dashed')
    }
    if(!is.null(yintercept)){
      p = p+geom_hline(yintercept = yintercept, color='black', linetype='dashed')
    }
    if(!axis){
      p = p+theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank())
    }
    if(plotborder){
      p = p+theme(panel.border = element_rect(color='black',linewidth = 1,fill=NA))+xlab('')+ylab('')
    }
  }
  
  p = p+guides(color=guide_legend(override.aes = list(size=3)))
  return(p)
}

plotkmeans = function(frac,meta, returnheat=T){
  plotcol = c('#440154','#984ea3','#3288bd','#35b779','#fde725','red')
  ctlevel = c('Non-invasive','Invasive','CAF','T','B','Mye','NF','Endo','Gland','Mus')
  
  dt = frac %>% gather('ct','frac',-coreid) %>% left_join(meta)
  dt_aggr = aggregate(list(frac=dt$frac),by=as.list(dt[,c('ct','cluster')]),mean)
  
  if(returnheat){
    # heat 
    dt_aggr1 = dt_aggr %>% spread(ct,frac)
    rownames(dt_aggr1)=dt_aggr1$cluster %>% as.character()
    dt_aggr1=dt_aggr1[,c('cluster',ctlevel)]
    p = ComplexHeatmap::pheatmap(dt_aggr1[,-1],color=plotcol,cluster_cols = F,cluster_rows = F)
  }else{
    dt_prop = aggregate(list(prop=(dt$frac>0.05)),by=as.list(dt[,c('ct','cluster')]),mean)

    dt_plot = left_join(dt_aggr,dt_prop)
    dt_plot$cluster = dt_plot$cluster %>% factor(levels = levels(meta$cluster))
    dt_plot$ct = dt_plot$ct %>% factor(levels = ctlevel)
    p = dt_plot %>%
      ggplot(aes(ct,cluster))+
      theme1+rotate_x_text()+
      scale_fill_gradientn(colors = plotcol)+
      geom_point(aes(fill=frac,size=prop),shape=21)
  }
  return(p)
}

cellfrac_all = lapply(grep('T',unique(allmeta_m$sid),value=T), function(slideid){
  rst = read_rds(paste0('stp3_tme/data/Cellular_Neighborhood/cellfrac_50/',slideid,'_cellfrac.rds'))
}) %>% Reduce(rbind,.)
rownames(cellfrac_all) = cellfrac_all$coreid

km_test = read_rds('stp3_tme/data/Cellular_Neighborhood/cellfrac_50/kmeans_15.rds')
allmeta_m = read_rds('stp5_batch2process/data/allcell/allmerge_allmeta.rds')

clevels = c('Normal epi','Normal epi & Invasive epi','Normal epi & NF','Invasive epi','CAF-epi niche','CAF','TLS','NF','Endo & Fibro','Gland','Muscle')
clevels1 = c('n','n-inv','n-nf','inv','ce','caf','tls','nf','en-f','gld','mus')


# summary #####
cnrst = read_rds('stp3_tme/data/Cellular_Neighborhood/cellfrac_50/clusterrst_all.rds')
propplot2(cnrst,coltmp,'cluster')+scale_y_continuous(expand = c(0,0))
ggsave('stp3_tme/plot/CN_primtype/cluster_prop_bar.pdf',width = 6,height = 6)

clusprop = getprop(cnrst$cluster,cnrst[,c('sampleID','stage')]) %>% gather('cluster','prop',all_of(names(coltmp)))
colnames(clusprop) = c('sampleID','stage','cluster','prop')

clusprop %>% 
  ggplot(aes(stage,prop))+
  theme1+NoLegend()+
  xlab('')+ylab('Cluster proportion')+
  scale_fill_manual(values = coltmp)+scale_color_manual(values = coltmp)+
  scale_y_continuous(expand = expansion(c(0,0.1)))+
  # stat_compare_means(comparisons = makecomparison(4))+
  geom_pwc(method = 'dunn_test')+
  ggbeeswarm::geom_quasirandom(aes(color=cluster))+
  geom_boxplot(aes(fill=cluster),outlier.color = NA,alpha=0.5)+
  facet_wrap(~cluster,scales = 'free')
ggsave('stp3_tme/plot/CN_primtype/cluster_prop_box.pdf',width = 12,height = 11)

{# test
  dunnrst = lapply(unique(clusprop$cluster),function(ct){
    print(as.character(ct))
    dt1= clusprop %>% filter(cluster==ct)
    dt1$stage = dt1$stage %>% droplevels()
    if(length(unique(dt1$stage))<3){
      stg1 = unique(dt1$stage)[1]
      stg2 = unique(dt1$stage)[2]
      rst = wilcox.test(dt1$prop[dt1$stage==stg1],dt1$prop[dt1$stage==stg2])
      rst = data.frame(Comparison=paste0(stg1,' - ',stg2),P.unadj=rst$p.value,celltype=ct)
    }else{
      rst = FSA::dunnTest(prop~stage,data=dt1)$res
      rst$celltype=ct
      rst = rst[,c(1,3,5)]
    }
    return(rst)
  }) %>% Reduce(rbind,.)
}





# cluster heatmap #####
cellfrac_all = lapply(grep('T',unique(allmeta_m$sid),value=T), function(slideid){
  rst = read_rds(paste0('stp3_tme/data/Cellular_Neighborhood/cellfrac_50/',slideid,'_cellfrac.rds'))
}) %>% Reduce(rbind,.)
rownames(cellfrac_all) = cellfrac_all$coreid
cnrst = read_rds('stp3_tme/data/Cellular_Neighborhood/cellfrac_50/clusterrst_all.rds')

pdf('stp3_tme/plot/CN_primtype/cluster_heat.pdf',width = 5,height = 7.5)
plotkmeans(cellfrac_all,cnrst)
dev.off()




# spat plot #####
cnrst = read_rds('stp3_tme/data/Cellular_Neighborhood/cellfrac_50/clusterrst_all.rds')

## T1S2 
{
  dttmp = cnrst %>% filter(sid=='T1S2')
  plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,ptlevels = rev(names(coltmp)))
  ggsave('stp3_tme/plot/CN_primtype/T1S2_clustspat.png',width = 20,height = 10)
  plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,ptlevels = rev(names(coltmp)),
              highlight = dttmp$coreid[dttmp$cluster%in%c('CAF','CAF-epi niche','Invasive epi')])
  ggsave('stp3_tme/plot/CN_primtype/T1S2_clustspat_highlight.png',width = 20,height = 10)
  
  p = plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,highlight = dttmp$coreid[dttmp$cluster%in%c('CAF','CAF-epi niche','Invasive epi')])/
    plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,highlight = dttmp$coreid[dttmp$cluster=='Invasive epi'])/
    plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,highlight = dttmp$coreid[dttmp$cluster=='CAF'])/
    plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,highlight = dttmp$coreid[dttmp$cluster=='CAF-epi niche'])
  ggsave('stp3_tme/plot/CN_primtype/T1S2_clustspat_split.png',p,width = 20,height = 30,limitsize=F)
}

## T1S3 
{
  dttmp = cnrst %>% filter(sid=='T1S3')
  plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,ptlevels = rev(names(coltmp)))
  ggsave('stp3_tme/plot/CN_primtype/T1S3_clustspat.png',width = 20,height = 10)
  plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,ptlevels = rev(names(coltmp)),
              highlight = dttmp$coreid[dttmp$cluster%in%c('CAF','CAF-epi niche','Invasive epi')])
  ggsave('stp3_tme/plot/CN_primtype/T1S3_clustspat_highlight.png',width = 20,height = 10)
  
  p = plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,highlight = dttmp$coreid[dttmp$cluster%in%c('CAF','CAF-epi niche','Invasive epi')])/
    plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,highlight = dttmp$coreid[dttmp$cluster=='Invasive epi'])/
    plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,highlight = dttmp$coreid[dttmp$cluster=='CAF'])/
    plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,highlight = dttmp$coreid[dttmp$cluster=='CAF-epi niche'])
  ggsave('stp3_tme/plot/CN_primtype/T1S3_clustspat_split.png',p,width = 20,height = 30,limitsize=F)
}

## T3S1 
{
  dttmp = cnrst %>% filter(sid=='T3S1')
  plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,ptlevels = rev(names(coltmp)))
  ggsave('stp3_tme/plot/CN_primtype/T3S1_clustspat.png',width = 20,height = 10)
  plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,ptlevels = rev(names(coltmp)),
              highlight = dttmp$coreid[dttmp$cluster%in%c('CAF','CAF-epi niche','Invasive epi')])
  ggsave('stp3_tme/plot/CN_primtype/T3S1_clustspat_highlight.png',width = 20,height = 10)
  
  p = plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,highlight = dttmp$coreid[dttmp$cluster%in%c('CAF','CAF-epi niche','Invasive epi')])/
    plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,highlight = dttmp$coreid[dttmp$cluster=='Invasive epi'])/
    plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,highlight = dttmp$coreid[dttmp$cluster=='CAF'])/
    plotcluster(dttmp,'cluster',reversex = T,color = coltmp,bgcol = 'white',gridline = T,highlight = dttmp$coreid[dttmp$cluster=='CAF-epi niche'])
  ggsave('stp3_tme/plot/CN_primtype/T3S1_clustspat_split.png',p,width = 20,height = 30,limitsize=F)
}




