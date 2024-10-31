source('~/xenium/header.R')
source('~/xenium/xenium_function.R')



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

allmeta = read_rds('data/xenium_seurat_processed/allcell_meta.rds')



primcol = brewer.pal(10, 'Paired')
primcol[7] = '#999999'
names(primcol) = c('B','Epi','Endo','Fib','Gland','Mye','Mus','Mast','Pla','T')

primcol2 = c('#1f78b4','#ffc089','#B5EFB5','#793c1b',
             '#6a3d9a','#333333','#ffff33','#df928e')
names(primcol2) = c('B','Endo','Fib','Gland','Mye','Mus','T','Epi')

group_colors <- c( 
  "Basal"="#0339f8",
  "Proliferative"="#FFA500",
  "Differentiated"="#76bbfdb3", # "#76bbfd"
  "Invasive"="#de0c62",
  
  "CD4-Treg" = "#E3BE00", 
  "CD4-Tfh"='#C5DEBA',
  "CD4-Tn/m"="#E59CC4", 
  
  "CD8-Tex"  = "#00441B",
  "CD8-Teff" = "#1B9E77",  
  "CD8-Tn/m" = '#984EA3', 
  
  "DC" = '#d78f4b',
  "Mast" = '#C6B5A0',
  "Mono" = '#E59CC4',
  "Macro-other" ='#b581b2' ,
  "Macro-LYVE1"='#743117',   
  "Macro-SPP1" ='#1e66ac',  
  
  'NF'= "#BC9DCC", 
  "CAF" ="#666666", #'#b91c23' 
  "Epithelial cells"="#ebc2bb",
  "other" = "lightgrey","Others" = "#E6E6E6"
)

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



# aggregate subtype #####
ctnew = lapply(list(fibmeta,tcmeta,myemeta),function(x){
  return(x[,c('cellID','cellType_merge')])
}) %>% Reduce(rbind,.)

ctnew$cellType_merge = ctnew$cellType_merge %>% 
  mapvalues(c('NF-PI16','CD8T-eff','CD8T-n/m','CD8T-ex','CD4T-n/m','CD4T-reg','CD4T-fh','Mac','Mac-2','Mac-3'),
            c('NF','Teff','CD8-naive','Tex','CD4-naive','Treg','Tfh','Macro-LYVE1','Macro-SPP1','Macro-IL1B')) %>% as.character()
allmeta1 = allmeta %>% mutate(cellType_merge=ct3)
allmeta1[ctnew$cellID,'cellType_merge'] = ctnew$cellType_merge
allmeta1$celltype = allmeta1$cellType_merge
allmeta1$cellType_merge = ifelse(allmeta1$cellType_merge%in%names(immcol),allmeta1$cellType_merge,'Others')
allmeta1$cellType_merge[allmeta1$cellType_merge=='CAF' & allmeta1$stage1%in%c('NOR','LGIN')] = 'Others'






# immune spatial #####
objlist = setNames(list('xe1_1','xe1_2','xe1_3','xe2_1','xe2_2','xe2_3','xe3_1','xe3_2',
                        'xe4_1','xe4_2','xe5_1','xe5_2_split','xe5_3_split','xe6_1','xe6_2','xe6_3'),
                   c('T1S1','T1S2','T1S3','T2S1','T2S2','T2S3','T3S1','T3S2',
                     'T4S1','T4S2','T5S1','T5S2','T5S3','T6S1','T6S2','T6S3'))
pfull = lapply(names(objlist), function(x){
  
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
  
  # psplit = lapply(unique(dt$sp_stg),function(y){
  #   dt1 = subset(dt,sp_stg==y)
  #   p = xedimplot(dt1, groupby = 'cellType_merge',flip = flip,reversex = reversex,reversey = reversey,color = immcol,bgcol = 'black',pt.size = 1)+
  #     ggtitle(y)+theme(plot.title = element_text(size=18))
  #   return(p)
  # }) %>% `names<-`(unique(dt$sp_stg))
  # 
  fname=paste0('stp3_tme/plot/immuneSpatial/full/',x,'_immspaital.png')
  if(!file.exists(fname)){
    print(x)
    dt = get(objlist[[x]])
    dt = renewMeta(dt, allmeta1[,c('cellID','cellType_merge')])
    
    pfull = xedimplot(dt, groupby = 'cellType_merge',flip = flip,reversex = reversex,reversey = reversey,
                      color = immcol,bgcol = 'black',pt.size = 0.1,ptlevels = names(immcol))
    # return(pfull)
    
    rg = getcoordsrange(dt,dt$cellID,flip = flip)
    wd=(rg[2]-rg[1])/500+5
    ht=(rg[4]-rg[3])/500
    ggsave(fname,pfull,width = wd,height = ht,limitsize = F)
  }
}) %>% `names<-`(names(objlist))

psplit = lapply(names(objlist), function(x){
  print(x)
  dt = get(objlist[[x]])
  dt = renewMeta(dt, allmeta1[,c('cellID','cellType_merge')])
  
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
  
  psplit = lapply(unique(dt$sp_stg),function(y){
    dt1 = subset(dt,sp_stg==y)
    p = xedimplot(dt1, groupby = 'cellType_merge',flip = flip,reversex = reversex,reversey = reversey,
                  color = immcol,bgcol = 'black',pt.size = 0.1,ptlevels = names(immcol))+
      ggtitle(y)+theme(plot.title = element_text(size=18))
    return(p)
  }) %>% `names<-`(unique(dt$sp_stg))

  return(psplit)
}) %>% `names<-`(names(objlist))





# imm spatial 2 #####
hltcell1 = c('NF','CAF')
hltcell2 = c('Macro-SPP1','CD8-Tn/m','CD8-Teff','CD8-Tex','CD4-Tn/m','CD4-Tfh','CD4-Treg')
hltcell = c(hltcell1,hltcell2)

## roie t1s2 #####
{
  ## add niche
  cnrst = read_rds('stp3_tme/data/Cellular_Neighborhood/cellfrac_50/clusterrst_all.rds')
  dttmp_t1s2 = cnrst %>% filter(sid=='T1S2') %>% filter(y>=1500)
  plotcluster(dttmp_t1s2,groupby = 'celltype',axis = T,reversex = T,gridline = T,color = group_colors,bgcol = 'white',pt.size = 0.5)
  bgalp = 1
  xmin = range(dttmp_t1s2$x)[1]
  xmax = range(dttmp_t1s2$x)[2]
  ymin = range(dttmp_t1s2$y)[1]
  ymax = range(dttmp_t1s2$y)[2]
  gridgap = 100
  
  p=plotcluster(dttmp_t1s2,groupby = 'celltype',reversex = T,gridline = T,color = group_colors,highlight = dttmp_t1s2$coreid[dttmp_t1s2$celltype%in%hltcell],
              bgcol = 'white',pt.size = 0.5)+
    geom_point(data=dttmp_t1s2[dttmp_t1s2$cluster=='CAF-epi niche',],aes(x,y),color='#ff9408',size=0.3)+
    geom_point(data=dttmp_t1s2[dttmp_t1s2$cluster=='Normal epi & NF',],aes(x,y),color='#238b45',size=0.3)
  
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
  
  p=dttmp_t1s2 %>%
    ggplot(aes(x,y))+
    theme1+theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank())+
    theme(panel.background = element_rect(color='black',fill=NA))+
    xlab('')+ylab('')+coord_fixed()+scale_x_reverse()+
    scale_color_manual(values = group_colors)+
    geom_hline(yintercept = gly,linewidth=0.03,color='black')+
    geom_vline(xintercept = glx,linewidth=0.03,color='black')+
    geom_point(data=dttmp_t1s2_3,color='#d3d3d3',size=0.2,alpha=bgalp,shape=16)+
    geom_point(data=dttmp_t1s2_1,aes(color=celltype),size=0.2,alpha=bgalp,shape=16)+
    geom_point(data=dttmp_t1s2[dttmp_t1s2$cluster=='CAF-epi niche',],color='#ff9408',size=0.5,alpha=bgalp,shape=16)+
    geom_point(data=dttmp_t1s2[dttmp_t1s2$cluster=='Normal epi & NF',],color='#238b45',size=0.5,alpha=bgalp,shape=16)+
    geom_point(data=dttmp_t1s2_2,aes(color=celltype),size=1,shape=16)
  ggsave('stp3_tme/plot/immuneSpatial/withNiche/T1S2.png',p,width = 30,height = 15)
}

## roie t3s1 #####
{
  roie = c(6000,14000,0,3000)
  xedimplot(xe3_1,'ct3',flip = T,reversex = T,axis=T,gridline = T,xintercept = roie[1:2],yintercept = roie[3:4],zoomin = roie,color = primcol2)
  
  ## add niche
  cnrst = read_rds('stp3_tme/data/Cellular_Neighborhood/cellfrac_50/clusterrst_all.rds')
  dttmp_t3s1 = cnrst %>% filter(sid=='T3S1') %>%
    filter(x<=roie[2],x>=roie[1],y<=roie[4],y>=roie[3])
  plotcluster(dttmp_t3s1,groupby = 'celltype',gridline = T,color = group_colors,highlight = dttmp_t3s1$coreid[dttmp_t3s1$celltype%in%hltcell],
              bgcol = 'white',pt.size = 0.5)+
    geom_point(data=dttmp_t3s1[dttmp_t3s1$cluster=='CAF-epi niche',],aes(x,y),color='#ff9408',size=0.5)+
    geom_point(data=dttmp_t3s1[dttmp_t3s1$cluster=='Normal epi & NF',],aes(x,y),color='#238b45',size=0.5)
  ggsave('stp3_tme/plot//ECMscore_ESCC_T3S1.pdf',width = 15,height = 10)
}




# epi area #####
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











