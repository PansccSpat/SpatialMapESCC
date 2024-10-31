source('~/xenium/header.R')
source('~/xenium/xenium_function.R')
library(ComplexHeatmap)
library(trend)



# input #####
allmeta_m = read_rds('stp5_batch2process/data/allcell/allmerge_allmeta.rds')
epimt = read_rds('stp5_batch2process/data/epi_merge/epi_merge_meta.rds')
crds = read_rds('stp1_summary/data/allsample_coords.rds')


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




epicol = c('#0339f8','#ff9408','#75bbfd','#de0c62')
names(epicol) = c('Basal','Proliferation','Differentiation','Invasive')
epicol1 = epicol 
names(epicol1) = names(epicol1) %>% mapvalues(c('Proliferation','Differentiation'),c('Proliferative','Differentiated'))

primcol2 = c('#1f78b4','#ffc089','#B5EFB5','#793c1b',
             '#6a3d9a','#333333','#ffff33','#df928e')
names(primcol2) = c('B','Endo','Fib','Gland','Mye','Mus','T','Epi')

progcol = c('Cell adhesion'="#0339f8",
            'Cycling'="#ff9408",
            'Mucosal defense'="#75bbfd",
            'Inflammation'="#B3DE69",
            'DNA repair'="#a6eee6",
            'Oncogenic TFs'="#9467bd",
            'Angiogenesis'="#41ab5d",
            'EMT'="#ff0000")





# preprocessing #####
borderdist = read_rds('stp3_tme/data/epi_borderdist/allepi_borderdist.rds')
colnames(borderdist) = colnames(borderdist) %>% mapvalues('downDist','borderdist')





# plot by fov ####
## t1 lgin fov #####
{# t1 lgin fov
  dir.create('stp2_epi/plot/epiProgram_spatial2/T1_L/',recursive = T)
  dt13 = allmeta_m %>% filter(sid=='T1S3') %>% filter(primaryType=='Epi')
  roi = c(19800,20800,2650,3250)
  plotcluster(dt13,'celltype',color = epicol1,xintercept = roi[1:2],yintercept = roi[3:4],zoomin = roi,axis = T,pt.size = 1)
  dt13_plot = dt13 %>% filter(x<=roi[2],x>=roi[1],y<=roi[4],y>=roi[3])
  
  inepic = borderdist$cellID[borderdist$inEpi]
  dt13_plot = dt13_plot %>% filter(cellID%in%inepic)
  
  # generate contour line 
  distdt = borderdist %>% filter(cellID%in%dt13_plot$cellID)
  distdt$dst = distdt$downDist
  
  crd = dt13_plot[,c('x','y')]
  bdlist = list()
  for(i in 1:10){
    print(i)
    tmp = crd[distdt$cellID[distdt$relativeDist>0.1*(i-1)],]
    as = alphahull::ashape(tmp,alpha = 120)
    plot(as)
    print('wait')
    bdlist[[i]] = extractpoly(as)[[1]] %>% data.frame()
    colnames(bdlist[[i]]) = c('x','y')
    bdlist[[i]]$group = paste0('g',i)
  }
  bds = Reduce(rbind,bdlist)
  
  plotcluster(dt13_plot,groupby = 'celltype',color = epicol1, bgcol = 'white',pt.size = 0.8,gridline = T,gridcol = 'black',zoomin = roi)+
    geom_path(data=bds,aes(group=group),color='black',linetype='dashed')
  ggsave('stp2_epi/plot/epiProgram_spatial2/T1_L/celltype_spatial.pdf',width = 6,height = 6)
}




## t5 lgin fov #####
{# t5 lgin fov
  dir.create('stp2_epi/plot/epiProgram_spatial2/T5_L_v2/',recursive = T)
  dt52 = allmeta_m %>% filter(sid=='T5S2') %>% filter(primaryType=='Epi')
  roi = c(2500,3500,2930,3530)
  plotcluster(dt52,'celltype',color = epicol1,xintercept = roi[1:2],yintercept = roi[3:4],axis = T,zoomin = roi,pt.size = 1) %>% plotly::ggplotly()
  ancpt = c('T5S2_djllofgl-1')
  
  # rotate coordinates
  crd = dt52[,c('x','y')]
  rotdeg = -13
  rotrad = rotdeg/180*pi
  wt1 = sin(rotrad)
  wt2 = cos(rotrad)
  newcrd = crd
  newcrd$x = crd$x*wt2-crd$y*wt1
  newcrd$y = crd$x*wt1+crd$y*wt2
  
  dt52_new = dt52
  dt52_new[,c('x','y')] = newcrd
  roi = dt52_new[ancpt,c('x','y')] %>% unlist()
  roi = c(roi[1],roi[1]+1500,roi[2]-250,roi[2]+750)
  plotcluster(dt52_new,'celltype',color = epicol1,xintercept = roi[1:2],yintercept = roi[3:4],axis = T,zoomin = roi,pt.size = 1)
  dt52_plot = dt52_new %>% filter(x<=roi[2],x>=roi[1],y<=roi[4],y>=roi[3])
  
  # generate contour line 
  distdt = borderdist %>% filter(cellID%in%dt52_plot$cellID)
  distdt$dst = distdt$downDist
  
  crd = dt52_plot[,c('x','y')]
  bdlist = list()
  for(i in 1:10){
    print(i)
    tmp = crd[distdt$cellID[distdt$relativeDist>0.1*(i-1)],]
    as = alphahull::ashape(tmp,alpha = 120)
    plot(as)
    print('wait')
    bdlist[[i]] = extractpoly(as)[[1]] %>% data.frame()
    colnames(bdlist[[i]]) = c('x','y')
    bdlist[[i]]$group = paste0('g',i)
  }
  bds = Reduce(rbind,bdlist)
  
  plotcluster(dt52_plot,groupby = 'celltype',color = epicol1, bgcol = 'white',pt.size = 0.8,gridline = T,gridcol = 'black',zoomin = roi)+
    geom_path(data=bds,aes(group=group),color='black',linetype='dashed')
  ggsave('stp2_epi/plot/epiProgram_spatial2/T5_L_v2/celltype_spatial.pdf',width = 6,height = 6)
}







