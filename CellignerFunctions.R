SingleDimPlot<-function(data,dims,col.by = NULL,cols = NULL,pt.size = NULL,shape.by = NULL,
                        alpha.by = NULL,order = NULL,label = FALSE,repel = FALSE,label.size = 4,
                        cells.highlight = NULL,cols.highlight = '#DE2D26',sizes.highlight = 1,na.value = 'grey50',raster = NULL){
  if(is.null(pt.size) ){pt.size <- AutoPointSize(data = data, raster = raster)}
  rownames(data)<-NULL
  if ((nrow(x = data) > 1e5) & !isFALSE(raster)){
    message("Rasterizing points since number of points exceeds 100,000.",
            "\nTo disable this behavior set `raster=FALSE`")
  }
  if(is.null(raster)){raster<-(nrow(data)>1e5)}
  #raster <- raster %||% (nrow(x = data) > 1e5)
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  } else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    if(is.null(cols)){col.base="#C3C3C3"}else{
      col.base=cols[1]
    }
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = sizes.highlight,
      cols.highlight = cols.highlight,
      col.base = col.base,
      pt.size = pt.size
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(!is.na(x = data[, col.by]), data[, col.by]), ]
      }
    } else {
      order <- rev(x = c(
        order,
        setdiff(x = unique(x = data[, col.by]), y = order)
      ))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  } else {
    # col.index <- grep(pattern = col.by, x = colnames(x = data), fixed = TRUE)
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = '^\\d', x = col.by)) {
      # Do something for numbers
      col.by <- paste0('x', col.by)
    } else if (grepl(pattern = '-', x = col.by)) {
      # Do something for dashes
      col.by <- gsub(pattern = '-', replacement = '.', x = col.by)
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  if (!is.null(x = alpha.by) && !alpha.by %in% colnames(x = data)) {
    warning(
      "Cannot find alpha variable ",
      alpha.by,
      " in data, setting to NULL",
      call. = FALSE,
      immediate. = TRUE
    )
    alpha.by <- NULL
  }
  
  plot <- ggplot(data = data)
  plot <- if(isTRUE(x = raster)){
    plot + geom_scattermore(
      mapping = aes_string(
        x = dims[1],
        y = dims[2],
        color = paste0("`", col.by, "`"),
        shape = shape.by,
        alpha = alpha.by
      ),
      pointsize = pt.size
    )
  }else {
    plot + geom_point(
      mapping = aes_string(
        x = dims[1],
        y = dims[2],
        color = paste0("`", col.by, "`"),
        shape = shape.by,
        alpha = alpha.by
      ),
      size = pt.size
    )
  }
  plot <- plot +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(color = NULL, title = col.by) 
  
  if (label && !is.null(x = col.by)) {
    plot <- LabelClusters(
      plot = plot,
      id = col.by,
      repel = repel,
      size = label.size
    )
  }
  if (!is.null(cols)) {
    if (length(cols) == 1 && (is.numeric(cols) || cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_color_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 && (cols %in% c('alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped'))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), palette = cols)
      scale <- scale_color_manual(values = colors, na.value = na.value)
    } else {
      scale <- scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + theme_cowplot()
  return(plot)
}

cell_line_tumor_class <- function(x, dist_mat, ann_mat, k=25, decreasing = T) {
  
  names(x) <- rownames(dist_mat)
  
  x <- sort(x, decreasing = decreasing)
  tumor_type <- dplyr::filter(ann_mat, sampleID %in% names(x[1:k]))$lineage %>%
    table() %>%
    as.data.frame() %>%
    dplyr::arrange(dplyr::desc(Freq)) %>%
    head(1) %>%
    .[['.']] %>% as.character()
  
  return(tumor_type)
}

get_cell_line_tumor_class <- function(tumor_CL_cor, alignment,removeUnknown=TRUE) {
  if(removeUnknown){
    unknownT<-alignment[alignment$lineage=="unknown","sampleID"]
    tumor_CL_cor<-tumor_CL_cor[setdiff(rownames(tumor_CL_cor),unknownT),]
    tumor_CL_cor<-tumor_CL_cor[,setdiff(colnames(tumor_CL_cor),unknownT)]
  }
  cl_tumor_classes <- apply(tumor_CL_cor, 2, function(x) cell_line_tumor_class(x, tumor_CL_cor, alignment)) %>% 
    as.character()
  names(cl_tumor_classes) <- colnames(tumor_CL_cor)
  
  return(cl_tumor_classes)
}
cell_line_tumor_class_plot <- function(cl_tumor_classes, alignment, tumor_CL_cor, filename) {
  cl_tissue_type <- dplyr::filter(alignment, type=='CL')
  #cl_tissue_type[grep('rhabdomyosarcoma', cl_tissue_type$subtype),'tissue'] <- 'rhabdomyosarcoma'
  rownames(cl_tissue_type) <- cl_tissue_type$sampleID
  classification_freq <- table(cl_tumor_classes, cl_tissue_type[colnames(tumor_CL_cor),'lineage']) %>% as.data.frame()
  classification_freq <- reshape2::dcast(classification_freq, cl_tumor_classes ~ Var2, value.var = 'Freq') %>%
    tibble::column_to_rownames('cl_tumor_classes')
  print(setdiff(intersect(unique(dplyr::filter(alignment, type=='CL')$lineage),
                          unique(dplyr::filter(alignment, type=='tumour')$lineage)), 
                rownames(classification_freq)))
  
  esophagus_tumor <- rep(0, ncol(classification_freq))
  thyroid_tumor <- rep(0, ncol(classification_freq))
  classification_freq <- rbind(classification_freq,`esophagus`= esophagus_tumor,`thyroid`=thyroid_tumor) 
  common_types <- intersect(rownames(classification_freq), colnames(classification_freq))
  
  prop_agree <- sum(diag(as.matrix(classification_freq[common_types, common_types])))/sum(as.matrix(classification_freq[common_types, common_types]))
  print(prop_agree)
  for(i in 1:ncol(classification_freq)) {
    classification_freq[,i] <- classification_freq[,i]/sum(classification_freq[,i])
  }
  
  
  agreement <- diag(as.matrix(classification_freq[common_types, common_types]))
  agreement_CL <- agreement
  names(agreement_CL) <- common_types
  agreement_tumor <- agreement
  names(agreement_tumor) <- common_types
  
  agreement_CL <- base::sort(agreement_CL, decreasing=T)
  agreement_tumor <- base::sort(agreement_tumor, decreasing=T)
  
  
  
  
  classification_freq <- classification_freq[names(agreement_tumor), names(agreement_CL)]
  rownames(classification_freq) <- gsub("_", " ", rownames(classification_freq))
  colnames(classification_freq) <- gsub("_", " ", colnames(classification_freq))
  classification_freq<-classification_freq[setdiff(rownames(classification_freq),"unknown"),setdiff(colnames(classification_freq),"unknown")]
  pheatmap::pheatmap(classification_freq, 
                     border_color = heatmap_params$square_border_color, 
                     na_col= heatmap_params$na_color, 
                     cluster_rows = F, 
                     cluster_cols = F, 
                     main="", 
                     fontsize = heatmap_params$title_font_size,
                     fontsize_col = heatmap_params$column_font_size,
                     fontsize_row = heatmap_params$row_font_size,
                     width = 3.5,
                     height = 3,
                     fontface = heatmap_params$font_face,
                     angle_col=90, 
                     filename = filename,
                     color= heatmap_params$color_vector)
  return(classification_freq)
}
heatmap_params <- list(
  na_color = '#666666',
  title_font_size = 6,
  row_font_size = 6,
  column_font_size = 6,
  font_face = "plain",
  square_border_color = "white",
  color_palette = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'YlOrRd')), space='Lab'),
  color_vector = c('#e3e3e3', rev(grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'YlOrRd')), space='Lab')(100)))
)
calc_tumor_CL_cor <- function(Celligner_aligned_data, Celligner_info) {
  tumors_samples <- intersect(dplyr::filter(Celligner_info, type=='tumour')$sampleID,colnames(Celligner_aligned_data))
  cl_samples <- intersect(dplyr::filter(Celligner_info, type=='CL')$sampleID,colnames(Celligner_aligned_data))
  print(length(tumors_samples))
  print(length(cl_samples))
  tumor_CL_cor <- cor(Celligner_aligned_data[,tumors_samples], Celligner_aligned_data[,cl_samples],
                      use='pairwise')
  
  
  return(tumor_CL_cor)
}

