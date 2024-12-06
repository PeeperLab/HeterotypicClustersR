min_thres_circos <- function(vis_circos_obj, transparency = FALSE, args.circos.text = list(), link.visible = TRUE, ...) {
  if (!is.logical(transparency)) 
    stop("transparency should be a logical")
  
  if (!all(c("links_circle", "ligand_colors", "order", "gaps") %in% names(vis_circos_obj))) 
    stop("vis_circos_obj should contain the elements 'links_circle', 'ligand_colors', 'order', and 'gaps'")
  
  # Initialize transparency value
  transparency_val <- 0
  
  if (transparency) {
    threshold <- 0.1  # Set weight threshold
    transparency_val <- vis_circos_obj$links_circle %>%
      mutate(weight = (weight - min(weight)) / (max(weight) - min(weight))) %>%  # Apply sqrt transformation
      mutate(transparency = ifelse(weight < threshold, 0.9, 1 - weight)) %>%  # Fixed transparency for low weights
      .$transparency
  }
  
  default_params <- list(
    x = vis_circos_obj$links_circle,
    order = vis_circos_obj$order,
    grid.col = vis_circos_obj$ligand_colors,
    transparency = transparency_val,
    directional = 1,
    link.sort = TRUE,
    link.decreasing = FALSE,
    diffHeight = 0.005,
    direction.type = c("diffHeight", "arrows"),
    link.arr.type = "big.arrow",
    link.visible = link.visible,  # Override if necessary
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.075)
  )
  default_params[names(list(...))] <- list(...)
  
  circos_text_default_params <- list(facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
  circos_text_default_params[names(args.circos.text)] <- args.circos.text
  
  circos.par(gap.degree = vis_circos_obj$gaps)
  do.call(chordDiagram, default_params)
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    do.call(circos.text, c(list(x = CELL_META$xcenter, y = CELL_META$ylim[1], 
                                label = CELL_META$sector.index), circos_text_default_params))
  }, bg.border = NA)
  circos.clear()
}
