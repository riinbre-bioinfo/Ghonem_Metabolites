# script to create strip plots for the metabolite data, run ANOVA for the metabolite data and create bar plots for enzyme data

metabsStripPlot <- function(data_meta.df, data_log.df, xvar, shape, color, grouped = FALSE, shapes = FALSE, colors = FALSE, brewer = FALSE, legend = TRUE) {

  metabNames = names(data_log.df)
  data.df = cbind(data_meta.df, data_log.df)

  stripPlot <- function(metab, data, xvar, shape, color, grouped, shapes, colors, brewer) {  
    p = ggplot(data, aes_(x = as.name(xvar), y = as.name(metab), shape = as.name(shape), color = as.name(color))) + theme_bw()

    if (legend == FALSE) {
      p = p + theme(legend.position = "none")
    }
      
    if (is.vector(shapes)) {
      p = p + scale_shape_manual(values = shapes)
    }
  
    if (is.character(colors) && brewer == FALSE) {
      p = p + scale_color_manual(values = colors)
    } else if (is.character(colors) && brewer == TRUE) {
      p = p + scale_color_brewer(palette = colors)
    } else {
      p = p + scale_color_brewer(palette="Dark2")
    }
  
    if (grouped == TRUE) {
      p = p + 
        geom_boxplot(position = position_dodge(0.5)) + 
        geom_jitter(position = position_dodge(0.5), size = 3)
    } else {
      p = p + 
        geom_jitter(position = position_jitter(0.2), size = 3)
    }
  
    p = p + 
      stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = 0.5, size = 1.2) +
      stat_summary(fun.y=mean, geom = "point", color = "black", size = 4, pch = 1) +
      theme(text = element_text(face = "bold", size = 18), axis.text = element_text(face = "bold", size = 14))
  }
  
  lapply(metabNames, stripPlot, data = data.df, xvar = xvar, shape = shape, color = color, grouped = grouped, shapes = shapes, colors = colors, brewer = brewer)
}

metabANOVA <- function(data_meta.df, data_log.df, fixed, random) { # run ANOVA and multicomp Tukey on all metabolites
  metabNames = names(data_log.df)
  data.df = cbind(data_meta.df, data_log.df)
  runANOVA <- function(metab, data, fixed, random) {
    form = as.formula(paste(metab, paste(fixed, random, sep = "+"), sep = "~")) # construct formula
    data.aov = aov(form, data = data) # run anova
    args = list("Tukey") # next commands pass "Treatment" to mcp
    names(args) = fixed
    cmp = do.call(mcp, args)
    data.glht = glht(data.aov, linfct = cmp) # run glht
    summary(data.glht)
  }

  summaries = lapply(metabNames, runANOVA, data = data.df, fixed = fixed, random = random) # return list of anova summaries
  names(summaries) <- names(data_log.df)
  summaries
}

enzymesBarPlot <- function(enzyme, data.df, group) {

  df.summ <- data.df %>% group_by(get(group)) %>% summarize(Mean = mean(get(enzyme)), se = sd(get(enzyme))/sqrt(n()))
  names(df.summ)[1] <- group
  df.summ
    
#  p <- ggplot(df.summ, aes(x = get(group), y = Mean, fill = get(group))) +
#    geom_bar(stat = "identity") +
#    stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = 0.25)
  
#  p <- ggplot(df.summ, aes(x = TreatCohort, y = Mean, ymin = Min, ymax = Max, fill = TreatCohort)) +
#    geom_bar(stat = "identity")
  
#  enz_summaries <- lapply(enzymes, function(x) tapply(data.df[, x], data.df[, group], mean))
#  names(enz_summaries) <- enzymes
#  means.df <- rownames_to_column(as.data.frame(enz_summaries), var = group)
  
#  barPlot <- function(enzyme, data, xvar) {  
#    p = ggplot(data.df[, enzymes], aes_(x = as.name(xvar), y = as.name(enzyme))) + theme_bw()
    
#    p = p + geom_bar(stat("identity"))
    
#    p = p + stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = 0.25)
#  }
  
#  lapply(enzymes, barPlot, data = data.df, xvar = "TreatCohort")
  
#  means.df
}