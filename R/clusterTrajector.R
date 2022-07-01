#' Cluster the nonlinear trajector
#' @description
#'
#' @param dataset dataframe,microbe or lipid datset, row:sample id ,col: feature
#' @param phemeta dataframe,phenotye dataset
#' @param metavar character,the intersting feature
#' @param scale  optinal character, log or int
#' @param span   numerical value
#' @param na.rm  logical value indicating whether NA values should be stripped before the computation proceeds.
#' @param cut   numerical value
#' @param ...  urther arguments passed to or from other methods
#'
#' @return list
#' @export
#'
#' @examples
#'
#'
clusterTraj <- function(dataset, phemeta, metavar , scale = "log", span= 0.7, na.rm = T, cut,...){

  # rm the NA
  if(na.rm){
    phemeta <- phemeta[!is.na(phemeta[, metavar]), ,drop=F]
  }
  # to rm the outlier of metavar
  if(cut){
    phemeta <- cutQuan(phemeta, metavar, quan = quan)
  }

  # mathch the sample id
  id <- intersect(rownames(dataset), rownames(phemeta))
  dataset2 <- dataset[id, ]
  value <- phemeta[id,]
  # to scale the log10 value
  if(scale == "log"){
    datasetScale <- apply(dataset2, 2, function(x){scale(log10(x+1))})
  }else if(scale == "int"){
    datasetScale <- apply(dataset2, 2, function(x){INT(x)})
  }else{
    datasetScale <- dataset2
  }

  # to compute the mdian value on the same X
  datasetScale2 <- apply(datasetScale, 2,
                         function(x){tapply(x, value, median, na.rm=T)})

  # loess estimate
  value <- as.numeric(rownames(datasetScale2))
  datasetLoess <- apply(datasetScale2, 2,
                        function(x){predict(stats::loess(x~value, span = span))})

  # cluster
  # library(NbClust)
  # here paremete can to modify
  hc <- stats::hclust(dist(t(datasetLoess), method = "euclidean"), method = "ward.D2")
  # get the best clust
  cluster <- NbClust::NbClust(t(datasetLoess), distance = "euclidean", min.nc=2, max.nc=8,
                     method = "complete", index = "ch")

  bestC <- as.numeric(cluster$Best.nc[1])
  clustlable <- cutree(hc, bestC)

  hashtmp <- hash::hash(names(clustlable), as.numeric(clustlable))

  # plot
  datasetLoess <- as.data.frame(datasetLoess, check.names=F)
  datasetLoess$var <- sort(value)
  qdat <- reshape2::melt(datasetLoess, id.vars = "var")
  qdat$group <- sapply(qdat$variable, function(x){values(hashtmp[x])})


  finalTheme <- theme_set(theme_bw()) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    )
  fontTheme <- theme(
    axis.title=element_text(size=14, face = "bold", colour = "black"),
    text=element_text(size=12),
    legend.text = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    plot.title = element_text(size = 14, face = "bold", colour = "black"),
    legend.title = element_text(size = 14, face = "bold", colour = "black"),
    strip.text = element_text(face = "bold", size = 14)
  )

  figure <- ggplot2::ggplot(qdat, mapping = aes(x = var, y = value, group=variable))+
            geom_line()+fontTheme+finalTheme+facet_wrap(.~group, scale="free")

  return(list(figure, as.data.frame(clustlable)))
}

cutQuan <- function(phe, metavar, quan){
  phe <- phe[!is.na(phe[,metavar]), ,drop=F]
  value <- phe[,metavar]
  sortvalue <- sort(value)
  minquan <- sortvalue[round(length(sortvalue)*quan)]
  maxquan <- sortvalue[round(length(sortvalue)*(1-quan))]
  index <- which(value>minquan & value<maxquan)
  return(phe[index, ,drop=F])
}

INT <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}



