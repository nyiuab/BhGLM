
# construct main-effect variable matrix

make.main <- function (geno, model = c("Cockerham", "codominant", "additive", "dominant", "recessive", "overdominant"),
                       fill.missing = TRUE, ind.group = NULL, geno.order = TRUE,
                       loci.names = c("marker", "position"), imprint = TRUE, verbose = FALSE, ...)
{
  model <- model[1]
  if(!any(model == c("Cockerham","codominant","additive","dominant","recessive","overdominant","full"))){
    model <- c("Cockerham")
    warning("You give wrong model. use Cockerham.", call. = FALSE, immediate. = TRUE)
  }
  loci.names <- loci.names[1]
  if(!any(loci.names == c("marker","position"))) {
    loci.names <- c("marker")
    warning("You give wrong loci name. use marker names.", call. = FALSE, immediate. = TRUE)
  }

#for crosses qtl mapping
  if(length(class(geno)) == 2){
    require(qtlbim)
    if(class(geno)[2] != "cross") stop("not an object of class 'cross'")
    cross = geno
    for(i in 1:nchr(cross)) {
      geno.m = apply(cross$geno[[i]]$data, 2, max, na.rm = TRUE)
      if(class(cross)[1] == "f2") m = 3
      else m = 2
      if(any(geno.m < m)) {
        geno.0 = geno.m[geno.m < m]
        cross = drop.markers(cross, names(geno.0))
        warning("markers with incomplete genotypes removed!", call. = FALSE, immediate. = TRUE)
      }
    }
    n.chr = nchr(cross)
    n.ind = nind(cross)
    n.mar = nmar(cross)
    x = vector(mode="list", length = n.chr)
    names(x) = names(cross$geno)

    if( is.null(cross$geno[[1]]$prob) ) {
      warning("first run qb.genoprob(cross, step = 0)", call. = FALSE, immediate. = TRUE)
      cross = qb.genoprob(cross, step = 0)
    }
    else if(attr(cross$geno[[1]]$prob, "step") != 0) cross = qb.genoprob(cross, step = 0)
    grid = pull.loci(cross) # in order to deal with loci in marker intervals

    for(i in 1:n.chr)
    {
      if(class(grid[[i]])=="matrix") grid[[i]] = grid[[i]][1, ]
      grid[[i]] = round(grid[[i]],digits = 1)
      if(loci.names == "marker") names(grid[[i]]) = paste(names(grid[[i]]))
      if(loci.names == "position")
        names(grid[[i]]) = paste("c", names(grid)[i], ".", grid[[i]], sep = "")

      if(dim(cross$geno[[i]]$prob)[3] == 2 | class(cross$geno[[i]]) == "X"){
        x[[i]] = as.matrix(0.5 * (cross$geno[[i]]$prob[, , 2] - cross$geno[[i]]$prob[, , 1]))
        colnames(x[[i]]) = names(grid[[i]])
      }

      if(dim(cross$geno[[i]]$prob)[3] == 3 & class(cross$geno[[i]]) != "X"){
        if (model == "codominant") {
          x1 = as.matrix(cross$geno[[i]]$prob[, , 3])
          x2 = as.matrix(cross$geno[[i]]$prob[, , 2])
          colnames(x1) = paste(names(grid[[i]]), "r", sep = "")
          colnames(x2) = paste(names(grid[[i]]), "h", sep = "")
          x[[i]] = cbind(x1, x2)
        }
        if (model == "Cockerham") {
          x1 = as.matrix(cross$geno[[i]]$prob[, , 3] - cross$geno[[i]]$prob[, , 1])
          x2 = as.matrix(0.5 * (cross$geno[[i]]$prob[, , 2] - cross$geno[[i]]$prob[, , 1] - cross$geno[[i]]$prob[, , 3]))
          colnames(x1) = paste(names(grid[[i]]), "a", sep = "")
          colnames(x2) = paste(names(grid[[i]]), "d", sep = "")
          x[[i]] = cbind(x1, x2)
        }
        if (model == "additive") {
          x[[i]] = as.matrix(cross$geno[[i]]$prob[, , 3] - cross$geno[[i]]$prob[, , 1])
          colnames(x[[i]]) = names(grid[[i]])
        }
        if (model == "recessive") {
          x[[i]] = as.matrix(cross$geno[[i]]$prob[, , 1] + cross$geno[[i]]$prob[, , 2])
          colnames(x[[i]]) = names(grid[[i]])
        }
        if (model == "dominant") {
          x[[i]] = as.matrix(cross$geno[[i]]$prob[, , 3] + cross$geno[[i]]$prob[, , 2])
          colnames(x[[i]]) = names(grid[[i]])
        }
        if (model == "overdominant") {
          x[[i]] = as.matrix(cross$geno[[i]]$prob[, , 2])
          colnames(x[[i]]) = names(grid[[i]])
        }
      }

      if(dim(cross$geno[[i]]$prob)[3] == 4 & class(cross$geno[[i]]) != "X"){
        x1 = as.matrix(cross$geno[[i]]$prob[, , 4] - cross$geno[[i]]$prob[, , 1])
        x2 = as.matrix(0.5 * (cross$geno[[i]]$prob[, , 2] + cross$geno[[i]]$prob[, , 3] - cross$geno[[i]]$prob[, , 1] - cross$geno[[i]]$prob[, , 4]))
        x3 = as.matrix(0.5 * (cross$geno[[i]]$prob[, , 3] - cross$geno[[i]]$prob[, , 2]))
        colnames(x1) = paste(names(grid[[i]]), "a", sep = "")
        colnames(x2) = paste(names(grid[[i]]), "d", sep = "")
        colnames(x3) = paste(names(grid[[i]]), "i", sep = "")
        x[[i]] = cbind(x1, x2)
        if (imprint) x[[i]] = cbind(x[[i]], x3)
      }

      if (ncol(x[[i]]) > n.mar[i]) {
        o = order(rep(1:n.mar[i], as.integer(ncol(x[[i]])/n.mar[i])))
        x[[i]] = x[[i]][, o]
      }
    }

    x.main = x[[1]]
    if(nchr(cross) > 1)
      for(j in 2:nchr(cross)) x.main = data.frame(x.main, x[[j]])

    return(x.main)
  }

#for human association analysis
  if(length(class(geno)) == 1) {

    geno = as.matrix(geno)
    if(is.null(colnames(geno))) colnames(geno) = paste("m", 1:ncol(geno), sep = "")
    z = geno
    func = function(w){ # recode genotype as 1,2,3
      w = as.factor(w)
      levels(w) = as.character(1:length(levels(w)))
      w = as.numeric(w)
      return(w)
    }
    z = apply(z, 2, func)

    w = apply(z, 2, var, na.rm = TRUE)
    if(length(which(w == 0)) != 0){
      z = z[ , which(w != 0), drop = FALSE]
      d = which(w == 0)
      warning(length(d), " markers with only one genotype are removed!", call. = FALSE, immediate. = TRUE )
      if(verbose)
        for(j in 1:length(d)) cat(names(d)[j], "\n")
    }
    
    func = function(w){
      na.prop = length(which(is.na(w)))/length(w)
      return(na.prop)
    }
    na.prop = apply(z, 2, func)
    if(any(na.prop > 30/100)) {
      d <- which(na.prop > 30/100)
      warning(length(d)," markers have more than 30% genotypes!", call. = FALSE, immediate. = TRUE)  
      if(verbose)
        for(j in 1:length(d)) cat(names(d)[j], "\n")
    }

    geno.m = apply(z, 2, max, na.rm = TRUE)
    if(any(geno.m < 3)) {
      two = names(which(geno.m < 3))
      warning(length(geno.m[geno.m < 3])," markers have only two genotypes!", call. = FALSE, immediate. = TRUE)
      if(verbose)
        for(j in 1:length(two)) cat(two[j], "\n")
      if(model == "Cockerham")
        warning("Cockerham model cannot be used to markers with only two genotypes!", call. = FALSE, immediate. = TRUE)
      if(model == "codominant")
        warning("Codominant model cannot be used to markers with only two genotypes!", call. = FALSE, immediate. = TRUE)
    }
    if(any(geno.m > 3)) {
      z = z[ , which(geno.m <= 3), drop = FALSE]
      d = which(geno.m > 3)
      warning(length(d), " markers with more than three genotypes are removed!", call. = FALSE, immediate. = TRUE)
      if(verbose)
        for(j in 1:length(d)) cat(names(d)[j], "\n")
    }
    if(ncol(z) == 0) stop("All markers have been removed!")

    if(geno.order) { # order genotypes: 1=c, 2=h, 3=r
      func = function(x) {
        w = table(x)
        x = (w[1] < w[length(w)]) * (length(w) + 1 - x) + (w[1] >= w[length(w)]) * x
      }
      z = apply(z, 2, func)
    }

    geno.m = apply(z, 2, max, na.rm = TRUE)
    z2 = z[ , which(geno.m == 2), drop = FALSE]
    z3 = z[ , which(geno.m == 3), drop = FALSE]
    x = NULL
    if(ncol(z3) > 0){
      if(model == "codominant") {
        x1 = 1 * (z3 == 3)
        x2 = 1 * (z3 == 2)
        colnames(x1) = paste(colnames(z3), "r", sep = "")
        colnames(x2) = paste(colnames(z3), "h", sep = "")
        x = cbind(x1, x2)
        o = order(rep(1:ncol(z3), 2))
        x = x[ , o]
      }
      if(model == "Cockerham") {
        x1 = z3 - 2
        x2 = 0.5 - abs(x1)
        colnames(x1) = paste(colnames(z3), "a", sep = "")
        colnames(x2) = paste(colnames(z3), "d", sep = "")
        x = cbind(x1, x2)
        o = order(rep(1:ncol(z3), 2))
        x = x[ , o]
      }
      if(model == "recessive") x = 1 * (z3 == 3)
    }
    x = cbind(z2 - 1, x)
    if(model == "additive") x = z - 1
    if(model == "dominant") x = 1 * (z != 1)
    if(model == "overdominant") x = 1 * (z == 2)


    if(fill.missing){
      if(!is.null(ind.group)) {
        ind.group = as.factor(ind.group)
        if(nrow(geno) != length(ind.group)) {
          warning("geno and ind.group have different obs. Cannot use group information!", call. = FALSE, immediate. = TRUE)
          ind.group = NULL
        }
      }
      if(is.null(ind.group)) {
        func = function(w){
          na.index = which(is.na(w))
          w[na.index] = mean(w, na.rm = TRUE)
          return(w)
        }
        x = apply(x, 2, func)
      }
      if(!is.null(ind.group)){
        func = function(w){
          group.means = tapply(w, ind.group, mean, na.rm = TRUE)
          group.means = ifelse(is.na(group.means), mean(w, na.rm = TRUE), group.means)
          for(k in 1:length(group.means)){
            na.index = which(is.na(w) & ind.group == names(group.means)[k])
            w[na.index] = mean(w[ind.group == names(group.means)[k]], na.rm = TRUE)
          }
          return(w)
        }
        x = apply(x, 2, func)
      }
    }

#    if(ncol(x)>ncol(z)) {
#      o = order(rep(1:ncol(z), as.integer(ncol(x)/ncol(z))))
#      x = x[,o]
#    }
    x.main = data.frame(x)
    return(x.main)
  }

}