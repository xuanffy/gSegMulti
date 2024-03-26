# Multiple Change-point Search Using Graph-based Method with Seeded Binary Segmentation
gsegsbs = function(s = 1,
                          e,
                          alpha = 0.05,
                          cutoff = 0.05,
                          g.max = 30,
                          MinLen = 7,
                          g.exp = 0.5,
                          dist.type = 'euclidean',
                          decay.a = sqrt(0.5),
                          keep.dist = FALSE) {
  main = FALSE
  stat = 'g'
  argchar = 'generalized'
  gfun = function(m) {
    floor(m ^ g.exp)
  }
  n = e - s + 1
  if (n < MinLen) {
    return(NULL)
  }
  if (!exists('FTM', '.GlobalEnv')) {
    main = TRUE
    FTM <<- matrix(c(s, e), ncol = 2, byrow = TRUE)
    all.pvalue <<- numeric(dim(FTM)[1])
    all.position <<- numeric(dim(FTM)[1])
    for (k in 2:ceiling(log(e, base = 1 / decay.a))) {
      lk = e * decay.a ^ (k - 1)
      if (lk >= MinLen) {
        nk = 2 * ceiling((1 / decay.a) ^ (k - 1)) - 1
        sk = (e - lk) / (nk - 1)
        Ik = matrix(c(floor((0:(
          nk - 1
        )) * sk), ceiling((0:(
          nk - 1
        )) * sk + lk)),
        ncol = 2,
        byrow = FALSE)
        FTM <<- rbind(FTM, Ik)
      }
    }
    FTM[, 1] <<- ifelse(FTM[, 1] >= 1, FTM[, 1], 1)
    FTM[, 2] <<- ifelse(FTM[, 2] <= e, FTM[, 2], e)
    
    ydist <<- as.matrix(dist(y, method = dist.type))
    for (l in 1:dim(FTM)[1]) {
      m = FTM[l, 2] - FTM[l, 1] + 1
      r = gseg1(
        m,
        mstree(as.dist(ydist[FTM[l, 1]:FTM[l, 2], FTM[l, 1]:FTM[l, 2]]), ngmax =
                 min(gfun(m), g.max)),
        statistics = stat,
        n0 = cutoff * m,
        n1 = (1 - cutoff) * m
      )
      all.pvalue[l] <<-
        eval(parse(text = paste('r$pval.appr$', argchar, sep = '')))
      all.position[l] <<-
        eval(parse(text = paste(
          'r$scanZ$', argchar, '$tauhat', sep = ''
        )))
    }
    
  }
  
  Mse = which((FTM[, 1] >= s) & FTM[, 2] <= e)
  pvalue = all.pvalue[Mse]
  r = gseg1(
    n,
    mstree(as.dist(ydist[s:e, s:e]), ngmax = min(gfun(n), g.max)),
    statistics = stat,
    n0 = cutoff * n,
    n1 = (1 - cutoff) * n
  )
  bestindex = which.min(pvalue)
  if (min(pvalue) > eval(parse(text = paste('r$pval.appr$', argchar, sep =
                                            '')))) {
    if (eval(parse(text = paste('r$pval.appr$', argchar, sep = ''))) < alpha) {
      newtau = eval(parse(text = paste(
        'r$scanZ$', argchar, '$tauhat', sep = ''
      )))
      tauhat <<- c(tauhat, newtau + s - 1)
      gsegsbs(
        s = s,
        e = newtau + s - 1,
        alpha = alpha,
        cutoff = cutoff,
        g.max = g.max,
        MinLen = MinLen,
        decay.a = decay.a,
        g.exp = g.exp
      )
      gsegsbs(
        s = newtau + s,
        e = e,
        alpha = alpha,
        cutoff = cutoff,
        g.max = g.max,
        MinLen = MinLen,
        decay.a = decay.a,
        g.exp = g.exp
      )
    }
  } else{
    if (min(pvalue) < alpha) {
      newtau = all.position[Mse[bestindex]]
      tauhat <<- c(tauhat, newtau + FTM[Mse[bestindex], 1] - 1)
      gsegsbs(
        s = s,
        e = newtau + FTM[Mse[bestindex], 1] - 1,
        alpha = alpha,
        cutoff = cutoff,
        g.max = g.max,
        MinLen = MinLen,
        decay.a = decay.a,
        g.exp = g.exp
      )
      gsegsbs(
        s = newtau + FTM[Mse[bestindex], 1],
        e = e,
        alpha = alpha,
        cutoff = cutoff,
        g.max = g.max,
        MinLen = MinLen,
        decay.a = decay.a,
        g.exp = g.exp
      )
    }
  }
  if (main == TRUE) {
    rm(FTM, all.pvalue, all.position, pos = '.GlobalEnv')
    if (keep.dist == FALSE) {
      rm(ydist, pos = '.GlobalEnv')
    }
  }
  return(NULL)
}
