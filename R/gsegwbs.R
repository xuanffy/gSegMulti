# Multiple Change-point Search Using Graph-based Method with Wild Binary Segmentation
gsegwbs = function(s = 1,
                       e,
                       alpha = 0.05,
                       L = 100,
                       cutoff = 0.05,
                       g.max = 30,
                       MinLen = 7,
                       surp.error = FALSE,
                       g.exp = 0.5,
                       dist.type = 'euclidean',
                       keep.dist = FALSE) {
  main = FALSE
  graph = 'mst'
  stat = 'g'
  argchar = 'generalized'
  gfun = function(m) {
    floor(m ^ g.exp)
  }
  n = e - s + 1
  if (n < MinLen) {
    return(NULL)
  } else{
    len = e - s + 1
    FTM = matrix(0, ncol = 2, nrow = 0)
    i = len - MinLen + 1
    for (l in s:(e - MinLen + 1)) {
      FTM = rbind(FTM, matrix(c(rep(l, i), (l + MinLen - 1):e), ncol = 2, nrow =
                                i))
      i = i - 1
    }
    if (L < ((len - MinLen + 1) * (len - MinLen + 2) / 2)) {
      FTM = FTM[sample(dim(FTM)[1], L, replace = F), ]
    }
  }
  all.pvalue = numeric(dim(FTM)[1])
  all.position = numeric(dim(FTM)[1])
  if (!exists('ydist', '.GlobalEnv')) {
    main = TRUE
    ydist <<- as.matrix(dist(y, method = dist.type))
  }
  if (graph == 'mst') {
    for (l in 1:dim(FTM)[1]) {
      m = FTM[l, 2] - FTM[l, 1] + 1
      if (surp.error) {
        r = try(gseg1(
          m,
          mstree(as.dist(ydist[FTM[l, 1]:FTM[l, 2], FTM[l, 1]:FTM[l, 2]]), ngmax =
                   min(gfun(m), g.max)),
          statistics = stat,
          n0 = cutoff * m,
          n1 = (1 - cutoff) * m
        ))
        if (inherits(r, "try-error")) {
          all.pvalue[l] = 1
          all.position[l] = 1
        } else{
          all.pvalue[l] = eval(parse(text = paste('r$pval.appr$', argchar, sep = '')))
          all.position[l] = eval(parse(text = paste(
            'r$scanZ$', argchar, '$tauhat', sep = ''
          )))
        }
      } else{
        r = gseg1(
          m,
          mstree(as.dist(ydist[FTM[l, 1]:FTM[l, 2], FTM[l, 1]:FTM[l, 2]]), ngmax =
                   min(gfun(m), g.max)),
          statistics = stat,
          n0 = cutoff * m,
          n1 = (1 - cutoff) * m
        )
        all.pvalue[l] = eval(parse(text = paste('r$pval.appr$', argchar, sep =
                                                  '')))
        all.position[l] = eval(parse(text = paste(
          'r$scanZ$', argchar, '$tauhat', sep = ''
        )))
      }
    }
  }
  
  if (graph == 'mst') {
    if (surp.error) {
      r = try(gseg1(
        n,
        mstree(as.dist(ydist[s:e, s:e]), ngmax = min(gfun(n), g.max)),
        statistics = stat,
        n0 = cutoff * n,
        n1 = (1 - cutoff) * n
      ))
      if (inherits(r, "try-error")) {
        lastpv = 1
        lasttau = 1
      } else{
        lastpv = eval(parse(text = paste('r$pval.appr$', argchar, sep = '')))
        lasttau = eval(parse(text = paste(
          'r$scanZ$', argchar, '$tauhat', sep = ''
        )))
      }
    } else{
      r = gseg1(
        n,
        mstree(as.dist(ydist[s:e, s:e]), ngmax = min(gfun(n), g.max)),
        statistics = stat,
        n0 = cutoff * n,
        n1 = (1 - cutoff) * n
      )
      lastpv = eval(parse(text = paste('r$pval.appr$', argchar, sep = '')))
      lasttau = eval(parse(text = paste(
        'r$scanZ$', argchar, '$tauhat', sep = ''
      )))
    }
  }
  
  
  bestindex = which.min(all.pvalue)
  if (min(all.pvalue) > lastpv) {
    if (lastpv < alpha) {
      newtau = lasttau
      tauhat <<- c(tauhat, newtau + s - 1)
      
      gsegwbs(
        s = s,
        e = newtau + s - 1,
        alpha = alpha,
        L = L,
        cutoff = cutoff,
        g.max = g.max,
        MinLen = MinLen,
        surp.error = surp.error,
        g.exp = g.exp
      )
      gsegwbs(
        s = newtau + s,
        e = e,
        alpha = alpha,
        L = L,
        cutoff = cutoff,
        g.max = g.max,
        MinLen = MinLen,
        surp.error = surp.error,
        g.exp = g.exp
      )
    }
  } else{
    if (min(all.pvalue) < alpha) {
      newtau = all.position[bestindex]
      tauhat <<- c(tauhat, newtau + FTM[bestindex, 1] - 1)
      gsegwbs(
        s = s,
        e = newtau + FTM[bestindex, 1] - 1,
        alpha = alpha,
        L = L,
        cutoff = cutoff,
        g.max = g.max,
        MinLen = MinLen,
        surp.error = surp.error,
        g.exp = g.exp
      )
      gsegwbs(
        s = newtau + FTM[bestindex, 1],
        e = e,
        alpha = alpha,
        L = L,
        cutoff = cutoff,
        g.max = g.max,
        MinLen = MinLen,
        surp.error = surp.error,
        g.exp = g.exp
      )
    }
  }
  if (main == TRUE & keep.dist == FALSE) {
    rm(ydist, pos = '.GlobalEnv')
  }
  return(NULL)
}
