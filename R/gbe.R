# Change-point Pruning by Stepwise Method with Extended Pseudo BIC (ep-BIC)
gbe = function(y,
                    tauhat,
                    g.max = 5,
                    detail = FALSE,
                    penalty = FALSE,
                    dist.type = 'euclidean') {
  m = length(tauhat)
  if (m <= 1) {
    return(tauhat)
  }
  stat = 'g'
  argchar = 'generalized'
  
  max.iter = 1000
  if (penalty == TRUE) {
    pen = function(tau, T) {
      m = length(tau)
      return(-(2 * m) * log(T))
    }
  } else{
    pen = function(tau, T) {
      0
    }
  }
  tauhat = sort(tauhat)
  n = dim(y)[1]
  
  mergetau = numeric(0)
  gofstat = numeric(0)
  gofoverall = numeric(0)
  tauhat = c(0, tauhat, n)
  if (!exists('ydist', '.GlobalEnv')) {
    ydist <<- as.matrix(dist(y, method = dist.type))
  }
  for (i in 1:m) {
    r = g.tests(
      mstree(as.dist(ydist[(tauhat[i] + 1):tauhat[i + 2], (tauhat[i] + 1):tauhat[i +
                                                                                   2]]), min(floor(
                                                                                     sqrt(tauhat[i + 2] - tauhat[i])
                                                                                   ), g.max)),
      1:(tauhat[i + 1] - tauhat[i]),
      (tauhat[i + 1] - tauhat[i] + 1):(tauhat[i + 2] - tauhat[i]),
      stat
    )
    gofstat[i] = eval(parse(text = paste(
      'r$', argchar, '$test.statistic', sep = ''
    )))
  }
  gofoverall[1] = sum(gofstat) + pen(tauhat[c(-1, -length(tauhat))], n)
  change = numeric(0)
  maxind = NA
  if (m >= 3) {
    for (i in 1:min(c(m - 2, max.iter))) {
      oldpen = pen(tauhat[c(-1, -length(tauhat))], n)
      if (is.na(maxind)) {
        scanrange = 1:(length(tauhat) - 2)
      }
      if (1 %in% scanrange) {
        r = g.tests(
          mstree(as.dist(ydist[1:tauhat[4], 1:tauhat[4]]), min(floor(
            sqrt(tauhat[4])
          ), g.max)),
          1:tauhat[3],
          (tauhat[3] + 1):tauhat[4],
          test.type = stat,
          maxtype.kappa = maxtype.kappa
        )
        change[1] = eval(parse(text = paste(
          'r$', argchar, '$test.statistic', sep = ''
        ))) - sum(gofstat[1:2]) + pen(tauhat[c(-1, -2, -length(tauhat))], n) - oldpen
      }
      for (j in scanrange[!scanrange %in% c(1, (length(tauhat) - 2))]) {
        rl = g.tests(
          mstree(as.dist(ydist[(tauhat[j - 1] + 1):tauhat[j + 2], (tauhat[j - 1] +
                                                                     1):tauhat[j + 2]]), min(floor(
                                                                       sqrt(tauhat[j + 2] - tauhat[j - 1])
                                                                     ), g.max)),
          1:(tauhat[j] - tauhat[j - 1]),
          (tauhat[j] - tauhat[j - 1] + 1):(tauhat[j + 2] - tauhat[j - 1]),
          stat
        )
        rr = g.tests(
          mstree(as.dist(ydist[(tauhat[j] + 1):tauhat[j + 3], (tauhat[j] + 1):tauhat[j +
                                                                                       3]]), min(floor(
                                                                                         sqrt(tauhat[j + 3] - tauhat[j])
                                                                                       ), g.max)),
          1:(tauhat[j + 2] - tauhat[j]),
          (tauhat[j + 2] - tauhat[j] + 1):(tauhat[j + 3] - tauhat[j]),
          stat
        )
        change[j] = eval(parse(text = paste(
          'rl$', argchar, '$test.statistic', sep = ''
        ))) + eval(parse(text = paste(
          'rr$', argchar, '$test.statistic', sep = ''
        ))) - sum(gofstat[(j - 1):(j + 1)]) + pen(tauhat[c(-1, -length(tauhat), -(j +
                                                                                    1))], n) - oldpen
      }
      if (length(tauhat) - 2 %in% scanrange) {
        j = length(tauhat) - 2
        r = g.tests(
          mstree(as.dist(ydist[(tauhat[j - 1] + 1):n, (tauhat[j - 1] + 1):n]), min(floor(
            sqrt(n - tauhat[j - 1])
          ), g.max)),
          1:(tauhat[j] - tauhat[j - 1]),
          (tauhat[j] - tauhat[j - 1] + 1):(n - tauhat[j - 1])
        )
        change[j] = eval(parse(text = paste(
          'r$', argchar, '$test.statistic', sep = ''
        ))) - sum(gofstat[(j - 1):j]) + pen(tauhat[c(-1, -length(tauhat), -length(tauhat) +
                                                       1)], n) - oldpen
      }
      rm(j)
      maxind = which.max(change)
      mergetau[i] = tauhat[maxind + 1]
      if (maxind == 1) {
        gofstat = gofstat[-1]
        tauhat = tauhat[-2]
        r = g.tests(mstree(as.dist(ydist[1:tauhat[3], 1:tauhat[3]]), min(floor(
          sqrt(tauhat[3])
        ), g.max)),
        1:tauhat[2],
        (tauhat[2] + 1):tauhat[3],
        stat)
        gofstat[1] = eval(parse(text = paste(
          'r$', argchar, '$test.statistic', sep = ''
        )))
        scanrange = 1:2
      } else if (maxind == length(change)) {
        gofstat = gofstat[-length(change)]
        tauhat = tauhat[-length(change) - 1]
        r = g.tests(
          mstree(as.dist(ydist[(tauhat[maxind - 1] + 1):n, (tauhat[maxind - 1] + 1):n]), min(floor(
            sqrt(n - tauhat[maxind - 1])
          ), g.max)),
          1:(tauhat[maxind] - tauhat[maxind - 1]),
          (tauhat[maxind] - tauhat[maxind - 1] + 1):(n - tauhat[maxind - 1]),
          stat
        )
        gofstat[maxind - 1] = eval(parse(text = paste(
          'r$', argchar, '$test.statistic', sep = ''
        )))
        scanrange = (maxind - 2):(maxind - 1)
      } else{
        gofstat = gofstat[-maxind]
        tauhat = tauhat[-maxind - 1]
        
        rl = g.tests(
          mstree(as.dist(ydist[(tauhat[maxind - 1] + 1):tauhat[maxind + 1], (tauhat[maxind -
                                                                                      1] + 1):tauhat[maxind + 1]]), min(floor(
                                                                                        sqrt(tauhat[maxind + 1] - tauhat[maxind - 1])
                                                                                      ), g.max)),
          1:(tauhat[maxind] - tauhat[maxind - 1]),
          (tauhat[maxind] - tauhat[maxind - 1] + 1):(tauhat[maxind + 1] - tauhat[maxind -
                                                                                   1]),
          stat
        )
        rr = g.tests(
          mstree(as.dist(ydist[(tauhat[maxind] + 1):tauhat[maxind + 2], (tauhat[maxind] +
                                                                           1):tauhat[maxind + 2]]), min(floor(
                                                                             sqrt(tauhat[maxind + 2] - tauhat[maxind])
                                                                           ), g.max)),
          1:(tauhat[maxind + 1] - tauhat[maxind]),
          (tauhat[maxind + 1] - tauhat[maxind] + 1):(tauhat[maxind + 2] - tauhat[maxind]),
          stat
        )
        
        gofstat[maxind - 1] = eval(parse(text = paste(
          'rl$', argchar, '$test.statistic', sep = ''
        )))
        gofstat[maxind] = eval(parse(text = paste(
          'rr$', argchar, '$test.statistic', sep = ''
        )))
        if (maxind == 2) {
          scanrange = 1:3
        } else if (maxind == (length(change) - 1)) {
          scanrange = (maxind - 2):(maxind)
        } else{
          scanrange = (maxind - 2):(maxind + 1)
        }
      }
      change = change[-maxind]
      gofoverall[i + 1] = sum(gofstat) + pen(tauhat[c(-1, -length(tauhat))], n)
    }
  }
  
  if (i < max.iter) {
    oldpen = pen(tauhat[c(-1, -length(tauhat))], n)
    i = m - 1
    change = numeric(0)
    
    r = g.tests(mstree(as.dist(ydist[1:tauhat[4], 1:tauhat[4]]), min(floor(sqrt(
      tauhat[4]
    )), g.max)),
    1:tauhat[3],
    (tauhat[3] + 1):tauhat[4],
    stat)
    
    change[1] = eval(parse(text = paste(
      'r$', argchar, '$test.statistic', sep = ''
    ))) - sum(gofstat[1:2]) + pen(tauhat[c(-1, -2, -length(tauhat))], n) - oldpen
    j = length(tauhat) - 2
    
    r = g.tests(
      mstree(as.dist(ydist[(tauhat[j - 1] + 1):n, (tauhat[j - 1] + 1):n]), min(floor(
        sqrt(n - tauhat[j - 1])
      ), g.max)),
      1:(tauhat[j] - tauhat[j - 1]),
      (tauhat[j] - tauhat[j - 1] + 1):(n - tauhat[j - 1]),
      stat
    )
    
    change[j] = eval(parse(text = paste(
      'r$', argchar, '$test.statistic', sep = ''
    ))) - sum(gofstat[(j - 1):j]) + pen(tauhat[c(-1, -3, -length(tauhat))], n) -
      oldpen
    rm(j)
    maxind = which.max(change)
    mergetau[i] = tauhat[maxind + 1]
    if (maxind == 1) {
      gofstat = gofstat[-1]
      tauhat = tauhat[-2]
      
      r = g.tests(mstree(as.dist(ydist[1:tauhat[3], 1:tauhat[3]]), min(floor(
        sqrt(tauhat[3])
      ), g.max)),
      1:tauhat[2],
      (tauhat[2] + 1):tauhat[3],
      stat)
      
      gofstat[1] = eval(parse(text = paste(
        'r$', argchar, '$test.statistic', sep = ''
      )))
    } else{
      gofstat = gofstat[-length(change)]
      tauhat = tauhat[-length(change) - 1]
      
      r = g.tests(
        mstree(as.dist(ydist[(tauhat[maxind - 1] + 1):n, (tauhat[maxind - 1] + 1):n]), min(floor(
          sqrt(n - tauhat[maxind - 1])
        ), g.max)),
        1:(tauhat[maxind] - tauhat[maxind - 1]),
        (tauhat[maxind] - tauhat[maxind - 1] + 1):(n - tauhat[maxind - 1]),
        stat
      )
      
      gofstat[maxind - 1] = eval(parse(text = paste(
        'r$', argchar, '$test.statistic', sep = ''
      )))
    }
    gofoverall[i + 1] = sum(gofstat) + pen(tauhat[c(-1, -length(tauhat))], n)
    ind = which.max(gofoverall)
    
    tauhat = tauhat[2]
    if (ind != length(gofoverall)) {
      tauhat = c(tauhat, mergetau[(m - 1):ind])
    }
    if (detail == TRUE) {
      rm(ydist, pos = '.GlobalEnv')
      return(list(tauhat, mergetau, gofoverall))
    } else{
      rm(ydist, pos = '.GlobalEnv')
      return(tauhat)
    }
  } else{
    tauhat = tauhat[-c(1, length(tauhat))]
    ind = which.max(gofoverall)
    if (ind != length(gofoverall)) {
      tauhat = c(tauhat, mergetau[i:ind])
    }
    if (detail == TRUE) {
      rm(ydist, pos = '.GlobalEnv')
      return(list(tauhat, mergetau, gofoverall))
    } else{
      rm(ydist, pos = '.GlobalEnv')
      return(tauhat)
    }
  }
  
}
