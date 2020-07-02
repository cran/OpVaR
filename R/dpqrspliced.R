dspliced<-
  function(x, sevdist){
    sapply(x,function(x){
      ifelse(x<=sevdist$thresh, sevdist$weights[1]*sevdist$par[[2]][[1]](x),sevdist$weights[2]*sevdist$par[[1]][[1]](x))
    }
    )
  }

pspliced <-
  function(q, sevdist){
    sapply(q,function(q){
      ifelse(q<=sevdist$thresh, sevdist$weights[1]*sevdist$par[[2]][[2]](q), sevdist$weights[1]+sevdist$weights[2]*(sevdist$par[[1]][[2]](q)))
    }
    )
  }

qspliced<-
  function(p, sevdist){
    sapply(p,function(p){
      ifelse(p<=sevdist$weights[1],sevdist$par[[2]][[3]](p/sevdist$weights[1]),sevdist$par[[1]][[3]]((p-sevdist$weights[1])/sevdist$weights[2]))
    }
    )
  }

rspliced<-
  function (n,sevdist){
    qspliced(runif(n),sevdist)
  }