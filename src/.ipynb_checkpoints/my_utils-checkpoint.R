make_table_prop <- function(dd, xvar, fillvar) {
  tt = 100*table(dd[,xvar],dd[,fillvar]) / apply(table(dd[,xvar],dd[,fillvar]),1,sum)
  dt = data.frame(tt)
  colnames(dt) = c(xvar,fillvar,"prop")
  return(dt)
}


make_factor <- function(vec_name, vec_value, myfun=median) {
  tt = sapply(unique(vec_name), function(x) myfun(vec_value[vec_name==x]))
  tt = sort(tt)
  ddt = data.frame(num=tt)
  ddt$name=rownames(ddt)
  ddt$label = paste(ddt$name,"\n", paste0("median=",round(ddt$num,2)))
  return(ddt)
}