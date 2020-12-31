disProBar <- function(k, NN){
  # k is current loop number
  # NN is total loop number
  if(k%%(NN/10)==0){
    print(paste('==============running process: ',10*floor(k/(NN/10)),'% =============', sep = ''))
  }
}
