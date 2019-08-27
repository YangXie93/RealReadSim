

slidingWindowMaker <- function(vec){
  res = c()
  size = 1
  n = 1
  
  while(size <=  length(vec)/10000 ){
    size = size *10
  }
  if(size > 1){  
    for(i in seq(1,length(vec),by = size/2)){
      if(i + size -1 > length(vec)){
        x = length(vec)
      }
      else{
        x = (i +size -1)
      }
      res[n] = mean(vec[i:x ])
      n = n +1
    }
  }
  
  return(res)
}