single_predict <- function(rit,instance){
  class_names <- names(rit[["Class_priors"]])
  nb_class <- length(class_names)
  
  probs <- vector(mode="numeric",length=nb_class) + log10(rit[["Class_priors"]])
  names(probs) <- class_names
  
  model <- rit[["Model"]]
  isCat <- rit[["cat_attr"]]
  
  for(elem in model){
    coeff <- compute_sim(elem[["interaction"]],instance,isCat)
    probs <- probs + coeff*log10(elem[["prevalence"]])
    probs <- probs + {1-coeff}*{log10(1-elem[["prevalence"]])}
  }
  
  # Compute argmax
  response <- names(probs)[which.max(probs)]
  response
}

rit_predict <- function(rit,testset){
  # Apply min-max scaling
  sc <- rit[["scaling"]]
  mins <- sc[["Mins"]]
  spans <- sc[["Spans"]]
  
  sc_data <- mapply(function(attr,idx){
      if(is.ordered(attr))
        {as.vector(attr,mode="numeric") - mins[idx]} / spans[idx]
      else if(is.factor(attr))
        attr
      else
        {attr - mins[idx]} / spans[idx]
  },testset,seq_along(testset),SIMPLIFY=FALSE)
  sc_data <- as.data.frame(sc_data)
  
  # Predict
  response <- vector(mode="character",length=nrow(testset))
  for(i in 1:nrow(testset)){
    response[i] <- single_predict(rit,sc_data[i,])
  }
  response
}