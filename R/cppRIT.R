# Apply min-max scaling over the dataframe
# Ordered factors are reduced to numeric vectors and scaled
# Unordered factors are untouched
mm_scaling <- function(data){
  
  # Min values of features
  mins <- mapply(function(elem){
    if(is.ordered(elem)){
      1
    }
    else if(is.factor(elem)){
      NA
    }
    else{
      min(elem)
    }
  },data)
  
  # Spans of features
  span <- mapply(function(elem,idx){
    if(is.ordered(elem)){
      nlevels(elem) - 1
    }
    else if(is.factor(elem)){
      NA
    }
    else{
      max(elem) - mins[idx]
    }
  },data,seq_along(data))
  
  # Apply min-max scaling
  sc_data <- mapply(function(elem,idx){
    if(is.ordered(elem)){
      {as.vector(elem,mode="numeric") - mins[idx]} / span[idx]
    }
    else if(is.factor(elem)){
      elem
    }
    else{
      {elem - mins[idx]} / span[idx]
    }
  },data,seq_along(data),SIMPLIFY=FALSE)
  
  list(Mins=mins,Spans=span,data=as.data.frame(sc_data))
}

stopParameters <- function(data,classes,epsilon,epsilon_cat,n_trees,depth,min_inter_sz,branch){
  # check parameters
  if(!is.data.frame(data))
    stop("Data parameter must be a dataframe.")
  if(length(classes) != nrow(data))
    stop("Data and response vector must have the same number of rows.")
  if(nrow(data) == 0)
    stop("Data cannot be empty.")
  if(any(is.na(classes)))
    stop("Response vector cannot contain missing values.")
  if(epsilon < 0)
    stop("epsilon hp should be positive.")
  if(epsilon_cat < 0)
    stop("epsilon_cat hp should be positive.")
  if(depth < 0)
    stop("depth should be positive.")
  if(n_trees < 1)
    stop("n_trees should be at least 1.")
  if(min_inter_sz < 2)
    stop("min_inter_sz should be at least 2.")
  if(branch < 1)
    stop("branch should be at least 1.")
}

relaxed_RIT <- function(data,classes,theta,epsilon,epsilon_cat,n_trees=100L,depth=10L,min_inter_sz=2L,branch=5,es=TRUE){
  
  stopParameters(data,classes,epsilon,epsilon_cat,n_trees,depth,min_inter_sz,branch)
  
  # Class information
  classes <- factor(classes)
  class_names <- levels(classes)
  nb_class <- nlevels(classes)
  
  # Check parameter theta
  if(missing(theta)){
    theta <- vector(mode="numeric",length=nb_class) - 1
    names(theta) <- class_names
  }
  else{
    if(length(theta) == 1){ # theta is scalar
      theta <- vector(mode="numeric",length=nb_class) + theta
      names(theta) <- class_names
    }
    else if(!(all(class_names %in% names(theta)))){
      stop("Named vector theta doesn't contain prevalence thresholds for some classes.")
    }
  }
  
  # Apply min-max scaling to data
  mm <- mm_scaling(data)
  data <- mm[["data"]]
  mm[["data"]] <- NULL
  
  isCat <- sapply(data,is.factor) # Logical vector: if each feature/attr is categorical or not
  
  # Split dataset per class and order theta vector
  datas <- vector("list",length=nb_class)
  names(datas) <- class_names
  o_theta <- vector(mode="numeric",length=nb_class)
  names(o_theta) <- class_names
  for(cls in class_names){
    datas[[cls]] <- data[which(classes==cls),]
    o_theta[[cls]] <- theta[[cls]]
  }
  
  # Call main algo
  all_leaves <- cpp_Relaxed_RIT(datas, o_theta, isCat, epsilon, epsilon_cat, n_trees, depth, branch, min_inter_sz,es)
  
  # Filter interactions
  all_leaves <- Filter(function(inter){
    all_low <- all(inter[["prevalence"]] < o_theta)
    all_high <- all(inter[["prevalence"]] > o_theta)
    !all_low && !all_high 
  },all_leaves)
  
  # Compute class prior probabilities
  priors <- data.table(classes)[,.N,keyby=classes]
  priors[,"N"] <- priors[,N] / length(classes)
  priors_vec <- vector(mode="numeric",length=nb_class)
  names(priors_vec) <- class_names
  for(i in 1:nrow(priors)){
    priors_vec[[priors[i,classes]]] <- priors[i,N]
  }
  
  list(Model=all_leaves,Class_priors=priors_vec,cat_attr=isCat,scaling=mm)
}
