

#' similR
#'
#' Document similarity based on word vectorization
#'
#' @param toks tokens representing a document corpus
#' @param vec if NULL, will calculate the vector embedding of words. If this lengthy calculation needs to be skipped - pass the prepared matrix of vectors here (for instance, taken from the previous run of this function with keep_vec=TRUE).
#' @param window_weights weights of the window to use for co-occurrence of tokens.
#' @param word_vectors_size dimensionality of vector space
#' @param x_max max number of co-occurrences to use in the weighting function
#' @param n_iter number of GloVe iterations
#' @param ik initial number of clusters (the final one will be greater or equal than this)
#' @param clustering_algorithm as implemented in kmeans(), but note that 'Hartigan-Wong' fails for large data, so here we made 'MacQueen' the default.
#' @param clustering_itermax as implemented in kmeans()
#' @param similarity_method as implemented in quanteda::textstat_simil(), plus an extra method 'jsd' for Jensen-Shannon divergence.
#' @param keep_vec  whether to return the matrix of the word-vectors
#'
#' @return depending on keep_vec, either just the simil matrix/object, or a list with that and the matrix of the word-vectors
#' @export
#'
#' @examples
similR=function(toks, vec=NULL, window_weights=1/(1:5), word_vectors_size=30, x_max=10, n_iter=30, ik=100,
                clustering_algorithm=c("MacQueen", "Hartigan-Wong", "Lloyd", "Forgy"),
                clustering_itermax=1000,
                similarity_method=c("jsd", "correlation", "cosine", "jaccard", "ejaccard", "dice", "edice", "simple matching", "hamann", "faith"),
                keep_vec=FALSE

                ){

  tic=proc.time()[3]
  cooc=quanteda::fcm(toks,context='window', window=length(window_weights), weights=window_weights, count='weighted')

  #Vectorize using text2vec
  if(is.null(vec)){
    glove=text2vec::GlobalVectors$new(word_vectors_size = word_vectors_size, vocabulary = featnames(cooc), x_max = x_max)
    vec_main=text2vec::fit_transform(cooc,glove,n_iter=n_iter)
    vec_context=glove$components
    vec=vec_main+t(vec_context)
  }


  hotspots=xmeans(vec, ik=ik, algorithm=clustering_algorithm[1], iter.max = clustering_itermax)$cluster

  dfm=quanteda::dfm(quanteda::tokens_select(toks,featnames(cooc)))
  colnames(dfm)=hotspots
  dfm=quanteda::dfm_compress(dfm)
  dfm=quanteda::dfm_weight(dfm,scheme='prop')
  if(similarity_method[1]=='jsd'){
    simmat=1-dfm%>%as.matrix()%>%proxy::dist(method=jsd)
  }else{
    simmat=quanteda::textstat_simil(dfm,method=similarity_method[1])
  }
  cat('Done in', round(proc.time()[3]-tic),'sec.\n')
  if(keep_vec) return(list(simmat=simmat,vec=vec)) else return(simmat)

}









#' Jensen-Shannon divergence of discrete probability distibutions JSD(p||q)
#'
#' @param p first probability distribution
#' @param q second probability distribution
#'
#' @return JSD(p||q), a numeric value
#' @export
#' @details Assumes that both p and q are normalized (sum(p)=sum(q)=1) and that they contain no NA values. Calculated with base-2 log, so the bounds are [0,1].
jsd=function(p,q){

  # p=p/sum(p)
  # q=q/sum(q)

  m=0.5*(p+q)

  x=log(m/p)
  x[is.infinite(x)]=0
  x[is.na(x)]=0

  y=log(m/q)
  y[is.infinite(y)]=0
  y[is.na(y)]=0

  #Divide by log(2) so to normalize to the interval
  return(-0.5*(sum(p*x)+sum(q*y))/log(2))

}

#' Kullback-Leibler divergence between discrete probability distributions KLD(p||q)
#'
#' @param p first probability distribution ("to which"), a numeric vector, no NAs
#' @param q second probability distribution ("from which"), a numeric vector of same length as \code{p}, no NAs
#'
#' @return KLD(p||q), a numeric value
#' @export
#' @details Assumes that both p and q are normalized (sum(p)=sum(q)=1). Note that by definition the Kullback-Leibler divergence is asymmetric: KLD(p||q) is not equal to KLD(q||p).
kld=function(p,q){

  # p=p/sum(p)
  # q=q/sum(q)

  z=log(q/p)
  z[is.infinite(z)]=0
  z[is.na(z)]=0

  return(-sum(p*z))


}
