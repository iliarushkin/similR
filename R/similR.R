

#' similR
#'
#' Document similarity based on word vectorization
#'
#' @param toks1 tokens from a document corpus
#' @param toks2 tokens from another document corpus (if NULL, assumes toks2=toks1)
#' @param vec if NULL, will calculate the vector embedding of words. If this lengthy calculation needs to be skipped - pass the prepared matrix of vectors here (for instance, taken from the previous run of this function with keep_vec=TRUE).
#' @param window_weights weights of the window to use for co-occurrence of tokens.
#' @param word_vectors_size dimensionality of vector space
#' @param x_max max number of co-occurrences to use in the weighting function
#' @param n_iter number of GloVe iterations
#' @param ik initial number of clusters (if soft_ik is TRUE, the final one will be greater or equal than this)
#' @param soft_ik whether to use xmeans instead of kmeans for hotspots
#' @param clustering_algorithm as implemented in kmeans(), but note that 'Hartigan-Wong' fails for large data, so here we made 'MacQueen' the default.
#' @param clustering_itermax as implemented in kmeans()
#' @param similarity_method as implemented in quanteda::textstat_simil(), plus an extra method 'jsd' for Jensen-Shannon divergence.
#' @param keep_vec  whether to return the matrix of the word-vectors
#'
#' @return depending on keep_vec, either just the simil matrix/object, or a list with that and the matrix of the word-vectors
#' @export
#'
#' @examples
similR=function(toks1, toks2, vec=NULL, hotspots=NULL, window_weights=1/(1:5), word_vectors_size=300, x_max=10, n_iter=30, ik=100,
                soft_ik=FALSE,
                clustering_algorithm=c("MacQueen", "Hartigan-Wong", "Lloyd", "Forgy"),
                clustering_itermax=1000,
                similarity_method=c("cosine", "jaccard", "ejaccard", "dice", "edice", "simple matching", "hamann", "faith", "correlation", "jsd"),
                keep_vec=FALSE,
                keep_hotspots=FALSE,
                ...
                ){

  tic=proc.time()[3]

  if(!is.null(toks2)){
    quanteda::docnames(toks1)=paste0('IDtoks1_', quanteda::docnames(toks1))
    quanteda::docnames(toks2)=paste0('IDtoks2_', quanteda::docnames(toks2))
    toks1=add_tokens_with_empties(toks1,toks2)
  }

  cooc=quanteda::fcm(toks1,context='window', window=length(window_weights), weights=window_weights, count='weighted')

  if(is.null(hotspots)){
    #Vectorize using text2vec
    if(is.null(vec)){
      glove=text2vec::GlobalVectors$new(word_vectors_size = word_vectors_size, vocabulary = featnames(cooc), x_max = x_max, ...)
      vec_main=text2vec::fit_transform(cooc,glove,n_iter=n_iter)
      vec_context=glove$components
      vec=vec_main+t(vec_context)
    }

    cat('Computing hotspots...\n')
    if(soft_ik){
      hotspots=xmeans(vec, ik=ik, algorithm=clustering_algorithm[1], iter.max = clustering_itermax)$cluster
    }else{
      hotspots=kmeans(vec, centers=ik, algorithm=clustering_algorithm[1], iter.max = clustering_itermax)$cluster
    }
  }

  dfm=quanteda::dfm(quanteda::tokens_select(toks1,quanteda::featnames(cooc)))
  quanteda::featnames(dfm)=hotspots
  dfm=quanteda::dfm_compress(dfm)
  dfm=quanteda::dfm_weight(dfm,scheme='prop')

  cat('Similarity method:', similarity_method[1],'\n')
  if(is.null(toks2)){
    simmat=quanteda::textstat_simil(dfm, method=similarity_method[1])
  }else{
    dfm2=quanteda::dfm_subset(dfm, grepl('^IDtoks2_',quanteda::docnames(dfm)))
    dfm=quanteda::dfm_subset(dfm, grepl('^IDtoks1_',quanteda::docnames(dfm)))
    quanteda::docnames(dfm)=gsub('^IDtoks1_','',quanteda::docnames(dfm))
    quanteda::docnames(dfm2)=gsub('^IDtoks2_','',quanteda::docnames(dfm2))
    simmat=quanteda::textstat_proxy(y=dfm, x=dfm2, method=similarity_method[1])%>%t()%>%as.matrix()
  }


  # if(similarity_method[1]=='jsd'){
    # simmat=1-(dfm%>%as.matrix()%>%proxy::dist(method=jsd))/log(2)
  # }
  cat('Done in', round(proc.time()[3]-tic),'sec.\n')
  if(keep_hotspots){
    simmat=list(simmat=simmat, hotspots=hotspots)
  }else{
      if(keep_vec) simmat=list(simmat=simmat, vec=vec)
  }

  return(simmat)

}


add_tokens_with_empties=function(toks1, toks2){

  toks1=toks1%>%as.list()%>%lapply(function(x){x[x=='']=' '; x})%>%quanteda::as.tokens()
  toks2=toks2%>%as.list()%>%lapply(function(x){x[x=='']=' '; x})%>%quanteda::as.tokens()
  attr(toks2,'what')=attr(toks1,'what')
  (toks1+toks2)%>%quanteda::tokens_remove(' ', padding=TRUE)


}






#' Jensen-Shannon divergence of discrete probability distibutions JSD(p||q)
#'
#' @param p first probability distribution
#' @param q second probability distribution
#'
#' @return JSD(p||q), a numeric value
#' @export
#' @details Assumes that both p and q are normalized (sum(p)=sum(q)=1) and that they contain no NA values. Calculated with base-e log, so the bounds are [0,log(2)].
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


  return(-0.5*(sum(p*x)+sum(q*y)))

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
