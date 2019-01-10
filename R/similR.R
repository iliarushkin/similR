

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
#' @param similarity_method as implemented in quanteda::textstat_simil()
#' @param keep_vec  whether to return the matrix of the word-vectors
#'
#' @return depending on keep_vec, either just the simil matrix/object, or a list with that and the matrix of the word-vectors
#' @export
#'
#' @examples
similR=function(toks, vec=NULL, window_weights=1/(1:5), word_vectors_size=30, x_max=10, n_iter=30, ik=100,
                clustering_algorithm=c("MacQueen", "Hartigan-Wong", "Lloyd", "Forgy"),
                clustering_itermax=1000,
                similarity_method='cosine',
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
  simmat=quanteda::textstat_simil(dfm,method=similarity_method)
  cat('Done in', round(proc.time()[3]-tic),'sec.\n')
  if(keep_vec) return(list(simmat=simmat,vec=vec)) else return(simmat)

}
