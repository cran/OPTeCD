#' Optimal Partial Tetra-Allele Cross Designs for Prime Number of Lines
#'
#' @param v Any prime number (>3)
#'
#' @return This function generates a series of universally optimal family of designs using MOLS. The method starts with selecting any of the (N-1)/2  MOLS of a given order N (the number of lines) and retaining the first four rows and making crosses with the lines occurring in each column. The parameters of the developed class of design is total number of crosses (T) =N*(N-1)/2, number of blocks (b) = (N-1)/2, number of replications (r), block sizes (k) = N and  degree of fractionation (f)=4/(N-2)(N-3)
#' @export
#'
#' @examples
#' library(OPTeCD)
#' OPTeCD(5)
OPTeCD<-function(v){
  MOLS<-function(v){
    MOLS<-list()
    for(i in 1:(v-1)){
      mols<-NULL
      seq<-c(1:v)
      j=0
      repeat{
        mols<-rbind(mols,c(seq+j))
        j=j+i
        if(nrow(mols)==v){
          mols=mols%%v
          mols[mols==0]<-v
          MOLS<-append(MOLS,list(mols))
          break
        }
      }
    }
    return(MOLS)
  }
  ############obtain half of the MOLS
  mols=MOLS(v)[1:((v-1)/2)]
  ###Take four rows from the to and then transpose each
  mols1<-lapply(mols,function(mat)t(mat[1:4,]))
  ###Show like original
  Tet_format<-function(mat){
    tet<-cbind("(",mat[,1],"X",mat[,2],")","X","(",mat[,3],"X",mat[,4],")")
    prmatrix(tet,collab = rep("",(ncol(tet))),quote = FALSE)
  }
  ###########Naming for each block

  ######### Transform to the final design
  #######
  names(mols1)<-paste("     Block",seq_along(mols1))
  ##########
  # View the names and the corresponding matrix display
  message("PTac Design")
  cat("\n")
  for (i in 1:length(mols1)) {
    message(names(mols1)[i])
    Tet_format(mols1[[i]])
    cat("\n")
  }
  #######Parameters of the design
  parameter_list<-list("N"=v,"T"=v*(v-1)*0.5,"b"=(v-1)*0.5,"k"=v,"f"=4/((v-2)*(v-3)))
  return(parameter_list)
}

