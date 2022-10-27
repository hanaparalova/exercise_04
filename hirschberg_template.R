#' Align two sequences globally using Hirschberg's algorithm
#' 
#' @param seq1 DNAString object representing NT or AA sequence to align
#' @param seq2 DNAString object representing NT or AA sequence to align
#' @param align A list of DNAString objects with alignment of input sequences
#' @param match An integer value of a score for matching bases
#' @param mismatch An integer value of score for mismatching bases
#' @param gap An integer value of penalty for gap insertion

library(Biostrings)
library(seqinr)

NW<-function(s1,s2,match,mismatch,gap){
  l1<-length(s1[[1]])+1; l2<-length(s2[[1]])+1; mat<-matrix(data=0,nrow=l1,ncol=l2)
  for (i in 2:l2) {mat[1,i]<-(i-1)*gap}
  for (i in 2:l1) {mat[i,1]<-(i-1)*gap}
  for (i in 2:l1) {
    for (j in 2:l2){
      mUp<-mat[i-1,j]+gap; mLeft<-mat[i,j-1]+gap;
      if (s1[[1]][i-1]==s2[[1]][j-1]){
        mDiag<-mat[i-1,j-1]+match
      }else{mDiag<-mat[i-1,j-1]+mismatch}
      win<-which(c(mUp,mLeft,mDiag)==max(mUp,mLeft,mDiag))
      if (length(win)>1){win<-win[length(win)]}
      mat[i,j]<-c(mUp,mLeft,mDiag)[win]
    }
  }
  return(mat[l1,])
}

hirschberg_template <- function(seq1, seq2, align, match, mismatch, gap){
    
    first_align_row <- align[[1]] # initialize the first row of alignment
    second_align_row <- align[[2]] # initialize the second row of alignment
    l1 <- length(seq1[[1]])
    l2 <- length(seq2[[1]])
  
    if (l1==0) # length of seq1 is equal to zero
    {
        for (i in 1:l2)# for each character in seq2
        {
            first_align_row <- c(first_align_row,DNAString('-'))# add gap
            second_align_row <- c(second_align_row,seq2[[1]][i])# add character from seq2
        }
        align <- c(first_align_row, second_align_row)
        print(align)
    }
    else if (l2==0)# length of seq2 is equal to zero
    {
        for (j in 1:l1)# for each character in seq1
        {
            first_align_row <- c(first_align_row, seq1[[1]][j])# add character from seq1
            second_align_row <- c(second_align_row, DNAString('-'))# add gap
        }
        align <- c(first_align_row, second_align_row)
        print(align)
    }
    else if (l1==1 && l2==1)# length of seq1 and seq2 is equal to 1
    {
        first_align_row <- c(first_align_row, seq1[[1]][1])# add character from seq1
        second_align_row <- c(second_align_row,seq2[[1]][1])# add character from seq2
        align <- c(first_align_row, second_align_row)
        print(align)
    }
    else
    {
        x_len <- l1 # length of seq1
        x_mid <- round(l1/2) # half of the length of seq1
        y_len <- l2 # length of seq2 
        
        left_score <- NW(subseq(seq1,start=1,end=x_mid), seq2, match, mismatch, gap)# NW score for the first half of seq1 and the whole seq2
        right_score <- NW(reverse(subseq(seq1,start=x_mid+1,end=l1)), reverse(seq2), match, mismatch, gap)# NW score for the second half of seq1 and the whole seq2 (both are reversed)
        y_mid <- (which((left_score+rev(right_score))==max(left_score+rev(right_score))))-1 # index of division for seq2

        # The first half
        if (y_mid == 0) # index of division for seq2 is equal to 0
        {
            align <- hirschberg_template(subseq(seq1,start=1,end=x_mid), DNAStringSet(), align, match, mismatch, gap)  # call hirschberg function for the first half of seq1 and for an empty DNAString object
        }
        else
        {
            align <- hirschberg_template(subseq(seq1,start=1,end=x_mid), subseq(seq2, start=1, end=y_mid), align, match, mismatch, gap) # call hirschberg function for the first half of seq1 and for the first part of seq2
        }
        
        # The second half
        if ((x_mid + 1) > x_len) # seq1 cannot be further divided
        {
            align <- hirschberg_template(DNAStringSet(), subseq(seq2,start=y_mid,end=l2), align, match, mismatch, gap)# call hirschberg function for an empty DNAString object and the second half of seq2
        }
        else if ((y_mid + 1) > y_len) # seq2 cannot be further divided
        {
            align <- hirschberg_template(subseq(seq1,start=x_mid,end=l1), DNAStringSet(), align, match, mismatch, gap) # call hirschberg function for the second half of seq1 and for an empty DNAString object
        }
        else 
        {
            align <- hirschberg_template(subseq(seq1,start=x_mid,end=l1), subseq(seq2,start=y_mid,end=l2), align, match, mismatch, gap)# call hirschberg function for the second half of seq1 and the second part of seq2
        }
    }
    
    return(align)
}

s1<-DNAStringSet('AGTACGCA')
s2<-DNAStringSet('TATGC')
align<-c(DNAString(''),DNAString(''))
hirschberg_template(s1, s2, align, 2, -1, -2)
