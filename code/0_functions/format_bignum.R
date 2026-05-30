# function to format the cell number in nice, descriptive manner
format_bignum = function(n){
  case_when(
    n >= 1e9 ~ paste(round(n/1e9, 1), 'billion cells'),
    n >= 1e6 ~ paste(round(n/1e6, 1), 'million cells'),
    n >= 1e3 ~ paste(format(round(n), big.mark = ",", scientific = FALSE), "cells"),
    TRUE ~ as.character(n))
}
