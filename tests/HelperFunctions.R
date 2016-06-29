setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sourceCpp("cppGetMinPerCol.cpp")
GetMinPerCol(matrix(c(1,4,5,2), ncol = 2))
