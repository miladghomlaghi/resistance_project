transform_links <-function (links_options,user_input){
 
  binary_links <- rep(0,length(user_input))
  
  for (b in 1:length(user_input)){
    
  if (stringi::stri_cmp_eq(user_input[b],links_options[1])){
    binary_links[b] <- 1
  } else if (stringi::stri_cmp_eq(user_input[b],links_options[2])){
    binary_links[b] <- -1
  } else{
    binary_links[b] <- 0
  }
  }
return(binary_links)
}