radioImages <- function(inputId, images, values){
  radios <- lapply(
    seq_along(images),
    function(i) {
      id <- paste0(inputId, i)
      tagList(
        tags$input(
          type = "radio",
          name = inputId,
          id = id,
          class = "input-hidden",
          value = as.character(values[i]),
          checked=if(i==1){"checked"}
          
        ),
        tags$label(
          `for` = id,
          tags$img(
            src = images[i]
          )
        )
      )
    }
  )
  do.call(
    function(...) div(..., class = "shiny-input-radiogroup", id = inputId), 
    radios
  )
}