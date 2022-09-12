#Retrieved from https://rdrr.io/github/gadenbuie/regexhelp/src/R/shiny_modified_inputs.R

# Modified Text Input
#
# Standard [shiny::textInput()] with additional `width` parameter, added code
# font style for the input text and with `autocomplete`, `autocorrect`,
# `autocapitalize` and `spellcheck` set to `off` or `false`.
#
# @inheritParams shiny::textInput
# @param width Width of `shiny-input-container` div.
# @param ... Extra elements to be included in the `input-group` div.
# @family modified shiny inputs

textInput_capital <- function(
  inputId,
  label,
  value = "",
  width = NULL,
  placeholder = NULL,
  ...
) {
  `%AND%` <- getFromNamespace("%AND%", "shiny")
  value <- shiny::restoreInput(id = inputId, default = value)
  
  shiny::div(
    class = "input-group shiny-input-container",
    style = if (!is.null(width)) paste0("width: ", shiny::validateCssUnit(width), ";"),
    label %AND% shiny::tags$label(label, `for` = inputId),
    shiny::tags$input(
      id = inputId,
      type = "text",
      class = "form-control",
      value = value,
      style = 'font-family: "Monaco", "Inconsolata", monospace;',
      autocomplete = "off",
      autocorrect = "off",
      autocapitalize = "on",
      spellcheck = "false",
      placeholder = placeholder
    ),
    ...
  )
}