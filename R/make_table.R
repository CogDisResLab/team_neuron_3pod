make_table <- function(df, caption = NULL) {
  df <- df %>%
    dplyr::mutate(dplyr::across(
      .cols = dplyr::where(is.numeric),
      .fns = ~ round(.x, 3)
    ))
  
  DT::datatable(
    df,
    rownames = FALSE,
    caption = htmltools::tags$caption(style = 'text-align: left;', caption),
    options = list(
      scrollX = TRUE,
      scrollY = TRUE,
      paging = TRUE,
      fixedHeader = TRUE,
      pageLength = 10
    )
  )
}