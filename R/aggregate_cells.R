#' Convert array of quosure (e.g. c(col_a, col_b)) into character vector
#' 
#' @keywords internal
#'
#' @importFrom rlang quo_name
#' @importFrom rlang quo_squash
#' @importFrom purrr when
#' @importFrom magrittr equals
#'
#' @param v A array of quosures (e.g. c(col_a, col_b))
#'
#' @return A character vector
quo_names <- function(v) {
  
  v = rlang::quo_name(rlang::quo_squash(v))
  gsub('^c\\(|`|\\)$', '', v) %>% 
    strsplit(', ') %>% 
    unlist 
}

#' Remove class to abject
#'
#'
#' @param var A tibble
#' @param name A character name of the class
#'
#' @return A tibble with an additional attribute
drop_class = function(var, name) {
  class(var) <- class(var)[!class(var)%in%name]
  var
}

get_specific_annotation_columns = function(.data, .col){
  
  
  # Comply with CRAN NOTES
  . = NULL
  
  # Make col names
  .col = enquo(.col)
  
  # x-annotation df
  n_x = .data %>% dplyr::distinct_at(vars(!!.col)) %>% nrow
  
  # element wise columns
  .data %>%
    select(-!!.col) %>%
    colnames %>%
    map(
      ~
        .x %>%
        when(
          .data %>%
            distinct_at(vars(!!.col, .x)) %>%
            nrow %>%
            magrittr::equals(n_x) ~ (.),
          ~ NULL
        )
    ) %>%
    
    # Drop NULL
    {	(.)[lengths((.)) != 0]	} %>%
    unlist
  
}


subset = 		function(.data,	 .column)	{
  # Make col names
  .column = enquo(.column)
  
  # Check if column present
  if(quo_names(.column) %in% colnames(.data) %>% all %>% `!`)
    stop("nanny says: some of the .column specified do not exist in the input data frame.")
  
  .data %>%
    
    # Selecting the right columns
    select(	!!.column,	get_specific_annotation_columns(.data, !!.column)	) %>%
    distinct()
  
}

#' @export
aggregate_cells = function(.data, .sample) {
  
  .sample = enquo(.sample)
  
  .data %>%
    
    tidySingleCellExperiment::nest(data = -!!.sample) %>%
    mutate(data = map(data, ~ 
                        
                        # loop over assays
                        map2(
                          .x %>% assays %>% as.list() ,
                          .x %>% assays %>% names(),
                          
                          # Get counts
                          ~ .x %>%
                            Matrix::rowSums(na.rm = T) %>%
                            tibble::enframe(
                              name  = "transcript",
                              value = sprintf("abundance_%s", .y)
                            )
                        ) %>%
                        Reduce(function(...) full_join(..., by=c("transcript")), .)
                      
    )) %>%
    left_join(.data %>% tidySingleCellExperiment::as_tibble() %>% subset(!!.sample)) %>%
    tidySingleCellExperiment::unnest(data) %>%
    
    drop_class("tidySingleCellExperiment_nested") %>%
    
    tidybulk::as_SummarizedExperiment(
      .sample = !!.sample,
      .transcript = transcript,
      .abundance = !!as.symbol(sprintf("abundance_%s", names(assays(.data))[1]))
    ) 
  
}
