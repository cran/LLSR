###############################################################################
# Shared helpers for AQSearch* functions
###############################################################################

build_cas_search_db <- function(cas_db) {
  cas_db$CAS.NAME_lower <- tolower(cas_db$CAS.NAME)
  cas_db$CAS.COMMON_lower <- tolower(cas_db$CAS.COMMON)
  cas_db$CAS.CODE_lower <- tolower(cas_db$CAS.CODE)
  cas_db
}

search_component_names <- function(component_query, cas_search_db) {
  if (!is.valid(component_query)) return(character(0))
  query_lower <- tolower(component_query)
  match_idx <- grepl(query_lower, cas_search_db$CAS.NAME_lower, fixed = TRUE) |
    grepl(query_lower, cas_search_db$CAS.COMMON_lower, fixed = TRUE) |
    grepl(query_lower, cas_search_db$CAS.CODE_lower, fixed = TRUE)
  unique(cas_search_db$CAS.NAME[match_idx])
}

###############################################################################
