#' Global Names Dataset
#'
#' @description
#' A collection of 100 short, diverse names from around the world for use in
#' social network simulations.
#'
#' @format A character vector of 100 names from various cultures and regions.
#'
#' @details
#' Names are selected to be short (typically 3-7 characters) and represent
#' diverse cultural backgrounds including:
#' \itemize{
#'   \item European (Western, Eastern, Nordic)
#'   \item Asian (East Asian, South Asian, Southeast Asian, Middle Eastern)
#'   \item African
#'   \item Latin American
#'   \item Oceanian
#' }
#'
#' @examples
#' # View all names
#' GLOBAL_NAMES
#'
#' # Sample 10 random names
#' sample(GLOBAL_NAMES, 10)
#'
#' # Use in network simulation
#' simulate_edge_list(n_nodes = 20)
#'
#' @export
GLOBAL_NAMES <- c(

# European - Western
"Emma", "Liam", "Mia", "Noah", "Lea", "Max", "Zoe", "Leo", "Eva", "Tom",
# European - Eastern
"Ivan", "Olga", "Yuri", "Nina", "Igor", "Mila", "Sasha", "Daria", "Pavel", "Anya",
# Nordic
"Erik", "Freya", "Lars", "Astrid", "Sven", "Ingrid", "Odin", "Sigrid", "Bjorn", "Liv",
# East Asian
"Yuki", "Hiro", "Mei", "Jin", "Suki", "Kai", "Lin", "Wei", "Rin", "Chen",
# South Asian
"Arun", "Priya", "Raj", "Devi", "Amit", "Sita", "Ravi", "Maya", "Ajay", "Lata",
# Southeast Asian
"Linh", "Minh", "Anh", "Bao", "Mai", "Duc", "Lan", "Tuan", "Hoa", "Nam",
# Middle Eastern
"Ali", "Layla", "Omar", "Noor", "Amir", "Zara", "Karim", "Leila", "Tariq", "Hana",
# African
"Kofi", "Amara", "Kwame", "Zuri", "Juma", "Nia", "Tendai", "Amina", "Sekou", "Fatou",
# Latin American
"Juan", "Rosa", "Diego", "Luz", "Pablo", "Sol", "Luis", "Ana", "Carlos", "Lucia",
# Oceanian
"Koa", "Manu", "Tane", "Moana", "Ariki", "Kaia", "Tui", "Hemi", "Aroha", "Wiremu"
)

#' Get Global Names
#'
#' @description
#' Retrieve names from the global names dataset, optionally sampling a subset.
#'
#' @param n Integer or NULL. Number of names to return. If NULL, returns all names.
#'   Default: NULL.
#' @param seed Integer or NULL. Random seed for reproducibility when sampling.
#'   Default: NULL.
#'
#' @return Character vector of names.
#'
#' @examples
#' # Get all names
#' get_global_names()
#'
#' # Get 15 random names
#' get_global_names(n = 15, seed = 42)
#'
#' @export
get_global_names <- function(n = NULL, seed = NULL) {
if (!is.null(seed)) set.seed(seed)

if (is.null(n)) {
    return(GLOBAL_NAMES)
}

if (n > length(GLOBAL_NAMES)) {
    warning(sprintf(
    "Requested %d names but only %d available. Returning all names.",
    n, length(GLOBAL_NAMES)
    ))
    return(GLOBAL_NAMES)
}

sample(GLOBAL_NAMES, n)
}
