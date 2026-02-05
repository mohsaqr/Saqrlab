#' Global Names Dataset
#'
#' @description
#' A collection of 300 short, diverse names from around the world for use in
#' social network simulations.
#'
#' @format A character vector of 300 names from various cultures and regions.
#'
#' @details
#' Names are selected to be short (typically 3-7 characters) and represent
#' diverse cultural backgrounds including:
#' \itemize{
#'   \item European (Western, Eastern, Nordic)
#'   \item Arab / Middle Eastern
#'   \item Asian (East Asian, South Asian, Southeast Asian)
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

# European - Western (30)
"Emma", "Liam", "Mia", "Noah", "Lea", "Max", "Zoe", "Leo", "Eva", "Tom",
"Ella", "Jack", "Lily", "Hugo", "Sara", "Finn", "Ruby", "Adam", "Lucy", "Jake",
"Chloe", "Luke", "Grace", "Ethan", "Clara", "Oscar", "Julia", "Felix", "Marie", "Paul",

# European - Eastern (25)
"Ivan", "Yuri", "Nina", "Igor", "Mila", "Sasha", "Daria", "Pavel", "Anya", "Katya",
"Boris", "Tanya", "Dmitri", "Vera", "Alexei", "Irina", "Oleg", "Lena", "Vlad", "Nadia",
"Petra", "Milan", "Jiri", "Hana", "Marek",

# Nordic (25)
"Erik", "Freya", "Lars", "Astrid", "Sven", "Ingrid", "Odin", "Sigrid", "Bjorn", "Liv",
"Olaf", "Elsa", "Axel", "Ida", "Nils", "Saga", "Leif", "Maja", "Knut", "Elin",
"Tor", "Hedda", "Arne", "Greta", "Rune",

# Arab (40)
"Ali", "Layla", "Omar", "Noor", "Amir", "Zara", "Karim", "Leila", "Tariq", "Mariam",
"Yusuf", "Fatima", "Khalid", "Amira", "Samir", "Rania", "Hassan", "Dina", "Walid", "Salma",
"Faris", "Yasmin", "Rashid", "Mona", "Nabil", "Lina", "Jamal", "Huda", "Ziad", "Rana",
"Sami", "Reem", "Adel", "Noura", "Mazen", "Asma", "Faisal", "Dalal", "Bassam", "Aisha",

# Persian / Turkish (20)
"Cyrus", "Shirin", "Dara", "Parisa", "Reza", "Leyla", "Arman", "Soraya", "Babak", "Mina",
"Emre", "Elif", "Kaan", "Defne", "Cem", "Zeynep", "Deniz", "Aylin", "Burak", "Selin",

# East Asian (30)
"Yuki", "Hiro", "Mei", "Jin", "Suki", "Kai", "Lin", "Wei", "Rin", "Chen",
"Akira", "Sakura", "Kenji", "Yuna", "Taro", "Haruki", "Koji", "Miki", "Ryu", "Emi",
"Jing", "Xiao", "Ming", "Hua", "Feng", "Yan", "Jun", "Lei", "Bo", "Lan",

# Korean (10)
"Minho", "Jisoo", "Joon", "Seo", "Tae", "Sora", "Woo", "Minji", "Dae", "Eunji",

# South Asian (25)
"Arun", "Priya", "Raj", "Devi", "Amit", "Sita", "Ravi", "Maya", "Ajay", "Lata",
"Vijay", "Anita", "Sanjay", "Neha", "Rohit", "Pooja", "Arjun", "Kavita", "Rahul", "Sunita",
"Kiran", "Aditi", "Nikhil", "Shreya", "Vikram",

# Southeast Asian (25)
"Linh", "Minh", "Anh", "Bao", "Mai", "Duc", "Thao", "Tuan", "Hoa", "Nam",
"Ting", "Somchai", "Ploy", "Arya", "Dewi", "Rizal", "Putri", "Bagus", "Sari", "Wayan",
"Ayu", "Gede", "Nisa", "Fajar", "Indah",

# African (30)
"Kofi", "Amara", "Kwame", "Zuri", "Juma", "Nia", "Tendai", "Amina", "Sekou", "Fatou",
"Chidi", "Adaora", "Emeka", "Ngozi", "Obi", "Zainab", "Ayana", "Malik", "Imani", "Jelani",
"Asha", "Omari", "Khadija", "Bakari", "Eshe", "Jabari", "Sanaa", "Dayo", "Amani", "Kwesi",

# Latin American (25)
"Juan", "Rosa", "Diego", "Luz", "Pablo", "Sol", "Luis", "Ana", "Carlos", "Lucia",
"Miguel", "Elena", "Andres", "Maria", "Jose", "Sofia", "Pedro", "Camila", "Mateo", "Valentina",
"Rafael", "Isabel", "Sergio", "Carmen", "Alvaro",

# Oceanian / Pacific (15)
"Koa", "Manu", "Tane", "Moana", "Ariki", "Kaia", "Tui", "Hemi", "Aroha", "Wiremu",
"Leilani", "Keanu", "Mahina", "Malia", "Sina"
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
