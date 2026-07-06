# Get Global Names

Retrieve names from the global names dataset, with optional filtering by
region.

## Usage

``` r
get_global_names(n = NULL, regions = "all", seed = NULL)
```

## Arguments

- n:

  Integer or NULL. Number of names to return. If NULL, returns all
  matching names. Default: NULL.

- regions:

  Character vector. Regions to sample from. Can be:

  - `"all"`: All regions (default)

  - Specific regions: e.g., `"arab"`, `"east_asia"`, `"nordic"`

  - Shortcuts: e.g., `"europe"`, `"africa"`, `"asia"`

  - Multiple: e.g., `c("arab", "persian", "turkish")`

  Default: "all".

- seed:

  Integer or NULL. Random seed for reproducibility. Default: NULL.

## Value

Character vector of names.

## See also

[`list_name_regions`](https://pak.dynasite.org/Saqrlab/reference/list_name_regions.md)
for available regions,
[`GLOBAL_NAMES`](https://pak.dynasite.org/Saqrlab/reference/GLOBAL_NAMES.md)
for the full dataset.

## Examples

``` r
# Get 20 random names from all regions
get_global_names(20, seed = 42)
#>  [1] "Shu"        "Soraya"     "Rangi"      "Cuauhtemoc" "Gagik"     
#>  [6] "Anahera"    "Oyunbileg"  "Sem"        "Tevita"     "Eloise"    
#> [11] "Firuz"      "Joao"       "Oyun"       "Nursultan"  "Riad"      
#> [16] "Jisoo"      "Esther"     "Fernando"   "Bataar"     "Shoshana"  

# Get 10 names from Arab region
get_global_names(10, regions = "arab")
#>  [1] "Mona"  "Dina"  "Dalia" "Omar"  "Dalal" "Ahmed" "Noha"  "Rania" "Huda" 
#> [10] "Amr"  

# Get 15 names from any African region
get_global_names(15, regions = "africa")
#>  [1] "Palesa"  "Farai"   "Aminata" "Mpho"    "Akua"    "Zuri"    "Solange"
#>  [8] "Tapiwa"  "Fartun"  "Makeda"  "Sihem"   "Ngozi"   "Sylvie"  "Yassine"
#> [15] "Girma"  

# Get names from multiple regions
get_global_names(20, regions = c("east_asia", "south_asia"))
#>  [1] "Tahera" "Hamza"  "Rohit"  "Xiao"   "Surya"  "Varun"  "Arjun"  "Chen"  
#>  [9] "Rahul"  "Uma"    "Lin"    "Arun"   "Ganesh" "Haruki" "Jing"   "Yuki"  
#> [17] "Sabita" "Priya"  "Maya"   "Sung"  

# Get all Nordic names
get_global_names(regions = "nordic")
#>  [1] "Erik"      "Astrid"    "Oscar"     "Maja"      "Axel"      "Ella"     
#>  [7] "Viktor"    "Saga"      "Gustav"    "Wilma"     "William"   "Ebba"     
#> [13] "Linus"     "Agnes"     "Filip"     "Freja"     "Lars"      "Ingrid"   
#> [19] "Sven"      "Liv"       "Olaf"      "Nora"      "Leif"      "Thea"     
#> [25] "Magnus"    "Emma"      "Henrik"    "Sofie"     "Andreas"   "Ida"      
#> [31] "Kristian"  "Julie"     "Mikkel"    "Freya"     "Emil"      "Clara"    
#> [37] "Christian" "Emilie"    "Noah"      "Alma"      "Oliver"    "Lucas"    
#> [43] "Laura"     "Frederik"  "Anna"      "Mathias"   "Onni"      "Aino"     
#> [49] "Eetu"      "Venla"     "Lauri"     "Helmi"     "Eero"      "Siiri"    
#> [55] "Matias"    "Emilia"    "Aleksi"    "Sofia"     "Veeti"     "Aada"     
#> [61] "Elias"     "Ragnar"    "Sigrid"    "Bjorn"     "Gudrun"    "Eirik"    
#> [67] "Hildur"    "Odin"      "Frida"     "Thor"      "Helga"     "Gunnar"   
#> [73] "Unnur"     "Bjarni"    "Kristin"   "Arnar"     "Margret"  
```
