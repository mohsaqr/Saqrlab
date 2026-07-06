# Global Names Dataset

A collection of ~1000 culturally diverse names from around the world,
organized by region for use in social network simulations.

## Usage

``` r
GLOBAL_NAMES
```

## Format

A named list of character vectors, with each element containing names
from a specific region.

## Details

Regions are organized geographically and include:

**Europe:**

- `western_europe`: UK, Ireland, France, Germany, Netherlands, Belgium

- `eastern_europe`: Russia, Ukraine, Poland, Czech, Romania, Hungary,
  Bulgaria

- `southern_europe`: Spain, Italy, Portugal, Greece

- `nordic`: Sweden, Norway, Denmark, Finland, Iceland

**Middle East & Central Asia:**

- `arab`: Gulf states, Levant, Egypt

- `persian`: Iran, Afghanistan, Tajikistan

- `turkish`: Turkey, Azerbaijan

- `central_asia`: Kazakhstan, Uzbekistan, Kyrgyzstan, Turkmenistan

**Asia:**

- `south_asia`: India, Pakistan, Bangladesh, Sri Lanka, Nepal

- `east_asia`: China, Japan, Korea, Mongolia

- `southeast_asia`: Vietnam, Thailand, Indonesia, Philippines, Malaysia

**Africa:**

- `west_africa`: Nigeria, Ghana, Senegal, Mali

- `east_africa`: Kenya, Ethiopia, Tanzania, Somalia, Uganda

- `southern_africa`: South Africa, Zimbabwe, Botswana, Namibia

- `north_africa`: Morocco, Algeria, Tunisia, Libya

- `central_africa`: DRC, Cameroon, Congo

**Americas:**

- `north_america`: USA, Canada

- `latin_america`: Mexico, Central America, South America, Brazil

- `caribbean`: Jamaica, Haiti, Cuba, Trinidad, Dominican Republic

- `indigenous_americas`: Native American, First Nations, Mayan, Incan

**Oceania:**

- `australia_nz`: Australia, New Zealand

- `pacific_islands`: Hawaii, Samoa, Tonga, Fiji, PNG

- `indigenous_oceania`: Maori, Aboriginal Australian

**Other:**

- `jewish_hebrew`: Hebrew/Jewish names

- `armenian_georgian`: Armenia, Georgia

## See also

[`get_global_names`](https://pak.dynasite.org/Saqrlab/reference/get_global_names.md)
for sampling names,
[`list_name_regions`](https://pak.dynasite.org/Saqrlab/reference/list_name_regions.md)
for exploring available regions.

## Examples

``` r
# View all region names
names(GLOBAL_NAMES)
#>  [1] "western_europe"      "eastern_europe"      "southern_europe"    
#>  [4] "nordic"              "arab"                "persian"            
#>  [7] "turkish"             "central_asia"        "south_asia"         
#> [10] "east_asia"           "southeast_asia"      "west_africa"        
#> [13] "east_africa"         "southern_africa"     "north_africa"       
#> [16] "central_africa"      "north_america"       "latin_america"      
#> [19] "caribbean"           "indigenous_americas" "australia_nz"       
#> [22] "pacific_islands"     "indigenous_oceania"  "jewish_hebrew"      
#> [25] "armenian_georgian"  

# View names from a specific region
GLOBAL_NAMES$arab
#>  [1] "Ali"     "Fatima"  "Omar"    "Layla"   "Yusuf"   "Mariam"  "Khalid" 
#>  [8] "Noor"    "Sultan"  "Sheikha" "Rashid"  "Amal"    "Hamad"   "Moza"   
#> [15] "Faisal"  "Hessa"   "Hassan"  "Amira"   "Amir"    "Zara"    "Tariq"  
#> [22] "Salma"   "Faris"   "Rania"   "Karim"   "Leila"   "Walid"   "Dina"   
#> [29] "Nabil"   "Huda"    "Samir"   "Lina"    "Ahmed"   "Mona"    "Mahmoud"
#> [36] "Heba"    "Mostafa" "Yasmine" "Tarek"   "Dalia"   "Mohamed" "Noha"   
#> [43] "Amr"     "Aya"     "Khaled"  "Mai"     "Hazem"   "Nada"    "Ziad"   
#> [50] "Rana"    "Sami"    "Reem"    "Adel"    "Noura"   "Mazen"   "Asma"   
#> [57] "Jamal"   "Dalal"   "Bassam"  "Aisha"   "Rami"    "Hana"    "Tamer"  
#> [64] "Sara"   
GLOBAL_NAMES$east_asia
#>   [1] "Wei"          "Lin"          "Chen"         "Mei"          "Jing"        
#>   [6] "Xiao"         "Ming"         "Hua"          "Feng"         "Yan"         
#>  [11] "Jun"          "Lei"          "Bo"           "Lan"          "Tao"         
#>  [16] "Ying"         "Hao"          "Yue"          "Zhen"         "Xin"         
#>  [21] "Peng"         "Qian"         "Rui"          "Shu"          "Long"        
#>  [26] "Fang"         "Hong"         "Li"           "Chao"         "Na"          
#>  [31] "Gang"         "Ling"         "Yuki"         "Hiro"         "Akira"       
#>  [36] "Sakura"       "Kenji"        "Yuna"         "Taro"         "Haruki"      
#>  [41] "Koji"         "Miki"         "Ryu"          "Emi"          "Rin"         
#>  [46] "Sora"         "Kaito"        "Hana"         "Yuto"         "Aoi"         
#>  [51] "Sota"         "Mei"          "Ren"          "Yui"          "Haruto"      
#>  [56] "Koharu"       "Takumi"       "Himari"       "Kenta"        "Akari"       
#>  [61] "Daiki"        "Riko"         "Naoki"        "Ayumi"        "Minho"       
#>  [66] "Jisoo"        "Joon"         "Seo"          "Tae"          "Yuna"        
#>  [71] "Woo"          "Minji"        "Dae"          "Eunji"        "Sung"        
#>  [76] "Hana"         "Hyun"         "Jiwon"        "Seok"         "Nari"        
#>  [81] "Jiho"         "Somin"        "Junho"        "Yuri"         "Minsoo"      
#>  [86] "Soyeon"       "Donghyun"     "Chaeyoung"    "Bataar"       "Oyun"        
#>  [91] "Bold"         "Altai"        "Temuulen"     "Sarnai"       "Ganzorig"    
#>  [96] "Enkhtuya"     "Erdene"       "Narantsetseg" "Munkh"        "Tsetseg"     
#> [101] "Baatar"       "Oyunbileg"    "Chuluun"      "Bolormaa"    

# Count names per region
sapply(GLOBAL_NAMES, length)
#>      western_europe      eastern_europe     southern_europe              nordic 
#>                  64                  72                  72                  80 
#>                arab             persian             turkish        central_asia 
#>                  64                  48                  48                  48 
#>          south_asia           east_asia      southeast_asia         west_africa 
#>                  96                 104                  96                  64 
#>         east_africa     southern_africa        north_africa      central_africa 
#>                  56                  48                  48                  40 
#>       north_america       latin_america           caribbean indigenous_americas 
#>                  48                  72                  48                  48 
#>        australia_nz     pacific_islands  indigenous_oceania       jewish_hebrew 
#>                  40                  48                  48                  40 
#>   armenian_georgian 
#>                  32 

# Total names available
sum(sapply(GLOBAL_NAMES, length))
#> [1] 1472
```
