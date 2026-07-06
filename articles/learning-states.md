# Learning States Reference

``` r

library(Saqrlab)
#> 
#> Attaching package: 'Saqrlab'
#> The following object is masked from 'package:stats':
#> 
#>     simulate
```

## Overview

Saqrlab includes over 180 learning action verbs organized into 8
categories, designed for educational research simulations. This
reference documents all available learning states and how to use them.

## Available Categories

``` r

list_learning_categories()
#>           Category Count
#> 1    metacognitive    20
#> 2        cognitive    30
#> 3       behavioral    30
#> 4           social    30
#> 5     motivational    30
#> 6        affective    30
#> 7 group_regulation     9
#> 8              lms    30
#>                                                  Examples
#> 1         Plan, Monitor, Evaluate, Reflect, Regulate, ...
#> 2          Read, Study, Analyze, Summarize, Memorize, ...
#> 3       Practice, Annotate, Research, Review, Revise, ...
#> 4 Collaborate, Discuss, Seek_help, Question, Explain, ...
#> 5            Focus, Persist, Explore, Create, Strive, ...
#> 6        Enjoy, Appreciate, Value, Interest, Curious, ...
#> 7    Adapt, Cohesion, Consensus, Coregulate, Discuss, ...
#> 8             View, Access, Download, Upload, Submit, ...
```

## Complete Learning States Reference

### Metacognitive (20 verbs)

Self-regulation and awareness actions for planning, monitoring, and
evaluating one’s own learning.

``` r

get_learning_states("metacognitive")
#>  [1] "Plan"        "Monitor"     "Evaluate"    "Reflect"     "Regulate"   
#>  [6] "Adjust"      "Adapt"       "Check"       "Assess"      "Judge"      
#> [11] "Strategize"  "Prioritize"  "Set_goals"   "Track"       "Self_assess"
#> [16] "Calibrate"   "Diagnose"    "Forecast"    "Anticipate"  "Reconsider"
```

### Cognitive (30 verbs)

Mental processing actions for understanding, analyzing, and synthesizing
information.

``` r

get_learning_states("cognitive")
#>  [1] "Read"          "Study"         "Analyze"       "Summarize"    
#>  [5] "Memorize"      "Connect"       "Apply"         "Comprehend"   
#>  [9] "Synthesize"    "Compare"       "Contrast"      "Infer"        
#> [13] "Interpret"     "Elaborate"     "Encode"        "Retrieve"     
#> [17] "Process"       "Understand"    "Learn"         "Recognize"    
#> [21] "Recall"        "Integrate"     "Differentiate" "Abstract"     
#> [25] "Generalize"    "Classify"      "Categorize"    "Deduce"       
#> [29] "Reason"        "Conclude"
```

### Behavioral (30 verbs)

Observable study actions including reading, writing, practicing, and
note-taking.

``` r

get_learning_states("behavioral")
#>  [1] "Practice"  "Annotate"  "Research"  "Review"    "Revise"    "Test"     
#>  [7] "Write"     "Note"      "Highlight" "Underline" "Reread"    "Skim"     
#> [13] "Scan"      "Draft"     "Edit"      "Copy"      "Record"    "Complete" 
#> [19] "Submit"    "Attempt"   "Repeat"    "Drill"     "Exercise"  "Rehearse" 
#> [25] "Outline"   "Diagram"   "Map"       "List"      "Organize"  "Structure"
```

### Social (30 verbs)

Interaction and collaboration actions for discussing, helping, teaching,
and asking.

``` r

get_learning_states("social")
#>  [1] "Collaborate" "Discuss"     "Seek_help"   "Question"    "Explain"    
#>  [6] "Share"       "Teach"       "Tutor"       "Debate"      "Argue"      
#> [11] "Negotiate"   "Consult"     "Ask"         "Answer"      "Present"    
#> [16] "Participate" "Engage"      "Contribute"  "Support"     "Help"       
#> [21] "Clarify"     "Communicate" "Listen"      "Respond"     "Feedback"   
#> [26] "Critique"    "Peer_review" "Co_create"   "Brainstorm"  "Network"
```

### Motivational (30 verbs)

Effort and persistence actions for focusing, persisting, engaging, and
striving.

``` r

get_learning_states("motivational")
#>  [1] "Focus"       "Persist"     "Explore"     "Create"      "Strive"     
#>  [6] "Commit"      "Motivate"    "Endure"      "Overcome"    "Challenge"  
#> [11] "Aspire"      "Dedicate"    "Invest"      "Concentrate" "Attend"     
#> [16] "Sustain"     "Maintain"    "Initiate"    "Continue"    "Pursue"     
#> [21] "Drive"       "Hustle"      "Push"        "Achieve"     "Accomplish" 
#> [26] "Excel"       "Improve"     "Grow"        "Develop"     "Progress"
```

### Affective (30 verbs)

Emotional and attitudinal states including enjoying, coping, managing
stress, and curiosity.

``` r

get_learning_states("affective")
#>  [1] "Enjoy"      "Appreciate" "Value"      "Interest"   "Curious"   
#>  [6] "Worry"      "Stress"     "Relax"      "Cope"       "Manage"    
#> [11] "Calm"       "Frustrate"  "Satisfy"    "Excite"     "Bore"      
#> [16] "Confuse"    "Resolve"    "Embrace"    "Accept"     "Tolerate"  
#> [21] "Celebrate"  "Doubt"      "Confident"  "Anxious"    "Hopeful"   
#> [26] "Discourage" "Encourage"  "Inspire"    "Overwhelm"  "Relief"
```

### Group Regulation (9 verbs)

Socially shared regulation of learning (SSRL) actions for collaborative
planning, monitoring, and adapting.

``` r

get_learning_states("group_regulation")
#> [1] "Adapt"      "Cohesion"   "Consensus"  "Coregulate" "Discuss"   
#> [6] "Emotion"    "Monitor"    "Plan"       "Synthesis"
```

### LMS Actions (30 verbs)

Learning Management System interactions including viewing content,
accessing resources, and submitting assignments.

``` r

get_learning_states("lms")
#>  [1] "View"         "Access"       "Download"     "Upload"       "Submit"      
#>  [6] "Click"        "Navigate"     "Browse"       "Login"        "Logout"      
#> [11] "Post"         "Reply"        "Forum"        "Quiz"         "Assignment"  
#> [16] "Video"        "Resource"     "Grade"        "Attempt"      "Complete"    
#> [21] "Module"       "Page"         "File"         "Link"         "Course"      
#> [26] "Content"      "Discussion"   "Message"      "Announcement" "Calendar"
```

## Using Learning States

### Basic Retrieval

``` r

# Get all verbs from a category
metacognitive_verbs <- get_learning_states("metacognitive")
length(metacognitive_verbs)
#> [1] 20

# Get verbs from multiple categories
srl_verbs <- get_learning_states(c("metacognitive", "cognitive", "motivational"))
length(srl_verbs)
#> [1] 80

# Get all verbs
all_verbs <- get_learning_states("all")
length(all_verbs)
#> [1] 202
```

### Random Sampling

``` r

# Random 8 verbs from any category
random_8 <- get_learning_states(n = 8, seed = 42)
random_8
#> [1] "Reason"   "Edit"     "Satisfy"  "Rehearse" "Worry"    "Dedicate" "Calendar"
#> [8] "Initiate"

# Random 6 verbs from specific categories
random_cog <- get_learning_states(c("cognitive", "behavioral"), n = 6, seed = 42)
random_cog
#> [1] "Submit"     "Write"      "Read"       "Generalize" "Compare"   
#> [6] "Test"
```

### Smart Selection

The
[`select_states()`](https://pak.dynasite.org/Saqrlab/reference/select_states.md)
function intelligently selects states based on network size:

``` r

# Small network: auto-selects from one category
small_net <- select_states(5, seed = 42)
small_net
#> [1] "Regulate"   "Discourage" "Judge"      "Plan"       "Appreciate"

# Medium network: balanced across categories
medium_net <- select_states(10, seed = 42)
medium_net
#>  [1] "Sustain"    "Strive"     "Negotiate"  "Judge"      "Plan"      
#>  [6] "Create"     "Manage"     "Forecast"   "Tolerate"   "Brainstorm"

# Large network: wide variety
large_net <- select_states(20, seed = 42)
large_net
#>  [1] "Link"       "Maintain"   "Forecast"   "Strive"     "Commit"    
#>  [6] "Memorize"   "Categorize" "Evaluate"   "Create"     "Calm"      
#> [11] "Test"       "Develop"    "Reconsider" "Drive"      "Respond"   
#> [16] "Improve"    "Brainstorm" "Curious"    "Module"     "Accept"
```

### Biased Selection

``` r

# Prioritize metacognitive verbs
meta_focus <- select_states(
  n_states = 10,
  primary_categories = "metacognitive",
  secondary_categories = "cognitive",
  primary_ratio = 0.7,  # 70% metacognitive
  seed = 42
)
meta_focus
#>  [1] "Reconsider" "Judge"      "Diagnose"   "Reflect"    "Monitor"   
#>  [6] "Encode"     "Regulate"   "Abstract"   "Process"    "Plan"
```

## Integration with Simulation Functions

### With simulate_sequences()

``` r

# Auto-generates with learning states by default
sequences <- simulate_sequences(
  n_sequences = 50,
  seq_length = 15,
  n_states = 5,
  seed = 42
)

# Check states used
unique(unlist(sequences))
#> [1] "Appreciate" "Judge"      "Discourage" "Plan"       "Regulate"
```

### Specifying Categories

``` r

# Use specific categories
sequences <- simulate_sequences(
  n_sequences = 50,
  seq_length = 15,
  n_states = 6,
  categories = c("metacognitive", "cognitive"),
  seed = 42
)

unique(unlist(sequences))
#> [1] "Plan"     "Process"  "Retrieve" "Memorize" "Reason"   "Judge"
```

### With simulate_matrix()

``` r

# Matrix uses random category by default
mat <- simulate_matrix(n_nodes = 5, seed = 42)
rownames(mat)
#> [1] "Regulate" "Plan"     "Judge"    "Reflect"  "Monitor"

# Run again - different category selected
mat2 <- simulate_matrix(n_nodes = 5, seed = 123)
rownames(mat2)
#> [1] "Consensus" "Emotion"   "Synthesis" "Cohesion"  "Plan"
```

### With simulate_htna()

``` r

# HTNA uses different category per type
net <- simulate_htna(seed = 42)
lapply(net$node_types, head, 3)
#> $Metacognitive
#> [1] "Diagnose" "Regulate" "Plan"    
#> 
#> $Cognitive
#> [1] "Understand" "Classify"   "Process"   
#> 
#> $Behavioral
#> [1] "Write"   "Review"  "Outline"
#> 
#> $Social
#> [1] "Help"       "Critique"   "Contribute"
#> 
#> $Motivational
#> [1] "Overcome"   "Accomplish" "Improve"
```

## Direct Access to LEARNING_STATES

``` r

# Access the raw data
names(LEARNING_STATES)
#> [1] "metacognitive"    "cognitive"        "behavioral"       "social"          
#> [5] "motivational"     "affective"        "group_regulation" "lms"

# Get specific category
LEARNING_STATES$metacognitive
#>  [1] "Plan"        "Monitor"     "Evaluate"    "Reflect"     "Regulate"   
#>  [6] "Adjust"      "Adapt"       "Check"       "Assess"      "Judge"      
#> [11] "Strategize"  "Prioritize"  "Set_goals"   "Track"       "Self_assess"
#> [16] "Calibrate"   "Diagnose"    "Forecast"    "Anticipate"  "Reconsider"

# Count verbs per category
sapply(LEARNING_STATES, length)
#>    metacognitive        cognitive       behavioral           social 
#>               20               30               30               30 
#>     motivational        affective group_regulation              lms 
#>               30               30                9               30
```

## Global Names Dataset

Saqrlab also includes 300 diverse names for simulation:

``` r

# Access names
head(GLOBAL_NAMES, 20)
#> $western_europe
#>  [1] "Emma"     "Liam"     "Chloe"    "Jack"     "Lily"     "Oliver"  
#>  [7] "Grace"    "Harry"    "Aoife"    "Sean"     "Niamh"    "Finn"    
#> [13] "Saoirse"  "Ciaran"   "Roisin"   "Declan"   "Lea"      "Hugo"    
#> [19] "Manon"    "Lucas"    "Chloe"    "Jules"    "Camille"  "Louis"   
#> [25] "Eloise"   "Antoine"  "Amelie"   "Theo"     "Juliette" "Raphael" 
#> [31] "Margot"   "Adrien"   "Max"      "Mia"      "Felix"    "Lena"    
#> [37] "Leon"     "Anna"     "Paul"     "Clara"    "Lukas"    "Sophie"  
#> [43] "Jonas"    "Lea"      "Elias"    "Marie"    "Noah"     "Emilia"  
#> [49] "Daan"     "Sanne"    "Sem"      "Julia"    "Lotte"    "Ruben"   
#> [55] "Fien"     "Bram"     "Thijs"    "Eva"      "Lars"     "Lisa"    
#> [61] "Milan"    "Noor"     "Tim"      "Fleur"   
#> 
#> $eastern_europe
#>  [1] "Ivan"    "Olga"    "Dmitri"  "Natasha" "Alexei"  "Katya"   "Yuri"   
#>  [8] "Anya"    "Boris"   "Mila"    "Sasha"   "Vera"    "Vlad"    "Daria"  
#> [15] "Oleg"    "Irina"   "Nikita"  "Lena"    "Andrei"  "Tanya"   "Maxim"  
#> [22] "Yulia"   "Artem"   "Sveta"   "Jakub"   "Zuzanna" "Kacper"  "Maja"   
#> [29] "Antoni"  "Hanna"   "Piotr"   "Ewa"     "Michal"  "Anna"    "Tomek"  
#> [36] "Kasia"   "Bartek"  "Ola"     "Pawel"   "Magda"   "Jiri"    "Petra"  
#> [43] "Tomas"   "Hana"    "Marek"   "Eva"     "Milan"   "Jana"    "Ondrej" 
#> [50] "Lucie"   "Adam"    "Tereza"  "Filip"   "Klara"   "Jakub"   "Monika" 
#> [57] "Andrei"  "Elena"   "Mircea"  "Raluca"  "Bela"    "Eszter"  "Matyas" 
#> [64] "Zsofia"  "Ion"     "Maria"   "Vlad"    "Ana"     "Attila"  "Kata"   
#> [71] "Levente" "Reka"   
#> 
#> $southern_europe
#>  [1] "Pablo"      "Lucia"      "Diego"      "Sofia"      "Alvaro"    
#>  [6] "Carmen"     "Sergio"     "Marta"      "Javier"     "Ana"       
#> [11] "Carlos"     "Laura"      "Alejandro"  "Elena"      "Daniel"    
#> [16] "Paula"      "Raul"       "Ines"       "Adrian"     "Alba"      
#> [21] "Jorge"      "Nuria"      "Ivan"       "Cristina"   "Marco"     
#> [26] "Giulia"     "Luca"       "Francesca"  "Matteo"     "Chiara"    
#> [31] "Andrea"     "Sara"       "Alessandro" "Valentina"  "Lorenzo"   
#> [36] "Martina"    "Davide"     "Federica"   "Giuseppe"   "Elisa"     
#> [41] "Joao"       "Maria"      "Tiago"      "Ines"       "Diogo"     
#> [46] "Beatriz"    "Pedro"      "Ana"        "Miguel"     "Mariana"   
#> [51] "Rui"        "Catarina"   "Andre"      "Sofia"      "Bruno"     
#> [56] "Rita"       "Nikos"      "Elena"      "Dimitris"   "Maria"     
#> [61] "Kostas"     "Sofia"      "Yannis"     "Anna"       "Giorgos"   
#> [66] "Eleni"      "Christos"   "Katerina"   "Alexandros" "Ioanna"    
#> [71] "Stavros"    "Dimitra"   
#> 
#> $nordic
#>  [1] "Erik"      "Astrid"    "Oscar"     "Maja"      "Axel"      "Ella"     
#>  [7] "Viktor"    "Saga"      "Gustav"    "Wilma"     "William"   "Ebba"     
#> [13] "Linus"     "Agnes"     "Filip"     "Freja"     "Lars"      "Ingrid"   
#> [19] "Sven"      "Liv"       "Olaf"      "Nora"      "Leif"      "Thea"     
#> [25] "Magnus"    "Emma"      "Henrik"    "Sofie"     "Andreas"   "Ida"      
#> [31] "Kristian"  "Julie"     "Mikkel"    "Freya"     "Emil"      "Clara"    
#> [37] "Christian" "Emilie"    "Noah"      "Alma"      "Oliver"    "Freja"    
#> [43] "Lucas"     "Laura"     "Frederik"  "Anna"      "Mathias"   "Maja"     
#> [49] "Onni"      "Aino"      "Eetu"      "Venla"     "Lauri"     "Helmi"    
#> [55] "Eero"      "Siiri"     "Matias"    "Emilia"    "Aleksi"    "Sofia"    
#> [61] "Veeti"     "Aada"      "Elias"     "Ella"      "Ragnar"    "Sigrid"   
#> [67] "Bjorn"     "Gudrun"    "Eirik"     "Hildur"    "Odin"      "Frida"    
#> [73] "Thor"      "Helga"     "Gunnar"    "Unnur"     "Bjarni"    "Kristin"  
#> [79] "Arnar"     "Margret"  
#> 
#> $arab
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
#> 
#> $persian
#>  [1] "Cyrus"    "Shirin"   "Dara"     "Parisa"   "Reza"     "Maryam"  
#>  [7] "Arman"    "Soraya"   "Babak"    "Mina"     "Kaveh"    "Nazanin" 
#> [13] "Omid"     "Sara"     "Farhad"   "Azar"     "Behzad"   "Niloufar"
#> [19] "Kian"     "Mahsa"    "Arash"    "Pari"     "Dariush"  "Roxana"  
#> [25] "Amir"     "Leyla"    "Kamran"   "Shadi"    "Nima"     "Azadeh"  
#> [31] "Siavash"  "Setareh"  "Ahmad"    "Fatima"   "Farid"    "Zainab"  
#> [37] "Hamid"    "Maryam"   "Jawad"    "Hosna"    "Rustam"   "Gulnora" 
#> [43] "Firdavs"  "Munis"    "Jamshed"  "Dilbar"   "Firuz"    "Sitora"  
#> 
#> $turkish
#>  [1] "Emre"    "Elif"    "Kaan"    "Defne"   "Cem"     "Zeynep"  "Deniz"  
#>  [8] "Aylin"   "Burak"   "Selin"   "Mert"    "Ceren"   "Baris"   "Ebru"   
#> [15] "Arda"    "Melis"   "Kemal"   "Leyla"   "Ozan"    "Nehir"   "Alp"    
#> [22] "Yasemin" "Umut"    "Ezgi"    "Onur"    "Duygu"   "Tolga"   "Pinar"  
#> [29] "Serkan"  "Burcu"   "Hakan"   "Esra"    "Eldar"   "Gunel"   "Tural"  
#> [36] "Nigar"   "Orhan"   "Sevil"   "Farid"   "Aynur"   "Rauf"    "Gulnar" 
#> [43] "Ilkin"   "Samira"  "Nurlan"  "Leyla"   "Vugar"   "Arzu"   
#> 
#> $central_asia
#>  [1] "Aibek"     "Aida"      "Nurlan"    "Gulnara"   "Bekzat"    "Aigerim"  
#>  [7] "Ruslan"    "Dana"      "Nursultan" "Ainur"     "Almas"     "Dilnaz"   
#> [13] "Miras"     "Zhanna"    "Saken"     "Aliya"     "Bobur"     "Nilufar"  
#> [19] "Sardor"    "Malika"    "Jasur"     "Dilnoza"   "Nodir"     "Zarina"   
#> [25] "Bakyt"     "Madina"    "Azamat"    "Asel"      "Damir"     "Cholpon"  
#> [31] "Timur"     "Sabina"    "Merdan"    "Mahri"     "Oraz"      "Ogulgerek"
#> [37] "Serdar"    "Aylar"     "Dovlet"    "Jennet"    "Bahrom"    "Zebo"     
#> [43] "Farkhod"   "Nigina"    "Suhrab"    "Madina"    "Daler"     "Firuza"   
#> 
#> $south_asia
#>  [1] "Arun"    "Priya"   "Raj"     "Devi"    "Amit"    "Sita"    "Ravi"   
#>  [8] "Maya"    "Ajay"    "Lata"    "Vijay"   "Anita"   "Rohit"   "Pooja"  
#> [15] "Arjun"   "Kavita"  "Rahul"   "Sunita"  "Kiran"   "Aditi"   "Nikhil" 
#> [22] "Shreya"  "Vikram"  "Neha"    "Sanjay"  "Meera"   "Aditya"  "Divya"  
#> [29] "Varun"   "Isha"    "Dev"     "Tara"    "Surya"   "Lakshmi" "Karthik"
#> [36] "Deepa"   "Ganesh"  "Radha"   "Suresh"  "Vasuki"  "Rajan"   "Chitra" 
#> [43] "Ashok"   "Padma"   "Senthil" "Bhavani" "Mohan"   "Uma"     "Imran"  
#> [50] "Ayesha"  "Bilal"   "Sana"    "Hamza"   "Hira"    "Zain"    "Maryam" 
#> [57] "Usman"   "Fatima"  "Ali"     "Zoya"    "Hasan"   "Amina"   "Faisal" 
#> [64] "Rabia"   "Rafiq"   "Fatema"  "Rahim"   "Nasreen" "Shahid"  "Rupa"   
#> [71] "Kamal"   "Shanta"  "Tanvir"  "Nusrat"  "Mahbub"  "Tahera"  "Iqbal"  
#> [78] "Hasina"  "Jamal"   "Rima"    "Nimal"   "Kumari"  "Sandun"  "Dilani" 
#> [85] "Chamara" "Nimali"  "Lasith"  "Thilini" "Binod"   "Sushma"  "Rajan"  
#> [92] "Kamala"  "Sunil"   "Gita"    "Ramesh"  "Sabita" 
#> 
#> $east_asia
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
#> 
#> $southeast_asia
#>  [1] "Linh"     "Minh"     "Anh"      "Bao"      "Mai"      "Duc"     
#>  [7] "Thao"     "Tuan"     "Hoa"      "Nam"      "Lan"      "Hung"    
#> [13] "Hanh"     "Cuong"    "Ngoc"     "Phong"    "Trang"    "Huy"     
#> [19] "Nhi"      "Khanh"    "Vy"       "Long"     "Huong"    "Dung"    
#> [25] "Somchai"  "Ploy"     "Chai"     "Nong"     "Kiet"     "Mali"    
#> [31] "Somsak"   "Araya"    "Prem"     "Dao"      "Niran"    "Malai"   
#> [37] "Sakchai"  "Achara"   "Prasert"  "Supaporn" "Rizal"    "Putri"   
#> [43] "Bagus"    "Sari"     "Wayan"    "Ayu"      "Gede"     "Nisa"    
#> [49] "Dewi"     "Budi"     "Rani"     "Eko"      "Sinta"    "Agus"    
#> [55] "Ratna"    "Dimas"    "Adi"      "Fitri"    "Rudi"     "Wulan"   
#> [61] "Hendra"   "Lestari"  "Yusuf"    "Mega"     "Jose"     "Maria"   
#> [67] "Juan"     "Rosa"     "Carlo"    "Ana"      "Miguel"   "Luz"     
#> [73] "Antonio"  "Carmen"   "Ramon"    "Liza"     "Francis"  "Grace"   
#> [79] "Paolo"    "Joy"      "Fajar"    "Indah"    "Hafiz"    "Nurul"   
#> [85] "Amir"     "Siti"     "Ismail"   "Aishah"   "Zul"      "Farah"   
#> [91] "Azlan"    "Aisyah"   "Firdaus"  "Nur"      "Hakim"    "Amira"   
#> 
#> $west_africa
#>  [1] "Chidi"     "Adaora"    "Emeka"     "Ngozi"     "Obi"       "Amaka"    
#>  [7] "Chijioke"  "Nneka"     "Ikenna"    "Chinwe"    "Obinna"    "Uchenna"  
#> [13] "Nnamdi"    "Chinyere"  "Kelechi"   "Adaobi"    "Ade"       "Yemi"     
#> [19] "Tunde"     "Funke"     "Segun"     "Bola"      "Kunle"     "Nike"     
#> [25] "Dayo"      "Sade"      "Kayode"    "Toyin"     "Femi"      "Lola"     
#> [31] "Bayo"      "Ronke"     "Kofi"      "Ama"       "Kwame"     "Akua"     
#> [37] "Yaw"       "Abena"     "Kwesi"     "Efua"      "Kojo"      "Akosua"   
#> [43] "Kwadwo"    "Adwoa"     "Kweku"     "Afua"      "Kobi"      "Afia"     
#> [49] "Amadou"    "Fatou"     "Moussa"    "Aissatou"  "Oumar"     "Mariama"  
#> [55] "Sekou"     "Binta"     "Ibrahima"  "Aminata"   "Mamadou"   "Kadiatou" 
#> [61] "Boubacar"  "Fatoumata" "Cheikh"    "Oumou"    
#> 
#> $east_africa
#>  [1] "Juma"      "Amina"     "Bakari"    "Zuri"      "Mwangi"    "Wanjiku"  
#>  [7] "Kamau"     "Njeri"     "Ochieng"   "Auma"      "Otieno"    "Adhiambo" 
#> [13] "Kipchoge"  "Chebet"    "Korir"     "Jepkosgei" "Baraka"    "Neema"    
#> [19] "Hamisi"    "Saida"     "Rajabu"    "Rehema"    "Salum"     "Mwajuma"  
#> [25] "Abebe"     "Tigist"    "Haile"     "Makeda"    "Tadesse"   "Meron"    
#> [31] "Dawit"     "Sara"      "Yohannes"  "Bethlehem" "Girma"     "Selam"    
#> [37] "Tesfaye"   "Hana"      "Solomon"   "Rahel"     "Mohamed"   "Halima"   
#> [43] "Abdi"      "Amina"     "Hassan"    "Fartun"    "Omar"      "Sahra"    
#> [49] "Mugisha"   "Uwimana"   "Kato"      "Akello"    "Okello"    "Nkunda"   
#> [55] "Kagame"    "Umutesi"  
#> 
#> $southern_africa
#>  [1] "Themba"     "Lindiwe"    "Sipho"      "Nomzamo"    "Thabo"     
#>  [6] "Zanele"     "Mandla"     "Lerato"     "Sibusiso"   "Thandiwe"  
#> [11] "Bongani"    "Nolwazi"    "Siyabonga"  "Ayanda"     "Luyanda"   
#> [16] "Nomvula"    "Mpho"       "Kelebogile" "Neo"        "Goitseone" 
#> [21] "Tshepo"     "Dineo"      "Thato"      "Palesa"     "Tendai"    
#> [26] "Rudo"       "Tapiwa"     "Chipo"      "Tinashe"    "Tariro"    
#> [31] "Farai"      "Nyasha"     "Tawanda"    "Rutendo"    "Tatenda"   
#> [36] "Rumbidzai"  "Kudakwashe" "Ropafadzo"  "Tanaka"     "Rufaro"    
#> [41] "Mulenga"    "Mutinta"    "Chisomo"    "Thandie"    "Kgosi"     
#> [46] "Naledi"     "Amantle"    "Masego"    
#> 
#> $north_africa
#>  [1] "Youssef" "Fatima"  "Hassan"  "Amina"   "Said"    "Khadija" "Rachid" 
#>  [8] "Zahra"   "Mehdi"   "Salma"   "Amine"   "Houda"   "Yassine" "Imane"  
#> [15] "Hamza"   "Hajar"   "Malik"   "Samira"  "Kamel"   "Nadia"   "Hamid"  
#> [22] "Leila"   "Sofiane" "Salima"  "Riad"    "Yasmina" "Farid"   "Lamia"  
#> [29] "Nassim"  "Amira"   "Djamel"  "Sihem"   "Hedi"    "Sonia"   "Nizar"  
#> [36] "Amel"    "Fares"   "Meriem"  "Slim"    "Ines"    "Mukhtar" "Hana"   
#> [43] "Idris"   "Asma"    "Fathi"   "Marwa"   "Nuri"    "Salwa"  
#> 
#> $central_africa
#>  [1] "Patrice"    "Solange"    "Jean"       "Marie"      "Emmanuel"  
#>  [6] "Grace"      "David"      "Ruth"       "Fiston"     "Carine"    
#> [11] "Christian"  "Nadine"     "Patrick"    "Sylvie"     "Serge"     
#> [16] "Claudine"   "Samuel"     "Esther"     "Paul"       "Christiane"
#> [21] "Pierre"     "Nadege"     "Andre"      "Simone"     "Yves"      
#> [26] "Blanche"    "Roger"      "Jeanne"     "Fabrice"    "Mireille"  
#> [31] "Herve"      "Colette"    "Arsene"     "Prudence"   "Guy"       
#> [36] "Diane"      "Didier"     "Lydie"      "Landry"     "Ornella"   
#> 
#> $north_america
#>  [1] "James"     "Emily"     "Michael"   "Sarah"     "David"     "Jessica"  
#>  [7] "John"      "Ashley"    "Robert"    "Amanda"    "William"   "Jennifer" 
#> [13] "Chris"     "Nicole"    "Matt"      "Lauren"    "Ryan"      "Rachel"   
#> [19] "Tyler"     "Megan"     "Brian"     "Stephanie" "Kevin"     "Brittany" 
#> [25] "Andrew"    "Samantha"  "Justin"    "Elizabeth" "Brandon"   "Heather"  
#> [31] "Josh"      "Michelle"  "Liam"      "Emma"      "Noah"      "Olivia"   
#> [37] "Ethan"     "Ava"       "Mason"     "Sophia"    "Logan"     "Charlotte"
#> [43] "Jacob"     "Amelia"    "Lucas"     "Harper"    "Jack"      "Evelyn"   
#> 
#> $latin_america
#>  [1] "Juan"      "Rosa"      "Diego"     "Luz"       "Pablo"     "Sol"      
#>  [7] "Luis"      "Ana"       "Carlos"    "Lucia"     "Miguel"    "Elena"    
#> [13] "Jose"      "Sofia"     "Pedro"     "Camila"    "Alejandro" "Valentina"
#> [19] "Fernando"  "Daniela"   "Ricardo"   "Mariana"   "Eduardo"   "Fernanda" 
#> [25] "Mateo"     "Isabella"  "Rafael"    "Isabel"    "Sergio"    "Carmen"   
#> [31] "Alvaro"    "Paulina"   "Andres"    "Maria"     "Santiago"  "Martina"  
#> [37] "Gabriel"   "Florencia" "Nicolas"   "Julieta"   "Sebastian" "Catalina" 
#> [43] "Benjamin"  "Antonia"   "Martin"    "Emilia"    "Tomas"     "Isidora"  
#> [49] "Joao"      "Julia"     "Pedro"     "Fernanda"  "Lucas"     "Amanda"   
#> [55] "Gustavo"   "Beatriz"   "Rafael"    "Larissa"   "Thiago"    "Bruna"    
#> [61] "Matheus"   "Carolina"  "Vitor"     "Leticia"   "Bruno"     "Mariana"  
#> [67] "Felipe"    "Gabriela"  "Leonardo"  "Camila"    "Gabriel"   "Isabela"  
#> 
#> $caribbean
#>  [1] "Marlon"    "Keisha"    "Dwayne"    "Shanique"  "Andre"     "Natalie"  
#>  [7] "Wayne"     "Tanya"     "Usain"     "Shelly"    "Damian"    "Khadija"  
#> [13] "Tristan"   "Yolanda"   "Jermaine"  "Sasha"     "Jean"      "Marie"    
#> [19] "Pierre"    "Rose"      "Jacques"   "Nadine"    "Claude"    "Carole"   
#> [25] "Wyclef"    "Michaelle" "Stanley"   "Fabienne"  "Herby"     "Guerda"   
#> [31] "Frantz"    "Ketty"     "Alejandro" "Yolanda"   "Orlando"   "Marisol"  
#> [37] "Ramon"     "Yesenia"   "Jorge"     "Luz"       "Lazaro"    "Yanelis"  
#> [43] "Yunel"     "Dania"     "Pedro"     "Yuliesky"  "Miguel"    "Yamilet"  
#> 
#> $indigenous_americas
#>  [1] "Takoda"     "Winona"     "Koda"       "Aiyana"     "Chayton"   
#>  [6] "Aponi"      "Ahanu"      "Chenoa"     "Hinto"      "Halona"    
#> [11] "Shilah"     "Mika"       "Nayeli"     "Kimi"       "Huritt"    
#> [16] "Ayita"      "Sequoia"    "Cochise"    "Sacagawea"  "Tecumseh"  
#> [21] "Hiawatha"   "Pocahontas" "Dakota"     "Cheyenne"   "Taima"     
#> [26] "Kaya"       "Nanuq"      "Sedna"      "Amaruq"     "Siku"      
#> [31] "Tulok"      "Atka"       "Itzamna"    "Xochitl"    "Cuauhtemoc"
#> [36] "Citlali"    "Quetzal"    "Itzel"      "Tlaloc"     "Ixchel"    
#> [41] "Tupac"      "Qori"       "Inti"       "Sumaq"      "Amaru"     
#> [46] "Killa"      "Rumi"       "Wayra"

# Get random names
get_global_names(n = 10, seed = 42)
#>  [1] "Shu"        "Soraya"     "Rangi"      "Cuauhtemoc" "Gagik"     
#>  [6] "Anahera"    "Oyunbileg"  "Sem"        "Tevita"     "Eloise"
```

## Educational Research Applications

### Self-Regulated Learning (SRL)

``` r

# Zimmerman's SRL phases
forethought <- get_learning_states("metacognitive")
performance <- get_learning_states(c("cognitive", "behavioral"))
reflection <- c("Evaluate", "Reflect", "Assess", "Judge")

cat("Forethought:", head(forethought, 5), "...\n")
#> Forethought: Plan Monitor Evaluate Reflect Regulate ...
cat("Performance:", head(performance, 5), "...\n")
#> Performance: Read Study Analyze Summarize Memorize ...
cat("Reflection:", reflection, "\n")
#> Reflection: Evaluate Reflect Assess Judge
```

### Socially Shared Regulation of Learning (SSRL)

``` r

# SSRL components
ssrl_verbs <- get_learning_states("group_regulation")
ssrl_verbs
#> [1] "Adapt"      "Cohesion"   "Consensus"  "Coregulate" "Discuss"   
#> [6] "Emotion"    "Monitor"    "Plan"       "Synthesis"
```

### Learning Analytics (LMS Data)

``` r

# Common LMS actions
lms_verbs <- get_learning_states("lms")
head(lms_verbs, 15)
#>  [1] "View"       "Access"     "Download"   "Upload"     "Submit"    
#>  [6] "Click"      "Navigate"   "Browse"     "Login"      "Logout"    
#> [11] "Post"       "Reply"      "Forum"      "Quiz"       "Assignment"
```

### Collaborative Learning

``` r

# Social and group regulation
collab_verbs <- get_learning_states(c("social", "group_regulation"))
length(collab_verbs)
#> [1] 38
head(collab_verbs, 10)
#>  [1] "Collaborate" "Discuss"     "Seek_help"   "Question"    "Explain"    
#>  [6] "Share"       "Teach"       "Tutor"       "Debate"      "Argue"
```

## Custom State Names

You can always use your own state names:

``` r

# With simulate_sequences
sequences <- simulate_sequences(
  n_sequences = 50,
  seq_length = 15,
  n_states = 4,
  states = c("Explore", "Learn", "Practice", "Master"),
  seed = 42
)

unique(unlist(sequences))
#> [1] "Plan"       "Regulate"   "Discourage" "Judge"

# With simulate_matrix
mat <- simulate_matrix(
  n_nodes = 4,
  names = c("Phase1", "Phase2", "Phase3", "Phase4"),
  seed = 42
)

rownames(mat)
#> [1] "Phase1" "Phase2" "Phase3" "Phase4"
```

## Complete Reference Table

``` r

# Create comprehensive reference
ref_table <- data.frame(
  Category = names(LEARNING_STATES),
  Count = sapply(LEARNING_STATES, length),
  First_5 = sapply(LEARNING_STATES, function(x) {
    paste(head(x, 5), collapse = ", ")
  }),
  row.names = NULL
)

ref_table
#>           Category Count                                            First_5
#> 1    metacognitive    20         Plan, Monitor, Evaluate, Reflect, Regulate
#> 2        cognitive    30          Read, Study, Analyze, Summarize, Memorize
#> 3       behavioral    30       Practice, Annotate, Research, Review, Revise
#> 4           social    30 Collaborate, Discuss, Seek_help, Question, Explain
#> 5     motivational    30            Focus, Persist, Explore, Create, Strive
#> 6        affective    30        Enjoy, Appreciate, Value, Interest, Curious
#> 7 group_regulation     9    Adapt, Cohesion, Consensus, Coregulate, Discuss
#> 8              lms    30             View, Access, Download, Upload, Submit
```

## Summary

### Key Functions

| Function | Purpose |
|----|----|
| `LEARNING_STATES` | Raw dataset of all verbs |
| [`get_learning_states()`](https://pak.dynasite.org/Saqrlab/reference/get_learning_states.md) | Retrieve verbs by category |
| [`list_learning_categories()`](https://pak.dynasite.org/Saqrlab/reference/list_learning_categories.md) | Show category summary |
| [`select_states()`](https://pak.dynasite.org/Saqrlab/reference/select_states.md) | Intelligent selection |
| `GLOBAL_NAMES` | 300 diverse names |
| [`get_global_names()`](https://pak.dynasite.org/Saqrlab/reference/get_global_names.md) | Retrieve names |

### Category Summary

| Category         | Count   | Focus               |
|------------------|---------|---------------------|
| metacognitive    | 20      | Self-regulation     |
| cognitive        | 30      | Mental processing   |
| behavioral       | 30      | Observable actions  |
| social           | 30      | Interaction         |
| motivational     | 30      | Effort/persistence  |
| affective        | 30      | Emotions            |
| group_regulation | 9       | SSRL                |
| lms              | 30      | System interactions |
| **Total**        | **209** |                     |

### Usage Tips

1.  **Default behavior**: Most functions use learning states
    automatically
2.  **Categories parameter**: Control which categories to sample from
3.  **Custom names**: Always supported via `states` or `names` parameter
4.  **Reproducibility**: Use `seed` parameter for consistent results
5.  **HTNA/MLNA**: Each type automatically gets different category
