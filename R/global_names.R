#' Global Names Dataset
#'
#' @description
#' A collection of ~1000 culturally diverse names from around the world,
#' organized by region for use in social network simulations.
#'
#' @format A named list of character vectors, with each element containing
#'   names from a specific region.
#'
#' @details
#' Regions are organized geographically and include:
#'
#' \strong{Europe:}
#' \itemize{
#'   \item \code{western_europe}: UK, Ireland, France, Germany, Netherlands, Belgium
#'   \item \code{eastern_europe}: Russia, Ukraine, Poland, Czech, Romania, Hungary, Bulgaria
#'   \item \code{southern_europe}: Spain, Italy, Portugal, Greece
#'   \item \code{nordic}: Sweden, Norway, Denmark, Finland, Iceland
#' }
#'
#' \strong{Middle East & Central Asia:}
#' \itemize{
#'   \item \code{arab}: Gulf states, Levant, Egypt
#'   \item \code{persian}: Iran, Afghanistan, Tajikistan
#'   \item \code{turkish}: Turkey, Azerbaijan
#'   \item \code{central_asia}: Kazakhstan, Uzbekistan, Kyrgyzstan, Turkmenistan
#' }
#'
#' \strong{Asia:}
#' \itemize{
#'   \item \code{south_asia}: India, Pakistan, Bangladesh, Sri Lanka, Nepal
#'   \item \code{east_asia}: China, Japan, Korea, Mongolia
#'   \item \code{southeast_asia}: Vietnam, Thailand, Indonesia, Philippines, Malaysia
#' }
#'
#' \strong{Africa:}
#' \itemize{
#'   \item \code{west_africa}: Nigeria, Ghana, Senegal, Mali
#'   \item \code{east_africa}: Kenya, Ethiopia, Tanzania, Somalia, Uganda
#'   \item \code{southern_africa}: South Africa, Zimbabwe, Botswana, Namibia
#'   \item \code{north_africa}: Morocco, Algeria, Tunisia, Libya
#'   \item \code{central_africa}: DRC, Cameroon, Congo
#' }
#'
#' \strong{Americas:}
#' \itemize{
#'   \item \code{north_america}: USA, Canada
#'   \item \code{latin_america}: Mexico, Central America, South America, Brazil
#'   \item \code{caribbean}: Jamaica, Haiti, Cuba, Trinidad, Dominican Republic
#'   \item \code{indigenous_americas}: Native American, First Nations, Mayan, Incan
#' }
#'
#' \strong{Oceania:}
#' \itemize{
#'   \item \code{australia_nz}: Australia, New Zealand
#'   \item \code{pacific_islands}: Hawaii, Samoa, Tonga, Fiji, PNG
#'   \item \code{indigenous_oceania}: Maori, Aboriginal Australian
#' }
#'
#' \strong{Other:}
#' \itemize{
#'   \item \code{jewish_hebrew}: Hebrew/Jewish names
#'   \item \code{armenian_georgian}: Armenia, Georgia
#' }
#'
#' @examples
#' # View all region names
#' names(GLOBAL_NAMES)
#'
#' # View names from a specific region
#' GLOBAL_NAMES$arab
#' GLOBAL_NAMES$east_asia
#'
#' # Count names per region
#' sapply(GLOBAL_NAMES, length)
#'
#' # Total names available
#' sum(sapply(GLOBAL_NAMES, length))
#'
#' @seealso \code{\link{get_global_names}} for sampling names,
#'   \code{\link{list_name_regions}} for exploring available regions.
#'
#' @export
GLOBAL_NAMES <- list(

  # ══════════════════════════════════════════════════════════════
  # EUROPE
  # ══════════════════════════════════════════════════════════════

  western_europe = c(
    # UK/Ireland
    "Emma", "Liam", "Chloe", "Jack", "Lily", "Oliver", "Grace", "Harry",
    "Aoife", "Sean", "Niamh", "Finn", "Saoirse", "Ciaran", "Roisin", "Declan",
    # France
    "Lea", "Hugo", "Manon", "Lucas", "Chloe", "Jules", "Camille", "Louis",
    "Eloise", "Antoine", "Amelie", "Theo", "Juliette", "Raphael", "Margot", "Adrien",
    # Germany/Austria
    "Max", "Mia", "Felix", "Lena", "Leon", "Anna", "Paul", "Clara",
    "Lukas", "Sophie", "Jonas", "Lea", "Elias", "Marie", "Noah", "Emilia",
    # Netherlands/Belgium
    "Daan", "Sanne", "Sem", "Julia", "Lotte", "Ruben", "Fien", "Bram",
    "Thijs", "Eva", "Lars", "Lisa", "Milan", "Noor", "Tim", "Fleur"
  ),

  eastern_europe = c(
    # Russia/Ukraine/Belarus
    "Ivan", "Olga", "Dmitri", "Natasha", "Alexei", "Katya", "Yuri", "Anya",
    "Boris", "Mila", "Sasha", "Vera", "Vlad", "Daria", "Oleg", "Irina",
    "Nikita", "Lena", "Andrei", "Tanya", "Maxim", "Yulia", "Artem", "Sveta",
    # Poland
    "Jakub", "Zuzanna", "Kacper", "Maja", "Antoni", "Hanna", "Piotr", "Ewa",
    "Michal", "Anna", "Tomek", "Kasia", "Bartek", "Ola", "Pawel", "Magda",
    # Czech/Slovakia
    "Jiri", "Petra", "Tomas", "Hana", "Marek", "Eva", "Milan", "Jana",
    "Ondrej", "Lucie", "Adam", "Tereza", "Filip", "Klara", "Jakub", "Monika",
    # Romania/Bulgaria/Hungary
    "Andrei", "Elena", "Mircea", "Raluca", "Bela", "Eszter", "Matyas", "Zsofia",
    "Ion", "Maria", "Vlad", "Ana", "Attila", "Kata", "Levente", "Reka"
  ),

  southern_europe = c(
    # Spain
    "Pablo", "Lucia", "Diego", "Sofia", "Alvaro", "Carmen", "Sergio", "Marta",
    "Javier", "Ana", "Carlos", "Laura", "Alejandro", "Elena", "Daniel", "Paula",
    "Raul", "Ines", "Adrian", "Alba", "Jorge", "Nuria", "Ivan", "Cristina",
    # Italy
    "Marco", "Giulia", "Luca", "Francesca", "Matteo", "Chiara", "Andrea", "Sara",
    "Alessandro", "Valentina", "Lorenzo", "Martina", "Davide", "Federica", "Giuseppe", "Elisa",
    # Portugal
    "Joao", "Maria", "Tiago", "Ines", "Diogo", "Beatriz", "Pedro", "Ana",
    "Miguel", "Mariana", "Rui", "Catarina", "Andre", "Sofia", "Bruno", "Rita",
    # Greece
    "Nikos", "Elena", "Dimitris", "Maria", "Kostas", "Sofia", "Yannis", "Anna",
    "Giorgos", "Eleni", "Christos", "Katerina", "Alexandros", "Ioanna", "Stavros", "Dimitra"
  ),

  nordic = c(
    # Sweden
    "Erik", "Astrid", "Oscar", "Maja", "Axel", "Ella", "Viktor", "Saga",
    "Gustav", "Wilma", "William", "Ebba", "Linus", "Agnes", "Filip", "Freja",
    # Norway
    "Lars", "Ingrid", "Sven", "Liv", "Olaf", "Nora", "Leif", "Thea",
    "Magnus", "Emma", "Henrik", "Sofie", "Andreas", "Ida", "Kristian", "Julie",
    # Denmark
    "Mikkel", "Freya", "Emil", "Clara", "Christian", "Emilie", "Noah", "Alma",
    "Oliver", "Freja", "Lucas", "Laura", "Frederik", "Anna", "Mathias", "Maja",
    # Finland
    "Onni", "Aino", "Eetu", "Venla", "Lauri", "Helmi", "Eero", "Siiri",
    "Matias", "Emilia", "Aleksi", "Sofia", "Veeti", "Aada", "Elias", "Ella",
    # Iceland
    "Ragnar", "Sigrid", "Bjorn", "Gudrun", "Eirik", "Hildur", "Odin", "Frida",
    "Thor", "Helga", "Gunnar", "Unnur", "Bjarni", "Kristin", "Arnar", "Margret"
  ),

  # ══════════════════════════════════════════════════════════════
  # MIDDLE EAST & CENTRAL ASIA
  # ══════════════════════════════════════════════════════════════

  arab = c(
    # Gulf States
    "Ali", "Fatima", "Omar", "Layla", "Yusuf", "Mariam", "Khalid", "Noor",
    "Sultan", "Sheikha", "Rashid", "Amal", "Hamad", "Moza", "Faisal", "Hessa",
    # Levant (Syria, Lebanon, Jordan, Palestine)
    "Hassan", "Amira", "Amir", "Zara", "Tariq", "Salma", "Faris", "Rania",
    "Karim", "Leila", "Walid", "Dina", "Nabil", "Huda", "Samir", "Lina",
    # Egypt
    "Ahmed", "Mona", "Mahmoud", "Heba", "Mostafa", "Yasmine", "Tarek", "Dalia",
    "Mohamed", "Noha", "Amr", "Aya", "Khaled", "Mai", "Hazem", "Nada",
    # Common Arabic
    "Ziad", "Rana", "Sami", "Reem", "Adel", "Noura", "Mazen", "Asma",
    "Jamal", "Dalal", "Bassam", "Aisha", "Rami", "Hana", "Tamer", "Sara"
  ),

  persian = c(
    # Iran
    "Cyrus", "Shirin", "Dara", "Parisa", "Reza", "Maryam", "Arman", "Soraya",
    "Babak", "Mina", "Kaveh", "Nazanin", "Omid", "Sara", "Farhad", "Azar",
    "Behzad", "Niloufar", "Kian", "Mahsa", "Arash", "Pari", "Dariush", "Roxana",
    "Amir", "Leyla", "Kamran", "Shadi", "Nima", "Azadeh", "Siavash", "Setareh",
    # Afghanistan
    "Ahmad", "Fatima", "Farid", "Zainab", "Hamid", "Maryam", "Jawad", "Hosna",
    # Tajikistan
    "Rustam", "Gulnora", "Firdavs", "Munis", "Jamshed", "Dilbar", "Firuz", "Sitora"
  ),

  turkish = c(
    # Turkey
    "Emre", "Elif", "Kaan", "Defne", "Cem", "Zeynep", "Deniz", "Aylin",
    "Burak", "Selin", "Mert", "Ceren", "Baris", "Ebru", "Arda", "Melis",
    "Kemal", "Leyla", "Ozan", "Nehir", "Alp", "Yasemin", "Umut", "Ezgi",
    "Onur", "Duygu", "Tolga", "Pinar", "Serkan", "Burcu", "Hakan", "Esra",
    # Azerbaijan
    "Eldar", "Gunel", "Tural", "Nigar", "Orhan", "Sevil", "Farid", "Aynur",
    "Rauf", "Gulnar", "Ilkin", "Samira", "Nurlan", "Leyla", "Vugar", "Arzu"
  ),

  central_asia = c(
    # Kazakhstan
    "Aibek", "Aida", "Nurlan", "Gulnara", "Bekzat", "Aigerim", "Ruslan", "Dana",
    "Nursultan", "Ainur", "Almas", "Dilnaz", "Miras", "Zhanna", "Saken", "Aliya",
    # Uzbekistan
    "Bobur", "Nilufar", "Sardor", "Malika", "Jasur", "Dilnoza", "Nodir", "Zarina",
    # Kyrgyzstan
    "Bakyt", "Madina", "Azamat", "Asel", "Damir", "Cholpon", "Timur", "Sabina",
    # Turkmenistan
    "Merdan", "Mahri", "Oraz", "Ogulgerek", "Serdar", "Aylar", "Dovlet", "Jennet",
    # Tajikistan (additional)
    "Bahrom", "Zebo", "Farkhod", "Nigina", "Suhrab", "Madina", "Daler", "Firuza"
  ),

  # ══════════════════════════════════════════════════════════════
  # SOUTH ASIA
  # ══════════════════════════════════════════════════════════════

  south_asia = c(
    # India - Hindi Belt
    "Arun", "Priya", "Raj", "Devi", "Amit", "Sita", "Ravi", "Maya",
    "Ajay", "Lata", "Vijay", "Anita", "Rohit", "Pooja", "Arjun", "Kavita",
    "Rahul", "Sunita", "Kiran", "Aditi", "Nikhil", "Shreya", "Vikram", "Neha",
    "Sanjay", "Meera", "Aditya", "Divya", "Varun", "Isha", "Dev", "Tara",
    # India - South (Tamil, Telugu, etc.)
    "Surya", "Lakshmi", "Karthik", "Deepa", "Ganesh", "Radha", "Suresh", "Vasuki",
    "Rajan", "Chitra", "Ashok", "Padma", "Senthil", "Bhavani", "Mohan", "Uma",
    # Pakistan
    "Imran", "Ayesha", "Bilal", "Sana", "Hamza", "Hira", "Zain", "Maryam",
    "Usman", "Fatima", "Ali", "Zoya", "Hasan", "Amina", "Faisal", "Rabia",
    # Bangladesh
    "Rafiq", "Fatema", "Rahim", "Nasreen", "Shahid", "Rupa", "Kamal", "Shanta",
    "Tanvir", "Nusrat", "Mahbub", "Tahera", "Iqbal", "Hasina", "Jamal", "Rima",
    # Sri Lanka
    "Nimal", "Kumari", "Sandun", "Dilani", "Chamara", "Nimali", "Lasith", "Thilini",
    # Nepal
    "Binod", "Sushma", "Rajan", "Kamala", "Sunil", "Gita", "Ramesh", "Sabita"
  ),

  # ══════════════════════════════════════════════════════════════
  # EAST ASIA
  # ══════════════════════════════════════════════════════════════

  east_asia = c(
    # China
    "Wei", "Lin", "Chen", "Mei", "Jing", "Xiao", "Ming", "Hua",
    "Feng", "Yan", "Jun", "Lei", "Bo", "Lan", "Tao", "Ying",
    "Hao", "Yue", "Zhen", "Xin", "Peng", "Qian", "Rui", "Shu",
    "Long", "Fang", "Hong", "Li", "Chao", "Na", "Gang", "Ling",
    # Japan
    "Yuki", "Hiro", "Akira", "Sakura", "Kenji", "Yuna", "Taro", "Haruki",
    "Koji", "Miki", "Ryu", "Emi", "Rin", "Sora", "Kaito", "Hana",
    "Yuto", "Aoi", "Sota", "Mei", "Ren", "Yui", "Haruto", "Koharu",
    "Takumi", "Himari", "Kenta", "Akari", "Daiki", "Riko", "Naoki", "Ayumi",
    # Korea
    "Minho", "Jisoo", "Joon", "Seo", "Tae", "Yuna", "Woo", "Minji",
    "Dae", "Eunji", "Sung", "Hana", "Hyun", "Jiwon", "Seok", "Nari",
    "Jiho", "Somin", "Junho", "Yuri", "Minsoo", "Soyeon", "Donghyun", "Chaeyoung",
    # Mongolia
    "Bataar", "Oyun", "Bold", "Altai", "Temuulen", "Sarnai", "Ganzorig", "Enkhtuya",
    "Erdene", "Narantsetseg", "Munkh", "Tsetseg", "Baatar", "Oyunbileg", "Chuluun", "Bolormaa"
  ),

  # ══════════════════════════════════════════════════════════════
  # SOUTHEAST ASIA
  # ══════════════════════════════════════════════════════════════

  southeast_asia = c(
    # Vietnam
    "Linh", "Minh", "Anh", "Bao", "Mai", "Duc", "Thao", "Tuan",
    "Hoa", "Nam", "Lan", "Hung", "Hanh", "Cuong", "Ngoc", "Phong",
    "Trang", "Huy", "Nhi", "Khanh", "Vy", "Long", "Huong", "Dung",
    # Thailand
    "Somchai", "Ploy", "Chai", "Nong", "Kiet", "Mali", "Somsak", "Araya",
    "Prem", "Dao", "Niran", "Malai", "Sakchai", "Achara", "Prasert", "Supaporn",
    # Indonesia
    "Rizal", "Putri", "Bagus", "Sari", "Wayan", "Ayu", "Gede", "Nisa",
    "Dewi", "Budi", "Rani", "Eko", "Sinta", "Agus", "Ratna", "Dimas",
    "Adi", "Fitri", "Rudi", "Wulan", "Hendra", "Lestari", "Yusuf", "Mega",
    # Philippines
    "Jose", "Maria", "Juan", "Rosa", "Carlo", "Ana", "Miguel", "Luz",
    "Antonio", "Carmen", "Ramon", "Liza", "Francis", "Grace", "Paolo", "Joy",
    # Malaysia/Singapore
    "Fajar", "Indah", "Hafiz", "Nurul", "Amir", "Siti", "Ismail", "Aishah",
    "Zul", "Farah", "Azlan", "Aisyah", "Firdaus", "Nur", "Hakim", "Amira"
  ),

  # ══════════════════════════════════════════════════════════════
  # AFRICA
  # ══════════════════════════════════════════════════════════════

  west_africa = c(
    # Nigeria - Igbo
    "Chidi", "Adaora", "Emeka", "Ngozi", "Obi", "Amaka", "Chijioke", "Nneka",
    "Ikenna", "Chinwe", "Obinna", "Uchenna", "Nnamdi", "Chinyere", "Kelechi", "Adaobi",
    # Nigeria - Yoruba
    "Ade", "Yemi", "Tunde", "Funke", "Segun", "Bola", "Kunle", "Nike",
    "Dayo", "Sade", "Kayode", "Toyin", "Femi", "Lola", "Bayo", "Ronke",
    # Ghana
    "Kofi", "Ama", "Kwame", "Akua", "Yaw", "Abena", "Kwesi", "Efua",
    "Kojo", "Akosua", "Kwadwo", "Adwoa", "Kweku", "Afua", "Kobi", "Afia",
    # Senegal/Mali
    "Amadou", "Fatou", "Moussa", "Aissatou", "Oumar", "Mariama", "Sekou", "Binta",
    "Ibrahima", "Aminata", "Mamadou", "Kadiatou", "Boubacar", "Fatoumata", "Cheikh", "Oumou"
  ),

  east_africa = c(
    # Kenya
    "Juma", "Amina", "Bakari", "Zuri", "Mwangi", "Wanjiku", "Kamau", "Njeri",
    "Ochieng", "Auma", "Otieno", "Adhiambo", "Kipchoge", "Chebet", "Korir", "Jepkosgei",
    # Tanzania
    "Baraka", "Neema", "Hamisi", "Saida", "Rajabu", "Rehema", "Salum", "Mwajuma",
    # Ethiopia
    "Abebe", "Tigist", "Haile", "Makeda", "Tadesse", "Meron", "Dawit", "Sara",
    "Yohannes", "Bethlehem", "Girma", "Selam", "Tesfaye", "Hana", "Solomon", "Rahel",
    # Somalia
    "Mohamed", "Halima", "Abdi", "Amina", "Hassan", "Fartun", "Omar", "Sahra",
    # Uganda/Rwanda
    "Mugisha", "Uwimana", "Kato", "Akello", "Okello", "Nkunda", "Kagame", "Umutesi"
  ),

  southern_africa = c(
    # South Africa - Zulu/Xhosa
    "Themba", "Lindiwe", "Sipho", "Nomzamo", "Thabo", "Zanele", "Mandla", "Lerato",
    "Sibusiso", "Thandiwe", "Bongani", "Nolwazi", "Siyabonga", "Ayanda", "Luyanda", "Nomvula",
    # South Africa - Sotho/Tswana
    "Mpho", "Kelebogile", "Neo", "Goitseone", "Tshepo", "Dineo", "Thato", "Palesa",
    # Zimbabwe
    "Tendai", "Rudo", "Tapiwa", "Chipo", "Tinashe", "Tariro", "Farai", "Nyasha",
    "Tawanda", "Rutendo", "Tatenda", "Rumbidzai", "Kudakwashe", "Ropafadzo", "Tanaka", "Rufaro",
    # Zambia/Botswana/Namibia
    "Mulenga", "Mutinta", "Chisomo", "Thandie", "Kgosi", "Naledi", "Amantle", "Masego"
  ),

  north_africa = c(
    # Morocco
    "Youssef", "Fatima", "Hassan", "Amina", "Said", "Khadija", "Rachid", "Zahra",
    "Mehdi", "Salma", "Amine", "Houda", "Yassine", "Imane", "Hamza", "Hajar",
    # Algeria
    "Malik", "Samira", "Kamel", "Nadia", "Hamid", "Leila", "Sofiane", "Salima",
    "Riad", "Yasmina", "Farid", "Lamia", "Nassim", "Amira", "Djamel", "Sihem",
    # Tunisia
    "Hedi", "Sonia", "Nizar", "Amel", "Fares", "Meriem", "Slim", "Ines",
    # Libya
    "Mukhtar", "Hana", "Idris", "Asma", "Fathi", "Marwa", "Nuri", "Salwa"
  ),

  central_africa = c(
    # DRC
    "Patrice", "Solange", "Jean", "Marie", "Emmanuel", "Grace", "David", "Ruth",
    "Fiston", "Carine", "Christian", "Nadine", "Patrick", "Sylvie", "Serge", "Claudine",
    # Cameroon
    "Samuel", "Esther", "Paul", "Christiane", "Pierre", "Nadege", "Andre", "Simone",
    "Yves", "Blanche", "Roger", "Jeanne", "Fabrice", "Mireille", "Herve", "Colette",
    # Congo/Gabon
    "Arsene", "Prudence", "Guy", "Diane", "Didier", "Lydie", "Landry", "Ornella"
  ),

  # ══════════════════════════════════════════════════════════════
  # AMERICAS
  # ══════════════════════════════════════════════════════════════

  north_america = c(
    # USA - Common names
    "James", "Emily", "Michael", "Sarah", "David", "Jessica", "John", "Ashley",
    "Robert", "Amanda", "William", "Jennifer", "Chris", "Nicole", "Matt", "Lauren",
    "Ryan", "Rachel", "Tyler", "Megan", "Brian", "Stephanie", "Kevin", "Brittany",
    "Andrew", "Samantha", "Justin", "Elizabeth", "Brandon", "Heather", "Josh", "Michelle",
    # Canada
    "Liam", "Emma", "Noah", "Olivia", "Ethan", "Ava", "Mason", "Sophia",
    "Logan", "Charlotte", "Jacob", "Amelia", "Lucas", "Harper", "Jack", "Evelyn"
  ),

  latin_america = c(
    # Mexico/Central America
    "Juan", "Rosa", "Diego", "Luz", "Pablo", "Sol", "Luis", "Ana",
    "Carlos", "Lucia", "Miguel", "Elena", "Jose", "Sofia", "Pedro", "Camila",
    "Alejandro", "Valentina", "Fernando", "Daniela", "Ricardo", "Mariana", "Eduardo", "Fernanda",
    # South America - Spanish
    "Mateo", "Isabella", "Rafael", "Isabel", "Sergio", "Carmen", "Alvaro", "Paulina",
    "Andres", "Maria", "Santiago", "Martina", "Gabriel", "Florencia", "Nicolas", "Julieta",
    "Sebastian", "Catalina", "Benjamin", "Antonia", "Martin", "Emilia", "Tomas", "Isidora",
    # Brazil
    "Joao", "Julia", "Pedro", "Fernanda", "Lucas", "Amanda", "Gustavo", "Beatriz",
    "Rafael", "Larissa", "Thiago", "Bruna", "Matheus", "Carolina", "Vitor", "Leticia",
    "Bruno", "Mariana", "Felipe", "Gabriela", "Leonardo", "Camila", "Gabriel", "Isabela"
  ),

  caribbean = c(
    # Jamaica/Trinidad
    "Marlon", "Keisha", "Dwayne", "Shanique", "Andre", "Natalie", "Wayne", "Tanya",
    "Usain", "Shelly", "Damian", "Khadija", "Tristan", "Yolanda", "Jermaine", "Sasha",
    # Haiti
    "Jean", "Marie", "Pierre", "Rose", "Jacques", "Nadine", "Claude", "Carole",
    "Wyclef", "Michaelle", "Stanley", "Fabienne", "Herby", "Guerda", "Frantz", "Ketty",
    # Cuba/Dominican Republic
    "Alejandro", "Yolanda", "Orlando", "Marisol", "Ramon", "Yesenia", "Jorge", "Luz",
    "Lazaro", "Yanelis", "Yunel", "Dania", "Pedro", "Yuliesky", "Miguel", "Yamilet"
  ),

  indigenous_americas = c(
    # Native American
    "Takoda", "Winona", "Koda", "Aiyana", "Chayton", "Aponi", "Ahanu", "Chenoa",
    "Hinto", "Halona", "Shilah", "Mika", "Nayeli", "Kimi", "Huritt", "Ayita",
    "Sequoia", "Cochise", "Sacagawea", "Tecumseh", "Hiawatha", "Pocahontas", "Dakota", "Cheyenne",
    # First Nations Canada
    "Taima", "Kaya", "Nanuq", "Sedna", "Amaruq", "Siku", "Tulok", "Atka",
    # Mayan/Aztec
    "Itzamna", "Xochitl", "Cuauhtemoc", "Citlali", "Quetzal", "Itzel", "Tlaloc", "Ixchel",
    # Incan/Quechua
    "Tupac", "Qori", "Inti", "Sumaq", "Amaru", "Killa", "Rumi", "Wayra"
  ),

  # ══════════════════════════════════════════════════════════════
  # OCEANIA & PACIFIC
  # ══════════════════════════════════════════════════════════════

  australia_nz = c(
    # Australia - Common
    "Jack", "Charlotte", "Oliver", "Amelia", "Noah", "Olivia", "William", "Isla",
    "Cooper", "Ava", "Lachlan", "Mia", "Riley", "Grace", "Flynn", "Sophie",
    "Harrison", "Chloe", "Hunter", "Emily", "Archie", "Zoe", "Oscar", "Lily",
    # New Zealand - Common
    "Jack", "Charlotte", "Oliver", "Olivia", "Leo", "Isla", "Liam", "Amelia",
    "James", "Ella", "Mason", "Mila", "Hunter", "Ava", "Finn", "Harper"
  ),

  pacific_islands = c(
    # Hawaii
    "Kai", "Leilani", "Keanu", "Mahina", "Koa", "Malia", "Manu", "Moana",
    "Keoni", "Nalani", "Makoa", "Kalani", "Kawika", "Alana", "Ikaika", "Kailani",
    # Samoa
    "Tane", "Sina", "Sione", "Mele", "Pita", "Ana", "Losa", "Teuila",
    "Manu", "Lagi", "Tavita", "Sieni", "Ioane", "Sala", "Pule", "Naea",
    # Tonga/Fiji
    "Sione", "Mele", "Tevita", "Ana", "Viliami", "Salote", "Mosese", "Sera",
    "Jone", "Mere", "Waisea", "Adi", "Samu", "Vasiti", "Ratu", "Litia"
  ),

  indigenous_oceania = c(
    # Maori (New Zealand)
    "Ariki", "Aroha", "Tui", "Hemi", "Wiremu", "Mere", "Rawiri", "Anahera",
    "Nikau", "Manaia", "Tama", "Hine", "Kauri", "Ngaio", "Rangi", "Marama",
    "Tamati", "Kiri", "Ihaia", "Parehuia", "Matiu", "Anika", "Hoani", "Ripeka",
    # Aboriginal Australian
    "Jarrah", "Kirra", "Bindi", "Marley", "Jinda", "Allira", "Koora", "Talia",
    "Yarran", "Alinta", "Jedda", "Kalinda", "Mirri", "Lowanna", "Burnum", "Merindah",
    # Torres Strait
    "Kala", "Seri", "Waia", "Gizu", "Maino", "Nazareth", "Elia", "Cessa"
  ),

  # ══════════════════════════════════════════════════════════════
  # OTHER
  # ══════════════════════════════════════════════════════════════

  jewish_hebrew = c(
    # Hebrew/Jewish names
    "David", "Sarah", "Daniel", "Rachel", "Aaron", "Miriam", "Noah", "Hannah",
    "Eli", "Leah", "Isaac", "Rebecca", "Jacob", "Esther", "Benjamin", "Ruth",
    "Avi", "Shira", "Yosef", "Noa", "Moshe", "Tamar", "Ezra", "Naomi",
    "Gideon", "Aviva", "Asher", "Maya", "Levi", "Talia", "Ethan", "Rivka",
    "Oren", "Yael", "Gilad", "Liora", "Noam", "Shoshana", "Eitan", "Adina"
  ),

  armenian_georgian = c(
    # Armenia
    "Armen", "Ani", "Hayk", "Lilit", "Aram", "Nare", "Tigran", "Sona",
    "Levon", "Anahit", "Artur", "Gohar", "Vardan", "Mariam", "Gagik", "Lusine",
    # Georgia
    "Giorgi", "Nino", "Luka", "Mariam", "Nika", "Ana", "Davit", "Tamara",
    "Sandro", "Salome", "Giga", "Maka", "Irakli", "Eka", "Lado", "Keti"
  )
)


#' Region Shortcuts for GLOBAL_NAMES
#'
#' @description
#' Internal mapping of shortcut names to actual region names in GLOBAL_NAMES.
#'
#' @keywords internal
#' @noRd
.REGION_SHORTCUTS <- list(
  # Continent shortcuts
  europe = c("western_europe", "eastern_europe", "southern_europe", "nordic"),
  middle_east = c("arab", "persian", "turkish", "central_asia"),
  asia = c("south_asia", "east_asia", "southeast_asia"),
  africa = c("west_africa", "east_africa", "southern_africa", "north_africa", "central_africa"),
  americas = c("north_america", "latin_america", "caribbean", "indigenous_americas"),
  oceania = c("australia_nz", "pacific_islands", "indigenous_oceania"),

  # Sub-region shortcuts
  western = c("western_europe", "north_america"),
  eastern = c("eastern_europe", "east_asia"),
  nordic_baltic = c("nordic"),
  mediterranean = c("southern_europe", "north_africa", "arab"),
  sub_saharan = c("west_africa", "east_africa", "southern_africa", "central_africa"),
  pacific = c("pacific_islands", "indigenous_oceania"),
  indigenous = c("indigenous_americas", "indigenous_oceania")
)


#' List Available Name Regions
#'
#' @description
#' Display all available regions and shortcuts for selecting names from GLOBAL_NAMES.
#'
#' @return A data frame with region names, counts, and example names.
#'
#' @examples
#' list_name_regions()
#'
#' @export
list_name_regions <- function() {
  regions <- names(GLOBAL_NAMES)
  counts <- sapply(GLOBAL_NAMES, length)
  examples <- sapply(GLOBAL_NAMES, function(x) {
    paste(head(x, 4), collapse = ", ")
  })

  df <- data.frame(
    Region = regions,
    Count = counts,
    Examples = paste0(examples, ", ..."),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # Add total row
  total_row <- data.frame(
    Region = "TOTAL",
    Count = sum(counts),
    Examples = "",
    stringsAsFactors = FALSE
  )

  df <- rbind(df, total_row)

  # Print shortcuts info
  cat("Available region shortcuts:\n")
  cat("  europe, middle_east, asia, africa, americas, oceania\n")
  cat("  western, eastern, mediterranean, sub_saharan, pacific, indigenous\n\n")

  return(df)
}


#' Get Global Names
#'
#' @description
#' Retrieve names from the global names dataset, with optional filtering by region.
#'
#' @param n Integer or NULL. Number of names to return. If NULL, returns all
#'   matching names. Default: NULL.
#' @param regions Character vector. Regions to sample from. Can be:
#'   \itemize{
#'     \item \code{"all"}: All regions (default)
#'     \item Specific regions: e.g., \code{"arab"}, \code{"east_asia"}, \code{"nordic"}
#'     \item Shortcuts: e.g., \code{"europe"}, \code{"africa"}, \code{"asia"}
#'     \item Multiple: e.g., \code{c("arab", "persian", "turkish")}
#'   }
#'   Default: "all".
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#'
#' @return Character vector of names.
#'
#' @examples
#' # Get 20 random names from all regions
#' get_global_names(20, seed = 42)
#'
#' # Get 10 names from Arab region
#' get_global_names(10, regions = "arab")
#'
#' # Get 15 names from any African region
#' get_global_names(15, regions = "africa")
#'
#' # Get names from multiple regions
#' get_global_names(20, regions = c("east_asia", "south_asia"))
#'
#' # Get all Nordic names
#' get_global_names(regions = "nordic")
#'
#' @seealso \code{\link{list_name_regions}} for available regions,
#'   \code{\link{GLOBAL_NAMES}} for the full dataset.
#'
#' @export
get_global_names <- function(n = NULL, regions = "all", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

 # Expand regions
  if (length(regions) == 1 && regions == "all") {
    selected_regions <- names(GLOBAL_NAMES)
  } else {
    selected_regions <- character(0)
    for (r in regions) {
      if (r %in% names(.REGION_SHORTCUTS)) {
        # It's a shortcut - expand it
        selected_regions <- c(selected_regions, .REGION_SHORTCUTS[[r]])
      } else if (r %in% names(GLOBAL_NAMES)) {
        # It's a direct region name
        selected_regions <- c(selected_regions, r)
      } else {
        warning(sprintf("Unknown region: '%s'. Use list_name_regions() to see available regions.", r))
      }
    }
    selected_regions <- unique(selected_regions)
  }

  # Collect names from selected regions
  all_names <- unlist(GLOBAL_NAMES[selected_regions], use.names = FALSE)
  all_names <- unique(all_names)

  if (is.null(n)) {
    return(all_names)
  }

  if (n > length(all_names)) {
    warning(sprintf(
      "Requested %d names but only %d available in selected regions. Returning all.",
      n, length(all_names)
    ))
    return(sample(all_names))
  }

  sample(all_names, n)
}
