library(tibble)
library(dplyr)

# slovenian fishes

slovenian_fish <- tribble(
  ~species, ~family, ~order, ~class,
  # Petromyzontida
  "Petromyzon marinus", "Petromyzontidae", "Petromyzontiformes", "Petromyzontida",
  "Eudontomyzon vladykovi", "Petromyzontidae", "Petromyzontiformes", "Petromyzontida",
  "Lampetra planeri", "Petromyzontidae", "Petromyzontiformes", "Petromyzontida",
  
  # Acipenseriformes
  "Acipenser ruthenus", "Acipenseridae", "Acipenseriformes", "Actinopterygii",
  "Acipenser baerii", "Acipenseridae", "Acipenseriformes", "Actinopterygii",
  "Polyodon spathula", "Polyodontidae", "Acipenseriformes", "Actinopterygii",
  
  # Anguilliformes
  "Anguilla anguilla", "Anguillidae", "Anguilliformes", "Actinopterygii",
  
  # Clupeiformes
  "Alosa immaculata", "Clupeidae", "Clupeiformes", "Actinopterygii",
  
  # Cypriniformes: Cyprinidae
  "Mylopharyngodon piceus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Rutilus rutilus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Rutilus virgo", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Rutilus aula", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Leuciscus leuciscus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Squalius cephalus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Squalius squalus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Squalius janae", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Leuciscus idus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Alburnus alburnus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Alburnus arborella", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Alburnus sarmaticus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Phoxinus lumaireul", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Phoxinus phoxinus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Ctenopharyngodon idella", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Scardinius erythrophthalmus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Leuciscus aspius", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Tinca tinca", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Chondrostoma nasus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Protochondrostoma genei", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Chondrostoma soetta", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Gobio obtusirostris", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Romanogobio vladykovi", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Romanogobio kesslerii", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Romanogobio uranoscopus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Romanogobio benacensis", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Pseudorasbora parva", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Barbus barbus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Barbus plebejus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Barbus balcanicus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Barbus caninus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Rhodeus amarus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Vimba vimba", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Abramis brama", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Ballerus sapa", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Ballerus ballerus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Blicca bjoerkna", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Alburnoides bipunctatus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Pelecus cultratus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Carassius carassius", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Carassius auratus", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Carassius gibelio", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Cyprinus carpio", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Hypophthalmichthys molitrix", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  "Hypophthalmichthys nobilis", "Cyprinidae", "Cypriniformes", "Actinopterygii",
  
  # Cypriniformes: Cobitidae
  "Cobitis elongatoides", "Cobitidae", "Cypriniformes", "Actinopterygii",
  "Cobitis bilineata", "Cobitidae", "Cypriniformes", "Actinopterygii",
  "Cobitis elongata", "Cobitidae", "Cypriniformes", "Actinopterygii",
  "Sabanejewia balcanica", "Cobitidae", "Cypriniformes", "Actinopterygii",
  "Misgurnus fossilis", "Cobitidae", "Cypriniformes", "Actinopterygii",
  
  # Cypriniformes: Nemacheilidae
  "Barbatula barbatula", "Nemacheilidae", "Cypriniformes", "Actinopterygii",
  
  # Siluriformes
  "Silurus glanis", "Siluridae", "Siluriformes", "Actinopterygii",
  "Ameiurus nebulosus", "Ictaluridae", "Siluriformes", "Actinopterygii",
  "Ameiurus melas", "Ictaluridae", "Siluriformes", "Actinopterygii",
  "Clarias gariepinus", "Clariidae", "Siluriformes", "Actinopterygii",
  
  # Esociformes
  "Esox lucius", "Esocidae", "Esociformes", "Actinopterygii",
  "Umbra krameri", "Umbridae", "Esociformes", "Actinopterygii",
  
  # Salmoniformes
  "Salmo trutta", "Salmonidae", "Salmoniformes", "Actinopterygii",
  "Salmo marmoratus", "Salmonidae", "Salmoniformes", "Actinopterygii",
  "Oncorhynchus mykiss", "Salmonidae", "Salmoniformes", "Actinopterygii",
  "Oncorhynchus kisutch", "Salmonidae", "Salmoniformes", "Actinopterygii",
  "Salvelinus fontinalis", "Salmonidae", "Salmoniformes", "Actinopterygii",
  "Salvelinus umbla", "Salmonidae", "Salmoniformes", "Actinopterygii",
  "Hucho hucho", "Salmonidae", "Salmoniformes", "Actinopterygii",
  "Coregonus lavaretus", "Coregonidae", "Salmoniformes", "Actinopterygii",
  "Thymallus thymallus", "Thymallidae", "Salmoniformes", "Actinopterygii",
  "Thymallus arcticus", "Thymallidae", "Salmoniformes", "Actinopterygii",
  
  # Gadiformes
  "Lota lota", "Lotidae", "Gadiformes", "Actinopterygii",
  
  # Cyprinodontiformes
  "Gambusia holbrooki", "Poeciliidae", "Cyprinodontiformes", "Actinopterygii",
  "Aphanius fasciatus", "Cyprinodontidae", "Cyprinodontiformes", "Actinopterygii",
  
  # Gasterosteiformes
  "Gasterosteus gymnurus", "Gasterosteidae", "Gasterosteiformes", "Actinopterygii",
  "Gasterosteus aculeatus", "Gasterosteidae", "Gasterosteiformes", "Actinopterygii",
  
  # Scorpaeniformes
  "Cottus gobio", "Cottidae", "Scorpaeniformes", "Actinopterygii",
  "Cottus metae", "Cottidae", "Scorpaeniformes", "Actinopterygii",
  
  # Perciformes
  "Padogobius bonelli", "Gobiidae", "Perciformes", "Actinopterygii",
  "Knipowitschia panizzae", "Gobiidae", "Perciformes", "Actinopterygii",
  "Ponticola kessleri", "Gobiidae", "Perciformes", "Actinopterygii",
  "Neogobius fluviatilis", "Gobiidae", "Perciformes", "Actinopterygii",
  "Neogobius melanostomus", "Gobiidae", "Perciformes", "Actinopterygii",
  "Babka gymnotrachelus", "Gobiidae", "Perciformes", "Actinopterygii",
  "Perca fluviatilis", "Percidae", "Perciformes", "Actinopterygii",
  "Sander lucioperca", "Percidae", "Perciformes", "Actinopterygii",
  "Gymnocephalus cernuus", "Percidae", "Perciformes", "Actinopterygii",
  "Gymnocephalus baloni", "Percidae", "Perciformes", "Actinopterygii",
  "Gymnocephalus schraetser", "Percidae", "Perciformes", "Actinopterygii",
  "Zingel streber", "Percidae", "Perciformes", "Actinopterygii",
  "Zingel zingel", "Percidae", "Perciformes", "Actinopterygii",
  "Micropterus salmoides", "Centrarchidae", "Perciformes", "Actinopterygii",
  "Lepomis gibbosus", "Centrarchidae", "Perciformes", "Actinopterygii",
  "Oreochromis niloticus", "Cichlidae", "Perciformes", "Actinopterygii"
) %>%
  mutate(genus = sub(" .*", "", species)) %>%
  select(species, genus, family, order, class)

# yangtze reptiles

yz_reptiles <- tribble(
  ~species, ~family, ~order,
  # Testudiformes
  "Platysternon megacephalum", "Platysternidae", "Testudiformes",
  "Chinemys megalocephala", "Emydidae", "Testudiformes",
  "Chinemys reevesii", "Emydidae", "Testudiformes",
  "Cuora aurocapitata", "Emydidae", "Testudiformes",
  "Cuora flavomarginata", "Emydidae", "Testudiformes",
  "Cuora pani", "Emydidae", "Testudiformes",
  "Cuora yunnanensis", "Emydidae", "Testudiformes",
  "Geoemyda spengleri", "Emydidae", "Testudiformes",
  "Mauremys mutica", "Emydidae", "Testudiformes",
  "Ocadia sinensis", "Emydidae", "Testudiformes",
  "Sacalia bealei", "Emydidae", "Testudiformes",
  "Sacalia quadriocellata", "Emydidae", "Testudiformes",
  "Manouria impressa", "Testudinidae", "Testudiformes",
  "Caretta caretta", "Cheloniidae", "Testudiformes",
  "Palea steindachneri", "Trionychidae", "Testudiformes",
  "Pelochelys bibroni", "Trionychidae", "Testudiformes",
  "Pelochelys maculatus", "Trionychidae", "Testudiformes",
  "Pelodiscus sinensis", "Trionychidae", "Testudiformes",
  
  # Crocodiliformes
  "Alligator sinensis", "Alligatoridae", "Crocodiliformes",
  
  # Squamata: Agamidae
  "Acanthosaura lepidogaster", "Agamidae", "Squamata",
  "Calotes emma", "Agamidae", "Squamata",
  "Calotes versicolor", "Agamidae", "Squamata",
  "Japalura dymondi", "Agamidae", "Squamata",
  "Japalura flaviceps", "Agamidae", "Squamata",
  "Japalura grahami", "Agamidae", "Squamata",
  "Japalura micangshanensis", "Agamidae", "Squamata",
  "Japalura splendida", "Agamidae", "Squamata",
  "Japalura zhechuanensis", "Agamidae", "Squamata",
  "Japalura varcoae", "Agamidae", "Squamata",
  "Japalura yunnanensis", "Agamidae", "Squamata",
  "Phrynocephalus vlangalii", "Agamidae", "Squamata",
  
  # Squamata: Gekkonidae
  "Gekko chinensis", "Gekkonidae", "Squamata",
  "Gekko hokouensis", "Gekkonidae", "Squamata",
  "Gekko japonicus", "Gekkonidae", "Squamata",
  "Gekko scabridus", "Gekkonidae", "Squamata",
  "Gekko subpalmatus", "Gekkonidae", "Squamata",
  "Gekko taibaiensis", "Gekkonidae", "Squamata",
  "Hemidactylus bowringii", "Gekkonidae", "Squamata",
  "Hemiphyllodactylus yunnanensis", "Gekkonidae", "Squamata",
  
  # Squamata: Scincidae
  "Ateuchosaurus chinensis", "Scincidae", "Squamata",
  "Eumeces capito", "Scincidae", "Squamata",
  "Eumeces chinensis", "Scincidae", "Squamata",
  "Eumeces elegans", "Scincidae", "Squamata",
  "Eumeces liui", "Scincidae", "Squamata",
  "Eumeces tunganus", "Scincidae", "Squamata",
  "Scincella barbouri", "Scincidae", "Squamata",
  "Scincella doriae", "Scincidae", "Squamata",
  "Scincella modesta", "Scincidae", "Squamata",
  "Scincella monticola", "Scincidae", "Squamata",
  "Scincella potanini", "Scincidae", "Squamata",
  "Scincella reevesii", "Scincidae", "Squamata",
  "Scincella schmidti", "Scincidae", "Squamata",
  "Scincella tsinlingensis", "Scincidae", "Squamata",
  "Sphenomorphus incognitus", "Scincidae", "Squamata",
  "Sphenomorphus indicus", "Scincidae", "Squamata",
  "Tropidophorus hainanus", "Scincidae", "Squamata",
  
  # Squamata: Lacertidae
  "Eremias argus", "Lacertidae", "Squamata",
  "Platyplacopus intermedius", "Lacertidae", "Squamata",
  "Platyplacopus kuehnei", "Lacertidae", "Squamata",
  "Takydromus septentrionalis", "Lacertidae", "Squamata",
  "Takydromus sexlineatus", "Lacertidae", "Squamata",
  "Takydromus wolteri", "Lacertidae", "Squamata",
  
  # Squamata: Dibamidae
  "Dibamus bourreti", "Dibamidae", "Squamata",
  
  # Squamata: Anguidae
  "Ophisaurus gracilis", "Anguidae", "Squamata",
  "Ophisaurus harti", "Anguidae", "Squamata",
  
  # Squamata: Typhlopidae
  "Ramphotyphlops braminus", "Typhlopidae", "Squamata",
  
  # Squamata: Boidae
  "Python molurus", "Boidae", "Squamata",
  
  # Squamata: Xenopeltidae
  "Xenopeltis hainanensis", "Xenopeltidae", "Squamata",
  
  # Squamata: Colubridae
  "Achalinus ater", "Colubridae", "Squamata",
  "Achalinus jinggangensis", "Colubridae", "Squamata",
  "Achalinus meiguensis", "Colubridae", "Squamata",
  "Achalinus rufescens", "Colubridae", "Squamata",
  "Achalinus spinalis", "Colubridae", "Squamata",
  "Ahaetulla prasina", "Colubridae", "Squamata",
  "Amphiesma atemporalis", "Colubridae", "Squamata",
  "Amphiesma boulengeri", "Colubridae", "Squamata",
  "Amphiesma craspedogaster", "Colubridae", "Squamata",
  "Amphiesma johannis", "Colubridae", "Squamata",
  "Amphiesma metusia", "Colubridae", "Squamata",
  "Amphiesma modesta", "Colubridae", "Squamata",
  "Amphiesma octolineata", "Colubridae", "Squamata",
  "Amphiesma optata", "Colubridae", "Squamata",
  "Amphiesma popei", "Colubridae", "Squamata",
  "Amphiesma sauteri", "Colubridae", "Squamata",
  "Amphiesma stolata", "Colubridae", "Squamata",
  "Boiga kraepelini", "Colubridae", "Squamata",
  "Boiga multomaculata", "Colubridae", "Squamata",
  "Calamaria pavimentata", "Colubridae", "Squamata",
  "Calamaria septentrionalis", "Colubridae", "Squamata",
  "Coluber spinalis", "Colubridae", "Squamata",
  "Cyclophiops major", "Colubridae", "Squamata",
  "Dinodon flavozonatum", "Colubridae", "Squamata",
  "Dinodon rufozonatum", "Colubridae", "Squamata",
  "Elaphe anomala", "Colubridae", "Squamata",
  "Elaphe bimaculata", "Colubridae", "Squamata",
  "Elaphe carinata", "Colubridae", "Squamata",
  "Elaphe dione", "Colubridae", "Squamata",
  "Elaphe frenata", "Colubridae", "Squamata",
  "Elaphe mandarina", "Colubridae", "Squamata",
  "Elaphe perlacea", "Colubridae", "Squamata",
  "Elaphe porphyracea", "Colubridae", "Squamata",
  "Elaphe prasina", "Colubridae", "Squamata",
  "Elaphe rufodorsata", "Colubridae", "Squamata",
  "Elaphe taeniura", "Colubridae", "Squamata",
  "Enhydris chinensis", "Colubridae", "Squamata",
  "Enhydris plumbea", "Colubridae", "Squamata",
  "Lycodon fasciatus", "Colubridae", "Squamata",
  "Lycodon ruhstrati", "Colubridae", "Squamata",
  "Macropisthodon rudis", "Colubridae", "Squamata",
  "Oligodon chinensis", "Colubridae", "Squamata",
  "Oligodon cinereus", "Colubridae", "Squamata",
  "Oligodon formosanus", "Colubridae", "Squamata",
  "Oligodon kunmingensis", "Colubridae", "Squamata",
  "Oligodon lacroixi", "Colubridae", "Squamata",
  "Oligodon lungshenensis", "Colubridae", "Squamata",
  "Oligodon multizonatum", "Colubridae", "Squamata",
  "Oligodon ningshanensis", "Colubridae", "Squamata",
  "Oligodon ornatus", "Colubridae", "Squamata",
  "Opisthotropis balteata", "Colubridae", "Squamata",
  "Opisthotropis kautunensis", "Colubridae", "Squamata",
  "Opisthotropis lateralis", "Colubridae", "Squamata",
  "Opisthotropis latouchii", "Colubridae", "Squamata",
  "Pareas boulengeri", "Colubridae", "Squamata",
  "Pareas chinensis", "Colubridae", "Squamata",
  "Pareas hamptoni", "Colubridae", "Squamata",
  "Pareas stanleyi", "Colubridae", "Squamata",
  "Plagiopholis blakewayi", "Colubridae", "Squamata",
  "Plagiopholis styani", "Colubridae", "Squamata",
  "Psammodynastes pulverulentus", "Colubridae", "Squamata",
  "Pseudoxenodon bambusicola", "Colubridae", "Squamata",
  "Pseudoxenodon karlschmidti", "Colubridae", "Squamata",
  "Pseudoxenodon macrops", "Colubridae", "Squamata",
  "Pseudoxenodon stejnegeri", "Colubridae", "Squamata",
  "Ptyas korros", "Colubridae", "Squamata",
  "Ptyas mucosus", "Colubridae", "Squamata",
  "Rhabdophis leonardi", "Colubridae", "Squamata",
  "Rhabdophis nigrocinctus", "Colubridae", "Squamata",
  "Rhabdophis nuchalis", "Colubridae", "Squamata",
  "Rhabdophis pentasupralabialis", "Colubridae", "Squamata",
  "Rhabdophis subminiatus", "Colubridae", "Squamata",
  "Rhabdophis tigrinus", "Colubridae", "Squamata",
  "Rhynchophis boulengeri", "Colubridae", "Squamata",
  "Sibynophis chinensis", "Colubridae", "Squamata",
  "Sinonatrix aequifasciata", "Colubridae", "Squamata",
  "Sinonatrix annularis", "Colubridae", "Squamata",
  "Sinonatrix percarinata", "Colubridae", "Squamata",
  "Xenochrophis piscator", "Colubridae", "Squamata",
  "Zaocys dhumnades", "Colubridae", "Squamata",
  "Zaocys nigromarginatus", "Colubridae", "Squamata",
  
  # Squamata: Elapidae
  "Bungarus multicinctus", "Elapidae", "Squamata",
  "Calliophis kelloggi", "Elapidae", "Squamata",
  "Calliophis macclellandi", "Elapidae", "Squamata",
  "Naja naja", "Elapidae", "Squamata",
  "Ophiophagus hannah", "Elapidae", "Squamata",
  
  # Squamata: Viperidae
  "Azemiops feae", "Viperidae", "Squamata",
  "Deinagkistrodon acutus", "Viperidae", "Squamata",
  "Ermia mangshanensis", "Viperidae", "Squamata",
  "Gloydius brevicaudus", "Viperidae", "Squamata",
  "Gloydius strauchi", "Viperidae", "Squamata",
  "Ovophis monticola", "Viperidae", "Squamata",
  "Protobothrops jerdonii", "Viperidae", "Squamata",
  "Protobothrops mucrosquamatus", "Viperidae", "Squamata",
  "Trimeresurus albolabris", "Viperidae", "Squamata",
  "Trimeresurus stejnegeri", "Viperidae", "Squamata",
  "Trimeresurus xiangchengensis", "Viperidae", "Squamata",
  "Trimeresurus yunnanensis", "Viperidae", "Squamata"
) %>%
  mutate(
    class = "Reptilia",
    genus = sub(" .*", "", species)
  ) %>%
  select(species, genus, family, order, class)

# Sava River fishes

library(tibble)

sava_fish <- tribble(
  ~species, ~common_name, ~genus, ~family, ~order, ~nativeness, ~occurring_state, ~occurring_river,
  "Eudontomyzon mariae", "Ukrainian lamprey", "Eudontomyzon", "Petromyzontidae", "Petromyzontiformes", "native", "SI; HR; BA; RS", "Sava catchment; Slovenian tributaries (Sora, Ljubljanica, Mirna, Krka, Kolpa, Savinja, Sotla)",
  "Acipenser ruthenus", "sterlet", "Acipenser", "Acipenseridae", "Acipenseriformes", "native", "SI; HR; BA; RS", "Sava R. (middle/lower section; historically Slovenia)",
  "Salmo trutta", "brown trout", "Salmo", "Salmonidae", "Salmoniformes", "native", "SI; HR; BA; RS; ME", "Sava catchment; upper rhithron streams (e.g., Sora, Ljubljanica, Kolpa, Una, Vrbas, Bosna, Drina)",
  "Salmo marmoratus", "marble trout", "Salmo", "Salmonidae", "Salmoniformes", "alien", "Sava catchment (introduced)", "introduced into Sava catchment streams",
  "Oncorhynchus mykiss", "rainbow trout", "Oncorhynchus", "Salmonidae", "Salmoniformes", "alien", "SI; HR; BA; RS", "Sava catchment; stocked in upper rhithron streams",
  "Salvelinus fontinalis", "brook trout", "Salvelinus", "Salmonidae", "Salmoniformes", "alien", "SI; HR; BA; RS", "Sava catchment; upper rhithron streams",
  "Salvelinus alpinus", "Arctic charr", "Salvelinus", "Salmonidae", "Salmoniformes", "alien", "SI", "Sava catchment; mountain streams (introduced)",
  "Hucho hucho", "huchen / Danube salmon", "Hucho", "Salmonidae", "Salmoniformes", "native", "SI; HR; BA; RS; ME", "Sava, Kolpa, Una, Vrbas, Bosna, Drina, Lim",
  "Thymallus thymallus", "European grayling", "Thymallus", "Salmonidae", "Salmoniformes", "native", "SI; HR; BA; RS; ME", "Sava catchment; upper/middle rhithron streams",
  "Esox lucius", "northern pike", "Esox", "Esocidae", "Esociformes", "native", "SI; HR; BA; RS", "Sava R. and lowland tributaries; oxbows and side arms",
  "Abramis brama", "common bream", "Abramis", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R. and lowland tributaries; potamon/lower rhithron",
  "Abramis sapa", "white-eye bream", "Abramis", "Cyprinidae", "Cypriniformes", "native", "HR; BA; RS", "Sava R. and lowland tributaries",
  "Abramis ballerus", "blue bream", "Abramis", "Cyprinidae", "Cypriniformes", "native", "HR; BA; RS", "Sava R. and lowland tributaries",
  "Blicca bjoerkna", "white bream", "Blicca", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R. and lowland tributaries",
  "Vimba vimba", "vimba", "Vimba", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R. and tributaries",
  "Tinca tinca", "tench", "Tinca", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R.; lowland tributaries, oxbows",
  "Cyprinus carpio", "common carp", "Cyprinus", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R.; lowland rivers and reservoirs",
  "Carassius carassius", "crucian carp", "Carassius", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R.; lowland tributaries, oxbows",
  "Carassius gibelio", "Prussian carp / gibel carp", "Carassius", "Cyprinidae", "Cypriniformes", "alien", "SI; HR; BA; RS", "Sava R.; lowland tributaries, oxbows; highly invasive",
  "Carassius auratus", "goldfish", "Carassius", "Cyprinidae", "Cypriniformes", "alien", "HR; BA; RS", "Sava catchment; lowland waters",
  "Ctenopharyngodon idella", "white grass carp", "Ctenopharyngodon", "Cyprinidae", "Cypriniformes", "alien", "SI; HR; BA; RS", "Sava catchment; lowland rivers and reservoirs",
  "Scardinius erythrophthalmus", "rudd", "Scardinius", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R.; lowland tributaries",
  "Aspius aspius", "asp", "Aspius", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R. and larger tributaries",
  "Rutilus pigus", "Danubian roach", "Rutilus", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R. and tributaries",
  "Rutilus rutilus", "roach", "Rutilus", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R. and tributaries",
  "Rhodeus sericeus", "bitterling", "Rhodeus", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R.; lowland tributaries",
  "Alburnus alburnus", "bleak", "Alburnus", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R. and tributaries",
  "Alburnoides bipunctatus", "spirlin", "Alburnoides", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R.; middle rhithron streams",
  "Phoxinus phoxinus", "minnow", "Phoxinus", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava catchment; upper/middle rhithron streams",
  "Leuciscus souffia", "blageon", "Leuciscus", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava catchment; streams and rivers",
  "Leuciscus leuciscus", "dace", "Leuciscus", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava catchment; streams and rivers",
  "Squalius cephalus", "chub", "Squalius", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R.; middle rhithron streams and rivers",
  "Chondrostoma nasus", "nase", "Chondrostoma", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R.; middle/lower rhithron",
  "Idus idus", "orfe / ide", "Idus", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R. and tributaries",
  "Barbus barbus", "common barbel", "Barbus", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava R.; middle/lower rhithron",
  "Barbus balcanicus", "brook barbel", "Barbus", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava catchment; upper/middle rhithron",
  "Gobio gobio", "gudgeon", "Gobio", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava catchment; streams and rivers",
  "Gobio uranoscopus", "Danubian gudgeon", "Gobio", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava catchment; streams and rivers",
  "Gobio albipinnatus", "whitefin gudgeon", "Gobio", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava catchment; streams and rivers",
  "Gobio kessleri", "Kessler's gudgeon", "Gobio", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava catchment; streams and rivers",
  "Pseudorasbora parva", "topmouth gudgeon", "Pseudorasbora", "Cyprinidae", "Cypriniformes", "alien", "SI; HR; BA; RS", "Sava catchment; lowland and slow-flowing waters",
  "Hypophthalmichthys molitrix", "silver carp", "Hypophthalmichthys", "Cyprinidae", "Cypriniformes", "alien", "HR; BA; RS", "Sava R.; lowland rivers and reservoirs",
  "Hypophthalmichthys nobilis", "bighead carp", "Hypophthalmichthys", "Cyprinidae", "Cypriniformes", "alien", "HR; BA; RS", "Sava R.; lowland rivers and reservoirs",
  "Alburnus sarmaticus", "Danube bleak (recently recognised)", "Alburnus", "Cyprinidae", "Cypriniformes", "native (recently recognised)", "SI; HR; BA; RS", "Sava catchment; rivers",
  "Chalcalburnus chalcoides", "Danube bleak", "Chalcalburnus", "Cyprinidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava catchment; rivers",
  "Barbatula barbatula", "stone loach", "Barbatula", "Nemacheilidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava catchment; upper/middle rhithron streams",
  "Misgurnus fossilis", "weather loach", "Misgurnus", "Cobitidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava catchment; lowland waters, oxbows",
  "Cobitis elongata", "Balkan loach", "Cobitis", "Cobitidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava catchment; streams and rivers",
  "Cobitis elongatoides", "riffe loach", "Cobitis", "Cobitidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava catchment; streams and rivers",
  "Cobitis taenia", "spined loach", "Cobitis", "Cobitidae", "Cypriniformes", "native", "HR; BA; RS", "Sava catchment; lowland waters",
  "Sabanejewia aurata", "golden loach", "Sabanejewia", "Cobitidae", "Cypriniformes", "native", "SI; HR; BA; RS", "Sava catchment; streams and rivers",
  "Silurus glanis", "wels catfish", "Silurus", "Siluridae", "Siluriformes", "native", "SI; HR; BA; RS", "Sava R.; lowland rivers and reservoirs",
  "Ameiurus nebulosus", "brown bullhead", "Ameiurus", "Ictaluridae", "Siluriformes", "alien", "SI; HR; BA; RS", "Sava R.; lowland waters; invasive",
  "Lota lota", "burbot", "Lota", "Lotidae", "Gadiformes", "native", "SI; HR; BA; RS", "Sava R.; larger tributaries",
  "Perca fluviatilis", "Eurasian perch", "Perca", "Percidae", "Perciformes", "native", "SI; HR; BA; RS", "Sava R.; lowland rivers and reservoirs",
  "Gymnocephalus cernuus", "common ruffe", "Gymnocephalus", "Percidae", "Perciformes", "native", "SI; HR; BA; RS", "Sava R.; lowland rivers",
  "Gymnocephalus baloni", "Balon's ruffe", "Gymnocephalus", "Percidae", "Perciformes", "native", "HR; BA; RS", "Sava R.; lower rhithron/potamon",
  "Gymnocephalus schraetser", "striped ruffe", "Gymnocephalus", "Percidae", "Perciformes", "native", "HR; BA; RS", "Sava R.; lower rhithron/potamon",
  "Sander lucioperca", "zander / pikeperch", "Sander", "Percidae", "Perciformes", "native", "SI; HR; BA; RS", "Sava R.; lowland rivers and reservoirs",
  "Zingel zingel", "zingel", "Zingel", "Percidae", "Perciformes", "native", "HR; BA; RS", "Sava R.; lower rhithron/potamon",
  "Zingel streber", "streber", "Zingel", "Percidae", "Perciformes", "native", "SI; HR; BA; RS", "Sava R.; lower rhithron",
  "Lepomis gibbosus", "pumpkinseed", "Lepomis", "Centrarchidae", "Perciformes", "alien", "SI; HR; BA; RS", "Sava R.; lowland waters and reservoirs",
  "Micropterus salmoides", "largemouth bass", "Micropterus", "Centrarchidae", "Perciformes", "alien", "SI; HR; BA; RS", "Sava catchment; lowland waters and reservoirs",
  "Neogobius fluviatilis", "monkey goby", "Neogobius", "Gobiidae", "Perciformes", "alien", "SI; HR; BA; RS", "Sava R.; lower rhithron/potamon",
  "Neogobius melanostomus", "round goby", "Neogobius", "Gobiidae", "Perciformes", "alien", "HR; BA; RS", "Sava R.; lower rhithron/potamon",
  "Neogobius gymnotrachelus", "racer goby", "Neogobius", "Gobiidae", "Perciformes", "alien", "HR; BA; RS", "Sava R.; lower rhithron/potamon",
  "Neogobius kessleri", "bighead goby", "Neogobius", "Gobiidae", "Perciformes", "alien", "HR; BA; RS", "Sava R.; lower rhithron/potamon",
  "Proterorhinus marmoratus", "tubenose goby", "Proterorhinus", "Gobiidae", "Perciformes", "alien", "RS", "Sava R.; lower rhithron/potamon",
  "Umbra krameri", "mudminnow", "Umbra", "Umbridae", "Esociformes", "native", "HR; BA; RS", "Lonja, Zasavica, Gromizelj; lowland wetlands",
  "Cottus gobio", "European bullhead", "Cottus", "Cottidae", "Scorpaeniformes", "native", "SI; HR; BA; RS", "Sava catchment; upper/middle rhithron streams",
  "Cottus metae", "Cottus metae (recently recognised)", "Cottus", "Cottidae", "Scorpaeniformes", "native (recently recognised)", "SI; HR; BA; RS", "Sava catchment; upper/middle rhithron streams",
  "Pelecus cultratus", "Sichel", "Pelecus", "Leuciscidae", "Cypriniformes", "alien", "BA; RS", "Sava catchment; Drina"
)

# validate species name using fishbase
## sava fish species
sava_fish_fb <- sava_fish |>
  mutate(valid_name = rfishbase::validate_names(species))
## check for difference with species
length(unique(sava_fish_fb$species))-length(unique(sava_fish_fb$valid_name)) # no different in length
## find differences
setdiff(unique(sava_fish_fb$species), unique(sava_fish_fb$valid_name))
## Abramis sapa -> Ballerus sapa
## Abramis ballerus -> Ballerus ballerus
## Aspius aspius -> Leuciscus aspius
## Leuciscus souffia -> Telestes souffia
## Idus idus -> Leuciscus idus
## Gobio uranoscopus -> Romanogobio uranoscopus
## Gobio albipinnatus -> Romanogobio albipinnatus
## Gobio kessleri -> Romanogobio kessleri
## Chalcalburnus chalcoides -> Alburnus chalcoides
## Gymnocephalus cernuus -> Gymnocephalus cernua
## Neogobius gymnotrachelus -> Babka gymnotrachelus
## Neogobius kessleri -> Ponticola kessleri
unique(sava_fish_fb$valid_name)

## slovenian fish species
si_fish_fb <- slovenian_fish |>
  mutate(valid_name = rfishbase::validate_names(species))
## check number of FB valid names
length(unique(si_fish_fb$species))-length(unique(si_fish_fb$valid_name)) # 2 differences
## find differences
setdiff(unique(si_fish_fb$species), unique(si_fish_fb$valid_name))
## Rutilus aula -> Leucos aula
## Phoxinus lumaireul -> Phoxinus phoxinus (might not be correct)
## Gasterosteus gymnurus -> Gasterosteus aculeatus
## Gymnocephalus cernuus -> Gymnocephalus cernua

si_fish_fb |>
  mutate(species = case_when(species=="Phoxinus lumaireul"~"Phoxinus lumaireul",
                             species!=valid_name~valid_name,
                             TRUE~species)) |>
  dplyr::select(species:class) |>
  write.csv("data-raw/slovenian_fish_species.csv", row.names = FALSE)
