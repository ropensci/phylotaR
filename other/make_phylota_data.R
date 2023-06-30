#library(phylotaR)
devtools::load_all()

print(getwd())

# # MAMMALS
# mammals <- readTree(file='other/mammals.tre')
# txnyms <- searchTxnyms(mammals, cache=TRUE, infer=TRUE,
#                        clean=TRUE, parent='Mammalia')
# mammals <- setTxnyms(mammals, txnyms)
# summary(mammals)
# save(mammals, file="data/mammals.rda", compress="xz")

# # BIRDS
# birds <- readTree(file='other/birds.tre')
# txnyms <- searchTxnyms(birds, cache=TRUE, infer=TRUE,
#                        clean=TRUE, parent='Aves')
# birds <- setTxnyms(birds, txnyms)
# summary(birds)
# save(birds, file="data/birds.rda", compress="xz")

# # PLANTS
# plants <- ape::read.tree(file='other/plants.tre')
# plants$root.edge <- NULL
# plants$node.label <- NULL  # treeman cannot handle special characters in nodel label
# plants <- as(plants, 'TreeMan')
# txnyms <- searchTxnyms(plants, cache=TRUE, infer=TRUE,
#                        clean=TRUE, parent='Viridiplantae')
# plants <- setTxnyms(plants, txnyms)
# summary(plants)
# save(plants, file="data/plants.rda", compress="xz")

aotus <- read_phylota(file.path('demos', 'aotus'))
save(aotus, file = file.path('data', 'aotus.rda'), compress = 'xz')
rm(aotus)

sturgeons <- read_phylota(file.path('demos', 'acipenseridae'))
save(sturgeons, file = file.path('data', 'sturgeons.rda'), compress = 'xz')
rm(sturgeons)

dragonflies <- read_phylota(file.path('demos', 'anisoptera'))
save(dragonflies, file = file.path('data', 'dragonflies.rda'), compress = 'xz')
rm(dragonflies)

bromeliads <- read_phylota(file.path('demos', 'bromeliaceae'))
save(bromeliads, file = file.path('data', 'bromeliads.rda'), compress = 'xz')
rm(bromeliads)

cycads <- read_phylota(file.path('demos', 'cycadidae'))
save(cycads, file = file.path('data', 'cycads.rda'), compress = 'xz')
rm(cycads)

tardigrades <- read_phylota(file.path('demos', 'eutardigrada'))
save(tardigrades, file = file.path('data', 'tardigrades.rda'), compress = 'xz')
rm(tardigrades)

tinamous <- read_phylota(file.path('demos', 'tinamiformes'))
save(tinamous, file = file.path('data', 'tinamous.rda'), compress = 'xz')

yeasts <- read_phylota(file.path('demos', 'kazachstania'))
save(yeasts, file = file.path('data', 'yeasts.rda'), compress = 'xz')
rm(yeasts)
