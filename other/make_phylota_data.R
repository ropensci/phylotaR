library(phylotaR)

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
