library(slendr)

map <- world(xrange = c(-12, 56), yrange = c(-10, 55), crs = "EPSG:3035")

# Creating the regions that will define population boundaries
# expanded African population region:
exp_afr <- region("Expanded African Population", map, 
  polygon = list(c(-10, 29), c(-10, 36), c(-1.8, 35), c(0.4, 32),
  c(6, 28.7), c(13.5, 27), c(20, 27.2), c(34.3, 34.4), c(37,37),
  c(46, 34), c(49, 29.7), c(37, 30.5), c(37, 28.5), c(32, 24), c(35.3, 10.2),
  c(38.44, 11.4), c(55, 8.13), c(44, -6), c(38, 2.9), c(30.9, 5.2),
  c(28, 7), c(0.6, 12), c(-5.4, 26.5))
)

# Moroccan population region:
morocco <- region("Morocco", map,
  polygon = list(c(-10, 29), c(-10, 36), c(-10.8, 35), c(-1.8, 33), c(-5.4, 29))
)

# European region:
europe <- region("Europe", map,
  polygon = list(c(4.3, 42.3), c(-9.3, 44.3), c(-10.6, 35.7), c(-1.7, 35.8))
)

# Middle-Eastern region:
mideast <- region("Middle East", map,
  polygon = list(c(37,37), c(46, 34), c(49, 29.7), c(37, 30.5), c(37, 28.5),
  c(32, 24), c(30, 24), c(28.5, 31.4))
)

# Expanded Middle-Eastern region:
migration <- region("Migratory Region", map,
 polygon = list(c(43.3, 33.7), c(40, 40.6), c(30.5, 48), c(18.9, 50.5), c(-3.2, 49.5),
 c(-10.5, 45.2), c(1.95, 40.1), c(13.7, 45.3), c(35.2, 35.6), c(33.6, 30.4))
)

# Setting up the population boundaries and dynamics
# ancestral African population (120 kya)
# note the removal 89998 ya - this is meant to simulate the extinction of the central 
# African populations due to the shift to arid conditions
afr <- population("AFR", parent = "ancestor", time = 1, N = 12000, remove = 30002, map = map, 
center = c(15, 18), radius = 600e3, dispersal = 100e3)

# African population range expansion (120-90 kya)
afr <- afr %>% expand_range(by = 3500e3, start = 1, end = 30001, polygon = exp_afr, 
snapshots = 30, lock = TRUE)

# Moroccan population split (90  kya)
mor <- population("MOR", parent = afr, time = 30001, N = 8000, dispersal = 100e3, map, 
polygon = morocco)

# Middle-eastern population split (90 kya)
mer = population("MER", parent = afr, time = 90000, N = 5000, dispersal = 1000e3, map, 
polygon = mideast)

# Expansion of Moroccan population into European region (90-85 kya)
mor <- mor %>% expand_range(by = 2500e3, start = 30001, end = 35001, snapshots = 30, 
lock = TRUE, polygon = join(morocco, europe))

# Expansion of Middle-Eastern population into European region (45-30 kya)
mer <- mer %>% expand_range(by = 5500e3, start = 75001, end = 90001, snapshots = 30, 
lock = TRUE, polygon = join(europe, join(migration, mideast)))

# Setting up the gene flow rates
# The only gene-flow we're setting up for now is from the Iberian population to the
# migrated European population (30-20 kya) -- this is subject to change
gf <- list(gene_flow(from = mor, to = mer, rate = 0.5, start = 90001, end = 120000),
list(gene_flow(from = mer, to = mor, rate = 0.5, start = 90001, end = 120000))

# Compiling the model
model_dir <- paste0(tempfile(), "_arabidopsis-test")
model <- compile_model(populations = list(afr, mor, mer),gene_flow = gf,generation_time = 1,
resolution = 10e3, competition = 1, mating = 1,dispersal = 1,path = model_dir)
  
# Running the simulation on slim
# Change method parameter to "batch" to get a tsv.gz file to produce animated GIF from
# Otherwise, if this isn't required, keep the following unchanged to open the model on SLiMgui:
slim(model, sequence_length = 1e6, recombination_rate = 8e-10, save_locations = TRUE, 
method = "gui", random_seed = 314159)

present_samples <- schedule_sampling(model, times = 120001, list(mor, 50),
locations = list(c(-7, 42.3), c(-6.7, 38.3), c(-3.7, 41), c(-0.11, 41.5), c(-1.1, 38.9), 
c(-5.5, 37)))