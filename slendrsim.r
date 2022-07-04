library(slendr)

# define a world map
map <- world(
  xrange = c(-12, 56),
  yrange = c(-10, 55),
  crs = "EPSG:3035")

# Creating the regions that will define population boundaries
# expanded African population region:
exp_afr <- region(
  "Expanded African Population", map, 
  polygon = list(c(-10, 29), c(-10, 36), c(-1.8, 35), c(0.4, 32),
  c(6, 28.7), c(13.5, 27), c(20, 27.2), c(34.3, 34.4), c(37,37),
  c(46, 34), c(49, 29.7), c(37, 30.5), c(37, 28.5), c(32, 24), c(35.3, 10.2),
  c(38.44, 11.4), c(55, 8.13), c(44, -6), c(38, 2.9), c(30.9, 5.2),
  c(28, 7), c(0.6, 12), c(-5.4, 26.5))
)

# Moroccan population region:
morocco <- region(
  "Morocco", map,
  polygon = list(c(-10, 29), c(-10, 36), c(-10.8, 35), c(-1.8, 33), c(-5.4, 29))
)

# European region:
eur <- region(
  "Europe", map,
  polygon = list(c(4.3, 42.3), c(-9.3, 44.3), c(-10.6, 35.7), c(-1.7, 35.8))
)

# Middle-Eastern region:
mideast = region(
  "Middle East", map,
  polygon = list(c(37,37), c(46, 34), c(49, 29.7), c(37, 30.5), c(37, 28.5), c(32, 24), c(30, 24), c(28.5, 31.4))
)

# Setting up the population boundaries and dynamics
# ancestral African population (120 kya)
# note the removal 89998 ya - this is meant to simulate the extinction of the central African populations
# due to the shift to arid conditions
afr <- population(
  "AFR", parent = "ancestor", time = 120000, N = 4000, remove = 
  map = map, center = c(15, 18), radius = 500e3
)

# African population range expansion (120-90 kya)
afr <- afr %>% expand_range(by = 3500e3, start = 120000, end = 90000,
polygon = exp_afr, snapshots = 30)

# Moroccan population split (90  kya)
mor <- population(
  "MOR", parent = afr, time = 90000, N = 3000,
  map, polygon = morocco
)

# Middle-eastern population split (90 kya)
mer = population(
  "MER", parent = afr, time = 90000, N = 4000,
  map, polygon = mideast
)

# Expansion of Moroccan population into European region (90-80 kya)
mor <- mor %>% expand_range(
  by = 2500e3, start = 90000, end = 80000,
  snapshots = 30, polygon = join(mor, eur)
)

# Migration from the Middle East to Europe (~90-30 kya)
mig <- population(
  "MIG", parent = mer, time = 89999, N = 3000, center = c(41, 31.7), radius = 500e3
) %>% move(
  trajectory = list(c(40, 35), c(15.2, 48.2), c(3.8, 46.8), c(-4.3, 40.1)),
  start = 89999, end = 30000, snapshots = 30
)

# Setting up the gene flow rates
# The only gene-flow we're setting up for now is from the Iberian population to the
# migrated European population (30-20 kya) -- this is subject to change
gf <- gene_flow(from = mor, to = mig, rate = 0.5, start = 30000, end = 20000)

# Compiling the model
model_dir <- paste0(tempfile(), "_arabidopsis-test")
model <- compile_model(populations = list(afr, mor, mer, mig),gene_flow = gf,generation_time = 1,resolution = 10e3,
  competition = 130e3, mating = 100e3,dispersal = 70e3,path = model_dir)
  
# Running the simulation on slim
# Change method parameter to "batch" to get a tsv.gz file to produce animated GIF from
# Otherwise, if this isn't required, keep the following unchanged to open the model on SLiMgui:
slim(model, sequence_length = 1, recombination_rate = 0, save_locations = TRUE, method = "gui", random_seed = 314159)


