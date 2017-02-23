library(qtl)
library(mpMap)

set.seed(1234)
bp_map <- sim.map(len = sample(50:150, 6), n.mar = sample(30:70, 6), 
                  include.x = FALSE)

bp_cross <- sim.cross(bp_map, type = "f2")

m4_ped <- sim.mpped(4, 3, 1, 6, 6)

m4_cross <- sim.mpcross(bp_map, m4_ped)
