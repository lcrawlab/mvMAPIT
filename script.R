library(mvMAPIT)
mvmapit(t(simulated_data$genotype[1:100,1:10]),
        t(simulated_data$trait[1:100,]),
        cores = 4,
        logLevel = "DEBUG"
)
