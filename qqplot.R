library(dplyr)
library(ggplot2)
library(qqplotr)

# Read the file into a tibble
df <- mvmapit_normal_null$pvalues
dp <- list(rate=log(10))
di <- "exp"
de <- FALSE # enabling the detrend option

# TODO: change causal snp facet labels
gg <- df %>% ggplot(mapping = aes(
    sample = -log10(p)
)) +
    theme_bw() +
    stat_qq_band(distribution = di,
                 dparams = dp,
                 detrend = de,
                 alpha = 0.5) +
    stat_qq_line(distribution = di, dparams = dp, detrend = de) +
    stat_qq_point(distribution = di, dparams = dp, detrend = de) +
    theme(legend.position = "none") +
    labs(x = bquote("Theoretical Quantiles " -log[10](p)),
         y = bquote("Sample Quantiles " -log[10](p)))

# # TODO: do we need to produce a different file format?
# ggsave(paste0(args$file_path, ".png"), plot = gg, width = 6, height = 4, dpi = 300)

#null just combine all data --> a lot of data points