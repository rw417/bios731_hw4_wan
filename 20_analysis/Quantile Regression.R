source(paste(here::here(), "/10_source/my_rq.R", sep=""))

# Load the dataset
cannabis <- readRDS(paste(here::here(), "02_cleaned_data/cannabis-2.rds", sep = "/"))

# Run quantile regression
beta_t0.25 <- my_rq(
  y=cannabis$t_mmr1,
  X=cannabis[,c("p_change", "h_hr", "i_composite_score")],
  tau=0.25
)

beta_t0.5 <- my_rq(
  y=cannabis$t_mmr1,
  X=cannabis[,c("p_change", "h_hr", "i_composite_score")],
  tau=0.5
)

beta_t0.75 <- my_rq(
  y=cannabis$t_mmr1,
  X=cannabis[,c("p_change", "h_hr", "i_composite_score")],
  tau=0.75
)

# Run linear regression
lm_model <- lm(t_mmr1 ~ p_change + h_hr + i_composite_score, data=cannabis)
