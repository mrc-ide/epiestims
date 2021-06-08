vary_si <- readRDS("results/vary_si_raw.rds")

summary_epsilon <- map_depth(
  vary_si, 3, function(x){
    if (inherits(x, "error")) NULL
    else summarise_epsilon(x)
  }
)

summary_epsilon <- map(summary_epsilon, function(incid_eps) {
  map_dfr(
    incid_eps, function(tmax_eps) bind_rows(tmax_eps, .id = "sim"),
    .id = "tmax"
  )
 }
)
sim_params <- readRDS("results/vary_si_sim_params.rds")
params <- as.list(sim_params)
params <- append(
  params, list(est_epsilon = summary_epsilon)
)

result <- pmap_dfr(
  params, function(rt_ref, epsilon, si_mu_variant,
                   si_std_variant, est_epsilon) {
    x <- data.frame(
      rt_ref = rt_ref, true_epsilon = epsilon,
      si_mu_variant = si_mu_variant,
      si_std_variant = si_std_variant
    )
    cbind(x, est_epsilon)
  }
)

##result <- readRDS("results/vary_si_processed.rds")
params <- as.list(sim_params)
params <- append(params, list(est_epsilon = vary_si))

eps_error <- pmap_dfr(
 params, function(rt_ref, epsilon, si_mu_variant, si_std_variant,
                  est_epsilon) {

    ## est_epsilon is list of length tmax_all
    map_dfr(est_epsilon, function(est_sim) {
      ## est_sim is list of length 100 i.e. the number of simulations
      out <- map_dfr(est_sim, function(x) {
        ## x is the actual result we are after
        summarise_epsilon_error(x, epsilon)
      }, .id = "sim"
      )
      pars <- data.frame(
        rt_ref = rt_ref, true_epsilon = epsilon,
        si_mu_variant = si_mu_variant,
      si_std_variant = si_std_variant
      )
      cbind(pars, out)
    }, .id = "tmax"
    )
  }
)


ggplot(eps_error, aes(tmax, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_linerange(alpha = 0.3, position = position_dodge(width = 0.1)) +
  facet_grid(si_mu_variant~rt_ref) +
  ylab("Estimated - True value") +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())


qntls <- tidyr::gather(eps_error, key = qntl, value = val, `2.5%`:`97.5%`)
cri_95 <- qntls[qntls$qntl %in% c("2.5%", "50%", "97.5%"), ]

p <- ggplot(cri_95) +
  geom_boxplot(aes(tmax, val, fill = qntl)) +
  facet_grid(si_mu_variant~rt_ref, labeller = label_context) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
  ylab("Estimated - True epsilon") +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())




ggplot(result, aes(tmax, `50%`)) + geom_point() +
  facet_wrap(~true_epsilon, scales = "free_y")


## Summary
by_tmax <- group_by(result, tmax, true_epsilon, rt_ref) %>%
  summarise(
    n = sum(true_epsilon > `2.5%` & true_epsilon < `97.5%`),
    total = n()
  ) %>% ungroup()

x <- filter(result, tmax == 10, si_mu_variant == result$si_mu_variant[1])

p <- ggplot(by_tmax, aes(tmax, n / total, col = factor(rt_ref))) +
  geom_point() +
  facet_wrap(~true_epsilon) +
  ylim(0, 1) +
  theme_minimal() +
  theme(legend.position = "top")


ggplot(result) +
  geom_linerange(aes(x = sim, ymin = `2.5%`, ymax = `97.5%`, col = tmax)) +
  facet_wrap(~true_epsilon, scales = "free_y")


ggplot(result) +
  geom_boxplot(aes(tmax, `50%`)) +
  facet_wrap(~true_epsilon, scales = "free_y")
