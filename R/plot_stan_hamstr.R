# Methods -------

#' Title
#'
#' @param object a hamstr_fit object
#' @param type One of "default", "age_models", "hier_acc_rates",
#'                               "acc_mean_prior_post", "mem_prior_post"
#' @inheritParams plot_hamstr  
#' @return A ggplot object
#'
#' @examples
#' @export
#' @method plot hamstr_fit
plot.hamstr_fit <- function(object,
                            type = c("default",
                              "age_models", "hier_acc_rates",
                              "acc_mean_prior_post", "mem_prior_post"
                              ),
                            summarise = TRUE,
                            ...){
  
  type <- match.arg(type)
  
  switch(type,
         default = plot_hamstr(object, summarise = summarise, ...),
         age_models = plot_hamstr(object, summarise = summarise,
                                  plot_diagnostics  = FALSE, ...),
         hier_acc_rates = plot_hierarchical_acc_rate(object),
         acc_mean_prior_post = plot_acc_mean_prior_posterior(object),
         mem_prior_post = plot_memory_prior_posterior(object))
  
}


# Functions ------

#' Plot an hamstr_fit object
#'
#' @param hamstr_fit The object returned from \code{stan_hamstr}.
#'
#' @param n.iter The number of iterations of the model to plot, defaults to
#'   1000.
#' @param summarise logical TRUE or FALSE. Plot the realisations as a summarised
#'  "ribbon" showing 50% and 95% intervals (faster), or as a spaghetti plot 
#'  showing individual realisations. Defaults to TRUE (ribbon).
#' @param plot_diagnostics logical, include diagnostic plots: traceplot of
#'   log-posterior, hierarchical accumulations rates, memory parameter. Defaults
#'   to TRUE.
#' @description Plots the HAMStR modelled age ~ depth relationship together with
#'   the depths, ages, and age uncertainties in the observed data. A random
#'   sample of size \code{n.iter} of the iterations of the posterior
#'   distribution are plotted as grey lines. The observed data are plotted as
#'   points with +- 2*se error bars.
#'   
#' @return A ggplot2 object
#' @export
#' @importFrom ggpubr ggarrange
#' @importFrom rstan extract
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error,
#'   K = c(10, 10), nu = 6,
#'   acc_mean_prior = 20,
#'   mem_mean = 0.7, mem_strength = 4,
#'   inflate_errors = 0,
#'   iter = 2000, chains = 3)
#'   
#' # With age models summarised as a ribbon. Faster than spaghetti plots.
#' plot_hamstr(fit)
#' 
#' # With age models as spaghetti plots. Can see individual realisations, but slower to plot.
#' plot_hamstr(fit, summarise = FALSE)
#' }
plot_hamstr <- function(hamstr_fit, summarise = TRUE, n.iter = 1000, plot_diagnostics = TRUE) {

  #summarise <- match.arg(summarise)

  if (summarise == TRUE){
    p.fit <- plot_summary_age_models(hamstr_fit)
  } else if (summarise == FALSE){
    p.fit <- plot_age_models(hamstr_fit, n.iter = n.iter)
  }

  if (plot_diagnostics == FALSE) return(p.fit)

  if (plot_diagnostics){
    p.mem <- plot_memory_prior_posterior(hamstr_fit)
    p.acc <- plot_hierarchical_acc_rate(hamstr_fit)
    }

  t.lp <- rstan::traceplot(hamstr_fit$fit, pars = c("lp__"), include = TRUE) +
    ggplot2::theme(legend.position = "top") +
    ggplot2::labs(x = "Iteration")

  ggpubr::ggarrange(
    ggpubr::ggarrange(t.lp, p.acc, p.mem, ncol = 3, widths = c(3,3,2)),
    p.fit,
    nrow = 2, heights = c(1, 2))
}

#' Plot Summary of Posterior Age Models
#'
#' @inheritParams plot_hamstr
#'
#' @return A ggplot2 object
#' @keywords internal
#' @importFrom readr parse_number
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error,
#'   K = c(10, 10), nu = 6,
#'   acc_mean_prior = 20,
#'   mem_mean = 0.7, mem_strength = 4,
#'   inflate_errors = 0,
#'   iter = 2000, chains = 3)
#'   
#' plot_summary_age_models(fit)
#' }
plot_summary_age_models <- function(hamstr_fit){
  
  age_summary <- summarise_age_models(hamstr_fit)
  
  obs_ages <- data.frame(
    depth = hamstr_fit$data$depth,
    age = hamstr_fit$data$obs_age,
    err = hamstr_fit$data$obs_err)
  
  obs_ages <- dplyr::mutate(obs_ages,
                            age_upr = age + 2*err,
                            age_lwr = age - 2*err)
  
  
  infl_errs <- rstan::summary(hamstr_fit$fit, par = "obs_err_infl")$summary %>% 
    tibble::as_tibble(., rownames = "par") %>% 
    dplyr::mutate(dat_idx = readr::parse_number(par))
  
  p.age.sum <- age_summary %>% 
    ggplot2::ggplot(ggplot2::aes(x = depth, y = mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymax = `2.5%`, ymin = `97.5%`), fill = "Lightgrey") +
    ggplot2::geom_ribbon(ggplot2::aes(ymax = `75%`, ymin = `25%`), fill = "Darkgrey") +
    ggplot2::geom_line() +
    ggplot2::geom_line(ggplot2::aes(y = `50%`), colour = "Green") +
    ggplot2::labs(x = "Depth", y = "Age") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  
  
  if (hamstr_fit$data$inflate_errors == 1){
    obs_ages <- obs_ages %>% 
      dplyr::mutate(infl_err = infl_errs$mean,
                    age_lwr_infl = age + 2*infl_err,
                    age_upr_infl = age - 2*infl_err)
    
    p.age.sum <- p.age.sum +
      ggplot2::geom_linerange(
        data = obs_ages,
        ggplot2::aes(x = depth, ymax = age_upr_infl, ymin = age_lwr_infl),
        group = NA,
        colour = "Red",
        alpha = 0.5, inherit.aes = F)
  }
  
  p.age.sum <- p.age.sum +
    ggplot2::geom_linerange(data = obs_ages,
                            ggplot2::aes(x = depth, 
                                         ymax = age_upr, ymin = age_lwr), inherit.aes = FALSE,
                            colour = "Blue", size = 1.25) +
    ggplot2::geom_point(data = obs_ages, ggplot2::aes(y = age),
                        colour = "Blue")
  
  
  p.age.sum <- add_subdivisions(p.age.sum, hamstr_fit)
  
  p.age.sum
}

#' Plot Age Models as Spaghetti Plot
#'
#' @inheritParams plot_hamstr 
#' 
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom readr parse_number
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error,
#'   K = c(10, 10), nu = 6,
#'   acc_mean_prior = 20,
#'   mem_mean = 0.7, mem_strength = 4,
#'   inflate_errors = 0,
#'   iter = 2000, chains = 3)
#'   
#' plot_age_models(fit)
#' }
plot_age_models <- function(hamstr_fit, n.iter = 1000){
  
  
  posterior_ages <- get_posterior_ages(hamstr_fit)
  
  obs_ages <- dplyr::tibble(
    depth = hamstr_fit$data$depth,
    age = hamstr_fit$data$obs_age,
    err = hamstr_fit$data$obs_err)
  
  obs_ages <- dplyr::mutate(obs_ages,
                            age_upr = .data$age + 2*.data$err,
                            age_lwr = .data$age - 2*.data$err)
  
  infl_errs <- rstan::summary(hamstr_fit$fit, par = "obs_err_infl")$summary %>% 
    tibble::as_tibble(., rownames = "par") %>% 
    dplyr::mutate(dat_idx = readr::parse_number(.data$par))
  
  p.fit <- posterior_ages %>%
    dplyr::filter(.data$iter %in% sample(unique(.data$iter), n.iter, replace = FALSE)) %>%
    ggplot2::ggplot(ggplot2::aes(x = depth, y = age, group = iter))
  
  
  p.fit <- p.fit +
    ggplot2::geom_line(alpha = 0.5 / sqrt(n.iter))
  
  if (hamstr_fit$data$inflate_errors == 1){
    obs_ages <- obs_ages %>% 
      dplyr::mutate(infl_err = infl_errs$mean,
                    age_lwr_infl = .data$age + 2*.data$infl_err,
                    age_upr_infl = .data$age - 2*.data$infl_err)
    
    p.fit <- p.fit +
      ggplot2::geom_linerange(
        data = obs_ages,
        ggplot2::aes(x = depth, ymax = age_upr_infl, ymin = age_lwr_infl),
        group = NA,
        colour = "Red",
        alpha = 0.5, inherit.aes = F)
  }
  
  p.fit <- p.fit +
    ggplot2::geom_linerange(
      data = obs_ages,
      ggplot2::aes(ymax = age_upr, ymin = age_lwr),
      group = NA,
      colour = "Blue",
      size = 1.2,
      alpha = 1) +
    ggplot2::geom_point(
      data = obs_ages,
      ggplot2::aes(y = age),
      group = NA,
      colour = "Blue",
      #size = 1.01,
      alpha = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(x = "Depth", y = "Age")
  
  
  # add subdivisions
  p.fit <- add_subdivisions(p.fit, hamstr_fit)
  
  return(p.fit)
  
}


#' Plot a Prior and Posterior
#'
#' @param prior 
#' @param posterior 
#'
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
plot_prior_posterior_hist <- function(prior, posterior){
  clrs <- c("Posterior" = "Blue", "Prior" = "Red")
  ggplot2::ggplot() +
    ggplot2::geom_histogram(data = posterior,
                   ggplot2::aes(x = x, ggplot2::after_stat(density),
                       fill = "Posterior"),
                   alpha = 0.5, bins = 100) +
    ggplot2::geom_line(data = prior, ggplot2::aes(x = x, y = d, colour = "Prior")) +
    ggplot2::facet_wrap(~par, scales = "free") +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::scale_colour_manual(values = clrs) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = c(0.5))))+
    ggplot2::labs(
      x = "Value",
      y = "Density",
      colour = "",
      fill = ""
    ) +
    ggplot2::theme_bw() 
}


#' Plot the Prior and Posterior Distributions of the Inflation Factor Parameters
#'
#' @return A ggplot2 object
#' @import rstan 
#' @import ggplot2
#' @importFrom readr parse_number
#' @importFrom rlang .data
#' @inheritParams plot_hamstr
#' @keywords internal
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error,
#'   K = c(10, 10), nu = 6,
#'   acc_mean_prior = 20,
#'   mem_mean = 0.7, mem_strength = 4,
#'   inflate_errors = 0,
#'   iter = 2000, chains = 3)
#'   
#' plot_infl_prior_posterior(fit)
#' }
plot_infl_prior_posterior <- function(hamstr_fit){
  
  clrs <- c("Posterior" = "Blue", "Prior" = "Red")
  
  hamstr_dat <- hamstr_fit$data
  
  infl_mean_shape_post <-
    tibble::tibble(infl_mean = as.vector(rstan::extract(hamstr_fit$fit, "infl_mean")[[1]]),
           infl_shape = as.vector(rstan::extract(hamstr_fit$fit, "infl_shape")[[1]])) %>% 
    dplyr::mutate(iter = 1:dplyr::n())
  
  
  max_x_shape <- with(hamstr_dat, {
    infl_shape_prior_upr <- stats::qgamma(c(0.99), shape = infl_shape_shape, rate =  infl_shape_shape / infl_shape_mean)
  
     max(c(infl_shape_prior_upr, infl_mean_shape_post$infl_shape))
  }) 
  
  max_x_mean <- with(hamstr_dat, {
    infl_mean_prior_upr <- stats::qnorm(c(0.99), 0, infl_sigma_sd)
    max(c(infl_mean_prior_upr, infl_mean_shape_post$infl_mean))
  })
  
  
  infl_fac <- rstan::extract(hamstr_fit$fit, "infl")[[1]] %>% 
    tibble::as_tibble() %>% 
    tidyr::gather() %>% 
    dplyr::mutate(key = readr::parse_number(.data$key))
  
 
  p.infl.fac <- rstan::stan_plot(hamstr_fit$fit, pars = "infl")
  
 
  
  infl_prior_shape <-
    tibble::tibble(x = seq(0, max_x_shape, length.out = 1000)) %>%
    dplyr::mutate(
      #infl_mean = 2*dnorm(x, 0,  sd = hamstr_dat$infl_sigma_sd),
      d = stats::dgamma(x-1, hamstr_dat$infl_shape_shape,  rate = hamstr_dat$infl_shape_shape / hamstr_dat$infl_shape_mean),
      par = "infl_shape")
  
  infl_prior_mean <-
    tibble::tibble(x = seq(0, max_x_mean, length.out = 1000)) %>%
    dplyr::mutate(
      d = 2*stats::dnorm(.data$x, 0,  sd = hamstr_dat$infl_sigma_sd),
      par = "infl_mean")
  
  infl_priors <- dplyr::bind_rows(infl_prior_mean, infl_prior_shape)
  
  infl_mean_shape_post_long <- infl_mean_shape_post %>% 
    tidyr::gather(.data$par, .data$x, -.data$iter)
 
  p.pars <- plot_prior_posterior_hist(infl_priors, infl_mean_shape_post_long)
  
  
  infl_mean_shape_post <- infl_mean_shape_post %>% 
    dplyr::mutate(q99 = stats::qgamma(0.75, .data$infl_shape,
                                      rate = .data$infl_shape/.data$infl_mean))
  
  
  infl_pars_prior_dist <- infl_mean_shape_post %>% 
    stats::filter(iter %in% sample.int(dplyr::n(), 10)) %>% 
    tidyr::crossing(., tibble::tibble(x = exp(seq(log(0.01), log(stats::quantile(infl_mean_shape_post$q99, prob = 0.95)), length.out = 100)))) %>% 
    dplyr::mutate(d = stats::dgamma(x, shape = .data$infl_shape, rate = .data$infl_shape / .data$infl_mean),
           #d = dgamma(x, shape = infl_shape, rate = infl_shape / 1),
           par = "Modelled prior for infl_fac")
  
  
  p.priors <- infl_pars_prior_dist  %>% 
    ggplot2::ggplot(ggplot2::aes(x = x, y = d, group = iter)) +
    ggplot2::geom_line(alpha = 1,#/sqrt(100),
              colour = "Red") +
    ggplot2::theme_bw() +
    #facet_wrap(~par+iter, scales = "free") +
    ggplot2::labs(y = "Density", x = "Value") 
  
  p <- ggpubr::ggarrange(plotlist = list(p.pars, p.priors,
                                   p.infl.fac), ncol = 2)
  
  return(p)
  
}


#' Plot Mean Accumulation Rate Prior and Posterior Distributions
#' @inheritParams plot_hamstr
#' 
#' @import ggplot2
#' @importFrom rlang .data
#' @return A ggplot2 object
#' @keywords internal
#' @examples 
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error,
#'   K = c(10, 10), nu = 6,
#'   acc_mean_prior = 20,
#'   mem_mean = 0.7, mem_strength = 4,
#'   inflate_errors = 0,
#'   iter = 2000, chains = 3)
#'   
#' plot_acc_mean_prior_posterior(fit)
#' }
plot_acc_mean_prior_posterior <- function(hamstr_fit) {
  clrs <- c("Posterior" = "Blue", "Prior" = "Red")
  
  hamstr_dat <- hamstr_fit$data
  
  prior_mean <- hamstr_dat$acc_mean_prior
  
  acc_prior_rng <- stats::qnorm(c(0.99), mean = 0, sd = 10 * prior_mean)
  
  acc_prior <-
    tibble::tibble(acc_rate = seq(0, acc_prior_rng[1], length.out = 1000)) %>%
    dplyr::mutate(
      density = 2 * stats::dnorm(.data$acc_rate, 0, 10 * prior_mean),
      density = ifelse(.data$acc_rate <= 0, 0, .data$density)
    )
  
  acc_post <-
    tibble::tibble(alpha = as.vector(rstan::extract(hamstr_fit$fit, "alpha[1]")[[1]]))
  
  p <- acc_prior %>%
    ggplot2::ggplot(ggplot2::aes(x = acc_rate, y = density)) +
    # plot the posterior first
    ggplot2::geom_histogram(
      data = acc_post,
      ggplot2::aes(x = alpha, ggplot2::after_stat(density), fill = "Posterior"),
      inherit.aes = FALSE,
      alpha = 0.5,
      # set the colour for the outline of the bins but don't include in colour 
      # legend
      colour = clrs["Posterior"],
      bins = 100
    ) +
    ggplot2::geom_line(ggplot2::aes(colour = "Prior")) +
    ggplot2::labs(
      x = "Mean accumulation rate",
      y = "Density",
      colour = "",
      fill = ""
    ) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::scale_colour_manual(values = clrs) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = c(0.5)))) +
    
    ggplot2::theme_bw()
  
  return(p)
  
}


#' Plot Memory Prior and Posterior
#'
#' @inheritParams plot_hamstr
#'
#' @return A ggplot2 object
#' @import ggplot2
#' @importFrom rlang .data
#' 
#' @examples
#' @keywords internal
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error,
#'   K = c(10, 10), nu = 6,
#'   acc_mean_prior = 20,
#'   mem_mean = 0.7, mem_strength = 4,
#'   inflate_errors = 0,
#'   iter = 2000, chains = 3)
#'   
#' plot_memory_prior_posterior(fit)
#' }
plot_memory_prior_posterior <- function(hamstr_fit){
  # memory prior
  mem.prior <- tibble::tibble(mem = seq(0, 1, length.out = 1000)) %>%
    dplyr::mutate(mem.dens = stats::dbeta(.data$mem, shape1 = hamstr_fit$data$mem_alpha,
                            shape2 = hamstr_fit$data$mem_beta))

  w <- rstan::extract(hamstr_fit$fit, "w")$w
  ifelse(is.matrix(w), w <- w[,ncol(w)],  w <- as.vector(w))

  mem.post <- tibble::tibble(w = w,
                     R = as.vector(rstan::extract(hamstr_fit$fit, "R")$R))


  p.mem <- mem.prior %>%
    ggplot2::ggplot(ggplot2::aes(x = mem, y = mem.dens)) +
    ggplot2::geom_density(data = mem.post, ggplot2::aes(x = R, fill = "at 1 cm 'R'"),
                 inherit.aes = FALSE, show.legend = TRUE) +
    ggplot2::geom_density(data = mem.post, ggplot2::aes(x = w, fill = "between\nsections 'w'"),
                 inherit.aes = FALSE, show.legend = TRUE) +
    ggplot2::geom_line(colour = "Red") +
    ggplot2::scale_x_continuous("Memory [correlation]", limits = c(0, 1)) +
    ggplot2::scale_y_continuous("") +
    ggplot2::scale_fill_discrete("") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")

  return(p.mem)
}

#' Add subdivision tickmarks 
#'
#' @param gg 
#' @inheritParams plot_hamstr
#'
#' @return A ggplot2 object
#' 
#' @import ggplot2
#'
#' @keywords internal
add_subdivisions <- function(gg, hamstr_fit){

  tick_dat <- hierarchical_depths(hamstr_fit$data)

  for (x in seq_along(tick_dat)){

    df <- data.frame(x = tick_dat[[x]])

    lnth <- length(tick_dat) - (x-1)

    gg <- gg + ggplot2::geom_rug(data = df, ggplot2::aes(x = x),
                        inherit.aes = F, sides = "top",
                        length = ggplot2::unit(0.01*lnth, "npc"))

  }

  return(gg)
}




#' Plot the hierarchical accumulation rate parameters
#'
#' @inheritParams plot_hamstr
#'
#' @return ggplot2 object
#' @keywords internal
#' @import ggplot2
#' @importFrom readr parse_number
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error,
#'   K = c(10, 10), nu = 6,
#'   acc_mean_prior = 20,
#'   mem_mean = 0.7, mem_strength = 4,
#'   inflate_errors = 0,
#'   iter = 2000, chains = 3)
#'   
#' plot_hierarchical_acc_rate(fit)
#' }
plot_hierarchical_acc_rate <- function(hamstr_fit){

  idx <- tibble::as_tibble(alpha_indices(hamstr_fit$data$K)[1:3]) %>%
    dplyr::mutate(alpha_idx = (alpha_idx))

  a3 <- rstan::summary(hamstr_fit$fit, pars = "alpha")$summary
  
  alph <- tibble::as_tibble(a3, rownames = "par") %>%
    dplyr::mutate(alpha_idx = readr::parse_number(par)) %>%
    dplyr::left_join(idx, .) %>%
    dplyr::mutate(lvl = factor(lvl))

  # for each unit at each level in hierarchy get max and min depth 
  alph$depth1 <- c(min(hamstr_fit$data$modelled_depths),
                   unlist(sapply((
                     hierarchical_depths(hamstr_fit$data)
                   ),
                   function(x) {
                    utils::head(x, -1)
                   })))
  
  alph$depth2 <- c(max(hamstr_fit$data$modelled_depths),
                   unlist(sapply((
                     hierarchical_depths(hamstr_fit$data)
                   ),
                   function(x) {
                     utils::tail(x, -1)
                   })))
  
  alph2 <- alph %>% 
    dplyr::select(lvl, alpha_idx, depth1, depth2, mean) %>% 
    dplyr::group_by(lvl) %>% 
    tidyr::gather(type, depth, -mean, -lvl, -alpha_idx) %>% 
    dplyr::select(lvl, alpha_idx, depth, mean) %>% 
    dplyr::arrange(lvl, alpha_idx, depth, mean)
  

  gg <- alph2 %>%
    ggplot2::ggplot(ggplot2::aes(x = depth, y = mean, colour = lvl)) +
    ggplot2::geom_path() +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(y = "Accummulation rate [age/depth]", x = "Depth",
         colour = "Hierarchical\nlevel") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "top")

  return(gg)
}

