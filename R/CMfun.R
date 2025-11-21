
#' @name CMfigures
#' @title Figures for conditioning model results
#'
#' @description Functions used by the markdown report to generate summary figures from the age-structured conditoning model
#' @param stanfit Output from [sample_CM()]
#' @param vars Character vector for variable names (see `names(stanfit@sim$samples[[1]])`). Regex and partial matching supported because it is passed to the `pattern` argument of [`grepl()`]
#' @param inc_warmup Logical, whether to include warmup MCMC samples
#' @returns
#' - `CM_trace()` returns a ggplot showing the MCMC trace plot (aka wormplot)
#' @importFrom dplyr bind_rows bind_cols mutate
#' @export
CM_trace <- function(stanfit, vars, inc_warmup = FALSE) {
  if (missing(vars)) vars <- names(stanfit@sim$samples[[1]])

  val <- lapply(1:stanfit@sim$chains, function(i) {
    x <- stanfit@sim$samples[[i]]
    vars_regex <- paste0(vars, collapse = "|")

    x_pars <- x[grepl(vars_regex, names(x))]

    if (length(x_pars)) {

      if (inc_warmup) {
        warmup <- 1
      } else {
        warmup <- stanfit@stan_args[[i]]$warmup
      }
      thin <- stanfit@stan_args[[i]]$thin
      iter <- stanfit@stan_args[[i]]$iter

      it <- seq(1, iter, by = thin)

      x[grepl(vars_regex, names(x))] %>%
        bind_cols() %>%
        as.data.frame() %>%
        mutate(Chain = i, Iteration = it) %>%
        filter(.data$Iteration > warmup) %>%
        reshape2::melt(id.vars = c("Chain", "Iteration"))
    } else {
      data.frame()
    }

  }) %>%
    bind_rows()

  if (nrow(val)) {
    g <- ggplot(val, aes(.data$Iteration, .data$value, colour = factor(.data$Chain))) +
      geom_line() +
      facet_wrap(vars(variable), scales = "free_y") +
      labs(colour = "Chain")
    g
  }
}


#' @rdname CMfigures
#' @returns
#' - `CM_pairs()` returns output from [graphics::pairs()], a matrix of scatterplots of MCMC posterior samples
#' @importFrom stats cor
#' @importFrom graphics pairs rect strwidth text
#' @export
CM_pairs <- function(stanfit, vars, inc_warmup = FALSE) {
  panel.hist <- function(x, ...) {
    usr <- par("usr"); on.exit(par(usr = usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
  }

  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr"); on.exit(par(usr = usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = 1.3)
  }

  val <- lapply(1:stanfit@sim$chains, function(i) {
    x <- stanfit@sim$samples[[i]]
    vars_str <- paste0(vars, collapse = "|")
    xout <- x[grepl(vars_str, names(x))] %>%
      bind_cols()

    if (nrow(xout)) {

      if (inc_warmup) {
        warmup <- 1
      } else {
        warmup <- stanfit@stan_args[[i]]$warmup
      }

      thin <- stanfit@stan_args[[i]]$thin
      iter <- stanfit@stan_args[[i]]$iter

      it <- data.frame(iter = seq(1, iter, by = thin))
      it$n <- 1:nrow(it)

      xout <- xout[it$n[it$iter > warmup], ]
    }
    return(xout)
  }) %>%
    bind_rows()

  if (nrow(val)) {
    pairs(x = val, pch = 19, cex = 0.2, diag.panel = panel.hist, upper.panel = panel.cor, gap = 0)
  }

}


#' @rdname CMfigures
#' @param report List, output of state variables from individual MCMC samples, obtained with `get_report()`
#' @param d List of data variables, obtained with `get_CMdata()`
#' @param year Vector of years
#' @returns
#' - `CM_fit_esc()` returns base graphics with fit to total escapement time series
#' @export
#' @importFrom graphics points polygon
CM_fit_esc <- function(report, d, year) {
  obs <- d$obsescape

  if (missing(year)) year <- 1:length(obs)

  esc <- exp(sapply(report, getElement, "logpredesc")) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975))

  plot(NULL, NULL, type = "n", pch = 16, xlim = range(year), ylim = c(0, 1.1) * range(esc, obs, na.rm = TRUE), xlab = "Year", ylab = "Escapement")
  polygon(c(rev(year), year), c(rev(esc[1, ]), esc[3, ]), col = alpha("grey", 0.5), border = NA)
  lines(year, esc[2, ], lwd = 1.5)

  points(year, obs, pch = 16)
  lines(year, obs, lty = 3)

  invisible()
}

#' @rdname CMfigures
#' @returns
#' - `CM_fit_pHOS()` returns base graphics with fit to pHOS (census) observations
#' @export
CM_fit_pHOS <- function(report, d, year) {
  obs <- d$obs_pHOS

  if (!is.null(obs) && sum(obs, na.rm = TRUE)) {

    if (missing(year)) year <- 1:length(obs)

    pred <- sapply(report, getElement, "pHOScensus") %>%
      apply(1, quantile, probs = c(0.025, 0.5, 0.975))

    plot(NULL, NULL, type = "n", pch = 16, xlim = range(year), ylim = c(0, 1.1) * range(pred, obs, na.rm = TRUE), xlab = "Year", ylab = "pHOS census")
    polygon(c(rev(year), year), c(rev(pred[1, ]), pred[3, ]), col = alpha("grey", 0.5), border = NA)
    lines(year, pred[2, ], lwd = 1.5)

    points(year, obs, pch = 16)
    lines(year, obs, lty = 3)
  }

  invisible()
}

CM_data <- function(obs, year, ylab) {

  if (!is.null(obs) && sum(obs)) {
    if (missing(year)) year <- 1:length(obs)
    plot(year, obs, type = "o", pch = 16, xlim = range(year),
         ylim = c(0, 1.1) * range(obs, na.rm = TRUE), xlab = "Year", ylab = ylab)
  }

  invisible()
}

CM_CWTrel <- function(obs, year1, rs_names) {
  if (missing(rs_names)) rs_names <- seq(1, ncol(obs))

  dat <- obs %>%
    structure(dimnames = list(Year = 1:nrow(obs) + year1 - 1, `Release Strategy` = rs_names)) %>%
    reshape2::melt() %>%
    mutate(`Release Strategy` = factor(`Release Strategy`, rs_names))

  g <- ggplot(dat, aes(.data$Year, .data$value, colour = .data$`Release Strategy`)) +
    geom_point() +
    geom_line() +
    labs(y = "CWT release") +
    theme(legend.position = "bottom")

  if (length(unique(dat$`Release Strategy`)) == 1) {
    g <- g +
      guides(colour = "none") +
      scale_colour_manual(values = "black")
  }

  g
}


#' @rdname CMfigures
#' @param year1 Numeric, first year of model
#' @param rs_names Character vector of hatchery release strategies
#' @returns
#' - `CM_fit_CWTesc()` returns ggplot of fit to CWT escapement at age
#' @export
CM_fit_CWTesc <- function(report, d, year1 = 1, rs_names) {

  if (missing(rs_names)) rs_names <- seq(1, d$n_r)

  ebrood <- sapply(report, getElement, "ebrood", simplify = "array") %>%
    apply(1:3, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    mutate(Age = paste("Age", Var3)) %>%
    mutate(`Release Strategy` = factor(rs_names[Var4], rs_names)) %>%
    reshape2::dcast(Age + Year + `Release Strategy` ~ Var1, value.var = "value")

  dat <- d$cwtesc

  nyears <- dim(dat)[1]
  nage <- dim(dat)[2]
  n_r <- dim(dat)[3]

  cwtesc <- dat %>%
    structure(dimnames = list(Year = year1 + seq(1, nyears) - 1, Age = seq(1, nage), `Release Strategy` = rs_names)) %>%
    reshape2::melt() %>%
    mutate(Age = paste("Age", Age), `Release Strategy` = factor(`Release Strategy`, rs_names)) %>%
    dplyr::filter(value > 0)

  g <- ggplot(ebrood, aes(.data$Year, colour = .data$`Release Strategy`, fill = .data$`Release Strategy`)) +
    geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), colour = NA, alpha = ifelse(n_r > 1, 0.25, 0.5)) +
    geom_line(aes(y = .data$`50%`), linewidth = 0.75) +
    facet_wrap(vars(Age), scales = "free_y") +
    geom_point(data = cwtesc, aes(y = .data$value)) +
    geom_line(data = cwtesc, aes(y = .data$value), linetype = 3) +
    labs(x = "Brood year", y = "CWT Escapement") +
    theme(legend.position = "bottom")

  if (n_r == 1) {
    g <- g +
      guides(colour = "none", fill = "none") +
      scale_colour_manual(values = "black") +
      scale_fill_manual(values = "grey")
  }

  g
}

#' @rdname CMfigures
#' @param PT Logical, whether to plot preterminal catch, otherwise (plot terminal catch)
#' @returns
#' - `CM_fit_CWTcatch()` returns ggplot of fit to CWT catch at age
#' @export
CM_fit_CWTcatch <- function(report, d, PT = TRUE, year1 = 1, rs_names) {
  dat <- d[[ifelse(PT, "cwtcatPT", "cwtcatT")]]

  if (sum(dat)) {

    nyears <- dim(dat)[1]
    nage <- dim(dat)[2]
    n_r <- dim(dat)[3]

    if (missing(rs_names)) rs_names <- seq(1, n_r)

    cwtcat <- dat %>%
      structure(dimnames = list(Year = year1 + seq(1, nyears) - 1, Age = seq(1, nage), `Release Strategy` = rs_names)) %>%
      reshape2::melt() %>%
      mutate(Age = paste("Age", Age), `Release Strategy` = factor(`Release Strategy`, rs_names)) %>%
      dplyr::filter(value > 0)

    cbrood <- sapply(report, getElement, ifelse(PT, "cbroodPT", "cbroodT"), simplify = "array") %>%
      apply(1:3, quantile, probs = c(0.025, 0.5, 0.975)) %>%
      reshape2::melt() %>%
      mutate(Year = Var2 + year1 - 1) %>%
      mutate(Age = paste("Age", Var3)) %>%
      mutate(`Release Strategy` = factor(rs_names[Var4], rs_names)) %>%
      reshape2::dcast(Age + Year + `Release Strategy` ~ Var1, value.var = "value")

    g <- ggplot(cbrood, aes(.data$Year, colour = .data$`Release Strategy`, fill = .data$`Release Strategy`)) +
      geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), colour = NA, alpha = ifelse(n_r > 1, 0.25, 0.5)) +
      geom_line(aes(y = .data$`50%`), linewidth = 0.75) +
      facet_wrap(vars(Age), scales = "free_y") +
      geom_point(data = cwtcat, aes(y = .data$value)) +
      geom_line(data = cwtcat, aes(y = .data$value), linetype = 3) +
      labs(x = "Brood year", y = ifelse(PT, "CWT preterminal catch", "CWT terminal catch")) +
      theme(legend.position = "bottom")

    if (n_r == 1) {
      g <- g +
        guides(colour = "none", fill = "none") +
        scale_colour_manual(values = "black") +
        scale_fill_manual(values = "grey")
    }

    g
  }
}

#' @rdname CMfigures
#' @param r Integer, the release strategy for the figure (only if `annual = FALSE`)
#' @param brood Logical, whether to show results by brood year or return year (FALSE)
#' @param annual Logical, whether to show panel figure by individual year (TRUE) or a single time series figure
#' @returns
#' - `CM_maturity()` returns ggplot of estimated maturity at age
#' @export
#' @importFrom dplyr rename
CM_maturity <- function(report, d, year1 = 1, r = 1, brood = TRUE, annual = FALSE, rs_names) {
  n_r <- d$n_r
  if (missing(rs_names)) rs_names <- seq(1, n_r)

  if (annual) {

    if (brood) {
      matt <- sapply(report, function(i) {
        sapply(1:n_r, function(rr) CY2BY(i[["matt"]][, , rr]), simplify = 'array')
      }, simplify = 'array')
    } else {
      matt <- sapply(report, getElement, "matt", simplify = 'array')
    }
    matt_q <- apply(matt, 1:3, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
      reshape2::melt() %>%
      mutate(Year = Var2 + year1 - 1) %>%
      rename(Age = Var3) %>%
      mutate(`Release Strategy` = factor(rs_names[Var4], rs_names)) %>%
      reshape2::dcast(Age + Var2 + Year + `Release Strategy` ~ Var1) %>%
      filter(!is.na(`50%`))

    g <- matt_q %>%
      ggplot(aes(Age, .data$`50%`, colour = .data$`Release Strategy`, fill = .data$`Release Strategy`)) +
      geom_line() +
      geom_point() +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), colour = NA, alpha = ifelse(n_r > 1, 0.25, 0.5)) +
      facet_wrap(vars(Year)) +
      theme(
        panel.spacing = unit(0, "in"),
        legend.position = "bottom",
        strip.background = element_blank()
      ) +
      labs(x = "Age", y = "Proportion mature", colour = "Release strategy", fill = "Release strategy",
           title = ifelse(brood, "Brood year", "Return year"))

    if (n_r == 1) {
      g <- g +
        guides(colour = "none", fill = "none") +
        scale_colour_manual(values = "black") +
        scale_fill_manual(values = "grey")
    }
  } else {

    bmatt <- data.frame(Age = 1:d$Nages, value = d$bmatt) %>%
      dplyr::filter(Age > 1)

    if (brood) {
      matt <- sapply(report, function(i) CY2BY(i[["matt"]][, , r]), simplify = 'array')
    } else {
      matt <- sapply(report, function(i) i[["matt"]][, , r], simplify = 'array')
    }
    matt_q <- apply(matt, 1:2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
      reshape2::melt() %>%
      mutate(Year = Var2 + year1 - 1) %>%
      rename(Age = Var3) %>%
      reshape2::dcast(Age + Var2 + Year ~ Var1) %>%
      filter(!is.na(`50%`))

    g <- matt_q %>%
      dplyr::filter(Age > 1) %>%
      ggplot(aes(Year, .data$`50%`, fill = factor(.data$Age), colour = factor(.data$Age))) +
      geom_line() +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
      geom_hline(data = bmatt, linetype = 2, aes(yintercept = .data$value, colour = factor(.data$Age))) +
      labs(x = ifelse(brood, "Brood year", "Return year"), y = "Proportion mature", colour = "Age", fill = "Age")

    if (n_r > 1) {
      g <- g +
        ggtitle(paste("Release strategy:", rs_names[r]))

    }
  }

  g
}


#' @rdname CMfigures
#' @param type Character, indicates type of variable to plot
#' @returns
#' - `CM_vul()` returns ggplot of estimated fishery vulnerability at age
#' @export
#' @importFrom dplyr rename
CM_vul <- function(report, type = c("vulPT", "vulT")) {
  type <- match.arg(type)
  vul <- sapply(report, getElement, type) # age x sim

  if (sum(vul)) {
    vul_q <- apply(vul, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
      reshape2::melt() %>%
      rename(Age = Var2) %>%
      reshape2::dcast(list("Age", "Var1"), value.var = "value")

    g <- vul_q %>%
      ggplot(aes(.data$Age, .data$`50%`)) +
      geom_line() +
      geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), alpha = 0.2) +
      labs(x = "Age", y = ifelse(type == "vulPT", "Preterminal vulnerability", "Terminal vulnerability"))
    g
  }
}

#' @rdname CMfigures
#' @returns
#' - `CM_SRR()` returns ggplot of estimated stock-recruit relationship (density-dependent juvenile production from egg production) with average
#' relationship and realized annual values. Years correspond to years of egg production (assumes juvenile production for the following brood year).
#' @export
CM_SRR <- function(report, year1 = 1) {
  egg <- sapply(report, getElement, "egg") %>% apply(1, median)
  smolt <- sapply(report, function(x) x$N[-1, 1, 1]) %>% apply(1, median)
  epred <- seq(0, 1.1 * max(egg), length.out = 50)
  spred <- sapply(report, function(x) x$alpha * epred * exp(-x$beta * epred)) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

  year <- year1 + seq(1, length(egg)) - 1

  df <- data.frame(
    year = year,
    egg = egg,
    smolt = smolt
  )

  df_med <- data.frame(
    egg = epred,
    smolt = spred[2, ]
  )

  df_poly1 <- data.frame(
    egg = epred,
    smolt = spred[1, ]
  )

  df_poly2 <- data.frame(
    egg = rev(epred),
    smolt = rev(spred[3, ])
  )

  g <- ggplot(df, aes(.data$egg, .data$smolt)) +
    geom_point(shape = 1) +
    geom_line(data = df_med) +
    geom_polygon(data = rbind(df_poly1, df_poly2), fill = "grey", alpha = 0.5) +
    labs(x = "Egg production", y = "Smolt production") +
    expand_limits(x = 0, y = 0)

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    g <- g + ggrepel::geom_text_repel(aes(label = .data$year))
  }
  g
}

.CM_prod <- function(report, d) {
  sapply(1:length(report), function(x) {
    sapply(1:d$Ldyr, function(y) {

      mo <- report[[x]]$mo[y, ]
      matt <- report[[x]]$matt[y, , d$r_matt]

      lo <- numeric(d$Nages)
      lo[1] <- 1
      for (a in 2:d$Nages) {
        lo[a] <- lo[a-1] * exp(-mo[a-1]) * (1 - matt[a-1]) # unfished juvenile survival
      }
      epro <- sum(lo * d$ssum * d$fec * matt)

      return(report[[x]]$alpha * epro)
    })
  })
}

#' @rdname CMfigures
#' @returns
#' - `CM_prod()` returns ggplot of realized productivity in the absence of fishery harvest, annual values are based on natural mortality and maturity at age
#' @export
CM_prod <- function(report, d, year1 = 1) {

  prod <- .CM_prod(report, d)

  prod_q <- apply(prod, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1"))

  g <- prod_q %>%
    ggplot(aes(Year, .data$`50%`)) +
    geom_line() +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
    labs(x = "Year", y = "Productivity")
  g

}

.CM_Srep <- function(report, d, type = c("spawner", "egg")) {
  type <- match.arg(type)

  sapply(1:length(report), function(x) {
    sapply(1:d$Ldyr, function(y) {

      mo <- report[[x]]$mo[y, ]
      matt <- report[[x]]$matt[y, , d$r_matt]

      lo <- numeric(d$Nages)
      lo[1] <- 1
      for (a in 2:d$Nages) {
        lo[a] <- lo[a-1] * exp(-mo[a-1]) * (1 - matt[a-1]) # unfished juvenile survival
      }
      epro <- sum(lo * d$ssum * d$fec * matt)
      spro <- sum(lo * d$ssum * matt)

      beta_s <- report[[x]]$beta * epro / spro
      alpha_s <- report[[x]]$alpha * epro

      Srep <- log(alpha_s)/beta_s

      if (type == "spawner") {
        return(Srep)
      } else {
        Erep <- Srep * epro / spro
        return(Erep)
      }
    })
  })
}

#' @rdname CMfigures
#' @returns
#' - `CM_Srep()` returns ggplot of realized spawner or egg production at replacement, annual values are based on natural mortality and maturity at age
#' @export
CM_Srep <- function(report, d, year1 = 1, type = c("spawner", "egg")) {
  type <- match.arg(type)

  Srep <- .CM_Srep(report, d, type)

  Srep_q <- apply(Srep, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1"))

  g <- Srep_q %>%
    ggplot(aes(Year, .data$`50%`)) +
    geom_line() +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
    labs(x = "Year", y = ifelse(type == "spawner", expression(S[rep]), expression(E[rep])))

  g
}


.CM_SMSY <- function(report, d, type = c("spawner", "egg")) {
  type <- match.arg(type)

  sapply(1:length(report), function(x) {
    sapply(1:d$Ldyr, function(y) {

      mo <- report[[x]]$mo[y, ]
      matt <- report[[x]]$matt[y, , d$r_matt]
      vulT <- report[[x]]$vulT
      FT <- report[[x]]$FT

      lo <- numeric(d$Nages)
      lo[1] <- 1
      for (a in 2:d$Nages) {
        lo[a] <- lo[a-1] * exp(-mo[a-1]) * (1 - matt[a-1]) # unfished juvenile survival
      }
      epro <- sum(lo * d$ssum * d$fec * matt)
      alpha_s <- report[[x]]$alpha * epro

      if (!sum(FT) && type == "spawner") {
        spro <- sum(lo * d$ssum * matt)
        beta_s <- report[[x]]$beta * epro / spro
        if (alpha_s > 1) {
          val <- calc_Smsy_Ricker(log(alpha_s), beta_s)
        } else {
          val <- 0
        }
      } else {
        vulPT <- report[[x]]$vulPT
        if (!sum(FT)) vulT[] <- 1
        SRRpars <- data.frame("kappa" = alpha_s, "Smax" = 1/report[[x]]$beta, "phi" = epro, "SRrel" = "Ricker")
        ref <- calc_MSY(Mjuv = c(mo, length(mo)), fec = d$fec, p_female = d$ssum, rel_F = c(0, 1), vulPT = vulPT, vulT = vulT, p_mature = matt,
                        s_enroute = 1, SRRpars = SRRpars)

        if (type == "spawner") {
          val <- ref["Spawners_MSY"]
        } else {
          val <- ref["Egg_MSY"]
        }
      }
      as.numeric(val)
    })
  })

}

CM_SMSY <- function(report, d, year1 = 1, type = c("spawner", "egg")) {
  type <- match.arg(type)

  SMSY <- .CM_SMSY(report, d, type)

  SMSY_q <- apply(SMSY, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1"))

  g <- SMSY_q %>%
    ggplot(aes(Year, .data$`50%`)) +
    geom_line() +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
    labs(x = "Year", y = ifelse(type == "spawner", expression(S[MSY]), expression(E[MSY])))

  g
}

.CM_Sgen <- function(report, d) {

  sapply(1:length(report), function(x) {
    sapply(1:d$Ldyr, function(y) {

      mo <- report[[x]]$mo[y, ]
      matt <- report[[x]]$matt[y, , d$r_matt]
      vulT <- report[[x]]$vulT
      FT <- report[[x]]$FT

      lo <- numeric(d$Nages)
      lo[1] <- 1
      for (a in 2:d$Nages) {
        lo[a] <- lo[a-1] * exp(-mo[a-1]) * (1 - matt[a-1]) # unfished juvenile survival
      }
      epro <- sum(lo * d$ssum * d$fec * matt)
      alpha_s <- report[[x]]$alpha * epro

      if (!sum(FT)) {
        spro <- sum(lo * d$ssum * matt)
        beta_s <- report[[x]]$beta * epro / spro
        if (alpha_s > 1) {
          val <- calc_Sgen_Ricker(log(alpha_s), beta_s)
        } else {
          val <- 0
        }
      } else {
        vulPT <- report[[x]]$vulPT
        SRRpars <- data.frame("kappa" = alpha_s, "Smax" = 1/report[[x]]$beta, "phi" = epro, "SRrel" = "Ricker")
        ref <- calc_MSY(Mjuv = c(mo, length(mo)), fec = d$fec, p_female = d$ssum, rel_F = c(0, 1), vulPT = vulPT, vulT = vulT, p_mature = matt,
                        s_enroute = 1, SRRpars = SRRpars)
        SMSY <- ref["Spawners_MSY"]

        val <- calc_Sgen(Mjuv = c(mo, length(mo)), fec = d$fec, p_female = d$ssum, rel_F = c(0, 1), vulPT = vulPT, vulT = vulT, p_mature = matt,
                         s_enroute = 1, SRRpars = SRRpars, SMSY = SMSY)
      }
      as.numeric(val)
    })
  })

}

CM_Sgen <- function(report, d, year1 = 1) {
  Sgen <- .CM_Sgen(report, d)

  Sgen_q <- apply(Sgen, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1"))

  g <- Sgen_q %>%
    ggplot(aes(Year, .data$`50%`)) +
    geom_line() +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
    labs(x = "Year", y = expression(S[gen]))

  g
}


.CM_statevarage <- function(report, year1 = 1, ci = TRUE, var, ylab, xlab = "Year", scales = "free_y") {
  arr <- sapply(report, getElement, var, simplify = "array")

  if (sum(arr, na.rm = TRUE)) {
    if (length(dim(arr)) == 4) { # Include NO/HO dimension

      df <- arr %>%
        apply(1:3, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
        reshape2::melt() %>%
        mutate(Year = Var2 + year1 - 1, Origin = ifelse(Var4 == 1, "Natural", "Hatchery")) %>%
        mutate(Age = paste("Age", Var3)) %>%
        reshape2::dcast(Year + Age + Origin ~ Var1, value.var = "value")

      g <- ggplot(df, aes(.data$Year, .data$`50%`, fill = .data$Origin, colour = .data$Origin)) +
        geom_line() +
        facet_wrap(vars(Age), scales = scales) +
        labs(x = xlab, y = ylab) +
        expand_limits(y = 0)

    } else {

      df <- arr %>%
        apply(1:2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
        reshape2::melt() %>%
        mutate(Year = Var2 + year1 - 1) %>%
        mutate(Age = paste("Age", Var3)) %>%
        reshape2::dcast(Year + Age ~ Var1, value.var = "value")

      g <- ggplot(df, aes(.data$Year, .data$`50%`)) +
        geom_line() +
        facet_wrap(vars(Age), scales = scales) +
        labs(x = xlab, y = ylab) +
        expand_limits(y = 0)
    }

    if (ci) g <- g + geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), alpha = 0.2)

    return(g)
  }
  invisible()
}

CM_ts_origin <- function(report, year1 = 1, ci = TRUE, var = "Spawners", ylab = var, xlab = "Year") {
  var <- match.arg(var)
  if (var == "Spawners") {
    arr <- sapply(report, getElement, "syear", simplify = "array") %>%
      apply(c(1, 3, 4), sum)
  }

  df <- arr %>%
    apply(1:2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1, Origin = ifelse(Var3 == 1, "Natural", "Hatchery")) %>%
    reshape2::dcast(Year + Origin ~ Var1, value.var = "value")

  g <- ggplot(df, aes(.data$Year, .data$`50%`, colour = .data$Origin, fill = .data$Origin)) +
    geom_line() +
    labs(x = xlab, y = ylab) +
    expand_limits(y = 0)

  if (ci) g <- g + geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), alpha = 0.2)
  g
}


#' @rdname CMfigures
#' @param ci Logical whether to show posterior intervals in addition to posterior median
#' @returns
#' - `CM_M()` returns ggplot of estimated natural mortality time series by age (marine stage)
#' @export
CM_M <- function(report, year1 = 1, ci = TRUE) {
  df <- sapply(report, getElement, "mo", simplify = "array") %>%
    apply(1:2, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    mutate(Age = paste("Age", Var3)) %>%
    reshape2::dcast(Year + Age ~ Var1, value.var = "value")

  g <- ggplot(df, aes(.data$Year, .data$`50%`)) +
    geom_line() +
    facet_wrap(vars(Age), scales = "free_y") +
    labs(x = "Year", y = "Natural mortality") +
    expand_limits(y = 0)

  if (ci) g <- g + geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), fill = alpha("grey", 0.5))
  g
}


#' @rdname CMfigures
#' @param surv Logical, whether to plot survival (values between 0 - 1) or instantaneous mortality rates
#' @returns
#' - `CM_Megg()` returns ggplot of egg-juvenile mortality time series
#' @export
CM_Megg <- function(report, year1 = 1, ci = TRUE, surv = FALSE) {
  megg <- sapply(report, getElement, 'megg')

  if (surv) megg <- exp(-megg)

  df <- megg %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1"), value.var = "value")

  g <- ggplot(df, aes(.data$Year, .data$`50%`)) +
    geom_line() +
    labs(x = "Year", y = ifelse(surv, "Egg-smolt survival", "Egg-smolt mortality")) +
    expand_limits(y = 0)

  if (ci) g <- g + geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), fill = alpha("grey", 0.5))
  g
}


#' @rdname CMfigures
#' @returns
#' - `CM_Njuv()` returns ggplot of juvenile abundance
#' @export
CM_Njuv <- function(report, year1 = 1, ci = TRUE) {
  .CM_statevarage(report, year1, ci, "N", "Juvenile abundance") +
    theme(legend.position = "bottom")
}

#' @rdname CMfigures
#' @returns
#' - `CM_recr()` returns ggplot of recruitment (mature return)
#' @export
CM_recr <- function(report, year1 = 1, ci = TRUE) {
  .CM_statevarage(report, year1, ci, "recr", "Recruitment") +
    theme(legend.position = "bottom")
}

#' @rdname CMfigures
#' @returns
#' - `CM_esc()` returns ggplot of escapement (after terminal harvest)
#' @export
CM_esc <- function(report, year1 = 1, ci = TRUE) {
  .CM_statevarage(report, year1, ci, "escyear", "Escapement") +
    theme(legend.position = "bottom")
}

#' @rdname CMfigures
#' @returns
#' - `CM_F()` returns ggplot of instantaneous fishing mortality
#' @export
CM_F <- function(report, PT = TRUE, year1 = 1, ci = TRUE) {
  .CM_ts(
    report, year1, ci,
    var = ifelse(PT, "FPT", "FT"),
    ylab = ifelse(PT, "Preterminal fishing mortality", "Terminal fishing mortality")
  )
}

.CM_ts <- function(report, year1 = 1, ci = TRUE, var, ylab, xlab = "Year") {
  ts <- sapply(report, getElement, var) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1"), value.var = "value")

  if (sum(ts$`50%`, na.rm = TRUE)) {
    g <- ggplot(ts, aes(.data$Year, .data$`50%`)) +
      geom_line() +
      labs(x = xlab, y = ylab) +
      expand_limits(y = 0)

    if (ci) g <- g + geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), fill = alpha("grey", 0.5))
    return(g)
  }
  invisible()
}

#' @rdname CMfigures
#' @returns
#' - `CM_surv()` returns ggplot of natural survival (converting from instantaneous units of natural mortality)
#' @export
CM_surv <- function(report, year1 = 1, ci = TRUE) {
  df <- exp(-sapply(report, getElement, "mo", simplify = "array")) %>%
    apply(1:2, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    mutate(Age = paste("Age", Var3)) %>%
    reshape2::dcast(Year + Age ~ Var1, value.var = "value")

  g <- ggplot(df, aes(.data$Year, .data$`50%`)) +
    geom_line() +
    facet_wrap(vars(Age), scales = "free_y") +
    labs(x = "Year", y = "Natural survival") +
    expand_limits(y = 0)

  if (ci) g <- g + geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), fill = alpha("grey", 0.5))
  g
}


#' @rdname CMfigures
#' @returns
#' - `CM_wt()` returns ggplot of annual deviations in egg-juvenile mortality from the Ricker function
#' @export
CM_wt <- function(stanfit, year1 = 1, ci = TRUE) {

  wt <- try(rstan::extract(stanfit, "wt")$wt, silent = TRUE)

  if (!is.character(wt)) {
    ts <- wt %>%
      apply(2, quantile, probs = c(0.025, 0.5, 0.975)) %>%
      reshape2::melt() %>%
      mutate(Year = Var2 + year1 - 1) %>%
      reshape2::dcast(list("Year", "Var1"), value.var = "value")

    g <- ggplot(ts, aes(.data$Year, .data$`50%`)) +
      geom_line() +
      geom_point() +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(x = "Year", y = "Egg mortality deviation")

    if (ci) g <- g + geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), fill = alpha("grey", 0.5))
    g
  }

}

#' @rdname CMfigures
#' @returns
#' - `CM_surv2()` returns ggplot of annual survival to age 2, which includes age-1 mortality (marine life stage) for both natural and hatchery origin fish.
#' Hatchery fish experience additional mortality specified by release mortality.
#' @importFrom dplyr n
#' @export
CM_surv2 <- function(report, year1 = 1, ci = TRUE, ylab = "Survival to age 2") {

  NO <- exp(-sapply(report, function(x) x$mo[, 1])) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    t() %>%
    as.data.frame() %>%
    mutate(type = "Natural", Year = year1 + 1:n() - 1)

  HO <- exp(-sapply(report, getElement, "moplot")) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    t() %>%
    as.data.frame() %>%
    mutate(type = "Hatchery", Year = year1 + 1:n() - 1)

  df <- rbind(NO, HO)

  g <- ggplot(df, aes(.data$Year, .data$`50%`)) +
    geom_line(aes(colour = .data$type)) +
    labs(x = "Year", y = ylab, fill = "Origin", colour = "Origin") +
    expand_limits(y = 0) +
    theme(legend.position = "bottom")

  if (ci) g <- g + geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`, fill = .data$type), alpha = 0.25)
  g
}

#' @rdname CMfigures
#' @returns
#' - `CM_wt()` returns ggplot of annual deviations in age 1 natural mortality (first year in marine life stage, deviations from time series average)
#' @export
CM_wto <- function(stanfit, year1 = 1, ci = TRUE) {

  wto <- try(rstan::extract(stanfit, "wto")$wto, silent = TRUE)

  if (!is.character(wto)) {
    ts <- wto %>%
      apply(2, quantile, probs = c(0.025, 0.5, 0.975)) %>%
      reshape2::melt() %>%
      mutate(Year = Var2 + year1 - 1) %>%
      reshape2::dcast(list("Year", "Var1"), value.var = "value")

    g <- ggplot(ts, aes(.data$Year, .data$`50%`)) +
      geom_line() +
      geom_point() +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(x = "Year", y = "Age 1 mortality deviation")

    if (ci) g <- g + geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), fill = alpha("grey", 0.5))
    g
  }

}

#CM_fitness <- function(report, year) {
#  .CM_ts(
#    report, year1, ci,
#    var = "fitness",
#    ylab = "Fitness"
#  )
#}

# Convert calendar year matrix to brood year matrix
CY2BY <- function(x) {

  nt <- nrow(x)
  na <- ncol(x)

  xbrood <- matrix(NA_real_, nt, na)
  for (t in 1:nt) {
    for (a in 1:na) {
      if (t+a-1 <= nt) {
        xbrood[t, a] <- x[t+a-1, a]
      }
    }
  }
  return(xbrood)
}

# Calculate adult equivalent
calc_AEQ <- function(report, brood = TRUE) {

  nt <- dim(report$matt)[1]
  na <- dim(report$matt)[2]
  n_r <- dim(report$matt)[3]

  Msurv <- array(NA_real_, c(nt, na))
  if (brood) {
    matt <- sapply(1:n_r, function(r) CY2BY(report$matt[, , r]), simplify = "array")
    Msurv[, -na] <- CY2BY(report$mo)
  } else {
    matt <- report$matt
    Msurv[, -na] <- report$mo
  }
  Msurv[, na] <- 1

  AEQ <- array(NA_real_, dim(matt))
  AEQ[, na, ] <- 1

  for (t in seq(nt, 2) - 1) {
    for (a in seq(na, 2) - 1) {
      AEQ[t, a, ] <- matt[t, a, ] + (1 - matt[t, a, ]) * Msurv[t, a+1] * AEQ[t, a+1, ]
    }
  }

  return(AEQ)
}


#' @rdname CMfigures
#' @param at_age Logical, whether to make figure by individual age
#' @returns
#' - `CM_ER()` returns ggplot of exploitation rate either by individual age or aggregate values using adult equivalents
#' @export
CM_ER <- function(report, brood = TRUE, type = c("PT", "T", "all"), year1 = 1, ci = TRUE, at_age = TRUE, r = 1) {

  type <- match.arg(type)

  if (at_age) {

    ER_list <- lapply(report, function(i) {
      name <- ifelse(type == "PT", "survPT", "survT")
      CYER <- 1 - i[[name]]
      if (brood) {
        ER <- CY2BY(CYER)
      } else {
        ER <- CYER
      }
      return(list(ER = ER))
    })

    g <- .CM_statevarage(
      ER_list,
      year1,
      ci,
      "ER",
      ylab = switch(
        type,
        "PT" = "Preterminal exploitation rate",
        "T" = "Terminal exploitation rate"
      ),
      xlab = ifelse(brood, "Brood year", "Return year"),
      scales = "fixed"
    )

  } else {

    ER_list <- lapply(report, function(i) {

      esc <- apply(i$escyear, 1:2, sum)
      morts_PT <- apply(i$cyearPT, 1:2, sum)
      morts_T <- apply(i$cyearT, 1:2, sum)
      AEQ_PT <- calc_AEQ(i, brood = brood)[, , r]
      AEQ_T <- array(1, dim(esc))

      if (brood) {
        esc <- CY2BY(esc)
        morts_PT <- CY2BY(morts_PT)
        morts_T <- CY2BY(morts_T)
      }

      denom <- rowSums(morts_PT * AEQ_PT + morts_T * AEQ_T + esc)

      if (type == "PT") {
        num <- rowSums(morts_PT * AEQ_PT)
      } else if (type == "T") {
        num <- rowSums(morts_T * AEQ_T)
      } else {
        num <- rowSums(morts_PT * AEQ_PT + morts_T * AEQ_T)
      }

      list(ER = num/denom)
    })

    g <- .CM_ts(
      ER_list, year1, ci,
      var = "ER",
      xlab = ifelse(brood, "Brood year", "Return year"),
      ylab = switch(
        type,
        "PT" = "Preterminal exploitation rate",
        "T" = "Terminal exploitation rate",
        "all" = "Total exploitation rate"
      )
    )
  }

  g
}

#' @rdname CMfigures
#' @returns
#' - `CM_CWT_ER()` returns ggplot of CWT exploitation rate (by release strategy)
#' @export
#' @importFrom dplyr left_join
CM_CWT_ER <- function(report, brood = TRUE, type = c("PT", "T", "all"), year1 = 1, ci = TRUE, rs_names) {

  type <- match.arg(type)

  ER <- sapply(report, function(i) {

    if (brood) {
      esc <- i$ebrood
      morts_PT <- i$cbroodPT
      morts_T <- i$cbroodT
    } else {
      esc <- i$ecwt
      morts_PT <- i$ccwtPT
      morts_T <- i$ccwtT
    }

    AEQ_PT <- calc_AEQ(i, brood = brood)
    AEQ_T <- array(1, dim(esc))

    denom <- apply(morts_PT * AEQ_PT + morts_T * AEQ_T + esc, c(1, 3), sum) # year x rs

    if (type == "PT") {
      num <- apply(morts_PT * AEQ_PT, c(1, 3), sum)
    } else if (type == "T") {
      num <- apply(morts_T * AEQ_T, c(1, 3), sum)
    } else {
      num <- apply(morts_PT * AEQ_PT + morts_T * AEQ_T, c(1, 3), sum)
    }

    return(num/denom)
  }, simplify = "array")

  if (missing(rs_names)) rs_names <- seq(1, dim(ER)[2])
  n_r <- length(rs_names)

  ts <- apply(ER, 1:2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1, `Release Strategy` = factor(rs_names[Var3], rs_names)) %>%
    reshape2::dcast(Year + `Release Strategy` ~ Var1, value.var = "value")

  if (sum(ts$`50%`, na.rm = TRUE)) {

    has_rel <- sapply(1:n_r, function(r) {
      rowSums(CY2BY(report[[1]]$ecwt[, , r]), na.rm = TRUE) > 0
    }, simplify = "array") %>%
      as.data.frame()
    colnames(has_rel) <- rs_names
    has_rel$Year <- year1 + seq(1, nrow(has_rel)) - 1

    has_rel <- has_rel %>%
      reshape2::melt(id.var = "Year", variable.name = "Release Strategy") %>%
      filter(value == TRUE)

    g <- left_join(has_rel, ts, by = c("Year", "Release Strategy")) %>%
      ggplot(aes(.data$Year, .data$`50%`, colour = .data$`Release Strategy`)) +
      geom_line() +
      geom_point() +
      labs(x = ifelse(brood, "Brood year", "Return year"),
           y = switch(
             type,
             "PT" = "Preterminal exploitation rate",
             "T" = "Terminal exploitation rate",
             "all" = "Total exploitation rate"
           )) +
      expand_limits(y = 0) +
      theme(legend.position = "bottom")

    if (ci) {
      g <- g +
        geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`, fill = `Release Strategy`),
                    alpha = ifelse(n_r > 1, 0.25, 0.5))
    }

    if (n_r == 1) {
      g <- g +
        guides(colour = "none", fill = "none") +
        scale_colour_manual(values = "black") +
        scale_fill_manual(values = "grey")
    }
  }

  g
}


#' @rdname CMfigures
#' @param x Matrix of covariates by year x covariate
#' @param names Character of covariate names
#' @param b Matrix of fixed effect coefficients by simulation x covariate. If missing only the covariates (`x`) are plotted,
#' otherwise, the dot product `sum(x * b)` is calculated by individual simulation and quantiles are plotted
#' @param ylab Character y axis label
#' @returns
#' - `CM_covariate()` returns ggplot of mortality covariates
#' @export
#' @importFrom dplyr left_join
CM_covariate <- function(x, names, year1 = 1, b, ylab = "Covariate") {

  if (sum(x)) {

    if (!is.matrix(x)) x <- matrix(x, ncol = 1)
    if (missing(names)) names <- paste("Covariate", 1:ncol(x))

    if (missing(b)) { # Plot covariate

      g <- structure(x, dimnames = list(Year = 1:nrow(x) + year1 - 1, Covariate = names)) %>%
        reshape2::melt() %>%
        ggplot(aes(.data$Year, .data$value, colour = .data$Covariate)) +
        geom_line() +
        geom_point() +
        labs(y = ylab, colour = NULL) +
        expand_limits(y = 0) +
        theme(legend.position = "bottom")

    } else { # Plot M = covariate x coefficient

      if (!is.matrix(b)) b <- matrix(b, ncol = 1) # simulation x cov
      nsim <- nrow(b)
      ncov <- ncol(b)
      nt <- nrow(x)

      Mcov <- array(NA_real_, c(nsim, ncov, nt))

      i <- expand.grid(
        x = 1:nsim,
        j = 1:ncov,
        t = 1:nt
      ) %>% as.matrix()

      i_b <- i[, c("x", "j")]
      i_x <- i[, c("t", "j")]

      Mcov[i] <- b[i_b] * x[i_x]

      g <- Mcov %>%
        apply(2:3, quantile, c(0.025, 0.5, 0.975)) %>%
        reshape2::melt() %>%
        mutate(Covariate = names[Var2], Year = Var3 + year1 - 1) %>%
        reshape2::dcast(Year + Covariate ~ Var1, value.var = "value") %>%
        ggplot(aes(.data$Year, .data$`50%`, colour = .data$Covariate, fill = .data$Covariate)) +
        geom_line() +
        geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), alpha = 0.2) +
        labs(y = ylab, colour = NULL, fill = NULL) +
        expand_limits(y = 0) +
        theme(legend.position = "bottom")
    }

    return(g)

  }

}

#' Conditioning model markdown report
#'
#' Generate a markdown report to plot time series and MCMC posteriors of estimates from the conditioning model. See [get_report()] for the
#' various plotting functions used in the report.
#'
#' @param stanfit Output from [sample_CM()]
#' @param year Optional vector of calendar years
#' @param cov1_names Optional character vector for names of covariates that predict age-1 natural mortality
#' @param cov_names Optional character vector for names of covariates that predict age-2+ natural mortality
#' @param rs_names Optional character vector for names of hatchery release strategies
#' @param name Optional character string for the model name to include in the report, e.g., model run number
#' @param filename Character string for the name of the markdown and HTML files
#' @param dir The directory in which the markdown and HTML files will be saved.
#' @param open_file Logical, whether the HTML document is opened after it is rendered
#' @param render_args List of arguments to pass to [rmarkdown::render()]
#' @param ... Additional arguments (not used)
#' @details Report excludes MCMC values from warmup iterations
#' @returns Returns invisibly the output of [rmarkdown::render()], typically the path of the output file
#' @seealso [fit_CM()] [get_report()]
#' @export
report_CM <- function(stanfit, year, cov1_names, cov_names, rs_names,
                      name, filename = "CM", dir = tempdir(), open_file = TRUE, render_args = list(), ...) {

  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    warning("Install ggrepel package to label years for stock-recruit figure.")
  }

  report <- get_report(stanfit, inc_warmup = FALSE)
  fit <- stanfit@.MISC$CMfit
  d <- get_CMdata(fit)

  if (missing(year)) year <- 1:d$Ldyr
  if (missing(name)) name <- "CM fit"
  year1 <- year[1]

  ####### Function arguments for rmarkdown::render
  rmd <- system.file("include", "CMreport.Rmd", package = "salmonMSE") %>% readLines()
  rmd_split <- split(rmd, 1:length(rmd))

  name_ind <- grep("NAME", rmd)
  rmd_split[[name_ind]] <- paste("#", name, "{.tabset}")

  filename_rmd <- paste0(filename, ".Rmd")

  render_args$input <- file.path(dir, filename_rmd)
  if (is.null(render_args$quiet)) render_args$quiet <- TRUE

  # Generate markdown report
  if (!dir.exists(dir)) {
    message("Creating directory: ", dir)
    dir.create(dir)
  }
  write(unlist(rmd_split), file = file.path(dir, filename_rmd))

  # Rendering markdown file
  message("Rendering markdown file: ", file.path(dir, filename_rmd))
  output_filename <- do.call(rmarkdown::render, render_args)
  message("Rendered file: ", output_filename)

  if (open_file) browseURL(output_filename)
  invisible(output_filename)
}
