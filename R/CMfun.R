
#' @importFrom dplyr bind_rows bind_cols mutate
stan_trace <- function(stanfit, vars) {
  val <- lapply(1:stanfit@sim$chains, function(i) {
    x <- stanfit@sim$samples[[i]]
    vars_regex <- paste0(vars, collapse = "|")

    x_pars <- x[grepl(vars_regex, names(x))]

    if (length(x_pars)) {
      x[grepl(vars_regex, names(x))] %>%
        bind_cols() %>%
        as.data.frame() %>%
        mutate(Chain = i, Iteration = 1:nrow(.)) %>%
        reshape2::melt(id.vars = c("Chain", "Iteration")) %>%
        dplyr::filter(.data$Iteration > stanfit@sim$warmup)
    } else {
      data.frame()
    }

  }) %>%
    bind_rows()

  if (nrow(val)) {
    g <- ggplot(val, aes(Iteration, value, colour = factor(Chain))) +
      geom_line() +
      facet_wrap(vars(variable), scales = "free_y") +
      labs(colour = "Chain")
    g
  }
}

#' @importFrom stats cor
#' @importFrom graphics pairs rect strwidth text
pairs_panel <- function(stanfit, vars) {
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
    xout[-seq(1, stanfit@sim$warmup), ]
  }) %>%
    bind_rows()

  if (nrow(val)) {
    pairs(x = val, pch = 19, cex = 0.2, diag.panel = panel.hist, upper.panel = panel.cor, gap = 0)
  }

}

#' @importFrom graphics points polygon
CM_fit_esc <- function(report, d, year) {
  obs <- d$obsescape
  esc <- exp(sapply(report, getElement, "logpredesc")) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975))

  plot(NULL, NULL, typ = "n", pch = 16, xlim = range(year), ylim = c(0, 1.1) * range(esc, obs), xlab = "Year", ylab = "Escapement")
  polygon(c(rev(year), year), c(rev(esc[1, ]), esc[3, ]), col = alpha("grey", 0.5), border = NA)
  lines(year, esc[2, ], lwd = 1.5)

  points(year, obs, pch = 16)
  lines(year, obs, lty = 3)

  invisible()
}

CM_data <- function(obs, year, ylab) {

  if (!is.null(obs) && sum(obs)) {
    plot(year, obs, typ = "o", pch = 16, xlim = range(year),
         ylim = c(0, 1.1) * range(obs), xlab = "Year", ylab = ylab)
  }

  invisible()
}

CM_fit_CWTesc <- function(report, d, year1 = 1) {
  ebrood <- sapply(report, getElement, "ebrood", simplify = "array") %>%
    apply(1:2, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    mutate(Age = paste("Age", Var3)) %>%
    reshape2::dcast(Age + Year ~ Var1, value.var = "value")

  dat <- d$cwtesc

  nyears <- nrow(dat)
  nage <- ncol(dat)

  cwtesc <- dat %>%
    structure(dimnames = list(Year = year1 + seq(1, nyears) - 1, Age = seq(1, nage))) %>%
    reshape2::melt() %>%
    mutate(Age = paste("Age", Age)) %>%
    dplyr::filter(value > 0)

  g <- ggplot(ebrood, aes(Year)) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = alpha("grey", 0.5)) +
    geom_line(aes(y = `50%`), linewidth = 0.75) +
    facet_wrap(vars(Age), scales = "free_y") +
    geom_point(data = cwtesc, aes(Year, value), inherit.aes = FALSE) +
    geom_line(data = cwtesc, aes(Year, value), linetype = 3, inherit.aes = FALSE) +
    labs(x = "Brood year", y = "CWT Escapement")
  g
}

CM_fit_CWTcatch <- function(report, d, PT = TRUE, year1 = 1) {
  dat <- d[[ifelse(PT, "cwtcatPT", "cwtcatT")]]

  if (sum(dat)) {
    nyears <- nrow(dat)
    nage <- ncol(dat)
    cwtcat <- dat %>%
      structure(dimnames = list(Year = year1 + seq(1, nyears) - 1, Age = seq(1, nage))) %>%
      reshape2::melt() %>%
      mutate(Age = paste("Age", Age)) %>%
      dplyr::filter(value > 0)

    cbrood <- sapply(report, getElement, ifelse(PT, "cbroodPT", "cbroodT"), simplify = "array") %>%
      apply(1:2, quantile, probs = c(0.025, 0.5, 0.975)) %>%
      reshape2::melt() %>%
      mutate(Year = Var2 + year1 - 1) %>%
      mutate(Age = paste("Age", Var3)) %>%
      reshape2::dcast(Age + Year ~ Var1, value.var = "value")

    g <- ggplot(cbrood, aes(Year)) +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = alpha("grey", 0.5)) +
      geom_line(aes(y = `50%`), linewidth = 0.75) +
      facet_wrap(vars(Age), scales = "free_y") +
      geom_point(data = cwtcat, aes(Year, value), inherit.aes = FALSE) +
      geom_line(data = cwtcat, aes(Year, value), linetype = 3, inherit.aes = FALSE) +
      labs(x = "Brood year", y = ifelse(PT, "CWT preterminal catch", "CWT terminal catch"))
    g
  }
}

#' @importFrom dplyr rename
CM_maturity <- function(report, d, year1 = 1, brood = TRUE, annual = FALSE) {
  if (brood) {
    matt <- sapply(report, function(i) CY2BY(i[["matt"]]), simplify = 'array')
  } else {
    matt <- sapply(report, getElement, "matt", simplify = 'array')
  }
  matt_q <- apply(matt, 1:2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    rename(Age = Var3) %>%
    reshape2::dcast(Age + Var2 + Year ~ Var1)

  if (annual) {
    g <- matt_q %>%
      ggplot(aes(Age, `50%`)) +
      geom_line() +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
      facet_wrap(vars(Year)) +
      theme(
        panel.spacing = unit(0, "in"),
        strip.background = element_blank()
      ) +
      labs(x = "Age", y = "Proportion mature", title = ifelse(brood, "Brood year", "Return year"))
  } else {
    bmatt <- data.frame(Age = 1:d$Nages, value = d$bmatt) %>%
      dplyr::filter(Age > 1)

    g <- matt_q %>%
      dplyr::filter(Age > 1) %>%
      ggplot(aes(Year, `50%`, fill = factor(Age), colour = factor(Age))) +
      geom_line() +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
      geom_hline(data = bmatt, linetype = 2, aes(yintercept = value, colour = factor(Age))) +
      labs(x = ifelse(brood, "Brood year", "Return year"), y = "Proportion mature", colour = "Age", fill = "Age")
  }

  g
}

CM_vul <- function(report, type = c("vulPT", "vulT")) {
  type <- match.arg(type)
  vul <- sapply(report, getElement, type) # age x sim

  if (sum(vul)) {
    vul_q <- apply(vul, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
      reshape2::melt() %>%
      rename(Age = Var2) %>%
      reshape2::dcast(list("Age", "Var1"), value.var = "value")

    g <- vul_q %>%
      ggplot(aes(Age, `50%`)) +
      geom_line() +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
      labs(x = "Age", y = ifelse(type == "vulPT", "Preterminal vulnerability", "Terminal vulnerability"))
    g
  }
}

CM_SRR <- function(report, year1 = 1, gg = TRUE) {
  egg <- sapply(report, getElement, "egg") %>% apply(1, median)
  smolt <- sapply(report, function(x) x$N[-1, 1, 1]) %>% apply(1, median)
  epred <- seq(0, 1.1 * max(egg), length.out = 50)
  spred <- sapply(report, function(x) x$alpha * epred * exp(-x$beta * epred)) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975))

  year <- year1 + seq(1, length(egg))

  #plot(egg, smolt, xlim = c(0, 1.1) * range(egg), ylim = c(0, 1.1) * range(smolt),
  #     xlab = "Egg production", ylab = "Smolt production")
  #lines(epred, spred[2, ])
  #polygon(c(rev(epred), epred), c(rev(spred[1, ]), spred[3, ]), col = alpha("grey", 0.5), border = NA)
  #invisible()

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

      g <- ggplot(df, aes(Year, fill = Origin)) +
        geom_line(aes(y = `50%`, colour = Origin)) +
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

      g <- ggplot(df, aes(Year)) +
        geom_line(aes(y = `50%`)) +
        facet_wrap(vars(Age), scales = scales) +
        labs(x = xlab, y = ylab) +
        expand_limits(y = 0)
    }

    if (ci) g <- g + geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2)
    g

  }
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

  g <- ggplot(df, aes(Year, colour = Origin, fill = Origin)) +
    geom_line(aes(y = `50%`)) +
    labs(x = xlab, y = ylab) +
    expand_limits(y = 0)

  if (ci) g <- g + geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2)
  g
}

CM_M <- function(report, year1 = 1, ci = TRUE) {
  df <- sapply(report, getElement, "mo", simplify = "array") %>%
    apply(1:2, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    mutate(Age = paste("Age", Var3)) %>%
    reshape2::dcast(Year + Age ~ Var1, value.var = "value")

  g <- ggplot(df, aes(Year)) +
    geom_line(aes(y = `50%`)) +
    facet_wrap(vars(Age), scales = "free_y") +
    labs(x = "Year", y = "Natural mortality") +
    expand_limits(y = 0)

  if (ci) g <- g + geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = alpha("grey", 0.5))
  g
}


CM_Megg <- function(report, year1 = 1, ci = TRUE, surv = FALSE) {
  megg <- sapply(report, getElement, 'megg')

  if (surv) megg <- exp(-megg)

  df <- megg %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1"), value.var = "value")

  g <- ggplot(df, aes(Year)) +
    geom_line(aes(y = `50%`)) +
    labs(x = "Year", y = ifelse(surv, "Egg-smolt survival", "Egg-smolt mortality")) +
    expand_limits(y = 0)

  if (ci) g <- g + geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = alpha("grey", 0.5))
  g
}


CM_Njuv <- function(report, year1 = 1, ci = TRUE) {
  .CM_statevarage(report, year1, ci, "N", "Juvenile abundance") +
    theme(legend.position = "bottom")
}

CM_recr <- function(report, year1 = 1, ci = TRUE) {
  .CM_statevarage(report, year1, ci, "recr", "Recruitment") +
    theme(legend.position = "bottom")
}

CM_esc <- function(report, year1 = 1, ci = TRUE) {
  .CM_statevarage(report, year1, ci, "escyear", "Escapement") +
    theme(legend.position = "bottom")
}

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
    g <- ggplot(ts, aes(Year)) +
      geom_line(aes(y = `50%`)) +
      labs(x = xlab, y = ylab) +
      expand_limits(y = 0)

    if (ci) g <- g + geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = alpha("grey", 0.5))
    g
  }
}

CM_surv <- function(report, year1 = 1, ci = TRUE) {
  df <- exp(-sapply(report, getElement, "mo", simplify = "array")) %>%
    apply(1:2, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    mutate(Age = paste("Age", Var3)) %>%
    reshape2::dcast(Year + Age ~ Var1, value.var = "value")

  g <- ggplot(df, aes(Year)) +
    geom_line(aes(y = `50%`)) +
    facet_wrap(vars(Age), scales = "free_y") +
    labs(x = "Year", y = "Natural survival") +
    expand_limits(y = 0)

  if (ci) g <- g + geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = alpha("grey", 0.5))
  g
}

CM_wt <- function(stanfit, year1 = 1, ci = TRUE) {

  wt <- try(rstan::extract(stanfit, "wt")$wt, silent = TRUE)

  if (!is.character(wt)) {
    ts <- wt %>%
      apply(2, quantile, probs = c(0.025, 0.5, 0.975)) %>%
      reshape2::melt() %>%
      mutate(Year = Var2 + year1 - 1) %>%
      reshape2::dcast(list("Year", "Var1"), value.var = "value")

    g <- ggplot(ts, aes(Year, `50%`)) +
      geom_line() +
      geom_point() +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(x = "Year", y = "Egg mortality deviation")

    if (ci) g <- g + geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = alpha("grey", 0.5))
    g
  }

}

CM_wto <- function(stanfit, year1 = 1, ci = TRUE) {

  wto <- try(rstan::extract(stanfit, "wto")$wto, silent = TRUE)

  if (!is.character(wto)) {
    ts <- wto %>%
      apply(2, quantile, probs = c(0.025, 0.5, 0.975)) %>%
      reshape2::melt() %>%
      mutate(Year = Var2 + year1 - 1) %>%
      reshape2::dcast(list("Year", "Var1"), value.var = "value")

    g <- ggplot(ts, aes(Year, `50%`)) +
      geom_line() +
      geom_point() +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(x = "Year", y = "Age 1 mortality deviation")

    if (ci) g <- g + geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = alpha("grey", 0.5))
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

  Msurv <- array(NA_real_, dim(report$matt))
  if (brood) {
    matt <- CY2BY(report$matt)
    Msurv[, -ncol(Msurv)] <- CY2BY(report$mo)
  } else {
    matt <- report$matt
    Msurv[, -ncol(Msurv)] <- report$mo
  }
  Msurv[, ncol(Msurv)] <- 1

  AEQ <- array(NA_real_, dim(Msurv))
  AEQ[, ncol(AEQ)] <- 1

  nt <- nrow(AEQ)
  na <- ncol(AEQ)
  for (t in seq(nt, 2) - 1) {
    for (a in seq(na, 2) - 1) {
      AEQ[t, a] <- matt[t,a] * (1 - matt[t, a]) * Msurv[t, a+1] * AEQ[t, a+1]
    }
  }

  return(AEQ)
}


CM_BYER <- function(report, type = c("PT", "T", "all"), year1 = 1, ci = TRUE, at_age = TRUE) {

  type <- match.arg(type)

  if (at_age) {

    BYER_list <- lapply(report, function(i) {
      name <- ifelse(type == "PT", "survPT", "survT")
      CYER <- 1 - i[[name]]
      BYER <- CY2BY(CYER)
      return(list(BYER = BYER))
    })

    g <- .CM_statevarage(
      BYER_list,
      year1,
      ci,
      "BYER",
      ylab = switch(
        type,
        "PT" = "Preterminal exploitation rate",
        "T" = "Terminal exploitation rate"
      ),
      xlab = "Brood year",
      scales = "fixed"
    )

  } else {

    BYER_list <- lapply(report, function(i) {

      esc <- apply(i$escyear, 1:2, sum)
      esc_brood <- CY2BY(esc)

      morts_PT <- apply(i$cyearPT, 1:2, sum)
      morts_brood_PT <- CY2BY(morts_PT)

      morts_T <- apply(i$cyearT, 1:2, sum)
      morts_brood_T <- CY2BY(morts_T)

      AEQ_PT <- calc_AEQ(i)
      AEQ_T <- array(1, dim(esc))

      denom <- rowSums(morts_brood_PT * AEQ_PT + morts_brood_T * AEQ_T + esc_brood)

      if (type == "PT") {
        num <- rowSums(morts_brood_PT * AEQ_PT)
      } else if (type == "T") {
        num <- rowSums(morts_brood_T * AEQ_T)
      } else {
        num <- rowSums(morts_brood_PT * AEQ_PT + morts_brood_T * AEQ_T)
      }

      list(BYER = num/denom)
    })

    g <- .CM_ts(
      BYER_list, year1, ci,
      var = "BYER",
      xlab = "Brood year",
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


CM_CYER <- function(report, type = c("PT", "T", "all"), year1 = 1, ci = TRUE, at_age = TRUE) {

  type <- match.arg(type)

  if (at_age) {

    CYER_list <- lapply(report, function(i) {
      name <- ifelse(type == "PT", "survPT", "survT")
      return(list(CYER = 1 - i[[name]]))
    })

    g <- .CM_statevarage(
      CYER_list,
      year1,
      ci,
      "CYER",
      ylab = switch(
        type,
        "PT" = "Preterminal exploitation rate",
        "T" = "Terminal exploitation rate"
      ),
      xlab = "Calendar year",
      scales = "fixed"
    )

  } else {

    CYER_list <- lapply(report, function(i) {

      esc <- apply(i$escyear, 1:2, sum)

      morts_PT <- apply(i$cyearPT, 1:2, sum)
      morts_T <- apply(i$cyearT, 1:2, sum)

      AEQ_PT <- calc_AEQ(i, brood = FALSE)
      AEQ_T <- array(1, dim(esc))

      denom <- rowSums(morts_PT * AEQ_PT + morts_T * AEQ_T + esc)

      if (type == "PT") {
        num <- rowSums(morts_PT * AEQ_PT)
      } else if (type == "T") {
        num <- rowSums(morts_T * AEQ_T)
      } else {
        num <- rowSums(morts_PT * AEQ_PT + morts_T * AEQ_T)
      }

      list(CYER = num/denom)
    })

    g <- .CM_ts(
      CYER_list, year1, ci,
      var = "CYER",
      xlab = "Calendar year",
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


CM_covariate <- function(x, names, year1 = 1, b, ylab = "Covariate") {

  if (sum(x)) {

    if (!is.matrix(x)) x <- matrix(x, ncol = 1)
    if (missing(names)) names <- paste("Covariate", 1:ncol(x))

    if (missing(b)) { # Plot covariate

      g <- structure(x, dimnames = list(Year = 1:nrow(x) + year1 - 1, Covariate = names)) %>%
        reshape2::melt() %>%
        ggplot(aes(Year, value, colour = Covariate)) +
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
        ggplot(aes(Year, `50%`, colour = Covariate, fill = Covariate)) +
        geom_line() +
        geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
        labs(y = ylab, colour = NULL, fill = NULL) +
        expand_limits(y = 0) +
        theme(legend.position = "bottom")
    }

    return(g)

  }

}

reportCM <- function(stanfit, year, cov1_names, cov_names,
                     name, filename = "CM", dir = tempdir(), open_file = TRUE, render_args = list(), ...) {

  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    warning("Install ggrepel package to label years for stock-recruit figure.")
  }

  report <- get_report(stanfit)
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
