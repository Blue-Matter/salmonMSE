
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
  lines(year, esc[2, ])
  points(year, obs, pch = 16)

  invisible()
}

CM_fit_CWTesc <- function(report, d, year1) {
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
    geom_line(aes(y = `50%`)) +
    facet_wrap(vars(Age), scales = "free_y") +
    geom_point(data = cwtesc, aes(Year, value), inherit.aes = FALSE) +
    labs(x = "Brood year", y = "CWT Escapement")
  g
}

CM_fit_CWTcatch <- function(report, d, PT = TRUE, year1) {
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
      geom_line(aes(y = `50%`)) +
      facet_wrap(vars(Age), scales = "free_y") +
      geom_point(data = cwtcat, aes(Year, value), inherit.aes = FALSE) +
      labs(x = "Brood year", y = ifelse(PT, "CWT preterminal catch", "CWT terminal catch"))
    g
  }
}

#' @importFrom dplyr rename
CM_maturity <- function(report, d, year1) {
  matt <- sapply(report, getElement, "matt", simplify = 'array')
  matt_q <- apply(matt, 1:2, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    rename(Age = Var3) %>%
    dplyr::filter(Age > 1)

  bmatt <- data.frame(Age = 1:6, value = d$bmatt) %>%
    dplyr::filter(Age > 1)

  g <- matt_q %>%
    reshape2::dcast(Age + Var2 + Year ~ Var1) %>%
    ggplot(aes(Year, `50%`, fill = factor(Age), colour = factor(Age))) +
    geom_line() +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
    geom_hline(data = bmatt, linetype = 2, aes(yintercept = value, colour = factor(Age))) +
    labs(x = "Year", y = "Proportion mature", colour = "Age", fill = "Age")
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

CM_SRR <- function(report) {
  egg <- sapply(report, getElement, "egg") %>% apply(1, median)
  smolt <- sapply(report, function(x) x$N[-1, 1, 1]) %>% apply(1, median)
  epred <- seq(0, 1.1 * max(egg), length.out = 50)
  spred <- sapply(report, function(x) x$alpha * epred * exp(-x$beta * epred)) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975))

  plot(egg, smolt, xlim = c(0, 1.1) * range(egg), ylim = c(0, 1.1) * range(smolt),
       xlab = "Egg production", ylab = "Smolt production")
  lines(epred, spred[2, ])
  polygon(c(rev(epred), epred), c(rev(spred[1, ]), spred[3, ]), col = alpha("grey", 0.5), border = NA)

  invisible()
}


.CM_statevarage <- function(report, year1, ci = TRUE, var, ylab) {
  df <- sapply(report, getElement, var, simplify = "array") %>%
    apply(1:3, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1, Origin = ifelse(Var4 == 1, "Natural", "Hatchery")) %>%
    mutate(Age = paste("Age", Var3)) %>%
    reshape2::dcast(Year + Age + Origin ~ Var1, value.var = "value")

  g <- ggplot(df, aes(Year, fill = Origin)) +
    geom_line(aes(y = `50%`, colour = Origin)) +
    facet_wrap(vars(Age), scales = "free_y") +
    labs(x = "Year", y = ylab) +
    expand_limits(y = 0)

  if (ci) g <- g + geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2)
  g
}

CM_M <- function(report, year1, ci = TRUE) {
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

CM_Njuv <- function(report, year1, ci = TRUE) {
  .CM_statevarage(report, year1, ci, "N", "Juvenile abundance") +
    theme(legend.position = "bottom")
}

CM_recr <- function(report, year1, ci = TRUE) {
  .CM_statevarage(report, year1, ci, "recr", "Recruitment") +
    theme(legend.position = "bottom")
}

CM_esc <- function(report, year1, ci = TRUE) {
  .CM_statevarage(report, year1, ci, "escyear", "Escapement") +
    theme(legend.position = "bottom")
}

CM_F <- function(report, PT = TRUE, year1, ci = TRUE) {
  .CM_ts(
    report, year1, ci,
    var = ifelse(PT, "FPT", "FT"),
    ylab = ifelse(PT, "Preterminal fishing mortality", "Terminal fishing mortality")
  )
}

.CM_ts <- function(report, year1, ci = TRUE, var, ylab) {
  ts <- sapply(report, getElement, var) %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1"), value.var = "value")

  if (sum(ts$`50%`)) {
    g <- ggplot(ts, aes(Year)) +
      geom_line(aes(y = `50%`)) +
      labs(x = "Year", y = ylab) +
      expand_limits(y = 0)

    if (ci) g <- g + geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = alpha("grey", 0.5))
    g
  }
}

CM_surv <- function(report, year1, ci = TRUE) {
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

CM_wt <- function(stanfit, year1, ci = TRUE) {

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

CM_wto <- function(stanfit, year1, ci = TRUE) {

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

reportCM <- function(stanfit, year,
                     name, filename = "CM", dir = tempdir(), open_file = TRUE, render_args = list(), ...) {

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
