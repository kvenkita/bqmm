# Standalone check of bqmm_corr_standata Zcov reconstruction (no Stan needed).
suppressMessages({ library(lme4); library(Matrix) })
invisible(lapply(c("family.R","priors.R","formula.R","standata.R"),
                 function(f) source(file.path("R", f))))
ok <- function(c, m) if (isTRUE(c)) cat("PASS:", m, "\n") else stop("FAIL: ", m)

set.seed(1); n <- 80
d <- data.frame(y = rnorm(n), x = rnorm(n), g = factor(rep(1:8, length.out = n)))
parsed <- bqmm_parse_formula(y ~ x + (1 + x | g), data = d)
pr <- bqmm_default_priors(NULL, parsed$y, ncol(parsed$X))
sd <- bqmm_corr_standata(parsed, 0.5, pr)

ok(sd$M == 2, "M = 2")
ok(sd$L == nlevels(d$g), "L = nlevels(g)")
ok(all(dim(sd$Zcov) == c(n, 2)), "Zcov dims n x 2")
ok(all(sd$Zcov[,1] == 1), "Zcov col1 all ones (intercept)")
ok(isTRUE(all.equal(sd$Zcov[,2], d$x)), "Zcov col2 == x")
ok(length(sd$level_id) == n && all(sd$level_id == as.integer(d$g)), "level_id matches g")

# multi-term should error
err <- tryCatch({ p2 <- bqmm_parse_formula(y ~ x + (1|g) + (1|g:x), d); FALSE },
                error = function(e) TRUE)
d2 <- data.frame(y=rnorm(n), x=rnorm(n), g=factor(rep(1:8,length.out=n)), h=factor(rep(1:4,each=n/4)))
p2 <- bqmm_parse_formula(y ~ x + (1|g) + (1|h), d2)
err2 <- tryCatch({ bqmm_corr_standata(p2, 0.5, bqmm_default_priors(NULL,p2$y,2)); FALSE },
                 error = function(e) TRUE)
ok(err2, "multi-term unstructured errors")
cat("\ncorr standata checks passed.\n")
