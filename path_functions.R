library(ggm)
library(igraph)
`ising` <- function(A, b) {
    v <- ncol(A)
    X <- data.matrix(expand.grid(rep(list(c(0, 1)), v)))
    p <- rep(0, 2^v)
    for (i in 1:(2^v)) {
        p[i] <- exp((1 / 2) * X[i, , drop = FALSE] %*% A %*% t(X[i, , drop = FALSE]) + X[i, ] %*% b)
    }
    p / sum(p)
}

`all_paths` <- function(G, a, b) {
    V <- names(V(G))
    pall <- all_simple_paths(G, a, b)
    pall <- lapply(pall, function(x) attr(x, "names"))
    for (pat in pall) {
        a <- pat[1]
        b <- pat[length(pat)]
        delta <- setdiff(pat, c(a, b))
        rest <- setdiff(V, c(a, delta, b))
        cat("Path: ", paste(pat, collapse = ","), "\n")
        cat("delta: ", paste(delta, collapse = ","), "\n")
        cat("rest: ", paste(rest, collapse = ","), "\n")
        cat("\n")
    }
    invisible(pall)
}

`path_prod_or` <- function(pall, A) {
    w1 <- rep(0, length(pall))
    w0 <- rep(0, length(pall))

    for (h in seq_along(pall)) {
        prodall <- prod(exp(A)[upper.tri(A)])
        pa <- pall[[h]]
        siz <- length(pa)
        ed <- c()
        for (i in 2:siz) {
            ed <- c(ed, exp(A[pa[i - 1], pa[i]]))
        }
        rst <- prodall / prod(ed)
        w0[h] <- rst
        w1[h] <- prod(ed)
        cat("Path: ", paste(pa, collapse = ","), "\n")
        cat(round(ed, 2), "omega: ", round(prod(ed), 2), "not-omega: ", round(rst, 2), "\n")
    }
    invisible(cbind(w0, w1))
}


`trivariate` <- function(pa, V, p) {
    v <- length(V)
    X <- data.frame(expand.grid(rep(list(c(0, 1)), v)))
    colnames(X) <- V
    a <- pa[1]
    b <- pa[length(pa)]
    delta <- setdiff(pa, c(a, b))
    rest <- setdiff(V, c(a, delta, b))
    Yd <- 0 + (apply(X[, delta, drop = FALSE] == 1, 1, all) & apply(X[, rest, drop = FALSE] == 0, 1, all))
    p_dab <- cbind(Yd, X[, c(a, b)], p = p)
    p_dab <- as.data.frame(xtabs(p ~ ., p_dab))
    p_dab
}

`bivariate` <- function(pa, V, p) {
    v <- length(V)
    X <- data.frame(expand.grid(rep(list(c(0, 1)), v)))
    colnames(X) <- V
    a <- pa[1]
    b <- pa[length(pa)]
    delta <- setdiff(pa, c(a, b))
    rest <- setdiff(V, c(a, delta, b))
    Zd <- 0 + cbind(apply(X[, c(a, b, delta)] == 1, 1, all), apply(X[, rest] == 1, 1, all))
    p_dd <- cbind(Zd, p)
    p_dd <- as.data.frame(xtabs(p ~ ., p_dd))
    p_dd
}

measures3 <- function(p_dab) {
    Psi <- p_dab$Freq[8] # prob. of P(Yd = 0, Xa = 1,Xb = 1)
    Theta <- p_dab$Freq[7] # prob. of P(Yd = 1, Xa = 1,Xb = 1)
    pc <- Psi / (Psi + Theta)
    c(pc = pc, Psi = Psi, Theta = Theta)
}

measures2 <- function(p_dd) {
    q <- p_dd$Freq
    OR <- q[1] * q[4] / (q[2] * q[3])
    Y <- (sqrt(OR) - 1) / (sqrt(OR) + 1)
    c(OR = OR, Y = Y)
}
# .............. all the calculations together .............
`path_calc` <- function(pat, V, p) {
    v <- length(V)
    if (v < 3) stop("The number of nodes is less then 3.")
    X <- data.frame(expand.grid(rep(list(c(0, 1)), v)))
    colnames(X) <- V
    a <- pat[1]
    b <- pat[length(pat)]
    delta <- setdiff(pat, c(a, b))
    rest <- setdiff(V, c(a, delta, b))
    Yd <- 0 + (apply(X[, delta] == 1, 1, all) & apply(X[, rest] == 0, 1, all))
    X3 <- cbind(Yd, X[, c(a, b)])
    # cells <- sapply(ggm::powerset(1:v, sort = FALSE), paste, collapse = "")
    # p <- c("p0", paste0("p", cells))
    p_dab <- cbind(X3, p)
    # o <- with(p_dab, order(p_dab[, 3], p_dab[, 2], p_dab[, 1]))
    dab <- as.data.frame(xtabs(p ~ ., p_dab))
    P <- dab$Freq
    # print((P[8] * P[2]) / (P[4] * P[6]))
    m_d <- glm(Freq ~ Yd * B * F, family = quasipoisson, data = dab)
    theta <- exp(coefficients(m_d))
    theta <- theta[c(1, 2, 3, 5, 4, 6, 7, 8)] # Inverse lex order
    num <- apply(X[, c(a, b, delta)] == 1, 1, all) & apply(X[, rest] == 0, 1, all)
    Psi <- cbind(X, p)[which(num == 1), v + 1] # probability that path is active
    Theta <- dab[7, 4] # probability of p(Yd = 0, Xa = 1,Xb = 1)
    pc <- Psi / (Psi + Theta)
    pc2 <- (P[c(2, 4, 6, 8)] / marg2(P, c(2, 3)))[4]
    browser()
    # Second decomposition
    Zd <- 0 + cbind(apply(X[, c(a, b, delta)] == 1, 1, all), apply(X[, rest] == 1, 1, all))
    p_dd <- cbind(Zd, p)
    dd <- as.data.frame(xtabs(p ~ ., p_dd))
    Q <- dd$Freq
    OR <- Q[1] * Q[4] / (Q[2] * Q[3])
    Y <- (sqrt(OR) - 1) / (sqrt(OR) + 1)
    list(Psi = Psi, Theta = Theta, pc = pc, OR = OR, Y = Y)
}

# ................. Utility functions.................

`kronpow` <- function(A, d) {
    H <- 1
    for (i in 1:d) {
        H <- H %x% A
    }
    H
}

`p2lam` <- function(p) {
    # Finds the loglinear parameters baseline coding for binary tables.
    d <- log2(length(p))
    M <- matrix(c(1, -1, 0, 1), 2, 2)
    Linv <- kronpow(M, d)
    lam <- Linv %*% log(p)
    lam
}

`marg2` <- function(p, m) {
    p <- as.matrix(p)
    d <- log2(nrow(p))
    a <- array(p, rep(2, d))
    matrix(margin.table(a, m), ncol = 1)
}

`logit_reg` <- function(p) {
    `logit` <- function(p) {
        # logit of a probability p.
        log(p) - log(1 - p)
    }
    p <- as.matrix(p)
    d <- log2(length(p))
    pc <- p[seq(2, 2^d, by = 2), ] / marg2(p, 2:d) # Conditional probabilities 1 | rest
    L <- matrix(c(1, 1, 0, 1), ncol = 2)
    L <- kronpow(L, log2(length(pc)))
    gamma <- solve(L) %*% logit(pc)
    gamma
}
`moebius` <- function(p) {
    p <- as.matrix(p)
    k <- nrow(p)
    d <- log2(k)
    m <- matrix(c(1, 1, 1, 0), 2, 2)
    X <- 1
    for (i in 1:d) {
        X <- X %x% m
    }
    X %*% p
}
# ..................................

# Two equivalente computations of (Yd, Xa, Xb)
ydelta <- function(pa, V) {
    v <- length(V)
    if (v < 3) stop("The number of nodes is less then 3.")
    X <- data.frame(expand.grid(rep(list(c(0, 1)), v)))
    colnames(X) <- V
    a <- pa[1]
    b <- pa[length(pa)]
    delta <- setdiff(pa, c(a, b))
    rest <- setdiff(V, c(a, delta, b))

    Ydelta <- matrix(0, nrow = 2^v, ncol = 2)
    for (i in 1:2^v) {
        Ydelta[i, ] <- c(
            all(X[i, delta] == 1),
            all(X[i, rest] == 0)
        )
    }

    Yd <- matrix(0, nrow = 2^v, ncol = 1)
    for (i in 1:2^v) {
        if (Ydelta[i, 1] == 1 & Ydelta[i, 2] == 1) {
            Yd[i] <- 1
        } else {
            Yd[i] <- 0
        }
    }
    data.frame(Yd, Xa = X[, a], Xb = X[, b])
}

ydelta2 <- function(pa, V, p) {
    v <- length(V)
    X <- data.frame(expand.grid(rep(list(c(0, 1)), v)))
    colnames(X) <- V
    a <- pa[1]
    b <- pa[length(pa)]
    delta <- setdiff(pa, c(a, b))
    rest <- setdiff(V, c(a, delta, b))
    Yd <- 0 + (apply(X[, delta, drop = FALSE] == 1, 1, all) & apply(X[, rest, drop = FALSE] == 0, 1, all))
    cbind(Yd, X[, c(a, b)])
}
