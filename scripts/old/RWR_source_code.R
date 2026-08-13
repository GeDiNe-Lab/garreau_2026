function (W, ind.positives, gamma = 0.6, tmax = 1000, eps = 1e-10, 
    norm = TRUE, seed.weights = NULL) 
{
    if (norm) {
        M <- Prob.norm(W)
    }
    else {
        M <- W
    }
    n <- nrow(M)
    p0 <- p <- numeric(n)
    names(p) <- names(p0) <- rownames(W)
    rm(W)
    n.positives <- length(ind.positives)
    if (n.positives == 0) {
        stop("RWR: number of core positives is equal to 0!")
    }
    if (is.null(seed.weights)) {
        p0[ind.positives] <- 1/n.positives
    }
    else {
        if (length(seed.weights) != n.positives) {
            stop("seed.weights must have the same length as ind.positives")
        }
        if (any(seed.weights < 0)) {
            stop("seed.weights must be non-negative")
        }
        p0[ind.positives] <- seed.weights/sum(seed.weights)
    }
    p <- p0
    for (t in 1:tmax) {
        pold <- p
        p <- ((1 - gamma) * as.vector(pold %*% M)) + gamma * 
            p0
        if (norm1(p - pold) < eps) 
            (break)()
    }
    return(list(p = p, ind.positives = ind.positives, n.iter = t))
}
