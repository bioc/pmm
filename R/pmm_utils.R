## Functions for calculating locfdr.

local.fdr <- function(z){

    ## Fitting f (by poisson regression)
    intervals <- seq(min(z), max(z), length = 120)  ## fixed breaks of 120
    x <- (intervals[-1] + intervals[-length(intervals)])/2
    y <- hist(z, breaks = intervals, plot = FALSE)$counts
    xx <- ns(x, df = 7)
    f <- glm(y ~ xx, family = "poisson")$fit  ## Poisson Regression


    ## Fitting f0 and pi0 (by MLE)
    N <- length(z)
    scale <- (quantile(z,probs = 0.75) -
                  quantile(z,probs = 0.25))/(2*qnorm(0.75))
    c <- 4.3 * exp(-0.26*log(N,10))
    mlests <- mle.fct(z, Alim = c(median(z), c*scale))
    if(!is.na(mlests["delta0"])){
        mlests = mle.fct(z, Alim = c(mlests["delta0"], c*mlests["sigma0"]),
            delta=mlests["delta0"], sigma=mlests["sigma0"])
    }
    if(is.na(mlests["delta0"])) warning("MLE for f0 failed")
    ## MLE estimates
    pi0 = mlests["pi0"]
    f0 = dnorm(x, mlests["delta0"], mlests["sigma0"])
    f0 = (sum(f) * f0)/sum(f0)

    ## locfdr
    locfdr = pmin((pi0 * f0)/f, 1)

    ## correct locfdr at the borders
    if(sum(x >= mlests["delta0"] & locfdr == 1) > 0){
        xle <- max(x[x >= mlests["delta0"] & locfdr == 1])
    }else{
        xle = mlests["delta0"]
    }
    if(sum(x <= mlests["delta0"] & locfdr == 1) > 0){
        xri <- min(x[x <= mlests["delta0"] & locfdr == 1])
    }else{
        xri = mlests["delta0"]
    }
    if(sum(x <= xle & x >= xri) > 0)
        locfdr[x <= xle & x >= xri] <- 1

    locfdr[x >= mlests["delta0"] - mlests["sigma0"] &
               x <= mlests["delta0"] + mlests["sigma0"]] <- 1
    locfdr <- approx(x, locfdr, z, rule = 2, ties="ordered")$y

    ## Return local false discovery rate
    return(locfdr)
}



mle.fct <- function(z, Alim, delta = 0, sigma = 1, eps = 1e-05, max.it = 50){

    ## Define for data range A0 (region with no-hits)
    A0.left <- Alim[1] - Alim[2]
    A0.right <- Alim[1] + Alim[2]

    ## Total number of genes
    N <- length(z)

    ## Genes in A0
    z0 <- z[z >= (A0.left) & z<= (A0.right)]
    N0 <- length(z0)

    ## Sufficient Statistic
    moments <- c(mean(z0), mean(z0^2))

    ## Initial value
    theta <- N0/N

    ## Estimation MLE
    for(i in 1:max.it){

        est <- c(delta/sigma^2, -1/(2 * sigma^2))

        ## A0
        a.left <- (A0.left - delta)/sigma
        a.right <- (A0.right - delta)/sigma
        H0 <- pnorm(a.right) - pnorm(a.left)

        ## Hp = integral over [a,b] of f(z) = z^p phi(z) dz
        fal <- dnorm(a.left)
        far <- dnorm(a.right)
        H1 <- fal - far
        H2 <- H0 + a.left * fal - a.right * far
        H3 <- (2 + a.left^2) * fal - (2 + a.right^2) * far
        H4 <-
            3 * H0 + (3 * a.left + a.left^3) * fal -
                (3 * a.right + a.right^3) * far
        H <- c(H0, H1, H2, H3, H4)

        ## Expectation and Covariance Matrix of sufficient statistic
        HI <- matrix(rep(0, 25), 5)
        for(j in 0:4) HI[j + 1, 0:(j + 1)] <- choose(j, 0:j)
        HI <- sigma^(0:4) * (HI * (delta/sigma)^(pmax(row(HI) - col(HI), 0)))
        E <- as.vector(HI %*% H)/H0
        V <- matrix(c(E[3] - E[2]^2, E[4] - E[2] * E[3],
                      E[4] - E[2] * E[3], E[5] - E[3]^2), nrow = 2)

        ## MLE solution
        m <- solve(V, moments - c(E[2], E[3]))
        esti <- est + m/(1 + 1/i^2)
        if(esti[2] > 0)
            esti <- est + 0.1*m/(1 + 1/i^2)
        if(is.na(esti[2])){
            break
        }else if(esti[2]>=0){
            break
        }
        delta <- - esti[1]/(2 * esti[2])
        sigma <- 1/sqrt(-2 * esti[2])
        if(sum((esti - est)^2)^0.5 < eps) break  ## if convergence
    }
    if(is.na(esti[2])){
        mle <- rep(NA,3)
    }else if(esti[2] >=0){
        mle <- rep(NA, 3)
    }else{
        a.left <- (A0.left - delta)/sigma
        a.right <- (A0.right - delta)/sigma
        H0 <- pnorm(a.right) - pnorm(a.left)
        pi0 <- theta/H0
        mle <- c(delta, sigma, pi0)
    }
    names(mle) <- c("delta0", "sigma0", "pi0")
    return(mle)
}
