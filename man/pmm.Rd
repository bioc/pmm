\name{pmm}
\alias{pmm}
\encoding{utf-8}
\title{Fitting the PMM}
\description{Fits the parallel mixed model.}
\usage{pmm(df.data, response, weight = "None", ignore = 3, simplify =
TRUE, gene.col = "GeneID", condition.col = "condition")}
\arguments{
    \item{df.data}{
        a data frame containing the variables for the model. Each row should
        correspond to one independent siRNA experiment. The data frame needs
        to have at least the following variables: GeneID, condition and
        a column with the measurements/readouts of the screens.}
    \item{response}{
        name of the column that contains the measurements/readouts of
        the screens.
    }
    \item{weight}{
        an optional vector of weights to be used in the fitting
        process of the linear mixed model. It should be a numeric
        vector. Default is a fit without weights.
    }
    \item{ignore}{number of minimal required sirna replicates for each
        gene. If a gene has less siRNA replicates it is ignored during the
        fitting process. Default is 3.
    }
    \item{simplify}{logical value that indicates whether the output of
        pmm should be simplified.
    }
    \item{gene.col}{name of the column that give a gene
        identifier. Default is "GeneID".
    }
    \item{condition.col}{name of the column that indicates the
        condition that was used for each measurement. Default is "condition".
    }
}
\details{
    The Parallel Mixed Model (PMM) is composed of a linear mixed model and
    an assessment of the local False Discovery Rate. The linear mixed
    model consists of a fixed effect for condition and of two random
    effects for gene g and for gene g within a condition c. We fit a
    linear mixed model by using \code{lmer} function from \code{lme4}
    R-package. To distinguish hit genes, PMM provides also an estimate of
    the local False Discovery Rate (FDR).
    \code{pmm} will only use the data of genes that have at least a
    certain number of siRNA replicates per condition. The number of
    ignored genes can be passed to \code{pmm} by the argument
    \code{ignore}. We recommend using at least 3 siRNA replicates per
    gene and condition in order to obtain a reliable fit.
}
\value{The simplified output of \code{pmm} is a matrix that contains the
    c_cg effects for each condition c and gene g, as well as an estimate
    for the local false discovery rate. A positive estimated c_cg effect
    means that the response was enhanced when the corresponding gene is
    knocked down. A negative effect means that the response was reduced.\cr
    The non-simplified output of \code{pmm} is a list of three
    components. The first component contains the simpilified output, i.e
    the matrix with the c_cg effects and fdr values, the second component
    contains the fit of the linear mixed model and the third component
    contains the a_g and b_cg values.
}
\examples{
 data(kinome)

 ## Fitting the parallel mixed model with weights
 fit1 <- pmm(kinome,"InfectionIndex","weight_library")
 head(fit1)

 ## Fitting the parallel mixed model without weights
 fit2 <- pmm(kinome,"InfectionIndex","None")
 head(fit2)

 ## Accessing the fit of the linear mixed model
 fit3 <- pmm(kinome,"InfectionIndex","weight_library",simplify=FALSE)
 identical(fit1,fit3[[1]])
 summary(fit3[[2]])

 ## NA-Handling
 kinome$InfectionIndex[kinome$GeneID == 10000 & kinome$condition ==
 "ADENO"] <- rep(NA,12)
 fit4 <- pmm(kinome,"InfectionIndex","weight_library",3)
 head(fit4)
}
\author{
Anna Drewek <adrewek@stat.math.ethz.ch>
}