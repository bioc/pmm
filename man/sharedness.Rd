\name{sharedness}
\alias{sharedness}
\encoding{utf-8}
\title{Sharedness Score}
\usage{
sharedness(fit, threshold = 0.2, na.action = "na.omit")
}
\arguments{
    \item{fit}{data frame returned by the \code{\link{pmm}} function.}
    \item{threshold}{threshold for the false discovery rate. Genes are
        counted as hits if they are below this threshold. Default is 0.2.}
    \item{na.action}{a function that indicates what happens if \code{fit}
        contains NAs. There are two options: "na.omit" (default) or "use". In
        the case of "na.omit" the sharedness score is applied to
        \code{na.omit(fit)} and in the other case the sharedness score is
        adapted for each gene to the number of conditions without NA.}
}
\value{The score returns a value between 0 and 1 for each gene. Score
    0 indicates that a gene is not shared among the condition and score 1
    that the gene is significant among all conditions.}
\description{
    The sharedness score computes the strength of sharedness of hit genes
    among the conditions.}
\details{
    The sharedness score is a combination of two quantities:\cr
    \deqn{sh_g = \frac{1}{2} \left( (1 - mean(fdr_{cg})) + \sum_{c}
        (fdr_{cg} < 1) \right)}{sh_g = 1/2 ( (1 - mean(fdr_cg)) + sum
        (fdr_cg < 1) )} \cr
    The first part defines the shift away from 1 and the second part
    describes how many pathogens support the shift (proportion
    of FDRs < 1).
}
\examples{
data(kinome)
fit <- pmm(kinome,"InfectionIndex","weight_library")
sh <- sharedness(fit, threshold = 0.2)
head(sh)

## NA-Handling
kinome$InfectionIndex[kinome$GeneID == 132158 & kinome$condition ==
"ADENO"] <- rep(NA,12)
fit <- pmm(kinome,"InfectionIndex","weight_library")
## Sharedness score for genes present in all conditions
sh <- sharedness(fit, threshold = 0.2, na.action = "na.omit")
head(sh)
## Sharedness score for all significant genes
sh <- sharedness(fit, threshold = 0.2, na.action = "use")
head(sh)
}
\author{
Anna Drewek <adrewek@stat.math.ethz.ch>
}