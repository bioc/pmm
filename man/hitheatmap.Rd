\name{hitheatmap}
\alias{hitheatmap}
\title{Visualization of the PMM results}
\usage{hitheatmap(fit, threshold = 0.2, sharedness.score = FALSE,
                  main = "", na.action = "use", \dots)}
\arguments{
    \item{fit}{data frame returned by the \code{\link{pmm}} function.}
    \item{threshold}{threshold for the false discovery rate. Genes are
        counted as hits if they are  below this threshold. Default is 0.2.}
    \item{sharedness.score}{logical value that indicates whether
        the sharedness score among the conditions should be additionally
        plotted. Default is FALSE.}
    \item{main}{the title at the top of the plot.}
    \item{na.action}{a function that indicates what happens if \code{fit}
        contains NAs. There are two options: "na.omit" or "use" (default). In
        the case of "na.omit" the hitheatmap is plotted for
        \code{na.omit(fit)} and in the other case the hitheatmap plots all
        data in fit.}
    \item{\dots}{further arguments passed to plot and par function.}
}
\value{A heat map}
\description{
    This function visualizes the results of PMM.
}
\details{
    The heat map represents the effects c_cg estimated by PMM. Red color
    indicates a positive c_cg coefficient, blue color a negative c_cg
    coefficient. The darker the color, the stronger is the c_cg effect. The
    heat map contains only the genes with false discovery rate below the given
    threshold for at least one condition. The yellow star indicates the hit
    genes in each condition. If sharedness.score = TRUE, an additional
    row is plotted. It represents the strength of sharedness for a gene
    among the conditions. The darker the color, the stronger is the
    sharedness effect. If na.omit = "use" then NAs are plotted in white
    color and marked by "NA".
}
\examples{
data(kinome)
fit1 <- pmm(kinome,"InfectionIndex","weight_library")

hitheatmap(fit1, threshold = 0.4)
hitheatmap(fit1, threshold = 0.2, main = "Results PMM")
hitheatmap(fit1, sharedness.score = TRUE)

## NA-Handling
kinome$InfectionIndex[kinome$GeneID == 3611 & kinome$condition ==
"ADENO"] <- rep(NA,12)
fit2 <- pmm(kinome,"InfectionIndex","weight_library")
hitheatmap(fit2, main = "Results PMM with NA")

## Using par options
hitheatmap(fit1, sharedness.score = TRUE, cex.main = 2,
                main = "My modified plot", col.main = "white",
                col.axis = "white", cex.axis = 0.8, bg = "black",
                mar = c(7,6,4,6))
}
\author{
Anna Drewek <adrewek@stat.math.ethz.ch>
}