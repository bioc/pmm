\name{pmm-package}
\alias{pmm-package}
\docType{package}
\encoding{utf-8}
\title{The PMM-Package}
\description{
    This package contains R functions for fitting the Parallel Mixed Model and
    analyzing its results.}
\details{
      The Parallel Mixed Model (PMM) approach is suitable for hit selection and
cross-comparison of RNAi screens generated in experiments that are performed in
parallel under several conditions. As an example, we could think of the
measurements or readouts from cells under RNAi knock-down, which are infected
with several pathogens or which are grown from different cell lines. PMM
simultaneously takes into account all the knock-down effects in order to gain
statistical power for the hit detection. As a special feature, PMM allows
incorporating RNAi weights that can be assigned according to the additional
information on the used RNAis or the screening quality.\cr
    The following functions are contained in this R-package:
    \tabular{cl}{
        \code{\link{pmm}} \tab fits the PMM  \cr
        \code{\link{hitheatmap}} \tab visualizes the results of PMM \cr
        \code{\link{sharedness}} \tab computes the sharedness
    }
}
\author{Anna Drewek <adrewek@stat.math.ethz.ch>}
\references{Rämö, P., Drewek, A., Arrieumerlou, C., Beerenwinkel, N., Ben-Tekaya
    H., Cardel, B., Casanova, A., Conde-Alvarez. R., Cossart, P., Csucs, G.,
    Eicher, S., Emmenlauer, M. Greber, U., Hardt, W.-D., Helenius, A., Kasper,
    C., Kaufmann, A., Kreibich, S., Kuebacher, A., Kunszt, P., Low, S.H.,
    Mercer, J., Mudrak, D., Muntwiler, S., Pelkmans, L., Pizarro-Cerda, J.,
    Podvinec, M., Pujadas, E., Rinn, B., Rouilly, V., Schmich F., Siebourg, J.,
    Snijder, B., Stebler, M., Studer, G., Szczurek, E., Truttmann, M., von
    Mering, C., Vonderheit, A., Yakimovich, A., Buehlmann, P. and Dehio, C.,
    \emph{Simultaneous analysis of large-scale RNAi screens for pathogen entry},
    BMC Genomics 15(1162): p.1471-2164. (2014)
}
