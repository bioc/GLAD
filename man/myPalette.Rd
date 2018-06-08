
\name{myPalette}
\alias{myPalette}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Microarray color palette}
\description{
This function returns a vector of color names corresponding to a range of colors specified in the arguments.
}
\usage{
myPalette(low = "white", high = c("green", "red"), mid=NULL, k =50)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{low}{Color for the lower end of the color palette, specified using any of the three kinds of R colors, i.e., either a color name (an element of \code{colors}), a hexadecimal string of the form \code{"#rrggbb"}, or an integer \code{i} meaning \code{palette()[i]}.}
  \item{high}{Color for the upper end of the color palette, specified
  using any of the three kinds of R colors, i.e., either a color name
  (an element of \code{colors}), a hexadecimal string of the form
  \code{"#rrggbb"}, or an integer \code{i} meaning \code{palette()[i]}.}
 \item{mid}{Color for the middle portion of the color palette, specified using any of the three kinds of R colors, i.e., either a color name (an element of \code{colors}), a hexadecimal string of the form \code{"#rrggbb"}, or an integer \code{i} meaning \code{palette()[i]}.}
  \item{k}{Number of colors in the palette.}
}

\value{
  A "character" vector of color names. This can be used to create a user-defined color palette for subsequent graphics by \code{palette}, in a \code{col=} specification in graphics functions, or in \code{\link{par}}.
}



\author{Sandrine Dudoit, Yee Hwa (Jean) Yang.}

\keyword{color}

\seealso{\code{palette}, \code{rgb},
\code{colors}, \code{col2rgb}, \code{\link{image}}, \code{\link{ColorBar}}, \code{\link{arrayPlot}}.}

\examples{

par(mfrow=c(1,4))
pal <- myPalette(low="red", high="green")
ColorBar(seq(-2,2, 0.2), col=pal, horizontal=FALSE, k=21)
pal <- myPalette(low="red", high="green", mid="yellow")
ColorBar(seq(-2,2, 0.2), col=pal, horizontal=FALSE, k=21)
pal <- myPalette()
ColorBar(seq(-2,2, 0.2), col=pal, horizontal=FALSE, k=21)
pal <- myPalette(low="purple", high="purple",mid="white")
ColorBar(seq(-2,2, 0.2), col=pal, horizontal=FALSE, k=21)

}


 
