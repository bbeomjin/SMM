\name{nottingham}
\alias{nottingham}
\title{
nottingham database
}
\description{
50 males and 50 females in neutral, frontal pose
Resolution: 358x463 to 468x536 monochrome as 200x200 monochrome
}
\usage{
nottingham
}
%- maybe also 'usage' for other objects documented here.
\format{
  A data frame with 100 rows and 40001 variables.
}

\details{
This is a collection of images useful for conducting experiments in psychology, primarily faces, though other submissions are welcome. They are free for research use.

Please note, there are no pictures of anyone famous here, at least, not outside the world of academia!

To cite PICS in a paper, please quote this URL, \url{pics.stir.ac.uk}
\itemize{
  \item gender. 1 is female, 0 is male.
  \item pixel. This is 200 x 200 images containing 40000 pixels(200 x 200)
}
}
\references{
Psychology, School of Natural Sciences
University of Stirling
\url{pics.stir.ac.uk}
}

\examples{
data(nottingham)

# sample of gender == 1
face_ind = which(nottingham$gender == 1)
plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", axes = FALSE)
rasterImage(matrix(as.numeric(nottingham[face_ind[1], -1]), ncol = 200), 0, 0, 1, 1)

# sample of gender == 0
face_ind = which(nottingham$gender == 0)
plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", axes = FALSE)
rasterImage(matrix(as.numeric(nottingham[face_ind[1], -1]), ncol = 200), 0, 0, 1, 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
