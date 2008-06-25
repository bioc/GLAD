--- glad.Rd	2008-05-07 15:08:11.000000000 -0400
+++ glad1.Rd	2008-05-07 15:32:50.000000000 -0400
@@ -5,7 +5,7 @@
 \title{Analysis of array CGH data}
 \description{
   This function allows the detection of breakpoints in genomic profiles
-  obtained by array CGH technology and affects a status (gain, normal
+  obtained by array CGH technology and attributes a status (gain, normal
   or lost) to each clone.
 }
 \usage{
@@ -24,7 +24,7 @@
 \arguments{
   \item{profileCGH}{Object of class \code{\link{profileCGH}}}
 
-  \item{mediancenter}{If \code{TRUE}, LogRatio are center on their median.}
+  \item{mediancenter}{If \code{TRUE}, LogRatio are centered on their median.}
   
   \item{smoothfunc}{Type of algorithm used to smooth \code{LogRatio} by a
     piecewise constant function. Choose either \code{lawsglad}, \code{\link[aws]{aws}} or
@@ -50,7 +50,7 @@
     stochastic penalty (see \code{\link[aws]{aws}} or
     \code{\link[aws]{laws}})}
   
-  \item{base}{If \code{TRUE}, the position of clone is the physical position onto
+  \item{base}{If \code{TRUE}, the position of clone is the physical position on
     the chromosome, otherwise the rank position is used.}
   
   \item{sigma}{Value to be passed to either argument \code{sigma2}   
@@ -58,7 +58,6 @@
     \code{\link[aws]{laws}}. If \code{NULL}, sigma is calculated from
     the data.}
   
-  
   \item{lambdabreak}{Penalty term (\eqn{\lambda'}) used during the 
     \bold{Optimization of the number of breakpoints} step.}
   
@@ -98,33 +97,31 @@
 \details{
 
   The function \code{glad} implements the methodology which
-  is described in the article : Analysis of array CGH data: from signal
+  is described in the article: Analysis of array CGH data: from signal
   ratio to gain and loss of DNA regions (Hupé et al., Bioinformatics 2004 20(18):3413-3422).
 
-  The principle of the GLAD algorithm:
+  The principles of the GLAD algorithm:
   First, the detection of breakpoints is based on the estimation of a
   piecewise constant function with the Adaptive Weights Smoothing (AWS)
   procedure (Polzehl and Spokoiny, 2002).
-  Thus, a procedure based on penalyzed maximum likelihood optimizes the number of
-  breakpoints allows the undesirable breakpoints to be removed.
+  Then, a procedure based on penalyzed maximum likelihood optimizes the number of
+  breakpoints and removes the undesirable breakpoints.
   Finally, based on the regions previously identified, a two-step
-  unsupervised classification (\bold{MSHR 
-    clustering by chromosome} and the \bold{HCSR
+  unsupervised classification (\bold{MSHR clustering by chromosome} and the \bold{HCSR 
     clustering throughout the genome}) with model selection criteria
   allows a status to be assigned for each region (gain, loss or
   normal).
 
   Main parameters to be tuned:
   \tabular{ll}{
-    \code{qlambda} \tab if you want the smoothing to fit some very local effect, choose a smaller \code{qlambda}. \cr
-    \code{bandwidth} \tab choose a bandwidth not to small otherwise you
-    will have a lot of little discontinuities. \cr
-    \code{lambdabreak} \tab More the parameter is high more the number
-    of undesirable breakpoints is high.\cr
-    \code{lambdacluster} \tab More the parameter is high more the regions
-    within a chromosome are supposed to belong to the same cluster. \cr
-    \code{lambdaclusterGen} \tab More the parameter is high more the
-    regions over the whole genome are supposed to belong to the same cluster.
+    \code{qlambda} \tab If you want the smoothing to fit some very local effects, choose a smaller \code{qlambda}. \cr
+    \code{bandwidth} \tab bandwidth shouldn't be too small otherwise several little discontinuities may appear. \cr
+    \code{lambdabreak} \tab The higher this parameter is, the higher the number
+    of undesirable breakpoints is.\cr
+    \code{lambdacluster} \tab The higher this parameter is, the higher is the number of regions
+    within a chromosome and supposed to belong to the same cluster. \cr
+    \code{lambdaclusterGen} \tab The higher this parameter is, the higher is the number of
+    regions over the whole genome that are supposed to belong to the same cluster.
     
     }
 
@@ -148,24 +145,24 @@
     \item{\bold{Breakpoints}}{The last position of a region with identical amount
       of DNA is flagged by 1 otherwise it is 0. Note that during the
       "Optimization of the number of breakpoints" step, removed
-      breakpoints are flagged by -1.}
+      breakpoints are flagged -1.}
     
 
     \item{\bold{Region}}{Each position between two breakpoints are labelled the
       same way with an integer value starting from one. The label is
-      incremented by one when a new breakpoints occurs or when moving to
+      incremented by one when a new breakpoint occurs or when moving to
       the next chromosome. The variable \code{region} is what we call MSHR.}
     
-    \item{\bold{Level}}{Each position with equal smoothing value are labelled the
+    \item{\bold{Level}}{Each position with equal smoothing value is labelled the
       same way with an integer value starting from one. The label is
-      incremented by one when a new level occurs or when moving to the
+      incremented by one when a new level is found or when moving to the
       next chromosome.}
     
-    \item{\bold{OutliersAws}}{Each AWS outliers are flagged by -1 or 1 
+    \item{\bold{OutliersAws}}{Each AWS outlier is flagged -1 or 1 
       otherwise  it is 0.}
     
 
-    \item{\bold{OutliersMad}}{Each MAD outliers are flagged by -1 (if it is
+    \item{\bold{OutliersMad}}{Each MAD outliers is flagged -1 (if it is
       in the \eqn{\alpha/2} lower tail of the distribution) or 1 (if it is
       in the \eqn{\alpha/2} upper tail of the distribution)
       otherwise  it is 0.}
