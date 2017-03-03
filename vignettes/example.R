library(normalregMix)
library(ggplot2)

# produce multiple plots on a single plot
# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

transition.probs <- matrix(c(0.8,0.15,0.2,0.85), ncol = 2)
beta <- matrix(c(0.7,0.3), ncol = 2)
mu = c(-1,1)
beta <- matrix(c(0.3,0.2,-0.2), ncol = 1)
mu = c(-1,1)
sigma = c(0.8)
theta2 <- list(beta = beta, mu = mu, sigma = sigma,
               transition.probs = transition.probs,
               initial.dist = c(1,rep(0,15)))
theta <- theta2

M <- ncol(theta$transition.probs)
s <- nrow(as.matrix(theta$beta))
is.beta.switching <- (ncol(as.matrix(theta$beta)) > 1)
is.sigma.switching <- (length(theta$sigma) > 1)
is.MSM <- TRUE
n = 300
# generates data
sample <- GenerateSample(theta, n = n,
                         is.MSM = is.MSM)
y <- sample$y
model <- EstimateMSAR(y = y, 
                      M = M, s = s,
                      is.beta.switching = is.beta.switching,
                      is.sigma.switching = is.sigma.switching,
                      is.MSM = is.MSM)
plot.actual <- DiagPlot(model, y = y) + ggtitle("MSM-AH(1) Model (Actual)") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"))
plot.estimated <- DiagPlot(sample$msar.model, y = y) + ggtitle("MSM-AH(1) Model (Estimated)") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"))
multiplot(plot.actual, plot.estimated, cols = 1)
(fisher.estimated <- model$fisher.estimated)
(theta.estimated <- model$theta)

new.y <- sample$y.sample + rnorm(n)
EstimateFisherInformation(theta = theta.estimated, y = new.y, y.lagged = sample$y.lagged, 
                          z.dependent = as.matrix(rep(0,n)), 
                          z.independent = as.matrix(rep(0,n)),
                          z.dependent.lagged = as.matrix(rep(0,n)), 
                          z.independent.lagged = as.matrix(rep(0,n)),
                          is.MSM = TRUE)
