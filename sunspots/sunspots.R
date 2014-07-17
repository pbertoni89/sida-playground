# Autoregressive fit of Wolf's sunspot numbers
#
# George Udny Yule, On a Method of Investigating Periodicities in Disturbed Series,
# with Special Reference to Wolfer's Sunspot Numbers
# Philosophical Transactions of the Royal Society of London, 1927
#
# Run as follows:
# R --silent --no-save -f sunspots.R

require(datasets)            # This package contains many standard datasets,
                             # among them Wolf's numbers from the 18th to the 20th century
x <- sunspot.year[50:225]    # Wolf's numbers from 1749 to 1924, as in Yule's paper
x <- x - mean(x)             # "Detrend", that is center around 0

# Fit an AR(2) model. aic=FALSE and order.max=2 instruct the routine not to choose
# the model order automatically, but to let it be 2. 
# method="ols" means "use ordinary least squares"
model <- ar(x, aic=FALSE, order.max=2, method="ols")

print( round(as.vector(model$ar), digits=2) )  # Print the coefficients
resid_var <- var(model$resid, na.rm=TRUE)      # The first two residuals are not available;
                                               # na.rm=TRUE tells var() to ignore them
print( round(resid_var, digits=2) )            # Print the variance of residuals

# this script DOESN'T calculate OLS neither NLLS !
