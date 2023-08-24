# Compute stats on LogicalRulesData using Joe's code

##First read in the arguments listed at the command line
args=(commandArgs(FALSE))
a <- args[length(args)-1]
b <- args[length(args)]
a <- sub("-","", a)
b <- sub("-","", b)

#rm(list = ls()) # Clear the global workspace
cat("\014")     # Clear the screen
graphics.off()

# Uncomment these if you haven't already installed them
# install.packages("sft")      # Install Houpt et al.'s SFT package

# Load all of the installed packages
library(here)
library(sft)      # sft
library(graphics)

# Relative filepath to the data directory
data_dir <- here("Rdata")

# Subject information
dataPrefix <- a # String at the beginning of data file
subjectNumber <- b

# Read data
datafilename <- paste(paste("R_analysis", dataPrefix, subjectNumber, sep = "_", collapse = NULL), ".dat", sep = "")
fullfilepath <- file.path(data_dir, datafilename)
data <- read.table(fullfilepath, sep = "\t", header = FALSE, row.names = NULL)
names(data) = c("Subject", "Condition", "RT", "Correct", "Channel1", "Channel2")

# Separate out target category items
target <- data[data$Channel1 > 0 & data$Channel2 > 0, ]

hh <- target$RT[target$Channel1 == 2 & target$Channel2 == 2]
hl <- target$RT[target$Channel1 == 2 & target$Channel2 == 1]
lh <- target$RT[target$Channel1 == 1 & target$Channel2 == 2]
ll <- target$RT[target$Channel1 == 1 & target$Channel2 == 1]

# Test SIC
sicresults <- sic(hh, hl, lh, ll)

# Set up time
mint = 0
maxt = 5000
dt = 10
t <- seq(mint, maxt, dt)

# Get SIC function
f_sic <- sicresults$SIC(t)

# Plot SIC function
plot(t, f_sic, type = "l", main = paste("SIC", "Subject", subjectNumber),
     xlab = "time", ylab = "SIC(t)")
lines(t, rep(0, 1, length(t)))

# Print stochastic dominance test
cat("\014")     # Clear the screen
print(sicresults$Dominance)
print("If stochastic dominance assumption is met:")
print("First 4 results should be significant, Last 4 should not")

# Print sic test results
sicresults$SICtest$positive

# Print sic test results
sicresults$SICtest$negative

# Print the MIC
sicresults$MIC