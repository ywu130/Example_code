curve_fitting
================

## plot the mean+-sem

``` r
# Load the required libraries
library(ggplot2)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ lubridate 1.9.3     ✔ tibble    3.2.1
    ## ✔ purrr     1.0.2     ✔ tidyr     1.3.0
    ## ✔ readr     2.1.4

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
# Calculate the mean and standard error for each group at each time point
summary_data <- data %>%
  reframe(
    c1_mean = rowMeans(.[, 2:7], na.rm = TRUE),
    c1_sem = apply(.[, 2:7], 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))),
    c2_mean = rowMeans(.[, 8:13], na.rm = TRUE),
    c2_sem = apply(.[, 8:13], 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))),
    n_mean = rowMeans(.[, 14:19], na.rm = TRUE),
    n_sem = apply(.[, 14:19], 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  )



# Convert row names to a column
summary_data$Days <- rownames(summary_data)
rownames(summary_data) <- NULL

# Reshape the data into long format
summary_data_long <- tidyr::pivot_longer(summary_data, 
                                         cols = -Days, 
                                         names_to = c("Group",".value"), 
                                         names_pattern = "(.*)_(.*)")

# Plot the dot plot
ggplot(summary_data_long, aes(x = Days, y = mean, color = Group)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - sem, ymax =  mean + sem), width = 0.2) +
  labs(x = "Days", y = "Cell number(10^5)", color = "Group") +
  theme_classic()
```

![](curve_fitting_files/figure-gfm/plot%20mean+-sem-1.png)<!-- -->

## fitting the curves

``` r
library(reshape2)
```

    ## 
    ## Attaching package: 'reshape2'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     smiths

``` r
# Convert the 'days' column to a factor
data$days <- factor(data$days)

##need to add column name before reshape
#change col names
column_names <- colnames(data)  # Get the current column names as a vector

# Change column names for each group
column_names[2:7] <- "c1"
column_names[8:13] <- "c2"
column_names[14:19] <- "n"

colnames(data) <- column_names  # Assign the updated column names back to the dataframe



# Reshape the data from wide to long format
df_long <- data %>%
  pivot_longer(cols = -days, names_to = "group", values_to = "values") %>%
  drop_na()

# Convert the 'group' column to a factor
df_long$group <- factor(df_long$group, levels = c("c1", "c2", "n"))

# Plot the dot plot
ggplot(df_long, aes(x = days, y = values, color = group)) +
  geom_point() +
  scale_color_manual(values = c("c1" = "red", "c2" = "blue", "n" = "green")) +
  labs(x = "Days", y = "Cell number", color = "Group") +
  theme_minimal()
```

![](curve_fitting_files/figure-gfm/fitting%20the%20curve-1.png)<!-- -->

``` r
# Assuming your dataframe is called 'df_long'

# Fit logistic regression models for each group
groups <- unique(df_long$group)

#change days into numeric (for model fitting)
df_long$days  <- as.numeric(df_long$days)


for (value in groups) {
  
  subset_data <- subset(df_long, group == value)
  
  # Fit logistic growth model
  model <- nls(values ~ SSlogis(days, Asym, xmid, scal), data = subset_data)
  
  # Print model summary
  cat("Group:", value, "\n")
  print(summary(model))
  cat("\n")
}
```

    ## Group: c1 
    ## 
    ## Formula: values ~ SSlogis(days, Asym, xmid, scal)
    ## 
    ## Parameters:
    ##      Estimate Std. Error t value Pr(>|t|)    
    ## Asym 37.16348    1.32818  27.981  < 2e-16 ***
    ## xmid  3.50103    0.11827  29.601  < 2e-16 ***
    ## scal  0.70171    0.09138   7.679 1.58e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.864 on 21 degrees of freedom
    ## 
    ## Number of iterations to convergence: 6 
    ## Achieved convergence tolerance: 5.098e-06
    ## 
    ## 
    ## Group: c2 
    ## 
    ## Formula: values ~ SSlogis(days, Asym, xmid, scal)
    ## 
    ## Parameters:
    ##      Estimate Std. Error t value Pr(>|t|)    
    ## Asym  33.3783     1.5527  21.496  < 2e-16 ***
    ## xmid   3.2011     0.1697  18.864 6.72e-16 ***
    ## scal   0.7160     0.1361   5.261 2.15e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 4.077 on 24 degrees of freedom
    ## 
    ## Number of iterations to convergence: 6 
    ## Achieved convergence tolerance: 3.219e-06
    ## 
    ## 
    ## Group: n 
    ## 
    ## Formula: values ~ SSlogis(days, Asym, xmid, scal)
    ## 
    ## Parameters:
    ##      Estimate Std. Error t value Pr(>|t|)    
    ## Asym  42.1459     2.0795  20.267  < 2e-16 ***
    ## xmid   3.9661     0.1616  24.541  < 2e-16 ***
    ## scal   1.0532     0.1252   8.413 6.79e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.925 on 26 degrees of freedom
    ## 
    ## Number of iterations to convergence: 4 
    ## Achieved convergence tolerance: 1.38e-06

``` r
# Set colors for each group
group_colors <- c("c1" = "red", "c2" = "blue", "n" = "green")

# Increase the plot margins
par(mar = c(5, 5, 4, 8))  # Adjust the values as needed

# Plot the original data points and fitted curves for each group
plot(NULL, xlim = c(min(df_long$days), max(df_long$days)), ylim = c(0, max(df_long$values)),
     xlab = "Days", ylab = "Cell number(10^5)", main = "Logistic Growth Model")

for (group in groups) {
  group_data <- df_long[df_long$group == group, ]
  
  # Fit logistic growth model for the current group
  model <- nls(values ~ SSlogis(days, Asym, xmid, scal), data = group_data)
  
  if (nrow(group_data) > 0) {
    # Generate a sequence of x values for plotting
    x <- seq(min(group_data$days), max(group_data$days), length.out = 100)
    
    # Compute the corresponding y values using the logistic growth model parameters
    y <- predict(model, newdata = data.frame(days = x))
    
    # Plot the original data points
    points(group_data$days, group_data$values, pch = 16, col = adjustcolor(group_colors[group], alpha = 0.2))
    
    
    
    # Plot the fitted curve
    lines(x, y, type = "l", col = group_colors[group], lwd = 2)
  }
}

# Add legend to the right and outside the plot
par(xpd = TRUE)  # Allow plotting outside the plot area
legend("topright", legend = groups, col = group_colors, pch = 16, lwd = 2, bty = "n",
       xjust = 1, yjust = 1, inset = c(-0.2, 0), xpd = TRUE)
```

![](curve_fitting_files/figure-gfm/original%20data%20and%20fitting%20curve-1.png)<!-- -->

``` r
par(xpd = FALSE)  # Reset the xpd parameter

# Reset the plot margins
par(mar = c(5, 4, 4, 2.5))  # Adjust the values as needed
```

``` r
# Install the "colorblindr" package if not already installed
# install.packages("colorblindr")

# Load the required library
#library(colorblindr)


file_name <- "model_fitting.png"

# Open the PNG graphics device
png(file_name)

# Define color-blind safe colors
group_colors <- c("c1" = "#E69F00", "c2" = "#56B4E9", "n" = "#009E73")

# Define dot symbols for each group
group_symbols <- c("c1" = 1, "c2" = 2, "n" = 3)

# Increase the plot margins
par(mar = c(5, 5, 4, 8))  # Adjust the values as needed

# Plot the original data points and fitted curves for each group
plot(NULL, xlim = c(min(df_long$days), max(df_long$days)), ylim = c(0, max(df_long$values)),
     xlab = "Days", ylab = "Cell number(10^5)", main = "Logistic Growth Model")

for (group in groups) {
  group_data <- df_long[df_long$group == group, ]
  
  # Fit logistic growth model for the current group
  model <- nls(values ~ SSlogis(days, Asym, xmid, scal), data = group_data)
  
  if (nrow(group_data) > 0) {
    # Generate a sequence of x values for plotting
    x <- seq(min(group_data$days), max(group_data$days), length.out = 100)
    
    # Compute the corresponding y values using the logistic growth model parameters
    y <- predict(model, newdata = data.frame(days = x))
    
    # Plot the original data points with the assigned dot symbol
    points(group_data$days, group_data$values, pch = group_symbols[group], col = adjustcolor(group_colors[group], alpha = 0.2))
    
    # Plot the fitted curve
    lines(x, y, type = "l", col = group_colors[group], lwd = 2)
  }
}

# Add legend to the right and outside the plot
par(xpd = TRUE)  # Allow plotting outside the plot area
legend("topright", legend = groups, col = group_colors, pch = group_symbols, lwd = 2, bty = "n",
       xjust = 1, yjust = 1, inset = c(-0.2, 0), xpd = TRUE)
par(xpd = FALSE)  # Reset the xpd parameter

# Reset the plot margins
par(mar = c(5, 4, 4, 2.5))  # Adjust the values as needed

dev.off()
```

    ## quartz_off_screen 
    ##                 2

### Play around for image style

``` r
# Define a color palette for the groups
group_palette <- c("c1" = "red", "c2" = "blue", "n" = "green")

# Create the base R plot
plot(summary_data_long$Days, summary_data_long$mean, type = "p", col = group_palette[summary_data_long$Group], pch = 16,
     xlab = "Days", ylab = "Cell number(10^5)", main = "Logistic Growth Model")

# Add error bars
apply(summary_data_long, 1, function(row) {
  x <- as.numeric(row["Days"])
  y <- as.numeric(row["mean"])
  sem <- as.numeric(row["sem"]) 
  group <- row["Group"]
  
  segments(x, y - sem, x, y + sem, col = group_palette[group], lwd = 2)
})
```

    ## NULL

``` r
# Add legend
legend("topright", legend = unique(summary_data_long$Group), col = group_palette[unique(summary_data_long$Group)], pch = 16)
```

![](curve_fitting_files/figure-gfm/sem%20from%20r%20base%20plot-1.png)<!-- -->

``` r
# Set the file name and path for the exported image
file_name <- "with_r_square_fitting_datapoints.png"

# Open the PNG graphics device
png(file_name)

# Increase the plot margins
par(mar = c(5, 5, 4, 8))  # Adjust the values as needed

# Plot the original data points and fitted curves for each group
plot(NULL, xlim = c(min(df_long$days), max(df_long$days)), ylim = c(0, max(df_long$values)),
     xlab = "Days", ylab = "Cell number(10^5)", main = "Logistic Growth Model")

# Initialize the y-coordinate for the R-squared labels
text_y <- c(40.8,37.8,34.8) # 40.8 is the maximum y value, and assign the 

for (i in 1:length(groups)) {
  group <- groups[i]
  group_data <- df_long[df_long$group == group, ]
  
  # Fit logistic growth model for the current group
  model <- tryCatch({
    nls(values ~ SSlogis(days, Asym, xmid, scal), data = group_data)
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(model)) {
    # Generate a sequence of x values for plotting
    x <- seq(min(group_data$days), max(group_data$days), length.out = 100)
    
    # Compute the corresponding y values using the logistic growth model parameters
    y <- predict(model, newdata = data.frame(days = x))
    
    # Plot the original data points with the assigned dot symbol
    points(group_data$days, group_data$values, pch = group_symbols[group], col = adjustcolor(group_colors[group], alpha = 0.5))
    
    # Plot the fitted curve
    lines(x, y, type = "l", col = group_colors[group], lwd = 2)
    
    # Calculate R-squared
    ssr <- sum((group_data$values - predict(model))^2)
    sst <- sum((group_data$values - mean(group_data$values))^2)
    r_squared <- 1 - ssr / sst
    
    # Add group information and R-squared to the left corner of the graph for each group
    text(x = min(group_data$days), y = text_y[i], 
         labels = paste0(group, ": R^2 = ", round(r_squared, 2)),
         pos = 4, adj = c(0, 1), col = group_colors[group])
    
  }
}

# Add legend to the right and outside the plot
par(xpd = TRUE)  # Allow plotting outside the plot area
legend("topright", legend = groups, col = group_colors, pch = group_symbols, lwd = 2, bty = "n",
       xjust = 1, yjust = 1, inset = c(-0.2, 0), xpd = TRUE)
par(xpd = FALSE)  # Reset the xpd parameter

# Reset the plot margins
par(mar = c(5, 4, 4, 2.5))  # Adjust the values as needed

# Close the PNG graphics device and save the plot
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## the final image employed

## Image exported

``` r
# Install the "colorblindr" package if not already installed
# install.packages("colorblindr")

# Load the required library
#library(colorblindr)
# Set the file name and path for the exported image
file_name <- "with_r_square_fitting_and_mean.png"

# Open the PNG graphics device
png(file_name,res = 120)



# Define color-blind safe colors
group_colors <- c("c1" = "#E69F00", "c2" = "#56B4E9", "n" = "#009E73")

# Define dot symbols for each group
group_symbols <- c("c1" = 1, "c2" = 16, "n" = 8)

# Increase the plot margins
par(mar = c(5, 4, 4, 2.5))  # Adjust the values as needed

# Plot the original data points and fitted curves for each group
plot(NULL, xlim = c(min(df_long$days), max(df_long$days)), ylim = c(0, max(df_long$values)),
     xlab = "Days", ylab = "Cell number(10^5)", main = "Logistic Growth Model")

# Initialize the y-coordinate for the R-squared labels
text_y <- c(40.8,37.8,34.8) # 40.8 is the maximum y value, and assign the

for (i in 1:length(groups)) {
  group <- groups[i]
  group_data <- df_long[df_long$group == group, ]
  
  
  # Fit logistic growth model for the current group
  model <- nls(values ~ SSlogis(days, Asym, xmid, scal), data = group_data)
  
  if (nrow(group_data) > 0) {
    # Generate a sequence of x values for plotting
    x <- seq(min(group_data$days), max(group_data$days), length.out = 100)
    
    # Compute the corresponding y values using the logistic growth model parameters
    y <- predict(model, newdata = data.frame(days = x))
    
    # # Plot the original data points with the assigned dot symbol
    # points(group_data$days, group_data$values, pch = group_symbols[group], col = adjustcolor(group_colors[group], alpha = 0.2))
    # 
    # Create the base R plot
    subset_summary <- summary_data_long[summary_data_long$Group == group, ]
    points(subset_summary$Days, subset_summary$mean, pch = group_symbols[group],col = adjustcolor(group_colors[group], alpha = 0.5))
    
    
    
    
    # Plot the fitted curve
    lines(x, y, type = "l", col = group_colors[group], lwd = 2)
    
    
    # Calculate R-squared
    ssr <- sum((group_data$values - predict(model))^2)
    sst <- sum((group_data$values - mean(group_data$values))^2)
    r_squared <- 1 - ssr / sst
    
    # Add group information and R-squared to the left corner of the graph for each group
    text(x = min(group_data$days), y = text_y[i], 
         labels = paste0(group, ": R^2 = ", round(r_squared, 2)),
         pos = 4, adj = c(0, 1), col = group_colors[group])
    
    
    
    
  }
}

# Add legend to the right and outside the plot
par(xpd = TRUE)  # Allow plotting outside the plot area
legend("topright", legend = groups, col = group_colors, pch = group_symbols, lwd = 2, bty = "n",
       xjust = 1, yjust = 1, inset = c(-0.3, 0), xpd = TRUE)
par(xpd = FALSE)  # Reset the xpd parameter

# Reset the plot margins
par(mar = c(5, 4, 4, 2.5))  # Adjust the values as needed

# Close the PNG graphics device and save the plot
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## Model parameters

``` r
library(nlme)
```

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

``` r
#model <- gnls(values ~ SSlogis(days, Asym, xmid, scal), data = df_long, params=list(Asym ~ group))


#Soy=as.data.frame(Soybean)
fm1 <- gnls(values ~ SSlogis(days, Asym, xmid, scal), df_long)
summary(fm1)
```

    ## Generalized nonlinear least squares fit
    ##   Model: values ~ SSlogis(days, Asym, xmid, scal) 
    ##   Data: df_long 
    ##        AIC      BIC    logLik
    ##   438.1674 447.6955 -215.0837
    ## 
    ## Coefficients:
    ##         Value Std.Error  t-value p-value
    ## Asym 37.10295 1.0340273 35.88198       0
    ## xmid  3.53568 0.0945682 37.38759       0
    ## scal  0.84576 0.0757271 11.16856       0
    ## 
    ##  Correlation: 
    ##      Asym  xmid 
    ## xmid 0.718      
    ## scal 0.644 0.424
    ## 
    ## Standardized residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -3.5878918 -0.4693795  0.1835391  0.5112357  2.2815654 
    ## 
    ## Residual standard error: 3.628134 
    ## Degrees of freedom: 80 total; 77 residual

``` r
fm1$coefficients
```

    ##       Asym       xmid       scal 
    ## 37.1029503  3.5356776  0.8457624

``` r
fm2 <- gnls(values ~ SSlogis(days, Asym, xmid, scal), df_long, 
            start= list(Asym = c(37, 37,37),
                        xmid= c(3.54, 3.54,3.54),
                        scal= c(0.8457,0.8457,0.8457)),
            params = list(Asym ~ group - 1, 
                          xmid ~ group - 1, 
                          scal ~ group - 1))
summary(fm2)
```

    ## Generalized nonlinear least squares fit
    ##   Model: values ~ SSlogis(days, Asym, xmid, scal) 
    ##   Data: df_long 
    ##        AIC      BIC    logLik
    ##   430.5931 454.4133 -205.2965
    ## 
    ## Coefficients:
    ##                 Value Std.Error   t-value p-value
    ## Asym.groupc1 37.16346 1.5504936 23.968792       0
    ## Asym.groupc2 33.37835 1.2733485 26.213049       0
    ## Asym.groupn  42.14592 2.3767738 17.732407       0
    ## xmid.groupc1  3.50103 0.1380689 25.357103       0
    ## xmid.groupc2  3.20112 0.1391587 23.003384       0
    ## xmid.groupn   3.96608 0.1847128 21.471578       0
    ## scal.groupc1  0.70171 0.1066762  6.577922       0
    ## scal.groupc2  0.71599 0.1116151  6.414792       0
    ## scal.groupn   1.05318 0.1430755  7.360974       0
    ## 
    ##  Correlation: 
    ##              Asym.1 Asym.2 Asym.g xmd.g1 xmd.g2 xmd.gr scl.g1 scl.g2
    ## Asym.groupc2 0.000                                                  
    ## Asym.groupn  0.000  0.000                                           
    ## xmid.groupc1 0.658  0.000  0.000                                    
    ## xmid.groupc2 0.000  0.602  0.000  0.000                             
    ## xmid.groupn  0.000  0.000  0.859  0.000  0.000                      
    ## scal.groupc1 0.594  0.000  0.000  0.394  0.000  0.000               
    ## scal.groupc2 0.000  0.556  0.000  0.000  0.310  0.000  0.000        
    ## scal.groupn  0.000  0.000  0.757  0.000  0.000  0.632  0.000  0.000 
    ## 
    ## Standardized residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -2.9118775 -0.5507306  0.1457179  0.4927899  2.4301872 
    ## 
    ## Residual standard error: 3.343243 
    ## Degrees of freedom: 80 total; 71 residual

``` r
fm2$coefficients
```

    ## Asym.groupc1 Asym.groupc2  Asym.groupn xmid.groupc1 xmid.groupc2  xmid.groupn 
    ##   37.1634591   33.3783465   42.1459197    3.5010285    3.2011214    3.9660762 
    ## scal.groupc1 scal.groupc2  scal.groupn 
    ##    0.7017078    0.7159874    1.0531752

``` r
summary(anova(fm1, fm2))
```

    ##      call               Model            df            AIC       
    ##  Length:2           Min.   :1.00   Min.   : 4.0   Min.   :430.6  
    ##  Class :character   1st Qu.:1.25   1st Qu.: 5.5   1st Qu.:432.5  
    ##  Mode  :character   Median :1.50   Median : 7.0   Median :434.4  
    ##                     Mean   :1.50   Mean   : 7.0   Mean   :434.4  
    ##                     3rd Qu.:1.75   3rd Qu.: 8.5   3rd Qu.:436.3  
    ##                     Max.   :2.00   Max.   :10.0   Max.   :438.2  
    ##                                                                  
    ##       BIC            logLik           Test      L.Ratio         p-value        
    ##  Min.   :447.7   Min.   :-215.1         :1   Min.   :19.57   Min.   :0.003296  
    ##  1st Qu.:449.4   1st Qu.:-212.6   1 vs 2:1   1st Qu.:19.57   1st Qu.:0.003296  
    ##  Median :451.1   Median :-210.2              Median :19.57   Median :0.003296  
    ##  Mean   :451.1   Mean   :-210.2              Mean   :19.57   Mean   :0.003296  
    ##  3rd Qu.:452.7   3rd Qu.:-207.7              3rd Qu.:19.57   3rd Qu.:0.003296  
    ##  Max.   :454.4   Max.   :-205.3              Max.   :19.57   Max.   :0.003296  
    ##                                              NA's   :1       NA's   :1

``` r
library(emmeans)
emmeans(fm2, pairwise~group, param="Asym")
```

    ## $emmeans
    ##  group emmean   SE df lower.CL upper.CL
    ##  c1      37.2 1.55 70     34.1     40.3
    ##  c2      33.4 1.27 70     30.8     35.9
    ##  n       42.1 2.38 70     37.4     46.9
    ## 
    ## Confidence level used: 0.95 
    ## 
    ## $contrasts
    ##  contrast estimate   SE df t.ratio p.value
    ##  c1 - c2      3.79 2.01 70   1.887  0.1501
    ##  c1 - n      -4.98 2.84 70  -1.756  0.1921
    ##  c2 - n      -8.77 2.70 70  -3.252  0.0050
    ## 
    ## P value adjustment: tukey method for comparing a family of 3 estimates

``` r
emmeans(fm2, pairwise~group, param="xmid")
```

    ## $emmeans
    ##  group emmean    SE df lower.CL upper.CL
    ##  c1      3.50 0.138 70     3.23     3.78
    ##  c2      3.20 0.139 70     2.92     3.48
    ##  n       3.97 0.185 70     3.60     4.33
    ## 
    ## Confidence level used: 0.95 
    ## 
    ## $contrasts
    ##  contrast estimate    SE df t.ratio p.value
    ##  c1 - c2     0.300 0.196 70   1.530  0.2833
    ##  c1 - n     -0.465 0.231 70  -2.017  0.1157
    ##  c2 - n     -0.765 0.231 70  -3.308  0.0042
    ## 
    ## P value adjustment: tukey method for comparing a family of 3 estimates
