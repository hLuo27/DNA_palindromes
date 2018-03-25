library(dplyr)
library(ggplot2)
library(reshape2)
library(knitr)

observed_dat <- read.table('hcmv.data', header = TRUE)

# Part 1: Random Scatter

num_basepairs=229354

#Randomly choose 296 palindrome sites from the DNA sequence of 229,354 bases (without replacement)
random_dat <- data.frame(random_location = sample(1:num_basepairs,296, replace = FALSE))
random_dat = random_dat %>% arrange(random_location)

dat <- bind_cols(data.frame(palindrome = 1:296),observed_dat,random_dat)
head(dat, 5)

melted_data = melt(dat,id="palindrome")

random_scatter_plot <- ggplot(data = melted_data) +
  geom_line(aes(x = palindrome, y = value, color = variable)) +
  labs(x = "Palindrome Index", y = "Location") +
  ggtitle("Comparing Locations by Random Scatter")
random_scatter_plot

# Part 2: Counts

# Analyzing counts (See if it follows Poisson)
bin_width = 3000
bins = seq(from = 1, to = num_basepairs, by = bin_width) 
num_intervals = round(229534/bin_width)

binned_data = data.frame(interval = cut(observed_dat$location,breaks = bins,labels = 1:(length(bins)-1)))
binned_data = binned_data %>% group_by(interval) %>% summarise(num_palindromes = n())
head(binned_data,5)

#For X~Poisson(lambda_1) 
lambda_1 = 296/num_intervals #Number of palindromes divided by number of intervals

#Creates summary df with number of intervals containing a specified number of palindromes
palindrome_summary <- binned_data %>% group_by(num_palindromes) %>% summarise(observed_num_intervals = n())
palindrome_summary = palindrome_summary %>% mutate(expected_num_intervals = num_intervals*(ppois(num_palindromes,lambda_1)-ppois(num_palindromes-1,lambda_1)))
palindrome_summary = palindrome_summary %>% add_row(num_palindromes = '7+', observed_num_intervals = sum(palindrome_summary$observed_num_intervals[7:10]),
                              expected_num_intervals = sum(palindrome_summary$expected_num_intervals[7:10])) %>% slice(c(1:6,11))
palindrome_summary
melt_counts = melt(palindrome_summary, id = 'num_palindromes')
counts_plot <- ggplot(data = melt_counts, aes(num_palindromes, value)) +
  geom_bar(aes(fill = variable),position = "dodge",stat="identity")+
  labs(x = "Number of Palindromes", y = "Number of Intervals") +
  ggtitle("Number of Palindromes per Interval")
counts_plot

#Perform goodness-of-fit test
observed_test_statistic = sum(((palindrome_summary$observed_num_intervals-palindrome_summary$expected_num_intervals)^2)/palindrome_summary$expected_num_intervals)
observed_test_statistic
pchisq(observed_test_statistic,df=5,lower.tail=FALSE)

# Part 3: Spacing
spacing_dat <- dat %>% select(palindrome,location) %>% mutate(distance = c(diff(dat$location),NA)) 
# Note that the difference in row i is distance between location of ith an (i+1)th palindromes
head(spacing_dat,5)

bins = seq(0,6000,100)
binned_spacing_data = data.frame(interval = cut(spacing_dat$distance[-nrow(spacing_dat)],breaks = bins,labels = seq(0,5900,100)))
binned_spacing_data = binned_spacing_data %>% group_by(interval) %>% summarise(count = n())

lambda_2 = 296/num_basepairs
interval_widths = as.numeric(levels(binned_spacing_data$interval)[binned_spacing_data$interval])
binned_spacing_data = binned_spacing_data %>% mutate(expected_count = 296*(exp(-1*lambda_2*interval_widths)-exp(-1*lambda_2*(interval_widths+100))))

binned_spacing_data

melt_bins = melt(binned_spacing_data,id = "interval")
spacing_plot = ggplot(data = melt_bins, aes(x=interval,y=value)) +
  geom_bar(aes(fill = variable),position = "dodge",stat = "identity") +
  scale_x_discrete(name ="Length of Spacing", breaks = seq(0,6000,500)) +
  labs(y = "Number of Occurences") + 
  ggtitle("Distance between Consecutive Palindromes")
spacing_plot

# Now analyzing spacing between pairs of palindromes (See if follows Gamma(2,lambda_3))
pair_dat <- dat %>% mutate(distance = c(diff(dat$location,lag = 2),0,0)) # Difference in row i is distance between location of ith an (i+2)nd palindromes
head(pair_dat,5)

bins = seq(0,6000,100)
binned_pairs_data = data.frame(interval = cut(pair_dat$distance[-nrow(pair_dat)],breaks = bins, labels = seq(0,5900,100)))
binned_pairs_data = binned_pairs_data %>% group_by(interval) %>% summarise(count = n())
interval_widths = as.numeric(levels(binned_pairs_data$interval)[binned_pairs_data$interval])
lambda_3 = 296/num_basepairs #Same as lambda_2
binned_pairs_data = binned_pairs_data %>% mutate(expected_count = 296*(pgamma(interval_widths+100,shape=2,rate=lambda_3)-pgamma(interval_widths,shape=2,rate=lambda_3)))
binned_pairs_data

melt_pairs = melt(binned_pairs_data,id = "interval")
pairs_plot = ggplot(data = melt_pairs, aes(x=interval,y=value)) +
  geom_bar(aes(fill=variable),position="dodge",stat="identity") +
  scale_x_discrete(name ="Length of Spacing", breaks = seq(0,6000,500)) +
  labs(y = "Number of Occurences") + 
  ggtitle("Distance between Palindromes (2 Apart)")
pairs_plot

# Now analyzing spacing between triplets of palindromes (See if follows Gamma(3,lambda_3))
triplet_dat <- dat %>% mutate(distance = c(diff(dat$location,lag = 3),rep(0,times=3))) # Difference in row i is distance between location of ith an (i+3)rd palindromes
head(triplet_dat,5)

bins = seq(0,8000,100)
binned_triplet_data = data.frame(interval = cut(triplet_dat$distance[-nrow(triplet_dat)],breaks = bins, labels = seq(0,7900,100)))
binned_triplet_data = binned_triplet_data %>% group_by(interval) %>% summarise(count = n())
interval_widths = as.numeric(levels(binned_triplet_data$interval)[binned_triplet_data$interval])
binned_triplet_data = binned_triplet_data %>% mutate(expected_count = 296*(pgamma(interval_widths+100,shape=3,rate=lambda_3)-pgamma(interval_widths,shape=3,rate=lambda_3)))
binned_triplet_data

melt_triplets = melt(binned_triplet_data,id = "interval")
triplet_plot = ggplot(data = melt_triplets, aes(x=interval,y=value)) +
  geom_bar(aes(fill=variable),position="dodge",stat="identity") +
  scale_x_discrete(name ="Length of Spacing", breaks = seq(0,8000,1000)) +
  labs(y = "Number of Occurences") + 
  ggtitle("Distance between Palindromes (3 Apart)")
triplet_plot

# Part 4: Max
k = max(binned_data$num_palindromes) #Highest number of hits in any interval

#Perform hypothesis test
ppois(k-1,lambda = lambda_1) #is the probability of having fewer than k palindromes in one interval
1-(ppois(k-1,lambda = lambda_1)^num_intervals) # Probability that the max count over m intervals is greater than or equal to k