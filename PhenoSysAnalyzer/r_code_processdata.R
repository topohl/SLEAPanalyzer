## Read in data from Folders and xlsx data. The data is organized in subfolders. 
#  This will read in the data of different xlsx sheets and merge them based on the name of the folder. 
#  Also, the data will be sorted and formatted
#  Later, the data points will be plotted.

library(readxl)      # load readxl package for reading excel files
library(dplyr)       # load dplyr package for data manipulation functions
library(purrr)       # load purrr package for functional programming
library(stringr)     # load stringr package for string manipulation
library(lubridate)   # load lubridate package for date and time functions
library(f)

# set the working directory to the parent directory containing the subfolders
setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/RFID/B2")

# get a list of all subfolders in the directory
subfolders <- list.dirs(".", recursive = FALSE)

# function to read and process a single Excel file
process_file <- function(file_path) {
  read_excel(file_path) %>%         # read the excel file
    mutate(filename = basename(file_path)) %>%     # add a column with the filename
    mutate(Change = str_extract(filename, "CC\\d")) %>%   # extract the Change number from the filename
    select(-filename) # remove the filename column
}

# read all Excel files in the subfolders and combine them into a single data frame
data <- map_dfr(subfolders, ~ list.files(path = .x, pattern = "E9_SIS_B2_CC\\d_ActivityIndex.xlsx", full.names = TRUE) %>% 
                  map(process_file) %>% 
                  bind_rows())

# separate the Animal column into AnimalNumber and Cage columns
data <- separate(data, Animal, c("AnimalNum", "Cage"), sep = "_", remove = FALSE)

# convert the DateTime column to a datetime format
data$DateTime <- as.POSIXct(data$DateTime, format = "%d.%m.%Y %H:%M")

data <- data %>%
  rename(ActivityIndex = ActivyIndex)


# create a Phase column based on the time of day
data <- data %>%
  mutate(Phase = ifelse(
    format(DateTime, "%H:%M") >= "18:30" | format(DateTime, "%H:%M") < "06:30",
    "Active",
    "Inactive"
  ),
  PriorActive = ifelse(
    format(DateTime, "%H:%M") >= "16:30" & format(DateTime, "%H:%M") <= "18:30",
    "TRUE",
    "FALSE"
  ),
  Hour = floor(difftime(DateTime, first(DateTime), units = "hours")),
  RecentChange = Hour <= 2,
  Group = ifelse(AnimalNum %in% c("OR126", "OR127", "OR128", "OR129"), "CON", "SIS")
  )

# remove intermediate variables
rm(subfolders)



# Create a new column to represent the consecutive Active phases for each animal
data$ConsecActive <- with(data, ave(Phase, AnimalNum, FUN=function(x) {
  cumsum(c(0, diff(ifelse(x == "Active", 1, 0))) == 1)
}))
data$ConsecActive <- as.numeric(data$ConsecActive)

# Aggregate data by AnimalNum, ConsecActive, and Group
hourly_data <- aggregate(ActivityIndex ~ AnimalNum + ConsecActive + Group + Phase + Hour + RecentChange + PriorActive, data = data, FUN = mean)

# Aggregate data by AnimalNum, ConsecActive, and Group for each night
nightly_data <- aggregate(ActivityIndex ~ AnimalNum + ConsecActive + Group + Phase + RecentChange + PriorActive, data = data, FUN = mean)



################################################################################
## Line plot of Active and Inactive Phase based on Group
#  This is used to compare activity patterns during the active and inactive phase

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

group_cols<-c("#1e3791","#00ac8c")

# Define the phases to loop through
phases <- c("Active", "Inactive")

# Create an empty list to store the plots
plots <- list()

# Loop through each phase and create a plot
for (phase in phases) {
  
  # Filter the data for the current phase
  phase_data <- nightly_data %>% 
    filter(Phase == phase) %>% 
    filter(RecentChange == FALSE) # Only keep rows where RecentChange is FALSE
  
  # Convert ConsecActive to factor
  phase_data <- phase_data %>% 
    mutate(ConsecActive = as.factor(ConsecActive))
    
  # Create the plot
  plot <- ggplot(phase_data, aes(x = fct_inorder(ConsecActive), ActivityIndex, group = Group, colour = Group)) +
    geom_line(stat = "summary") +
    stat_summary(aes(fill = Group), fun = mean,
                 fun.min = function(x) mean(x) - sd(x), 
                 fun.max = function(x) mean(x) + sd(x), 
                 geom = "ribbon", 
                 alpha = 0.2, colour = NA) +
    scale_color_manual(name = NULL, values = group_cols) + 
    scale_fill_manual(name = NULL, values = group_cols) +
    theme_minimal_hgrid(12, rel_small = 1) +
    labs(title = bquote(~bold(.(paste("Activity Index")))), subtitle = paste("", phase, "Phase"),
         y = "Activity [a.u.]", x = "Day/Night", caption = "") +
    theme(plot.title = element_text(hjust = 0.5, face = "plain"),
          plot.subtitle = element_text(hjust = 0.5, face = "plain"),
          legend.position = "top",
          legend.justification = "right",
          legend.text = element_text(size = 9),
          legend.box.spacing = unit(1, "pt"),
          axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1),
          axis.title.x = element_blank())
  
  # Add the plot to the list of plots
  plots[[phase]] <- plot
  
}

# Arrange the plots using grid.arrange
grid.arrange(grobs = plots, nrow = 2)


################################################################################
## Compare Means between Groups during 2 hours preceding the Active Phase
#  This is used to quantify activity before the onset of the active phase

library(cowplot)
## Compare Means between Groups during 2 hours preceding the Active Phase

# Filter the data for PriorActive == TRUE
prior_active_data <- nightly_data %>% 
  filter(PriorActive == TRUE)

# Calculate the means and standard deviations of ActivityIndex for each group
means <- prior_active_data %>% 
  group_by(Group, AnimalNum) %>% 
  summarise(ActivityIndex = mean(ActivityIndex),
            sd_ActivityIndex = sd(ActivityIndex))

# Normality test for CTRL group
con_norm <- shapiro.test(means$ActivityIndex[means$Group == "CON"])
con_norm

# Normality test for SIS group
sis_norm <- shapiro.test(means$ActivityIndex[means$Group == "SIS"])
sis_norm

# Determine which test to use based on normality of both groups
if (con_norm$p.value >= 0.05 & sis_norm$p.value >= 0.05) {
  # If normal, use t-test
  t_test <- t.test(ActivityIndex ~ Group, data = prior_active_data)
  test_name <- "t.test"
  test_pval <- t_test$p.value
  # Determine the significance level for t-test
  sig_levels <- sprintf("%.3f", test_pval)
} else {
  # If not normal, use Wilcoxon rank-sum test
  wilcox_test <- wilcox.test(ActivityIndex ~ Group, data = prior_active_data)
  test_name <- "wilcox.test"
  test_pval <- wilcox_test$p.value
  # Determine the significance level for Wilcoxon rank-sum test
  sig_levels <- sprintf("%.3f", p.adjust(test_pval, method = "BH"))
}

# Create plot for individual values
p_indiv <- ggplot(prior_active_data, aes(x = Group, y = ActivityIndex, color = Group)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5) + # Add jittered dots
  scale_color_manual(values = group_cols) + # Set color scheme
  scale_y_continuous(expand=c(0.1,0.1)) +
  labs(title = "Activity", x = "Group", y = "Activity Index [a.u.]") + # Add labels
  theme_minimal_hgrid(12, rel_small = 1) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(hjust = 0.5, face = "plain"),
        axis.text.x = element_text(angle = 0, vjust = 1, size = 11, hjust = 0.45),
        axis.title.x = element_blank()) # Use a white and black theme

# Create plot for means
p_means <- ggplot(means, aes(x = Group, y = ActivityIndex, fill = Group, colour = Group)) +
  geom_jitter(size=4, alpha=0.5, width=0.2, shape=16) + # Add dodged columns for mean values
  geom_errorbar(aes(ymin = ActivityIndex - sd_ActivityIndex, ymax = ActivityIndex + sd_ActivityIndex), 
                position = position_dodge(width = 0.9), width = 0.2) + # Add error bars
  stat_summary(fun.min=function(z) {quantile(z,0.25)}, fun.max=function(z) {quantile(z,0.75)}, fun=mean, color="black", size=1, shape=16) + # Add indication for mean and error bars
  scale_fill_manual(values = group_cols) + # Set fill colors
  scale_color_manual(values = group_cols) +
  scale_y_continuous(expand=c(0.1,0.1)) +
  labs(title = "Mean Activity", x = "Group", y = "Activity Index [a.u.]") + # Add labels
  theme_minimal_hgrid(12, rel_small = 1) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(hjust = 0.5, face = "plain"),
        axis.text.x = element_text(angle = 0, vjust = 1, size = 11, hjust = 0.45),
        axis.title.x = element_blank()) # Use a white and black theme

# Add significance asterisks
p_means <- p_means + geom_signif(comparisons = list(c("CON", "SIS")), 
                                 test = test_name,
                                 map_signif_level = TRUE, 
                                 textsize = 4,
                                 annotations = sig_levels,
                                 y_position = max(means$ActivityIndex)+5)


# Arrange two plots side by side
plot_grid(p_indiv, p_means, ncol = 2, align = "v", labels = c("A", "B"))





