##load necessary packages
library(lubridate) #work on dates
library(dplyr)
library(ggpubr)
library(viridisLite)
library(ggbeeswarm)
library(lmtest)
library(emmeans)
library(multcomp)
library(multcompView)
library(rstatix)
library(stringr)
library(vegan)
library(RColorBrewer)
library(iNEXT)

##basic dataset formatting for simulation analysis

MKE_transect<-read.csv("MKE_transect.csv")
#filter out surveys with uncertain transects
MKE_transect_only<-filter(MKE_transect, !Description=="")
MKE_transect_only$Taxon[MKE_transect_only$Taxon=="Celastrina"]<-"Celastrina ladon/neglecta"
#remove non species level id
sp<-taxo[!grepl("sp.",taxo$BAMONA.Name), ]
MKE_transect_only<-filter(MKE_transect_only, Taxon %in% sp$BAMONA.Name)
# Convert 'collection_date' to Date type
MKE_transect_only$teventDate<-as.Date(MKE_transect_only$eventDate,format=format)
# Calculate week number
MKE_transect_only$week_number <- format(MKE_transect_only$teventDate, "%U")
MKE_transect_only$week_number<-as.numeric(MKE_transect_only$week_number)
week_numbers <- unique(MKE_transect_only$week_number)
#get sampling matrix 
MKE_transect_only %>% distinct(Taxon)
MKE_transect_only %>% count(Taxon,sort=TRUE) #table like S1--number of transects each taxon appears in 
MKE_transect_only %>% distinct(eventDate, sort=TRUE)
MKE_transect_only %>% filter(eventDate > 2023-06-07)
###standardize dates
week_interval<-weeks(1)
biweek_interval<-weeks(2)
month_interval<-months(1)
MKE_transect_only<-MKE_transect_only %>% mutate(Site %in% c("KK River Trail North","KK River Trail South")) %>% separate(eventDate, into = c("start_datetime", "end_datetime"), sep = " - ")
MKE_transect_only$timestamp <- parse_date_time(MKE_transect_only$start_datetime, orders = c("ymd HM", "ymd", "mdy"))
MKE_transect_only$timestamp<-ymd(MKE_transect_only$timestamp)
class(MKE_transect_only$timestamp)
#create presence absence column
MKE_transect_only$Presence_Absence <- "1" #all presence values in dataframe
#trim for matrix
df<-data.frame(MKE_transect_only$Site,MKE_transect_only$Taxon,MKE_transect_only$Description,MKE_transect_only$timestamp,MKE_transect_only$Presence_Absence)
colnames(df)<-c("Site","Taxon","Transect","timestamp","Presence_Absence")
df_dedupe<-distinct(df) #removes multiples from different habitats, all cabbage whites and does not impact transect values
result_matrix <- df_dedupe %>%
  pivot_wider(
    id_cols = c(Transect, Site, timestamp),
    names_from = Taxon,
    values_from = Presence_Absence,
   ) %>%
      as.matrix()
result_matrix <- replace(result_matrix,is.na(result_matrix), 0)
# Print the result
print(result_matrix)
reference_table <- expand.grid(
  Site = unique(result_df$Site),
  Transect = unique(result_df$Transect),
  timestamp = unique(result_df$timestamp)
)
# Merge with the reference table to include all transects
sampled_data <- merge(reference_table, result_df, by = c("Site", "Transect", "timestamp"), all = TRUE)
write.csv(sampled_data, "data_matrix.csv")
data_matrix<-read.csv("data_matrix.csv",header=TRUE)
##fix dates
data_matrix <- replace(data_matrix,is.na(data_matrix), 0)
format<-"%m/%d/%Y"
# Convert 'collection_date' to Date type
data_matrix$timestamp<-as.Date(data_matrix$timestamp,format=format)
# Calculate week number
data_matrix$week_number <- format(data_matrix$timestamp, "%U")
data_matrix$week_number<-as.numeric(data_matrix$week_number)
week_numbers <- unique(data_matrix$week_number)

####simulations
###INTENSITY
  # Example #test with just the loops
  set.seed(123)
  ###new function for recreating routes with T1-4: short length, different number of surveys
  # Function to create subsamples at the specified time interval
  create_intensity_subsample <- function(df, intensity) {
    # conditional for size of intensity
    if (intensity <= 0) {
      stop("Number of samples must be a positive integer.")
    }
    
    # Identify unique transects
    unique_transects <- unique(df$Transect)
    # Identify unique week numbers
    unique_weeks <- unique(df$week_number)
    # If num_samples is greater than the number of transects, set it to the maximum possible value
      if (intensity > length(unique_weeks)) {
        intensity <- length(unique_weeks)
      }
    
    
    # Initialize a dataframe to store subsamples
    subsamples <- data.frame()
    
      # Select a random week number for each transect
      selected_rows <- lapply(unique_transects, function(transect) {
        df_sub <- df[df$Transect == transect, ]
        # Ensure that the intensity does not exceed the number of rows in df_sub
        if (intensity > nrow(df_sub)) {
          intensity <- nrow(df_sub)
        }
        random_row <- df_sub[sample(nrow(df_sub), size = intensity,replace=FALSE), ]
        return(random_row)
      })
      
      # Combine the selected rows into a single dataframe
      selected_rows_df <- do.call(rbind, selected_rows)
      
      subsamples <- rbind(subsamples, selected_rows_df)
    
    
    return(subsamples)
  }


###keep num of surveys the same (just one per week) just increase plot size per each week
  #test with MU site 4, repeat for each site
  MU<-site_list$`Marquette University`
  #KK site 3
  KKS<-site_list$`KK River Trail South`
  EB<-site_list$Estabrook #no T4
  KKN<-site_list$`KK River Trail North`
  EB$week_number<-as.character(EB$week_number)
  KKS$week_number<-as.character(KKS$week_number)
  # Function to create subsamples at the specified time interval
  create_intensity_subsample <- function(df, intensity) {
    # conditional for size of intensity
    if (intensity <= 0) {
      stop("Number of samples must be a positive integer.")
    }
    # If num_samples is greater than the number of rows, set it to the maximum possible value
    
      # Filter data for the specified time interval
  week <- df[df$week_number %in% week_number, ]  
    # Initialize a dataframe to store subsamples
    subsamples <- data.frame()
    subsample<-week[sample(nrow(week), size = intensity, replace = FALSE), ]
    subsample<-subsample[order(subsample$week_number),]
    subsamples <- rbind(subsamples, subsample)

    
    return(subsamples)
  }

  # List to store subsamples
subsamples_by_week<-list()
week_numbers <- unique(MU$week_number)#change for site
for(week_number in week_numbers){
  # Identify unique week numbers
    subsamples_by_intensity<-list()
    num_subsamples=1000
    intensities<-c(1:4)#change for number of transects in each site
    for(intensity in intensities){
      #list
      subsamples_intensity_sizes<-list()
      # Create subsamples at the specified time interval
      # Loop through subsamples
      for (i in 1:num_subsamples) {   
        subsample_intensity_sizes <- create_intensity_subsample(MU, intensity)#change for site
        # Bind the subsamples within the time interval 
        subsamples_intensity_sizes[[paste("Subsample",i)]] <- subsample_intensity_sizes
      }
      subsamples_by_intensity[[paste("Intensity",intensity)]] <-subsamples_intensity_sizes
    } 
    
    #Store subsamples by week for this site
    subsamples_by_week[[paste("Week", week_number)]] <- subsamples_by_intensity
  }

#remove duplicates within each site subsamples
unique_subsamples_by_week<-list()
# Function to sort rows within each data frame
# Function to keep data frames with the same rows regardless of order
unique_dataframes_two_levels_nested_lists <- function(lst) {
  # Sort rows within data frames in nested lists
  sorted_nested_lists <- sort_rows_within_nested_lists(lst)
  
  # Convert each data frame to a character matrix to make them comparable
  char_matrices <- lapply(sorted_nested_lists, function(outer_inner_list) {
    lapply(outer_inner_list, function(inner_list) {
      lapply(inner_list, as.matrix)
    })
  })
  
  # Convert character matrices to characters to make them hashable
  char_strings <- lapply(char_matrices, function(outer_inner_list) {
    lapply(outer_inner_list, function(inner_list) {
      lapply(inner_list, function(mat) {
        apply(mat, 1, paste, collapse = "")
      })
    })
  })
  
  # Find unique character strings
  unique_indices <- lapply(char_strings, function(outer_inner_list) {
    lapply(outer_inner_list, function(inner_list) {
      which(!duplicated(inner_list))
    })
  })
  
  # Return unique data frames within nested lists
  unique_nested_lists <- lapply(seq_along(unique_indices), function(i) {
    lapply(seq_along(unique_indices[[i]]), function(j) {
      lapply(unique_indices[[i]][[j]], function(unique_index) {
        if (length(unique_index) == 0) {
          NULL
        } else {
          sorted_nested_lists[[i]][[j]][[unique_index]]
        }
      })
    })
  })
  
  return(unique_nested_lists)
}

# Keep only unique data frames within two-level nested lists
unique_subsamples_by_week <- unique_dataframes_two_levels_nested_lists(subsamples_by_week)
#transfer names
names(unique_subsamples_by_week)<-names(subsamples_by_week)
intensities<-names(subsamples_by_week$`Week 25`) #intensity values
names(unique_subsamples_by_week$`Week 35`)<-intensities
names(unique_subsamples_by_week$`Week 34`)<-intensities
names(unique_subsamples_by_week$`Week 33`)<-intensities
names(unique_subsamples_by_week$`Week 32`)<-intensities
names(unique_subsamples_by_week$`Week 31`)<-intensities
names(unique_subsamples_by_week$`Week 29`)<-intensities
names(unique_subsamples_by_week$`Week 28`)<-intensities
names(unique_subsamples_by_week$`Week 27`)<-intensities
names(unique_subsamples_by_week$`Week 26`)<-intensities
names(unique_subsamples_by_week$`Week 23`)<-intensities
names(unique_subsamples_by_week$`Week 37`)<-intensities
names(unique_subsamples_by_week$`Week 25`)<-intensities
names(unique_subsamples_by_week$`Week 30`)<-intensities

##all named, calculate richness for each week each plot size (intensity)
# Create an empty list to store the results--deduplicated for T1-Tn, repeat for each site
results_intensity_plot_s4<-list()
results_intensity_plot_s3<-list()
results_intensity_plot_s1<-list()
results_intensity_plot_s2<-list()


# Define the process_data_frame function
process_data_frame <- function(df, list_name,week_name) {
  # Process the data frame here
  col_sums <- colSums(df[, 5:20])
  presence <- ifelse(col_sums > 0, 1, 0)
  rich_sum <- sum(presence)
  
  # Create a data frame containing the list name and sum
  richness_estimate <- data.frame(week_number=week_name,List_Name = list_name, Rich_Sum = rich_sum)
  
  return(richness_estimate)
}

# Loop through all the lists and data frames 
for (i in 1:length(unique_subsamples_by_week)) {
  current_list <- unique_subsamples_by_week[[i]]
  for (j in 1:length(current_list)) {
    intensity_list<-current_list[[j]]
    for(k in 1:length(intensity_list)){
      df <- intensity_list[[k]]
      list_name <- names(current_list)[j]  # Get the name of the data frame within the list
      
      week_name<-names(unique_subsamples_by_week)[i] #Get name of transect list
      # Use lapply to apply the process_data_frame function to each data frame within the intensity list
      intensity_results <- lapply(intensity_list, function(df){
        list_name = list_name
        week_name=week_name
        result <- process_data_frame(df, list_name, week_name)
        return(result)
      })
      
    }
    results_intensity_plot_s4 <- c(results_intensity_plot_s4, intensity_results)##change for each site
  }
}

######run through above for each site before combining
# Combine the results into a single data frame
result_df_intensity_plot_s4 <- do.call(rbind, results_intensity_plot_s4)
result_df_intensity_plot_s3 <- do.call(rbind, results_intensity_plot_s3)
result_df_intensity_plot_s1 <- do.call(rbind, results_intensity_plot_s1)
result_df_intensity_plot_s2 <- do.call(rbind, results_intensity_plot_s2)
# Print the result data frame
print(result_df_intensity_plot_s2)
result_df_intensity_plot_s1<-result_df_intensity_plot_s1%>% group_by(week_number,List_Name)%>%
  mutate(average=mean(Rich_Sum)) #calculate average for each Site/Intensity
#best_rich_est<-distinct(filter(result_df_intensity,List_Name.x %in% "12")) #filter to highest...can just use Intensity 12 for all
best_rich_est#get dataframe of richness per transect
site4_best_rich_est<-filter(best_rich_est, Site_Name %in% "Site 4")
site3_best_rich_est<-filter(best_rich_est, Site_Name %in% "Site 3")
site1_best_rich_est<-filter(best_rich_est, Site_Name %in% "Site 1")
site2_best_rich_est<-filter(best_rich_est, Site_Name %in% "Site 2")

# Calculate the richness ratio
result_df_intensity_plot_s2 <- result_df_intensity_plot_s2 %>%
  mutate(RichnessRatio = Rich_Sum / site2_best_rich_est$average)
result_df_intensity_plot_s1 <- result_df_intensity_plot_s1 %>%
  mutate(RichnessRatio = Rich_Sum / site1_best_rich_est$average)
result_df_intensity_plot_s3 <- result_df_intensity_plot_s3 %>%
  mutate(RichnessRatio = Rich_Sum / site3_best_rich_est$average)
result_df_intensity_plot_s4 <- result_df_intensity_plot_s4 %>%
  mutate(RichnessRatio = Rich_Sum / site4_best_rich_est$average)

# Fit a negative binomial regression model
nb_model <- glm.nb(RichnessRatio~List_Name,data=result_df_intensity_plot_s2)
# Perform the Likelihood Ratio Test
lr_test <- lrtest(poisson, nb_model)
# Print the test results
print(lr_test)

#overdispersion found
qpoisson<-glm(RichnessRatio~List_Name+Site ,family="quasipoisson",data=result_df_intensity_plot)
summary(qpoisson)
marginal = emmeans(qpoisson,~
                     List_Name+Site)
cld<-cld(marginal,
         alpha=0.05,
         Letters=letters,  ### Use lower-case letters for .group
         adjust="tukey")
##combo plot
result_df_intensity_plot_s1$Site<-"Site 1"
result_df_intensity_plot_s2$Site<-"Site 2"
result_df_intensity_plot_s3$Site<-"Site 3"
result_df_intensity_plot_s4$Site<-"Site 4"
result_df_intensity_plot<-rbind(result_df_intensity_plot_s1,result_df_intensity_plot_s2,result_df_intensity_plot_s3,result_df_intensity_plot_s4)
result_df_intensity_plot$List_Name[result_df_intensity_plot$List_Name=="Intensity 1"]<-"500m"
result_df_intensity_plot$List_Name<-factor(result_df_intensity_plot$List_Name,levels=c("500m", "1km", "1.5km", "2km"))
colnames(result_df_intensity_plot)[2]<-"Length"

# Calculate sample sizes
sample_sizes <- result_df_intensity_plot %>%
  group_by(Site,List_Name) %>%
  summarise(Sample_Size = n())
sample_sizes$Value<-0.1
##get medians and standard deviation for curves to compare intensity
plot_data <- result_df_intensity_plot %>%
  group_by(intensity,Site_Name,Length) %>%
  summarize(sd_y = sd(RichnessRatio),median_y=median(RichnessRatio))



###transects with different numbers of surveys
# Function to create subsamples at the specified time interval
  create_intensity_subsample <- function(df, intensity) {
    # conditional for size of intensity
    if (intensity <= 0) {
      stop("Number of samples must be a positive integer.")
    }
      # If num_samples is greater than the number of rows, set it to the maximum possible value
    intensity <- min(intensity, nrow(df))
    # Initialize a dataframe to store subsamples
    subsamples <- data.frame()
      subsample<-df[sample(nrow(df), size = intensity, replace = FALSE), ]
      subsample<-subsample[order(subsample$week_number),]
        subsamples <- rbind(subsamples, subsample)



    return(subsamples)
 }
  #subsample_list <- list()
  data_matrix_clean<-filter(data_matrix,!X==100)##remove EB T4, no need to change full route bc it did not add any new species
  # Split the data by site
  site_list <- split(data_matrix_clean, data_matrix_clean$Site)
  
  # Initialize a list to store subsamples within each site
  subsamples_by_site <- list()
  
  # Loop through sites
  for (site_id in 1:length(site_list)) {
    site_data <- site_list[[site_id]]
    
    # Split the site data by transect
    transect_list <- split(site_data, site_data$Transect)
    
    # Initialize a list to store subsamples within each transect
    subsamples_by_transect <- list()
    
    # Loop through transects within the site
    for (transect_id in 1:length(transect_list)) {
      transect_data <- transect_list[[transect_id]]
      #Intensity for subsampling (e.g., 1-8 times)
      # Adjust 'time_interval' as needed
      subsamples_by_intensity<-list()
      intensities <- c(1,2,3,4,6,8,9,10,12,14,15,16,18,20,22,24)  
      num_subsamples=1000
      for(intensity in intensities){
        #list
        subsamples_intensity_interval<-list()
        # Create subsamples at the specified time interval
        # Loop through subsamples
        for (i in 1:num_subsamples) {   
          subsample_intensity_interval <- create_intensity_subsample(transect_data, intensity)
          # Bind the subsamples within the time interval 
          subsamples_intensity_interval[[paste("Subsample",i)]] <- subsample_intensity_interval
          }
        subsamples_by_intensity[[paste("Intensity",intensity)]] <-subsamples_intensity_interval
        }
      subsamples_by_transect[[paste("_Site", site_id, "_Transect", transect_id)]] <- subsamples_by_intensity
      }
    #Store subsamples by transect for this site
    subsamples_by_site[[paste("Site", site_id)]] <- subsamples_by_transect
  }

#remove duplicates within each site
# Keep only unique data frames within three-level nested lists
#each site individually for per transect
unique_subsamples_by_site1<-list()
unique_subsamples_by_site1 <- unique_dataframes_two_levels_nested_lists(subsamples_by_site$`Site 1`)
unique_subsamples_by_site2 <- unique_dataframes_two_levels_nested_lists(subsamples_by_site$`Site 2`)
unique_subsamples_by_site3 <- unique_dataframes_two_levels_nested_lists(subsamples_by_site$`Site 3`)
unique_subsamples_by_site4 <- unique_dataframes_two_levels_nested_lists(subsamples_by_site$`Site 4`)
##assign names to new unique lists
Transects<-names(subsamples_by_site$`Site 4`)
intensities #intensity values
names(unique_subsamples_by_site4)<-Transects
names(unique_subsamples_by_site4$`_Site 4 _Transect 4`) <-intensities
#list
unique_all<-list(unique_subsamples_by_site1,unique_subsamples_by_site2,unique_subsamples_by_site3,unique_subsamples_by_site4)
names(unique_all)<-c("Site 1", "Site 2", "Site 3", "Site 4")
##calculate richness
# # Create an empty list to store the results, one site at a time
results_intensity_transect_s4 <- list()
unique_all_s4<-unique_all$`Site 4`#
# # Define the process_data_frame function
process_data_frame <- function(df, list_name,transect_name) {
  #   # Process the data frame here
  #   # Process the data frame here
  col_sums <- colSums(df[, 5:20])
  presence <- ifelse(col_sums > 0, 1, 0)
  rich_sum <- sum(presence)
  #
  #   # Create a data frame containing the list name and sum
  richness_estimate <- data.frame(Transect_Name=transect_name,List_Name = list_name, Rich_Sum = rich_sum)
  #
  return(richness_estimate)
}
#
# # Loop through all the lists and data frames
for (i in 1:length(unique_all_s4)) {
  current_list <- unique_all_s4[[i]]
  for (j in 1:length(current_list)) {
    intensity_list<-current_list[[j]]
    for(k in 1:length(intensity_list)){
      df <- intensity_list[[k]]
      list_name <- names(current_list)[j]  # Get the name of the data frame within the list
      #
      transect_name<-names(unique_all_s4)[i] #Get name of transect list
      #       # Use lapply to apply the process_data_frame function to each data frame within the intensity list
      intensity_results <- lapply(intensity_list, function(df){
        list_name = list_name
        transect_name=transect_name
        result <- process_data_frame(df, list_name, transect_name)
        return(result)
      })
      #
    }
    results_intensity_transect_s4 <- c(results_intensity_transect_s4, intensity_results)
  }
}

# # Combine the results into a single data frame
result_df_s4 <- do.call(rbind, results_intensity_transect_s4)
#calculate for each df
result_df_s4<-result_df_s4%>% group_by(Transect_Name,List_Name)%>%
  mutate(average=mean(Rich_Sum)) #calculate average for each Transect/Intensity
#or combo df and add sites
result_df_s1$Site_Name<-"Site 1"
result_df_s2$Site_Name<-"Site 2"
result_df_s3$Site_Name<-"Site 3"
result_df_s4$Site_Name<-"Site 4"

result_df_all<-rbind(result_df_s1,result_df_s2,result_df_s3,result_df_s4)
best_rich_est#get dataframe of richness per site!
result_df_all <- result_df_all%>%
  left_join(best_rich_est, by = "Site_Name","List_Name")
result_df_all<-result_df_all[-c(6:7,9:10)]
# Calculate the richness ratio
result_df_all <- result_df_all %>%
  mutate(RichnessRatio = Rich_Sum.x / average.y)
####full route
#combine transect data for full route
  head(data_matrix)
site_matrix<-data_matrix_clean %>%
    group_by(Site,week_number)%>%
  summarise(
    CP=sum(Colias.philodice),
    PR=sum(Pieris.rapae),
    DP=sum(Danaus.plexippus),
    PP=sum(Papilio.polyxenes),
    Pop=sum(Polites.peckius),
    EV=sum(Euphyes.vestris),
    EC=sum(Epargyreus.clarus),
    VA=sum(Vanessa.atalanta),
    CLN=sum(Celastrina.ladon.neglecta),
    CeP=sum(Cercyonis.pegala),
    HP=sum(Hylephila.phyleus),
    EB=sum(Erynnis.baptisiae),
    VV=sum(Vanessa.virginiensis),
    PI=sum(Polygonia.interrogationis),
    AN=sum(Ancyloxypha.numitor),
    PG=sum(Papilio.glaucus)
  )
  
  # Example #test with just the loops
  set.seed(123)
  
  # Function to create subsamples at the specified time interval
  create_intensity_subsample_full <- function(df, intensity) {
    # conditional for size of intensity
    if (intensity <= 0) {
      stop("Number of samples must be a positive integer.")
    }
    # If num_samples is greater than the number of rows, set it to the maximum possible value
    intensity <- min(intensity, nrow(df))
    # Initialize a dataframe to store subsamples
    subsamples <- data.frame()
    subsample<-df[sample(nrow(df), size = intensity, replace = FALSE), ]
    subsample<-subsample[order(subsample$week_number),]
    subsamples <- rbind(subsamples, subsample)
    
    
    
    return(subsamples)
  } 
  
  # List to store subsamples
  subsample_list_full <- list()
  
  # Split the data by site
  site_list <- split(site_matrix, site_matrix$Site)
  
  # Initialize a list to store subsamples within each site
  subsamples_by_site_full <- list()
  
  # Loop through sites
  for (site_id in 1:length(site_list)) {
    site_data <- site_list[[site_id]]
    
    
      # Intensity for subsampling (e.g., 1-8 times)
      # Adjust 'time_interval' as needed
      subsamples_by_intensity_full<-list()
      intensities <- c(1,2,3,4,5,6,8,10,12)  # Change to 4 for one sample every month
      num_subsamples=100
      for(intensity in intensities){
        #list
        subsamples_intensity_interval_full<-list()
        # Create subsamples at the specified time interval
        # Loop through subsamples
        for (i in 1:num_subsamples) {   
          subsample_intensity_interval_full <- create_intensity_subsample_full(site_data, intensity)
          # Bind the subsamples within the time interval 
          subsamples_intensity_interval_full[[paste("Subsample",i)]] <- subsample_intensity_interval_full
        }
        subsamples_by_intensity_full[[paste("Intensity",intensity)]] <-subsamples_intensity_interval_full
      }
      subsamples_by_site_full[[paste("_Site", site_id)]] <- subsamples_by_intensity_full}
####intensity by length--full routes
#remove duplicates within each site
unique_subsamples_by_site_full<-list()
# Function to sort rows within each data frame
# Function to keep data frames with the same rows regardless of order
unique_dataframes_two_levels_nested_lists <- function(lst) {
  # Sort rows within data frames in nested lists
  sorted_nested_lists <- sort_rows_within_nested_lists(lst)
  
  # Convert each data frame to a character matrix to make them comparable
  char_matrices <- lapply(sorted_nested_lists, function(outer_inner_list) {
    lapply(outer_inner_list, function(inner_list) {
      lapply(inner_list, as.matrix)
    })
  })
  
  # Convert character matrices to characters to make them hashable
  char_strings <- lapply(char_matrices, function(outer_inner_list) {
    lapply(outer_inner_list, function(inner_list) {
      lapply(inner_list, function(mat) {
        apply(mat, 1, paste, collapse = "")
      })
    })
  })
  
  # Find unique character strings
  unique_indices <- lapply(char_strings, function(outer_inner_list) {
    lapply(outer_inner_list, function(inner_list) {
      which(!duplicated(inner_list))
    })
  })
  
  # Return unique data frames within nested lists
  unique_nested_lists <- lapply(seq_along(unique_indices), function(i) {
    lapply(seq_along(unique_indices[[i]]), function(j) {
      lapply(unique_indices[[i]][[j]], function(unique_index) {
        if (length(unique_index) == 0) {
          NULL
        } else {
          sorted_nested_lists[[i]][[j]][[unique_index]]
        }
      })
    })
  })
  
  return(unique_nested_lists)
}
# Function to sort rows within each data frame
sort_rows_within_dataframe <- function(df) {
  df_sorted <- df[do.call(order, df), ]
  return(df_sorted)
}

# Function to apply sorting to data frames within nested lists
sort_rows_within_nested_lists <- function(lst) {
  sorted_outer_list <- lapply(lst, function(outer_inner_list) {
    lapply(outer_inner_list, function(inner_list) {
      lapply(inner_list, sort_rows_within_dataframe)
    })
  })
  return(sorted_outer_list)
}

# Keep only unique data frames within two-level nested lists
unique_subsamples_by_site_full <- unique_dataframes_two_levels_nested_lists(subsamples_by_site_full)
##note this function is the same for all double nested lists

names(unique_subsamples_by_site_full)<-c("Site 1", "Site 2", "Site 3", "Site 4")
intensities #intensity values
names(unique_subsamples_by_site_full$`Site 1`)<-intensities
names(unique_subsamples_by_site_full$`Site 2`)<-intensities
names(unique_subsamples_by_site_full$`Site 3`)<-intensities
names(unique_subsamples_by_site_full$`Site 4`)<-intensities

# Create an empty list to store the results--deduplicated full routes
results_intensity_full <- list()

# Define the process_data_frame function
process_data_frame_full <- function(df, list_name,site_name) {
  # Process the data frame here
  col_sums <- colSums(df[, 3:18])
  presence <- ifelse(col_sums > 0, 1, 0)
  rich_sum <- sum(presence)
  
  # Create a data frame containing the list name and sum
  richness_estimate <- data.frame(Site_Name=site_name,List_Name = list_name, Rich_Sum = rich_sum)
  
  return(richness_estimate)
}

# Loop through all the lists and data frames 
for (i in 1:length(unique_subsamples_by_site_full)) {
  current_list <- unique_subsamples_by_site_full[[i]]
  for (j in 1:length(current_list)) {
    intensity_list<-current_list[[j]]
    for(k in 1:length(intensity_list)){
      df <- intensity_list[[k]]
      list_name <- names(current_list)[j]  # Get the name of the data frame within the list
      
      site_name<-names(unique_subsamples_by_site_full)[i] #Get name of transect list
      # Use lapply to apply the process_data_frame function to each data frame within the intensity list
      intensity_results <- lapply(intensity_list, function(df){
        list_name = list_name
        site_name=site_name
        result <- process_data_frame_full(df, list_name, site_name)
        return(result)
      })
      
    }
    results_intensity_full <- c(results_intensity_full, intensity_results)
  }
}


# Combine the results into a single data frame
result_df_intensity_full <- do.call(rbind, results_intensity_full)

# Print the result data frame
print(result_df_intensity_full)
result_df_intensity_full<-result_df_intensity_full%>% group_by(Site_Name,List_Name)%>%
  mutate(average=mean(Rich_Sum)) #calculate average for each Site/Intensity
best_rich_est<-distinct(filter(result_df_intensity_full,List_Name %in% "12")) #filter to highest...can just use Intensity 12 for all
best_rich_est#get dataframe of richness per transect
result_df_intensity_full <- result_df_intensity_full%>%
  left_join(best_rich_est, by = "Site_Name","List_Name")
result_df_intensity_full<-result_df_intensity_full[-c(5:6)]
# Calculate the richness ratio
result_df_intensity_full <- result_df_intensity_full %>%
  mutate(RichnessRatio = Rich_Sum.x / average.y)
#factor variables for plot
result_df_intensity$List_Name.x<-factor(result_df_intensity$List_Name.x,levels=unique(result_df_intensity$List_Name.x))
result_df_intensity_full$List_Name<-factor(result_df_intensity_full$List_Name,levels=unique(result_df_intensity_full$List_Name))
names<-names(result_df_intensity_mixed)
names(result_df_intensity_full)<-names
result_df_intensity_mixed$List_Name<-factor(result_df_intensity_mixed$List_Name,levels=unique(result_df_intensity_mixed$List_Name))
result_df_intensity_full <- colnames(result_df_intensity_full)[colnames(result_df_intensity_full)  == "Length.x"] <- "RichnessRatio"

#compare the intensity vs. time question
result_df_all$Length<-"500m"
result_df_intensity_full$Length<- ifelse(result_df_intensity_full$Site_Name=="Site 2","1km",ifelse(result_df_intensity_full$Site_Name == "Site 4", "2km", "1.5km"))

result_df_intensity_full$List_Name.x<-as.character(result_df_intensity_full$List_Name.x)         ##reassign intensity values in full to match the transects
# filter mixed dataset for comparison
transect_comp<-filter(result_df_all,List_Name.x %in% c("1","2","3","4","5","6","7","8","9","10","11","12"))
#match intensity names, redo for each intensity case
transect_comp<-transect_comp %>%
  mutate(intensity = case_when(
    Site_Name %in% c("Site 1","Site 3","Site 4") & List_Name.x %in% "2" ~ "NA",
    Site_Name %in% "Site 2" & List_Name.x %in% "2" ~ "1",
    Site_Name %in% c("Site 1","Site 3") & List_Name.x %in% "3" ~ "1",
    Site_Name %in% c("Site 2","Site 4") & List_Name.x %in% "3" ~ "NA",
    Site_Name %in% c("Site 1","Site 3") & List_Name.x %in% "4" ~ "NA",
    Site_Name %in% "Site 2" & List_Name.x %in% "4" ~ "2",
    Site_Name %in% "Site 4" & List_Name.x %in% "4" ~ "1",
    Site_Name %in% c("Site 1","Site 3") & List_Name.x %in% "6" ~ "2",
    Site_Name %in% "Site 2" & List_Name.x %in% "6" ~ "3",
    Site_Name %in% "Site 4" & List_Name.x %in% "6" ~ "NA",
    Site_Name %in% c("Site 1","Site 3") & List_Name.x %in% "8" ~ "NA",
    Site_Name %in% "Site 2" & List_Name.x %in% "8" ~ "4",
    Site_Name %in% "Site 4" & List_Name.x %in% "8" ~ "2",
    Site_Name %in% c("Site 1","Site 3") & List_Name.x %in% "9" ~ "3",
    Site_Name %in% "Site 2" & List_Name.x %in% "9" ~ "NA",
    Site_Name %in% "Site 4" & List_Name.x %in% "9" ~ "NA",
    Site_Name %in% c("Site 1","Site 3") & List_Name.x %in% "12" ~ "4",
    Site_Name %in% "Site 2" & List_Name.x %in% "12" ~ "6",
    Site_Name %in% "Site 4" & List_Name.x %in% "12" ~ "3",
    Site_Name %in% c("Site 1","Site 3") & List_Name.x %in% "15" ~ "5",
    Site_Name %in% c("Site 2","Site 4") & List_Name.x %in% "15" ~ "NA",
    Site_Name %in% "Site 4" & List_Name.x %in% "16" ~ "4",
    Site_Name %in% "Site 3" & List_Name.x %in% "18" ~ "6",
    Site_Name %in% c("Site 1","Site 2","Site 4") & List_Name.x %in% "18" ~ "NA",
    Site_Name %in% "Site 4" & List_Name.x %in% "20" ~ "5",
    Site_Name %in% c("Site 1","Site 2","Site 3") & List_Name.x %in% "20" ~ "NA",
    Site_Name %in% "Site 4" & List_Name.x %in% "24" ~ "6",
    Site_Name %in% "Site 3" & List_Name.x %in% "24" ~ "8",
    Site_Name %in% c("Site 1","Site 2") & List_Name.x %in% "20" ~ "NA",
     ))
transect_comp_filter<-filter(transect_comp,!intensity %in% "NA")#remove NA values for filtered dataset
###combined datasets
result_df_intensity_full<-result_df_intensity_full %>% mutate(intensity=List_Name.x) #add column for comparison

intensity_comp<-full_join(transect_comp,result_df_intensity_full)
qpoisson<-glm(RichnessRatio~Length+List_Name.x+Site_Name,family="quasipoisson",data=intensity_comp)
summary(qpoisson)

##calculate for curves
# Calculate median and standard deviation
sd_data <- intensity_comp %>%
  group_by(List_Name.x,Site_Name,Length) %>%
  summarize(sd_y = sd(RichnessRatio),median_y=median(RichnessRatio))
sd_data_sites<-intensity_comp%>%
  group_by(intensity,Length) %>%
  summarize(sd_y = sd(RichnessRatio),median_y=median(RichnessRatio))
  sd_data_filter <- sd_data %>%
    filter((Site_Name %in% "Site 1" & List_Name.x %in% c("1", "2", "3", "4", "5")) |
             (Site_Name %in% "Site 2" & List_Name.x %in% c("1", "2", "3", "4", "5", "6")) |
             (Site_Name %in% "Site 3" & List_Name.x %in% c("1", "2", "3", "4", "5", "6", "8")) |
             (Site_Name %in% "Site 4" & List_Name.x %in% c("1", "2", "3", "4", "5", "6", "8", "9", "10", "12")))
#####Species accumulation plots
#specaccum all the transects
MUT1_sa <- specaccum(MUT1, method = "collector")
MUT2_sa <- specaccum(MUT2, method = "collector")
MUT3_sa <- specaccum(MUT3, method = "collector")
MUT4_sa <- specaccum(MUT4, method = "collector")
#create dataframe for ggplot
data_MUT1 <- data.frame(Sites=MUT1_sa$sites, Richness=MUT1_sa$richness)
data_MUT1['Transect']='1'
data_MUT1<- tibble::rownames_to_column(data_MUT1, "rn")
data_MUT2 <- data.frame(Sites=MUT2_sa$sites, Richness=MUT2_sa$richness)
data_MUT2['Transect']='2'
data_MUT2<- tibble::rownames_to_column(data_MUT2, "rn")
data_MUT3 <- data.frame(Sites=MUT3_sa$sites, Richness=MUT3_sa$richness)
data_MUT3['Transect']='3'
data_MUT3<- tibble::rownames_to_column(data_MUT3, "rn")
data_MUT4 <- data.frame(Sites=MUT4_sa$sites, Richness=MUT4_sa$richness)
data_MUT4['Transect']='4'
data_MUT4<- tibble::rownames_to_column(data_MUT4, "rn")

data_MU<-rbind(data_MUT1,data_MUT2,data_MUT3,data_MUT4)
data_MU

#################diversity plot and rarefactions
#produce table with counts per species per date per site--only for Site 4
Site4_matrix<-filter(site_matrix, Site %in% "Marquette University")
Site4_matrix<-as.data.frame(Site4_matrix)
Site4_matrix$week_number<-c("X23","X25","X26","X27","X28","X29","X31","X32","X33","X34","X35","X37")
Site4_matrix
Site4_matrix_biweekly<-filter(Site4_matrix,week_number %in% c("X25","X27","x29","X31","X33","X35","X37"))
rownames(Site4_matrix_biweekly)<-Site4_matrix_biweekly$week_number
Site4_matrix_biweekly<-Site4_matrix_biweekly[3:18]
Site4_matrix
Site4_matrix_biweekly<-data.matrix(Site4_matrix_biweekly)
Site4_matrix_biweekly<-ifelse(Site4_matrix_biweekly >0,1,0)
write.csv(Site4_matrix,"Site4.csv")
Site4_matrix<-read.csv("Site4.csv",header=TRUE)
Site4_matrix<-as.data.frame(Site4_matrix)
Site4_matrix_monthly<-filter(Site4_matrix,week_number %in% c("X25","X29","X33","X37"))
rownames(Site4_matrix_monthly)<-Site4_matrix_monthly$week_number
Site4_matrix_monthly<-Site4_matrix_monthly[3:18]
Site4_matrix
Site4_matrix_monthly<-data.matrix(Site4_matrix_monthly)
Site4_matrix_monthly<-ifelse(Site4_matrix_monthly >0,1,0)
Site4_matrix_bme<-filter(Site4_matrix,week_number %in% c("X25","X31","X37"))##beginning, middle, and end surveys
rownames(Site4_matrix_bme)<-Site4_matrix_bme$week_number
Site4_matrix_bme<-Site4_matrix_bme[3:18]
Site4_matrix_bme<-data.matrix(Site4_matrix_bme)
Site4_matrix_bme<-ifelse(Site4_matrix_bme >0,1,0)
Site4_matrix_be<-filter(Site4_matrix,week_number %in% c("X25","X37"))#beginning and end surveys
rownames(Site4_matrix_be)<-Site4_matrix_be$week_number
Site4_matrix_be<-Site4_matrix_be[3:18]
Site4_matrix_be<-data.matrix(Site4_matrix_be)
Site4_matrix_be<-ifelse(Site4_matrix_be >0,1,0)
Site4_matrix_mid<-filter(Site4_matrix,week_number %in% "X31")#midpoint survey
rownames(Site4_matrix_mid)<-Site4_matrix_mid$week_number
Site4_matrix_mid<-Site4_matrix_mid[3:18]
Site4_matrix_mid<-data.matrix(Site4_matrix_mid)
Site4_matrix_mid<-ifelse(Site4_matrix_mid >0,1,0)
rownames(Site4_matrix)<-Site4_matrix$week_number
Site4_matrix<-Site4_matrix[3:18]
Site4_matrix<-data.matrix(Site4_matrix)
Site4_matrix<-ifelse(Site4_matrix >0,1,0)
Site4<-list(Site4_matrix,Site4_matrix_biweekly,Site4_matrix_monthly,Site4_matrix_bme,Site4_matrix_be,Site4_matrix_mid)#has to be list regardless
list_names<-c("Weekly","Biweekly","Monthly","BME","B/E","Midpoint")
Site4_t<-lapply(Site4,t)
names(Site4_t)<-list_names#has to be list regardless
out.raw <- iNEXT(Site4_t, q = 0, datatype="incidence_raw") iNEXT function for rarefaction
out.raw.rich<-filter(out.raw$AsyEst,Diversity %in% "Species richness")
###small plot (transect 2) from same site
Site4_t2_matrix<-filter(data_matrix, Site %in% "Marquette University")
Site4_t2_matrix<-filter(Site4_t2_matrix, 
                        Transect %in% "T2")

Site4_t2_matrix<-as.data.frame(Site4_t2_matrix)
Site4_t2_matrix$week_number<-c("X23","X25","X26","X27","X28","X29","X31","X32","X33","X34","X35","X37")
Site4_t2_matrix
Site4_t2_matrix_biweekly<-filter(Site4_t2_matrix,week_number %in% c("X25","X27","x29","X31","X33","X35","X37"))
rownames(Site4_t2_matrix_biweekly)<-Site4_t2_matrix_biweekly$week_number
Site4_t2_matrix_biweekly<-Site4_t2_matrix_biweekly[3:18]
Site4_t2_matrix
Site4_t2_matrix_biweekly<-data.matrix(Site4_t2_matrix_biweekly)
Site4_t2_matrix_biweekly<-ifelse(Site4_t2_matrix_biweekly >0,1,0)
write.csv(Site4_t2_matrix,"Site4.csv")
Site4_t2_matrix<-read.csv("Site4.csv",header=TRUE)
Site4_t2_matrix<-as.data.frame(Site4_t2_matrix)
Site4_t2_matrix_monthly<-filter(Site4_t2_matrix,week_number %in% c("X25","X29","X33","X37"))
rownames(Site4_t2_matrix_monthly)<-Site4_t2_matrix_monthly$week_number
Site4_t2_matrix_monthly<-Site4_t2_matrix_monthly[3:18]
Site4_t2_matrix
Site4_t2_matrix_monthly<-data.matrix(Site4_t2_matrix_monthly)
Site4_t2_matrix_monthly<-ifelse(Site4_t2_matrix_monthly >0,1,0)
Site4_t2_matrix_bme<-filter(Site4_t2_matrix,week_number %in% c("X25","X31","X37"))
rownames(Site4_t2_matrix_bme)<-Site4_t2_matrix_bme$week_number
Site4_t2_matrix_bme<-Site4_t2_matrix_bme[3:18]
Site4_t2_matrix_bme<-data.matrix(Site4_t2_matrix_bme)
Site4_t2_matrix_bme<-ifelse(Site4_t2_matrix_bme >0,1,0)
Site4_t2_matrix_be<-filter(Site4_t2_matrix,week_number %in% c("X25","X37"))
rownames(Site4_t2_matrix_be)<-Site4_t2_matrix_be$week_number
Site4_t2_matrix_be<-Site4_t2_matrix_be[3:18]
Site4_t2_matrix_be<-data.matrix(Site4_t2_matrix_be)
Site4_t2_matrix_be<-ifelse(Site4_t2_matrix_be >0,1,0)
Site4_t2_matrix_mid<-filter(Site4_t2_matrix,week_number %in% "X31")
rownames(Site4_t2_matrix_mid)<-Site4_t2_matrix_mid$week_number
Site4_t2_matrix_mid<-Site4_t2_matrix_mid[3:18]
Site4_t2_matrix_mid<-data.matrix(Site4_t2_matrix_mid)
Site4_t2_matrix_mid<-ifelse(Site4_t2_matrix_mid >0,1,0)
rownames(Site4_t2_matrix)<-Site4_t2_matrix$week_number
Site4_t2_matrix<-Site4_t2_matrix[3:18]
Site4_t2_matrix<-data.matrix(Site4_t2_matrix)
Site4_t2_matrix<-ifelse(Site4_t2_matrix >0,1,0)
Site4_t2<-list(Site4_t2_matrix,Site4_t2_matrix_biweekly,Site4_t2_matrix_monthly,Site4_t2_matrix_bme,Site4_t2_matrix_be,Site4_t2_matrix_mid)#has to be list regardless
list_names<-c("Weekly","Biweekly","Monthly","BME","B/E","Midpoint")
Site4_t2_t<-lapply(Site4_t2,t)
names(Site4_t2_t)<-list_names#has to be list regardless
out.raw <- iNEXT(Site4_t2_t, q = 0, datatype="incidence_raw") #test

# Re-order the levels of the categorical variable
out.raw$AsyEst$Site <- factor(out.raw$AsyEst$Site, levels = c("Weekly", "Biweekly", "Monthly","BME","B/E","Midpoint"))
out.raw$DataInfo$site <- factor(out.raw$DataInfo$site, levels = c("Weekly", "Biweekly", "Monthly","BME","B/E","Midpoint"))

out.raw.rich<-filter(out.raw$AsyEst,Diversity %in% "Species richness")
plot+theme_classic(base_size=18)+theme(legend.position="bottom", legend.title=element_blank()) +scale_fill_viridis_d(option="D",aesthetics =c("color","fill"))+scale_color_viridis_d(option="D",aesthetics =c("color","fill"))+geom_hline(data=out.raw.rich,aes(yintercept = Observed, color=Site),linetype = "dashed",show.legend=TRUE)# Add asymptote line

###############effect of week on species presence, seasonality
head(result_df_all) ##get general richness and subset for equal transect data per each week
season<-filter(result_df_intensity_plot,Length
                   %in% "500m")
#glm for week number on richness estimates by transect
poisson<-glm(RichnessRatio~week_number+Site,data=seasons,family="poisson") #model tests effect of week number on richness ratio blocking sites
glmnb<-glm.nb(RichnessRatio~week_number+Site,data=seasons)
lr_test <- lrtest(poisson, glmnb)
# Print the test results
print(lr_test)
qpoisson<-glm(RichnessRatio~week_number+Site,data=season,family="quasipoisson") #model tests effect of week number on richness ratio blocking sites

summary(glmnb)
summary(qpoisson)
summary(glm)
summary(poisson)
marginal = emmeans(qpoisson,~
                     week_number|Site)
pairs(marginal,
      adjust="tukey")

cld<-cld(marginal,
         alpha=0.05,
         Letters=letters,  ### Use lower-case letters for .group
         adjust="tukey")

season$week_number<-factor(season$week_number,levels=c("Week 23","Week 25","Week 26","Week 27","Week 28","Week 29","Week 30","Week 31","Week 32","Week 33","Week 34","Week 35","Week 37"))

#Figures

####basic plot for phenology
ggplot(MKE_transect_only, aes(x = week_number, y = Count, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Abundance Counts of Different Taxa Over Different Dates",
       x = "Date",
       y = "Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
##increased area sampling and higher richness plot one site boxplot
ggplot(result_df_intensity_plot_s2,aes(x=List_Name,y=RichnessRatio,fill=List_Name))+scale_fill_viridis_d(option="D")+geom_boxplot(width=0.25,outlier.size = -.5, color="black",lwd=.75, alpha = 0.5) +
  geom_violin(trim = FALSE, width = 2, alpha = 0.5, position = position_dodge(width = 0.5), size = 0.75, color = "black")+
  ylim(0,1.0)+
  theme_minimal() +xlab("Sampling Intensity") +ylab("Similarity to Best Species Richness Estimate")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+guides(fill = guide_legend(override.aes = list(alpha = 1,color="black"))) + rremove("legend.title") + rremove("legend") #geom_text(data = cld, aes(x = List_Name, y = Value, label = .group),size = 3,position = position_dodge(width = 1))
##increased area same num of surveys all sites boxplot
ggplot(result_df_intensity_plot,aes(x=Length,y=RichnessRatio,fill=Length))+facet_wrap(~Site) +scale_fill_viridis_d(option="D")+geom_boxplot(width=0.25,outlier.size = -.5, color="black",lwd=.75, alpha = 0.5) +
  geom_violin(trim = FALSE, width = 2, alpha = 0.5, position = position_dodge(width = 0.5), size = 0.75, color = "black")+
  ylim(0,1.0)+
  theme_minimal() +xlab("Transect Length") +ylab("Similarity to Best Species Richness Estimate")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+guides(fill = guide_legend(override.aes = list(alpha = 1,color="black"))) +theme(text = element_text(size = 14))
#### intensity curves
ggplot(sd_data, aes(x = intensity, y = median_y,color=Length,fill=Length, group=Length)) +
  # Add data points
  facet_wrap(~Site_Name) +scale_color_viridis_d(option="D")+scale_fill_viridis_d(option="D")+geom_smooth(method="loess",lwd=1,se=FALSE) +geom_ribbon(data = sd_data, aes(ymin = median_y - sd_y, ymax = median_y + sd_y), alpha = 0.3) +theme_minimal() +xlab("Intensity") +ylab("Similarity to Best Species Richness Estimate")+theme(text = element_text(size = 16)) +coord_cartesian(ylim = c(0, 1.0))
ylim(0,1.0)
##boxplot for full routes and increasing surveying
result_df_intensity_event<-filter(result_df_intensity_full, List_Name.x %in% c("1","2","3","4","5","6","7","8"))
ggplot(result_df_intensity_event,aes(x=List_Name.x,y=RichnessRatio,fill=List_Name.x))+scale_fill_viridis_d(option="D")+geom_boxplot(width=0.25,outlier.size = -.5, color="black",lwd=.75, alpha = 0.5) +
  geom_violin(trim = FALSE, width = 2, alpha = 0.5, position = position_dodge(width = 0.5), size = 0.75, color = "black")+
  ylim(0,1.0)+
  theme_minimal() +xlab("Number of Survey Events") +ylab("Similarity to Best Species Richness Estimate")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+guides(fill = guide_legend(override.aes = list(alpha = 1,color="black"))) + rremove("legend.title") + rremove("legend") 
###saved
ggsave("Intensity_event.pdf",plot=last_plot())
##2km all sites--randomly select surveys
ggplot(result_df_intensity_full,aes(x=List_Name.x,y=RichnessRatio,fill=List_Name.x))+ facet_wrap(~Site_Name) +scale_fill_viridis_d(option="D")+geom_boxplot(width=0.25,outlier.size = -.5, color="black",lwd=.75, alpha = 0.5) +geom_violin(trim=FALSE,width=2,alpha=0.5, position = position_dodge(width=.5),size=.75,color="black") +ylim(0,1.0)+
  theme_minimal() +xlab("Sampling Intensity") +ylab("Similarity to Best Species Richness Estimate")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+guides(fill = guide_legend(override.aes = list(alpha = 1,color="black"))) + rremove("legend.title") + rremove("legend")#+geom_text(data = cld, aes(x = List_Name.x, y = emmean, label = cld$.group),size = 3,position=position_dodge(width=.75))
###rarefaction ggiNEXT
plot<-ggiNEXT(out.raw,type=1)
plot+theme_classic(base_size=16)+theme(legend.position="bottom", legend.title=element_blank()) +scale_fill_viridis_d(option="D",aesthetics =c("color","fill"))+scale_color_viridis_d(option="D",aesthetics =c("color","fill"))+geom_hline(data=out.raw.rich,aes(yintercept = Observed,color=Site),linetype = "dashed",show.legend=FALSE)# Add asymptote line
###weekly comparison plot
ggplot(season,aes(x=week_number,y=RichnessRatio,color=Site))+scale_color_viridis_d(option="D")+geom_boxplot(width=0.75,outlier.size = -.5, color="black",lwd=.75, alpha = 0.5) +geom_point(position="jitter",size=1)+
  theme_minimal() +xlab("Survey by Week") +ylab("Similarity to Best Species Richness Estimate")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))+guides(fill = guide_legend(override.aes = list(alpha = 1,color="black"))) + rremove("legend.title") +ylim(0,1.0)+theme(text = element_text(size = 14))
ggsave("Season_comp.pdf",plot=last_plot())
