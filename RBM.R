#Rules based model


##
## read in the background-corrected microbe data and sample qc metrics
##

overview <- read.csv("./sample_overviews.csv")
microbe_reports_bgc <- read.csv("./combined-reports.csv")

#Filter for significant taxa after background correction
sig_microbe_reports <- microbe_reports_bgc[microbe_reports_bgc$p_val < .01,] 

#Filter for bacteria and eukaryota
sig_microbe_reports %>% 
        filter(category == "bacteria" | category == "eukaryota") -> microbe_reports

samplenames <- unique(microbe_reports$sample_name)


##
## Rules based model
##

# function to get maximum drop-off from a single sample (identify the outliers)
get_difs <- function(l){
        difs <- list()
        for(i in seq_along(l[1:length(l)-1])){
                difs[[i]] <- l[i] - l[i+1]
        }
        return(unlist(difs))
}

# function to apply the RBM
apply_rbm <- function(df, overview, id, plot_results = TRUE, output_directory = "./outputs/RBM/"){
        #log_values <- log10(df + 1)                # log-transform rpm values
        log_values <- df                            # DON'T APPLY THE LOG-TRANSFORMATION
        values <- sort(log_values[,id][!is.na(log_values[,id])], decreasing = TRUE)
        if(length(values) > 15){
                values <-  values[1:15]                   # restrict analysis to the top 15 organisms   
        }
        x <- get_difs(values)                       # identify the index of the maximum drop-off in RPM values
        if(length(x) == 0){
                top_orgs <- values
        }else{
                top_orgs <- values[1:which.max(x)]
        }
        
        # plot the results for visual interogation
        if(plot_results){
                pdf(paste(output_directory, "RBM_",colnames(log_values)[id],".pdf", sep=""), height = 6, width = 6)
                par(mar = c(2, 10, 2, 2) + 0.2)
                title = paste(colnames(log_values)[id], overview[overview$sample_name == colnames(log_values)[id],c("initial_input_mass_ng")])
                barplot(values, las=2, horiz=TRUE, main = title, 
                        col= c(rep("red",length(top_orgs)), rep("blue", length(values)-length(top_orgs))), cex.names = .6, cex.main = .6)#, xlim=c(0,5))
                dev.off()
        }
        
        return(top_orgs)
}




##
## Apply the RBM to background-corrected data
##

result_info <- c()   # list of all results that will be output

# loop through all the samples and apply the RBM
for(sn in samplenames){
        print(sn)
        
        merged_report <- microbe_reports[microbe_reports$sample_name == sn,]
        
        merged_report %>%
                dplyr::filter(genus_tax_id > 0) %>%  # require genus taxid (removes "uncultured bacteria")
                dplyr::filter(nr_count > 0) %>%      # filter out taxons with 0 NR reads
                group_by(genus_tax_id) %>%           # group by genus, select the top species within the genus
                top_n(1, abs(nt_rpm)) %>% 
                as.data.frame() -> filtered_report
        
        if(dim(filtered_report)[1] > 0){
                mat <- acast(filtered_report, name~sample_name, value.var="nt_rpm")
                result_info <- c(result_info, sn)
                
                mass <- overview[overview$sample_name == sn,c("initial_input_mass_ng")]
                result_info <- c(result_info, mass)
                
                result <- apply_rbm(mat, overview, 1)
                result_info <- c(result_info, names(result))
                result_info <- c(result_info, "--")
                
        }
        
}

sink("./outputs/rbm_results.txt")
writeLines(unlist(lapply(result_info, paste, collapse=" ")))
sink()

