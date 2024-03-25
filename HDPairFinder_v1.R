# HDPairFinder
# This R script reads in the mzML or mzXML files
#               outputs csv files containing H-/D-labeled pairs
#               notes: The syntax of light and heavy refer to H-labeled and D-labeled compounds in the following code 
# Tingting Zhao, December 22, 2022
# Copyright @ The University of British Columbia

#------------------------------ Set the working directory ---------------------------------------
working_directory <- "D:/2022-07-05-HDPairFinder"

#------------------------------ Choose modules ----------------------------------------
run_pairPicking <- TRUE # TRUE: run the first module, extracting the H-/D-labeled pairs; 
                        # FALSE: skip the first module
run_alignment <- TRUE  # TRUE: run the second module, aligning the pairs across multiple samples
                        # FALSE: skip the second module
run_gapFilling <- TRUE # TRUE: run the third module, retrieving the missing pair in each sample
                        # FALSE: skip the third module
run_annotation <- TRUE   # TRUE: run the fourth module, annotating the compound
                        # FALSE: skip the fourth module

#-----------------------------     Set the parameters     ----------------------------------------
## 1) parameters for pair picking ##############
heavy_mz_tol <- 20  # m/z tolerance (in ppm) of mass difference between the  
                  # theoretical heavy m/z and experiential heavy m/z
rt_diff <- c(-0.2, 0.1)   # acceptable retention time difference (D-labeled - H-labeled): -0.2 ~ 0.1 min 
int_ratio <- c(0.4, 1.4)  # acceptable intensity ratio range (D-labeled / H-labeled): 0.4 ~ 1.4
cc_threshold <- 0.7       # threshold for cross correlation
run_inSourceFrag <- TRUE  # TRUE: use ISFrag to determine in-source fragments and remove them
                          # FALSE: do not identify in source fragments
                          # note: ISfrag takes more than half of the running time, if the users want to speed up the calculation, you can turn off ISFrag
## 2) parameters for alignment ##############
align_mz_tol <- 50    # precursor m/z tolerance (ppm)
align_rt_tol <- 0.5   # retention time tolerance (min)

## 3) parameters for gap filling ##############
gap_mz_tol <- 20      # mz_tol for gap filling (ppm)
gap_rt_tol <- 0.5     # rt_tol for gap filling (min)
int_threshold <- 1000  # intensity threshold for light and heavy compounds 

## 4) parameters for annotation ##############
anno_mz_tol <- 30     # relative mass tolerance (ppm) for the annotation

#-----------------------------  Functions -----------------------------------------------
## 1) extract EIC function ########
# rawlcms is the the result after applying xcmsRaw [  rawlcms <- xcmsRaw(lcmsFiles[num]) ]
# rt in min
EIC_matrix <- function(rawlcms, mz, rt, mz_tol=0.01, rt_tol=30){
        mzRange <- c(mz-mz_tol, mz + mz_tol)
        rtRange <- c(rt*60-rt_tol, rt*60+rt_tol)
        rawEIC <- rawEIC(rawlcms, mzrange= mzRange, rtrange = rtRange)
        eic_matrix <- cbind(rawlcms@scantime[rawEIC$scan],
                            rawEIC$intensity)
        return(eic_matrix)
}

##2) Peak smoothing function ########
peak_smooth <- function(x,level=2){
        n <- level
        if(length(x) < 2*n){
                return(x)
        } else if(length(unique(x))==1){
                return(x)
        } else{
                y <- vector(length=length(x))
                for(i in 1:n){
                        y[i] <- sum(c((n-i+2):(n+1),n:1)*x[1:(i+n)])/sum(c((n-i+2):(n+1),n:1))
                }
                for(i in (n+1):(length(y)-n)){
                        y[i] <-  sum(c(1:(n+1),n:1)*x[(i-n):(i+n)])/sum(c(1:(n+1),n:1))
                }
                for(i in (length(y)-n+1):length(y)){
                        y[i] <- sum(c(1:n,(n+1):(n+i-length(x)+1))*x[(i-n):length(x)])/sum(c(1:n,(n+1):(n+i-length(x)+1)))
                }
                return(y)
        }
}

#----------------------------- Main program       ----------------------------------------
#--------- 1st module: extraction of H-/D-labeled compounds ---------------
#### 1.1 Pick the peaks ####
run_peakPicking <- TRUE # TRUE: run the peak picking 
                        # FALSE: use customized skip the peak picking
setwd(working_directory)
remove_isfrag_pairs <- TRUE # remove the in source fragment

if (run_peakPicking){
        library(xcms)
        
        lcmsFiles <- list.files(pattern = ".mzML")
        
        for (i in 1:length(lcmsFiles)) {
                
                lcmsFile <- lcmsFiles[i]
                lc_ms <- xcmsSet(lcmsFile,
                                 method = "centWave", # centWave method to extract the feature
                                 ppm=20, # region of interest: 50 ppm HX recommend 
                                 peakwidth= c(10,120), 
                                 
                                 mzdiff= 0.005, # minimum difference in m/z for peaks with overlapping retention time
                                 
                                 snthresh=5, 
                                 prefilter= c(3, 500), # prefilter = c(k, I): mass trace are only retained if they contains 
                                 # at least k peaks with intensity > I
                                 integrate=1,
                                 
                                 noise=1000 # optional argument: peaks with intensity smaller than noise are omitted from ROI
                )
                
              
                # peak information: m/z, rt, peak area, peak height
                peaks <- as.data.frame(lc_ms@peaks)
                #peaks$rt <- peaks$rt/60 # convert retention time unit from second to minute
                # compoundName <- paste(round(peaks$mz,4),round(peaks$rt,2), sep = " / " ) # use two column to generate a new column (calculation or function)
                compoundName <- paste0("S",i) # use two column to generate a new column (calculation or function)
                # ft <- as.data.frame(cbind(compoundName, round(peaks$mz,4), round(peaks$rt,2), peaks$maxo)) 
                # into: peak area
                # maxo: peak height
                head(peaks[,c(1,4,5,6,9)])
                 ft <- as.data.frame(cbind(compoundName, peaks[,c(1,4,5,6,9)]))
                 colnames(ft)[1] <- "feature_index"
                 colnames(ft)[6] <- "Intensity"
                 for(j in 1:nrow(ft)){
                         ft$feature_index[j] <- paste0(ft$feature_index[j],"_",j)
                 }
                 # convert retention time from seconds to min
                 ft$rt <- ft$rt/60
                 ft$rtmin <- ft$rtmin/60
                 ft$rtmax <- ft$rtmax/60
                #colnames(ft) <- c("ID", "Precursor Mass", "Retention Time","Intensity") 
                sname <- strsplit(lcmsFile,"\\.")[[1]][1]
                write.csv(ft, paste0(sname,"_Raw.csv"), row.names = F)
                
        }
        Sys.time()
        # ft <- ft[ft$Intensity >=500,]
}

# identify in source fragments
if (run_inSourceFrag) {
        # use package
        library(ISFrag)
        files <- list.files(pattern = "Raw.csv")
        MSfiles <- list.files(pattern = ".mzML")
        parent_dir <- getwd()
        for (i in 1:length(files)) {
                # change the colnames of feature table 
                # note: 1) some column names have changed, remember to change it back in the output
                #       "feature_index"   
                #       "mz"   ---> "mz"
                #       "rt"   ---> "rt"
                #       "Intensity"
                #       2) retention time: the unit minute has been changed to second,
                #                          remember to change it back
                ft <- read.csv(files[i])
                # Sciex csv table contains some N/A, remove them 
                if (!run_peakPicking){
                        ft <- ft[ft$rt != "N/A",]
                }
                # convert retention time back to seconds
                customFT <- ft
                customFT$rt <- customFT$rt*60
                customFT$rtmin <- customFT$rtmin*60
                customFT$rtmax <- customFT$rtmax*60
                
                # MS2 assignment 
                # create a fold to save the corresponding raw mzML
                MS2directory_name <- strsplit(MSfiles[i],".mzML")[[1]]
                dir.create(MS2directory_name)
                MS2_dir <- paste0(parent_dir,"/",MS2directory_name)
                file.copy(MSfiles[i], MS2_dir)
                
                # this step takes 30 seconds
                featureTable <- ms2.assignment(MS2directory = MS2_dir, customFT = customFT)
                
                #Identification of ISF Features 
                # Identify level 3 in-source fragments.
                featureTable <- featureTable[,-1]
                level3 <- find.level3(MS1directory = MS2_dir, 
                                      MS1.files = MSfiles[i], 
                                      featureTable = featureTable, 
                                      type = "single")
                
                # Identify level 2 in-source fragments.
                level2 <- find.level2(ISFtable = level3)
                
                # Identify level 1 in-source fragments.
                level1 <- find.level1(ISF_putative = level2)
                results <- get.ISFrag.results(ISF_List = level1, featureTable = featureTable)
                isf_featuerTable <- cbind(ft$feature_index, results$FeatureTable)
                colnames(isf_featuerTable)[1] <- "feature_index"
                result <- isf_featuerTable[,c(1,2,3,6,12,13,14)]
                result$rt <- result$rt/60
                #result <- read.csv(list.files(pattern="ISFrag.csv")[1])
                #result <- result[,-8]
                for(k in 1:nrow(result)){
                        
                        if(result$Num_Level2[k] == 0 & result$Num_Level1[k] == 0){next}
                        
                        result[k,]
                        levels <- strsplit(result$ISF_level[k], ";")[[1]]
                        levels <- unique(levels)
                        for(l in 1:length(levels)){
                                
                                level <- levels[l]
                        
                                level <- paste0("S",i,"_",level)
                                levels[l] <- level
                        }
                        result$ISF_level[k] <- paste0(levels, collapse = ";")
                        
                }
                setwd(parent_dir)
                newName <- strsplit(files[i], "\\Raw.csv")[[1]][1]
                write.csv(result, paste0(newName,"ISFrag.csv"), row.names = F)
        }
}


#### 1.2 Pick the pairs #####
# remove_C isotope_pairs <- TRUE # 
weight_mz <- 0.5 # weight of precursor m/z when calculating the pair picking score
weight_rt <- 0.25     # weight of retention time when calculating the pair picking score
weight_ratio <- 0.25  # weight of intensity ratio when calculating the pair picking score
score_threshold <- 0.2 # threshold of the overall score.
# Candidate with a score larger than the threshold will be kept in the table of pairs
theoretical_ratio <- 1  # theoretical intensity ratio of the pairs; default: 1 
tagMz_cutoff <- 1.0 # m/z cutoff

rt_cutoff <- 1 # retention time cutoff
ratio_cutoff <- 1 # intensity ratio cutoff

tags_mz <- c(2.0126, 4.0251, 6.0377, 8.0502) 
# theoretical m/z difference of different number of tags
# one tag: 2.0126
# two tags: 4.0251
# three tags: 6.0377
# four tags: 8.0502
evaluation <- TRUE # whether evaluate the pairs
message(Sys.time())
if (run_pairPicking) {
        message("------ Start to pick the pairs ------")
        library(xcms)
        if(run_inSourceFrag){
                files= list.files(pattern = "ISFrag.csv")
        }else{
                files <- list.files(pattern = "Raw.csv")       
        }
        
        lcmsFiles <- list.files(pattern = ".mzML")
        
        for (num in 1:length(files)) {
                
                ##### find pairs ####
                original_file <- files[num]
                
                if(run_inSourceFrag){
                        sample_name <- strsplit(files[num], "ISFrag.csv")[[1]][1]
                }else{
                        sample_name <- strsplit(files[num], "Raw.csv")[[1]][1]     
                }
                
                mess <- paste0("------ Pick pairs for sample ", sample_name, "------")
                message(mess)
                ft <- read.csv(original_file)
                # exclude features w/o retention time information  
                if(!run_peakPicking){ft <- ft[ft$rt != "N/A",]} 

                ft$mz <- as.numeric(ft$mz)
                ft$rt <- as.numeric(ft$rt)
                ft$Intensity <- as.numeric(ft$Intensity)
                ft <- ft[order(ft$mz),] # feature table is sorted by increasing precursor m/z
                
                colname <- c("ID_light", "RT_light","mz_light","Int_light", "ID_heavy",
                             "RT_heavy","mz_heavy", "Int_heavy", "num_tag","mz_score","rt_score",
                             "ratio_score" ,"score")
                record <- as.data.frame(matrix(ncol = length(colname)))
                colnames(record) <- colname
                
                
                d <- 0
                N <- nrow(ft) - 1
                if(run_inSourceFrag){
                        # if run the in-source fragment identification, add a "ISFragment" column
                        record$ISFragment <- NA
                }
                individual_record <- record
                # for each metabolic feature, search for possible heavy compounds
                for (i in 1:N) {
                        #i<-which(ft$feature_index == "S1_1871")
                        #i<-13130
                        t_rt <- as.numeric(ft$rt[i])
                        t_mz <- ft$mz[i]
                        t_int <- ft$Intensity[i]
                        ft_sub <- ft[(i+1):nrow(ft),] # candidates whose m/z is larger than the target feature
                        
                        if(run_inSourceFrag){
                                # check if this light mz is in_source fragment
                                if(ft$Num_Level2[i] != 0 | ft$Num_Level1[i] !=0){
                                        t_in_source <- paste0("light:", ft$ISF_level[i])
                                }else{
                                        t_in_source <- ""
                                }
                                
                        }
                        # first filter by the rt, intensity ratio with wide tolerance 
                        ft_sub <- ft_sub[ft_sub$rt - t_rt <= rt_cutoff * rt_diff[2],]
                        ft_sub <- ft_sub[ft_sub$rt - t_rt >= rt_cutoff * rt_diff[1],]
                        ft_sub <- ft_sub[ ft_sub$Intensity/t_int  >= ratio_cutoff * int_ratio[1],]
                        ft_sub <- ft_sub[ ft_sub$Intensity/t_int <= ratio_cutoff * int_ratio[2],]
                        # if no candidate is found after primary filtering, skip to the next feature
                        if (nrow(ft_sub) == 0) {
                                next
                        }
                        
                        # record highest score with different number of tags: one, two, three, four 
                        cname <- c("ID_heavy","RT_heavy","mz_heavy", 
                                   "Int_heavy", "num_tag","mz_score",
                                   "rt_score","ratio_score" ,"score",
                                   "pearson","cross_cor" )
                        tag_df <- as.data.frame(matrix(ncol = length(cname), nrow = length(tags_mz)))
                        colnames(tag_df) <- cname
                        
                        for (k in 1:length(tags_mz)) {
                                
                                # filter by the m/z with wide m/z tolerance 
                                ft_sub2 <- ft_sub[abs(ft_sub$mz - t_mz - tags_mz[k]) <= tagMz_cutoff * (t_mz+tags_mz[k]) * heavy_mz_tol*0.000001,]
                                
                                if (nrow(ft_sub2) == 0) {
                                        tag_df$score[k] <- -500
                                        next
                                }
                                
                                mz_score <- 1 - abs(ft_sub2$mz - t_mz - tags_mz[k])/(t_mz+tags_mz[k])/(heavy_mz_tol*0.000001)
                                rt_score <- vector()
                                for (m in 1:nrow(ft_sub2)){
                                        if(ft_sub2$rt[m] < t_rt){
                                                # retention time shift earlier
                                                rt_score_ind <- 1 - abs(ft_sub2$rt[m] - t_rt)/(abs(rt_diff[1]))
                                        }else{
                                                rt_score_ind <- 1 - abs(ft_sub2$rt[m] - t_rt)/rt_diff[2]
                                        }
                                        rt_score <- c(rt_score, rt_score_ind)
                                }
                                
                                ratio_score <- vector()
                                for (m in 1:nrow(ft_sub2)){
                                        if(ft_sub2$Intensity[m] < t_int){
                                                # intensity lower than light intensity
                                                ratio_score_ind <- 1 - abs(ft_sub2$Intensity[m]/t_int - theoretical_ratio)/(1-int_ratio[1])
                                        }else{
                                                ratio_score_ind <- 1 - abs(ft_sub2$Intensity[m]/t_int - theoretical_ratio)/(int_ratio[2]-1)
                                        }
                                        ratio_score <- c(ratio_score, ratio_score_ind)
                                }
                                
                                score <- weight_rt * rt_score + weight_mz * mz_score + weight_ratio * ratio_score
                                index <- which.max(score)
                                
                                tag_df$ID_heavy[k] <- ft_sub2$feature_index[index]
                                tag_df$RT_heavy[k] <- ft_sub2$rt[index]
                                tag_df$mz_heavy[k] <- ft_sub2$mz[index]
                                tag_df$Int_heavy[k] <- ft_sub2$Intensity[index]
                                tag_df$num_tag[k] <- k
                                tag_df$mz_score[k] <- mz_score[index]
                                tag_df$rt_score[k] <- rt_score[index]
                                tag_df$ratio_score[k] <- ratio_score[index]
                                tag_df$score[k] <- score[index]
                                
                        }
                        
                        # if no heavy compounds are found, skip to the next feature
                        if (max(tag_df$score) == -500) {
                                next
                        }
                        
                        tag_index <- which.max(tag_df$score)
                        
                        # if multiple max, choose the one with highest m/z score 
                        if (length(tag_index) >1) {
                                index <- which(tag_df$mz_score[tag_index] == max(tag_df$mz_score[tag_index]))
                                tag_index <- tag_index[index]
                                
                                # if multiple max, choose the one with highest rt score 
                                if (length(tag_index) > 1) {
                                        index <- which(tag_df$rt_score[tag_index] == max(tag_df$rt_score[tag_index]))
                                        tag_index <- tag_index[index]
                                        
                                        if (length(tag_index) > 1) {
                                                index <- which(tag_df$ratio_score[tag_index] == max(tag_df$ratio_score[tag_index]))
                                                tag_index <- tag_index[index]
                                        }
                                }
                                
                        }
                        
                        individual_record$ID_light[1] <- ft$feature_index[i]
                        individual_record$RT_light[1] <- ft$rt[i]
                        individual_record$mz_light[1] <- ft$mz[i]
                        individual_record$Int_light[1] <- round(ft$Intensity[i],0)
                        
                        individual_record$ID_heavy[1] <- tag_df$ID_heavy[tag_index]
                        individual_record$RT_heavy[1] <- tag_df$RT_heavy[tag_index]
                        individual_record$mz_heavy[1] <- tag_df$mz_heavy[tag_index]
                        individual_record$Int_heavy[1] <- round(tag_df$Int_heavy[tag_index],0)
                        individual_record$num_tag[1] <- tag_df$num_tag[tag_index] 
                        individual_record$mz_score[1] <- tag_df$mz_score[tag_index]
                        individual_record$rt_score[1] <- tag_df$rt_score[tag_index]
                        individual_record$ratio_score[1] <- tag_df$ratio_score[tag_index]
                        individual_record$score[1] <- tag_df$score[tag_index]
                        
                        if(run_inSourceFrag){
                                # check if the target is also in_source fragment;
                                h_index <- which(ft$feature_index == tag_df$ID_heavy[tag_index])
                                
                                if(ft$Num_Level2[h_index] != 0 | ft$Num_Level1[h_index] != 0){
                                        heavy_in_source <- paste0("heavy:", ft$ISF_level[h_index])
                                }else{heavy_in_source = ""}
                                
                                if(t_in_source == "" | heavy_in_source == ""){
                                        individual_record$ISFragment[1] <- paste0(t_in_source,heavy_in_source)
                                }else{individual_record$ISFragment[1] <- paste0(t_in_source," " ,heavy_in_source) }
                                
                        }
                        
                        record <- rbind(record, individual_record)
                        # update the process 
                        #if((i*100)%/%nrow(ft) > d) {
                        #        d <- d+1
                        #        message(paste0(d,"% Completed"))
                        #}
                }
                record <- record[-1,]
                
                if(run_inSourceFrag){
                        filename <- strsplit(original_file,"_ISFrag")[[1]][1]
                }else{
                        filename <- strsplit(original_file,"_Raw")[[1]][1]     
                }
                
                #hist(record$score, breaks = 50, main = filename, xlim = c(-1,1))
                #write.csv(record, paste0(filename, "_v1_pairs_all.csv"), row.names = F)
                # record <- record[record$score >= score_threshold,]
                hist(record$score, breaks = 50)
                #write.csv(record, paste0(filename, "_pairs.csv"), row.names = F)
                
                rawlcms <- xcmsRaw(lcmsFiles[num])
                #### evaluate the quality of pairs using cross-correlation #####
                if (evaluation) {
                        message("        Evaluate quality of pairs with EIC: rt window 60 s ")
                        pair_tb <- record
                        if (!run_peakPicking){
                                pair_tb$corrected_lint <- NA
                                pair_tb$corrected_hint <- NA
                        }
 
                        parent_dir <- getwd()
                        eic_dir <- paste0(parent_dir, paste0("/EIC_of_pairs_for_",sample_name))
                        dir.create(eic_dir)
                        setwd(eic_dir)
                        for(i in 1:nrow(pair_tb)) {
                                
                                # extract EIC for light, RT range: 60 seconds
                                light_mz <- pair_tb$mz_light[i]
                                
                                light_rt <- pair_tb$RT_light[i]
                                #light_mzRange <- c(light_mz-0.01, light_mz+0.01)
                                #light_rtRange <- c(light_rt*60 - 30, light_rt*60 + 30)
                                
                                light_eic_matrix <- EIC_matrix(rawlcms, mz= light_mz, rt=light_rt)
                                light_int_v <- peak_smooth(light_eic_matrix[,2])
                                #light_raweic <- rawEIC(rawlcms, mzrange= light_mzRange, rtrange = light_rtRange)
                                #light_eic_matrix <- cbind(rawlcms@scantime[light_raweic$scan],
                                #                          light_raweic$intensity)
                                
                                heavy_mz <- pair_tb$mz_heavy[i]
                                heavy_rt <- pair_tb$RT_heavy[i]
                                #heavy_mzRange <- c(heavy_mz-0.01, heavy_mz+0.01)
                                #heavy_rtRange <- light_rtRange
                                heavy_eic_matrix <- EIC_matrix(rawlcms, mz = heavy_mz, rt = light_rt )
                                heavy_int_v <- peak_smooth(heavy_eic_matrix[,2])
                                #heavy_raweic <- rawEIC(rawlcms, mzrange= heavy_mzRange, rtrange = heavy_rtRange)
                                #heavy_eic_matrix_test <- cbind(rawlcms@scantime[heavy_raweic$scan],
                                #                          heavy_raweic$intensity)
                                #plot(heavy_eic_matrix_test[,1],heavy_eic_matrix_test[,2],
                                #     type="l")
                                #plot(heavy_eic_matrix[,1], heavy_eic_matrix[,2], type="l")
                                
                                if (!run_peakPicking){
                                        pair_tb$corrected_lint[i] <- max(light_eic_matrix[,2])
                                        pair_tb$corrected_hint[i] <- max(heavy_eic_matrix[,2])
                                }
                                
                                if(max(light_eic_matrix[,2]) == 0 || max(heavy_eic_matrix[,2]) ==0 ){
                                        pair_tb$pearson[i] <- 0
                                        pair_tb$cross_cor[i] <- 0
                                        next
                                }
                                
                                pearson <- round(cor(light_int_v, heavy_int_v, method = "pearson"),2)
                                cross_cor <-  ccf(light_int_v, heavy_int_v ) # cross correlation
                                ccvalue <- round(max(cross_cor$acf),2)
                                pair_tb$pearson[i] <- pearson
                                pair_tb$cross_cor[i] <- ccvalue
                                if(ccvalue >=cc_threshold ){
                                        # plot EIC
                                        png(paste0("cc_",ccvalue," l_mz_",round(light_mz,4),"_rt_",round(light_rt,2), "_h_mz",round(heavy_mz,4),"_rt_",round(heavy_rt,2),".png"))
                                        #par(mfrow=c(2,1))
                                        #plot(light_eic_matrix[,1]/60,
                                        #     light_eic_matrix[,2], type="l",
                                        #     xlab="retetion time (min)", ylab="intensity",
                                        #     main=paste0("l_mz:",light_mz,"_rt:",round(light_rt,2)," cross_cor:", ccvalue))
                                        
                                        #plot(heavy_eic_matrix[,1]/60,
                                        #     heavy_eic_matrix[,2], type="l",
                                        #     xlab="retention time (min)", ylab="intensity",
                                        #     main=paste0("h_mz:",heavy_mz,"_rt:",round(heavy_rt,2)," pearson:",pearson))
                                        
                                        #par(mfrow=c(2,1))
                                        #plot(light_eic_matrix[,1]/60,
                                        #     light_eic_matrix[,2], type="l",
                                        #     xlim=c(light_rt-1, light_rt+1),
                                        #     ylim=c(0, max(light_eic_matrix[,2],heavy_eic_matrix[,2])),
                                        #     xlab="retetion time (min)", ylab="intensity",
                                        #     main=paste0("lmz:",light_mz,"hmz:",heavy_mz,"_rt:",round(light_rt,2)," cross_cor:", ccvalue))
                                        #lines(heavy_eic_matrix[,1]/60,heavy_eic_matrix[,2], col="red")
                                        #dev.off()
                                        
                                        #plot smoothed peak
                                        plot(light_eic_matrix[,1]/60,
                                             light_int_v, type="l",
                                             xlim=c(light_rt-1, light_rt+1),
                                             ylim=c(0, max(light_eic_matrix[,2],heavy_eic_matrix[,2])),
                                             xlab="retetion time (min)", ylab="intensity",
                                             main=paste0("lmz:",light_mz,"hmz:",heavy_mz,"_rt:",round(light_rt,2)," cross_cor:", ccvalue))
                                        lines(heavy_eic_matrix[,1]/60,
                                              heavy_int_v, col="red")
                                        dev.off()
                                }

                        
                        }
                        setwd(parent_dir)
                        #write.csv(pair_tb,paste0(filename,"_pairs_before_cc_filter.csv"), row.names = FALSE)
                        
                        if (!run_peakPicking){
                                pair_tb <- pair_tb[pair_tb$corrected_lint >= 1000,]
                                pair_tb <- pair_tb[pair_tb$corrected_hint >= 1000,]
                        }
                        
                        #pair_tb <- read.csv(paste0(filename,"_pairs_no_int_ccf_valuation.csv"))
                        pair_tb <- pair_tb[pair_tb$cross_cor>= (cc_threshold - 0.00001),]
                        
                        write.csv(pair_tb,paste0(filename,"_pairs.csv"), row.names = FALSE)
                        #hist(pair_tb$RT_heavy - pair_tb$RT_light, breaks = 20)
                        #dens <- density(pair_tb$RT_heavy - pair_tb$RT_light)
                        #plot(dens$x, dens$y, type="l")
                        #median(pair_tb$RT_heavy - pair_tb$RT_light)
                        
                        #test <- pair_tb[pair_tb$cross_cor>=0.75,]
                        #hist(test$RT_heavy - test$RT_light, breaks = 20)
                        #dens <- density(test$RT_heavy - test$RT_light)
                        #plot(dens$x, dens$y, type="l")
                        #median(test$RT_heavy - test$RT_light)
                }
                
  ### 1.3 data cleaning of in source fragments, carbon isotope, salt adduct ####
                #pair_tb <- read.csv(list.files(pattern = "_pairs.csv")[1])
                pair_tb <- cbind(0, pair_tb)
                colnames(pair_tb)[1] <- c("pair_ID")
                pair_tb$pair_ID <- paste0("P", num,"_" ,1:nrow(pair_tb))
                ############ 1) in-source-fragment identification  ##############
                if(run_inSourceFrag){
                        
                        pair_tb$isFrag <- ""
                        for (i in 1:nrow(pair_tb)){
                                
                                isfrag <- pair_tb$ISFragment[i]
                                if(isfrag == ""){next}
                                if(grepl('light',isfrag, fixed = T) & grepl("heavy", isfrag, fixed = T)){
                                        
                
                                        # in this pair: both heavy and light are in-source fragment
                                        # extract the parent ion of light
                                        light_parent <- strsplit(isfrag," ")[[1]][1]
                                        light_parent <- strsplit(light_parent, ":")[[1]][2]
                                        light_parent <- strsplit(light_parent,";")[[1]]
                                        light_parent_v <- vector()
                                        for(j in 1:length(light_parent)){
                                                light_parent_v <- c(light_parent_v,strsplit(light_parent[j], "<-")[[1]][1])
                                        }
                                        # find all parent index in the pair table
                                        light_parent_index <- vector()
                                        for(j in 1:length(light_parent_v)){
                                                index <- which(grepl(light_parent_v[j], pair_tb$ID_light,fixed = T) == TRUE)
                                                if(length(index) != 0){
                                                        light_parent_index <- c(light_parent_index, index)
                                                }
                                        }
                                        # extract the parent ion of heavy
                                        heavy_parent <- strsplit(isfrag," ")[[1]][2]
                                        heavy_parent <- strsplit(heavy_parent, ":")[[1]][2]
                                        heavy_parent <- strsplit(heavy_parent,";")[[1]]
                                        heavy_parent_v <- vector()
                                        for(j in 1:length(heavy_parent)){
                                                heavy_parent_v <- c(heavy_parent_v,strsplit(heavy_parent[j], "<-")[[1]][1])
                                        }
                                        
                                        # find all parent index in the pair table
                                        heavy_parent_index <- vector()
                                        for(j in 1:length(heavy_parent_v)){
                                                index <- which(grepl(heavy_parent_v[j], pair_tb$ID_heavy,fixed = T) == TRUE)
                                                if(length(index) != 0){
                                                        heavy_parent_index <- c(heavy_parent_index, index)
                                                }
                                        }
                                        
                                        # check if any overlap between light_parent_index and heavy_parent_index
                                        overlap_index <- intersect(heavy_parent_index,light_parent_index)
                                        
                                        if(length(overlap_index) != 0){
                                                pair_tb$isFrag[i] <- paste0("ISFrag_of_",pair_tb$pair_ID[overlap_index][1])
                                        }
                                }
                        }
                }
                
                

                #  do the isotope removal and adduct determination directly for this file
                ######## 2) carbon isotope identification #######
                remove_index <- vector()
                pair_tb$Isotope <- ""
                pair_tb$Adduct <- ""
                pair_tb$Comment <- ""
                for (i in 1:nrow(pair_tb)){
                        
                        # if this pair has been determined as other isotope, skip to next
                        if(pair_tb$Isotope[i] != ""){
                                next
                        }
                        light_rt <- pair_tb$RT_light[i]
                        light_mz <- pair_tb$mz_light[i]
                        light_int <- pair_tb$Int_light[i]
                        num_tag <- pair_tb$num_tag[i]
                        
                        #light_mzRange <- c(light_mz-0.01, light_mz+0.01)
                        #light_rtRange <- c(light_rt*60 - 30, light_rt*60 + 30)
                        light_eic_matrix <- EIC_matrix(rawlcms, mz=light_mz, rt=light_rt)
                        light_int_v1 <- peak_smooth(light_eic_matrix[,2])
                        
                        rt_index <- which(abs(pair_tb$RT_light - light_rt)<=0.1)
                        rt_index <- rt_index[rt_index!=i]
                        tag_match_index <- which(pair_tb$num_tag[rt_index] == num_tag)
                        rt_index <- rt_index[tag_match_index]
                                
                        mz_index <- which(abs(pair_tb$mz_light[rt_index] - light_mz-1.0034) <= 0.01)
                        
                        # no rt and mz match the first isotope
                        if(length(mz_index) ==0){
                                pair_tb$Isotope[i] <- "[M+0]"
                                next
                        }
                        mz_index <- which.min(abs(pair_tb$mz_light[rt_index] - light_mz-1.0034))
                        final_index1 <- rt_index[mz_index]
                        # check the intensity for the first isotope #
                        light_int1 <- pair_tb$Int_light[final_index1] 
                        pair_tb$Isotope[i] <- "[M+0]"
                        if(light_int1 < light_int){
                                iso1_mz <- pair_tb$mz_light[final_index1]
                                #iso1_mzRange <- c(iso1_mz -0.01, iso1_mz+0.01)
                                iso1_eic_matrix <- EIC_matrix(rawlcms, iso1_mz, light_rt)
                                
                                PPC <-   cor(light_int_v1, iso1_eic_matrix[,2], method = "pearson")
                                if(PPC < 0.8){
                                        next
                                }
                                pair_tb$Isotope[final_index1] <- "[M+1]"
                                if(pair_tb$Comment[final_index1] != ""){
                                        pair_tb$Comment[final_index1] <- paste0(pair_tb$Comment[final_index1],";isotope_of_",pair_tb$pair_ID[i]) 
                                }else{
                                        pair_tb$Comment[final_index1] <- paste0("isotope_of_",pair_tb$pair_ID[i])
                                }

                                remove_index <- c(remove_index, final_index1)
                                
                                # check the second isotope #
                                mz_index <- which(abs(pair_tb$mz_light[rt_index] - light_mz-2.0067) <= 0.01)
                                if(length(mz_index) ==0){next}
                                mz_index <- which.min(abs(pair_tb$mz_light[rt_index] - light_mz-2.0067))
                                final_index2 <- rt_index[mz_index]
                                
                                # check the intensity for second isotope #
                                light_int2 <- pair_tb$Int_light[final_index2] 
                                if(light_int2 < light_int1){
                                        iso2_mz <- pair_tb$mz_light[final_index2]
                                        #iso2_mzRange <- c(iso2_mz -0.01, iso2_mz+0.01)
                                        iso2_eic_matrix <- EIC_matrix(rawlcms, iso2_mz, light_rt)
                                        light_int_v2 <- peak_smooth(iso2_eic_matrix[,2])
                                        PPC <-   cor(light_int_v1,light_int_v2 , method = "pearson")
                                        if(PPC < 0.8){next}
                                        pair_tb$Isotope[final_index2] <- "[M+2]"
                                        if(pair_tb$Comment[final_index2] != ""){
                                                pair_tb$Comment[final_index2] <- paste0(pair_tb$Comment[final_index2],";isotope_of_",pair_tb$pair_ID[i])
                                        }else{
                                                pair_tb$Comment[final_index2] <- paste0("isotope_of_",pair_tb$pair_ID[i])
                                        }
          
                                        remove_index <- c(remove_index, final_index2)
                                        
                                        # check the third isotope #
                                        mz_index <- which(abs(pair_tb$mz_light[rt_index] - light_mz-3.0101) <= 0.01)
                                        if(length(mz_index) == 0){next}
                                        mz_index <- which.min(abs(pair_tb$mz_light[rt_index] - light_mz-3.0101))
                                        final_index3 <- rt_index[mz_index]
                                        
                                        # check the intensity for third isotope #
                                        light_int3 <- pair_tb$Int_light[final_index3]
                                        if(light_int3 < light_int2){
                                                iso3_mz <- pair_tb$mz_light[final_index3]
                                                #iso3_mzRange <- c(iso3_mz -0.01, iso3_mz+0.01)
                                                iso3_eic_matrix <- EIC_matrix(rawlcms, iso3_mz, light_rt)
                                                light_int_v3 <- peak_smooth(iso3_eic_matrix[,2])
                                                PPC <-   cor(light_int_v1, light_int_v3, method = "pearson")
                                                if(PPC < 0.8){next}
                                                pair_tb$Isotope[final_index3] <- "[M+3]"
                                                if(pair_tb$Comment[final_index3]!=""){
                                                        pair_tb$Comment[final_index3] <- paste0(pair_tb$Comment[final_index3],";isotope_of_",pair_tb$pair_ID[i])
                                                }else{
                                                        pair_tb$Comment[final_index3] <- paste0("isotope_of_",pair_tb$pair_ID[i])
                                                }
                                                remove_index <- c(remove_index, final_index3)
                                        }else{next}
                                        
                                }else{next}
                        }else{next}
                }
                #pair_tb <- pair_tb[-remove_index,]
                #write.csv(pair_tb,paste0(filename,"_pairs_isotope.csv"), row.names = FALSE)
                
                ##### 3) salt adduct identification ####### 
                # H+: 1.007276 
                # Na+:  22.989221
                NaAdduct <- 22.989221- 1.007276
                # K+: 38.963158
                KAdduct <- 38.963158 - 1.007276
                # NH4+: 18.033826
                NH4Adduct <- 18.033826 - 1.007276
                # -H2O: 
                H2OAdduct <- 18.010565
                
                for (i in 1:nrow(pair_tb)){
                        light_rt <- pair_tb$RT_light[i]
                        light_mz <- pair_tb$mz_light[i]
                        num_tag <- pair_tb$num_tag[i]
                        #light_mzRange <- c(light_mz-0.01, light_mz+0.01)
                        #light_rtRange <- c(light_rt*60 - 30, light_rt*60 + 30)
                        # extract EIC for light, RT range: 60 seconds
                        light_eic_matrix <- EIC_matrix(rawlcms, light_mz, light_rt)
                        light_int_v <- light_eic_matrix[,2]
                        
                        
                        rt_index <- which(abs(pair_tb$RT_light - light_rt)<=0.1)
                        rt_index <- rt_index[rt_index!=i]
                        tag_match_index <- which(pair_tb$num_tag[rt_index] == num_tag)
                        rt_index <- rt_index[tag_match_index]
                        
                        ############# Na adduct ##############
                        mz_index <- which(abs(pair_tb$mz_light[rt_index] - light_mz - NaAdduct) <= 0.01)
                        if(length(mz_index) !=0){
                                mz_index <- which.min(abs(pair_tb$mz_light[rt_index]- light_mz  - NaAdduct))
                                final_index_Na <- rt_index[mz_index]
                                mz_Na <- pair_tb$mz_light[final_index_Na]
                                Na_eic_matrix <- EIC_matrix(rawlcms, mz_Na, light_rt)
                                Na_int_v <- Na_eic_matrix[,2]
                                Na_pearson <- round(cor(light_int_v, Na_int_v, method = "pearson"),2)
                                if(Na_pearson >= 0.8) {
                                        #pair_tb$Adduct[i] <- paste0(pair_tb$Adduct[i],"[M+H]+_of_", 
                                        #                            pair_tb$ID_light[final_index_Na],"/",
                                        #                            pair_tb$ID_heavy[final_index_Na]," ")
                                        if(pair_tb$Adduct[final_index_Na] != ""){
                                                pair_tb$Adduct[final_index_Na] <- paste0(pair_tb$Adduct[final_index_Na], 
                                                                                         ";[M+Na]+_of_",pair_tb$pair_ID[i])
                                        }else{pair_tb$Adduct[final_index_Na] <- paste0("[M+Na]+_of_",pair_tb$pair_ID[i])}
                                }
                        }
                        ############ K adduct ###########
                        mz_index <- which(abs(pair_tb$mz_light[rt_index] - light_mz - KAdduct) <= 0.01)
                        if(length(mz_index) !=0){
                                mz_index <- which.min(abs(pair_tb$mz_light[rt_index]- light_mz  - KAdduct))
                                final_index_K <- rt_index[mz_index]
                                mz_K <- pair_tb$mz_light[final_index_K]
                                K_eic_matrix <- EIC_matrix(rawlcms, mz_K, light_rt)
                                K_int_v <- K_eic_matrix[,2]
                                K_pearson <- round(cor(light_int_v, K_int_v, method = "pearson"),2)
                                if(K_pearson >= 0.8) {
                                        #pair_tb$Adduct[i] <- paste0(pair_tb$Adduct[i],"[M+H]+_of_", 
                                        #                            pair_tb$ID_light[final_index_K],"/",
                                        #                            pair_tb$ID_heavy[final_index_K]," ")
                                        if(pair_tb$Adduct[final_index_K] != ""){
                                                pair_tb$Adduct[final_index_K] <- paste0(pair_tb$Adduct[final_index_K], 
                                                                                         ";[M+K]+_of_",pair_tb$pair_ID[i])
                                        }else{pair_tb$Adduct[final_index_K] <- paste0("[M+K]+_of_",pair_tb$pair_ID[i])}
                                }
                        }
                        
                        ############ NH4 adduct ##########
                        mz_index <- which(abs(pair_tb$mz_light[rt_index] - light_mz - NH4Adduct) <= 0.01)
                        if(length(mz_index) !=0){
                                mz_index <- which.min(abs(pair_tb$mz_light[rt_index] - light_mz- NH4Adduct))
                                final_index_NH4 <- rt_index[mz_index]
                                mz_NH4 <- pair_tb$mz_light[final_index_NH4]
                                #mzRange_NH4 <- c(mz_NH4-0.01, mz_NH4+0.01)
                                # extract EIC for adduct, RT range: 40 seconds
                                NH4_eic_matrix <- EIC_matrix(rawlcms, mz_NH4 ,light_rt)
                                NH4_int_v <- NH4_eic_matrix[,2]
                                NH4_pearson <- round(cor(light_int_v, NH4_int_v, method = "pearson"),2)
                                if(NH4_pearson >= 0.8) {
                                        #pair_tb$Adduct[i] <- paste0(pair_tb$Adduct[i],"[M+H]+_of_", 
                                        #                            pair_tb$ID_light[final_index_NH4]," ")
                                        if(pair_tb$Adduct[final_index_NH4] != ""){
                                                pair_tb$Adduct[final_index_NH4] <- paste0(pair_tb$Adduct[final_index_NH4], 
                                                                                         ";[M+NH4]+_of_",pair_tb$pair_ID[i])
                                        }else{pair_tb$Adduct[final_index_NH4] <- paste0("[M+NH4]+_of_",pair_tb$pair_ID[i])}
                                }
                        }
                        ############ H2O adduct ##########
                        mz_index <- which(abs(pair_tb$mz_light[rt_index] + H2OAdduct -light_mz) <= 0.01)
                        
                        if(length(mz_index) !=0){
                                mz_index <- which.min(abs(pair_tb$mz_light[rt_index] + H2OAdduct - light_mz))
                                final_index_H2O <- rt_index[mz_index]
                                mz_H2O <- pair_tb$mz_light[final_index_H2O]
                                #mzRange_H2O <- c(mz_H2O-0.01, mz_H2O+0.01)
                                # extract EIC for adduct, RT range: 40 seconds
                                H2O_eic_matrix <- EIC_matrix(rawlcms, mz_H2O,light_rt)
                                H2O_int_v <- H2O_eic_matrix[,2]
                                H2O_pearson <- round(cor(light_int_v, H2O_int_v, method = "pearson"),2)
                                if(H2O_pearson >= 0.8) {
                                        #pair_tb$Adduct[i] <- paste0(pair_tb$Adduct[i],"[M+H]+_of_", 
                                        #                                 pair_tb$ID_light[final_index_H2O]," ")
                                        if(pair_tb$Adduct[final_index_H2O] != ""){
                                                pair_tb$Adduct[final_index_H2O] <- paste0(pair_tb$Adduct[final_index_H2O], 
                                                                                          ";[M+H-H2O]+_of_",pair_tb$pair_ID[i])
                                        }else{pair_tb$Adduct[final_index_H2O] <- paste0("[M+H-H2O]+_of_",pair_tb$pair_ID[i])}
                                }
                        }
                        
                }
                original_isfrag <- which(colnames(pair_tb) == "ISFragment")
                if(length(original_isfrag) != 0){
                        pair_tb <- pair_tb[,-original_isfrag]
                }
                write.csv(pair_tb,paste0(filename,"_pairs_isfrag_isotope_adduct.csv"), row.names = FALSE)
                #write.csv(pair_tb,paste0(filename,".csv"), row.names = FALSE)
        }
}

#------ 2rd module: alignment ---------------
align_weight_mz <- 0.5 # weight of precursor m/z when calculating the alignment score
align_weight_rt <- 0.5 # weight of retention time when calculating the alignment score

# filter by the number of tags, m/z (ppm), retention time and intensity ratio
# if multiple candidates, calculate the overall score
int_record <- "light" # determine how to record the intensity in the alignment table
                      # light: light compound 
                      # heavy: heavy compound
                      # average: average of light and heavy compound

if (run_alignment) {
        files <- list.files(pattern = "isfrag_isotope_adduct.csv")
        #files <- list.files(pattern = "Raw.csv")
        N <- length(files) # record the number of files
        
        # load each csv file in file_1, file_2 ... file_N
        sum_pairs <- 0
        in_source <- 0
        for (i in 1:N) {
                
                file <- read.csv(files[i])
                
                file <- file[file$Isotope == "[M+0]",]
                if(remove_isfrag_pairs){
                        # only keep none isfrag
                        # if no isfrag found: the table could be NA for all instead of ''
                        remove_index <- which(grepl('ISFrag_of_',file$isFrag) == TRUE)
                        if(length(remove_index) != 0) {
                                file <- file[-remove_index,]  # file[-empty,] will remove everthing
                        }
                        in_source <- in_source + length(remove_index)
                }
                assign(paste0("ft_",i), file)
                message(paste(files[i], nrow(file)))
                sum_pairs <- sum_pairs + nrow(file)
                
        }
        message(paste0("There are ", sum_pairs, " pairs for ", N, " samples"))
        
        # generate a table to record alignment results
        if(!remove_isfrag_pairs){
                label <- c("ID_light","RT_light","mz_light",
                           "ID_heavy","mz_heavy","num_tag",
                           "Adduct","ISFrag","pairID")
        }else{
                label <- c("ID_light","RT_light","mz_light",
                           "ID_heavy","mz_heavy","num_tag",
                           "Adduct","pairID")
        }

        for (i in 1:length(files)) {
                label <- c(label, strsplit(files[i], "_pairs")[[1]][1])
        }
        record <- as.data.frame(matrix(ncol = length(label)))
        colnames(record) <- label
        individual <- record # individual is used to generate new row
        
        # choose the source column of the intensity
        if (int_record == "light") {intCol <-  which(colnames(ft_1) == "Int_light")}
        if (int_record == "heavy") {intCol <-  which(colnames(ft_1) == "Int_heavy")}
        if (int_record == "average") {intCol <- c(which(colnames(ft_1) == "Int_light"),
                                                  which(colnames(ft_1) == "Int_heavy"))}
        
        # fill the table with the first file, leave int_2, int_3 ... int_n as NA
        for (i in 1:nrow(ft_1)) {
                
                individual$ID_light <- ft_1$ID_light[i]
                individual$RT_light <- ft_1$RT_light[i]
                individual$mz_light <- ft_1$mz_light[i]
                individual$ID_heavy <- ft_1$ID_heavy[i]
                individual$mz_heavy <- ft_1$mz_heavy[i]
                individual$num_tag <- ft_1$num_tag[i]
                if(!remove_isfrag_pairs){
                        individual$ISFrag <- ft_1$isFrag[i]
                }
                individual$Adduct <- ft_1$Adduct[i]
                individual$pairID <- ft_1$pair_ID[i]
                
                if ( length(intCol) == 1) {
                        individual[,length(label)-N+1] <- round(ft_1[i, intCol],0)
                }else{
                        individual[,length(label)-N+1] <- round((ft_1[i, intCol[1]] + ft_1[i, intCol[2]])/2,0)
                }
                
                record <- rbind(record, individual)
        }
        record <- record[-1,]
        
        if(N==1){
                write.csv(record, paste0("alignment_for_",length(files),"_files.csv"), row.names=F)
                break
        }
        
        # loop all remaining files
        for (i in 2:length(files)) {
                
                ft <- get(paste0("ft_", i))
                # loop each metabolic feature in this i-th feature table based on three factors:
                ## RT, precursor m/z, number of tags
                for (j in 1:nrow(ft)) {
                        
                        t_pairID <- ft$pair_ID[j]
                        t_ID <- ft$ID_light[j]
                        t_RT <- ft$RT_light[j]
                        t_preMz <- ft$mz_light[j]
                        t_tag <- ft$num_tag[j]
                        heavy_ID <- ft$ID_heavy[j]
                        heavy_mz <- ft$mz_heavy[j]
                        adduct_form <- ft$Adduct[j]
                        isfrag <- ft$isFrag[j]
                        
                        if ( length(intCol) == 1) {
                                t_Int <- round(ft[j, intCol],0)
                        }else{
                                t_Int <- round((ft[j, intCol[1]] + ft[j, intCol[2]])/2,0)
                        }
                        
                        tag_index <- which(record$num_tag == t_tag)
                        
                        if (length(tag_index) != 0) {
                                
                                rt_index <- which((abs(record[tag_index,]$RT_light - t_RT)) <= align_rt_tol)
                                
                                if(length(rt_index) != 0) {
                                        mz_index <- which((abs(record[tag_index[rt_index],]$mz_light - t_preMz)) <= (align_mz_tol * t_preMz * 0.000001))
                                        
                                        if (length(mz_index) != 0) {
                                                
                                                # if the lose filter has multiple candidates
                                                if(length(mz_index) > 1) {
                                                        mzScore <- align_weight_mz*(1 - abs(record[tag_index[rt_index[mz_index]],]$mz_light - t_preMz)/(align_mz_tol * t_preMz * 0.000001))
                                                        rtScore <- align_weight_rt*(1 - abs(record[tag_index[rt_index[mz_index]],]$RT_light - t_RT)/align_rt_tol)
                                                        Score <- mzScore + rtScore
                                                        score_index <- which(Score  == max(Score))
                                                        mz_index <- mz_index[score_index]
                                                }
                                               
                                                # merge the this pair with the pair existing in the pair table
                                                final_index <- tag_index[rt_index[mz_index]]
                                                record[final_index,length(label)-N+i]  <- round(t_Int,0)
                                                record$pairID[final_index] <- paste0(record$pairID[final_index],";",t_pairID) 
                                               
                                                next
                                        }
                                }
                        }
                        
                        # no match based on RT, pre_mz, or tag, 
                        # directly add this feature to the end of the record table, including all the information
                        
                        individual <- as.data.frame(matrix(ncol = ncol(record)))
                        colnames(individual) <- colnames(record)
                        
                        individual$ID_light <- t_ID
                        individual$RT_light <- t_RT
                        individual$mz_light <- t_preMz
                        individual$num_tag <- t_tag
                        individual$ID_heavy <- heavy_ID
                        individual$mz_heavy <- heavy_mz
                        individual$Adduct <- adduct_form
                        if(!remove_isfrag_pairs){
                                individual$ISFrag <- isfrag
                        }
                        individual$pairID <- t_pairID
                        individual[,length(label)-N+i] <- round(t_Int,0)
                        
                        record <- rbind(record, individual)
                }
                
        }
        record <- cbind(0, record)
        colnames(record)[1] <- "AlignmentID"
        record$AlignmentID <- 1:nrow(record)
        # post-curation between the pairs #
        ID_v <- strsplit(record$pairID,";")
        #1.  for in-source-fragment
        if(!remove_isfrag_pairs){
                for(j in 1:nrow(record)){
                        if(record$ISFrag[j] =="" | is.na(record$ISFrag[j])){next}
                        isfrag <- strsplit(record$ISFrag[j],"_of_")[[1]]
                        
                        for(h in 1:length(ID_v)){
                                if(length(which(ID_v[[h]] == isfrag[2])) !=0){
                                        index <- h
                                        record$ISFrag[j] <- paste0("ISFrag_of_",index)
                                        break
                                }
                        }
                }
        }
        #2. adduct form
        for(j in 1:nrow(record)){
                
                if(record$Adduct[j] =="" | is.na(record$Adduct[j])){next}
                adduct_v <- strsplit(record$Adduct[j],";")[[1]]
                adduct_form <- vector()
                for(k in 1:length(adduct_v)){
                        
                        adduct <- strsplit(adduct_v[k],  "_of_")[[1]]
                        
                        for(h in 1:length(ID_v)){
                                if(length(which(ID_v[[h]] == adduct[2])) !=0){
                                        index <- h
                                        # found, then record, and skip to next adduct information 
                                        adduct_form <- c(adduct_form, paste0(adduct[1],"_of_",index))
                                        break
                                }
                        }
                        
                }
                record$Adduct[j] <- paste(adduct_form, collapse = ";")
        }
        
        # sort the table in the order of increasing pre_mz
        record <- record[order(record$AlignmentID),]
        write.csv(record, paste0("alignment_for_",length(files),"_files.csv"), row.names=F)
}

#------- 3rd module: missing value retrieval (gap-filling) ---------------
# corThreshold <- 0.8 # peak-peak correlation  >= 0.8, this pair's peak will be record in the gap
# gap_ratio_tol <- 0.4 # ratio tolerance of intensity ratio of the light and heavy compound 

if (run_gapFilling) {
        
        library(xcms)
        # fill the gap from one raw LC-MS to one raw LC-MS file
        lcmsFiles <- list.files(pattern = ".mzML")
        
        N <- length(lcmsFiles)
        
        # ------------- fill the missing value in the alignment result 
        dir.create("EIC_of_missingPair") # create a new folder to save the EIC
        cur_dire <- getwd()
        EIC_dire <- paste0(cur_dire,"/EIC_of_missingPair" ) 
        
        align <- read.csv(paste0("alignment_for_",length(lcmsFiles),"_files.csv"))
        align$Fill <- 0
        for(i in 1:nrow(align)){
                int_v <- vector()
                j<-1
                while(j <= N){
                        int_v <- c(int_v, align[i,ncol(align)- N -1 + j] )
                        j <- j+1
                }
                align$Fill[i] <- length(which(int_v != "NA") == TRUE)
        }
        
        align$recovered <- 0 # number of recovered after gap filling
        
        #align <- record
        
        corRecordName <- c("RT","light_mz","heavy_mz","light_int","heavy_int","correlation")
        corRecord <- as.data.frame(matrix(ncol=length(corRecordName)))
        colnames(corRecord) <- corRecordName
        
        
        corRecord_indi <- corRecord
        
        ## read the MS1 information in each mzML file, 
        #  corresponding to column 9 to 8+N is the intensity of pairs
        for (i in 1:N) {
                
                # because different file may have different RT range, 
                # we can not use xcmsRaw and EICRaw function to get the EIC (some error message from the function)
                lcms <- readMSData(lcmsFiles[i], msLevel. = 1, centroided. = TRUE)
                
                # record the retention time
                rt_v <- vector()
                for (scan in 1:length(lcms)) {
                        rt_v <- c(rt_v, lcms[[scan]]@rt)
                }
                #assign(paste0("lcms",i), lcmsRaw )
                #assign(paste0("rt_v",i), rt_v)
                num_col <- ncol(align) - N - 2 + i 
                
                # check which pair is missing, row by row
                miss_index <- which(is.na(align[,num_col]))
                for (j in 1:length(miss_index)) {
                        
                        k<-miss_index[j]
                        t_rt <- align$RT_light[k]*60
                        t_mz_light <- align$mz_light[k]
                        t_mz_heavy <- align$mz_heavy[k]
                        rtRange <- c(t_rt - gap_rt_tol*60, t_rt + gap_rt_tol*60)
                        
                        #/ get the EIC of light compound and heavy compound, and determine to record or not /# 
                        # rt_light +/- gap_rt_tol*60 (12s)
                        # 1) intensity > 1000; intensity ratio in (0.4,14); correlation > threshold (0.7)
                        
                        ## if the retention time is out of the scan range, skip to next file
                        if( rtRange[1] >= min(rt_v) && rtRange[2] <= max(rt_v) ) { 
                                
                                # first filter the MS1 spectra by retention time
                                scanIndex_rt <- which(abs(rt_v - t_rt) <= (gap_rt_tol*60))
                                if (length(scanIndex_rt) == 0) {next}
                                
                                # find the EIC for light and heavy mz
                                int_light_v <- vector()
                                int_heavy_v <- vector()
                                for ( j in scanIndex_rt) {
                                        
                                        # light compound
                                        mz_index_light <- which( abs(lcms[[j]]@mz - t_mz_light) <= (t_mz_light * gap_mz_tol*0.000001 +0.000000001))
                                        if (length(mz_index_light) == 0) {
                                                int_light_v <- c(int_light_v, 0)
                                        }
                                        # if multiple candidate, choose m/z closest one
                                        if (length(mz_index_light) >= 1) {
                                                int_light <- lcms[[j]]@intensity[mz_index_light]
                                                int_light_v <- c(int_light_v, max(int_light)) 
                                        }
                                        
                                        # heavy compound
                                        mz_index_heavy <- which( abs(lcms[[j]]@mz - t_mz_heavy) <= (t_mz_heavy * gap_mz_tol*0.000001 +0.000000001))
                                        if (length(mz_index_heavy) == 0) {
                                                int_heavy_v <- c(int_heavy_v, 0)
                                        }
                                        # if multiple candidate, choose m/z closest one
                                        if (length(mz_index_heavy) >= 1) {
                                                int_heavy <- lcms[[j]]@intensity[mz_index_heavy]
                                                int_heavy_v <- c(int_heavy_v, max(int_heavy)) 
                                        }
                                        
                                }
                                
                                # ensemble light_EIC, heavy_EIC
                                intMax_light <- max(int_light_v)
                                plot_matrix_light <- cbind(rt_v[scanIndex_rt]/60, peak_smooth(int_light_v))
                                
                                intMax_heavy <- max(int_heavy_v)
                                plot_matrix_heavy <- cbind(rt_v[scanIndex_rt]/60,peak_smooth(int_heavy_v))
                                
                                ratio_int<- intMax_heavy/intMax_light
                                
                               
                                # determine whether to record this pair: int_threshold + ratio + peak correlation
                                if (intMax_light >= int_threshold && intMax_heavy >= int_threshold) {
                                              # if int ratio out of the acceptable range, exclude 
                                              if(ratio_int > int_ratio[2] | ratio_int < int_ratio[1]){next}
                                        
                                                ccvalue <- ccf(plot_matrix_light[,2], plot_matrix_heavy[,2])
                                                cross_cor <- round(max(ccvalue$acf),2)
                                                
                                                # correlation lower than threshold, skip to next missing value
                                                if(cross_cor < cc_threshold) {next}
                                                
                                                ##  draw the EIC ##
                                                #fileName <- strsplit(lcmsFiles[i],"\\.")[[1]][1]
                                                setwd(EIC_dire)
                                                filename <- paste0("pair",align$AlignmentID[k],
                                                                   "_file",i,
                                                                   "_cor",round(cross_cor,2),
                                                                   "_ratio", round(ratio_int,2))
                                                png(paste0(filename,".png"))
                                                par(mfrow=c(1,1))
                                                
                                                plot(plot_matrix_light, type = "l",
                                                     xlab = "retention time (min)", ylab = "intensity",
                                                     xlim = c(t_rt/60 -1, t_rt/60 +1),
                                                     ylim=c(0, max(plot_matrix_heavy[,2], plot_matrix_light[,2])),
                                                     main = filename)
                                                lines(plot_matrix_heavy, type="l",col="red")
                                                dev.off()
                                                setwd(cur_dire)
                                                
                                                corRecord_indi$RT <- round(t_rt,2)
                                                corRecord_indi$light_mz <- t_mz_light
                                                corRecord_indi$heavy_mz <- t_mz_heavy
                                                corRecord_indi$light_int <- intMax_light
                                                corRecord_indi$heavy_int <- intMax_heavy
                                                corRecord_indi$correlation <- cross_cor
                                                
                                                corRecord <- rbind(corRecord, corRecord_indi)
                                                
                                                # record the recovered pair when the peak correlation > 0.8
                                                
                                                if (int_record == "light"){
                                                                align[k, num_col]  <- round(intMax_light,0) + 0.1
                                                }
                                                        
                                                if (int_record == "heavy"){
                                                                align[k, num_col]  <- round(intMax_heavy,0) + 0.1
                                                }
                                                        
                                                if (int_record == "average"){
                                                        align[k,num_col]  <- round((intMax_light + intMax_heavy)/2,0) + 0.1
                                                }
                                                align$Fill[k] <- align$Fill[k]+1
                                                align$recovered[k] <- align$recovered[k] +1

                                }
                        }
                }
                
        }
        setwd(EIC_dire)
        write.csv(corRecord, "peakCorelation_of_missed_pairs.csv", row.names = F)
        setwd(cur_dire)
        par(mfrow = c(1,1))
        
        highCor <- corRecord[corRecord$correlation >= cc_threshold,]
        hist(corRecord$correlation, breaks = 100, 
             xlim = c(0,1),
             xlab="correaltion",
             main = paste0("missing pairs:",nrow(corRecord), "_high_cor_num:",nrow(highCor)))
        write.csv(align, paste0("alignment_after_gap_filling.csv"), row.names = F)
}
sum(align$recovered)

#------ 4th module: Annotation -------------------------------------
if(run_alignment){
        filename_to_annotate <- list.files(pattern = "alignment_after_gap_filling.csv")
}else{
        filename_to_annotate <- list.files(pattern = "isotope_adduct.csv")
}
if (run_annotation) {
        lightLabel <- 11.9995 + 2*1.0073
        heavyLabel <- 11.9995 + 2*2.0136
        H <- 1.007276
        
        #filename_to_annotate <- "alignment_after_gap_filling.csv"
        #filename_to_annotate <- "alignment_for_32_files.csv"
        #filename_to_annotate <- "10_features.csv"
        # top_N <- 1 # determine how many candidates to record. 
        # 1: top 1; 
        # 2: all candidate

        gf_ft <- read.csv(filename_to_annotate)
        gf_ft$unlab_light_mz <- NA
        gf_ft$unlab_heavy_mz <- NA
        gf_ft$putative_annotation <- NA
        gf_ft$candidate_exact_mass <- NA
        gf_ft$candidate_formula <- NA
        gf_ft$mass_error <- NA
        
        ##### load the AMINES library ####
        #db <- read.csv("C:/Users/User/Desktop/2022-07-05-MSPairFinder/MSPairFinder_database/all_four/25138_primary_13176_secondary_38314_total.csv") 
        db <- read.csv("AMINES_library.csv")
        db <- db[!is.na(db$Name),]

        Sys.time()
        for (i in 1:nrow(gf_ft)) {
                
                gf_ft$unlab_light_mz[i] <- gf_ft$mz_light[i] - gf_ft$num_tag[i] * lightLabel - H
                gf_ft$unlab_heavy_mz[i] <- gf_ft$mz_heavy[i] - gf_ft$num_tag[i] * heavyLabel - H
                t_mass <- ( gf_ft$unlab_light_mz[i] + gf_ft$unlab_heavy_mz[i])/2
                
                # filter based on the neutral mass tolerance
                sub_db <- db[abs(db$Monoisotopic_Neutral_Mass - t_mass) <=  t_mass * anno_mz_tol * 0.000001,]
                
                sub_db <- sub_db[!is.na(sub_db$Name),]
                 if (nrow(sub_db) != 0) {
                        sub_db$mass_error <- round(abs(sub_db$Monoisotopic_Neutral_Mass - t_mass)/t_mass * 1000000,0)
                        sub_db <- sub_db[order(sub_db$mass_error,decreasing = F),]
                        
                        gf_ft$putative_annotation[i] <- paste0(sub_db$Name,collapse = ";")
                        gf_ft$candidate_exact_mass[i] <- paste0(round(sub_db$Monoisotopic_Neutral_Mass,4), collapse = ";")
                        gf_ft$candidate_formula[i] <- paste0(sub_db$Formula, collapse = ";")
                        gf_ft$mass_error[i] <- paste0(sub_db$mass_error, collapse = ";")
                }
        }
        Sys.time()
        annotated_feature  <- gf_ft[gf_ft$putative_annotation != "",]

        num_annotated <- nrow(annotated_feature)
        num_unique_annotated <- length(unique(annotated_feature$putative_annotation))
        #write.csv(gf_ft, paste0(nrow(annotated_feature),"_putative_annotation_for_",nrow(gf_ft),"_features","_mass_accuracy_", anno_mz_tol,"ppm.csv"), row.names = F)
        write.csv(gf_ft, "putative_annotation.csv", row.names = F)
}

