

merged_data$class<-'class'
merged_data$class[merged_data$topos_ID%in%
                    c(6218,6260,6428,6434,6260, 6809,
                      7235,7427,7430,7616,7853,
                      8009,8039,8186,9212,9611,
                      9620,9788,9968)] <- "NFBL"
merged_data$class[merged_data$topos_ID%in%
                    c(7820,7841, 8213, 8228,9986)] <- "IFFL"
merged_data$class[merged_data$topos_ID%in%
                    c(5837,6431,8180,8798,10580)]<-"NFBL+IFFL"
merged_data$class[merged_data$topos_ID%in%
                    c(6071,6665,7406)]<-"NFBL+DNFBL"
plot_freq_match(merged_data,save=F)
plot_binary_topology(topos[14416,])
