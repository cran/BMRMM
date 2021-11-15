########################################
### Treating ISI as a Discrete State ###
########################################

# This function is called when 'isi_type' is 'Discrete' 
# We divide each ISI by 'isi_unit' and treat each ISI as a block of a new state
# For example, if ISI=20 and isi_unit=5, then we will add 4 consecutive "ISI" states  

#############
### Input ###
#############

# data <- original data set
# isi_unit <- specifies the unit of each block of ISI

##############
### Output ###
##############

# aug_data <- augmented data set by adding ISI as a new state

add_isi_as_state <- function(data,isi_unit) {
  aug_data <- NULL
  isi_ind <- ncol(data)
  cur_state_ind <- ncol(data)-1
  prev_state_ind <- ncol(data)-2 
  num_states <- length(unique(c(data[,prev_state_ind],data[,cur_state_ind])))
  new_state <- num_states+1
  for(row in 1:nrow(data)) {
    num_block <- floor(data[row,isi_ind]/isi_unit)
    if(row != nrow(data) && data[row,1]==data[row+1,1] && num_block>0) {
      attr <- as.numeric(data[row,1:(prev_state_ind-1)])
      prev_state <- data[row,prev_state_ind]
      cur_state <- data[row,cur_state_ind]
      aug_data <- rbind(aug_data,c(attr,prev_state,new_state))
      if(num_block>1) {
        aug_data <- rbind(aug_data,c(attr,new_state,new_state))
      }
      aug_data <- rbind(aug_data,c(attr,new_state,cur_state))
    } else {
      aug_data <- rbind(aug_data,as.numeric(data[row,1:cur_state_ind]))
    }
  }
  aug_data <- data.frame(aug_data)
  colnames(aug_data) <- colnames(data)[1:cur_state_ind]
  return(aug_data)
}