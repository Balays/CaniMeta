

ps_files <- list.files('EMU/ps/unfiltered', '*.PS.rds', full.names = T)

## merged phyloseq

ps_all <- read_rds(ps_files[2])

## import unfiltered EMU output phyloseqs

ps_files <- ps_files[-2]

merge_PS.rds_files <- function(ps_files) {
  
  ps_m <- read_rds(ps_files[1])
  
  for(i in ps_files[-1]) {
    
    ps   <- read_rds(ps_files[i])
    
    ps_m <- merge_phyloseq(ps, ps_m)  
    
  }
  
  return(ps_m)
}


  