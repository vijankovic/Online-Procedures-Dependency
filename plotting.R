################
### Plotting ###
################

library(ggplot2)
library(ggpubr)

# convert format for 'plotting' #
convert_results <- function(data) {
  fwer <- data$fwer_mat
  pwr  <- data$pwr_mat
  # remove first column
  fwer <- fwer[,-1]
  pwr  <- pwr[,-1]
  df   <- as.data.frame(list(as.factor(rep(rownames(fwer),each = ncol(fwer))),as.numeric(rep(colnames(fwer),times = nrow(fwer))),as.numeric(as.vector(t(fwer))),as.numeric(as.vector(t(pwr)))))
  colnames(df) <- c("Method","pi1","FWER","Power")
  return(df)
}


# generate plot #
plotting <- function(data) {
    p1 <- ggplot(data, aes(pi1, FWER, group = Method, color = Method)) +
      geom_line(aes(linetype = Method)) +
      geom_point() +
      geom_hline(yintercept=0.05, color = "black") +
      xlab(expression(pi[1])) +
      ylab("FWER") +
      scale_x_continuous(breaks = seq(0.1,0.9,0.1))
    p2 <- ggplot(data, aes(pi1, Power, group = Method, color = Method)) +
      geom_line(aes(linetype = Method)) +
      geom_point() +
      xlab(expression(pi[1])) +
      ylab("Power") +
      scale_x_continuous(breaks = seq(0.1,0.9,0.1))
    ggarrange(p1,p2,ncol = 1,nrow = 2,common.legend = TRUE,legend = "bottom")
}

# save plot #
save_plot <- function(name, plot) {
  ggsave(plot = plot, filename = paste(name,"jpg",sep = "."),device = "jpg",height = 1500,width = 2000,units = "px")
}
