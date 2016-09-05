export_plots=function(x,n=20){

  if (requireNamespace("grDevices", quietly = TRUE)) {
    pdf(file="Plots.pdf")
    beta_plot(x,n)
    betaprime_plot(x,n)
    betasecond_plot(x,n)
    curvature_plot(x,n)
    radius_plot(x,n)
    dev.off()
  }  else {
    print("You need to install and load the grDevices package to perform this function")
  }


}
