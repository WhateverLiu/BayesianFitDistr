library(plotly)


param = optRst$param
d = solve_d(param, eps= 1e-8, maxit = 100)
abcd = as.data.frame(t(rbind(param[1:3, , drop = F], d)))
abcd$MDR = targetMDRs
colnames(abcd) = c('a', 'b', 'c', 'd', 'MDR')


brks_transformed <- pretty(log(abcd$d), n = 5)
brks_untransformed <- sprintf("%.1f", exp(brks_transformed))
fig <- plot_ly(data = abcd, x = ~a, y = ~b, z = ~c, 
               type = 'scatter3d', mode = 'markers', 
               marker = list(
                 size = ~MDR*100,
                 color = ~log(d),
                 line = list(width = 0),
                 colorbar = list( 
                   title = "d",
                   tickmode = "array",
                   ticktext = brks_untransformed,
                   tickvals = brks_transformed)),
               text = ~paste('a:', a, '<br>b:', b, '<br>c:', c, '<br>d:',  d, '<br>MDR:', MDR),
               hoverinfo = 'text') %>%
  layout(title = list(text = paste0("<br>Zero-inflated limited TRB parameter estimates"),
                      font = list(size = 30)),
         scene = list(  
           xaxis = list(title = list(text = "a", font = list(size = 20)), tickfont = list(size = 15)),
           yaxis = list(title = list(text = "b", font = list(size = 20)), tickfont = list(size = 15)), 
           zaxis = list(title = list(text = "c", font = list(size = 20)), tickfont = list(size = 15)),
           aspectratio = list(x = 1, y = 1, z = 1))
         )









