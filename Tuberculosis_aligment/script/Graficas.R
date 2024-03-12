library(ggplot2)
graph_volharina <- ggplot(resumenharina, aes(x=harina, y=vol)) + ## la grafica usara a los 2 niveles de harina como eje x, y usara a los del volumen como Y
  labs(title= "Volumen de Co2 producido harinas")+ ## Anotaciones que queremos hacer en nuestra grafica, titulo y nombre de ambos ejes
  geom_bar(position=position_dodge(), stat="identity", fill = c("green", "purple")) + ## Generamos nuestra grafica de barras, dodge evita que las barras se sobrepngan y fill indica el color de las barras
  geom_errorbar(aes(ymin=vol-ci, ymax=vol+ci), width=.2,position=position_dodge(.9))+ ## Bigotes de la grafica, tomaras el valor del volumen y le restaras el ci para generar el valor minimo, lo mismo para el maximo, solo que se suma el ci
  xlab("Marcas de harina")+
  ylab("Volumen de di?xido de carbono")+
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"))
graph_volharina