#############################################################################
#
# PRACTICA 1
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
## ENTREGA EL 22 OCTUBRE 23:59
## Se requiere la entrega de un Jupyter Notebook con los códigos más las imágenes y las respuestas a las preguntas
## Adjuntar en la entrega el PDF final y el archivo con los genes
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/datos-trabajoR(2)", followlocation= TRUE)
data = read.table(file = "datos-trabajoR(2).txt", head = T)


# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data) # esto mide las dimensiones
head(data) #vemos las primeras linesas de la tabla
tail(data) # vemos las últimas líneas de la tabla

# Hacemos un primer histograma para explorar los datos
hist(data) 
hist(data,col = "gray", main ="GSE5583 - Histogram") # con esto hemos cambiado el titulo del histograma y también el color

# Transformamos los datos con un logaritmo 
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
#Los datos en nuestro histograman cambian en base al logaritmo. Para sivualizar mejor nuestros datos
data_log=log2(data)
hist(data_log)


# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
# ¿Qué es un boxplot?
#Estos parámetros sirven para cambiar los colores del boxplot. Utilizamos un vector para poder acceder a los elementos que queremos.En este caso los 3 primeros son WT y los 3 siguientes son KO.
#Le ponemos título con main. 
#las=2 nos sirve para poner los títulos de los ejes en vertical
#Un boxplot es una gráfica que nos muestra una serie de datos en base a sus cuartiles. En él vemos los cuartiles, la mediana y valores anómalos

boxplot(data_log, col=c("turquoise", "turquoise", "turquoise", "pink", "pink", "pink"),
        main="GSE5583 - boxplot", las=2)

boxplot(data_log)

# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación
# de los valores de expresión. ¿Es correcta la separación?
# Sí, es correcta
hc =hclust(as.dist(1-cor(data_log)))
plot(hc, main="GSE5583 - Hierarchical Clustering")


#######################################
# Análisis de Expresión Diferencial 
#######################################

# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
# Hemos generado de los datos las columnas separadas de ambos grupos WT y KO. 
# Con head hemos generado una matrix de los datos llamados wt.

wt<-data[,1:3]
ko<-data[,4:6]
class(wt)
head(wt)

# Calcula las medias de las muestras para cada condición. Usa apply
# 1 indica hacer la media de las filas, 2 indica hacer la media de las columnas

wt.mean =apply(wt, 1, mean)
ko.mean=apply(ko, 1, mean)
head(wt.mean)
head(ko.mean)


# ¿Cuál es la media más alta?
# La media más alta está en KO

max(wt.mean)
max(ko.mean)
 
# Ahora hacemos un scatter plot (gráfico de dispersión)

plot(ko.mean ~ wt.mean)
plot(ko.mean ~ wt.mean, 
     col="hot pink")
# Cambiamos el nombre de los ejes y le ponemos título
plot(ko.mean ~ wt.mean, 
  col="hotpink", 
  xlab ="WT", ylab = "KO", 
  main = "GSE5583 - Scatter Plot")

# Añadir una línea diagonal con abline


# ¿Eres capaz de añadirle un grid?


# Calculamos la diferencia entre las medias de las condiciones


# Hacemos un histograma de las diferencias de medias


# Calculamos la significancia estadística con un t-test.
# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué?
# ¿Cuántas valores tiene cada muestra?


# Ahora comprobamos que hemos hecho TODOS los cálculos


# Hacemos un histograma de los p-values.
# ¿Qué pasa si le ponemos con una transformación de -log10?


# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística


# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?


# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)


# Ahora el filtro de p-value


# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios?


# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo


# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés?


# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap?
# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors


# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)

# Hacemos nuestro heatmap


# Lo guardamos en un archivo PDF


# Guardamos los genes diferencialmente expresados y filtrados en un fichero

