# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/datos-trabajoR.txt", followlocation= TRUE)
data = read.table(file = "datos-trabajoR.txt", head = T)
getwd #Antes de cargar los datos comprobamos el directorio

data #Comprobamos que cargamos los dtos correctamente

#1. Carga los datos y exáminalos en R. Emplea las funciones head(), summary(), dim() y str(). ¿Cuántas variables hay? ¿Cuántos tratamientos? Hay 2 variables y 5 tratamientos.
dim(data) # esto mide las dimensiones
head(data) #vemos las primeras linesas de la tabla, hay 2 variables
tail(data) # vemos las últimas líneas de la tabla
summary (data)
str (data)

#2. Haz un boxplot para nuestros datos. Uno para cada variable. Colorea a Variable 1 y a Variable 2 de forma diferente (guarda esos colores para las siguientes gráficas)
boxplot(Variable1~Tratamiento, data=data, main="Variable1boxplot", col="pink")
boxplot(Variable2~Tratamiento, data=data, main="Variable2boxplot", col="blue")


#3. Haz un gráfico de dispersión con las dos variables. Cada tratamiento debe de ir de un color distinto. ¡Como en la siguienteimagen! Pista: usa col=datos$Tratamiento
# Crear un gráfico de dispersión con colores según la variable Tratamiento.

plot(data$Variable1, data$Variable2, col = data$Tratamiento, pch = 19, main = "Gráfico de Dispersión", xlab = "Variable1", ylab = "Variable2")

#4. Ponle leyenda al gráfico del apartado anterior. En el margen inferior derecho. Pista: investiga sobre legend()
#Agregar una leyenda en el margen inferior derecho.
legend("bottomright", legend = c("Tratamiento 1", "Tratamiento 2", "Tratamiento 3", "Tratamiento 4", "Tratamiento 5"), fill = c("black", "red", "green", "light blue", "blue"))

#5. Haz un histograma para cada variable. Recuerda mantener los colores.
# Creamos un histograma para Variable1 con colores según Tratamiento.
par(mfrow = c(1, 2))  # Dividimos la ventana gráfica en 1 fila y 2 columnas

#Histograma para Variable1
hist(data$Variable1, col = "pink", main = "Histograma de Variable1",
     xlab = "Variable1", ylab = "Frecuencia")

# Histograma para Variable2
hist(data$Variable2, col = "blue", main = "Histograma de Variable2",
     xlab = "Variable2", ylab = "Frecuencia")

#6. Haz un factor en la columna tratamiento y guárdalo en una variable. Pista: factor(factor$Tratamiento)
# Crear un factor en la columna Tratamiento y lo guardamos en una nueva variable llamada TratamientoFactor.
data$TratamientoFactor <- factor(data$Tratamiento)

# Verificar que se ha creado el nuevo factor.
head(data)  # Muestra las primeras filas del data frame con la nueva columna TratamientoFactor.

#7. Calcula la media y la desviación estándar para cada tratamiento. Recomendación: es más fácil si usas aggregate() o tapply().
#aggregate(Variable~factor,datos,función) (USADO)
#tapply(datos$Variable,factor,función)

resultados_aggregate <- aggregate(. ~ TratamientoFactor, data = data, FUN = function(x) c(Mean = mean(x), SD = sd(x)))

# Mostrar los resultados
print(resultados_aggregate)

#8. Averigua cuántos elementos tiene cada tratamiento. Recomendación: es más fácil si usas table() con el factor. 
# Tiene 10 elementos por cada tratamiento
# Calcular la frecuencia de elementos para cada tratamiento
frecuencia_tratamientos <- table(data$TratamientoFactor)

# Mostrar los resultados
print(frecuencia_tratamientos)

#9. Extrae los datos para el tratamiento 1 y el tratamiento 4 y guárdalos cada uno en una variable diferente.
# Crear un subconjunto de datos para el tratamiento 1
tratamiento1_data <- data[data$TratamientoFactor == 1, ]
tratamiento1_data #vemos si los datos son correctos


# Crear un subconjunto de datos para el tratamiento 4
tratamiento4_data <- data[data$TratamientoFactor == 4, ]
tratamiento4_data #vemos si los datos son correctos

#10. Nuestra hipótesis nula es que las medias de tratamiento 1 y tratamiento 4 para la Variable 1 son iguales. ¿Puedes comprobarlo? Para ello, necesitarás comprobar primero si los datos se distribuyen de forma normal. En función del resultado de la prueba de normalidad, ¿qué test usarías? ** En general, asumimos que las muestras son independientes, pero ¿son sus varianzas iguales? Actúa de acuerdo a tus resultados.

# Realizar un test de normalidad (Shapiro-Wilk) para Tratamiento 1
shapiro_test_1 <- shapiro.test(tratamiento1_data$Variable1)

# Realizar un test de normalidad (Shapiro-Wilk) para Tratamiento 4
shapiro_test_4 <- shapiro.test(tratamiento4_data$Variable1)

# Mostrar los resultados de la prueba de normalidad
cat("Prueba de normalidad para Tratamiento 1 (Variable 1):\n")
print(shapiro_test_1)

cat("\nPrueba de normalidad para Tratamiento 4 (Variable 1):\n")
print(shapiro_test_4)

#Si el valor p (p-value) obtenido en ambas pruebas es mayor que un nivel de significancia (0.05), podrías concluir que los datos siguen una distribución normal. Sin embargo, si el valor p es menor que el nivel de significancia, los datos no seguirían una distribución normal.
#Por lo tanto, los datos tienen una distribución normal.

# Realizar un t-test para comparar las medias de Tratamiento 1 y Tratamiento 4 para Variable 1
t_test_result <- t.test(tratamiento1_data$Variable1, tratamiento4_data$Variable1, var.equal = TRUE)

# Mostrar los resultados del t-test
cat("Prueba t para medias iguales:\n")
print(t_test_result)

#Si el valor p es menor que un nivel de significancia (usualmente 0.05), se rechaza la hipótesis nula (las medias no  son iguales)

#F test to compare two variances

fisher_test_result <- var.test(tratamiento1_data$Variable1, tratamiento4_data$Variable1)

# Mostrar los resultados del test de Fisher
cat("Prueba de igualdad de varianzas (Fisher's F-Test) entre Tratamiento 1 y Tratamiento 4 para Variable 1:\n")
print(fisher_test_result)

#Si el valor p es menor que un nivel de significancia (usualmente 0.05), se rechaza la hipótesis nula (las varianzas no son iguales)
