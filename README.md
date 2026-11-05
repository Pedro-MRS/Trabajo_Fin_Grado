# Trabajo Fin Grado
## Semejanza funcional y asociación espacial en un matorral gipsófilo
Comprender la relación entre la asociación espacial de especies y la respuesta funcional de las especies a las condiciones ambientales del sitio y a las interacciones bióticas, es importante para conocer los procesos que determinan la estructuración de una comunidad. Se tienen en cuenta procesos como la competencia, que puede producir una segregación espacial de especies con disimilitud de rasgos; y el filtrado ambiental, que puede producir una agregación espacial de especies con similitud de rasgos. Además, la influencia de estos procesos varía con la escala. Sin embargo, en una comunidad de matorral mediterráneo, como la de este estudio, el mecanismo principal es la facilitación, donde existe convergencia funcional de los rasgos que promueven la coexistencia entre especies. Evaluamos estas premisas en un matorral enano mediterráneo cuantificando el grado de asociación espacial en dos escalas diferentes (0,25 m y 2 m), y evaluando su relación con el grado de disimilitud funcional entre especies a partir cinco rasgos importantes como el área foliar especifica, el contenido de materia seca de la hoja, espesor foliar, masa de la semilla y la altura máxima. Nuestra hipótesis es que el filtrado ambiental sea el factor más significativo en el ensamblaje de esta comunidad y que las especies coexistentes presenten similitud entre sus rasgos.
Los resultados muestran que las especies de esta comunidad tienden a agregarse como respuesta a las condiciones ambientales, aunque también una importante proporción de especies se segregan por diferentes preferencias ambientales; y que la relación entre el grado de asociación espacial de especies y el grado de disimilitud de rasgos se correlacionan significativamente y negativamente en ambas escalas. Debido a la influencia del filtrado ambiental en el ensamblaje de esta comunidad, y sugiere la presencia de procesos facilitativos.

## 🛠️ Herramientas bioinformáticas necesarias

Para poder ejecutar los scripts de este repositorio de manera local, es necesario contar con las siguientes herramientas instaladas:

### Entorno de R
* **R** (v4.0 o superior)
* Librerías requeridas: `spatstat`, `idar`, `ggplot2`, `ecespa`, `ecodist`, `FD`.

## 📂 Estructura del Repositorio

El proyecto está organizado de la siguiente manera:
* `pipeline/`: Carpeta destinada a almacenar los diferentes scripts para el procesamiento de los datos.
  * `1_modelo_homogeneo.R`: evaluación del grado de asociación espacial del modelo Poisson homogéneo y representación gráfica de asociaciones.
  * `2_modelo_inhomogeneo.R`: definición de funciones para automatizar los procesos de evaluación del grado de asociación espacial del modelo Poisson inhomogéneo y su representación gráfica.
  * `3_gráficos_inhomogeneo.R`: representación gráfica de asociación espacial del modelo inhomogéneo y estimación de la relación entre la asociación espacial con la disimilitud de rasgos del modelo inhomogéneo mediante el test de Mantel entre las matrices de distancia correspondientes junto con su representación gráfica.
  * `4_mantel_homogeneo.R`: test de Mantel para la estimación de la relación entre la asociación espacial con la disimilitud de rasgos del modelo homogéneo entre las matrices de distancia correspondientes junto con su representación gráfica.
  * `PCAdapt/`: Carpeta donde se guardan los archivos generados por PCAdapt con los resultados de los outliers obtenidos.
  * `iHS/`: Carpeta donde se guardan los archivos generados por iHS con los resultados obtenidos para el análisis de seleccion positiva.
