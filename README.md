# Modelamiento Bayesiano de redes sociales online de influencia y su impacto en la formación de la opinión pública: Repositorio.

Este trabajo contiene todo el código fuente y los datos necesarios para replicar los resultados del trabajo "Modelamiento Bayesiano de redes sociales online de influencia y su impacto en la formación de la opinión pública" de Samuel Hernando Sánchez Gutiérrez.

Algunos de los archivos `.R` exportan objetos de tipo `.RData` que otros cargan. También pueden exportar archivos `.RData` con la función `save`para leerlos inmediatamente con la función `load`, ésto se hace para evitar tener que correr mútiples veces segmentos de código que pueden tomar horas en ejecutarse. Dichos segmentos se encuentran comentados en algunos casos para evitar correr accidentalmente una tarea de gran magnitud, por tanto, se le recomienda revisar el código comentado pues puede ser importante ejecutarlo para generar los archivos `.RData` que se utlizan más adelante o en otros archivos.

A continuación se listan los archivos indicando a qué sección del documento corresponden y sus dependencias (archivos que se deben ejecutar primero para tener los `.RData` necesarios):

* [`network_example.R`](/network_example.R): Ejemplo aplicado de la sección 1.3.2.4. No tiene dependencias.
* [`gamma_densities_plots.R`](/gamma_densities_plots.R): Generación del gráfico de la figura 2-4. No tiene dependencias.
* [`reading_data.r`](/reading_data.r) : Carga de la red social de la sección 3.1.1. Genera la figura 3-1 y la tabla 3-1. No tiene dependencias. 
* [`jags_implementation.R`](/jags_implementation.R): Ajuste del modelo de influencia a la red social de la Reforma Tributaria y el análisis de los resultados. Corresponde a las secciones 3.1.1.1 a 3.1.1.4. Depende de [`reading_data.r`](/reading_data.r).
* [`aditional_inference_results.R`](/aditional_inference_results.R): Genera la figura 3-7 y la tabla 3-4 de la sección 3.1.1.3. Depende de [`jags_implementation.R`](/jags_implementation.R).
* [`jags_implementation_null_model.R`](/jags_implementation_null_model.R): Ajuste del modelo de la sección 3.1.1.5. Depende de [`reading_data.r`](/reading_data.r).
* [`simulations_for_influence_model.R`](/simulations_for_influence_model.R): Simulación para el modelo de influencia contemplada en la sección 3.1.2. depende de [`jags_implementation.R`](/jags_implementation.R).
* [`simulate_from_influence_model.R`](/simulate_from_influence_model.R): Genera todos los resultados presentados en la sección 3.2 (Estudio de simulación de difusiones en una red de influencia). No tiene dependencias.
* [`cascade.cpp`](/cascade.cpp): Exporta la función que se utiliza para simular las difusiones en [`simulate_from_influence_model.R`](/simulate_from_influence_model.R). Este archivo no tiene dependencias ni debe ser ejecutado por sí solo.
* [`optimize.cpp`](/optimize.cpp): Exporta la función que halla k* en las simulaciones ejecutadas en [`simulate_from_influence_model.R`](/simulate_from_influence_model.R). Este archivo no tiene dependencias ni debe ser ejecutado por sí solo.
