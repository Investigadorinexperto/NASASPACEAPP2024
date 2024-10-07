import pandas as pd
import math
import numpy as np

#Aquí ingreso los datos en bruto provenientes de la NASA y luego filtro las filas que se repiten para agilizar
datos = pd.read_csv("Datos_estrellas.csv")
datos = datos.drop_duplicates(subset=["hostname", "ra", "dec", "sy_dist"])
datos["sy_dist_ya"] = datos["sy_dist"]*3.26

nombres_estrellas = datos["hostname"].tolist()

estrella_origen = "Kepler-220"

#la idea es analizar cada posibilidad 

#buscar las coordenadas de estrella_origen
coordenadas_estrella_origen = datos[datos["hostname"].str.lower() == estrella_origen.lower()].copy()

#print(f"las coordenadas: {coordenadas_estrella_origen}")

def esfericas_a_cartesianas(radio,a,b): #(radio,a,b) 
    #Recordar que se están calculando en parsecs, por lo que se debe reconvertir a años luz mediante
    x = radio*math.sin(a)*math.cos(b)
    y = radio*math.cos(a)*math.cos(b)
    z = radio*math.sin(b)
    
    return(x,y,z)

estrella_origen_distancia = coordenadas_estrella_origen["sy_dist_ya"].values[0]
estrella_origen_RA = math.radians(coordenadas_estrella_origen["ra"].values[0])
estrella_origen_dec = math.radians(coordenadas_estrella_origen["dec"].values[0])

estrella_origen_coordenadas = esfericas_a_cartesianas(estrella_origen_distancia,estrella_origen_RA, math.pi/2 - estrella_origen_dec)

datos_estrellas = datos[["hostname","sy_dist_ya","ra","dec"]].copy()

distancias_en_ya = datos["sy_dist_ya"].tolist()
estrellas_RA = datos["ra"].tolist()
estrellas_dec = datos["dec"].tolist()

estrellas_RA_radianes = [math.radians(ra) for ra in estrellas_RA]
estrellas_dec_radianes = [math.radians(dec) for dec in estrellas_dec]

coord_x_estrellas = []
coord_y_estrellas = []
coord_z_estrellas = []

for i in range(len(nombres_estrellas)):
    coordenada_cartesiana = esfericas_a_cartesianas(distancias_en_ya[i],estrellas_RA_radianes[i],math.pi/2 - estrellas_dec_radianes[i])
    
    coord_x_estrellas.append(coordenada_cartesiana[0])
    coord_y_estrellas.append(coordenada_cartesiana[1])
    coord_z_estrellas.append(coordenada_cartesiana[2])
    
tabla_coordenadas_estrellas = pd.DataFrame({
    "nombre" : nombres_estrellas,
    "x" : coord_x_estrellas,
    "y" : coord_y_estrellas,
    "z" : coord_z_estrellas
})

#ya se han convertido a cartesianas las coordenadas cilíndricas a cartesianas, por lo que falta filtrar 
#aquellos puntos que se encuentren en una región cilíndrica con extremos en el punto "Sol" (0,0,0) y en
#estrella_origen-220 con un radio de 1 (ya)

centro_base_1 = np.array((0, 0, 0))
centro_base_2 = np.array(estrella_origen_coordenadas)
radio_cilindro = 50 ###################################################################################################

#es un vector que va desde el sol en el centro hasta Kepler en las coordendas de la estrella origen
vector_director = centro_base_2-centro_base_1
long_cilindro = np.linalg.norm(vector_director)

vector_director_normalizado = vector_director/long_cilindro

puntos_xyz = tabla_coordenadas_estrellas[["x", "y", "z"]].values
productos_cruzados = np.cross(vector_director_normalizado, puntos_xyz)
distancias_cilindro = np.linalg.norm(productos_cruzados, axis=1)

puntos_en_cilindro = tabla_coordenadas_estrellas[distancias_cilindro < radio_cilindro]

print(f"Estos son los puntos dentro del cilindro\n\n{puntos_en_cilindro}")

coordenadas_puntos = puntos_en_cilindro[['x', 'y', 'z']].values.tolist()

distancia_a_estrella_origen = []
for index, row in puntos_en_cilindro.iterrows():
    # Calcular la distancia
    distancia = np.linalg.norm(np.array([row['x'], row['y'], row['z']]) - estrella_origen_coordenadas)
    distancia_a_estrella_origen.append(distancia)

# Aquí creo una nueva tabla que compendia el nombre del sistema, sus coordenadas en forma de tupla y 
tabla_pts_en_cilindro_y_dists = pd.DataFrame({
    'hostname': puntos_en_cilindro['nombre'].values,  
    'coordenadas': puntos_en_cilindro[['x', 'y', 'z']].values.tolist(),
    f"distancia a {estrella_origen}": distancia_a_estrella_origen
})

tabla_pts_en_cilindro_y_dists= tabla_pts_en_cilindro_y_dists.sort_values(by=f"distancia a {estrella_origen}")
print(tabla_pts_en_cilindro_y_dists)

tabla_pts_en_cilindro_y_dists.to_csv("resultados_puntos_en_cilindro.csv")

print(f"coordenadas de {estrella_origen}: {estrella_origen_coordenadas}")