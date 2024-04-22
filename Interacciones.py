import pandas as pd
from Bio import PDB
from Bio.PDB import *
import numpy as np
import yaml
import math
import argparse
import sys
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel
from scipy.spatial.distance import cdist
import mdtraj as md
import dataframe_image as dfi

## Objetos ##

# Definición Clase Receptor
class Receptor:

    def __init__(self, pdb_name,model):
        self.receptor_name = pdb_name
        self.active_site = pd.DataFrame(columns=['Serial', 'Pos', 'Residue', 'Atom', 'X' , 'Y' ,'Z'])
        self.hetatom = pd.DataFrame(columns=['Serial', 'Pos', 'Residue', 'Atom', 'X' , 'Y' ,'Z'])
        self.waters = pd.DataFrame(columns=['Serial', 'Pos', 'Residue', 'Atom', 'X' , 'Y' ,'Z'])
        for residue in model.get_residues():
            for atom in residue:
                Res_name = residue.get_resname()
                Res_id = residue.get_id()[1]
                atom_name = atom.get_name()
                Coor = list(atom.get_coord())
                Serial = atom.get_serial_number()
                if residue.get_id()[0] == ' ':
                    self.active_site.loc[len(self.active_site.index)] = [Serial ,Res_id, Res_name, atom_name ,round(float(Coor[0]),3),round(float(Coor[1]),3),round(float(Coor[2]),3)]
                elif residue.get_id()[0] == 'W':
                    self.waters.loc[len(self.waters.index)] = [Serial ,Res_id, Res_name, atom_name ,round(float(Coor[0]),3),round(float(Coor[1]),3),round(float(Coor[2]),3)]
                else:
                    self.hetatom.loc[len(self.hetatom.index)] = [Serial ,Res_id, Res_name, atom_name ,round(float(Coor[0]),3),round(float(Coor[1]),3),round(float(Coor[2]),3)]
    
    def DF_Receptor(self):
        print(self.active_site)
    
    def DF_H2O(self):
        print(self.waters)

    def DF_Ligands(self):
        print(self.hetatom)

    def Ligands(self):
        if self.hetatom.empty:
            print('Sin Ligandos')
            return(list(self.hetatom['Residue'].unique()))
        else:
            return(list(self.hetatom['Residue'].unique()))

    def Receptor_CA(self):
        return(self.active_site[self.active_site['Atom'] == 'CA'].copy())

    def Generate_ligand(self, lig):
        new_ligand = Ligand(self.hetatom[self.hetatom['Residue'] == lig])
        return new_ligand
    
    def Mass_center(self):
        Res = self.active_site['Pos'].unique()
        active_site_cm = pd.DataFrame(columns=['Pos', 'X' , 'Y' ,'Z'])
        for r in Res:
            Sub_Set = self.active_site.loc[self.active_site['Pos'] == int(r)]
            coordenadas = Sub_Set[['X' , 'Y' ,'Z']].to_numpy()
            cm = centro_de_masa(coordenadas)
            active_site_cm.loc[len(active_site_cm.index)] = [str(r),cm[0] ,cm[1],cm[2] ]
        return(active_site_cm)

    def Gen_Active_Site(self,Pos):
        AS = (self.active_site[self.active_site['Pos'].isin(Pos)])
        return(AS)

    def residue_num(self):
        return(len(self.active_site['Pos'].unique()))
    
    def ligand_num(self):
        return(len(self.hetatom['Pos'].unique()))

    def pdb_name(self):
        print(self.receptor_name)

    def PDB_Active_Site(self,Pos , li , pdb,folder):
        # Cargar el archivo PDB
        traj = md.load(f'PDB/{pdb}.pdb')
        # Seleccionar los residuos especificados
        seleccion = traj.atom_slice(traj.topology.select("resid " + " ".join(map(str, Pos))))

        # Guardar la selección en un nuevo archivo PDB
        ruta_nuevo_pdb = f'{folder}/AS_{pdb}_{li}.pdb'
        seleccion.save(ruta_nuevo_pdb)

class Ligand:
    def __init__(self, lig_name , pdb_file,folder):
        with open(f'PDB/{pdb_file}.pdb', 'r') as f:
            lineas = f.readlines()
            lineas_filtradas = [linea for linea in lineas if lig_name in linea]
        with open(f'{folder}/{lig_name}.pdb', 'w') as f:
            f.writelines(lineas_filtradas)
        mol = next(pybel.readfile("pdb", f'{folder}/{lig_name}.pdb'))
        # Agregar hidrógenos implícitos
        mol.addh()
        # Guardar la molécula con hidrógenos implícitos en un nuevo archivo PDB
        mol.write("pdb", f'{folder}/{lig_name}_H.pdb', overwrite=True)
        number_hydrogens_in_pdb(f'{folder}/{lig_name}_H.pdb',f'{folder}/{lig_name}_H.pdb')
        ## Load File PDB ###
        pdb_parser_lig = PDBParser()
        structure_lig = pdb_parser_lig.get_structure('pdb', f'{folder}/{lig_name}_H.pdb')
        Cad = []
        for a in structure_lig[0]:
            Cad.append(a.id)
        model_lig = structure_lig[0][Cad[0]]
        self.Ligand_DF = pd.DataFrame(columns=['Serial', 'Residue', 'Atom', 'X' , 'Y' ,'Z'])
        for residue in model_lig.get_residues():
            for atom in residue:
                Res_name = lig_name
                Res_id = residue.get_id()[1]
                atom_name = atom.get_name()
                Coor = list(atom.get_coord())
                Serial = atom.get_serial_number()
                self.Ligand_DF.loc[len(self.Ligand_DF.index)] = [Serial , Res_name, atom_name ,round(float(Coor[0]),3),round(float(Coor[1]),3),round(float(Coor[2]),3)]
       
        
    def Hot_points(self):
        dadores = {}
        aceptores = []

        ## Busqueda de dadores (OH / NH / FH)##

        ## Busqueda OH
        # Filtrar los átomos de oxígeno y de hidrógeno
        oxygen_atoms = self.Ligand_DF[self.Ligand_DF['Atom'].str.startswith('O')]
        hydrogen_atoms = self.Ligand_DF[self.Ligand_DF['Atom'].str.startswith('H')]
        dadores = busqueda_dadores(oxygen_atoms,hydrogen_atoms,dadores,0.97)

        ## Busqueda NH
        nitrogen_atoms = self.Ligand_DF[self.Ligand_DF['Atom'].str.startswith('N')]
        hydrogen_atoms = self.Ligand_DF[self.Ligand_DF['Atom'].str.startswith('H')]
        dadores = busqueda_dadores(nitrogen_atoms,hydrogen_atoms,dadores,1.05)

        ## Busqueda fH
        fluor_atoms = self.Ligand_DF[self.Ligand_DF['Atom'].str.startswith('N')]
        hydrogen_atoms = self.Ligand_DF[self.Ligand_DF['Atom'].str.startswith('H')]
        dadores = busqueda_dadores(fluor_atoms,hydrogen_atoms,dadores,0.93)

        ## Aceptores ##
        oxygen_atoms = self.Ligand_DF[self.Ligand_DF['Atom'].str.startswith('O')]
        aceptores_O = list(oxygen_atoms['Atom'].unique())
        nitrogen_atoms = self.Ligand_DF[self.Ligand_DF['Atom'].str.startswith('N')]
        aceptores_N =   list(nitrogen_atoms['Atom'].unique())
        fluor_atoms = self.Ligand_DF[self.Ligand_DF['Atom'].str.startswith('N')]
        aceptores_F =   list(nitrogen_atoms['Atom'].unique())
        aceptores = aceptores_O + aceptores_N
        return(aceptores,dadores)
        
    def find_aromatic(self,pdb,folder):
        
        mol = Chem.MolFromPDBFile('{}'.format(f'{folder}/{pdb}.pdb'))
        # Get the aromatic rings
        sssr = Chem.GetSSSR(mol)

        aromatic_rings = [ring for ring in sssr if len(ring) == 6]  # Filter out rings with fewer than 6 atoms
        # Obtener los átomos participantes de los anillos
        participating_rings = []

        for ring in aromatic_rings:
            participating_rings.append(list(ring))
        
        aromaticos = []
        centers = []
        for j in range(0,len(participating_rings)):
            r = []
            for k in range(0,len(participating_rings[j])):
                aro = (self.Ligand_DF.iloc[participating_rings[j][k],2])
                if any(char.isdigit() for char in aro):
                    r.append(aro.strip())
                else:
                    r.append(str(aro)+str((self.Ligand_DF.iloc[pos,1]).strip()))
            aromaticos.append(r)
            Aromatic_Ring = self.Ligand_DF.loc[self.Ligand_DF['Atom'].isin(r)][['X' , 'Y' , 'Z']]
            centers.append(center_aromatic_ring(np.array(Aromatic_Ring)))
        
        

        return(aromaticos , centers)


    def Mass_center(self):
        Res = self.ligand['Pos'].unique()
        for r in Res:
            Sub_Set = self.ligand.loc[self.ligand['Pos'] == int(r)]
            coordenadas = Sub_Set[['X' , 'Y' ,'Z']].to_numpy()
            cm = centro_de_masa(coordenadas)
        return(cm)


    def DF_Ligand(self):
        print(self.Ligand_DF)
        
    
   


## Funciones ##

def table_to_jpg(DF,folder,out_put_name):
    
    if DF.shape[0] < 100:
                  
        df_styled = DF.style.background_gradient()
        dfi.export(df_styled, f'{folder}/{out_put_name}.jpg')
    else:
        # Define el máximo de filas por grupo
        max_rows_per_group = 100
        # Calcula el número total de grupos necesarios
        num_groups = (len(DF) + max_rows_per_group - 1) // max_rows_per_group

        # Divide el DataFrame en grupos más pequeños
        
        for i in range(num_groups):
            start_idx = i * max_rows_per_group
            end_idx = (i + 1) * max_rows_per_group
            df_group = DF.iloc[start_idx:end_idx, :]
            df_group = df_group.style.background_gradient()
            dfi.export(df_group, f'{folder}/{out_put_name}_{i}.jpg')



def create_folders(carpetas):
    
    # Crear las carpetas una por una
    for carpeta in carpetas:
        if not os.path.exists(carpeta):  # Verificar si la carpeta no existe
            os.makedirs(carpeta)  # Crear la carpeta


def carga_variables():
    # Cargo Variables Generales #
    with open(r'Interacciones_variables.yml') as file:
        Interaciones = yaml.load(file, Loader=yaml.FullLoader)

    Distances_Hidrogen_Bonds =float(Interaciones['distancias']['Distances_Hidrogen_Bonds'])
    Distances_Aromatic = float(Interaciones['distancias']['Distances_Aromatic'])
    Distancia_Hidrofobica = float(Interaciones['distancias']['Distances_Hidrofobica'])
    
    Aceptores_Prot = Interaciones['acceptors']
    Dadores_Prot = Interaciones['donors']
    Aceptot_antecedent = Interaciones['acceptors_antecedent']
    Special_case = Interaciones['special']
    return(Distances_Hidrogen_Bonds,Distances_Aromatic,Distancia_Hidrofobica,Aceptores_Prot,Dadores_Prot,Aceptot_antecedent,Special_case)

def centro_de_masa(coordenadas):
    """
    Calcula el centro de masa para un conjunto de coordenadas (x, y, z).

    Parámetros:
    - coordenadas: array numpy de shape (n, 3) donde n es el número de coordenadas
                   y cada fila contiene las coordenadas (x, y, z) de un punto.

    Devuelve:
    - Un array numpy de 3 elementos que contiene las coordenadas del centro de masa (x, y, z).
    """

    # Calcula el centro de masa ponderado por la masa
    masa_total = coordenadas.shape[0]
    centro = np.sum(coordenadas, axis=0) / masa_total

    return centro

def number_hydrogens_in_pdb(input_pdb_file, output_pdb_file):
    with open(input_pdb_file, 'r') as f:
        pdb_content = f.readlines()

    with open(output_pdb_file, 'w') as f:
        hydrogen_count = 0
        for line in pdb_content:
            if line.startswith("HETATM") and line[12:16].strip() == "H":
                hydrogen_count += 1
                line = line[:12] + "{:>3}".format("H" + str(hydrogen_count)) + line[15:]
            f.write(line)

def busqueda_dadores(oxygen_atoms,hydrogen_atoms,dador,dist):
    # Calcular las distancias entre cada par de átomos
    distances = []
    for i, oxygen_atom in oxygen_atoms.iterrows():
        for j, hydrogen_atom in hydrogen_atoms.iterrows():
            distance = np.linalg.norm(oxygen_atom[['X', 'Y', 'Z']] - hydrogen_atom[['X', 'Y', 'Z']])
            distances.append((oxygen_atom['Atom'], hydrogen_atom['Atom'], distance))

    # Filtrar las distancias menores a 1 Å
    distances_filtered = [(oxygen, hydrogen, distance) for oxygen, hydrogen, distance in distances if distance < dist]

    # Mostrar los pares de átomos con distancia menor a X Å
    for oxygen, hydrogen, distance in distances_filtered:
        #print(f'Distancia entre {oxygen} y {hydrogen}: {distance:.3f} Å')
        dador[oxygen] = [hydrogen , distance]
    return(dador)

def active_site(df1 , df2 , pdb ,lig,folder):
        
        # Convertir las coordenadas XYZ de ambos DataFrames en una matriz numpy
        coords_df1 = df1[['X', 'Y', 'Z']].values
        coords_df2 = df2[['X', 'Y', 'Z']].values

        # Calcular las distancias euclidianas entre todos los puntos de df1 y df2
        distances = cdist(coords_df1, coords_df2)

        # Encontrar el índice del punto más cercano en df2 para cada punto en df1
        indices_punto_mas_cercano = distances.argmin(axis=1)

        # Crear un nuevo DataFrame con la información de los puntos más cercanos
        residue_atom_df1 = df1[['Atom']]
        residue_atom_df2 = df2[['Pos','Residue', 'Atom']]
        distances_df = pd.DataFrame(distances.min(axis=1), columns=['Distance'])

        output_df = pd.concat([residue_atom_df1, residue_atom_df2.iloc[indices_punto_mas_cercano].reset_index(drop=True), distances_df], axis=1)

        # Filtrar los pares de puntos que están a menos de 8 Å de distancia
        output_df = output_df[output_df['Distance'] < 8]

        output_df.to_csv(f'{folder}/Active_Site_Points_{pdb}_{lig}.csv')

        return(output_df['Pos'].unique())

def center_aromatic_ring(Aromatic_Ring):
	x,y,z = [],[],[]

	for  j in range(0,len(Aromatic_Ring)):
		x.append(float(Aromatic_Ring[j][0]))
		y.append(float(Aromatic_Ring[j][1]))
		z.append(float(Aromatic_Ring[j][2]))
		
	CD1 = (x[1],y[1],z[1])
	CE1 = (x[3],y[3],z[3])

	vector_1 = (np.add(CD1, CE1))
	
	CD2 = (x[2],y[2],z[2])
	CE2 = (x[4],y[4],z[4])
		
	vector_2 = (np.add(CD2, CE2))
	
	center = (np.add(vector_2/2, vector_1/2))/2
	return(np.round(center,3))

def Busqueda_Interacciones(caso,points,Ligando,DF_Active,Interacciones):
    Sub_Set = Ligando.loc[Ligando['Atom'] == points][['X', 'Y', 'Z']]
    coords_df1 = Sub_Set.values
    coords_df2 = DF_Active[['X', 'Y', 'Z']].values
    # Calcular las distancias euclidianas entre todos los puntos de df1 y df2
    distances = cdist(coords_df1, coords_df2)
    indices_punto_mas_cercano = distances.argmin(axis=1)
    # Crear un DataFrame temporal con los resultados
    temp_df = DF_Active.iloc[indices_punto_mas_cercano]
    # Crear un DataFrame con los resultados que queremos añadir a Interacciones
    interacciones_df = pd.DataFrame({
        'Interaccion': caso,
        'At Ligand': points,
        'Pos': temp_df['Pos'].values,
        'Residue': temp_df['Residue'].values,
        'Atom': temp_df['Atom'].values,
        'Dist': distances[0][indices_punto_mas_cercano]
    })
    # Añadir este DataFrame a Interacciones
    Interacciones = pd.concat([Interacciones, interacciones_df], ignore_index=True)
    return(Interacciones)

def get_aromatic_coord(Res,AA):
    Aromatic_Ring = []
    if (Res == 'TYR') or (Res == 'PHE'):
            Coordenada = (AA.loc[AA['Atom'] == "CG", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CD1", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CD2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CE1", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CE2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CZ", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
    elif ((Res) == 'TRP') :
            Coordenada = (AA.loc[AA['Atom'] == "CE3", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CD2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CZ3", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CE2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CH2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CZ2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])    
    center = center_aromatic_ring(Aromatic_Ring)
    return(center[0],center[1],center[2])

def center_aromatic_ring(Aromatic_Ring):
	x,y,z = [],[],[]

	for  j in range(0,len(Aromatic_Ring)):
		x.append(float(Aromatic_Ring[j][0]))
		y.append(float(Aromatic_Ring[j][1]))
		z.append(float(Aromatic_Ring[j][2]))
		
	CD1 = (x[1],y[1],z[1])
	CE1 = (x[3],y[3],z[3])

	vector_1 = (np.add(CD1, CE1))
	
	CD2 = (x[2],y[2],z[2])
	CE2 = (x[4],y[4],z[4])
		
	vector_2 = (np.add(CD2, CE2))
	
	center = (np.add(vector_2/2, vector_1/2))/2
	return(np.round(center,3))


def Busqueda_Aromatica(DF_Active,Centros_Lig,Distances_Aromatic,Interacciones):
    Center_Receptor = {}
    Aromaticos = ['TYR' , 'PHE' , 'TRP']
    Res_index = []
    Sub_Set = DF_Active[DF_Active['Residue'].isin(Aromaticos)]
    Pos = list(Sub_Set['Pos'].unique())
    for p in Pos:
        Res = Sub_Set.iloc[0,2]
        AA = Sub_Set[Sub_Set['Pos'] == p]
        Res_index.append(Res)
        coord = get_aromatic_coord(Res,AA)
        Center_Receptor[p] = coord
    values_array = np.array(list(Center_Receptor.values()))
    for k in range(0,len(Centros_Lig)):
        # Calcular las distancias entre el punto y los puntos del conjunto
        distancias = cdist(np.expand_dims(Centros_Lig[k], axis=0), values_array)

        # Encontrar la distancia mínima
        distancia_minima = np.min(distancias)
        indice_punto_cercano = np.argmin(distancias)
        if distancia_minima < Distances_Aromatic:
            Interacciones.loc[len(Interacciones.index)] = 'Aromatica' , f'Anillo {k+1}' ,Pos[indice_punto_cercano] ,Res_index[indice_punto_cercano] ,'dummy' ,distancia_minima
    return(Interacciones)

def Distancia_Interaccion(Interacciones,DF_Active,Ligand_DF,Distances_Hidrogen_Bonds,Acep_antecedent):
    Data_Frame_PH = pd.DataFrame(columns=['Interaccion' , 'At Ligand' , 'Pos' ,  'Residue' , 'Atom' , 'Dist' , 'Angle'])
    Candidatos_PH = Interacciones[Interacciones['Dist'] <= float(Distances_Hidrogen_Bonds)]
    ### Angulos P.H. ###
    ### Corregir bien que es dador y aceptor aggregar aromaticos
    for j in range(0,Candidatos_PH.shape[0]):
        ## Aceptores Lig ##
        if Candidatos_PH.iloc[j,0] == 'Aceptor':
            ## ^H^ -> AnteAcep  - Aceptor - Dador
            # 1 Aceptor Lig
            Coord_Acep_Lig = (np.array(Ligand_DF[Ligand_DF['Atom'] == Candidatos_PH.iloc[j,1]][['X','Y','Z']]))
            # 2 Antecesor Lig
            Coord_Ante_Lig = Antecesor_Ligando(Ligand_DF[Ligand_DF['Atom'].str.startswith('H')],Coord_Acep_Lig)
            # 3 Dador
            filtered_df = DF_Active[(DF_Active['Pos'] == Candidatos_PH.iloc[j,2]) & (DF_Active['Atom'] == Candidatos_PH.iloc[j,4])]
            Coord_Dador = np.array(filtered_df[['X','Y','Z']])
            angle = (calcular_angulo(Coord_Dador, Coord_Acep_Lig, Coord_Ante_Lig))
            Data_Frame_PH.loc[len(Data_Frame_PH.index)]= Candidatos_PH.iloc[j,0] , Candidatos_PH.iloc[j,1] , Candidatos_PH.iloc[j,2] , Candidatos_PH.iloc[j,3] , Candidatos_PH.iloc[j,4] ,Candidatos_PH.iloc[j,5] ,angle
        if Candidatos_PH.iloc[j,0] == 'Dador':
            ## ^H^ -> AnteAcep  - Aceptor - Dador
            # Dador
            Coord_Dador = (np.array(Ligand_DF[Ligand_DF['Atom'] == Candidatos_PH.iloc[j,1]][['X','Y','Z']]))
            # 1 Aceptor 
            Coord_Acep = DF_Active[(DF_Active['Pos'] == Candidatos_PH.iloc[j,2]) & (DF_Active['Atom'] == Candidatos_PH.iloc[j,4])]
            Coord_Acep = np.array(Coord_Acep[['X','Y','Z']])
            # 2 Antecesor depende el caso
            if Candidatos_PH.iloc[j,4] == 'O':
                filtered_df = DF_Active[(DF_Active['Pos'] == Candidatos_PH.iloc[j,2]) & (DF_Active['Atom'] == 'C')]
                Coord_Ant = np.array(filtered_df[['X','Y','Z']])
                angle = (calcular_angulo(Coord_Dador, Coord_Acep, Coord_Ant))
                Data_Frame_PH.loc[len(Data_Frame_PH.index)]= Candidatos_PH.iloc[j,0] , Candidatos_PH.iloc[j,1] , Candidatos_PH.iloc[j,2] , Candidatos_PH.iloc[j,3] , Candidatos_PH.iloc[j,4] ,Candidatos_PH.iloc[j,5] ,angle
            if Candidatos_PH.iloc[j,4] == 'N':
                filtered_df = DF_Active[(DF_Active['Pos'] == Candidatos_PH.iloc[j,2]) & (DF_Active['Atom'] == 'CA')]
                Coord_Ant = np.array(filtered_df[['X','Y','Z']])
                angle = (calcular_angulo(Coord_Dador, Coord_Acep, Coord_Ant))
                Data_Frame_PH.loc[len(Data_Frame_PH.index)]= Candidatos_PH.iloc[j,0] , Candidatos_PH.iloc[j,1] , Candidatos_PH.iloc[j,2] , Candidatos_PH.iloc[j,3] , Candidatos_PH.iloc[j,4] ,Candidatos_PH.iloc[j,5] ,angle
            if Candidatos_PH.iloc[j,3] in list(Acep_antecedent.keys()):
                try:
                    Atomo_Cercano = (Acep_antecedent[Candidatos_PH.iloc[j,3]][Candidatos_PH.iloc[j,3]])
                    filtered_df = DF_Active[(DF_Active['Pos'] == Candidatos_PH.iloc[j,2]) & (DF_Active['Atom'] == Atomo_Cercano)]
                    Coord_Ant = np.array(filtered_df[['X','Y','Z']])
                    angle = (calcular_angulo(Coord_Dador, Coord_Acep, Coord_Ant))
                    Data_Frame_PH.loc[len(Data_Frame_PH.index)]= Candidatos_PH.iloc[j,0] , Candidatos_PH.iloc[j,1] , Candidatos_PH.iloc[j,2] , Candidatos_PH.iloc[j,3] , Candidatos_PH.iloc[j,4] ,Candidatos_PH.iloc[j,5] ,angle
                except KeyError:
                    pass
    return(Data_Frame_PH)

def calcular_angulo(donante, aceptor, antecedente_aceptor):
    donante = np.squeeze(donante)
    aceptor = np.squeeze(aceptor)
    antecedente_aceptor = np.squeeze(antecedente_aceptor)
    # Convertir coordenadas a vectores
    vec_DA = donante - aceptor
    vec_AA = antecedente_aceptor - aceptor

    # Calcular el ángulo utilizando el producto punto
    cos_theta = np.dot(vec_DA, vec_AA) / (np.linalg.norm(vec_DA) * np.linalg.norm(vec_AA))
    angulo_rad = np.arccos(cos_theta)

    # Convertir de radianes a grados
    angulo_grados = np.degrees(angulo_rad)
    return angulo_grados

def Distancia_Aromaticas(Interacciones,DF_Active,Ligand_DF,aro_rings , centers,Data_Frame_PH):
    Puntos_Planos = {'TYR' : ['CD1' ,'CD2'] ,'PHE' : ['CD1' ,'CD2']}
    Candidatos_Aro = Interacciones[Interacciones['Atom'] == 'dummy']
    #get_aromatic_coord(Res,AA) # AA dataframe del TRP
    for j in range(0,Candidatos_Aro.shape[0]):
        Puntos_Plano = []
        Puntos_Plano_Lig = []
        # Busco centro de aa
        filtered_df = DF_Active[(DF_Active['Pos'] == Candidatos_Aro.iloc[j,2])]
        coord_centro = get_aromatic_coord(Candidatos_Aro.iloc[j,3],filtered_df)
        Puntos_Plano.append(np.array(coord_centro))
        for k in Puntos_Planos[Candidatos_Aro.iloc[j,3]]:
            filtered_df = DF_Active[(DF_Active['Pos'] == Candidatos_Aro.iloc[j,2]) & (DF_Active['Atom'] == k)]
            Coord_Ant = np.array(filtered_df[['X','Y','Z']])
            Puntos_Plano.append(Coord_Ant)
        Vector_Prot = calcular_vector_normal(Puntos_Plano[0],Puntos_Plano[1],Puntos_Plano[2])
        # Q anillo busco
        Pos_Anillo = int(Candidatos_Aro.iloc[j,1].split(' ')[1]) - 1
        Puntos_Plano_Lig.append(np.array(centers[Pos_Anillo]))
        for l in range(0,2):
            Coord_Dador = (np.array(Ligand_DF[Ligand_DF['Atom'] == aro_rings[Pos_Anillo][l]][['X','Y','Z']]))
            Puntos_Plano_Lig.append(Coord_Dador)
        Vector_Lig = calcular_vector_normal(Puntos_Plano_Lig[0],Puntos_Plano_Lig[1],Puntos_Plano_Lig[2])
        Angulo_Planos = Angulo_entre_planos(Vector_Prot,Vector_Lig)
        Data_Frame_PH.loc[len(Data_Frame_PH.index)] = Candidatos_Aro.iloc[j,0] , Candidatos_Aro.iloc[j,1] , Candidatos_Aro.iloc[j,2] , Candidatos_Aro.iloc[j,3] , Candidatos_Aro.iloc[j,4] ,Candidatos_Aro.iloc[j,5] ,Angulo_Planos
    return(Data_Frame_PH)


def Angulo_entre_planos(normal_plane1,normal_plane2):
    normal_plane1 = np.squeeze(normal_plane1)
    normal_plane2 = np.squeeze(normal_plane2)

    # Calcula el coseno del ángulo entre los vectores normales
    cos_angle = np.dot(normal_plane1, normal_plane2) / (np.linalg.norm(normal_plane1) * np.linalg.norm(normal_plane2))

    # Calcula el ángulo en radianes
    angle_radians = np.arccos(cos_angle)

    # Convierte el ángulo a grados
    angle_degrees = np.degrees(angle_radians)

    # Si el ángulo es mayor a 90 grados, calcular el ángulo complementario
    if angle_degrees > 90:
        angle_degrees = 180 - angle_degrees

    return angle_degrees



def calcular_vector_normal(punto1, punto2, punto3):
    # Calcula los vectores formados por los puntos
    vector1 = punto2 - punto1
    vector2 = punto3 - punto1
    
    # Calcula el producto vectorial de los dos vectores
    vector_normal = np.cross(vector1, vector2)
    
    return vector_normal


def Antecesor_Ligando(h_atoms,Coord_Dador_Lig):
    # Calcular la distancia euclidiana entre el punto dado y todos los puntos "H"
    h_atoms = h_atoms.copy()
    h_atoms.loc[:, 'Distancia'] = np.linalg.norm(h_atoms[['X', 'Y', 'Z']].values - Coord_Dador_Lig, axis=1)
    # Encontrar el punto "H" más cercano
    h_mas_cercano = h_atoms.loc[h_atoms['Distancia'].idxmin()]
    return(np.array(h_mas_cercano[['X','Y','Z']]))




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Cargar pdb de receptor y cadena de interes')

    # Agregar argumentos
    parser.add_argument('-p', '--PDB', type=str, help='Receptor PDB')
    parser.add_argument('-c', '--Chain', type=str, help='Chain')

    # Parsear los argumentos de la línea de comandos
    args = parser.parse_args()

    if not (args.PDB and args.Chain):
        parser.print_help()
        sys.exit(1)

    # Acceder a las variables ingresadas
    pdb = args.PDB
    chain = args.Chain.upper()

    folder = pdb.upper()
    
    create_folders([folder])


    # Cargo Variables ##
    Distances_Hidrogen_Bonds,Distances_Aromatic,Distancia_Hidrofobica,Aceptores_Prot,Dadores_Prot,Aceptot_antecedent,Special_case = carga_variables()
    
    ## Load File PDB ###
    File = f'PDB/{pdb}.pdb'
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure('pdb', File)
    model = structure[0][chain]

    
    ### Objeto Receptor ##
    Estructura_1 = Receptor(pdb, model)
    CM_DF = Estructura_1.Mass_center()
    #Estructura_1.DF_Receptor() # imprime prot

    ## Precencia de ligandos ##
    Lig = Estructura_1.Ligands()

    if not Lig:
        print("Estructura {} no posee ligando".format(pdb))
        sys.exit()

    for li in Lig:
        print(f'{pdb} - {li}')    
        Ligando = Ligand(li, pdb,folder)
        Lig_Aceptor , Lig_Dador = Ligando.Hot_points()
        aro_rings , centers = (Ligando.find_aromatic(li,folder))

        ## Busqueda Sitio Activo ligando ##

        ## Ver si usar el CA ##
        Pos = active_site(Ligando.Ligand_DF , Estructura_1.active_site , pdb ,li,folder) # Ligando , receptor
        Estructura_1.PDB_Active_Site(Pos , li,pdb,folder)
        DF_Active = (Estructura_1.Gen_Active_Site(Pos))

        ## Calculos Intereacciones ##

        # Crear el DataFrame Interacciones
        Interacciones = pd.DataFrame(columns=['Interaccion', 'At Ligand', 'Pos', 'Residue', 'Atom', 'Dist'])

        # Aceptores #
        for points in Lig_Aceptor:
            Interacciones = Busqueda_Interacciones('Aceptor',points,Ligando.Ligand_DF,DF_Active,Interacciones)

        # Dadores #
        for points in Lig_Dador.keys():
            Interacciones = Busqueda_Interacciones('Dador',points,Ligando.Ligand_DF,DF_Active,Interacciones)

        # Aromaticos #
        Interacciones = Busqueda_Aromatica(DF_Active,centers,Distances_Aromatic,Interacciones)


        ## Validar P.H. ##
        Data_Frame_PH = Distancia_Interaccion(Interacciones,DF_Active,Ligando.Ligand_DF,Distances_Hidrogen_Bonds,Aceptot_antecedent)
        Data_Frame_PH = Distancia_Aromaticas(Interacciones,DF_Active,Ligando.Ligand_DF,aro_rings , centers,Data_Frame_PH)

        ## Salida ##
        print(Data_Frame_PH)
        Data_Frame_PH.to_csv(f'{folder}/{pdb}_interaciones_{li}.csv')
        table_to_jpg(Data_Frame_PH,folder,f'{pdb}_interaciones_{li}')


