from typing import (
    TYPE_CHECKING,
)

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )
import os
import re
import datetime
import numpy as np
import math
from pathlib import Path
from itertools import islice

from nomad.config import config
from nomad.datamodel.metainfo.workflow import Workflow
from nomad.parsing.parser import MatchingParser
from runschema.run import Run, Program
from runschema.calculation import Calculation
from graphene_parser.schema_packages.schema_package import GrapheneCalculation, Dimensions_Graphene, ChemReactions_Graphene, Concentrations_Graphene

configuration = config.get_plugin_entry_point(
    'graphene_parser.parsers:parser_entry_point'
)

def DetailedParser(filepath, archive):
    
    if os.path.exists(str(filepath.parent) + r'/status_graphene.csv'):
        with open(str(filepath.parent) + r'/status_graphene.csv') as status_file:
            time_run = archive.m_setdefault("run.time_run")
            time_run.cpu1_start = 0
            calc = archive.m_setdefault("run.calculation")
        
            for i, line in enumerate(status_file):
                line = line.strip("\n")
                parts = line.split(",")
                _cpu = 0
                _time = 0
                _step = 0
                if parts[0] == None:
                    continue
                if re.search(r"cpu", line.lower()):
                    try:
                        _cpu = float(parts[1])
                        time_run.cpu1_end = _cpu
                    except:
                        _cpu = float('nan')
                        time_run.cpu1_end = _cpu
                if 'kmc time' in line.lower():
                    try:
                        _time = float(parts[1])    
                        calc.time = _time
                    except:
                        _time = float('nan')
                        calc.time = _time
                if 'step' in line.lower():
                    try:
                        _step = int(float(parts[1]))
                        calc.step = _step
                    except:
                        _step = -1
                        calc.step = _step
    

    if os.path.exists(str(filepath.parent) + r'/occurrence_res_graphene.csv'):
        with open(str(filepath.parent) + r"/occurrence_res_graphene.csv") as occurrence_file:
            occurence_array = []
            residence_time_array = []
            for i, line in enumerate(occurrence_file):
                if re.search("Oc", line):
                    continue
                parts = line.split(",")
                occurence_array.append(int(parts[1]))
                if len(parts)>=3:
                    residence_time_array.append(float(parts[2]))
                else:
                    pass
    if os.path.exists(str(filepath.parent) + r'/input_graphene.yml'):    
        with open(str(filepath.parent) + r'/input_graphene.yml') as file:
            dim = calc.m_create(Dimensions_Graphene)
            j = 0
            for i, line in enumerate(file):
                parts  = line.split(": ")
               
                if 'x_dim'  in line.lower():
                    dim.x = float(parts[1])
                if 'y_dim' in line.lower():
                    dim.y = float(parts[1])
                if 'p_ch4' in line.lower():
                    calc.pressure_CH4 = float(parts[1])
                if 'p_h2' in line.lower():
                    calc.pressure_H2 = float(parts[1])        
                if re.search(r'\-', line):
                    chem_reactions = calc.m_create(ChemReactions_Graphene)
                    parts[0] = parts[0].lstrip('- ').rstrip(' ')
                    chem_reactions.name = parts[0]
                    chem_reactions.barrier = float(parts[1])
                    
                    chem_reactions.occurences = occurence_array[i-9]
                    if len(residence_time_array) > 0:
                        chem_reactions.residence_time = residence_time_array[i-9]
                    else:
                        pass
    if os.path.exists(str(filepath.parent) + r'/last_step_graphene.csv'):
        with open(str(filepath.parent) + r'/last_step_graphene.csv') as last_step_file:
            species_array = []
            for j, x in enumerate(last_step_file):
                pass
            bonds = []
            flag = []
            coord_x_final = []
            coord_y_final = []
            coordinates_final =  np.zeros((j,3))
            
            last_step_file.seek(0)    
            for i, line in enumerate(last_step_file):
                parts = line.strip("\n").split(",")
                if re.search(r'species', line):
                    continue
                if re.search(r'\d', parts[0]):
                    coord_x_final.append(float(parts[0].strip('"').replace('[', '')))
                else:
                    coord_x_final.append(-1)
                if re.search(r'\d', parts[1]):
                    coord_y_final.append(float(parts[1].strip('"').replace(']', '')))
                else:
                    coord_y_final.append(-1)

                species_array.append(parts[2])
 
                if re.search(r'\d', parts[3]):
                    bonds.append(int(parts[3].strip('"').replace('[', '').replace(']', '')))  
                else:
                    bonds.append(-1)
                if re.search(r'\d', parts[4]):
                   flag.append(int(parts[4].strip('"').replace('[', '').replace(']', '')))  
                else:
                   flag.append(-1)
            
            
            calc.bonds = np.array(bonds)
            calc.flag = np.array(flag)
            coordinates_final[:, 0] = np.array(coord_x_final)
            coordinates_final[:, 1] = np.array(coord_y_final)
            
            calc.species_coordinates = coordinates_final
            calc.species = species_array

                  
    if os.path.exists(str(filepath.parent) + r'/concentration_graphene.csv'):
        with open(str(filepath.parent) + r'/concentration_graphene.csv') as conc_file:
            for x, line in enumerate(conc_file):
                if x == 0:
                    first_line_parts = line.strip("\n").split(",")
                rows = x
            conc_array = np.zeros((rows, len(first_line_parts)-1))
            time_array = []
            conc_file.seek(0) 
            for j, line in enumerate(conc_file):
                if re.search(r'time', line.lower()):
                    continue
                
                parts = line.strip("\n").split(",")
                parts = [float(x) for x in parts]
                time_array.append(parts[0])
                parts = parts[1:]            
                
                # j-1 because of first line
                conc_array[j-1] = parts
            
            if rows > 0:
                calc.concentration_time = np.array(time_array)
                for i in range(len(first_line_parts)-1):
                    conc = calc.m_create(Concentrations_Graphene)   
                    conc.name = first_line_parts[i+1]
                    conc.concentration = conc_array[:,i]
            else:
                pass
    '''
    if os.path.exists(str(filepath.parent) + r'/concentration_graphene_00'):
        file_list = [f'conc_{x}_file' for x in range(100)]        
        row_array = [0 for a in range(100)]
        file_index_array = [f'{j:02d}' for j in range(100)]
        for i in file_index_array:    
            with open(str(filepath.parent) + f'/concentration_graphene_{i}') as file_list[int(i)]:
                for j, line in enumerate(file_list[int(i)]):
                    if j == 0 and i == '00':
                        first_line_parts = line.strip("\n").split(",")
                        continue 
                    else:
                        row_array[int(i)] = j
                    
                file_list[int(i)].seek(0)    
        
        _molecules_name = []                
        _molecules_name = first_line_parts[1:]
        calc.molecules_name = _molecules_name
        time_array = np.zeros(math.floor(sum(row_array)/1000 + 1))
        conc_array = np.zeros((len(time_array), len(first_line_parts) - 1))
        l = 0
        k = 0
        for i in file_index_array:        
            with open(str(filepath.parent) + f'/concentration_graphene_{i}') as file_list[int(i)]:     
                
                for j, line in enumerate(file_list[int(i)]):
                    k += 1
                    if re.search(r'time', line.lower()):
                        continue
                    
                    parts = line.strip("\n").split(",")
                    parts = [float(x) for x in parts]
                    
                    if k % 1000 == 0:
                            
                        time_array[l] = parts[0]
                        parts = parts[1:]
                        conc_array[l] = parts
                        l += 1
                    else:
                        pass
                
        sec_concentration = calc.m_create(Concentrations_Graphene)            
        sec_concentration.concentration = conc_array
        sec_concentration.concentration_time = time_array                            
        '''
    if os.path.exists(str(filepath.parent) + r'/growth_graphene.csv'):
        with open(str(filepath.parent) + r'/growth_graphene.csv') as growth_file:
            mean_radius_growth = []
            mean_radius_growth_time = []
            for i, line in enumerate(growth_file):
                if re.search(r'time', line.lower()):
                    continue
                parts = line.strip('"').split(',')
                try:
                    mean_radius_growth.append(float(parts[0]))
                except:
                    mean_radius_growth.append(-1)
                try:
                    mean_radius_growth_time.append(float(parts[1]))
                except:
                    mean_radius_growth_time.append(-1)
            calc.mean_radius_growth = np.array(mean_radius_growth)
            calc.mean_radius_growth_time = np.array(mean_radius_growth_time)

    if os.path.exists(str(filepath.parent) + r'/properties_graphene.csv'):
        with open(str(filepath.parent) + r'/properties_graphene.csv') as prop_file:
            dim = calc.m_create(Dimensions_Graphene)
            for i, line in enumerate(prop_file):
                line = line.strip("\n")
                parts = line.split(",")
                if parts[0] == None:
                    continue
                if re.search(r'mean', parts[0].lower()):
                    dim.mean_radius = float(parts[1])    
                if re.search(r'ra', parts[0].lower()):
                    dim.ra = float(parts[1])
                if re.search(r'rq', parts[0].lower()):
                    dim.rq = float(parts[1])
                if re.search(r'hy', parts[0].lower()):
                    dim.hydro_edge = float(parts[1])
                if re.search(r'defe', parts[0].lower()):
                    dim.defects = float(parts[1])    

class NewParser(MatchingParser):
    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] = None,
    ) -> None:
        logger.info('NewParser.parse', parameter=configuration.parameter)

        archive.workflow2 = Workflow(name='test')
        mainfile = Path(mainfile)
        DetailedParser(mainfile, archive)