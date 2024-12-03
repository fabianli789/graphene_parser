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

import numpy as np
from nomad.config import config
from nomad.datamodel.data import Schema
from nomad.metainfo import Quantity, SchemaPackage, Section, MSection, SubSection
from runschema.run import Run
from runschema.calculation import Calculation

configuration = config.get_plugin_entry_point(
    'graphene_parser.schema_packages:schema_package_entry_point'
)

m_package = SchemaPackage()


class Dimensions_Graphene(MSection):
    m_def = Section(validate=False)
    x = Quantity(type=float, description='x-dimension of 2d lattice in nm')
    y = Quantity(type=float, description='y-dimension of 2d lattice in nm')
    mean_radius = Quantity(type=float, description='Mean Radius of Graphene sample in nm')
    ra = Quantity(type=float, description='Roughness RA in nm')
    rq = Quantity(type=float, description='Roughness RQ in nm')
    hydro_edge = Quantity(type=float, description='Share of Hydrogenated edges in %')
    defects = Quantity(type=float, description='Share of defects in %')

class Concentrations_Graphene(MSection):
    m_def = Section(validate=False)

    name = Quantity(type=str, description="""names of the molecules. "e" means already attached
                                             to graphene edge.""")
    concentration = Quantity(type=np.float64, shape = ['*'], description="""concentration of molecule
                                                                            at specific time.""")

class ChemReactions_Graphene(MSection):
    m_def  = Section(validate=False)
    name = Quantity(type=str, description = 'name of chem. reaction')
    barrier = Quantity(type=float, shape=[], description='energetic barrier in eV')
    occurences = Quantity(type=int, shape=[], description ='number of occurences of this chem. reaction')
    residence_time = Quantity(type=float, shape=[], description =  'time of each chem. reaction')


class GrapheneCalculation(Calculation):
    m_def = Section(validate=False, extends_base_section=False)
    graphene_dimensional_properties = SubSection(sub_section=Dimensions_Graphene.m_def, repeats=False)

    graphene_chem_reactions = SubSection(sub_section=ChemReactions_Graphene.m_def, repeats=True)
    graphene_concentrations = SubSection(sub_section=Concentrations_Graphene.m_def, repeats=True)
    volume_fraction = Quantity(type=float, description='volume if SEI got pressed together')
    porosity = Quantity(type=float, description='share of porous volume wrt to total SEI volume')

    concentration_time = Quantity(type=np.float64, shape=['*'], description="""time evolution for the
                                                                            concentration of molecules,
                                                                            corresponding to
                                                                            "concentration"-array.""")

    pressure_CH4 = Quantity(type=float, description='Pressure of CH4 in Torr')
    pressure_H2 = Quantity(type=float, description='Pressure of H2 in Torr')
    bonds  = Quantity(type=np.int32, shape=['*'], description = 'index of paired atom. -1 means unavailable. Same order as "cartesian_site_coordinates" - array.')
    flag = Quantity(type=np.int32, shape=['*'], description = '# of bonds in graphene flake. Same order as "cartesian_site_coordinates" - array.')
    mean_radius_growth = Quantity(type=np.float64, shape=['*'], description = 'change of mean radius  in nm over time. See mean_radius_growth_time for time steps.')
    mean_radius_growth_time = Quantity(type=np.float64, shape=['*'], description = 'time steps for mean radius growth. Same order of array as "mean_radius_growth".')
    species_coordinates = Quantity(type=np.float64, shape=['*', 3], description='2D cartesian coordinates of the species in nm.')
    species = Quantity(type=str, shape=['*'], description='Type of species, same array length as calculation.species_coordinates.')


m_package.__init_metainfo__()
