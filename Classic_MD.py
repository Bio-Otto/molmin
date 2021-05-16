##########################################################################
# IMPORTS
import os

import simtk
from simtk.openmm import app
import simtk.openmm as mm
from simtk.unit import femtosecond, picosecond, nanometer, kelvin, angstrom, atmospheres
from sys import stdout
from apply_pdbfixer import fix_pdb
from simtk.openmm import *


class Classic_MD_Engine:
    def __init__(self, pdb_path, protein_ff=None, water_ff=None, time_step=2.0, nonbondedCutoff=12.0,
                 water_padding=5, total_Steps=250000, temp=310,
                 friction_cofficient=5.0, minimize=True, minimize_steps=5000, just_minimize=False,
                 equilibrate=True, equilibration_step=500, report_interval=500, last_pdb_filename='last.pdb',
                 output_directory=os.getcwd()):

        global platform
        print("Simulation parameters preparing for the start ...")
        speed = 0
        for i in range(Platform.getNumPlatforms()):
            p = Platform.getPlatform(i)
            # print(p.getName(), p.getSpeed())
            if p.getSpeed() > speed:
                platform = p
                speed = p.getSpeed()

        if platform.getName() == 'CUDA' or platform.getName() == 'OpenCL':
            platform.setPropertyDefaultValue('Precision', 'mixed')
            print('Set precision for platform', platform.getName(), 'to mixed')

        self.pdb_path = pdb_path
        self.protein_ff = protein_ff + '.xml'

        if protein_ff == 'charmm36':
            print("Protein FF: %s" % self.protein_ff)
            self.water_ff = '%s/%s.xml' % (protein_ff, water_ff)
            print("Water FF: %s" % self.water_ff)
        else:
            self.water_ff = '%s.xml' % water_ff
            print("Protein FF: %s" % self.protein_ff)
            print("Water FF: %s" % self.water_ff)

        self.time_step = time_step * femtosecond
        self.nonbondedCutoff = nonbondedCutoff * angstrom

        self.water_padding = water_padding * angstrom
        self.total_Steps = total_Steps
        self.temp = temp * kelvin

        self.friction_cofficient = friction_cofficient / picosecond
        self.minimize = minimize
        self.minimize_steps = minimize_steps
        self.just_minimize = just_minimize
        self.equilibrate = equilibrate
        self.equilibration_step = equilibration_step
        self.report_interval = report_interval
        self.last_pdb_filename = last_pdb_filename

        self.output_directory = output_directory

        ## SOLUTION WATER MODEL
        if len((self.water_ff.split('/')[-1]).split('.')[0]) <= 5:
            self.water_model = (self.water_ff.split('/')[-1]).split('.')[0]
            print("Water Model for Solution 1: %s" % self.water_model)

        elif (self.water_ff.split('/')[-1]).split('.')[0] == 'tip4pew':
            self.water_model = (self.water_ff.split('/')[-1]).split('.')[0]
            print("Water Model for Solution 2: %s" % self.water_model)

        else:
            self.water_model = (self.water_ff.split('/')[-1]).split('.')[0][0:5]
            print("Water Model for Solution 3: %s" % self.water_model)

        ## FIXING PDB
        print('pdb file fixing and preparing for simulation ...')
        fixed_pdb_name = fix_pdb(self.pdb_path, self.output_directory)

        print('Loading pdb to simulation engine ...')
        pdb = app.PDBFile(fixed_pdb_name)

        # box = pdb.topology.getUnitCellDimensions()

        print('Modeller of pdb file is preparing ...')
        modeller = mm.app.Modeller(pdb.topology, pdb.positions)
        # modeller.topology.setUnitCellDimensions(box)

        print('Forcefield parameters loading to the simulation system ...')
        forcefield = app.ForceField(self.protein_ff, self.water_ff, 'mg.xml')

        print('Adding missing hydrogens to the model ...')
        modeller.addHydrogens(forcefield)

        print('Adding solvent (both water and ions) to the model to fill a rectangular box ...')
        modeller.addSolvent(forcefield, model=self.water_model, padding=self.water_padding)

        print('Constructing an OpenMM System')
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME,
                                         nonbondedCutoff=self.nonbondedCutoff, constraints=None,
                                         rigidWater=True,
                                         ewaldErrorTolerance=0.00001)

        # self.system.addForce(mm.AndersenThermostat(310 * unit.kelvin, self.friction_cofficient))

        system.addForce(mm.MonteCarloBarostat(1.0 * atmospheres, self.temp, 25))

        print('Creating a LangevinIntegrator.')
        integrator = mm.LangevinIntegrator(self.temp, self.friction_cofficient, self.time_step)

        simulation = app.Simulation(modeller.topology, system, integrator, platform)
        simulation.context.setPositions(modeller.positions)

        simulation.context.computeVirtualSites()

        print('Minimizing...')
        minimize_pdb_path = os.path.join(self.output_directory, 'minimized.pdb')
        simulation.minimizeEnergy()
        print("Minimization done, the energy is", simulation.context.getState(getEnergy=True).getPotentialEnergy())
        positions = simulation.context.getState(getPositions=True).getPositions()
        print("Minimized geometry is written to 'minimized.pdb'")

        app.PDBFile.writeModel(modeller.topology, positions, open(minimize_pdb_path, 'w'), keepIds=True)
