# I had to split this up into two scripts because of a bug in gmsh.initialize().
# This script creates the setup for swap mesh.

import sys

from pyelmer import elmer
from pyelmer.gmsh_objects import Model, Shape
from pyelmer.gmsh_utils import *
from pyelmer.execute import run_elmer_grid, run_elmer_solver
from pyelmer.post import scan_logfile

def geometry(crystal_length, name):
    """Gmsh setup of simple CZ growth model."""
    model = Model(name)
    
    # bodies
    rect1 = factory.addRectangle(0, 0, 0, 0.1, 0.1)
    melt = Shape(model, 2, 'melt', [rect1])
    melt.mesh_size = 0.01

    rect2 = factory.addRectangle(0, 0.1, 0, 0.03, crystal_length)
    crystal = Shape(model, 2, 'crystal', [rect2])
    crystal.mesh_size = 0.005
    crystal.set_interface(melt)

    # boundaries
    if_crys_melt = Shape(model, 1, 'if_crys_melt', melt.get_interface(crystal))
    bnd_melt_bt = Shape(model, 1, 'bnd_melt_bt', [melt.bottom_boundary])
    bnd_melt_side = Shape(model, 1, 'bnd_melt_side', [melt.right_boundary])
    bnd_melt_top = Shape(model, 1, 'bnd_melt_top', melt.get_boundaries_in_box([0.03, 0.1], [0.1, 0.1]))

    bnd_crys_top = Shape(model, 1, 'bnd_crys_top', [crystal.top_boundary])
    bnd_crys_side = Shape(model, 1, 'bnd_crys_side', [crystal.right_boundary])
    bnd_sym_ax = Shape(model, 1, 'bnd_sym_ax', model.symmetry_axis)

    # physical groups, mesh
    model.synchronize()
    model.make_physical()
    model.set_const_mesh_sizes()
    model.generate_mesh()

    # model.show()
    model.write_msh(f'./{name}.msh')
    print(model)
    return model

def sif(model, v_pull=40):  # v_pull in mm/min
    """Elmer setup for simple cz growth model."""

    sim = elmer.load_simulation('axi-symmetric_transient')

    # materials
    tin_solid = elmer.load_material('tin_solid', sim)
    tin_liquid = elmer.load_material('tin_liquid', sim)

    # solver, equations
    solver_heat = elmer.load_solver('HeatSolver', sim)
    solver_phasechange = elmer.load_solver('TransientPhaseChange', sim, './elmer_input.yml')
    solver_mesh = elmer.load_solver('MeshUpdate', sim, './elmer_input.yml')
    elmer.load_solver('ResultOutputSolver', sim)
    eqn_main = elmer.Equation(sim, 'eqn_heat', [solver_heat, solver_mesh])
    eqn_phase = elmer.Equation(sim, 'eqn_phase', [solver_phasechange])

    # initial conditions
    ic_melt = elmer.InitialCondition(sim, 'ic_melt', {'Temperature': 510})
    ic_crys = elmer.InitialCondition(sim, 'ic_crys', {'Temperature': 500})
    ic_phase = elmer.InitialCondition(sim, 'ic_phase', {'Temperature': 505,
                                                        'PhaseSurface': 'Real 0.0'})
    # bodies
    melt = elmer.Body(sim, 'melt', [model['melt'].ph_id])
    melt.material = tin_liquid
    melt.equation = eqn_main
    melt.initial_condition = ic_melt

    crystal = elmer.Body(sim, 'crystal', [model['crystal'].ph_id])
    crystal.material = tin_solid
    crystal.equation = eqn_main
    crystal.initial_condition = ic_crys

    phase_if = elmer.Body(sim, 'phase_if', [model['if_crys_melt'].ph_id])
    phase_if.material = tin_solid
    phase_if.equation = eqn_phase
    phase_if.initial_condition = ic_phase

    # boundaries
    bc_phase = elmer.Boundary(sim, 'bc_phase', [model['if_crys_melt'].ph_id])
    bc_phase.normal_target_body = crystal
    bc_phase.smart_heater = True
    bc_phase.smart_heater_T = 505
    bc_phase.phase_change_transient = True
    bc_phase.phase_change_vel = v_pull / 6e4  # mm/min to m/s
    bc_phase.material = tin_solid
    bc_phase.phase_change_body = phase_if

    bc_melt_bt = elmer.Boundary(sim, 'bc_melt_bt', [model['bnd_melt_bt'].ph_id])
    bc_melt_bt.fixed_heatflux = 100  # just a number, is set by heat control anyway
    bc_melt_bt.data.update({'Smart Heater Control': 'Logical True'})
    bc_melt_bt.mesh_update = [0, 0]

    bc_melt_other = elmer.Boundary(sim, 'bc_melt_other', [model['bnd_melt_top'].ph_id, model['bnd_melt_side'].ph_id])
    bc_melt_other.mesh_update = [0, 0]

    bc_crys_side = elmer.Boundary(sim, 'bc_crys_side', [model['bnd_crys_side'].ph_id])
    bc_crys_side.mesh_update = [0, None]

    bc_crys_top = elmer.Boundary(sim, 'bc_crys_top', [model['bnd_crys_top'].ph_id])
    bc_crys_top.fixed_temperature = 500
    bc_crys_top.mesh_update = [0, f'Variable Time\n    real MATC "{v_pull / 6e4 }*tx"']

    bc_sym_ax = elmer.Boundary(sim, 'bc_sim_ax', [model['bnd_sym_ax'].ph_id])
    bc_sym_ax.mesh_update = [0, None]

    # export
    sim.write_startinfo('./')
    sim.write_sif('./')

def run_simple_simulation():
    """Create geometry, run ElmerGrid and ElmerSolver."""
    model = geometry(0.1, 'cz-simple')
    sif(model)
    print('Starting ElmerGrid...')
    run_elmer_grid('./', 'cz-simple.msh')
    print('Starting ElmerSolver...')
    run_elmer_solver('./')
    err, warn, stats = scan_logfile('./')
    print('Errors:', err)
    print('Warnings:', warn)
    print('Statistics', stats)


if __name__ == "__main__":    
    if len(sys.argv) == 4:  # this is called from setup_swapmesh.py
        # This is a workaround to execute the functions in here using
        # a python call instead of importing it. It is required because
        # of a Gmsh bug on Windows, that seems to destroy the path
        # environment / cannot be initialized twice.
        crys_len = float(sys.argv[1])
        name = sys.argv[2]
        v_pull = float(sys.argv[3])
        model = geometry(crys_len, name)
        sif(model, v_pull)
    else:  # this file is executed directly
        run_simple_simulation()

