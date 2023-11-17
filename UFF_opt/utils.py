from openmm import app, unit
import openmm as mm
def create_supercell(pdb_path, scaling_factors,supercell_path):
    """

    The functions generate many atoms with the same name, when openmm read the function again, it will be delete,so I write generate_supercell

    Input:
        pdb_path: pdb file contains the unit cell
        scaling_factors: a tuple or other iterative object (can use scaling_factors[idx])
        supercell_path: write a pdb file including the supercell
    Output:
        None, just write a new function.
    """
    pdbfile = app.PDBFile(pdb_path)

    # Get the original topology and positions from the PDB file
    original_topology = pdbfile.topology
    original_positions = pdbfile.getPositions()

    # Get the original box vectors
    original_box_vectors = original_topology.getPeriodicBoxVectors()

    # Scale the box vectors
    scaled_box_vectors = [
        original_box_vectors[0] * scaling_factors[0],
        original_box_vectors[1] * scaling_factors[1],
        original_box_vectors[2] * scaling_factors[2]
    ]
    # Convert the scaled box vectors to Vec3 objects
    scaled_box_vectors_plain = [vec.value_in_unit(unit.nanometer) for vec in scaled_box_vectors]
    original_box_vectors_plain = [vec.value_in_unit(unit.nanometer) for vec in original_box_vectors]
    original_positions = [vec.value_in_unit(unit.nanometer) for vec in original_positions]
    # Create a new topology with the scaled box vectors and more atoms
    new_topology = app.Topology()
    new_positions = []

    for chain in original_topology.chains():
        new_chain = new_topology.addChain()
        for residue in chain.residues():
            new_residue = new_topology.addResidue(residue.name, new_chain)
            for atom in residue.atoms():
                for i in range(scaling_factors[0]):
                    for j in range(scaling_factors[1]):
                        for k in range(scaling_factors[2]):
                            new_atom = new_topology.addAtom(atom.name, atom.element, new_residue)
                            vec = original_positions[atom.index] + i * original_box_vectors_plain[0] + \
                                            j * original_box_vectors_plain[1] + k * original_box_vectors_plain[2]
                            new_positions.append(vec*10)
    # the new_topology default unit is nanometer, so we can input nanometer value
    new_topology.setPeriodicBoxVectors(scaled_box_vectors_plain)
    print(len(new_positions)) 
    # Here new_positions should be Quantity or vec in unit of angstrom 
    app.PDBFile.writeFile(new_topology, new_positions, open(supercell_path, 'w'))

def generate_supercell(pdb_path, loading_path, scaling_factors,supercell_path):
    """

    The functions generate many atoms with the same name, when openmm read the function again, it will be delete,so I write generate_supercell

    Input:
        pdb_path: pdb file contains the unit cell
        scaling_factors: a tuple or other iterative object (can use scaling_factors[idx])
        supercell_path: write a pdb file including the supercell
    Output:
        None, just write a new function.
    """
    pdbfile = app.PDBFile(pdb_path)

    # Get the original topology and positions from the PDB file
    original_topology = pdbfile.topology
    original_positions = pdbfile.getPositions()

    # Get the original box vectors
    original_box_vectors = original_topology.getPeriodicBoxVectors()

    # Scale the box vectors
    scaled_box_vectors = [
        original_box_vectors[0] * scaling_factors[0],
        original_box_vectors[1] * scaling_factors[1],
        original_box_vectors[2] * scaling_factors[2]
    ]
    # Convert the scaled box vectors to Vec3 objects
    scaled_box_vectors_plain = [vec.value_in_unit(unit.nanometer) for vec in scaled_box_vectors]
    original_box_vectors_plain = [vec.value_in_unit(unit.nanometer) for vec in original_box_vectors]
    original_positions = [vec.value_in_unit(unit.nanometer) for vec in original_positions]
    # Create a new topology with the scaled box vectors and more atoms
    new_topology = app.Topology()
    new_positions = []

    for chain in original_topology.chains():
        new_chain = new_topology.addChain()
        for residue in chain.residues():
            new_residue = new_topology.addResidue(residue.name, new_chain)
            for atom in residue.atoms():
                for i in range(scaling_factors[0]):
                    for j in range(scaling_factors[1]):
                        for k in range(scaling_factors[2]):
                            new_atom = new_topology.addAtom(atom.name, atom.element, new_residue)
                            vec = original_positions[atom.index] + i * original_box_vectors_plain[0] + \
                                            j * original_box_vectors_plain[1] + k * original_box_vectors_plain[2]
                            new_positions.append(vec*10)
    # the new_topology default unit is nanometer, so we can input nanometer value
    new_topology.setPeriodicBoxVectors(scaled_box_vectors_plain)
    print(len(new_positions)) 
    # Here new_positions should be Quantity or vec in unit of angstrom 
    app.PDBFile.writeFile(new_topology, new_positions, open(supercell_path, 'w'))

def cutoff_topology(topo):
    """
    Input: 
        openmm.app.topology.Topology object, simulation.topology or pdb.topology

    Output:
        MOF framework (MOL residue) without bondings
        GAS molecule (GAS residue) without bondings
    """
    subset_topologies = []

    # Iterate over chains and create a new topology for each chain
    # Iterate over chains and create a new topology for each unique residue type within the chain
    for chain in topo.chains():
        unique_residues = set(residue.name for residue in chain.residues())
        
        for unique_residue_name in unique_residues:
            subset_topology = app.Topology()
            
            # Add the chain to the subset topology
            new_chain = subset_topology.addChain(id=chain.id)
            
            # Copy residues and atoms for the unique residue type
            atom_map = {}  # Keep track of the mapping between original atoms and new atoms
            for residue in chain.residues():
                if residue.name == unique_residue_name:
                    new_residue = subset_topology.addResidue(residue.name, new_chain, id=residue.id)
                    for atom in residue.atoms():
                        new_atom = subset_topology.addAtom(atom.name, atom.element, new_residue)
                        atom_map[atom] = new_atom
            
            # Copy bonds for the unique residue type
            for bond in topo.bonds():
                atom1, atom2 = bond
                if atom1 in atom_map and atom2 in atom_map:
                    subset_topology.addBond(atom_map[atom1], atom_map[atom2])  # Add the bond to the subset topology
            subset_topology.setPeriodicBoxVectors(topo.getPeriodicBoxVectors())
            subset_topologies.append(subset_topology)
    return subset_topologies

if __name__ == '__main__':
    pdb=app.PDBFile("MIL120_loading.pdb") #the working directory must be the UFF_opt directory
    frame_top, gas_top = cutoff_topology(pdb.topology)