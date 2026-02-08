# slight modification from https://github.com/ibivu/hydrophobic_patches/
# @author: jeff-lafrence, gorantlal
import random
from typing import Optional, TextIO

import networkx as nx
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, Selection
from Bio.PDB.ResidueDepth import get_surface
from scipy.spatial import cKDTree

from ordaos_patterns.bio.protein_language import threetoonedict
from ordaos_patterns.bio.stride_dictionary import return_stride_dictionary


class ResiduePatch():
    """
    Parse pisite files to dict
    Attributes
    ----------
    ids : list
        residue ids -> ('protein', model_id, chain_id, ('', residue_index, ''))
    dssp: dict
        dssp dict -> (amino acid, secondary structure, relative ASA, phi, psi, residue yes/no)
    dssp_dict_keys: dict
        dssp dict keys
    """

    def __init__(self, ids, dssp, dssp_dict_keys):
        self.ids = ids
        self.dssp = dssp
        self.dssp_dict_keys = dssp_dict_keys

    def get_ids(self):
        return self.ids

    def size(self):
        largest_hydro_sasa = 0
        for i in self.ids:
            try:
                largest_hydro_sasa += self.dssp[i[-2:]][2]
            except KeyError:
                if i == self.ids[-1] or i[-1] == self.ids[0]:
                    print(f'{i} is first or last residue, likely incompletely resolved, missing from DSSP dictionary')
                else:
                    print(f'{i} missing from DSSP disctionary, likely incompletely resolved or not split off to chain file')
        return round(largest_hydro_sasa,2)

    def residues(self):
        return [self.dssp[i[-2:]][0] for i in self.ids]

    def patch_length(self):
        return len(self.ids)

    def residue_on_surface(self):
        """
        Get all residues accessible from the surface
        Return
        ------
        list
            residue ids
        """
        return [x for x in self.dssp_dict_keys if x[1][0] == x[1][2] == ' ' and self.dssp[x][2] > 0]

    def check_pisite(self, pisite_dict):
        """
        Get all interaction sites
        Attributes
        ----------
        pisite_dict: pisite dict
        Return
        ------
        int
            interaction sites
        """
        interaction_sites = 0
        for i in self.ids:
            id = i[-1][1]
            if self.dssp[i[-2:]][0] == pisite_dict[id][1] and int(pisite_dict[id][2]) >= 1:
                interaction_sites += 1

        return interaction_sites

    def random_patch_ppis(self, pisite_dict):
        """
        Get all interaction sites
        Attributes
        ----------
        pisite_dict: pisite dict
        Return
        ------
        int
            fraction of interaction site in a random patch
        """
        residues_in_patch = self.patch_length()

        random_selected = [0]*(len(self.residue_on_surface())-residues_in_patch) + [1]*residues_in_patch

        random.shuffle(random_selected)

        random_patch_size = 0

        for i,j in enumerate(self.residue_on_surface()):
            id = j[-1][1]
            if self.dssp[j][0] == pisite_dict[id][1] and random_selected[i] == 1 and int(pisite_dict[id][2]) == 1:
                random_patch_size += 1
        return random_patch_size


class ProteinPatch():
    """
    Parse pisite files to dict
    ...
    Attributes
    ----------
    id : str
        protein id
    pdb_file: str
        pdb file
    cif_file: str
        cif file
    residues_in_patch: list
        list of allowed residues in patch
    r: float
        search radius
    msms: str
        msms command
    """

    def __init__(self, id, pdb_file: Optional[TextIO] = None, cif_file: Optional[TextIO] = None):
        if pdb_file and cif_file:
            raise ValueError('Either pdb_file or cif_file must be specified but not both.')

        residues_in_patch = {'ALA', 'PRO', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'VAL', 'GLY'}
        r = 1.25
        msms = 'msms -density 1.5'

        if pdb_file:
            structure = PDBParser(get_header=True, QUIET=True).get_structure(id, pdb_file)
        elif cif_file:
            structure = MMCIFParser(QUIET=True, auth_residues=True).get_structure(id, cif_file)
        else:
            raise ValueError("input not 'cif' or 'pdb'")

        self.model = structure[0]
        self.r = r
        if not msms.startswith('msms'):
            msms = 'msms ' + msms
        self.msms = msms
        self.residues_in_patch = residues_in_patch
        self.dssp_dict, self.dssp_dict_keys = return_stride_dictionary(structure)
        self.G = self.dot_cloud_graph()
        self.G = self.patch_network(self.G)
        self.patches = self.create_patches()
        self.largest_hydro_patch = self.largest_patch()

    def dot_cloud_graph(self):
        """
        Create a dotted cloud of the proteins surface and label the dots
        hydrophobic or hydrophilic
        ...
        Return
        ------
        networkx.Graph
            Graph with nodes labeled as hydrophobic or not
        """
        surface_points = get_surface(self.model, MSMS=self.msms)
        residue_list = list(Selection.unfold_entities(self.model, 'R'))
        center_vectices = [self._sidechain_center(r.get_atoms()) for r in residue_list]
        T = cKDTree(center_vectices)
        closest_residues = T.query(surface_points, k=1)[1]
        G = nx.Graph()
        for node, coordinates in enumerate(surface_points):
            G.add_node(node)
            G.nodes[node]['selected'] = 0
            closest_residue = residue_list[closest_residues[node]]
            if closest_residue.get_resname() in self.residues_in_patch:
                G.nodes[node]['selected'] = 1
            G.nodes[node]['surface_vector_pos'] = coordinates
            G.nodes[node]['closest_residue_id'] = closest_residue.get_full_id()
            G.nodes[node]['closest_residue_aa'] = closest_residue.get_resname()
        return G

    def patch_network(self, G):
        """
        Create a edges between hydrophobic nodes if they are within r
        hydrophobic or hydrophilic
        ...
        Attributes
        ----------
        G : networkx.Graph
            Graph with nodes labeled as hydrophobic or not
        Return
        ------
        networkx.Graph
            Graph with nodes labeled as hydrophobic
            and edges between nodes within range r
        """
        node_list = [i for i in G.nodes if G.nodes[i]['selected']]
        x = [G.nodes[i]['surface_vector_pos']
             for i in G.nodes if G.nodes[i]['selected']]
        T = cKDTree(x)
        pairs = T.query_pairs(self.r)
        G.add_edges_from([(node_list[x[0]],node_list[x[1]]) for x in pairs])
        return G

    def largest_patch(self):
        """
        get largest patch
        Return
        ------
        asa of largest patch
        residue indexes of largest patch
        percent hydrophobic SASA
        """
        largest_patch = max(self.patches, key=(lambda x: x.size()))
        largest_patch_asa = largest_patch.size()
        largest_patch_residues = sorted(
            [(res[3][0] + str(res[3][1]) + res[3][2]).strip()
              for res in largest_patch.get_ids()])
        total_sasa = sum(value[2] for value in self.dssp_dict.values())
        one_letter_hydro_res = [threetoonedict[res] for res in list(self.residues_in_patch)]
        total_hydro = sum(value[2] for value in self.dssp_dict.values() if value[0] in one_letter_hydro_res)
        percent_hydro = round(total_hydro / total_sasa * 100,2)
        percent_hydro_patch = round(largest_patch_asa / total_sasa * 100,2)
        return [largest_patch_asa, largest_patch_residues, percent_hydro, percent_hydro_patch]

    def largest_patch_evaluator(
        self, minimum_asa: float = 1000, maximum_asa: float = 3000, minimum_percent: float = 3, maximum_percent: float = 8
    ):
        largest_patch_asa, _largest_patch_residues, _percent_hydro, percent_hydro_patch = self.largest_patch()
        h_asa_score = self.get_score(largest_patch_asa, minimum_asa, maximum_asa)
        percent_sasa_hydro_score = self.get_score(percent_hydro_patch, minimum_percent, maximum_percent)

        eval_scores = [h_asa_score, percent_sasa_hydro_score]
        return eval_scores

    def get_score(self, input_value: float, eval_floor: int, eval_ceiling: int) -> float:
        if input_value <= eval_floor:
            return 1
        if input_value >= eval_ceiling:
            return 0
        return float(1-(input_value-eval_floor)/(eval_ceiling-eval_floor))

    def create_patches(self):
        """
        Get the components of the graph
        Return
        ------
        list
            list of Graph components where an item is a patch
        """
        patched_G = self.patch_network(self.G)
        components = nx.connected_components(patched_G)
        patch_dict = []
        for component in components:
            if len(component) <= 1:
                continue
            residue_ids_in_patch = list({self.G.nodes[i]['closest_residue_id'] for i in component})
            patch_dict.append(ResiduePatch(residue_ids_in_patch, self.dssp_dict, self.dssp_dict_keys))
        return sorted(patch_dict, key=(lambda x: x.size()), reverse=True)

    def _sidechain_center(self, atoms):
        vectors = [atom.get_vector().get_array() for atom in atoms]
        center = np.array(vectors).mean(axis=0)
        return center
