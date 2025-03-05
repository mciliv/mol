from rdkit import Chem
from rdkit.Chem import SDWriter

def sdf(smiles, overwrite=False, directory_path=compound_dir()):
    destination = directory_path / Path(smiles).with_suffix(".sdf")
    if overwrite or not destination.exists():
        try:
            mol = Chem.MolFromSmiles(smiles)
            Chem.Kekulize(mol)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            Chem.RemoveHs(mol)
            
            mol.SetProp("SMILES", smiles) 
            
            with open(destination, "w") as file:
                with SDWriter(file) as writer:
                    writer.write(mol)
        except Exception as e:
            logging.exception(f"\n{name}\n{e}\n")
    return destination


def is_valid_sdf_with_molecule(filepath):
    try:
        supplier = Chem.SDMolSupplier(str(filepath))
        for mol in supplier:
            if mol:
                return True
    finally:
        return False


import os
from openbabel import openbabel
from pxr import Usd, UsdGeom, UsdShade, Sdf

def smiles_to_usdz(smiles: str, output_usdz: str):
    # Step 1: Convert SMILES to 3D coordinates using Open Babel
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats("smi", "mol")
    
    molecule = openbabel.OBMol()
    ob_conversion.ReadString(molecule, smiles)
    
    # Generate 3D coordinates
    builder = openbabel.OBBuilder()
    builder.Build(molecule)
    
    # Save as MOL file
    mol_filename = "temp.mol"
    ob_conversion.WriteFile(molecule, mol_filename)
    
    # Step 2: Export MOL to OBJ (or another 3D format)
    obj_filename = "temp.obj"
    os.system(f"obabel {mol_filename} -O {obj_filename}")
    
    # Step 3: Convert OBJ to USDZ
    stage = Usd.Stage.CreateNew(output_usdz)
    root = UsdGeom.Xform.Define(stage, "/Root")
    
    # Create a mesh for the OBJ file
    mesh = UsdGeom.Mesh.Define(stage, "/Root/Molecule")
    mesh.GetPrim().GetReferences().AddReference(obj_filename)
    
    # Optionally add a simple material
    material = UsdShade.Material.Define(stage, "/Root/Material")
    shader = UsdShade.Shader.Define(stage, "/Root/Material/Shader")
    shader.CreateIdAttr("UsdPreviewSurface")
    material.CreateSurfaceOutput().ConnectToSource(shader, "out")
    mesh.GetMaterialBind().ConnectToSource(material)
    
    # Save the USDZ file
    stage.GetRootLayer().Save()
    
    # Cleanup temporary files
    os.remove(mol_filename)
    os.remove(obj_filename)
