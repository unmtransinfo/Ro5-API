#ro5 compute
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import Crippen

# https://rdkit.org/docs/api-docs.html


#  The Ro5 main computation:
#  Compute Rule of Five from a Smiles string.
#  Return dict with fields:
#   smiles, mwt, hbd, hba, logp, violations
def ro5_compute(smiles: str):

    #mwt - avg molecular weight ; mwt <= 500
    #hbd - hydrogen bond donor count ; hbd <= 5
    #hba - hydrogen bond acceptor count ;hba <= 10
    #LogP - predicted octanol-water partition coefficient ; LogP < 5.0
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"smiles": smiles, "error": "Invalid SMILES"}), 400

    
    mwt = Descriptors.MolWt(mol)
    hbd = int(Lipinski.NumHDonors(mol))
    hba = int(Lipinski.NumHAcceptors(mol))
    logp = Crippen.MolLogP(mol)

    

    violations = sum ([mwt>500,hbd>5, hba>10, logp>=5.0])

    return {
        "smiles": smiles,
        "mwt": round(mwt, 3),
        "hbd": hbd,
        "hba": hba,
        "logp": round(logp, 3),
        "violations": int(violations),
    }

