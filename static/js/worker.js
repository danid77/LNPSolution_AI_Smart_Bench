self.importScripts('https://unpkg.com/@rdkit/rdkit/Code/MinimalLib/dist/RDKit_minimal.js');

let RDKit;

self.onmessage = async function (e) {
  const smiles = e.data;
  
  if (!RDKit) {
    RDKit = await initRDKitModule();
  }

  try {
    const mol = RDKit.get_mol(smiles);
    const molBlock = mol.get_molblock();
    self.postMessage({ success: true, molBlock });
  } catch (err) {
    self.postMessage({ success: false, error: err.message });
  }
};