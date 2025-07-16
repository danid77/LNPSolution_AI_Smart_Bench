function smilesTo3D(smiles, callback) {
    fetch(window.generate3dSdfUrl + "?smiles=" + encodeURIComponent(smiles))
      .then(res => res.text())
      .then(sdf => callback(sdf, smiles))  // SMILES도 함께 넘김
      .catch(err => console.error("3D 생성 실패", err));
}
  
function show3DInModal(sdfText, smiles) {
    const viewerDiv = document.getElementById("viewer3d");
    if (!viewerDiv) return;
    viewerDiv.innerHTML = "";

    const viewer = $3Dmol.createViewer("viewer3d", { backgroundColor: "white" });
    viewer.addModel(sdfText, "sdf");
    viewer.setStyle({}, { stick: {}, sphere: { scale: 0.3 } });
    viewer.zoomTo();
    viewer.render();

    // SMILES 표시
    document.getElementById("modalSmiles").textContent = smiles;

    new bootstrap.Modal(document.getElementById("structureModal3D")).show();
}

function attachSmilesClickHandlers() {
    document.querySelectorAll(".smiles-canvas").forEach(canvas => {
        canvas.addEventListener("click", function () {
            const smiles = this.getAttribute("data-smiles");
            smilesTo3D(smiles, show3DInModal);
        });
    });
}