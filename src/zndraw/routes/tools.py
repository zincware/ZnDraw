"""Tool endpoints for molecule rendering and other utilities."""

import base64
import typing as t

from fastapi import APIRouter
from pydantic import BaseModel

router = APIRouter(prefix="/v1/tools", tags=["tools"])


class ConvertMoleculeToImageRequest(BaseModel):
    """Request model for molecule image conversion."""

    type: t.Literal["smiles", "inchi"]
    data: str
    dark: bool = False


class ConvertMoleculeToImageResponse(BaseModel):
    """Response model for molecule image conversion."""

    image: str
    status: str


@router.post("/rdkit-img", response_model=ConvertMoleculeToImageResponse)
async def rdkit_image(
    body: ConvertMoleculeToImageRequest,
) -> ConvertMoleculeToImageResponse:
    """Convert SMILES or InChI to a molecule image using RDKit.

    Returns a base64-encoded PNG image of the molecule.
    """
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D

    if body.type == "smiles":
        mol = Chem.MolFromSmiles(body.data)
    else:
        mol = Chem.MolFromInchi(body.data)

    if mol is None:
        return ConvertMoleculeToImageResponse(
            image="",
            status=f"Invalid {body.type} string",
        )

    Chem.rdDepictor.Compute2DCoords(mol)  # type: ignore[attr-defined]

    drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
    opts = drawer.drawOptions()
    opts.clearBackground = False

    if body.dark:
        opts.setSymbolColour((0.9, 0.9, 0.9))
        opts.updateAtomPalette({6: (0.9, 0.9, 0.9)})

    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    png_data = drawer.GetDrawingText()

    img_base64 = base64.b64encode(png_data).decode("utf-8")

    return ConvertMoleculeToImageResponse(
        image=f"data:image/png;base64,{img_base64}",
        status="success",
    )
