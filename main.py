import argparse
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import os
from PIL import Image
import io
import logging
import tkinter as tk
from tkinter import ttk, messagebox

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

images_dir = Path("./images").resolve()
try:
    images_dir.mkdir(exist_ok=True, parents=True)
except Exception as e:
    logger.error(f"Error creating images directory: {str(e)}")
    raise

def smiles_to_image(smiles, size=300, format="png"):
    if smiles == '':
        return
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return
        safe_smiles = smiles.replace('/', '_').replace('\\', '_').replace(':', '_')
        path = images_dir / f"{safe_smiles}_{size}.{format}"
        if path.is_file():
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        if format == "svg":
            try:
                drawer = rdMolDraw2D.MolDraw2DSVG(int(size), int(size))
                drawer.drawOptions().clearBackground = False
                drawer.drawOptions().backgroundColour = (1, 1, 1)  # white
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
                with open(path, 'w', encoding='utf-8') as f:
                    f.write(svg)
                return
            except Exception:
                return
        elif format == "jpg":
            img = Draw.MolToImage(mol, size=(int(size), int(size)))
            img = img.convert('RGB')
            buffer = io.BytesIO()
            img.save(buffer, format='JPEG', quality=95)
            with open(path, 'wb') as f:
                f.write(buffer.getvalue())
            return
        else:
            img = Draw.MolToImage(mol, size=(int(size), int(size)))
            buffer = io.BytesIO()
            img.save(buffer, format='PNG')
            with open(path, 'wb') as f:
                f.write(buffer.getvalue())
            return
    except Exception:
        return


def terminal_main():
    print("SMILES Molecule Image Generator")
    smiles = input("Enter SMILES string: ").strip()
    if not smiles:
        return
    try:
        size_input = input("Enter image size (px, default 200): ").strip()
        size = int(size_input) if size_input else 200
        if size < 50 or size > 1000:
            size = 200
    except ValueError:
        size = 200
    fmt = input("Enter image format (png/jpg/svg, default png): ").strip().lower() or "png"
    if fmt not in ["png", "jpg", "svg"]:
        fmt = "png"
    smiles_to_image(smiles, size, fmt)
    safe_smiles = smiles.replace('/', '_').replace('\\', '_').replace(':', '_')
    path = images_dir / f"{safe_smiles}_{size}.{fmt}"
    print(f"Image saved: {path}")

if __name__ == "__main__":
    terminal_main()
