from fastapi import FastAPI
from fastapi.responses import FileResponse, HTMLResponse, Response
from fastapi.staticfiles import StaticFiles
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import datetime
import os
from PIL import Image
import io
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI()

images_dir = Path("./images").resolve()
try:
    images_dir.mkdir(exist_ok=True, parents=True)
    logger.info(f"Images directory created/verified at: {images_dir}")
except Exception as e:
    logger.error(f"Error creating images directory: {str(e)}")
    raise

app.mount("/images", StaticFiles(directory=str(images_dir)), name="images")

@app.get("/", response_class=HTMLResponse)
async def get_index():
    with open("templates/index.html", "r", encoding="utf-8") as f:
        return f.read()

@app.get("/smilesRender")
async def get_image(smiles="", size=300, format="png"):
    if smiles == '':
        return {"error": "Please provide a valid SMILES"}

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES string"}

        safe_smiles = smiles.replace('/', '_').replace('\\', '_').replace(':', '_')
        path = images_dir / f"{safe_smiles}_{size}.{format}"
        
        logger.info(f"Processing SMILES: {smiles}")
        logger.info(f"Target path: {path}")

        if path.is_file():
            logger.info(f"File exists, returning cached version: {path}")
            if format == 'svg':
                with open(path, 'r', encoding='utf-8') as f:
                    content = f.read()
                return Response(content=content, media_type="image/svg+xml")
            else:
                with open(path, 'rb') as f:
                    content = f.read()
                return Response(content=content, media_type=f"image/{format}")

        path.parent.mkdir(parents=True, exist_ok=True)

        if format == "svg":
            try:
                drawer = rdMolDraw2D.MolDraw2DSVG(int(size), int(size))
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
                
                with open(path, 'w', encoding='utf-8') as f:
                    f.write(svg)
                logger.info(f"SVG file saved: {path}")
                
                return Response(content=svg, media_type="image/svg+xml")
            except Exception as e:
                logger.error(f"生成 SVG 时出错: {str(e)}")
                return {"error": f"生成 SVG 时出错: {str(e)}"}
        elif format == "jpg":
            img = Draw.MolToImage(mol, size=(int(size), int(size)))
            img = img.convert('RGB')
            buffer = io.BytesIO()
            img.save(buffer, format='JPEG', quality=95)
            
            with open(path, 'wb') as f:
                f.write(buffer.getvalue())
            logger.info(f"JPG file saved: {path}")
            
            return Response(content=buffer.getvalue(), media_type="image/jpeg")
        else:
            img = Draw.MolToImage(mol, size=(int(size), int(size)))
            buffer = io.BytesIO()
            img.save(buffer, format='PNG')
            
            with open(path, 'wb') as f:
                f.write(buffer.getvalue())
            logger.info(f"PNG file saved: {path}")
            
            return Response(content=buffer.getvalue(), media_type="image/png")

    except Exception as e:
        logger.error(f"Error generating images: {str(e)}")
        return {"error": f"Error generating images: {str(e)}"}
