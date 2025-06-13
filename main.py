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

# HTML模板
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>smilesDrawer</title>
    <style>
        body {
            max-width: 800px;
            margin: 0 auto;
            padding: 20px;
        }
        .container {
            padding: 20px;
        }
        h1 {
            text-align: center;
        }
        .input-group {
            margin: 20px 0;
            text-align: center;
        }
        .input-group input[type="text"] {
            width: 350px;
        }
        .controls input[type="number"] {
            width: 50px;
        }
        .result {
            margin-top: 20px;
            text-align: center;
        }
        .controls {
            margin: 10px 0;
            text-align: center;
            display: flex;
            justify-content: center;
            gap: 20px;
        }
        #molecule {
            display: none;
            max-width: 100%;
            margin: 0 auto;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>smilesDrawer</h1>
        <div class="input-group">
            <input type="text" id="smiles" placeholder="SMILES">
            <button onclick="getMolecule()">View</button>
        </div>
        <div class="controls">
            <div>
                <label>Image size:</label>
                <input type="number" id="imageSize" value="200" min="100" max="1000" step="50">
                <span>px</span>
            </div>
            <div>
                <label>Image format:</label>
                <select id="imageFormat">
                    <option value="png">PNG</option>
                    <option value="jpg">JPG</option>
                    <option value="svg">SVG</option>
                </select>
            </div>
        </div>
        <div class="result">
            <img id="molecule">
            <div id="path-info" class="path-info"></div>
        </div>
        <div id="error" class="error"></div>
    </div>

    <script>
        function getMolecule() {
            const smiles = document.getElementById('smiles').value;
            if (!smiles) {
                document.getElementById('error').textContent = 'Input SMILES';
                return;
            }
            
            const img = document.getElementById('molecule');
            const error = document.getElementById('error');
            const pathInfo = document.getElementById('path-info');
            const imageSize = document.getElementById('imageSize').value;
            const imageFormat = document.getElementById('imageFormat').value;
            
            img.style.display = 'none';
            error.textContent = '';
            pathInfo.textContent = '';
            
            fetch(`/smilesRender?smiles=${encodeURIComponent(smiles)}&size=${imageSize}&format=${imageFormat}`)
                .then(response => {
                    if (!response.ok) {
                        throw new Error('Failed to get molecular structure');
                    }
                    return response.blob();
                })
                .then(blob => {
                    const url = URL.createObjectURL(blob);
                    img.src = url;
                    img.style.display = 'block';
                    const safeSmiles = smiles.replace(/[\/\\:]/g, '_');
                    pathInfo.textContent = `Saved to: images/${safeSmiles}_${imageSize}.${imageFormat}`;
                })
                .catch(err => {
                    error.textContent = err.message;
                });
        }

        document.getElementById('smiles').addEventListener('keypress', function(e) {
            if (e.key === 'Enter') {
                getMolecule();
            }
        });
    </script>
</body>
</html>"""

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
    return HTML_TEMPLATE

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
