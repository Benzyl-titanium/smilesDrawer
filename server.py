from fastapi import FastAPI, HTTPException
from fastapi.responses import Response, HTMLResponse
from fastapi.middleware.cors import CORSMiddleware
from rdkit import Chem
from rdkit.Chem import Draw
import io
import logging
import time
import base64

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/", response_class=HTMLResponse)
async def get_index():
    return """
    <!DOCTYPE html>
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
            .error {
                color: red;
                text-align: center;
                margin-top: 10px;
            }
            .download-btn {
                display: none;
                margin-top: 10px;
                padding: 5px 15px;
                background-color: #4CAF50;
                color: white;
                border: none;
                border-radius: 4px;
                cursor: pointer;
            }
            .download-btn:hover {
                background-color: #45a049;
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
                <button id="downloadBtn" class="download-btn" onclick="downloadImage()">Download Image</button>
            </div>
            <div id="error" class="error"></div>
        </div>

        <script>
            let currentImageBlob = null;
            let currentFormat = 'png';

            function getMolecule() {
                const smiles = document.getElementById('smiles').value;
                if (!smiles) {
                    document.getElementById('error').textContent = 'Input SMILES';
                    return;
                }
                
                const img = document.getElementById('molecule');
                const error = document.getElementById('error');
                const downloadBtn = document.getElementById('downloadBtn');
                const imageSize = document.getElementById('imageSize').value;
                const imageFormat = document.getElementById('imageFormat').value;
                currentFormat = imageFormat;
                
                img.style.display = 'none';
                downloadBtn.style.display = 'none';
                error.textContent = '';
                
                fetch(`/generate/${smiles}?size=${imageSize}&format=${imageFormat}`)
                    .then(response => {
                        if (!response.ok) {
                            throw new Error('Failed to get molecular structure');
                        }
                        return response.blob();
                    })
                    .then(blob => {
                        currentImageBlob = blob;
                        const url = URL.createObjectURL(blob);
                        img.src = url;
                        img.style.display = 'block';
                        downloadBtn.style.display = 'inline-block';
                    })
                    .catch(err => {
                        error.textContent = err.message;
                    });
            }

            function downloadImage() {
                if (currentImageBlob) {
                    const url = URL.createObjectURL(currentImageBlob);
                    const a = document.createElement('a');
                    a.href = url;
                    const smiles = document.getElementById('smiles').value;
                    const imageSize = document.getElementById('imageSize').value;
                    const safeSmiles = smiles.replace(/[\/\\:]/g, '_');
                    a.download = `${safeSmiles}_${imageSize}.${currentFormat}`;
                    document.body.appendChild(a);
                    a.click();
                    document.body.removeChild(a);
                    URL.revokeObjectURL(url);
                }
            }

            document.getElementById('smiles').addEventListener('keypress', function(e) {
                if (e.key === 'Enter') {
                    getMolecule();
                }
            });
        </script>
    </body>
    </html>
    """

@app.get("/generate/{smiles}")
async def generate_image(smiles: str, size: int = 300, format: str = "png"):
    try:
        logger.info(f"Received request for SMILES: {smiles}, size: {size}, format: {format}")
        
        # Validate SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Invalid SMILES string: {smiles}")
            raise HTTPException(status_code=400, detail="Invalid SMILES string")

        # Generate image
        img = Draw.MolToImage(mol, size=(int(size), int(size)))
        
        # Convert to bytes
        img_byte_arr = io.BytesIO()
        
        # Save in selected format
        if format.lower() == "svg":
            drawer = rdMolDraw2D.MolDraw2DSVG(int(size), int(size))
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            return Response(content=svg, media_type="image/svg+xml")
        elif format.lower() == "jpg":
            img = img.convert('RGB')
            img.save(img_byte_arr, format='JPEG', quality=95)
            media_type = "image/jpeg"
        else:  # default to PNG
            img.save(img_byte_arr, format='PNG')
            media_type = "image/png"

        img_byte_arr.seek(0)

        # Return image
        return Response(
            content=img_byte_arr.getvalue(),
            media_type=media_type
        )

    except Exception as e:
        logger.error(f"Error generating image: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    try:
        print("Starting server...")
        print("Please access the website at: http://localhost:8000")
        print("For external access, use your server's IP address")
        uvicorn.run(
            app,
            host="0.0.0.0",
            port=8000,
            log_level="info",
            access_log=True
        )
    except Exception as e:
        print(f"Error starting server: {str(e)}")
        raise 