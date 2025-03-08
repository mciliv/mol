const express = require('express');
const cors = require('cors');
const fs = require('fs');
const path = require('path');
const { spawn } = require('child_process');

const app = express();
const PORT = 3000;
const SDF_DIR = path.join(__dirname, 'sdf_files');

if (!fs.existsSync(SDF_DIR)) {
    fs.mkdirSync(SDF_DIR, { recursive: true });
}

app.use(cors());
app.use(express.json());

app.post('/generate-sdf', (req, res) => {
    const { smiles, overwrite = false } = req.body;

    if (!smiles) {
        return res.status(400).json({ error: "SMILES string is required" });
    }

    const sdfPath = path.join(SDF_DIR, `${smiles}.sdf`);

    if (fs.existsSync(sdfPath) && !overwrite) {
        return res.json({ sdfPath: `/sdf_files/${smiles}.sdf`, message: "File already exists" });
    }

    const pythonProcess = spawn('poetry' ['run', 'python', 'src/convert/sdf.py', smiles, "--dir", SDF_DIR]);

    pythonProcess.stdout.on('data', (data) => {
        console.log(`Python Output: ${data.toString().trim()}`);
    });

    pythonProcess.stderr.on('data', (data) => {
        console.error(`Error: ${data.toString()}`);
    });

    pythonProcess.on('close', (code) => {
        if (code === 0) {
            res.json({ sdfPath: `/sdf_files/${smiles}.sdf`, message: "File generated" });
        } else {
            res.status(500).json({ error: "SDF generation failed" });
        }
    });
});

app.use('/sdf_files', express.static(SDF_DIR));

app.listen(PORT, () => {
    console.log(`🚀 Server running on http://localhost:${PORT}`);
});

