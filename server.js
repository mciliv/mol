import express from 'express';
import cors from 'cors';
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import { spawn } from 'child_process';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const app = express();
const PORT = 3000;
const SDF_DIR = path.join(__dirname, 'sdf_files');

if (!fs.existsSync(SDF_DIR)) {
    fs.mkdirSync(SDF_DIR, { recursive: true });
}

app.use(cors());
app.use(express.json());
app.use('/sdf_files', express.static(SDF_DIR));

app.get('/', (req, res) => {
    res.sendFile(path.join(__dirname, 'index.html'));
});

app.post('/generate-sdfs', async (req, res) => {
    const { smiles, overwrite = false } = req.body;
    const sdfPaths = [];
    const errors = [];

    const sdfPromises = smiles.map(s => {
        if (!s) return Promise.resolve();

        const sdfPath = `/sdf_files/${s}.sdf`;
        const fullPath = path.join(SDF_DIR, `${s}.sdf`);

        if (fs.existsSync(fullPath) && !overwrite) {
            sdfPaths.push(sdfPath);
            return Promise.resolve();
        }

        const args = ['sdf.py', s, '--dir', SDF_DIR];
        if (overwrite) args.push('--overwrite');

        return new Promise((resolve, reject) => {
            const pythonProcess = spawn('python', args);

            pythonProcess.stdout.on('data', data => console.log(`Python Output: ${data.toString().trim()}`));
            pythonProcess.stderr.on('data', data => console.error(`Error: ${data.toString().trim()}`));

            pythonProcess.on('close', code => {
                if (code === 0) {
                    sdfPaths.push(sdfPath);
                    resolve();
                } else {
                    const errorMsg = `SDF generation failed for ${s}`;
                    console.error(errorMsg);
                    errors.push(errorMsg);
                    reject(new Error(errorMsg));
                }
            });
        });
    });

    await Promise.allSettled(sdfPromises);

    if (errors.length > 0) {
        return res.status(500).json({ errors });
    }

    res.json({ sdfPaths });
});

app.listen(PORT, () => console.log(`Node server running on http://localhost:${PORT}`));

app.use((req, res, next) => {
    console.log(`Incoming request: ${req.method} ${req.url}`);
    next();
});