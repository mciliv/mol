const express = require('express');
const cors = require('cors');
const fs = require('fs');
const path = require('path');
const { spawn } = require('child_process');
const { env } = require('process');

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

        return new Promise((resolve, reject) => {
            const { command, args } = sdf(s, overwrite);
            const pythonProcess = spawn(command, args);

            pythonProcess.stdout.on('data', data =>
                console.log(`Python Output: ${data.toString().trim()}`));
            pythonProcess.stderr.on('data', data =>
                console.error(`Error: ${data.toString().trim()}`));

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

function sdf(s, overwrite) {
    let command = "python"
    let args = ['sdf.py', s, '--dir', SDF_DIR];
    if (overwrite) args.push('--overwrite');
    if (env.PY_DEBUG) {
        command = "debugpy"
        args = ["--listen", "5678", "--wait-for-client"] + args;
    }
    return { command, args };
}

app.listen(PORT, () => console.log(`Node server running on http://localhost:${PORT}`));

app.use((req, res, next) => {
    console.log(`Incoming request: ${req.method} ${req.url}`);
    next();
});

module.exports = app;