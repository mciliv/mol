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

app.use((req, res, next) => {
    console.log(`Incoming request: ${req.method} ${req.url}`);
    next();
});

app.use(cors());

app.use(express.json({ limit: '50mb' }));

app.use(express.static(__dirname));
app.use('/sdf_files', express.static(SDF_DIR));
app.use('/favicon.ico', express.static(path.join(__dirname, 'favicon.ico')));

app.get('/', (req, res) => {
    res.sendFile(path.join(__dirname, 'index.html'));
});

app.post('/analyze-image', async (req, res) => {
    try {
        const { image: base64ImageData, coordinates } = req.body;
        if (!base64ImageData) {
            return res.status(400).json({ error: 'No image data provided' });
        }

        const headers = {
            'Authorization': `Bearer ${process.env.OPENAI_API_KEY}`,
            'Content-Type': 'application/json'
        };

        const payload = {
            model: 'gpt-4o-mini',
            messages: [
                {
                    role: 'user',
                    content: [
                        {
                            type: 'text',
                            text: `What do you see in this image? I am curious about the materials around us - from everyday objects to complex structures. I want to know what everything is made out of. First identify one main object in the image (around these coordinates: X: ${coordinates?.x}, Y: ${coordinates?.y}.) and then analyze any chemical compounds and answer with the compounds listed as SMILES in a json array dont have json markdown-- nothing else, only an array. X:0 Y:0 is the top left corner of the image. The image is 1000x1000 pixels.`
                        },
                        {
                            type: 'image_url',
                            image_url: {
                                url: `data:image/jpeg;base64,${base64ImageData}`
                            }
                        }
                    ]
                }
            ]
        };

        const response = await fetch('https://api.openai.com/v1/chat/completions', {
            method: 'POST',
            headers: headers,
            body: JSON.stringify(payload)
        });

        if (!response.ok) {
            throw new Error(`OpenAI API error: ${response.statusText}`);
        }

        const result = await response.json();
        console.log(result.choices[0].message.content);

        res.json({
            analysis: result.choices[0].message.content
        });

    } catch (error) {
        console.error('Error:', error);
        res.status(500).json({ error: error.message });
    }
});

app.post('/generate-sdfs', async (req, res) => {
    const { smiles, overwrite = false } = req.body;
    const sdfPaths = [];
    const errors = [];

    const sdfPromises = smiles.map(s => {
        console.log(`Processing SMILES: ${s}`);
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
    return { command, args };
}

app.listen(PORT, () => console.log(`Node server running on http://localhost:${PORT}`));

module.exports = app;