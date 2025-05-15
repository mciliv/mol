document.addEventListener('DOMContentLoaded', () => {
    const videoElement = document.getElementById('video-feed');
    const permissionMessage = document.querySelector('.permission-message');
    const snapshotsContainer = document.createElement('div');
    snapshotsContainer.style.cssText = 'position: fixed; right: 20px; top: 20px; display: flex; flex-direction: column; gap: 10px; z-index: 1000;';
    document.body.appendChild(snapshotsContainer);

    // Create camera switch button
    const switchCameraBtn = document.createElement('button');
    switchCameraBtn.textContent = '🔄 Switch Camera';
    switchCameraBtn.style.cssText = 'position: fixed; left: 50%; top: 20px; transform: translateX(-50%); padding: 12px 20px; background: rgba(0, 0, 0, 0.7); color: white; border: none; border-radius: 25px; font-size: 16px; cursor: pointer; z-index: 1000; transition: background 0.3s;';
    document.body.appendChild(switchCameraBtn);

    let currentFacingMode = 'user';
    let currentStream = null;

    // Check if the browser supports getUserMedia
    if (!navigator.mediaDevices || !navigator.mediaDevices.getUserMedia) {
        permissionMessage.textContent = 'Camera access is not supported in your browser';
        permissionMessage.style.display = 'block';
        return;
    }

    // Function to get video constraints
    function getVideoConstraints() {
        return {
            video: {
                width: { min: 1280, ideal: 1920, max: 3840 },
                height: { min: 720, ideal: 1080, max: 2160 },
                facingMode: currentFacingMode
            }
        };
    }

    // Function to switch camera
    async function switchCamera() {
        try {
            // Stop current stream
            if (currentStream) {
                currentStream.getTracks().forEach(track => track.stop());
            }

            // Toggle facing mode
            currentFacingMode = currentFacingMode === 'user' ? 'environment' : 'user';

            // Get new stream
            const stream = await navigator.mediaDevices.getUserMedia(getVideoConstraints());
            videoElement.srcObject = stream;
            currentStream = stream;
            permissionMessage.style.display = 'none';
        } catch (error) {
            console.error('Error switching camera:', error);
            permissionMessage.textContent = 'Error switching camera';
            permissionMessage.style.display = 'block';
        }
    }

    // Add click event to switch camera button
    switchCameraBtn.addEventListener('click', switchCamera);

    // Handle click/tap events on the video feed
    videoElement.addEventListener('click', handleInteraction);
    videoElement.addEventListener('touchstart', (e) => {
        e.preventDefault(); // Prevent default touch behavior
        handleInteraction(e);
    });

    async function analyzeImage(imageBase64, x, y) {
        try {
            const response = await fetch('/analyze-image', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    image: imageBase64,
                    coordinates: {
                        x: Math.round(x),
                        y: Math.round(y)
                    }
                })
            });

            const data = await response.json();
            if (data.error) {
                throw new Error(data.error);
            }
            return data.analysis;
        } catch (error) {
            console.error('Error analyzing image:', error);
            return 'Error analyzing image';
        }
    }

    function handleInteraction(event) {
        const rect = videoElement.getBoundingClientRect();
        let x, y;

        if (event.type === 'touchstart') {
            const touch = event.touches[0];
            x = touch.clientX - rect.left;
            y = touch.clientY - rect.top;
        } else {
            x = event.clientX - rect.left;
            y = event.clientY - rect.top;
        }

        // Log the coordinates
        console.log(`Interaction at coordinates - X: ${Math.round(x)}, Y: ${Math.round(y)}`);

        // Create and show the feedback box
        const box = document.createElement('div');
        box.style.position = 'absolute';
        box.style.width = '100px';
        box.style.height = '100px';
        box.style.border = '2px solid red';
        box.style.left = `${x - 50}px`; // Center the box on x coordinate
        box.style.top = `${y - 50}px`; // Center the box on y coordinate
        box.style.pointerEvents = 'none'; // Prevent the box from interfering with clicks
        box.style.transition = 'opacity 0.3s ease-out';
        box.style.display = 'flex';
        box.style.justifyContent = 'center';
        box.style.alignItems = 'center';

        // Add coordinates text
        const coordsText = document.createElement('div');
        coordsText.textContent = `${Math.round(x)},${Math.round(y)}`;
        coordsText.style.color = 'red';
        coordsText.style.fontWeight = 'bold';
        coordsText.style.fontSize = '16px';
        box.appendChild(coordsText);

        videoElement.parentElement.appendChild(box);

        // Create a canvas to capture the full video feed
        const canvas = document.createElement('canvas');
        const context = canvas.getContext('2d');
        
        // Set canvas size to match video dimensions
        canvas.width = videoElement.videoWidth;
        canvas.height = videoElement.videoHeight;

        // Draw the entire video frame to the canvas
        context.drawImage(videoElement, 0, 0, canvas.width, canvas.height);

        // Get base64 image data
        const imageBase64 = canvas.toDataURL('image/jpeg', 0.8).split(',')[1];

        // Create snapshot container with fixed display size
        const snapshot = document.createElement('div');
        snapshot.style.cssText = 'width: 160px; background: white; border: 2px solid #333; border-radius: 5px; overflow: hidden; transition: opacity 0.3s ease-out;';
        
        // Create a display canvas that maintains aspect ratio
        const displayCanvas = document.createElement('canvas');
        displayCanvas.width = 160;
        displayCanvas.height = 90;
        const displayContext = displayCanvas.getContext('2d');
        displayContext.drawImage(canvas, 0, 0, canvas.width, canvas.height, 0, 0, 160, 90);
        
        // Create analysis result container
        const analysisContainer = document.createElement('div');
        analysisContainer.style.cssText = 'padding: 8px; font-size: 12px; color: #333;';
        analysisContainer.textContent = 'Analyzing...';
        
        snapshot.appendChild(displayCanvas);
        snapshot.appendChild(analysisContainer);
        snapshotsContainer.appendChild(snapshot);

        // Analyze the image
        analyzeImage(imageBase64, x, y).then(result => {
            analysisContainer.textContent = result;
        }).then(result => generateSDFs(result));

        // Remove the box and snapshot after 10 seconds with fade effect
        setTimeout(() => {
            box.style.opacity = '0';
            snapshot.style.opacity = '0';
            setTimeout(() => {
                box.remove();
                snapshot.remove();
            }, 300);
        }, 15000);
    }

    // Request camera access
    navigator.mediaDevices.getUserMedia(getVideoConstraints())
        .then(stream => {
            videoElement.srcObject = stream;
            currentStream = stream;
            permissionMessage.style.display = 'none';
        })
        .catch(error => {
            console.error('Error accessing camera:', error);
            permissionMessage.textContent = 'Camera access denied';
            permissionMessage.style.display = 'block';
        });
});

async function generateSDFs(smiles) {
    if (!smiles || smiles.length === 0) return;
    try {
        const response = await fetch('/generate-sdfs', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles, overwrite: false })
        });

        const data = await response.json();
        if (data.sdfPaths) {
            load3DmolGrid(data.sdfPaths);
        } else {
            console.error("Error:", data.error);
        }
    } catch (error) {
        console.error("Fetch error:", error);
    }
}

async function load3DmolGrid(sdfFiles) {
    const gldiv = document.getElementById('gldiv');
    gldiv.innerHTML = '';

    const viewers = [];

    for (const sdfFile of sdfFiles) {
        const viewerContainer = document.createElement('div');
        viewerContainer.className = 'mol-viewer-container';
        gldiv.appendChild(viewerContainer);

        const viewer = await render(sdfFile, viewerContainer);
        if (viewer) {
            viewers.push(viewer);
        }
    }

    function resizeAllViewers() {
        viewers.forEach(viewer => {
            viewer.resize()
            viewer.render();
        });
    }

    resizeAllViewers();
    window.addEventListener("resize", resizeAllViewers);
}

async function render(sdfFile, container) {
    try {
        console.log(`Fetching ${sdfFile}`);
        const response = await fetch(sdfFile, { cache: "no-store" });
        if (!response.ok) {
            throw new Error(`HTTP error ${response.status} for ${sdfFile}`);
        }
        const sdfData = await response.text();
        if (!sdfData.trim()) {
            throw new Error(`Empty SDF data for ${sdfFile}`);
        }
        if (!sdfData.includes("$$$$")) {
            throw new Error(`Invalid SDF format for ${sdfFile}`);
        }
        console.log(`SDF data for ${sdfFile}:`, sdfData.substring(0, 100));
        const viewer = $3Dmol.createViewer(container);
        viewer.addModel(sdfData, "sdf");
        viewer.setStyle({}, { sphere: {} });
        viewer.zoomTo();
        viewer.render();
        return viewer;
    } catch (error) {
        console.error(`Failed to load ${sdfFile})`, error);
        container.innerText = `Error: ${error.message}`;
        container.style.color = 'red';
        return null;
    }
}
