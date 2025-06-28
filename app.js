document.addEventListener('DOMContentLoaded', () => {
    const videoElement = document.getElementById('video-feed');
    const snapshotsContainer = document.querySelector('.snapshots-container');
    const permissionMessage = document.querySelector('.permission-message');
    snapshotsContainer.style.display = 'flex';
    
    
    let currentFacingMode = 'user';
    let currentStream = null;
    
    
    // Check if the browser supports getUserMedia
    if (!navigator.mediaDevices || !navigator.mediaDevices.getUserMedia) {
        permissionMessage.textContent = 'Camera access is not supported in your browser';
        permissionMessage.style.display = 'block';
        return;
    }
    
    
    async function configureCameras() {
        try {
            const devices = await navigator.mediaDevices.enumerateDevices();
            const videoDevices = devices.filter(device => device.kind === 'videoinput');
            console.log('Found video devices:', videoDevices.length);
            
            
            if (videoDevices.length > 1) {
                const switchCameraBtn = document.createElement('button');
                switchCameraBtn.textContent = 'Switch Camera';
                switchCameraBtn.classList.add('switch-camera-btn');
                snapshotsContainer.appendChild(switchCameraBtn);
                switchCameraBtn.style.display = 'block';
                async function switchCamera() {
                    try {
                        // Stop current stream
                        if (currentStream) {
                            const permissionMessage = document.querySelector('.permission-message');
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
                switchCameraBtn.addEventListener('click', switchCamera);

            }
        } catch (error) {
                console.error('Error checking ameras:', error);
            }
        }
        
        function getVideoConstraints() {
            return {
                video: {
                    width: { min: 1280, ideal: 1920, max: 3840 },
                    height: { min: 720, ideal: 1080, max: 2160 },
                    facingMode: currentFacingMode
                }
            };
        }
        
        // Request camera access
        navigator.mediaDevices.getUserMedia(getVideoConstraints())
        .then(stream => {
            videoElement.srcObject = stream;
            permissionMessage.style.display = 'none';
            configureCameras();
        })
        .catch(error => {
            console.error('Error accessing camera:', error);
            permissionMessage.textContent = 'Camera access denied';
            permissionMessage.style.display = 'block';
        });
        
        videoElement.addEventListener('click', handleInteraction);
        videoElement.addEventListener('touchstart', (e) => {
            e.preventDefault(); // Prevent default touch behavior
            handleInteraction(e);
        });
        
        async function listMolecules(imageBase64, croppedImageBase64, x, y) {
            if (!imageBase64) {
                throw new Error(`${imageBase64} is not a valid base64 image string`);
            }
            
            try {
                const response = await fetch('/list-molecules', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                    },
                    body: JSON.stringify({ 
                        imageBase64, 
                        croppedImageBase64,
                        x: Math.round(x), 
                        y: Math.round(y) 
                    })
                });
                
                const data = await response.json();
                
                if (!response.ok) {
                    throw new Error(`Error analyzing image (HTTP ${response.status}): ${data.error || response.statusText}`);
                }
                
                return data.output;
            } catch (error) {
                console.error('Error:', error);
                throw error;
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
            box.className = 'feedback-box';
            box.style.left = `${x - 50}px`; // Center the box on x coordinate
            box.style.top = `${y - 50}px`; // Center the box on y coordinate

            // Add coordinates text
            const coordsText = document.createElement('div');
            coordsText.className = 'coords-text';
            coordsText.textContent = `${Math.round(x)},${Math.round(y)}`;
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

            // Create a canvas for the cropped image
            const cropCanvas = document.createElement('canvas');
            const cropContext = cropCanvas.getContext('2d');
            cropCanvas.width = 100;
            cropCanvas.height = 100;
            
            // Draw the cropped region
            cropContext.drawImage(
                canvas,
                x - 50, y - 50, 100, 100,  // Source coordinates and size
                0, 0, 100, 100             // Destination coordinates and size
            );
            
            // Get cropped image data
            const croppedImageBase64 = cropCanvas.toDataURL('image/jpeg', 0.8).split(',')[1];

            // Create snapshot container
            const snapshot = document.createElement('div');
            snapshot.className = 'snapshot';
            
            // Create a display canvas that maintains aspect ratio
            const displayCanvas = document.createElement('canvas');
            displayCanvas.width = 160;
            displayCanvas.height = 90;
            const displayContext = displayCanvas.getContext('2d');
            displayContext.drawImage(canvas, 0, 0, canvas.width, canvas.height, 0, 0, 160, 90);
            
            // Create analysis result container
            const analysisContainer = document.createElement('div');
            analysisContainer.className = 'analysis-container';
            
            snapshot.appendChild(displayCanvas);
            snapshot.appendChild(analysisContainer);
            snapshotsContainer.appendChild(snapshot);

            // Chain the analysis, UI updates, and 3D structure generation
            listMolecules(imageBase64, croppedImageBase64, x, y)
                .then(result => {
                    if (Array.isArray(result)) {
                        return generateSDFs(result);
                    }
                    if (result && result.object && Array.isArray(result.smiles)) {
                        objectInput.value = result.object;
                        return generateSDFs(result.smiles);
                    }
                    analysisContainer.textContent = result;
                    analysisContainer.style.color = 'blue';
                })
                .catch(error => {
                    console.error('Analysis error:', error);
                    analysisContainer.textContent = error.message || 'Error analyzing image';
                    analysisContainer.style.color = 'red';
                });

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
});

// Client-side API functions
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
