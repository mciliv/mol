document.addEventListener("DOMContentLoaded", () => {
  /* ---------- DOM shortcuts ---------- */
  const video = document.getElementById("video-feed");
  const snapshots = document.querySelector(".snapshots-container");
  const msgBox = document.querySelector(".permission-message");
  const objectInput = document.getElementById("object-input");
  const appContainer = document.querySelector(".app-container");
  const instructionText = document.querySelector(".instruction-text");
  const cameraMode = document.getElementById("camera-mode");
  const photoMode = document.getElementById("photo-mode");
  const cameraContainer = document.querySelector(".camera-container");
  const photoOptions = document.getElementById("photo-options");
  const photoUpload = document.getElementById("photo-upload");
  const photoUrl = document.getElementById("photo-url");
  const urlAnalyze = document.getElementById("url-analyze");

  // Mode switching functionality
  function updateInputMode() {
    if (cameraMode.checked) {
      cameraContainer.style.display = "flex";
      photoOptions.style.display = "none";
    } else {
      cameraContainer.style.display = "none";
      photoOptions.style.display = "flex";
    }
  }

  cameraMode.addEventListener("change", updateInputMode);
  photoMode.addEventListener("change", updateInputMode);
  
  // Initialize mode
  updateInputMode();
  
  // Text input handling - always active
  objectInput.addEventListener("keyup", async (e) => {
    if (e.key !== "Enter") return;
    const object = objectInput.value.trim();
    if (!object) return;

    // Show loading state
    const loadingMsg = createLoadingMessage(`🔍 Analyzing "${object}"...`);
    snapshots.appendChild(loadingMsg);

    try {
      const res = await fetch("/object-molecules", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ object }),
      });
      
      loadingMsg.remove();
      
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const { output } = await res.json();

      // Enhanced result display
      processAnalysisResult(output, snapshots, "📝", object, true);
    } catch (err) {
      loadingMsg.remove();
      
      const h3 = document.createElement("h3");
      h3.textContent = `⚠ Error analyzing "${object}": ${err.message}`;
      h3.style.color = "red";
      snapshots.appendChild(h3);
    }
    objectInput.value = "";
  });

  // Store event handlers for re-attachment
  const photoUploadHandler = async (e) => {
    const file = e.target.files[0];
    if (!file) return;

    // Validate file type
    if (!file.type.startsWith('image/')) {
      alert('Please select an image file');
      return;
    }

    // Display the uploaded image for interactive clicking
    displayUploadedImage(file);
  };

  const photoUrlHandler = (e) => {
    if (e.key === "Enter") {
      urlAnalyzeHandler();
    }
  };

  const urlAnalyzeHandler = () => {
    const url = photoUrl.value.trim();
    if (!url) {
      alert('Please enter an image URL');
      return;
    }
    
    // Basic URL validation
    try {
      new URL(url);
    } catch {
      alert('Please enter a valid URL');
      return;
    }
    
    analyzeImageFromUrl(url);
  };

  // Photo upload functionality
  photoUpload.addEventListener("change", photoUploadHandler);

  // Display uploaded image for interactive clicking
  async function displayUploadedImage(file) {
    const photoOptions = document.getElementById('photo-options');
    
    // Clear existing content
    photoOptions.innerHTML = '';
    
    // Create image display container (similar to camera container)
    const imageContainer = document.createElement('div');
    imageContainer.className = 'uploaded-image-container';
    imageContainer.style.cssText = `
      position: relative;
      flex: 1;
      background: #000;
      display: flex;
      align-items: center;
      justify-content: center;
      cursor: crosshair;
      min-height: 60vh;
      overflow: hidden;
      -webkit-overflow-scrolling: touch;
    `;
    
    // Create image element
    const img = document.createElement('img');
    img.style.cssText = `
      width: 100%;
      height: 100%;
      object-fit: contain;
      background: #000;
      position: relative;
      -webkit-transform: translateZ(0);
      transform: translateZ(0);
    `;
    
    // Create crosshair (same as camera)
    const crosshair = document.createElement('div');
    crosshair.className = 'crosshair';
    crosshair.style.cssText = `
      position: absolute;
      top: 50%;
      left: 50%;
      transform: translate(-50%, -50%);
      width: 20px;
      height: 20px;
      border: 2px solid rgba(255, 255, 255, 0.5);
      border-radius: 50%;
      pointer-events: none;
      z-index: 10;
    `;
    
    // Add crosshair lines
    const beforeLine = document.createElement('div');
    beforeLine.style.cssText = `
      position: absolute;
      top: -10px;
      left: 50%;
      transform: translateX(-50%);
      width: 2px;
      height: 20px;
      background: rgba(255, 255, 255, 0.5);
    `;
    
    const afterLine = document.createElement('div');
    afterLine.style.cssText = `
      position: absolute;
      top: 50%;
      left: -10px;
      transform: translateY(-50%);
      width: 20px;
      height: 2px;
      background: rgba(255, 255, 255, 0.5);
    `;
    
    crosshair.appendChild(beforeLine);
    crosshair.appendChild(afterLine);
    
    // Create instruction text
    const instructionText = document.createElement('div');
    instructionText.className = 'instruction-text';
    instructionText.textContent = 'Click on image parts or type descriptions above';
    instructionText.style.cssText = `
      position: absolute;
      bottom: 60px;
      left: 50%;
      transform: translateX(-50%);
      background: rgba(0, 0, 0, 0.85);
      padding: 12px 20px;
      font-size: 14px;
      font-weight: 500;
      color: rgba(255, 255, 255, 0.95);
      pointer-events: none;
      z-index: 10;
      letter-spacing: 0.01em;
    `;
    
    // Create close button to return to upload options
    const closeButton = document.createElement('button');
    closeButton.textContent = '×';
    closeButton.style.cssText = `
      position: absolute;
      top: 20px;
      right: 20px;
      background: rgba(0, 0, 0, 0.8);
      border: none;
      color: white;
      font-size: 24px;
      width: 40px;
      height: 40px;
      border-radius: 50%;
      cursor: pointer;
      z-index: 20;
      display: flex;
      align-items: center;
      justify-content: center;
      transition: background 0.2s ease;
    `;
    
    closeButton.addEventListener('mouseenter', () => {
      closeButton.style.background = 'rgba(255, 255, 255, 0.2)';
    });
    
    closeButton.addEventListener('mouseleave', () => {
      closeButton.style.background = 'rgba(0, 0, 0, 0.8)';
    });
    
    closeButton.addEventListener('click', () => {
      // Return to upload options
      photoOptions.innerHTML = `
        <div class="upload-option">
          <input type="file" id="photo-upload" accept="image/*" aria-label="Upload photo for molecular analysis">
          <label for="photo-upload" class="upload-label">
            <svg class="photo-icon" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
              <rect x="3" y="3" width="18" height="18" rx="2" ry="2"/>
              <circle cx="8.5" cy="8.5" r="1.5"/>
              <polyline points="21,15 16,10 5,21"/>
            </svg>
            <svg class="upload-text" width="80" height="20" viewBox="0 0 80 20">
              <text x="40" y="14" text-anchor="middle" font-family="system-ui, sans-serif" font-size="11" font-weight="500" fill="currentColor">Upload Photo</text>
            </svg>
          </label>
          
          <div class="url-input-container">
            <input type="url" id="photo-url" placeholder="Paste image URL..." aria-label="Enter image URL for molecular analysis">
            <button type="button" id="url-analyze" class="url-button">
              <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                <path d="M10 13a5 5 0 0 0 7.54.54l3-3a5 5 0 0 0-7.07-7.07l-1.72 1.71"/>
                <path d="M14 11a5 5 0 0 0-7.54-.54l-3 3a5 5 0 0 0 7.07 7.07l1.72-1.71"/>
              </svg>
              <svg class="url-text" width="50" height="16" viewBox="0 0 50 16">
                <text x="25" y="12" text-anchor="middle" font-family="system-ui, sans-serif" font-size="10" font-weight="500" fill="currentColor">Analyze</text>
              </svg>
            </button>
          </div>
        </div>
      `;
      
      // Re-attach event listeners
      const newPhotoUpload = document.getElementById('photo-upload');
      const newPhotoUrl = document.getElementById('photo-url');
      const newUrlAnalyze = document.getElementById('url-analyze');
      
      // Re-attach with the stored handlers
      newPhotoUpload.addEventListener('change', photoUploadHandler);
      newPhotoUrl.addEventListener('keyup', photoUrlHandler);
      newUrlAnalyze.addEventListener('click', urlAnalyzeHandler);
    });
    
    // Load image
    const reader = new FileReader();
    reader.onload = (e) => {
      img.src = e.target.result;
      
      // Store the base64 data for later use
      img.dataset.base64 = e.target.result.split(',')[1];
      
      // Add click handler for image analysis
      imageContainer.addEventListener('click', async (evt) => {
        await handleImageClick(evt, img);
      });
      
      // Add touch handler for mobile
      imageContainer.addEventListener('touchstart', (e) => {
        e.preventDefault();
        handleImageClick(e.touches[0], img);
      });
    };
    
    reader.readAsDataURL(file);
    
    // Assemble the container
    imageContainer.appendChild(img);
    imageContainer.appendChild(crosshair);
    imageContainer.appendChild(instructionText);
    imageContainer.appendChild(closeButton);
    
    // Add to photo options
    photoOptions.appendChild(imageContainer);
  }

  // Handle clicks on uploaded image (similar to camera interaction)
  async function handleImageClick(evt, img) {
    // Hide instruction text on first click
    const instructionText = document.querySelector('.uploaded-image-container .instruction-text');
    if (instructionText) {
      instructionText.style.opacity = "0";
      instructionText.style.transition = "opacity 0.5s ease";
    }

    // Get click coordinates relative to the image
    const rect = img.getBoundingClientRect();
    const clickX = evt.clientX - rect.left;
    const clickY = evt.clientY - rect.top;
    
    // Convert to relative coordinates (0-1)
    const relativeX = clickX / rect.width;
    const relativeY = clickY / rect.height;
    
    // Get the base64 data
    const imageBase64 = img.dataset.base64;
    
    // Create a canvas to crop the image around the click point
    const canvas = document.createElement('canvas');
    const ctx = canvas.getContext('2d');
    
    // Load the image into canvas
    const tempImg = new Image();
    tempImg.onload = () => {
      canvas.width = tempImg.width;
      canvas.height = tempImg.height;
      ctx.drawImage(tempImg, 0, 0);
      
      // Calculate crop area (similar to camera logic)
      const cropSize = Math.min(tempImg.width, tempImg.height) * 0.1; // 10% of image size
      const cropX = Math.max(0, relativeX * tempImg.width - cropSize / 2);
      const cropY = Math.max(0, relativeY * tempImg.height - cropSize / 2);
      
      // Create cropped canvas
      const cropCanvas = document.createElement('canvas');
      cropCanvas.width = cropSize;
      cropCanvas.height = cropSize;
      const cropCtx = cropCanvas.getContext('2d');
      
      cropCtx.drawImage(
        canvas, 
        cropX, cropY, cropSize, cropSize,
        0, 0, cropSize, cropSize
      );
      
      const croppedBase64 = cropCanvas.toDataURL('image/jpeg', 0.9).split(',')[1];
      
      // Show loading state
      const loadingMsg = createLoadingMessage("🔍 Analyzing...");
      snapshots.appendChild(loadingMsg);
      
      // Send for analysis
      fetch("/image-molecules", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          imageBase64,
          croppedImageBase64: croppedBase64,
          x: relativeX * tempImg.width,
          y: relativeY * tempImg.height,
        }),
      })
      .then(res => {
        loadingMsg.remove();
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        return res.json();
      })
      .then(({ output }) => {
        const objectName = output.object || "Uploaded image";
        processAnalysisResult(output, snapshots, "📁", objectName);
      })
      .catch(err => {
        loadingMsg.remove();
        const errorMsg = document.createElement("h3");
        errorMsg.textContent = `⚠ Error: ${err.message}`;
        errorMsg.style.color = "red";
        snapshots.appendChild(errorMsg);
      });
    };
    
    tempImg.src = `data:image/jpeg;base64,${imageBase64}`;
  }

  // Photo URL functionality
  async function analyzeImageFromUrl(url) {
    try {
      // Fetch and convert image to base64
      const imageBase64 = await urlToBase64(url);
      
      // Create a blob from the base64 data to display the image
      const byteCharacters = atob(imageBase64);
      const byteNumbers = new Array(byteCharacters.length);
      for (let i = 0; i < byteCharacters.length; i++) {
        byteNumbers[i] = byteCharacters.charCodeAt(i);
      }
      const byteArray = new Uint8Array(byteNumbers);
      const blob = new Blob([byteArray], { type: 'image/jpeg' });
      const file = new File([blob], 'url-image.jpg', { type: 'image/jpeg' });
      
      // Display the image for interactive clicking
      displayUploadedImage(file);
      
      // Clear URL input
      photoUrl.value = '';

    } catch (err) {
      const errorMsg = document.createElement("h3");
      errorMsg.textContent = `⚠ Error loading image from URL: ${err.message}`;
      errorMsg.style.color = "red";
      snapshots.appendChild(errorMsg);
    }
  }

  // URL input event listeners
  urlAnalyze.addEventListener("click", urlAnalyzeHandler);
  photoUrl.addEventListener("keyup", photoUrlHandler);

  // Helper function to convert file to base64
  function fileToBase64(file) {
    return new Promise((resolve, reject) => {
      const reader = new FileReader();
      reader.onload = () => {
        const result = reader.result.split(',')[1]; // Remove data:image/...;base64, prefix
        resolve(result);
      };
      reader.onerror = reject;
      reader.readAsDataURL(file);
    });
  }

  // Helper function to convert URL to base64
  async function urlToBase64(url) {
    return new Promise((resolve, reject) => {
      const img = new Image();
      img.crossOrigin = 'anonymous'; // Try to enable CORS
      
      img.onload = () => {
        const canvas = document.createElement('canvas');
        const ctx = canvas.getContext('2d');
        
        canvas.width = img.width;
        canvas.height = img.height;
        
        ctx.drawImage(img, 0, 0);
        
        try {
          const dataURL = canvas.toDataURL('image/jpeg', 0.9);
          const base64 = dataURL.split(',')[1];
          resolve(base64);
        } catch (err) {
          reject(new Error('Failed to convert image to base64'));
        }
      };
      
      img.onerror = () => {
        reject(new Error('Failed to load image from URL. Check if URL is valid and accessible.'));
      };
      
      img.src = url;
    });
  }

  /* ---------- Helper Functions ---------- */
  
  // Create loading message element
  function createLoadingMessage(text) {
    const loadingMsg = document.createElement("h3");
    loadingMsg.textContent = text;
    loadingMsg.style.fontStyle = "italic";
    loadingMsg.style.opacity = "0.7";
    return loadingMsg;
  }
  
  // Generate result message for molecular analysis
  function createResultMessage(icon, objectName, smilesCount, useQuotes = false) {
    const name = useQuotes ? `"${objectName}"` : objectName;
    const plural = smilesCount !== 1 ? 's' : '';
    return `${icon} ${name} → ${smilesCount} molecule${plural} found`;
  }

  // Handle result processing for all analysis types
  function processAnalysisResult(output, container, icon, objectName, useQuotes = false) {
    const smilesCount = output.smiles ? output.smiles.length : 0;
    const result = document.createElement("h3");
    result.textContent = createResultMessage(icon, objectName, smilesCount, useQuotes);
    container.appendChild(result);

    if (output.smiles && output.smiles.length > 0) {
      generateSDFs(output.smiles, objectName);
    }
    
    return result;
  }

  /* ---------- Camera setup ---------- */
  console.log("🐛 Starting camera setup...");
  console.log("Navigator object:", {
    mediaDevices: !!navigator.mediaDevices,
    getUserMedia: !!navigator.mediaDevices?.getUserMedia,
    userAgent: navigator.userAgent,
    location: location.href
  });
  
  // Simple check - if we get here, the browser supports cameras
  if (!navigator.mediaDevices || !navigator.mediaDevices.getUserMedia) {
    console.error("❌ Camera API not available");
    msgBox.hidden = false;
    msgBox.textContent = "Camera API not supported in this browser";
    return;
  }
  
  console.log("✅ Camera API is available, proceeding...");

  // Check for HTTPS requirement (especially important for mobile)
  const isSecureContext = window.isSecureContext || location.protocol === 'https:';
  if (!isSecureContext && !location.hostname.includes('localhost')) {
    msgBox.hidden = false;
    msgBox.textContent = "⚠️ Camera requires HTTPS on mobile devices. Please use HTTPS or localhost.";
    return;
  }

  // Mobile detection for debugging
  const isMobile = /Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(navigator.userAgent);
  console.log(`Device type: ${isMobile ? 'Mobile' : 'Desktop'}`);
  console.log(`User agent: ${navigator.userAgent}`);
  console.log(`Secure context: ${isSecureContext}`);

  let facingMode = isMobile ? "environment" : "user";
  let currentStream = null;



  // Simple constraints - only two levels needed
  const simpleConstraints = () => ({ video: { facingMode } });
  const basicConstraints = () => ({ video: true });

  // Simplified camera start
  async function startCamera() {
    currentStream?.getTracks().forEach(t => t.stop());

    try {
      // Try simple constraints first, then basic fallback
      let stream;
      try {
        stream = await navigator.mediaDevices.getUserMedia(simpleConstraints());
        console.log("Using simple constraints");
      } catch (err) {
        console.log("Simple failed, trying basic:", err.message);
        stream = await navigator.mediaDevices.getUserMedia(basicConstraints());
        console.log("Using basic constraints");
      }
      
      currentStream = stream;
      video.srcObject = stream;
      await video.play();
      msgBox.hidden = true;
      console.log("✅ Camera started");
      
    } catch (err) {
      console.error("Camera error:", err);
      msgBox.hidden = false;
      msgBox.textContent = `📷 Camera error: ${err.message}`;
    }
  }



  async function setupSwitchCamera() {
    const devices = await navigator.mediaDevices.enumerateDevices();
    if (devices.filter((d) => d.kind === "videoinput").length < 2) return;

    const btn = document.getElementById("switch-camera-btn");
    btn.style.display = "flex";
    btn.onclick = () => {
      facingMode = facingMode === "user" ? "environment" : "user";
      console.log(`Switching to ${facingMode} camera`);
      startCamera();
    };
  }

  /* ---------- Enhanced interaction handler ---------- */
  video.addEventListener("click", handleInteraction);
  video.addEventListener("touchstart", (e) => {
    e.preventDefault();
    handleInteraction(e.touches[0]);
  });

  async function handleInteraction(evt) {
    if (!cameraMode.checked) return; // Only handle interaction if in camera mode
    
    // Hide instruction text on first click
    if (instructionText) {
      instructionText.style.opacity = "0";
      instructionText.style.transition = "opacity 0.5s ease";
    }

    // Capture high-quality frame
    const canvas = document.createElement("canvas");
    canvas.width = video.videoWidth;
    canvas.height = video.videoHeight;
    canvas.getContext("2d").drawImage(video, 0, 0, canvas.width, canvas.height);
    const imageBase64 = canvas.toDataURL("image/jpeg", 0.9).split(",")[1];

    // Scale coordinates from display size to actual video resolution
    const scaleX = video.videoWidth / video.clientWidth;
    const scaleY = video.videoHeight / video.clientHeight;
    const clickX = Math.round(evt.clientX - video.getBoundingClientRect().left);
    const clickY = Math.round(evt.clientY - video.getBoundingClientRect().top);
    const actualX = Math.round(clickX * scaleX);
    const actualY = Math.round(clickY * scaleY);
    console.log(`Scaled coordinates: (${actualX}, ${actualY}) from video ${video.videoWidth}x${video.videoHeight}`);

    // Crop area around click using scaled coordinates
    const crop = document.createElement("canvas");
    crop.width = crop.height = 100;
    crop
      .getContext("2d")
      .drawImage(canvas, actualX - 50, actualY - 50, 100, 100, 0, 0, 100, 100);
    const croppedBase64 = crop.toDataURL("image/jpeg", 0.9).split(",")[1];

    // Show loading state
    const loadingMsg = createLoadingMessage("🔍 Analyzing...");
    snapshots.appendChild(loadingMsg);

    try {
      const res = await fetch("/image-molecules", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          imageBase64,
          croppedImageBase64: croppedBase64,
          x: actualX,
          y: actualY,
        }),
      });

      loadingMsg.remove();

      if (!res.ok) throw new Error(`HTTP ${res.status}`);

      const { output } = await res.json();
      
      // Display result using consistent method
      const objectName = output.object || "Unknown object";
      processAnalysisResult(output, snapshots, "🔍", objectName);

    } catch (err) {
      loadingMsg.remove();
      
      const errorMsg = document.createElement("h3");
      errorMsg.textContent = `⚠ Error: ${err.message}`;
      errorMsg.style.color = "red";
      snapshots.appendChild(errorMsg);
    }
  }

  /* ---------- 3D Molecule Rendering ---------- */
  
  // Get descriptive name for chemical string (SMILES or mineral formula)
  function getMoleculeName(chemical) {
    const moleculeNames = {
      // Simple molecules
      'O': 'Water',
      'CCO': 'Ethanol',
      'CC(=O)O': 'Acetic Acid',
      'C': 'Methane',
      'CO': 'Methanol',
      'C(CO)N': 'Ethanolamine',
      'C(C(=O)O)N': 'Glycine',
      'C(CC(=O)O)N': 'GABA',
      'C(CC(=O)O)C(=O)O': 'Succinic Acid',
      'C(C(=O)O)O': 'Glycolic Acid',
      
      // Aromatic compounds
      'C1=CC=CC=C1': 'Benzene',
      'CC1=CC=CC=C1': 'Toluene',
      'C1=CC=C(C=C1)O': 'Phenol',
      'C1=CC=C(C=C1)N': 'Aniline',
      
      // Common organics
      'CCCCCCCC=O': 'Octanal',
      'CCC(C)C(=O)OC(C)C': 'Isopropyl Isovalerate',
      'CCCOC(=O)N1CCCC1C(=O)OCCC': 'Polyurethane Monomer',
      'C1CCOC1': 'Tetrahydrofuran',
      'NCCCN': 'Propanediamine',
      
      // Complex molecules
      'CN1C=NC2=C1C(=O)N(C(=O)N2C)C': 'Caffeine',
      'CC(=O)OC1=CC=CC=C1C(=O)O': 'Aspirin',
      'C1NC(=O)N(C)C2=CC=CC=C12': 'N-Methylbenzimidazolone',
      'NC1=NC(=NC(=N1)Cl)Cl': 'Cyanuric Chloride',
      
      // Pharmaceuticals/Complex
      'C[C@H](NC(=O)[C@@H](N)Cc1c[nH]c2ccc(Cl)cc12)C(=O)O': 'Chlorotryptophan Derivative',
      'O=C(Nc1ccc(cc1)C(=O)O)C(F)(F)F': 'Trifluoroacetyl-p-aminobenzoic Acid',
      'C(C(C(=O)NC(C(=O)O)C(C)O)O)O': 'Threonine Derivative',
      
      // Minerals and Crystal Structures
      'CaCO3': 'Calcite (Calcium Carbonate)',
      'CaCO₃': 'Calcite (Calcium Carbonate)',
      'SiO2': 'Quartz (Silicon Dioxide)',
      'SiO₂': 'Quartz (Silicon Dioxide)',
      'Al2O3': 'Corundum (Aluminum Oxide)',
      'Al₂O₃': 'Corundum (Aluminum Oxide)',
      'FeS2': 'Pyrite (Iron Disulfide)',
      'FeS₂': 'Pyrite (Iron Disulfide)',
      'NaCl': 'Halite (Sodium Chloride)',
      'quartz': 'Quartz (Silicon Dioxide)',
      'calcite': 'Calcite (Calcium Carbonate)',
      'corundum': 'Corundum (Aluminum Oxide)',
      'pyrite': 'Pyrite (Iron Disulfide)',
      'halite': 'Halite (Sodium Chloride)',
      'salt': 'Halite (Sodium Chloride)'
    };
    
    return moleculeNames[chemical] || `Structure (${chemical.substring(0, 20)}${chemical.length > 20 ? '...' : ''})`;
  }

  async function generateSDFs(smiles, objectName) {
    if (!smiles || smiles.length === 0) return;
    
    console.log(`🧬 Generating SDFs for ${smiles.length} chemicals from ${objectName}...`);
    
    try {
      const response = await fetch('/generate-sdfs', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles, overwrite: false })
      });

      if (!response.ok) {
        throw new Error(`HTTP ${response.status}: ${response.statusText}`);
      }

      const data = await response.json();
      console.log(`✅ SDF Response:`, data);
      
      // Create column with all chemicals found
      const summary = {
        total: smiles.length,
        visualizable: data.sdfPaths ? data.sdfPaths.length : 0,
        skipped: data.skipped ? data.skipped.length : 0,
        errors: data.errors ? data.errors.length : 0
      };
      
      createObjectColumn(objectName, data.sdfPaths || [], smiles, null, summary, data.skipped || []);
      
    } catch (error) {
      console.error("❌ SDF generation error:", error);
      
      // Create an error column
      createObjectColumn(objectName, [], smiles, '🧬 Working on 3D structures...');
    }
  }

  async function createObjectColumn(objectName, sdfFiles, smiles = [], errorMessage = null, summary = null, skippedChemicals = []) {
    console.log(`🧬 Creating object column for "${objectName}" with ${sdfFiles.length} molecules`);
    
    const gldiv = document.getElementById('gldiv');
    
    // Create object column container
    const objectColumn = document.createElement("div");
    objectColumn.className = "object-column";
    
    // Create title container with close button
    const titleContainer = document.createElement("div");
    titleContainer.className = "object-title";
    
    const titleText = document.createElement("span");
    titleText.textContent = objectName;
    titleContainer.appendChild(titleText);
    
    const closeButton = document.createElement("button");
    closeButton.className = "column-close";
    closeButton.textContent = "×";
    closeButton.onclick = () => objectColumn.remove();
    titleContainer.appendChild(closeButton);
    
    objectColumn.appendChild(titleContainer);
    
    // Add summary information if available
    if (summary) {
      const summaryDiv = document.createElement("div");
      summaryDiv.className = "chemical-summary";
      summaryDiv.innerHTML = `
        <div>Total chemicals found: ${summary.total}</div>
        <div>3D visualizable: ${summary.visualizable}</div>
        ${summary.skipped > 0 ? `<div>Non-SMILES formats: ${summary.skipped}</div>` : ''}
        ${summary.errors > 0 ? `<div>Failed: ${summary.errors}</div>` : ''}
      `;
      objectColumn.appendChild(summaryDiv);
    }
    
    // Show skipped chemicals if any
    if (skippedChemicals.length > 0) {
      const skippedDiv = document.createElement("div");
      skippedDiv.className = "skipped-chemicals";
      skippedDiv.innerHTML = `
        <div class="skipped-title">🧪 Other chemicals found:</div>
        <div class="skipped-list">${skippedChemicals.join(', ')}</div>
        <div class="skipped-note">* These likely represent minerals/crystals that can't be shown in 3D molecular view</div>
      `;
      objectColumn.appendChild(skippedDiv);
    }
    
    if (errorMessage) {
      // Show error message
      const errorDiv = document.createElement("div");
      errorDiv.textContent = errorMessage;
      errorDiv.style.color = '#ffffff';
      errorDiv.style.textAlign = 'center';
      errorDiv.style.padding = '20px';
      objectColumn.appendChild(errorDiv);
    } else {
      // Load molecules vertically in this column
      const viewers = [];
      for (let i = 0; i < sdfFiles.length; i++) {
        const sdfFile = sdfFiles[i];
        const moleculeSmiles = smiles[i] || '';
        
        // Create molecule container with name
        const moleculeContainer = document.createElement("div");
        moleculeContainer.className = "molecule-container";
        
        // Add molecule name
        const moleculeName = document.createElement("div");
        moleculeName.className = "molecule-name";
        moleculeName.textContent = getMoleculeName(moleculeSmiles);
        moleculeContainer.appendChild(moleculeName);
        
        // Add 3D viewer container
        const container = document.createElement("div");
        container.className = "mol-viewer-container";
        moleculeContainer.appendChild(container);
        
        objectColumn.appendChild(moleculeContainer);
        
        const viewer = await render(sdfFile, container);
        if (viewer) viewers.push(viewer);
      }
      
      console.log(`✅ Successfully loaded ${viewers.length} 3D molecular viewers for ${objectName}`);
      viewers.forEach(viewer => { viewer.resize(); viewer.render(); });
    }
    
    // Add the complete column to the main container
    gldiv.appendChild(objectColumn);
  }

  async function render(sdfFile, container) {
    try {
      console.log(`🧬 Loading ${sdfFile}`);
      
      const response = await fetch(sdfFile);
      if (!response.ok) {
        throw new Error(`HTTP error ${response.status}`);
      }
      
      const sdfData = await response.text();
      if (!sdfData.trim()) {
        throw new Error(`Empty SDF data`);
      }
      if (!sdfData.includes("$$$$")) {
        throw new Error(`Invalid SDF format`);
      }
      
      const viewer = $3Dmol.createViewer(container);
      viewer.addModel(sdfData, "sdf");
      
      viewer.setBackgroundColor('#000000');
      
      // Use van der Waals radii for accurate representation
      viewer.setStyle({}, { 
        sphere: { 
          scale: 0.8  // Scale factor for van der Waals radii (Bohr radius)
        }
      });
      
      viewer.zoomTo();
      viewer.render();
      
      console.log(`✅ Loaded molecule successfully`);
      return viewer;
    } catch (error) {
      console.error(`❌ Failed to load molecule:`, error);
      container.textContent = `❌ Error loading molecule`;
      container.style.color = 'red';
      container.style.textAlign = 'center';
      container.style.padding = '20px';
      return null;
    }
  }

  /* ---------- Initialize ---------- */
  (async () => {
    console.log("🚀 Starting app initialization...");
    
    try {
      await startCamera();
      await setupSwitchCamera();
      console.log("✅ Camera setup complete");
    } catch (err) {
      console.error("❌ Camera setup failed:", err);
    }
  })();
});
