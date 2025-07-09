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
      
      const p = document.createElement("p");
      p.textContent = `⚠ Error analyzing "${object}": ${err.message}`;
      p.style.color = "red";
      snapshots.appendChild(p);
    }
    objectInput.value = "";
  });

  // Photo upload functionality
  photoUpload.addEventListener("change", async (e) => {
    const file = e.target.files[0];
    if (!file) return;

    // Validate file type
    if (!file.type.startsWith('image/')) {
      alert('Please select an image file');
      return;
    }

    // Show loading state
    const loadingMsg = createLoadingMessage(`🔍 Analyzing uploaded photo...`);
    snapshots.appendChild(loadingMsg);

    try {
      // Convert to base64
      const imageBase64 = await fileToBase64(file);
      
      // Use center of image as click point (50%, 50%)
      const mockClickX = 0.5;
      const mockClickY = 0.5;

      const res = await fetch("/image-molecules", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          imageBase64,
          croppedImageBase64: imageBase64, // Use full image as crop
          x: mockClickX,
          y: mockClickY,
        }),
      });

      loadingMsg.remove();

      if (!res.ok) throw new Error(`HTTP ${res.status}`);

      const { output } = await res.json();
      
      // Display result
      const objectName = output.object || "Uploaded image";
      processAnalysisResult(output, snapshots, "📁", objectName);

      // Clear URL input
      photoUrl.value = '';

    } catch (err) {
      loadingMsg.remove();
      
      const errorMsg = document.createElement("p");
      errorMsg.textContent = `⚠ Error analyzing photo: ${err.message}`;
      errorMsg.style.color = "red";
      snapshots.appendChild(errorMsg);
    }

    // Reset file input
    photoUpload.value = '';
  });

  // Photo URL functionality
  async function analyzeImageFromUrl(url) {
    try {
      // Show loading state
      const loadingMsg = createLoadingMessage(`🔍 Analyzing image from URL...`);
      snapshots.appendChild(loadingMsg);

      // Fetch and convert image to base64
      const imageBase64 = await urlToBase64(url);
      
      // Use center of image as click point
      const mockClickX = 0.5;
      const mockClickY = 0.5;

      const res = await fetch("/image-molecules", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          imageBase64,
          croppedImageBase64: imageBase64,
          x: mockClickX,
          y: mockClickY,
        }),
      });

      loadingMsg.remove();

      if (!res.ok) throw new Error(`HTTP ${res.status}`);

      const { output } = await res.json();
      
      // Display result
      const objectName = output.object || "Image from URL";
      processAnalysisResult(output, snapshots, "🔗", objectName);

      // Clear URL input
      photoUrl.value = '';

    } catch (err) {
      // Remove loading message if it exists
      const loadingMsgs = snapshots.querySelectorAll('p');
      loadingMsgs.forEach(msg => {
        if (msg.textContent.includes('Analyzing image from URL')) {
          msg.remove();
        }
      });
      
      const errorMsg = document.createElement("p");
      errorMsg.textContent = `⚠ Error loading image from URL: ${err.message}`;
      errorMsg.style.color = "red";
      snapshots.appendChild(errorMsg);
    }
  }

  // URL input event listeners
  urlAnalyze.addEventListener("click", () => {
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
  });

  photoUrl.addEventListener("keyup", (e) => {
    if (e.key === "Enter") {
      urlAnalyze.click();
    }
  });

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
    const loadingMsg = document.createElement("p");
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
    const result = document.createElement("p");
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
      
      const errorMsg = document.createElement("p");
      errorMsg.textContent = `⚠ Error: ${err.message}`;
      errorMsg.style.color = "red";
      snapshots.appendChild(errorMsg);
    }
  }

  /* ---------- 3D Molecule Rendering ---------- */
  
  // Get descriptive name for SMILES string
  function getMoleculeName(smiles) {
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
      'C(C(C(=O)NC(C(=O)O)C(C)O)O)O': 'Threonine Derivative'
    };
    
    return moleculeNames[smiles] || `Molecule (${smiles.substring(0, 20)}${smiles.length > 20 ? '...' : ''})`;
  }

  async function generateSDFs(smiles, objectName) {
    if (!smiles || smiles.length === 0) return;
    
    console.log(`🧬 Generating SDFs for ${smiles.length} molecules from ${objectName}...`);
    
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
      if (data.sdfPaths) {
        console.log(`✅ Loading 3D molecules for ${objectName} with ${data.sdfPaths.length} structures`);
        createObjectColumn(objectName, data.sdfPaths, smiles);
      } else {
        console.error("❌ No SDF paths in response:", data);
      }
    } catch (error) {
      console.error("❌ SDF generation error:", error);
      
      // Create an error column
      createObjectColumn(objectName, [], smiles, '🧬 Working on 3D structures...');
    }
  }

  async function createObjectColumn(objectName, sdfFiles, smiles = [], errorMessage = null) {
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
      
      // Use ball and stick representation for better visibility
      viewer.setStyle({}, {
        stick: { radius: 0.15, colorscheme: 'default' },
        sphere: { scale: 0.25, colorscheme: 'default' }
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
