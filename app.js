document.addEventListener("DOMContentLoaded", () => {
  /* ---------- DOM shortcuts ---------- */
  const video = document.getElementById("video-feed");
  const snapshots = document.querySelector(".snapshots-container");
  const msgBox = document.querySelector(".permission-message");
  const objectInput = document.getElementById("object-input");
  const appContainer = document.querySelector(".app-container");
  const instructionText = document.querySelector(".instruction-text");
  
  // Get mode toggle elements
  const imageMode = document.getElementById("image-mode");
  const textMode = document.getElementById("text-mode");

  // Handle mode switching
  function updateMode() {
    const body = document.body;
    if (textMode.checked) {
      body.classList.add("text-mode");
      body.classList.remove("image-mode");
      objectInput.focus();
    } else {
      body.classList.remove("text-mode");
      body.classList.add("image-mode");
    }
  }

  imageMode.addEventListener("change", updateMode);
  textMode.addEventListener("change", updateMode);
  
  // Initialize mode
  updateMode();

  // text input handling
  objectInput.addEventListener("keyup", async (e) => {
    if (e.key !== "Enter" || !textMode.checked) return;
    const object = objectInput.value.trim();
    if (!object) return;

    // Show loading state
    const loadingMsg = document.createElement("p");
    loadingMsg.textContent = `Analyzing "${object}"...`;
    loadingMsg.style.fontStyle = "italic";
    loadingMsg.style.opacity = "0.7";
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

      // Simple result display
      const smilesCount = output.smiles ? output.smiles.length : 0;
      const result = document.createElement("p");
      result.textContent = `📝 ${object} → ${smilesCount} structure${smilesCount !== 1 ? 's' : ''} found`;
      snapshots.appendChild(result);

      if (output.smiles && output.smiles.length > 0) {
        document.getElementById('gldiv').innerHTML = '';
        generateSDFs(output.smiles);
      }
    } catch (err) {
      loadingMsg.remove();
      
      const p = document.createElement("p");
      p.textContent = `⚠ Error analyzing "${object}": ${err.message}`;
      p.style.color = "red";
      snapshots.appendChild(p);
    }
    objectInput.value = "";
  });

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
    if (!imageMode.checked) return;
    
    const rect = video.getBoundingClientRect();
    const clickX = Math.round(evt.clientX - rect.left);
    const clickY = Math.round(evt.clientY - rect.top);
    console.log(`Click at display (${clickX}, ${clickY})`);

    // Simple click feedback
    const mark = document.createElement("div");
    mark.className = "feedback-box";
    mark.style.left = `${clickX - 30}px`;
    mark.style.top = `${clickY - 30}px`;
    video.appendChild(mark);

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

    // Simple snapshot display
    const snapshot = document.createElement("p");
    snapshot.textContent = "🔍 Analyzing...";
    snapshots.appendChild(snapshot);

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

      if (!res.ok) throw new Error(`HTTP ${res.status}`);

      const { output } = await res.json();
      
      // Simple image result display
      const objectName = output.object || "Unknown object";
      const smilesCount = output.smiles ? output.smiles.length : 0;
      
      snapshot.textContent = `🔍 ${objectName} → ${smilesCount} structure${smilesCount !== 1 ? 's' : ''} found`;

      if (output.smiles && output.smiles.length > 0) {
        document.getElementById('gldiv').innerHTML = '';
        generateSDFs(output.smiles);
      }

    } catch (err) {
      snapshot.textContent = `⚠ Error: ${err.message}`;
      snapshot.style.color = "red";
    }

    // Auto-remove feedback after delay
    setTimeout(() => {
      if (mark.parentElement) {
        mark.style.opacity = "0.3";
      }
    }, 2000);
  }

  /* ---------- 3D Molecule Rendering ---------- */
  async function generateSDFs(smiles) {
    if (!smiles || smiles.length === 0) return;
    
    console.log(`🧬 Generating SDFs for ${smiles.length} molecules...`);
    
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
        console.log(`✅ Loading 3D grid with ${data.sdfPaths.length} molecules`);
        load3DmolGrid(data.sdfPaths);
      } else {
        console.error("❌ No SDF paths in response:", data);
      }
    } catch (error) {
      console.error("❌ SDF generation error:", error);
      
      const gldiv = document.getElementById('gldiv');
      gldiv.textContent = '🧬 Working on 3D structures...';
    }
  }

  async function load3DmolGrid(sdfFiles) {
    console.log(`🧬 Loading 3D grid with ${sdfFiles.length} molecules`);
    
    const gldiv = document.getElementById('gldiv');
    gldiv.innerHTML = '';
    
    const viewers = [];
    for (const sdfFile of sdfFiles) {
      const container = document.createElement("div");
      container.className = "mol-viewer-container";
      gldiv.appendChild(container);
      const viewer = await render(sdfFile, container);
      if (viewer) viewers.push(viewer);
    }

    console.log(`✅ Successfully loaded ${viewers.length} 3D molecular viewers`);
    viewers.forEach(viewer => { viewer.resize(); viewer.render(); });
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
      
      // Set background color for the 3D viewer
      // Options: 'white', 'black', 'lightgray', 'darkgray', or hex colors like '#f8f9fa'
      viewer.setBackgroundColor('#000000');  // Black background to match app theme
      
      // Use van der Waals radii for most accurate representation
      viewer.setStyle({}, { 
        sphere: { 
          scale: 0.8  // Scale factor for van der Waals radii (0.8 gives good visual balance)
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
