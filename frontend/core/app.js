// app.js - Core molecular analysis application
import { paymentManager } from '../components/payment.js';
import { cameraManager } from '../components/camera.js';
import { cameraHandler } from '../components/camera-handler.js';
import { uiManager } from '../components/ui-utils.js';

class MolecularApp {
  constructor() {
    this.snapshots = null;
    this.objectInput = null;
    this.viewers = [];
  }

  // Initialize the application
  async initialize() {
    // Get DOM elements
    this.snapshots = document.querySelector(".snapshots-container");
    this.objectInput = document.getElementById("object-input");

    // Initialize modules
    uiManager.initialize();
    uiManager.setupDebuggingFunctions();
    uiManager.clearDevelopmentStates();
    uiManager.showMainApp();

    // Setup main event listeners
    this.setupEventListeners();

    // Initialize payment system
    await paymentManager.checkInitialPaymentSetup();
      
    // Initialize camera system
    await cameraManager.initialize();

    // Safari-specific camera permission request
    if (cameraManager.isSafari) {
      console.log('🍎 Safari detected - requesting camera permission early');
    setTimeout(() => {
        cameraManager.requestPermission();
      }, 1000);
  }
  
    console.log('✅ Molecular analysis app initialized');
  }

  // Setup main application event listeners
  setupEventListeners() {
    // Setup text analysis with debugging support
    this.setupTextAnalysis();

    // Setup camera handler event listeners
    cameraHandler.setupEventListeners();
      
    // Mode switching based on user interaction
    const video = document.getElementById("video-feed");
    if (video) {
      video.addEventListener("click", () => uiManager.switchToCameraMode());
      video.addEventListener("touchstart", () => uiManager.switchToCameraMode());
  }
  
    // Photo mode switching
    const photoUpload = document.getElementById("photo-upload");
    const photoUrl = document.getElementById("photo-url");
    const urlAnalyze = document.getElementById("url-analyze");
    
    if (photoUpload) {
      photoUpload.addEventListener("change", () => uiManager.switchToPhotoMode());
    }
    if (photoUrl) {
      photoUrl.addEventListener("focus", () => uiManager.switchToPhotoMode());
    }
    if (urlAnalyze) {
      urlAnalyze.addEventListener("click", () => uiManager.switchToPhotoMode());
    }
  
    // Text input clears mode selection
    this.objectInput.addEventListener("focus", () => uiManager.clearModeSelection());

    // Listen for image analysis completion
    document.addEventListener('imageAnalysisComplete', (e) => {
      const { output, icon, objectName, useQuotes, croppedImageData } = e.detail;
      this.processAnalysisResult(output, icon, objectName, useQuotes, croppedImageData);
    });
  }
  
  // Handle Enter key press for text analysis
  setupTextAnalysis() {
    this.objectInput.addEventListener("keyup", async (e) => {
      if (e.key !== "Enter") return;
      
      // 🔴 BREAKPOINT: Set breakpoint here to debug text analysis trigger
      console.log('🚀 Text analysis triggered from Enter key');
      console.log('📊 App state before analysis:', {
        isProcessing: this.isProcessing,
        hasPaymentSetup: this.hasPaymentSetup,
        inputValue: this.objectInput.value,
        paymentVisible: this.paymentPopdown.style.display !== 'none'
      });
      
      await this.handleTextAnalysis();
    });
  }

  // Main text analysis handler
  async handleTextAnalysis() {
    // 🔴 BREAKPOINT: Set breakpoint here to debug main analysis flow
    console.log('🔬 Starting handleTextAnalysis');
    console.log('📊 Current state:', {
      isProcessing: this.isProcessing,
      hasPaymentSetup: this.hasPaymentSetup,
      inputValue: this.objectInput.value
    });
    
    if (this.isProcessing) {
      console.log('⚠️ Already processing, skipping analysis');
      return;
    }

    const inputValue = this.objectInput.value.trim();
    console.log('📝 Input value:', inputValue);
    
    if (!inputValue) {
      console.log('❌ No input value, skipping analysis');
      return;
    }

    // 🔴 BREAKPOINT: Set breakpoint here to debug payment check
    if (!this.hasPaymentSetup) {
      console.log('💳 Payment not set up, showing simple message');
      
      // Create a simple message instead of analysis
      const messageColumn = uiManager.createColumn("See payment setup above", "payment-required");
      messageColumn.innerHTML = `
        <div class="molecule-container">
          <div class="molecule-info">
            <h3>Payment Required</h3>
            <p>See payment setup above</p>
            <div class="analysis-note">Complete payment setup to analyze molecules</div>
          </div>
        </div>
      `;
      
      // Clear the input
      this.objectInput.value = "";
      return;
    }

    // 🔴 BREAKPOINT: Set breakpoint here to debug processing start
    this.isProcessing = true;
    this.currentAnalysisType = 'text';
    console.log('🏁 Starting processing with type:', this.currentAnalysisType);
    
    try {
      this.hidePaymentPopdown();
      this.showProcessing();
      
      // 🔴 BREAKPOINT: Set breakpoint here to debug API call preparation
      console.log('🌐 Preparing API call for text analysis');
      const response = await fetch("/analyze-text", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ text: inputValue }),
      });

      // 🔴 BREAKPOINT: Set breakpoint here to debug API response
      console.log('📡 API response received:', {
        status: response.status,
        statusText: response.statusText,
        ok: response.ok
      });

      if (!response.ok) {
        const errorText = await response.text();
        console.error('❌ API error response:', errorText);
        throw new Error(`Failed to analyze text: ${response.status} ${errorText}`);
      }

      const result = await response.json();
      // 🔴 BREAKPOINT: Set breakpoint here to debug API result processing
      console.log('📋 API result:', result);
      
      this.lastAnalysis = result;
      this.displayResults(result);
      
      // Clear the input
      this.objectInput.value = "";
      
    } catch (error) {
      // 🔴 BREAKPOINT: Set breakpoint here to debug errors
      console.error('💥 Error in handleTextAnalysis:', error);
      this.handleError(error);
    } finally {
      // 🔴 BREAKPOINT: Set breakpoint here to debug cleanup
      console.log('🧹 Cleaning up processing state');
      this.isProcessing = false;
      this.hideProcessing();
    }
  }



  // Process analysis results and display molecules
  processAnalysisResult(output, icon, objectName, useQuotes = false, croppedImageData = null) {
    const chemicals = output.chemicals || [];

    // Handle description responses
    if (chemicals.length === 1 && chemicals[0].smiles && chemicals[0].smiles.startsWith("DESCRIPTION: ")) {
      const description = chemicals[0].smiles.replace("DESCRIPTION: ", "");
      this.generateSDFs([], objectName, description, null, croppedImageData);
      return;
    }

    // Handle molecular responses
    const smiles = chemicals.map((chem) => chem.smiles).filter(Boolean);
    if (smiles.length > 0) {
      this.generateSDFs(smiles, objectName, null, chemicals, croppedImageData);
    }
  }

  // Generate SDF files and create 3D visualizations
  async generateSDFs(smiles, objectName, description = null, chemicals = null, croppedImageData = null) {
    if (smiles.length === 0 && !description) {
      this.createClosableErrorMessage("No valid molecules found for visualization");
    return;
  }

    try {
      let sdfPaths = [];
      
      if (smiles.length > 0) {
        const response = await fetch("/generate-sdfs", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ smiles, overwrite: true }),
      });

        if (!response.ok) throw new Error(`SDF generation failed: ${response.status}`);
        const result = await response.json();
        sdfPaths = result.sdfPaths || [];
      }

      await this.createObjectColumn(
    objectName,
        sdfPaths,
        smiles,
        null,
        null,
        [],
        description,
        chemicals,
        croppedImageData
      );

    } catch (error) {
      console.error("SDF generation error:", error);
      this.createClosableErrorMessage(`Error generating 3D models: ${error.message}`);
    }
  }

  // Create object column with 3D molecular visualizations
  async createObjectColumn(objectName, sdfFiles, smiles = [], errorMessage = null, 
                          summary = null, skippedChemicals = [], description = null, 
                          chemicals = null, croppedImageData = null) {
    
    const gldiv = document.getElementById("gldiv");
    const objectColumn = document.createElement("div");
    objectColumn.className = "object-column";

    // Create header
    const header = document.createElement("div");
    header.className = "object-header";

    const titleContainer = document.createElement("div");
    titleContainer.className = "object-title-container";

    const icon = document.createElement("div");
    icon.className = "object-icon";
    icon.textContent = description ? "📖" : "🔬";
    
    const name = document.createElement("div");
    name.className = "object-name";
    name.textContent = objectName;

    const closeBtn = document.createElement("button");
    closeBtn.className = "object-close";
    closeBtn.innerHTML = '<img src="close.svg" alt="Close" width="16" height="16" />';
    closeBtn.onclick = () => {
      objectColumn.remove();
      this.updateScrollHandles();
    };

    titleContainer.appendChild(icon);
    titleContainer.appendChild(name);
    header.appendChild(titleContainer);
    header.appendChild(closeBtn);
    objectColumn.appendChild(header);

    // Add cropped image if provided
    if (croppedImageData) {
      const imageContainer = document.createElement("div");
      imageContainer.className = "cropped-image-container";

      const img = document.createElement("img");
      img.src = `data:image/jpeg;base64,${croppedImageData}`;
      img.className = "image-highlighted";
      img.alt = "Analysis region";

      imageContainer.appendChild(img);
      objectColumn.appendChild(imageContainer);
    }

    // Handle description vs molecules
    if (description) {
      const descDiv = document.createElement("div");
      descDiv.className = "description-content";
      descDiv.textContent = description;
      objectColumn.appendChild(descDiv);
    } else {
      // Render 3D molecules
      for (let i = 0; i < sdfFiles.length; i++) {
        const sdfFile = sdfFiles[i];
        const chemical = chemicals?.[i] || { smiles: smiles[i] };
        
        const container = document.createElement("div");
        container.className = "molecule-viewer";
        objectColumn.appendChild(container);

        const moleculeContainer = document.createElement("div");
        moleculeContainer.className = "molecule-container";

        const moleculeName = document.createElement("div");
        moleculeName.className = "molecule-name";

        const displayName = uiManager.getMoleculeName(chemical);
        moleculeName.textContent = displayName;

        // Add Wikipedia link
        const wikipediaLink = document.createElement("a");
        wikipediaLink.textContent = displayName;
        wikipediaLink.href = `https://en.wikipedia.org/wiki/${encodeURIComponent(displayName)}`;
        wikipediaLink.target = "_blank";
        wikipediaLink.rel = "noopener noreferrer";
        wikipediaLink.className = "wikipedia-link";

        moleculeName.appendChild(wikipediaLink);
        moleculeContainer.appendChild(moleculeName);
        objectColumn.appendChild(moleculeContainer);

        const viewer = await this.render(sdfFile, container);
        if (viewer) this.viewers.push(viewer);
      }

      this.viewers.forEach((viewer) => {
        viewer.resize();
        viewer.render();
      });
    }

    gldiv.appendChild(objectColumn);
    this.updateScrollHandles();
  }

  // Render 3D molecule visualization
  async render(sdfFile, container) {
    try {
      const urlParts = sdfFile.split("/");
      const filename = urlParts.pop();
      const encodedFilename = encodeURIComponent(filename);
      const encodedPath = urlParts.join("/") + "/" + encodedFilename;

      const response = await fetch(encodedPath);
      if (!response.ok) {
        throw new Error(`HTTP error ${response.status}`);
      }

      const sdfData = await response.text();
      const viewer = $3Dmol.createViewer(container, {
        defaultcolors: $3Dmol.rasmolElementColors,
      });

      viewer.addModel(sdfData, "sdf");
      viewer.setStyle({}, { sphere: { scale: 0.8 } });
      viewer.zoomTo();
      viewer.render();

      return viewer;
      
    } catch (error) {
      console.error(`Failed to load molecule:`, error);
      container.textContent = `Error loading molecule: ${error.message}`;
      container.className += " error-text";
      return null;
    }
  }

  // Create error message
  createClosableErrorMessage(message) {
    const errorDiv = uiManager.createErrorMessage(message, this.snapshots);
    this.updateScrollHandles();
    return errorDiv;
  }

  // Update scroll handles (placeholder)
  updateScrollHandles() {
    if (window.updateScrollHandles) {
      window.updateScrollHandles();
    }
  }

  // Cleanup resources
  cleanup() {
    cameraManager.cleanup();
    uiManager.cleanup();
    
    // Clean up 3D viewers
    this.viewers.forEach(viewer => {
      if (viewer && viewer.clear) {
        viewer.clear();
      }
    });
    this.viewers = [];
    }
  }

// Initialize app when DOM is ready
document.addEventListener("DOMContentLoaded", async () => {
  const app = new MolecularApp();
  await app.initialize();

  // Make app globally available for debugging
  window.molecularApp = app;
});
