// Simple App Logic - Following ui.mdc Guidelines
import { simplePaymentManager } from '../components/simple-payment.js';
import { cameraManager } from '../components/camera.js';
import { cameraHandler } from '../components/camera-handler.js';
import { uiManager } from '../components/ui-utils.js';

class MolecularApp {
  constructor() {
    this.snapshots = null;
    this.objectInput = null;
    this.viewers = [];
    this.isProcessing = false;
    this.hasPaymentSetup = false;
    this.currentAnalysisType = null;
    this.lastAnalysis = null;
  }

  async initialize() {
    this.snapshots = document.querySelector(".snapshots-area");
    this.objectInput = document.getElementById("object-input");

    uiManager.initialize();
    uiManager.setupDebuggingFunctions();
    uiManager.showMainApp();

    this.setupEventListeners();

    // Simple payment check
    simplePaymentManager.checkPaymentRequired();
    this.hasPaymentSetup = true; // Keep simple for now
    
    await cameraManager.initialize();

    if (cameraManager.isSafari && !cameraManager.hasStoredCameraPermission()) {
      setTimeout(() => {
        cameraManager.requestPermission();
      }, 1000);
    }
  }

  setupEventListeners() {
    this.setupTextAnalysis();
    cameraHandler.setupEventListeners();
    
    // Keyboard shortcut (Cmd+K / Ctrl+K)
    document.addEventListener('keydown', (event) => {
      if ((event.metaKey || event.ctrlKey) && event.key === 'k') {
        event.preventDefault();
        const textInput = document.getElementById('object-input');
        if (textInput) {
          textInput.focus();
          textInput.select();
        }
      }
    });
    
    // Platform-specific placeholder
    const textInput = document.getElementById('object-input');
    if (textInput) {
      const isMac = navigator.platform.toUpperCase().indexOf('MAC') >= 0;
      const shortcutKey = isMac ? '⌘K' : 'Ctrl+K';
      textInput.placeholder = `Type any molecule name (e.g., caffeine, aspirin, water)... (${shortcutKey} to focus)`;
    }

    document.addEventListener('imageAnalysisComplete', (e) => {
      const { output, icon, objectName, useQuotes, croppedImageData } = e.detail;
      this.processAnalysisResult(output, icon, objectName, useQuotes, croppedImageData);
    });
  }
   
  setupTextAnalysis() {
    this.objectInput.addEventListener("keyup", async (e) => {
      if (e.key !== "Enter") return;
      await this.handleTextAnalysis();
    });
  }

  async handleTextAnalysis() {
    if (this.isProcessing) return;

    const inputValue = this.objectInput.value.trim();
    if (!inputValue) return;

    this.isProcessing = true;
    this.currentAnalysisType = 'text';
    
    try {
      this.showProcessing();
      
      const response = await fetch("/analyze-text", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ text: inputValue }),
      });

      if (!response.ok) {
        const errorText = await response.text();
        throw new Error(`Failed to analyze text: ${response.status} ${errorText}`);
      }

      const result = await response.json();
      this.lastAnalysis = result;
      
      const objectName = this.objectInput.value.trim();
      this.processAnalysisResult(result, null, objectName, false, null);
      
      this.objectInput.value = "";
      
    } catch (error) {
      console.error('Analysis error:', error);
      this.handleError(error);
    } finally {
      this.isProcessing = false;
      this.hideProcessing();
    }
  }

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
    closeBtn.innerHTML = '<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M18 6L6 18M6 6l12 12"/></svg>';
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

  async render(sdfFile, container) {
    try {
      const response = await fetch(sdfFile);
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

  showProcessing() {
    // Simple processing state
    const objectInput = document.getElementById('object-input');
    if (objectInput) {
      objectInput.disabled = true;
      objectInput.placeholder = 'Analyzing...';
    }
  }

  hideProcessing() {
    const objectInput = document.getElementById('object-input');
    if (objectInput) {
      objectInput.disabled = false;
      const isMac = navigator.platform.toUpperCase().indexOf('MAC') >= 0;
      const shortcutKey = isMac ? '⌘K' : 'Ctrl+K';
      objectInput.placeholder = `Type any molecule name (e.g., caffeine, aspirin, water)... (${shortcutKey} to focus)`;
    }
  }

  handleError(error) {
    console.error('Error in molecular analysis:', error);
    this.createClosableErrorMessage(`Analysis failed: ${error.message}`);
  }

  createClosableErrorMessage(message) {
    const errorDiv = uiManager.createErrorMessage(message, this.snapshots);
    this.updateScrollHandles();
    return errorDiv;
  }

  updateScrollHandles() {
    if (window.updateScrollHandles) {
      window.updateScrollHandles();
    }
  }

  cleanup() {
    cameraManager.cleanup();
    uiManager.cleanup();
    
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
  window.app = app; // Make app globally available for debugging
  await app.initialize();

  window.molecularApp = app;
});
