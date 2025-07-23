// Simple App Logic - Following ui.mdc Guidelines
import { simplePaymentManager } from '../components/simple-payment.js';
import { cameraManager } from '../components/camera.js';
import { cameraHandler } from '../components/camera-handler.js';
import { uiManager } from '../components/ui-utils.js';
import { logger } from '../components/logger.js';

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
    this.snapshots = document.querySelector(".snapshots-container");
    this.objectInput = document.getElementById("object-input");

    uiManager.initialize();
    uiManager.setupDebuggingFunctions();
    uiManager.showMainApp();

    this.setupEventListeners();

    // Initialize sidebar payment system
    await simplePaymentManager.checkPaymentRequired();
    
    // Auto-enable dev mode for localhost
    if (location.hostname === 'localhost' || location.hostname === '127.0.0.1') {
      logger.info('Auto-enabling developer mode for localhost');
      this.hasPaymentSetup = true;
    }

    logger.info('Molecular analysis app initialized');
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
      textInput.placeholder = `Describe an object... (${shortcutKey} to focus)`;
    }

    document.addEventListener('imageAnalysisComplete', (e) => {
      const { output, icon, objectName, useQuotes, croppedImageData } = e.detail;
      this.processAnalysisResult(output, icon, objectName, useQuotes, croppedImageData);
    });

    // Camera mode selection handler
    const cameraMode = document.getElementById('camera-mode');
    if (cameraMode) {
      cameraMode.addEventListener('change', async (e) => {
        if (e.target.checked) {
          logger.cameraEvent('camera_mode_activated');
          try {
            await cameraManager.initialize();
            const permissionGranted = await cameraManager.requestPermission();
            if (!permissionGranted) {
              e.target.checked = false;
            }
          } catch (error) {
            logger.error('Camera initialization failed', error);
            e.target.checked = false;
          }
        }
      });
    }

    // Card management button handler
    const cardBtn = document.getElementById('card-icon-btn');
    if (cardBtn) {
      cardBtn.addEventListener('click', () => {
        logger.userAction('payment_sidebar_toggle');
        const paymentSection = document.getElementById('payment-section');
        if (paymentSection) {
          paymentSection.classList.toggle('collapsed');
        }
      });
    }
  }

  setupTextAnalysis() {
    const textInput = document.getElementById('object-input');
    if (textInput) {
      textInput.addEventListener('keydown', (event) => {
        if (event.key === 'Enter') {
          event.preventDefault();
          this.handleTextAnalysis();
        }
      });
    }
  }

  async handleTextAnalysis() {
    if (this.isProcessing) return;

    const inputValue = this.objectInput.value.trim();
    if (!inputValue) return;

    // Check payment setup for sidebar
    const paymentSetup = await this.checkPaymentSetupForSidebar();
    if (!paymentSetup) {
      logger.warn('Payment not set up, showing message');
      this.showError('Payment setup required. Complete setup in the sidebar on the right.');
      return;
    }

    this.isProcessing = true;
    this.currentAnalysisType = 'text';
    
    try {
      this.showProcessing();
      
      logger.analysisEvent('text_analysis_started', { input: inputValue });
      
      // Create AbortController for timeout
      const controller = new AbortController();
      const timeoutId = setTimeout(() => controller.abort(), 45000); // 45 second timeout
      
      const response = await fetch("/analyze-text", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ text: inputValue }),
        signal: controller.signal
      });

      clearTimeout(timeoutId);

      if (!response.ok) {
        const errorText = await response.text();
        throw new Error(`Failed to analyze text: ${response.status} ${errorText}`);
      }

      const result = await response.json();
      logger.analysisEvent('text_analysis_completed', { input: inputValue, result });

      this.lastAnalysis = {
        type: 'text',
        input: inputValue,
        result: result
      };

      this.processAnalysisResult(result, null, inputValue, false, null);
      
    } catch (error) {
      logger.error('Text analysis failed', error);
      
      // Enhanced error handling with user-friendly messages
      let errorMessage = 'Analysis failed: ';
      
      if (error.name === 'AbortError') {
        errorMessage = 'Analysis timed out. Please check your internet connection and try again.';
      } else if (error.message.includes('Failed to fetch') || error.message.includes('fetch')) {
        errorMessage = 'Network connection failed. Please check your internet connection and try again.';
      } else if (error.message.includes('Network connection failed')) {
        errorMessage = 'Unable to connect to analysis service. Please check your internet connection.';
      } else if (error.message.includes('Rate limit exceeded')) {
        errorMessage = 'Too many requests. Please wait a moment and try again.';
      } else if (error.message.includes('503') || error.message.includes('temporarily unavailable')) {
        errorMessage = 'Analysis service temporarily unavailable. Please try again in a few moments.';
      } else {
        errorMessage += error.message;
      }
      
      this.showError(errorMessage);
    } finally {
      this.hideProcessing();
      this.isProcessing = false;
    }
  }

  // Check payment setup for sidebar-based system
  async checkPaymentSetupForSidebar() {
    // Check if dev mode is enabled (localhost auto-enable)
    if (this.hasPaymentSetup === true) {
      logger.info('Developer mode active - bypassing payment check');
      return true;
    }
    
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      // Ensure payment section is visible
      const paymentSection = document.getElementById('payment-section');
      if (paymentSection) {
        paymentSection.classList.remove('hidden');
      }
      return false;
    }
    
    try {
      const response = await fetch('/validate-payment', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ device_token: deviceToken })
      });
      
      if (!response.ok) {
        localStorage.removeItem('molDeviceToken');
        localStorage.removeItem('molCardInfo');
        // Ensure payment section is visible
        const paymentSection = document.getElementById('payment-section');
        if (paymentSection) {
          paymentSection.classList.remove('hidden');
        }
        return false;
      }
      
      return true;
      
    } catch (error) {
      logger.error('Payment validation error', error);
      return true; // Fallback to allow analysis
    }
  }

  showProcessing() {
    // Create loading indicator
    const loadingDiv = document.createElement('div');
    loadingDiv.className = 'loading-indicator';
    loadingDiv.innerHTML = `
      <div class="spinner"></div>
      <div>Analyzing...</div>
    `;
    loadingDiv.style.cssText = `
      position: fixed;
      top: 50%;
      left: 50%;
      transform: translate(-50%, -50%);
      background: rgba(0, 0, 0, 0.9);
      color: white;
      padding: 20px;
      border-radius: 10px;
      text-align: center;
      z-index: 1000;
    `;
    
    document.body.appendChild(loadingDiv);
    this.loadingIndicator = loadingDiv;
  }

  hideProcessing() {
    if (this.loadingIndicator) {
      this.loadingIndicator.remove();
      this.loadingIndicator = null;
    }
  }

  showError(message) {
    logger.error('Application error', { message });
    // Create simple error display
    const errorDiv = document.createElement('div');
    errorDiv.className = 'error-message';
    errorDiv.textContent = message;
    errorDiv.style.cssText = `
      background: rgba(255, 68, 68, 0.1);
      border: 1px solid #ff4444;
      border-radius: 4px;
      padding: 12px;
      margin: 12px 16px;
      color: #ff4444;
      font-size: 12px;
      text-align: center;
    `;
    
    const resultsSection = document.querySelector('.results-section');
    if (resultsSection) {
      // Remove existing error messages
      const existingErrors = resultsSection.querySelectorAll('.error-message');
      existingErrors.forEach(err => err.remove());
      
      resultsSection.insertBefore(errorDiv, resultsSection.firstChild);
      
      // Auto-remove after 5 seconds
      setTimeout(() => errorDiv.remove(), 5000);
    }
  }

  async processAnalysisResult(output, icon, objectName, useQuotes, croppedImageData) {
    logger.analysisEvent('result_processed', { 
      type: this.currentAnalysisType, 
      objectName,
      hasCroppedData: !!croppedImageData 
    });

    try {
      if (!output || typeof output !== 'object') {
        throw new Error('Invalid analysis output');
      }

      let displayName = objectName;
      if (useQuotes && objectName) {
        displayName = `"${objectName}"`;
      }

      // Handle different output formats
      let chemicals = [];
      let description = null;
      let errorMessage = null;
      
      if (output.chemicals && Array.isArray(output.chemicals)) {
        chemicals = output.chemicals;
      } else if (output.description) {
        description = output.description;
      } else if (output.error) {
        errorMessage = output.error;
        description = `Error: ${output.error}`;
      } else if (output.summary) {
        description = output.summary;
      }

      // Generate SDFs for chemicals if we have them
      let sdfFiles = [];
      let smiles = [];
      
      if (chemicals.length > 0) {
        const sdfResult = await this.generateSDFs(chemicals);
        sdfFiles = sdfResult.sdfFiles;
        smiles = sdfResult.smiles;
      }

      // Create the display column
      await this.createObjectColumn(
        displayName || "Analysis Result", 
        sdfFiles, 
        smiles, 
        errorMessage, 
        output.summary || null, 
        output.skippedChemicals || [], 
        description,
        chemicals,
        croppedImageData
      );

    } catch (error) {
      logger.error('Failed to process analysis result', error);
      this.showError(`Failed to display results: ${error.message}`);
    }
  }

  async generateSDFs(chemicals) {
    const sdfFiles = [];
    const smiles = [];

    // Extract SMILES from all chemicals
    const smilesArray = chemicals.filter(c => c.smiles).map(c => c.smiles);
    
    if (smilesArray.length === 0) {
      return { sdfFiles, smiles };
    }

    try {
      // Call the correct endpoint with array of SMILES
      const response = await fetch("/generate-sdfs", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smiles: smilesArray }),
      });

      if (response.ok) {
        const result = await response.json();
        logger.info('SDF generation result:', result);
        
        // Process the returned SDF paths
        for (let i = 0; i < result.sdfPaths.length; i++) {
          const sdfPath = result.sdfPaths[i];
          smiles.push(smilesArray[i]);
          sdfFiles.push(sdfPath); // Use the path directly since it's already a URL
          logger.info(`Generated SDF: ${smilesArray[i]} -> ${sdfPath}`);
        }
        
        if (result.errors && result.errors.length > 0) {
          logger.warn('Some SDF generation failed:', result.errors);
        }
        
      } else {
        const errorText = await response.text();
        logger.error(`Failed to generate SDFs: ${response.status} ${errorText}`);
      }
    } catch (error) {
      logger.error(`Error generating SDFs:`, error);
    }

    return { sdfFiles, smiles };
  }

  async createObjectColumn(objectName, sdfFiles, smiles = [], errorMessage = null, 
                          summary = null, skippedChemicals = [], description = null, 
                          chemicals = null, croppedImageData = null) {

    const gldiv = document.getElementById("gldiv");
    if (!gldiv) {
      logger.error("No gldiv found for molecular display");
      return;
    }

    // Create object column
    const objectColumn = document.createElement("div");
    objectColumn.className = "object-column";

    // Create title with close button
    const titleContainer = document.createElement("div");
    titleContainer.className = "object-title";
    
    const titleText = document.createElement("span");
    titleText.textContent = objectName;
    titleContainer.appendChild(titleText);

    const closeButton = document.createElement("button");
    closeButton.innerHTML = "✕";
    closeButton.className = "close-button";
    closeButton.title = "Close";
    closeButton.onclick = () => {
      objectColumn.remove();
      this.updateScrollHandles();
    };
    titleContainer.appendChild(closeButton);
    
    objectColumn.appendChild(titleContainer);

    // Handle description vs molecules
    if (description) {
      const descDiv = document.createElement("div");
      descDiv.className = "description-content";
      descDiv.textContent = description;
      descDiv.style.cssText = `
        color: rgba(255, 255, 255, 0.7);
        font-size: 13px;
        line-height: 1.4;
        padding: 12px;
        background: rgba(255, 255, 255, 0.02);
        border-radius: 4px;
        max-height: 200px;
        overflow-y: auto;
      `;
      objectColumn.appendChild(descDiv);
    } else {
      // Render 3D molecules with SPHERE REPRESENTATION ONLY
      logger.info(`Rendering ${sdfFiles.length} molecules`);
      
      for (let i = 0; i < sdfFiles.length; i++) {
        const sdfFile = sdfFiles[i];
        const chemical = chemicals?.[i] || { smiles: smiles[i] };
        
        logger.info(`Rendering molecule ${i + 1}/${sdfFiles.length}: ${sdfFile}`);
        
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
        if (viewer) {
          this.viewers.push(viewer);
          logger.info(`Successfully rendered molecule: ${displayName}`);
        } else {
          logger.error(`Failed to render molecule: ${displayName}`);
        }
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
      logger.info(`Fetching SDF file: ${sdfFile}`);
      const response = await fetch(sdfFile);
      if (!response.ok) {
        throw new Error(`HTTP error ${response.status}: ${response.statusText}`);
      }

      const sdfData = await response.text();
      logger.info(`SDF data length: ${sdfData.length} characters`);
      
      if (!window.$3Dmol) {
        throw new Error('3DMol.js library not loaded');
      }
      
      const viewer = $3Dmol.createViewer(container, {
        defaultcolors: $3Dmol.rasmolElementColors,
      });

      viewer.addModel(sdfData, "sdf");
      // CRITICAL: Use ONLY sphere representation with van der Waals radii at 0.8 scale
      viewer.setStyle({}, { sphere: { scale: 0.8 } });
      viewer.zoomTo();
      viewer.render();
      
      logger.info('3D molecule viewer created successfully');
      return viewer;
      
    } catch (error) {
      logger.error(`Failed to load molecule:`, error);
      container.textContent = `Error loading molecule: ${error.message}`;
      container.className += " error-text";
      return null;
    }
  }

  updateScrollHandles() {
    // Simple implementation for horizontal scroll management
    const gldiv = document.getElementById("gldiv");
    if (gldiv) {
      const hasOverflow = gldiv.scrollWidth > gldiv.clientWidth;
      gldiv.style.overflowX = hasOverflow ? 'auto' : 'hidden';
    }
  }

  clearResults() {
    const gldiv = document.getElementById("gldiv");
    if (gldiv) {
      gldiv.innerHTML = "";
    }
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
