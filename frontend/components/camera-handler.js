// camera-handler.js - Consolidated camera and image handling functionality

import { paymentManager } from './payment.js';
import { uiManager } from './ui-utils.js';

class CameraHandler {
  constructor() {
    this.isMobile = /Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(navigator.userAgent);
  }

  // Handle photo upload analysis
  async handlePhotoUpload(e) {
    const file = e.target.files[0];
    if (!file) return;

    if (!file.type.startsWith("image/")) {
      alert("Please select an image file");
      return;
    }

    await this.displayUploadedImage(file);
  }

  // Handle URL-based image analysis
  async handleUrlAnalysis() {
    const photoUrl = document.getElementById("photo-url");
    const url = photoUrl.value.trim();
    
    if (!url) {
      alert("Please enter an image URL");
      return;
    }

    try {
      new URL(url);
    } catch {
      alert("Please enter a valid URL");
      return;
    }

    try {
      const imageBase64 = await uiManager.urlToBase64(url);
      const byteCharacters = atob(imageBase64);
      const byteNumbers = new Array(byteCharacters.length);
      
      for (let i = 0; i < byteCharacters.length; i++) {
        byteNumbers[i] = byteCharacters.charCodeAt(i);
      }
      
      const byteArray = new Uint8Array(byteNumbers);
      const blob = new Blob([byteArray], { type: "image/jpeg" });
      const file = new File([blob], "url-image.jpg", { type: "image/jpeg" });

      await this.displayUploadedImage(file);
      photoUrl.value = "";
      
    } catch (err) {
      this.createClosableErrorMessage(`Error loading image from URL: ${err.message}`);
    }
  }

  // Display uploaded image for interactive analysis
  async displayUploadedImage(file) {
    const photoOptions = document.getElementById("photo-options");
    photoOptions.innerHTML = "";

    const imageContainer = document.createElement("div");
    imageContainer.className = "uploaded-image-container";

    const img = document.createElement("img");
    
    // Add mobile reticle if needed
    if (this.isMobile) {
      const crosshair = document.createElement("div");
      crosshair.className = "crosshair";
      const beforeLine = document.createElement("div");
      beforeLine.className = "crosshair-line vertical";
      const afterLine = document.createElement("div");
      afterLine.className = "crosshair-line horizontal";
      crosshair.appendChild(beforeLine);
      crosshair.appendChild(afterLine);
      imageContainer.appendChild(crosshair);
    }

    const instructionText = document.createElement("div");
    instructionText.className = "instruction-text";
    instructionText.textContent = this.isMobile
      ? "Center object in circle & tap, or type name above"
      : "Click on object or type name above";

    const closeButton = document.createElement("button");
    closeButton.className = "close-button";
    closeButton.innerHTML = '<img src="../assets/close.svg" alt="Close" width="24" height="24" />';
    closeButton.onclick = () => {
      photoOptions.innerHTML = "";
      const template = document.getElementById("photo-upload-template");
      const clone = template.content.cloneNode(true);
      photoOptions.appendChild(clone);

      const newPhotoUpload = photoOptions.querySelector("#photo-upload");
      newPhotoUpload.addEventListener("change", (e) => this.handlePhotoUpload(e));
    };

    try {
      const imageBase64 = await uiManager.fileToBase64(file);
      img.src = `data:${file.type};base64,${imageBase64}`;
      img.dataset.base64 = imageBase64;
      img.addEventListener("click", (e) => this.handleImageClick(e, img));
      
    } catch (error) {
      this.createClosableErrorMessage(`Error processing image: ${error.message}`);
      return;
    }

    imageContainer.appendChild(img);
    imageContainer.appendChild(instructionText);
    imageContainer.appendChild(closeButton);
    photoOptions.appendChild(imageContainer);
  }

  // Handle click on uploaded image
  async handleImageClick(evt, img) {
    console.log('🖱️ Image click detected:', evt);
    
    // Temporarily bypass payment check for testing
    let paymentCheck = false;
    try {
      paymentCheck = await paymentManager.checkPaymentMethod();
      console.log('💳 Payment check result:', paymentCheck);
    } catch (error) {
      console.log('⚠️ Payment check failed, proceeding anyway:', error);
      paymentCheck = true; // Fallback to allow analysis
    }
    
    if (!paymentCheck) {
      console.log('🚫 Payment required - showing message');
      // Show simple message instead of popdown
      const messageColumn = uiManager.createColumn("See payment setup above", "payment-required");
      messageColumn.innerHTML = `
        <div class="molecule-container">
          <div class="molecule-info">
            <h3>Payment Required</h3>
            <p>See payment setup above</p>
            <div class="analysis-note">Complete payment setup to analyze molecules from images</div>
          </div>
        </div>
      `;
      return;
    }

    console.log('✅ Payment check passed, proceeding with analysis');
    
    const rect = img.getBoundingClientRect();
    const clickX = evt.clientX - rect.left;
    const clickY = evt.clientY - rect.top;

    const relativeX = clickX / rect.width;
    const relativeY = clickY / rect.height;

    console.log('📊 Click coordinates:', { 
      clickX, 
      clickY, 
      relativeX, 
      relativeY,
      actualX: relativeX * tempImg.width,
      actualY: relativeY * tempImg.height,
      imageWidth: tempImg.width,
      imageHeight: tempImg.height
    });

    const imageBase64 = img.dataset.base64;
    if (!imageBase64) {
      console.error('❌ No image data found');
      this.createClosableErrorMessage('No image data available for analysis');
      return;
    }

    // Create crop canvas
    const canvas = document.createElement("canvas");
    const ctx = canvas.getContext("2d");

    const tempImg = new Image();
    tempImg.onload = async () => {
      console.log('🖼️ Image loaded, processing crop');
      canvas.width = tempImg.width;
      canvas.height = tempImg.height;
      ctx.drawImage(tempImg, 0, 0);

      const cropSize = Math.min(tempImg.width, tempImg.height) * 0.1;
      const cropX = Math.max(0, relativeX * tempImg.width - cropSize / 2);
      const cropY = Math.max(0, relativeY * tempImg.height - cropSize / 2);

      const cropCanvas = document.createElement("canvas");
      cropCanvas.width = cropSize;
      cropCanvas.height = cropSize;
      const cropCtx = cropCanvas.getContext("2d");
      cropCtx.imageSmoothingEnabled = false;

      cropCtx.drawImage(canvas, cropX, cropY, cropSize, cropSize, 0, 0, cropSize, cropSize);

      const middleX = Math.floor(cropSize / 2);
      const middleY = Math.floor(cropSize / 2);
      const boxSize = Math.max(8, Math.floor(cropSize * 0.1));
      
      cropCtx.save();
      cropCtx.strokeStyle = "#ff0000";
      cropCtx.lineWidth = Math.max(2, Math.floor(cropSize * 0.02));
      cropCtx.strokeRect(middleX - boxSize / 2, middleY - boxSize / 2, boxSize, boxSize);
      cropCtx.restore();

      const croppedBase64 = cropCanvas.toDataURL("image/jpeg", 0.9).split(",")[1];
      console.log('✂️ Crop created, size:', cropSize);
      
      const loadingColumn = uiManager.createLoadingColumn("Analyzing...", croppedBase64);

      try {
        console.log('🌐 Sending analysis request to server');
        const response = await fetch("/image-molecules", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            imageBase64,
            croppedImageBase64: croppedBase64,
            x: relativeX * tempImg.width,
            y: relativeY * tempImg.height,
            cropMiddleX: middleX,
            cropMiddleY: middleY,
            cropSize: cropSize,
          }),
        });

        console.log('📡 Server response status:', response.status);
        if (!response.ok) {
          const errorText = await response.text();
          throw new Error(`HTTP ${response.status}: ${errorText}`);
        }
        
        const { output } = await response.json();
        console.log('✅ Analysis completed:', output);

        loadingColumn.remove();
        this.updateScrollHandles();

        const objectName = output.object || "Uploaded image";
        // Emit event for app to handle analysis result
        this.emitAnalysisResult(output, "Photo", objectName, false, croppedBase64);
        
        // Try to increment usage, but don't fail if it doesn't work
        try {
          await paymentManager.incrementUsage();
        } catch (usageError) {
          console.log('⚠️ Usage increment failed:', usageError);
        }
        
      } catch (err) {
        console.error('❌ Analysis error:', err);
        loadingColumn.remove();
        this.updateScrollHandles();
        this.createClosableErrorMessage(`Error: ${err.message}`);
      }
    };

    tempImg.onerror = () => {
      console.error('❌ Failed to load image for processing');
      this.createClosableErrorMessage('Failed to process image for analysis');
    };

    tempImg.src = `data:image/jpeg;base64,${imageBase64}`;
  }

  // Emit analysis result event for app to handle
  emitAnalysisResult(output, icon, objectName, useQuotes = false, croppedImageData = null) {
    const event = new CustomEvent('imageAnalysisComplete', {
      detail: { output, icon, objectName, useQuotes, croppedImageData }
    });
    document.dispatchEvent(event);
  }

  // Create error message
  createClosableErrorMessage(message) {
    const snapshots = document.querySelector(".snapshots-container");
    const errorDiv = uiManager.createErrorMessage(message, snapshots);
    this.updateScrollHandles();
    return errorDiv;
  }

  // Update scroll handles (placeholder)
  updateScrollHandles() {
    if (window.updateScrollHandles) {
      window.updateScrollHandles();
    }
  }

  // Setup event listeners for camera functionality
  setupEventListeners() {
    // Photo upload handling
    const photoUpload = document.getElementById("photo-upload");
    if (photoUpload) {
      photoUpload.addEventListener("change", (e) => this.handlePhotoUpload(e));
    }
  
    // URL analysis
    const photoUrl = document.getElementById("photo-url");
    const urlAnalyze = document.getElementById("url-analyze");
    if (photoUrl && urlAnalyze) {
      photoUrl.addEventListener("keyup", (e) => {
        if (e.key === "Enter") this.handleUrlAnalysis();
      });
      urlAnalyze.addEventListener("click", () => this.handleUrlAnalysis());
    }
  }
}

// Create and export singleton instance
export const cameraHandler = new CameraHandler(); 