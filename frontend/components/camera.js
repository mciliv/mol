// camera.js - Camera management and interaction module

import { paymentManager } from './payment.js';
import { uiManager } from './ui-utils.js';

class CameraManager {
  constructor() {
    this.video = null;
    this.currentStream = null;
    this.facingMode = 'environment';
    this.isMobile = /Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(navigator.userAgent);
    this.isSafari = /^((?!chrome|android).)*safari/i.test(navigator.userAgent);
    this.isIOS = /iPad|iPhone|iPod/.test(navigator.userAgent);
  }

  // Initialize camera system
  async initialize() {
    this.video = document.getElementById("video-feed");
    
    if (!navigator.mediaDevices || !navigator.mediaDevices.getUserMedia) {
      this.showError("Camera API not supported in this browser");
      return false;
    }

    const isSecureContext = window.isSecureContext || location.protocol === "https:";
    if (!isSecureContext && !location.hostname.includes("localhost")) {
      this.showError("Camera requires HTTPS on mobile devices. Please use HTTPS or localhost.");
      return false;
    }

    this.setupEventListeners();
    
    try {
      await this.startCamera();
      await this.setupSwitchCamera();
      this.updateReticleVisibility();
      return true;
    } catch (err) {
      console.error("Camera setup failed:", err);
      return false;
    }
  }

  // Setup camera event listeners
  setupEventListeners() {
    this.video.addEventListener("click", (e) => this.handleInteraction(e));
    this.video.addEventListener("touchstart", (e) => {
      e.preventDefault();
      this.handleInteraction(e.touches[0]);
    });

    // Safari-specific touch handling
    if (this.isSafari || this.isIOS) {
      this.video.addEventListener("touchend", (e) => e.preventDefault());
      
      let lastTouchEnd = 0;
      this.video.addEventListener("touchend", (e) => {
        const now = new Date().getTime();
        if (now - lastTouchEnd <= 300) {
          e.preventDefault();
        }
        lastTouchEnd = now;
      }, false);
    }
  }

  // Start camera stream
  async startCamera() {
    if (paymentManager.isPaymentRequired()) {
      console.log('Camera access blocked - payment setup required');
      return;
    }

    this.stopCurrentStream();

    try {
      let stream;

      if (this.isSafari || this.isIOS) {
        try {
          stream = await navigator.mediaDevices.getUserMedia(this.getSafariConstraints());
        } catch (err) {
          console.log("Safari constraints failed, trying basic:", err);
          stream = await navigator.mediaDevices.getUserMedia(this.getBasicConstraints());
        }
      } else {
        try {
          stream = await navigator.mediaDevices.getUserMedia(this.getSimpleConstraints());
        } catch (err) {
          stream = await navigator.mediaDevices.getUserMedia(this.getBasicConstraints());
        }
      }

      this.currentStream = stream;
      this.video.srcObject = stream;

      // Safari-specific video attributes
      this.video.setAttribute("playsinline", "true");
      this.video.setAttribute("webkit-playsinline", "true");
      this.video.setAttribute("x-webkit-airplay", "allow");

      await this.video.play();
      this.hideError();
      this.updateReticleVisibility();
      
    } catch (err) {
      console.error("Camera error:", err);
      this.handleCameraError(err);
    }
  }

  // Stop current camera stream
  stopCurrentStream() {
    if (this.currentStream) {
      this.currentStream.getTracks().forEach(track => track.stop());
      this.currentStream = null;
    }
  }

  // Get camera constraints for different browsers
  getSimpleConstraints() {
    return {
      video: {
        facingMode: this.facingMode,
        width: { ideal: 1280 },
        height: { ideal: 720 },
      },
    };
  }

  getBasicConstraints() {
    return {
      video: {
        width: { ideal: 1280 },
        height: { ideal: 720 },
      },
    };
  }

  getSafariConstraints() {
    return {
      video: {
        facingMode: this.facingMode,
        width: { min: 640, ideal: 1280, max: 1920 },
        height: { min: 480, ideal: 720, max: 1080 },
      },
    };
  }

  // Handle camera interaction (click/tap)
  async handleInteraction(evt) {
    if (!await paymentManager.checkPaymentMethod()) {
      // Show simple message instead of popdown
      const messageColumn = uiManager.createColumn("See payment setup above", "payment-required");
      messageColumn.innerHTML = `
        <div class="molecule-container">
          <div class="molecule-info">
            <h3>Payment Required</h3>
            <p>See payment setup above</p>
            <div class="analysis-note">Complete payment setup to analyze molecules via camera</div>
          </div>
        </div>
      `;
      return;
    }

    // Mobile reticle validation
    if (this.isMobile && !this.isWithinReticle(evt)) {
      this.showReticleFeedback();
      return;
    }

    this.showCropOutline(evt);
    
    try {
      // Capture and analyze the image
      const captureData = await this.captureAndAnalyze(evt);
      
      // Create loading column
      const loadingColumn = uiManager.createLoadingColumn("Analyzing...", captureData.croppedBase64);
      
      // Send to server for analysis
      const response = await fetch("/image-molecules", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          imageBase64: captureData.imageBase64,
          croppedImageBase64: captureData.croppedBase64,
          x: captureData.coordinates.x,
          y: captureData.coordinates.y,
          cropMiddleX: captureData.coordinates.cropMiddleX,
          cropMiddleY: captureData.coordinates.cropMiddleY,
          cropSize: captureData.coordinates.cropSize,
        }),
      });

      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      const { output } = await response.json();

      // Remove loading column
      loadingColumn.remove();
      
      // Emit analysis result event for app to handle
      const event = new CustomEvent('imageAnalysisComplete', {
        detail: { 
          output, 
          icon: "📷", 
          objectName: output.object || "Camera capture", 
          useQuotes: false, 
          croppedImageData: captureData.croppedBase64 
        }
      });
      document.dispatchEvent(event);
      
      // Increment usage
      await paymentManager.incrementUsage();
      
    } catch (error) {
      console.error('Camera analysis error:', error);
      this.createClosableErrorMessage(`Analysis failed: ${error.message}`);
    }
  }

  // Check if tap is within mobile reticle
  isWithinReticle(evt) {
    const rect = this.video.getBoundingClientRect();
    const centerX = rect.left + rect.width / 2;
    const centerY = rect.top + rect.height / 2;
    const tapX = evt.clientX;
    const tapY = evt.clientY;

    const reticleRadius = 50;
    const distanceFromCenter = Math.sqrt(
      Math.pow(tapX - centerX, 2) + Math.pow(tapY - centerY, 2)
    );

    return distanceFromCenter <= reticleRadius;
  }

  // Show reticle feedback animation
  showReticleFeedback() {
    const reticle = document.querySelector(".mobile-reticle");
    if (reticle) {
      reticle.style.animation = "reticlePulse 0.5s ease";
      setTimeout(() => {
        reticle.style.animation = "";
      }, 500);
    }
  }

  // Show crop outline at interaction point
  showCropOutline(evt) {
    const cropSize = 100;
    const outline = document.createElement("div");
    outline.className = "crop-outline";
    outline.style.width = cropSize + "px";
    outline.style.height = cropSize + "px";
    outline.style.left = evt.clientX - cropSize / 2 + "px";
    outline.style.top = evt.clientY - cropSize / 2 + "px";
    document.body.appendChild(outline);
    
    setTimeout(() => {
      outline.style.opacity = "0";
      setTimeout(() => outline.remove(), 200);
    }, 500);
  }

  // Capture image and prepare for analysis
  async captureAndAnalyze(evt) {
    const canvas = document.createElement("canvas");
    canvas.width = this.video.videoWidth;
    canvas.height = this.video.videoHeight;
    canvas.getContext("2d").drawImage(this.video, 0, 0, canvas.width, canvas.height);
    const imageBase64 = canvas.toDataURL("image/jpeg", 0.9).split(",")[1];

    const scaleX = this.video.videoWidth / this.video.clientWidth;
    const scaleY = this.video.videoHeight / this.video.clientHeight;
    const clickX = Math.round(evt.clientX - this.video.getBoundingClientRect().left);
    const clickY = Math.round(evt.clientY - this.video.getBoundingClientRect().top);
    const actualX = Math.round(clickX * scaleX);
    const actualY = Math.round(clickY * scaleY);

    const crop = document.createElement("canvas");
    crop.width = crop.height = 100;
    const cropCtx = crop.getContext("2d");
    cropCtx.imageSmoothingEnabled = false;
    cropCtx.drawImage(canvas, actualX - 50, actualY - 50, 100, 100, 0, 0, 100, 100);

    const middleX = Math.floor(100 / 2);
    const middleY = Math.floor(100 / 2);
    const boxSize = Math.max(8, Math.floor(100 * 0.1));
    
    cropCtx.save();
    cropCtx.strokeStyle = "#ff0000";
    cropCtx.lineWidth = Math.max(2, Math.floor(100 * 0.02));
    cropCtx.strokeRect(middleX - boxSize / 2, middleY - boxSize / 2, boxSize, boxSize);
    cropCtx.restore();

    const croppedBase64 = crop.toDataURL("image/jpeg", 0.9).split(",")[1];

    return {
      imageBase64,
      croppedBase64,
      coordinates: { x: actualX, y: actualY, cropMiddleX: middleX, cropMiddleY: middleY, cropSize: 100 }
    };
  }

  // Setup camera switching functionality
  async setupSwitchCamera() {
    const devices = await navigator.mediaDevices.enumerateDevices();
    const videoDevices = devices.filter(d => d.kind === "videoinput");
    
    if (videoDevices.length > 1) {
      const switchCameraContainer = document.querySelector(".switch-camera-container");
      const switchCameraBtn = document.getElementById("switch-camera-btn");
      
      if (switchCameraContainer && switchCameraBtn) {
        switchCameraContainer.style.display = "block";
        switchCameraBtn.onclick = () => {
          this.facingMode = this.facingMode === "user" ? "environment" : "user";
          this.startCamera();
        };
      }
    }
  }

  // Update mobile reticle visibility
  updateReticleVisibility() {
    if (this.isMobile) {
      this.addMobileTargetingReticle();
    } else {
      this.removeMobileTargetingReticle();
    }
  }

  // Add mobile targeting reticle
  addMobileTargetingReticle() {
    this.removeMobileTargetingReticle();

    const reticle = document.createElement("div");
    reticle.className = "mobile-reticle";
    const template = document.getElementById("mobile-reticle-template");
    const clone = template.content.cloneNode(true);
    reticle.appendChild(clone);

    this.video.appendChild(reticle);

    setTimeout(() => {
      reticle.style.animation = "reticlePulse 2s ease-in-out infinite";
    }, 1000);
  }

  // Remove mobile targeting reticle
  removeMobileTargetingReticle() {
    const existingReticle = document.querySelector(".mobile-reticle");
    if (existingReticle) {
      existingReticle.remove();
    }
  }

  // Handle camera errors
  handleCameraError(err) {
    let errorMessage;
    
    if (err.name === "NotAllowedError") {
      errorMessage = "Camera access denied. Please allow camera access in browser settings.";
    } else if (err.name === "NotFoundError") {
      errorMessage = "No camera found on this device.";
    } else if (this.isSafari && err.name === "NotSupportedError") {
      errorMessage = "Camera not supported in Safari. Try using Chrome or Firefox.";
    } else {
      errorMessage = `Camera error: ${err.message}`;
    }
    
    this.showError(errorMessage);
  }

  // Show error message
  showError(message) {
    const msgBox = document.querySelector(".permission-message");
    if (msgBox) {
      msgBox.hidden = false;
      msgBox.textContent = message;
    }
  }

  // Hide error message
  hideError() {
    const msgBox = document.querySelector(".permission-message");
    if (msgBox) {
      msgBox.hidden = true;
    }
  }

  // Create closable error message
  createClosableErrorMessage(message) {
    const snapshots = document.querySelector(".snapshots-container");
    const errorDiv = uiManager.createErrorMessage(message, snapshots);
    return errorDiv;
  }

  // Request camera permission (Safari)
  async requestPermission() {
    try {
      console.log('🎥 Requesting camera permission...');
      
      const stream = await navigator.mediaDevices.getUserMedia({ 
        video: { facingMode: 'environment' } 
      });
      
      stream.getTracks().forEach(track => track.stop());
      
      console.log('✅ Camera permission granted');
      return true;
    } catch (error) {
      console.error('❌ Camera permission denied:', error);
      
      if (error.name === 'NotAllowedError') {
        alert('Please allow camera access in Safari:\n1. Click Safari menu → Settings → Websites → Camera\n2. Set this site to "Allow"');
      }
      
      return false;
    }
  }

  // Cleanup resources
  cleanup() {
    this.stopCurrentStream();
    this.removeMobileTargetingReticle();
  }
}

// Export singleton instance
export const cameraManager = new CameraManager(); 