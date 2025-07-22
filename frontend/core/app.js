// Simple App Logic - Following ui.mdc Guidelines
import { simplePaymentManager } from '../components/simple-payment.js';

// Initialize app
document.addEventListener('DOMContentLoaded', () => {
  console.log('🚀 App starting...');
  
  // Check payment status
  simplePaymentManager.checkPaymentRequired();
  
  // Setup event listeners
  setupEventListeners();
  
  console.log('✅ App ready');
});

function setupEventListeners() {
  // Account button
  const accountBtn = document.getElementById('account-btn');
  if (accountBtn) {
    accountBtn.addEventListener('click', handleAccountClick);
  }
  
  // Main input
  const mainInput = document.getElementById('object-input');
  if (mainInput) {
    mainInput.addEventListener('keypress', handleInputKeypress);
  }
  
  // Mode toggles
  const cameraMode = document.getElementById('camera-mode');
  const photoMode = document.getElementById('photo-mode');
  
  if (cameraMode) {
    cameraMode.addEventListener('change', handleCameraMode);
  }
  
  if (photoMode) {
    photoMode.addEventListener('change', handlePhotoMode);
  }
  
  // URL analyze button
  const urlAnalyze = document.getElementById('url-analyze');
  if (urlAnalyze) {
    urlAnalyze.addEventListener('click', handleUrlAnalyze);
  }
  
  // Photo upload
  const photoUpload = document.getElementById('photo-upload');
  if (photoUpload) {
    photoUpload.addEventListener('change', handlePhotoUpload);
  }
}

// Account button handler
function handleAccountClick() {
  const cardInfo = localStorage.getItem('molCardInfo');
  
  if (!cardInfo) {
    simplePaymentManager.showPaymentSection();
  } else {
    // Show account info or manage payment
    console.log('💳 Account management not implemented yet');
  }
}

// Input keypress handler
function handleInputKeypress(event) {
  if (event.key === 'Enter') {
    const query = event.target.value.trim();
    if (query) {
      console.log('🔍 Analyzing:', query);
      analyzeInput(query);
    }
  }
}

// Camera mode handler
function handleCameraMode(event) {
  if (event.target.checked) {
    console.log('📸 Camera mode activated');
    // Uncheck photo mode
    const photoMode = document.getElementById('photo-mode');
    if (photoMode) photoMode.checked = false;
    
    // Initialize camera (placeholder)
    console.log('📷 Camera initialization not implemented yet');
  }
}

// Photo mode handler
function handlePhotoMode(event) {
  if (event.target.checked) {
    console.log('🖼️ Photo mode activated');
    // Uncheck camera mode
    const cameraMode = document.getElementById('camera-mode');
    if (cameraMode) cameraMode.checked = false;
  }
}

// URL analyze handler
function handleUrlAnalyze() {
  const urlInput = document.getElementById('photo-url');
  if (urlInput && urlInput.value.trim()) {
    console.log('🔗 Analyzing URL:', urlInput.value);
    analyzeUrl(urlInput.value);
  }
}

// Photo upload handler
function handlePhotoUpload(event) {
  const file = event.target.files[0];
  if (file) {
    console.log('📁 File uploaded:', file.name);
    analyzeFile(file);
  }
}

// Analysis functions (placeholders)
async function analyzeInput(query) {
  // Check payment
  if (!await simplePaymentManager.checkPaymentMethod()) {
    return;
  }
  
  console.log('🧪 Analyzing input:', query);
  // Placeholder for actual analysis
  
  // Increment usage
  simplePaymentManager.incrementUsage();
}

async function analyzeUrl(url) {
  // Check payment
  if (!await simplePaymentManager.checkPaymentMethod()) {
    return;
  }
  
  console.log('🧪 Analyzing URL:', url);
  // Placeholder for actual analysis
  
  // Increment usage
  simplePaymentManager.incrementUsage();
}

async function analyzeFile(file) {
  // Check payment
  if (!await simplePaymentManager.checkPaymentMethod()) {
    return;
  }
  
  console.log('🧪 Analyzing file:', file.name);
  // Placeholder for actual analysis
  
  // Increment usage
  simplePaymentManager.incrementUsage();
}

// Developer shortcut (for testing)
window.devSetup = () => {
  simplePaymentManager.setupDeveloperAccount();
};
