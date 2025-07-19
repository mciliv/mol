document.addEventListener("DOMContentLoaded", () => {
  // DEVELOPMENT TOGGLE: Programmatic bypass for payment requirements
  // Usage: Run `window.bypassPayment()` in console to toggle
  window.bypassPayment = function() {
    const mainInterface = document.getElementById('main-app-interface');
    const paymentPopdown = document.getElementById('payment-popdown');
    
    if (mainInterface.classList.contains('payment-required')) {
      console.log('🔧 Bypassing payment requirements');
      mainInterface.classList.remove('payment-required');
      mainInterface.style.opacity = '1';
      mainInterface.style.filter = 'none';
      mainInterface.style.pointerEvents = 'auto';
      if (paymentPopdown) {
        paymentPopdown.style.display = 'none';
      }
    } else {
      console.log('🔒 Restoring payment requirements');
      mainInterface.classList.add('payment-required');
      mainInterface.style.opacity = '';
      mainInterface.style.filter = '';
      mainInterface.style.pointerEvents = '';
      if (paymentPopdown) {
        paymentPopdown.style.display = 'block';
      }
    }
  };

  // SHOW BOTH SECTIONS: Payment section + clear app interface
  window.showBothSections = function() {
    console.log('🔧 Showing both payment and app interface clearly');
    const mainInterface = document.getElementById('main-app-interface');
    const paymentPopdown = document.getElementById('payment-popdown');
    
    // Make sure both sections are visible
    if (paymentPopdown) {
      paymentPopdown.style.display = 'block';
      paymentPopdown.style.visibility = 'visible';
    }
    
    if (mainInterface) {
      mainInterface.style.display = 'block';
      mainInterface.style.visibility = 'visible';
      mainInterface.style.opacity = '1';
      mainInterface.style.filter = 'none';
      mainInterface.style.pointerEvents = 'auto';
      mainInterface.classList.remove('payment-required');
    }
    
    console.log('✅ Both sections should now be clearly visible');
  };

  // SIMPLE APPROACH: Just show the app interface
  window.showApp = function() {
    console.log('🔧 Simple app show');
    const mainInterface = document.getElementById('main-app-interface');
    console.log('Main interface element:', mainInterface);
    
    if (mainInterface) {
      mainInterface.style.display = 'block';
      mainInterface.style.opacity = '1';
      mainInterface.style.filter = 'none';
      mainInterface.style.pointerEvents = 'auto';
      mainInterface.style.visibility = 'visible';
      mainInterface.classList.remove('payment-required');
      console.log('✅ App interface forced visible');
    } else {
      console.error('❌ Main interface element not found');
    }
  };

  // CHECK WHAT EXISTS
  window.checkElements = function() {
    console.log('📋 Element Check:');
    console.log('main-app-interface:', document.getElementById('main-app-interface'));
    console.log('payment-popdown:', document.getElementById('payment-popdown'));
    console.log('app-container:', document.querySelector('.app-container'));
    console.log('top-bar:', document.querySelector('.top-bar'));
  };
  // END DEVELOPMENT TOGGLE

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

  // ==================== ACCOUNT SETUP INTEGRATION ====================
  
  let stripe;
  let cardElement;
  let paymentRequest;
  let setupInProgress = false;
  
  // Check for existing payment setup first
  async function checkInitialPaymentSetup() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    // Always show main app interface
    showMainApp();
    
    if (!deviceToken || !cardInfo) {
      // No payment setup - show popdown
      showPaymentPopdown();
      return false;
    }
    
    try {
      // Quick validation with server
      const response = await fetch('/validate-payment', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ device_token: deviceToken })
      });
      
      if (!response.ok) {
        // Invalid token - clear and show setup
        localStorage.removeItem('molDeviceToken');
        localStorage.removeItem('molCardInfo');
        showPaymentPopdown();
        return false;
      }
      
      const result = await response.json();
      
      // Update account status display
      updateAccountStatus(result.user);
      
      return true;
      
    } catch (error) {
      console.error('Initial payment check failed:', error);
      // Network error - allow offline usage without popdown
      return true;
    }
  }
  
  // Show payment popdown (blocking access to app)
  function showPaymentPopdown() {
    const popdown = document.getElementById('payment-popdown');
    const mainInterface = document.getElementById('main-app-interface');
    
    // Show popdown and fade main interface
    popdown.style.display = 'block';
    mainInterface.classList.add('payment-required');
    
    // Trigger animation
    setTimeout(() => {
      popdown.classList.add('show');
    }, 10);
    
    // Initialize payment setup
    initializePaymentSetup();
  }
  
  // Hide payment popdown
  function hidePaymentPopdown() {
    const popdown = document.getElementById('payment-popdown');
    const mainInterface = document.getElementById('main-app-interface');
    
    // Hide popdown and restore main interface
    popdown.classList.remove('show');
    mainInterface.classList.remove('payment-required');
    
    setTimeout(() => {
      popdown.style.display = 'none';
    }, 400);
  }
  
  // Show main app interface
  function showMainApp() {
    document.getElementById('main-app-interface').style.display = 'block';
    
    // Initialize camera and app (but don't start camera if payment required)
    initializeMainApp();
  }
  
  // Update account status in top bar
  function updateAccountStatus(user) {
    const accountStatus = document.getElementById('account-status');
    const accountName = document.getElementById('account-name');
    
    if (user && user.name) {
      accountName.textContent = user.name;
    } else {
      accountName.textContent = 'Account';
    }
    
    accountStatus.style.display = 'flex';
  }
  
  // Initialize payment setup functionality
  async function initializePaymentSetup() {
    if (setupInProgress) return;
    setupInProgress = true;
    
    try {
      // Get Stripe configuration
      const response = await fetch('/stripe-config');
      const config = await response.json();
      
      stripe = Stripe(config.publishableKey);
      
      // Initialize Payment Request (Apple Pay, Google Pay, etc.)
      initializePaymentRequest();
      
      // Create card element
      const elements = stripe.elements();
      cardElement = elements.create('card', {
        style: {
          base: {
            fontSize: '16px',
            color: '#ffffff',
            fontFamily: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif',
            '::placeholder': {
              color: 'rgba(255, 255, 255, 0.6)'
            }
          },
          invalid: {
            color: '#ff6b6b',
            iconColor: '#ff6b6b'
          }
        },
        hidePostalCode: true
      });
      
      cardElement.mount('#card-element');
      
      // Handle card errors
      cardElement.on('change', function(event) {
        const displayError = document.getElementById('card-errors');
        if (event.error) {
          displayError.textContent = event.error.message;
          displayError.style.display = 'block';
        } else {
          displayError.textContent = '';
          displayError.style.display = 'none';
        }
      });
      
      // Setup form handler
      document.getElementById('card-setup-form').addEventListener('submit', handleCardSetup);
      document.getElementById('start-analyzing-btn').addEventListener('click', () => {
        hidePaymentPopdown();
      });
      
      // Setup autofill detection
      setupAutofillDetection();
      
    } catch (error) {
      console.error('Payment setup initialization failed:', error);
      // If payment setup fails, still allow access to the app
      showMainApp();
    }
  }
  
  // Initialize Payment Request for express payments
  function initializePaymentRequest() {
    paymentRequest = stripe.paymentRequest({
      country: 'US',
      currency: 'usd',
      total: {
        label: 'Molecular Analysis Setup',
        amount: 25, // $0.25 in cents
      },
      requestPayerName: true,
      requestPayerEmail: false,
    });

    paymentRequest.canMakePayment().then(function(result) {
      if (result) {
        const paymentRequestButton = paymentRequest.mount('#payment-request-button');
        document.getElementById('payment-request-button').style.display = 'block';
        document.getElementById('or-divider').style.display = 'block';
        document.getElementById('express-payment').style.display = 'block';
      }
    });

    paymentRequest.on('paymentmethod', async function(ev) {
      await handleExpressPayment(ev.paymentMethod, ev);
    });
  }
  
  // Setup autofill detection
  function setupAutofillDetection() {
    const autofillFields = document.querySelectorAll('.autofill-detector input');
    
    autofillFields.forEach(field => {
      field.addEventListener('input', handleAutofillDetection);
      field.addEventListener('change', handleAutofillDetection);
      
      setTimeout(() => {
        if (field.value) {
          handleAutofillDetection();
        }
      }, 100);
    });
  }
  
  // Handle autofill detection
  function handleAutofillDetection() {
    const nameField = document.querySelector('input[name="cc-name"]');
    const numberField = document.querySelector('input[name="cc-number"]');
    
    if (nameField && nameField.value) {
      const userNameField = document.getElementById('user-name');
      if (!userNameField.value) {
        userNameField.value = nameField.value;
      }
    }
    
    if (numberField && numberField.value) {
      const note = document.querySelector('.autofill-note');
      if (note) {
        note.innerHTML = `
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
            <path d="M9 12l2 2 4-4"/>
            <circle cx="12" cy="12" r="10"/>
          </svg>
          <span style="color: #00ff88;">Autofill detected - click below to use</span>
        `;
      }
    }
  }
  
  // Handle card setup form submission
  async function handleCardSetup(event) {
    event.preventDefault();
    
    const setupBtn = document.getElementById('setup-btn');
    const btnText = setupBtn.querySelector('.btn-text');
    const btnLoading = setupBtn.querySelector('.btn-loading');
    
    // Show loading state
    btnText.style.display = 'none';
    btnLoading.style.display = 'flex';
    setupBtn.disabled = true;
    
    try {
      // Create payment method with Stripe
      const { error, paymentMethod } = await stripe.createPaymentMethod({
        type: 'card',
        card: cardElement,
        billing_details: {
          name: document.getElementById('user-name').value || 'Molecular Analysis User',
        },
      });
      
      if (error) {
        throw new Error(error.message);
      }
      
      // Setup payment method with server
      const response = await fetch('/setup-payment-method', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          payment_method: paymentMethod.id,
          device_info: getDeviceFingerprint(),
          name: document.getElementById('user-name').value || null
        }),
      });
      
      if (!response.ok) {
        throw new Error('Setup failed. Please try again.');
      }
      
      const result = await response.json();
      
      // Store user info locally
      const deviceToken = generateDeviceToken();
      const cardInfo = {
        last4: paymentMethod.card.last4,
        brand: paymentMethod.card.brand,
        usage: 0,
        name: document.getElementById('user-name').value || null,
        setupDate: new Date().toISOString()
      };
      
      localStorage.setItem('molDeviceToken', deviceToken);
      localStorage.setItem('molCardInfo', JSON.stringify(cardInfo));
      
      // Show success state
      showSetupSuccess();
      
    } catch (error) {
      console.error('Setup error:', error);
      showSetupError(error.message);
      
      // Reset button state
      btnText.style.display = 'inline';
      btnLoading.style.display = 'none';
      setupBtn.disabled = false;
    }
  }
  
  // Handle express payments
  async function handleExpressPayment(paymentMethod, event) {
    try {
      const response = await fetch('/setup-payment-method', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          payment_method: paymentMethod.id,
          device_info: getDeviceFingerprint(),
          name: paymentMethod.billing_details.name || null,
          express_type: 'payment_request'
        }),
      });
      
      if (!response.ok) {
        throw new Error('Setup failed. Please try again.');
      }
      
      // Store user info locally
      const deviceToken = generateDeviceToken();
      const cardInfo = {
        last4: paymentMethod.card ? paymentMethod.card.last4 : '••••',
        brand: paymentMethod.card ? paymentMethod.card.brand : 'express',
        usage: 0,
        name: paymentMethod.billing_details.name || null,
        setupDate: new Date().toISOString(),
        type: 'express'
      };
      
      localStorage.setItem('molDeviceToken', deviceToken);
      localStorage.setItem('molCardInfo', JSON.stringify(cardInfo));
      
      event.complete('success');
      showSetupSuccess();
      
    } catch (error) {
      console.error('Express payment error:', error);
      event.complete('fail');
      showSetupError(error.message);
    }
  }
  
  // Show setup success
  function showSetupSuccess() {
    document.getElementById('card-setup-form').style.display = 'none';
    document.getElementById('express-payment').style.display = 'none';
    document.getElementById('setup-success').style.display = 'block';
  }
  
  // Show setup error
  function showSetupError(message) {
    const errorDiv = document.getElementById('card-errors');
    errorDiv.textContent = message;
    errorDiv.style.display = 'block';
    
    setTimeout(() => {
      errorDiv.style.display = 'none';
    }, 5000);
  }
  
  // Generate device fingerprint
  function getDeviceFingerprint() {
    const screen = `${window.screen.width}x${window.screen.height}`;
    const timezone = Intl.DateTimeFormat().resolvedOptions().timeZone;
    const language = navigator.language;
    const platform = navigator.platform;
    
    return btoa(`${screen}-${timezone}-${language}-${platform}`);
  }
  
  // Generate device token
  function generateDeviceToken() {
    const fingerprint = getDeviceFingerprint();
    const timestamp = Date.now();
    const random = Math.random().toString(36).substring(2);
    
    return btoa(`${fingerprint}-${timestamp}-${random}`).replace(/[+/=]/g, '');
  }
  
  // Initialize main app functionality  
  function initializeMainApp() {
    // Initialize camera and app logic, but camera won't start if payment popdown is showing
    initializeCameraAndApp();
  }
  
  // Start the payment check process
  checkInitialPaymentSetup();
  
  // ==================== ORIGINAL APP LOGIC (wrapped in function) ====================
  
  function initializeCameraAndApp() {

  // ==================== PAYMENT INTEGRATION ====================
  
  // Check if user has valid payment method (zero-overhead validation)
  async function checkPaymentMethod() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      // No payment method setup - show popdown
      showPaymentPopdown();
      return false;
    }
    
    try {
      // Validate with server (this ensures user still exists server-side)
      const response = await fetch('/validate-payment', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ device_token: deviceToken })
      });
      
      if (!response.ok) {
        // Server doesn't recognize this device - clear local data and require setup
        localStorage.removeItem('molDeviceToken');
        localStorage.removeItem('molCardInfo');
        showPaymentPopdown();
        return false;
      }
      
      const result = await response.json();
      
      // Update local usage counter with server data
      const localCardInfo = JSON.parse(cardInfo);
      localCardInfo.usage = result.user.usage;
      localStorage.setItem('molCardInfo', JSON.stringify(localCardInfo));
      
      return true;
      
    } catch (error) {
      console.error('Payment validation error:', error);
      // Network error - allow offline usage but show warning
      console.warn('Unable to verify payment method online, allowing offline usage');
      return true;
    }
  }
  
  // Increment usage counter after successful analysis
  async function incrementUsage() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    if (!deviceToken) return;
    
    try {
      const response = await fetch('/increment-usage', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ device_token: deviceToken })
      });
      
      if (response.ok) {
        const result = await response.json();
        
        // Update local counter
        const cardInfo = JSON.parse(localStorage.getItem('molCardInfo') || '{}');
        cardInfo.usage = result.usage;
        localStorage.setItem('molCardInfo', JSON.stringify(cardInfo));
        
        console.log(`📊 Analysis complete - Total usage: ${result.usage}`);
      }
    } catch (error) {
      console.error('Usage increment error:', error);
      // Increment locally as fallback
      const cardInfo = JSON.parse(localStorage.getItem('molCardInfo') || '{}');
      cardInfo.usage = (cardInfo.usage || 0) + 1;
      localStorage.setItem('molCardInfo', JSON.stringify(cardInfo));
    }
  }
  
  // Helper function to check if payment is required
  function isPaymentRequired() {
    const paymentPopdown = document.getElementById('payment-popdown');
    return paymentPopdown && paymentPopdown.classList.contains('show');
  }

  function updateInputMode() {
    // Show/hide camera based on checkbox state
    cameraContainer.style.display = cameraMode.checked ? "flex" : "none";

    // Show/hide photo options based on checkbox state
    photoOptions.style.display = photoMode.checked ? "flex" : "none";
  }

  // Auto-switch checkboxes based on user interaction
  function switchToCameraMode() {
    cameraMode.checked = true;
    photoMode.checked = false;
    updateInputMode();
  }

  function switchToPhotoMode() {
    photoMode.checked = true;
    cameraMode.checked = false;
    updateInputMode();
  }

  function clearModeSelection() {
    cameraMode.checked = false;
    photoMode.checked = false;
    updateInputMode();
  }

  // Event listeners for checkbox changes
  cameraMode.addEventListener("change", updateInputMode);
  photoMode.addEventListener("change", updateInputMode);

  // Auto-switch based on user interaction
  video.addEventListener("click", switchToCameraMode);
  video.addEventListener("touchstart", switchToCameraMode);

  photoUpload.addEventListener("change", switchToPhotoMode);
  photoUrl.addEventListener("focus", switchToPhotoMode);
  urlAnalyze.addEventListener("click", switchToPhotoMode);

  // Text input interaction clears mode selection
  objectInput.addEventListener("focus", clearModeSelection);

  updateInputMode();

  objectInput.addEventListener("keyup", async (e) => {
    if (e.key !== "Enter") return;
    const object = objectInput.value.trim();
    if (!object) return;

    // Show immediate loading feedback
    const loadingColumn = createLoadingColumn(`Analyzing "${object}"...`);

    try {
      const res = await fetch("/object-molecules", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ object }),
      });

      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const { output } = await res.json();

      // Remove loading column
      loadingColumn.remove();
      updateScrollHandles();

      processAnalysisResult(output, snapshots, "Text", object, true);
      
      // Increment usage counter after successful analysis
      incrementUsage();
    } catch (err) {
      // Remove loading column on error
      loadingColumn.remove();
      updateScrollHandles();
      createClosableErrorMessage(`Error analyzing "${object}": ${err.message}`);
    }
    objectInput.value = "";
  });

  const photoUploadHandler = async (e) => {
    const file = e.target.files[0];
    if (!file) return;

    if (!file.type.startsWith("image/")) {
      alert("Please select an image file");
      return;
    }

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
      alert("Please enter an image URL");
      return;
    }

    try {
      new URL(url);
    } catch {
      alert("Please enter a valid URL");
      return;
    }

    analyzeImageFromUrl(url);
  };

  photoUpload.addEventListener("change", photoUploadHandler);

  async function displayUploadedImage(file) {
    const photoOptions = document.getElementById("photo-options");

    photoOptions.innerHTML = "";

    const imageContainer = document.createElement("div");
    imageContainer.className = "uploaded-image-container";

    const img = document.createElement("img");

    // Only show crosshair on mobile devices
    const isMobile =
      /Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(
        navigator.userAgent,
      );
    let crosshair;
    if (isMobile) {
      crosshair = document.createElement("div");
      crosshair.className = "crosshair";
      const beforeLine = document.createElement("div");
      beforeLine.className = "crosshair-line vertical";
      const afterLine = document.createElement("div");
      afterLine.className = "crosshair-line horizontal";
      crosshair.appendChild(beforeLine);
      crosshair.appendChild(afterLine);
    }

    const instructionText = document.createElement("div");
    instructionText.className = "instruction-text";
    instructionText.textContent = isMobile
      ? "Center object in circle & tap, or type name above"
      : "Click on object or type name above";

    const closeButton = document.createElement("button");
    closeButton.innerHTML = '<img src="close.svg" alt="Close" width="16" height="16" />';
    closeButton.className = "close-button";

    closeButton.addEventListener("click", () => {
      const template = document.getElementById("photo-upload-template");
      const clone = template.content.cloneNode(true);
      photoOptions.innerHTML = "";
      photoOptions.appendChild(clone);

      const newPhotoUpload = document.getElementById("photo-upload");
      const newPhotoUrl = document.getElementById("photo-url");
      const newUrlAnalyze = document.getElementById("url-analyze");

      newPhotoUpload.addEventListener("change", photoUploadHandler);
      newPhotoUrl.addEventListener("keyup", photoUrlHandler);
      newUrlAnalyze.addEventListener("click", urlAnalyzeHandler);
    });

    const reader = new FileReader();
    reader.onload = (e) => {
      img.src = e.target.result;
      img.dataset.base64 = e.target.result.split(",")[1];

      imageContainer.addEventListener("click", async (evt) => {
        switchToPhotoMode();
        await handleImageClick(evt, img);
      });

      imageContainer.addEventListener("touchstart", (e) => {
        e.preventDefault();
        switchToPhotoMode();
        handleImageClick(e.touches[0], img);
      });
    };

    reader.readAsDataURL(file);

    imageContainer.appendChild(img);
    if (isMobile && crosshair) imageContainer.appendChild(crosshair);
    imageContainer.appendChild(instructionText);
    imageContainer.appendChild(closeButton);

    photoOptions.appendChild(imageContainer);
  }

  async function handleImageClick(evt, img) {
    // Check payment method before analysis
    if (!checkPaymentMethod()) {
      return;
    }

    // Keep instruction text visible - it's meant to describe the functionality

    const rect = img.getBoundingClientRect();
    const clickX = evt.clientX - rect.left;
    const clickY = evt.clientY - rect.top;

    const relativeX = clickX / rect.width;
    const relativeY = clickY / rect.height;

    const imageBase64 = img.dataset.base64;

    const canvas = document.createElement("canvas");
    const ctx = canvas.getContext("2d");

    const tempImg = new Image();
    tempImg.onload = () => {
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

      cropCtx.drawImage(
        canvas,
        cropX,
        cropY,
        cropSize,
        cropSize,
        0,
        0,
        cropSize,
        cropSize,
      );

      // Calculate the exact middle pixel of the cropped image
      const middleX = Math.floor(cropSize / 2);
      const middleY = Math.floor(cropSize / 2);

      // Draw a visible red square centered on the middle pixel
      // Scale the box size to be visible on screen (about 10% of crop size)
      const boxSize = Math.max(8, Math.floor(cropSize * 0.1));
      cropCtx.save();
      cropCtx.strokeStyle = "#ff0000";
      cropCtx.lineWidth = Math.max(2, Math.floor(cropSize * 0.02)); // Thicker line for visibility
      cropCtx.strokeRect(
        middleX - boxSize / 2,
        middleY - boxSize / 2,
        boxSize,
        boxSize,
      );
      cropCtx.restore();

      const croppedBase64 = cropCanvas
        .toDataURL("image/jpeg", 0.9)
        .split(",")[1];

      // Show immediate loading feedback
      const loadingColumn = createLoadingColumn("Analyzing...", croppedBase64);

      fetch("/image-molecules", {
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
      })
        .then((res) => {
          if (!res.ok) throw new Error(`HTTP ${res.status}`);
          return res.json();
        })
        .then(({ output }) => {
          // Remove loading column
          loadingColumn.remove();
          updateScrollHandles();

          const objectName = output.object || "Uploaded image";
          processAnalysisResult(
            output,
            snapshots,
            "Photo",
            objectName,
            false,
            croppedBase64,
          );
          
          // Increment usage counter after successful analysis
          incrementUsage();
        })
        .catch((err) => {
          // Remove loading column on error
          loadingColumn.remove();
          updateScrollHandles();
          createClosableErrorMessage(`Error: ${err.message}`);
        });
    };

    tempImg.src = `data:image/jpeg;base64,${imageBase64}`;
  }

  async function analyzeImageFromUrl(url) {
    // Check payment method before analysis
    if (!checkPaymentMethod()) {
      return;
    }

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
      const blob = new Blob([byteArray], { type: "image/jpeg" });
      const file = new File([blob], "url-image.jpg", { type: "image/jpeg" });

      // Display the image for interactive clicking
      displayUploadedImage(file);

      // Clear URL input
      photoUrl.value = "";
    } catch (err) {
      createClosableErrorMessage(
        `Error loading image from URL: ${err.message}`,
      );
    }
  }

  urlAnalyze.addEventListener("click", urlAnalyzeHandler);
  photoUrl.addEventListener("keyup", photoUrlHandler);

  function fileToBase64(file) {
    return new Promise((resolve, reject) => {
      const reader = new FileReader();
      reader.onload = () => {
        const result = reader.result.split(",")[1]; // Remove data:image/...;base64, prefix
        resolve(result);
      };
      reader.onerror = reject;
      reader.readAsDataURL(file);
    });
  }

  async function urlToBase64(url) {
    return new Promise((resolve, reject) => {
      const img = new Image();
      img.crossOrigin = "anonymous"; // Try to enable CORS

      img.onload = () => {
        const canvas = document.createElement("canvas");
        const ctx = canvas.getContext("2d");

        canvas.width = img.width;
        canvas.height = img.height;

        ctx.drawImage(img, 0, 0);

        try {
          const dataURL = canvas.toDataURL("image/jpeg", 0.9);
          const base64 = dataURL.split(",")[1];
          resolve(base64);
        } catch (err) {
          reject(new Error("Failed to convert image to base64"));
        }
      };

      img.onerror = () => {
        reject(
          new Error(
            "Failed to load image from URL. Check if URL is valid and accessible.",
          ),
        );
      };

      img.src = url;
    });
  }

  function createLoadingMessage(text) {
    const loadingMsg = document.createElement("h3");
    loadingMsg.textContent = text;
    loadingMsg.style.fontStyle = "italic";
    loadingMsg.style.opacity = "0.7";
    return loadingMsg;
  }

  function createLoadingColumn(loadingText, croppedImageData = null) {
    const gldiv = document.getElementById("gldiv");

    const loadingColumn = document.createElement("div");
    loadingColumn.className = "object-column loading-column";

    const titleContainer = document.createElement("div");
    titleContainer.className = "object-title";

    const titleText = document.createElement("span");
    titleText.textContent = loadingText;
    titleContainer.appendChild(titleText);

    loadingColumn.appendChild(titleContainer);

    // Add cropped image if available
    if (croppedImageData) {
      const croppedImageContainer = document.createElement("div");
      croppedImageContainer.className = "cropped-image-container";

      const croppedImage = document.createElement("img");
      croppedImage.src = `data:image/jpeg;base64,${croppedImageData}`;
      croppedImage.alt = "Cropped region";
      croppedImage.className = "cropped-image";
      croppedImage.style.border = "1px solid #ff0000";

      croppedImageContainer.appendChild(croppedImage);
      loadingColumn.appendChild(croppedImageContainer);
    }

    // Add loading indicator
    const loadingIndicator = document.createElement("div");
    loadingIndicator.className = "loading-indicator loading-message-italic";
    loadingIndicator.textContent = "Processing";
    loadingColumn.appendChild(loadingIndicator);

    gldiv.appendChild(loadingColumn);

    return loadingColumn;
  }

  function createClosableErrorMessage(message, container = snapshots) {
    const errorDiv = document.createElement("div");
    errorDiv.className = "error-message-container";

    errorDiv.innerHTML = `
      <div class="error-message-title">⚠️ Error</div>
      <div>${message}</div>
      <button class="error-close-btn" onclick="this.parentElement.remove(); updateScrollHandles();">×</button>
    `;

    container.appendChild(errorDiv);
    updateScrollHandles();
  }

  function createResultMessage(
    icon,
    objectName,
    smilesCount,
    useQuotes = false,
  ) {
    const name = useQuotes ? `"${objectName}"` : objectName;
    const plural = smilesCount !== 1 ? "s" : "";
    return `${icon} ${name} → ${smilesCount} molecule${plural} found`;
  }

  function processAnalysisResult(
    output,
    container,
    icon,
    objectName,
    useQuotes = false,
    croppedImageData = null,
  ) {
    // Handle new chemicals structure with names and SMILES
    const chemicals = output.chemicals || [];

    // Check if we have a description response
    if (
      chemicals.length === 1 &&
      chemicals[0].smiles &&
      chemicals[0].smiles.startsWith("DESCRIPTION: ")
    ) {
      const description = chemicals[0].smiles.replace("DESCRIPTION: ", "");
      // For description responses, just show the description in the object column header
      generateSDFs([], objectName, description, null, croppedImageData);
      return;
    }

    // Extract SMILES array from chemicals
    const smiles = chemicals.map((chem) => chem.smiles).filter(Boolean);

    // For molecule responses, generate SDFs and show in object column
    if (smiles.length > 0) {
      generateSDFs(smiles, objectName, null, chemicals, croppedImageData);
    }
  }

  if (!navigator.mediaDevices || !navigator.mediaDevices.getUserMedia) {
    console.error("Camera API not available");
    msgBox.hidden = false;
    msgBox.textContent = "Camera API not supported in this browser";
    return;
  }

  const isSecureContext =
    window.isSecureContext || location.protocol === "https:";
  if (!isSecureContext && !location.hostname.includes("localhost")) {
    msgBox.hidden = false;
    msgBox.textContent =
      "Camera requires HTTPS on mobile devices. Please use HTTPS or localhost.";
    return;
  }

  const isMobile =
    /Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(
      navigator.userAgent,
    );
  const isSafari = /^((?!chrome|android).)*safari/i.test(navigator.userAgent);
  const isIOS = /iPad|iPhone|iPod/.test(navigator.userAgent);
  let facingMode = isMobile ? "environment" : "user";
  let currentStream = null;

  const simpleConstraints = () => ({
    video: {
      facingMode,
      width: { ideal: 1280 },
      height: { ideal: 720 },
    },
  });
  const basicConstraints = () => ({
    video: {
      width: { ideal: 1280 },
      height: { ideal: 720 },
    },
  });
  const safariConstraints = () => ({
    video: {
      facingMode,
      width: { min: 640, ideal: 1280, max: 1920 },
      height: { min: 480, ideal: 720, max: 1080 },
    },
  });

  async function startCamera() {
    // Don't start camera if payment setup is required
    if (isPaymentRequired()) {
      console.log('Camera access blocked - payment setup required');
      return;
    }
    
    currentStream?.getTracks().forEach((t) => t.stop());

    try {
      let stream;

      // Safari-specific handling
      if (isSafari || isIOS) {
        try {
          stream =
            await navigator.mediaDevices.getUserMedia(safariConstraints());
        } catch (err) {
          console.log("Safari constraints failed, trying basic:", err);
          stream =
            await navigator.mediaDevices.getUserMedia(basicConstraints());
        }
      } else {
        try {
          stream =
            await navigator.mediaDevices.getUserMedia(simpleConstraints());
        } catch (err) {
          stream =
            await navigator.mediaDevices.getUserMedia(basicConstraints());
        }
      }

      currentStream = stream;
      video.srcObject = stream;

      // Safari-specific video attributes
      video.setAttribute("playsinline", "true");
      video.setAttribute("webkit-playsinline", "true");
      video.setAttribute("x-webkit-airplay", "allow");

      await video.play();
      msgBox.hidden = true;

      // Add mobile targeting reticle if on mobile, remove if on desktop
      if (isMobile) {
        addMobileTargetingReticle();
      } else {
        removeMobileTargetingReticle();
      }
    } catch (err) {
      console.error("Camera error:", err);
      msgBox.hidden = false;

      // Better error messages for Safari
      if (err.name === "NotAllowedError") {
        msgBox.textContent =
          "Camera access denied. Please allow camera access in Safari settings.";
      } else if (err.name === "NotFoundError") {
        msgBox.textContent = "No camera found on this device.";
      } else if (isSafari && err.name === "NotSupportedError") {
        msgBox.textContent =
          "Camera not supported in Safari. Try using Chrome or Firefox.";
      } else {
        msgBox.textContent = `Camera error: ${err.message}`;
      }
    }
  }

  function addMobileTargetingReticle() {
    // Remove existing reticle if any
    const existingReticle = document.querySelector(".mobile-reticle");
    if (existingReticle) {
      existingReticle.remove();
    }

    const reticle = document.createElement("div");
    reticle.className = "mobile-reticle";
    const template = document.getElementById("mobile-reticle-template");
    const clone = template.content.cloneNode(true);
    reticle.appendChild(clone);

    // Append to video element instead of camera container for proper positioning
    video.appendChild(reticle);

    // Add pulsing animation to draw attention
    setTimeout(() => {
      reticle.style.animation = "reticlePulse 2s ease-in-out infinite";
    }, 1000);
  }

  function removeMobileTargetingReticle() {
    const existingReticle = document.querySelector(".mobile-reticle");
    if (existingReticle) {
      existingReticle.remove();
    }
  }

  async function setupSwitchCamera() {
    const devices = await navigator.mediaDevices.enumerateDevices();
    const videoDevices = devices.filter((d) => d.kind === "videoinput");
    console.log("Found video devices:", videoDevices.length);

    if (videoDevices.length > 1) {
      const switchCameraContainer = document.querySelector(
        ".switch-camera-container",
      );
      const switchCameraBtn = document.getElementById("switch-camera-btn");
      if (switchCameraContainer && switchCameraBtn) {
        switchCameraContainer.style.display = "block";
        switchCameraBtn.onclick = () => {
          facingMode = facingMode === "user" ? "environment" : "user";
          startCamera();
        };
      }
    }
  }

  video.addEventListener("click", (e) => {
    switchToCameraMode();
    handleInteraction(e);
  });
  video.addEventListener("touchstart", (e) => {
    e.preventDefault();
    switchToCameraMode();
    handleInteraction(e.touches[0]);
  });

  // Safari-specific touch handling
  if (isSafari || isIOS) {
    video.addEventListener("touchend", (e) => {
      e.preventDefault();
    });

    // Prevent zoom on double tap
    let lastTouchEnd = 0;
    video.addEventListener(
      "touchend",
      (e) => {
        const now = new Date().getTime();
        if (now - lastTouchEnd <= 300) {
          e.preventDefault();
        }
        lastTouchEnd = now;
      },
      false,
    );
  }

  async function handleInteraction(evt) {
    if (!cameraMode.checked) {
      return;
    }

    // Check payment method before analysis
    if (!checkPaymentMethod()) {
      return;
    }

    // On mobile, check if tap is within reticle area
    if (isMobile) {
      const rect = video.getBoundingClientRect();
      const centerX = rect.left + rect.width / 2;
      const centerY = rect.top + rect.height / 2;
      const tapX = evt.clientX;
      const tapY = evt.clientY;

      // Larger reticle area for easier targeting (100px diameter circle)
      const reticleRadius = 50;
      const distanceFromCenter = Math.sqrt(
        Math.pow(tapX - centerX, 2) + Math.pow(tapY - centerY, 2),
      );

      if (distanceFromCenter > reticleRadius) {
        showReticleFeedback();
        return;
      }
    }

    // Show red crop outline at click/tap position
    showCropOutline(evt);

    // Keep instruction text visible - it's meant to describe the functionality

    const canvas = document.createElement("canvas");
    canvas.width = video.videoWidth;
    canvas.height = video.videoHeight;
    canvas.getContext("2d").drawImage(video, 0, 0, canvas.width, canvas.height);
    const imageBase64 = canvas.toDataURL("image/jpeg", 0.9).split(",")[1];

    const scaleX = video.videoWidth / video.clientWidth;
    const scaleY = video.videoHeight / video.clientHeight;
    const clickX = Math.round(evt.clientX - video.getBoundingClientRect().left);
    const clickY = Math.round(evt.clientY - video.getBoundingClientRect().top);
    const actualX = Math.round(clickX * scaleX);
    const actualY = Math.round(clickY * scaleY);

    const crop = document.createElement("canvas");
    crop.width = crop.height = 100;
    const cropCtx = crop.getContext("2d");
    cropCtx.imageSmoothingEnabled = false;
    cropCtx.drawImage(
      canvas,
      actualX - 50,
      actualY - 50,
      100,
      100,
      0,
      0,
      100,
      100,
    );

    // Calculate the exact middle pixel of the cropped image
    const middleX = Math.floor(100 / 2);
    const middleY = Math.floor(100 / 2);

    // Draw a visible red square centered on the middle pixel
    // Scale the box size to be visible on screen (about 10% of crop size)
    const boxSize = Math.max(8, Math.floor(100 * 0.1)); // 10px for 100px crop
    cropCtx.save();
    cropCtx.strokeStyle = "#ff0000";
    cropCtx.lineWidth = Math.max(2, Math.floor(100 * 0.02)); // Thicker line for visibility
    cropCtx.strokeRect(
      middleX - boxSize / 2,
      middleY - boxSize / 2,
      boxSize,
      boxSize,
    );
    cropCtx.restore();

    const croppedBase64 = crop.toDataURL("image/jpeg", 0.9).split(",")[1];

    // Show immediate loading feedback
    const loadingColumn = createLoadingColumn("Analyzing...", croppedBase64);

    try {
      const res = await fetch("/image-molecules", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          imageBase64,
          croppedImageBase64: croppedBase64,
          x: actualX,
          y: actualY,
          cropMiddleX: middleX,
          cropMiddleY: middleY,
          cropSize: 100,
        }),
      });

      if (!res.ok) throw new Error(`HTTP ${res.status}`);

      const { output } = await res.json();

      // Remove loading column
      loadingColumn.remove();
      updateScrollHandles();

      const objectName = output.object || "Unknown object";
      processAnalysisResult(
        output,
        snapshots,
        "Camera",
        objectName,
        false,
        croppedBase64,
      );
      
      // Increment usage counter after successful analysis
      incrementUsage();
    } catch (err) {
      console.error("API call failed:", err);
      // Remove loading column on error
      loadingColumn.remove();
      updateScrollHandles();
      createClosableErrorMessage(`Error: ${err.message}`);
    }
  }

  function showReticleFeedback() {
    const reticle = document.querySelector(".mobile-reticle");
    if (reticle) {
      reticle.style.animation = "reticlePulse 0.5s ease";
      setTimeout(() => {
        reticle.style.animation = "";
      }, 500);
    }
  }

  function showCropOutline(evt) {
    const rect = video.getBoundingClientRect();
    const cropSize = 100;
    const outline = document.createElement("div");
    outline.className = "crop-outline";
    // Only set dynamic positioning styles
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

  function getMoleculeName(chemical) {
    const moleculeNames = {
      O: "Water",
      CCO: "Ethanol",
      "CC(=O)O": "Acetic Acid",
      C: "Methane",
      CO: "Methanol",
      "C(CO)N": "Ethanolamine",
      "C(C(=O)O)N": "Glycine",
      "C(CC(=O)O)N": "GABA",
      "C(CC(=O)O)C(=O)O": "Succinic Acid",
      "C(C(=O)O)O": "Glycolic Acid",
      "C1=CC=CC=C1": "Benzene",
      "CC1=CC=CC=C1": "Toluene",
      "C1=CC=C(C=C1)O": "Phenol",
      "C1=CC=C(C=C1)N": "Aniline",
      "CCCCCCCC=O": "Octanal",
      "CCC(C)C(=O)OC(C)C": "Isopropyl Isovalerate",
      "CCCOC(=O)N1CCCC1C(=O)OCCC": "Polyurethane Monomer",
      C1CCOC1: "Tetrahydrofuran",
      NCCCN: "Propanediamine",
      "CN1C=NC2=C1C(=O)N(C(=O)N2C)C": "Caffeine",
      "CC(=O)OC1=CC=CC=C1C(=O)O": "Aspirin",
      "C1NC(=O)N(C)C2=CC=CC=C12": "N-Methylbenzimidazolone",
      "NC1=NC(=NC(=N1)Cl)Cl": "Cyanuric Chloride",
      "C[C@H](NC(=O)[C@@H](N)Cc1c[nH]c2ccc(Cl)cc12)C(=O)O":
        "Chlorotryptophan Derivative",
      "O=C(Nc1ccc(cc1)C(=O)O)C(F)(F)F": "Trifluoroacetyl-p-aminobenzoic Acid",
      "C(C(C(=O)NC(C(=O)O)C(C)O)O)O": "Threonine Derivative",
      CaCO3: "Calcite (Calcium Carbonate)",
      "CaCO₃": "Calcite (Calcium Carbonate)",
      SiO2: "Quartz (Silicon Dioxide)",
      "SiO₂": "Quartz (Silicon Dioxide)",
      Al2O3: "Corundum (Aluminum Oxide)",
      "Al₂O₃": "Corundum (Aluminum Oxide)",
      FeS2: "Pyrite (Iron Disulfide)",
      "FeS₂": "Pyrite (Iron Disulfide)",
      NaCl: "Halite (Sodium Chloride)",
      quartz: "Quartz (Silicon Dioxide)",
      calcite: "Calcite (Calcium Carbonate)",
      corundum: "Corundum (Aluminum Oxide)",
      pyrite: "Pyrite (Iron Disulfide)",
      halite: "Halite (Sodium Chloride)",
      salt: "Halite (Sodium Chloride)",
    };

    return (
      moleculeNames[chemical] ||
      `Structure (${chemical.substring(0, 20)}${chemical.length > 20 ? "..." : ""})`
    );
  }

  async function generateSDFs(
    smiles,
    objectName,
    description = null,
    chemicals = null,
    croppedImageData = null,
  ) {
    // Handle description responses
    if (description) {
      createObjectColumn(
        objectName,
        [],
        [],
        null,
        null,
        [],
        description,
        chemicals,
        croppedImageData,
      );
      return;
    }

    if (!smiles || smiles.length === 0) return;

    try {
      const response = await fetch("/generate-sdfs", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smiles, overwrite: false }),
      });

      if (!response.ok) {
        throw new Error(`HTTP ${response.status}: ${response.statusText}`);
      }

      const data = await response.json();

      const summary = {
        total: smiles.length,
        visualizable: data.sdfPaths ? data.sdfPaths.length : 0,
        skipped: data.skipped ? data.skipped.length : 0,
        errors: data.errors ? data.errors.length : 0,
      };

      createObjectColumn(
        objectName,
        data.sdfPaths || [],
        smiles,
        null,
        summary,
        data.skipped || [],
        null,
        chemicals,
        croppedImageData,
      );
    } catch (error) {
      console.error("SDF generation error:", error);
      createObjectColumn(
        objectName,
        [],
        smiles,
        "Working on 3D structures...",
        null,
        [],
        null,
        chemicals,
        croppedImageData,
      );
    }
  }

  async function createObjectColumn(
    objectName,
    sdfFiles,
    smiles = [],
    errorMessage = null,
    summary = null,
    skippedChemicals = [],
    description = null,
    chemicals = null,
    croppedImageData = null,
  ) {
    const gldiv = document.getElementById("gldiv");

    const objectColumn = document.createElement("div");
    objectColumn.className = "object-column";

    const titleContainer = document.createElement("div");
    titleContainer.className = "object-title";

    const titleText = document.createElement("span");
    // Show description in header if available, otherwise show object name
    titleText.textContent = description
      ? `${objectName}: ${description}`
      : objectName;
    titleContainer.appendChild(titleText);

    const closeButton = document.createElement("button");
    closeButton.className = "column-close";
    closeButton.innerHTML = '<img src="close.svg" alt="Close" width="16" height="16" />';
    closeButton.onclick = () => {
      objectColumn.remove();
      updateScrollHandles();
    };
    titleContainer.appendChild(closeButton);

    objectColumn.appendChild(titleContainer);

    // Add cropped image if available
    if (croppedImageData) {
      const croppedImageContainer = document.createElement("div");
      croppedImageContainer.className = "cropped-image-container";

      const croppedImage = document.createElement("img");
      croppedImage.src = `data:image/jpeg;base64,${croppedImageData}`;
      croppedImage.alt = "Cropped region";
      croppedImage.className = "cropped-image";
      croppedImage.style.border = "1px solid #ff0000";

      croppedImageContainer.appendChild(croppedImage);
      objectColumn.appendChild(croppedImageContainer);
    }

    if (summary) {
      const template = document.getElementById("chemical-summary-template");
      const clone = template.content.cloneNode(true);
      const summaryDiv = clone.querySelector(".chemical-summary");

      summaryDiv.querySelector(".summary-total").textContent = summary.total;
      summaryDiv.querySelector(".summary-visualizable").textContent =
        summary.visualizable;

      if (summary.skipped > 0) {
        const skippedDiv = summaryDiv.querySelector(".summary-skipped");
        skippedDiv.style.display = "block";
        skippedDiv.querySelector(".summary-skipped-count").textContent =
          summary.skipped;
      }

      if (summary.errors > 0) {
        const errorsDiv = summaryDiv.querySelector(".summary-errors");
        errorsDiv.style.display = "block";
        errorsDiv.querySelector(".summary-errors-count").textContent =
          summary.errors;
      }

      objectColumn.appendChild(summaryDiv);
    }

    if (skippedChemicals.length > 0) {
      const template = document.getElementById("skipped-chemicals-template");
      const clone = template.content.cloneNode(true);
      const skippedDiv = clone.querySelector(".skipped-chemicals");

      skippedDiv.querySelector(".skipped-list").textContent =
        skippedChemicals.join(", ");

      objectColumn.appendChild(skippedDiv);
    }

    if (errorMessage) {
      const errorDiv = document.createElement("div");
      errorDiv.textContent = errorMessage;
      errorDiv.style.color = "#ffffff";
      errorDiv.style.textAlign = "center";
      errorDiv.style.padding = "20px";
      objectColumn.appendChild(errorDiv);
    } else {
      const viewers = [];
      for (let i = 0; i < sdfFiles.length; i++) {
        const sdfFile = sdfFiles[i];
        const moleculeSmiles = smiles[i] || "";

        const moleculeContainer = document.createElement("div");
        moleculeContainer.className = "molecule-container";

        const container = document.createElement("div");
        container.className = "mol-viewer-container";
        moleculeContainer.appendChild(container);

        const moleculeName = document.createElement("div");
        moleculeName.className = "molecule-name";

        // Use chemical name from data if available, otherwise fall back to lookup
        let displayName;
        if (chemicals && chemicals[i] && chemicals[i].name) {
          displayName = chemicals[i].name;
        } else {
          displayName = getMoleculeName(moleculeSmiles);
        }

        // Create clickable Wikipedia link
        const wikipediaLink = document.createElement("a");
        wikipediaLink.textContent = displayName;
        wikipediaLink.href = `https://en.wikipedia.org/wiki/${encodeURIComponent(displayName)}`;
        wikipediaLink.target = "_blank";
        wikipediaLink.rel = "noopener noreferrer";
        wikipediaLink.className = "wikipedia-link";

        moleculeName.appendChild(wikipediaLink);
        moleculeContainer.appendChild(moleculeName);

        objectColumn.appendChild(moleculeContainer);

        const viewer = await render(sdfFile, container);
        if (viewer) viewers.push(viewer);
      }

      viewers.forEach((viewer) => {
        viewer.resize();
        viewer.render();
      });
    }

    gldiv.appendChild(objectColumn);

    // Update scroll handles after adding new column
    updateScrollHandles();
  }

  async function render(sdfFile, container) {
    try {
      console.log(`Attempting to fetch SDF file: ${sdfFile}`);
      const urlParts = sdfFile.split("/");
      const filename = urlParts.pop();
      const encodedFilename = encodeURIComponent(filename);
      const encodedPath = urlParts.join("/") + "/" + encodedFilename;

      console.log(`URL-encoded path: ${encodedPath}`);
      const response = await fetch(encodedPath);
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

      viewer.setBackgroundColor("#000000");

      viewer.setStyle(
        {},
        {
          sphere: {
            scale: 0.8,
          },
        },
      );

      viewer.zoomTo();
      viewer.render();

      return viewer;
    } catch (error) {
      console.error(`Failed to load molecule:`, error);
      container.textContent = `Error loading molecule: ${error.message}`;
      container.style.color = "red";
      container.style.textAlign = "center";
      container.style.padding = "20px";
      return null;
    }
  }

  (async () => {
    try {
      await startCamera();
      await setupSwitchCamera();

      // Ensure reticle is only shown on mobile
      if (!isMobile) {
        removeMobileTargetingReticle();
      }
    } catch (err) {
      console.error("Camera setup failed:", err);
    }
  })();

  // Add keyboard navigation for 3D viewer scrolling
  document.addEventListener("keydown", (e) => {
    const gldiv = document.getElementById("gldiv");
    if (!gldiv) return;

    const scrollAmount = 450; // Width of one column + gap

    switch (e.key) {
      case "ArrowLeft":
        e.preventDefault();
        gldiv.scrollLeft -= scrollAmount;
        updateScrollHandles();
        break;
      case "ArrowRight":
        e.preventDefault();
        gldiv.scrollLeft += scrollAmount;
        updateScrollHandles();
        break;
      case "Home":
        e.preventDefault();
        gldiv.scrollLeft = 0;
        updateScrollHandles();
        break;
      case "End":
        e.preventDefault();
        gldiv.scrollLeft = gldiv.scrollWidth;
        updateScrollHandles();
        break;
    }
  });

  // Create scroll handle buttons
  function createScrollHandles() {
    const gldiv = document.getElementById("gldiv");
    if (!gldiv) return;

    // Create left scroll button
    const leftHandle = document.createElement("button");
    leftHandle.className = "scroll-handle scroll-handle-left";
    leftHandle.innerHTML = '<img src="chevron-left.svg" alt="Scroll left" width="16" height="16" />';
    leftHandle.setAttribute("aria-label", "Scroll left");
    leftHandle.setAttribute("title", "Scroll left (or press Left arrow)");
    leftHandle.id = "scroll-left-btn";

    // Create right scroll button
    const rightHandle = document.createElement("button");
    rightHandle.className = "scroll-handle scroll-handle-right";
    rightHandle.innerHTML = '<img src="chevron-right.svg" alt="Scroll right" width="16" height="16" />';
    rightHandle.setAttribute("aria-label", "Scroll right");
    rightHandle.setAttribute("title", "Scroll right (or press Right arrow)");
    rightHandle.id = "scroll-right-btn";

    // Add click handlers with throttling
    let clickThrottle = false;
    leftHandle.addEventListener("click", (e) => {
      e.preventDefault();
      e.stopPropagation();
      if (clickThrottle) return;
      clickThrottle = true;
      setTimeout(() => (clickThrottle = false), 100);

      const scrollAmount = 450;
      gldiv.scrollLeft -= scrollAmount;
      updateScrollHandles();
    });

    rightHandle.addEventListener("click", (e) => {
      e.preventDefault();
      e.stopPropagation();
      if (clickThrottle) return;
      clickThrottle = true;
      setTimeout(() => (clickThrottle = false), 100);

      const scrollAmount = 450;
      gldiv.scrollLeft += scrollAmount;
      updateScrollHandles();
    });

    // Add keyboard support
    leftHandle.addEventListener("keydown", (e) => {
      if (e.key === "Enter" || e.key === " ") {
        e.preventDefault();
        leftHandle.click();
      }
    });

    rightHandle.addEventListener("keydown", (e) => {
      if (e.key === "Enter" || e.key === " ") {
        e.preventDefault();
        rightHandle.click();
      }
    });

    // Append to gldiv
    gldiv.appendChild(leftHandle);
    gldiv.appendChild(rightHandle);
  }

  // Update scroll handle visibility
  let updateScrollHandlesTimeout;
  function updateScrollHandles() {
    // Debounce rapid calls
    clearTimeout(updateScrollHandlesTimeout);
    updateScrollHandlesTimeout = setTimeout(() => {
      updateScrollHandlesInternal();
    }, 10);
  }

  function updateScrollHandlesInternal() {
    const gldiv = document.getElementById("gldiv");
    if (!gldiv) return;

    const leftHandle = document.getElementById("scroll-left-btn");
    const rightHandle = document.getElementById("scroll-right-btn");
    if (!leftHandle || !rightHandle) return;

    // Avoid unnecessary updates if scroll position hasn't changed
    const currentScrollLeft = gldiv.scrollLeft;
    const currentScrollWidth = gldiv.scrollWidth;
    const currentClientWidth = gldiv.clientWidth;

    if (
      updateScrollHandlesInternal.lastScrollLeft === currentScrollLeft &&
      updateScrollHandlesInternal.lastScrollWidth === currentScrollWidth &&
      updateScrollHandlesInternal.lastClientWidth === currentClientWidth
    ) {
      return;
    }

    updateScrollHandlesInternal.lastScrollLeft = currentScrollLeft;
    updateScrollHandlesInternal.lastScrollWidth = currentScrollWidth;
    updateScrollHandlesInternal.lastClientWidth = currentClientWidth;

    const canScrollLeft = gldiv.scrollLeft > 0;
    const canScrollRight =
      gldiv.scrollLeft < gldiv.scrollWidth - gldiv.clientWidth;

    // Only show handles if content overflows
    const hasOverflow = gldiv.scrollWidth > gldiv.clientWidth;

    const shouldShowLeft = hasOverflow && canScrollLeft;
    const shouldShowRight = hasOverflow && canScrollRight;

    // Only update if visibility state actually changes
    if (shouldShowLeft !== leftHandle.classList.contains("show")) {
      if (shouldShowLeft) {
        leftHandle.classList.add("show");
        leftHandle.setAttribute("tabindex", "0");
      } else {
        leftHandle.classList.remove("show");
        leftHandle.setAttribute("tabindex", "-1");
      }
    }

    if (shouldShowRight !== rightHandle.classList.contains("show")) {
      if (shouldShowRight) {
        rightHandle.classList.add("show");
        rightHandle.setAttribute("tabindex", "0");
      } else {
        rightHandle.classList.remove("show");
        rightHandle.setAttribute("tabindex", "-1");
      }
    }
  }

  // Create scroll handles after DOM is ready
  createScrollHandles();

  // Update scroll handles on scroll (throttled)
  let scrollTimeout;
  document.addEventListener(
    "scroll",
    (e) => {
      if (e.target.id === "gldiv") {
        clearTimeout(scrollTimeout);
        scrollTimeout = setTimeout(updateScrollHandles, 16); // ~60fps
      }
    },
    true,
  );

  // Update scroll handles on window resize (throttled)
  let resizeTimeout;
  window.addEventListener("resize", () => {
    clearTimeout(resizeTimeout);
    resizeTimeout = setTimeout(updateScrollHandles, 100); // Throttle resize events
  });

  // Initial update
  updateScrollHandles();
  
  } // End of initializeCameraAndApp function

  // CAMERA PERMISSIONS: Safari-specific improvements
  window.requestCameraPermission = async function() {
    try {
      console.log('🎥 Requesting camera permission...');
      
      // Request permission explicitly
      const stream = await navigator.mediaDevices.getUserMedia({ 
        video: { facingMode: 'environment' } 
      });
      
      // Stop the stream immediately after getting permission
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
  };

  // Auto-request camera permission on load (for Safari)
  if (/^((?!chrome|android).)*safari/i.test(navigator.userAgent)) {
    console.log('🍎 Safari detected - requesting camera permission early');
    setTimeout(() => {
      window.requestCameraPermission();
    }, 1000);
  }

});
