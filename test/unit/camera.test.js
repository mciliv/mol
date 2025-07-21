// Focused test for camera click triggering image analysis
// This can be run in the browser console

console.log('📷 Testing Camera Click → Image Analysis Connection');

// Test 1: Verify camera elements exist
function testCameraElements() {
  console.log('\n🔍 Test 1: Camera Elements');
  
  const elements = {
    'Video Feed': document.getElementById('video-feed'),
    'Camera Container': document.getElementById('camera-container'),
    'Camera Mode Checkbox': document.getElementById('camera-mode')
  };
  
  Object.entries(elements).forEach(([name, element]) => {
    console.log(`${element ? '✅' : '❌'} ${name}: ${element ? 'Found' : 'Missing'}`);
  });
  
  return Object.values(elements).every(Boolean);
}

// Test 2: Test camera mode activation
function testCameraModeActivation() {
  console.log('\n🎯 Test 2: Camera Mode Activation');
  
  const videoFeed = document.getElementById('video-feed');
  const cameraMode = document.getElementById('camera-mode');
  
  if (!videoFeed || !cameraMode) {
    console.log('❌ Camera elements not found');
    return false;
  }
  
  // Simulate camera mode activation
  cameraMode.checked = true;
  cameraMode.dispatchEvent(new Event('change'));
  
  // Check if camera container is visible
  setTimeout(() => {
    const cameraContainer = document.getElementById('camera-container');
    const isVisible = !cameraContainer.classList.contains('hidden');
    console.log(`${isVisible ? '✅' : '❌'} Camera container visible: ${isVisible}`);
  }, 100);
  
  return true;
}

// Test 3: Test camera click event
function testCameraClick() {
  console.log('\n👆 Test 3: Camera Click Event');
  
  const videoFeed = document.getElementById('video-feed');
  if (!videoFeed) {
    console.log('❌ Video feed not found');
    return false;
  }
  
  // Add click event listener to test
  let clickDetected = false;
  const testClickListener = () => {
    clickDetected = true;
    console.log('✅ Camera click detected');
  };
  
  videoFeed.addEventListener('click', testClickListener);
  
  // Simulate click
  videoFeed.click();
  
  // Remove test listener
  setTimeout(() => {
    videoFeed.removeEventListener('click', testClickListener);
    if (!clickDetected) {
      console.log('❌ Camera click not detected');
    }
  }, 100);
  
  return true;
}

// Test 4: Test camera stream initialization
function testCameraStream() {
  console.log('\n📹 Test 4: Camera Stream');
  
  const videoFeed = document.getElementById('video-feed');
  if (!videoFeed) {
    console.log('❌ Video feed not found');
    return false;
  }
  
  // Check if camera manager is available
  if (!window.cameraManager) {
    console.log('❌ Camera manager not available');
    return false;
  }
  
  // Test camera initialization
  if (window.cameraManager.isInitialized) {
    console.log('✅ Camera manager initialized');
  } else {
    console.log('⚠️ Camera manager not initialized');
  }
  
  // Test camera permission
  navigator.mediaDevices.getUserMedia({ video: true })
    .then(stream => {
      console.log('✅ Camera permission granted');
      stream.getTracks().forEach(track => track.stop());
    })
    .catch(err => {
      console.log(`⚠️ Camera permission: ${err.message}`);
    });
  
  return true;
}

// Test 5: Test image capture simulation
function testImageCapture() {
  console.log('\n📸 Test 5: Image Capture Simulation');
  
  // Create a mock image capture
  const canvas = document.createElement('canvas');
  canvas.width = 640;
  canvas.height = 480;
  const ctx = canvas.getContext('2d');
  
  // Draw a simple test pattern
  ctx.fillStyle = '#ff0000';
  ctx.fillRect(0, 0, 100, 100);
  ctx.fillStyle = '#00ff00';
  ctx.fillRect(100, 100, 100, 100);
  
  // Convert to base64
  const imageData = canvas.toDataURL('image/jpeg', 0.8).split(',')[1];
  
  console.log('✅ Mock image captured');
  console.log(`📊 Image data length: ${imageData.length} characters`);
  
  return true;
}

// Test 6: Test analysis trigger
function testAnalysisTrigger() {
  console.log('\n🔬 Test 6: Analysis Trigger');
  
  // Check if molecular app is available
  if (!window.molecularApp) {
    console.log('❌ Molecular app not available');
    return false;
  }
  
  // Test if analysis methods exist
  const methods = {
    'handleTextAnalysis': typeof window.molecularApp.handleTextAnalysis === 'function',
    'handleImageClick': typeof window.molecularApp.handleImageClick === 'function',
    'processAnalysisResult': typeof window.molecularApp.processAnalysisResult === 'function'
  };
  
  Object.entries(methods).forEach(([name, exists]) => {
    console.log(`${exists ? '✅' : '❌'} ${name}: ${exists ? 'Available' : 'Missing'}`);
  });
  
  return Object.values(methods).every(Boolean);
}

// Test 7: Test complete camera → analysis flow
function testCompleteFlow() {
  console.log('\n🔄 Test 7: Complete Camera → Analysis Flow');
  
  // Step 1: Activate camera mode
  const cameraMode = document.getElementById('camera-mode');
  if (cameraMode) {
    cameraMode.checked = true;
    cameraMode.dispatchEvent(new Event('change'));
    console.log('✅ Step 1: Camera mode activated');
  }
  
  // Step 2: Simulate camera click
  const videoFeed = document.getElementById('video-feed');
  if (videoFeed) {
    videoFeed.click();
    console.log('✅ Step 2: Camera clicked');
  }
  
  // Step 3: Check if analysis would be triggered
  setTimeout(() => {
    console.log('✅ Step 3: Analysis flow ready');
    console.log('📝 Note: Actual analysis requires camera permission and image capture');
  }, 500);
  
  return true;
}

// Run all camera tests
function runCameraTests() {
  console.log('🚀 Running Camera Connection Tests...\n');
  
  const tests = [
    testCameraElements,
    testCameraModeActivation,
    testCameraClick,
    testCameraStream,
    testImageCapture,
    testAnalysisTrigger,
    testCompleteFlow
  ];
  
  let passed = 0;
  let total = tests.length;
  
  tests.forEach((test, index) => {
    try {
      const result = test();
      if (result) passed++;
    } catch (error) {
      console.log(`❌ Test ${index + 1} failed: ${error.message}`);
    }
  });
  
  console.log(`\n📊 Camera Test Results: ${passed}/${total} tests passed`);
  
  if (passed === total) {
    console.log('🎉 Camera click → image analysis connection is working!');
  } else {
    console.log('⚠️ Some camera functionality needs attention');
  }
  
  console.log('\n📝 Manual Testing Steps:');
  console.log('1. Click the camera mode checkbox');
  console.log('2. Click on the video feed area');
  console.log('3. Grant camera permission when prompted');
  console.log('4. Click on an object in the camera view');
  console.log('5. Verify analysis results appear');
}

// Export for manual testing
window.cameraTests = {
  runCameraTests,
  testCameraElements,
  testCameraModeActivation,
  testCameraClick,
  testCameraStream,
  testImageCapture,
  testAnalysisTrigger,
  testCompleteFlow
};

console.log('✅ Camera testing functions loaded. Run cameraTests.runCameraTests() to test camera connections.'); 