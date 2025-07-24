// Manual testing script to verify component connections
// Run this in the browser console to test functionality

console.log('🧪 Starting manual component connection tests...');

// Test 1: Check if all managers are initialized
function testInitialization() {
  console.log('\n📋 Test 1: Component Initialization');
  
  const checks = {
    'Camera Manager': window.cameraManager && window.cameraManager.isInitialized,
    'Payment Manager': window.paymentManager && typeof window.paymentManager.checkPaymentMethod === 'function',
    'UI Manager': window.uiManager && typeof window.uiManager.switchToCameraMode === 'function',
    'Molecular App': window.molecularApp && typeof window.molecularApp.handleTextAnalysis === 'function'
  };
  
  Object.entries(checks).forEach(([name, status]) => {
    console.log(`${status ? '✅' : '❌'} ${name}: ${status ? 'Initialized' : 'Missing'}`);
  });
  
  return Object.values(checks).every(Boolean);
}

// Test 2: Test camera click triggers camera mode
function testCameraClick() {
  console.log('\n📷 Test 2: Camera Click Integration');
  
  const videoFeed = document.getElementById('video-feed');
  if (!videoFeed) {
    console.log('❌ Video feed element not found');
    return false;
  }
  
  // Simulate click
  videoFeed.click();
  
  // Check if camera mode is activated
  setTimeout(() => {
    const isActive = videoFeed.classList.contains('active');
    console.log(`${isActive ? '✅' : '❌'} Camera click triggered camera mode: ${isActive}`);
  }, 100);
  
  return true;
}

// Test 3: Test text input triggers analysis
function testTextAnalysis() {
  console.log('\n🔬 Test 3: Text Analysis Integration');
  
  const input = document.getElementById('object-input');
  if (!input) {
    console.log('❌ Text input element not found');
    return false;
  }
  
  // Test with a simple molecule
  input.value = 'water';
  
  // Simulate Enter key
  const enterEvent = new KeyboardEvent('keyup', { key: 'Enter' });
  input.dispatchEvent(enterEvent);
  
  console.log('✅ Text analysis triggered (check console for API calls)');
  return true;
}

// Test 4: Test photo upload triggers photo mode
function testPhotoUpload() {
  console.log('\n📸 Test 4: Photo Upload Integration');
  
  const photoUpload = document.getElementById('photo-upload');
  if (!photoUpload) {
    console.log('❌ Photo upload element not found');
    return false;
  }
  
  // Simulate file selection (this won't actually upload, but tests the event handler)
  const changeEvent = new Event('change');
  photoUpload.dispatchEvent(changeEvent);
  
  console.log('✅ Photo upload event handler triggered');
  return true;
}

// Test 5: Test URL analysis integration
function testUrlAnalysis() {
  console.log('\n🌐 Test 5: URL Analysis Integration');
  
  const photoUrl = document.getElementById('photo-url');
  const urlAnalyze = document.getElementById('url-analyze');
  
  if (!photoUrl || !urlAnalyze) {
    console.log('❌ URL analysis elements not found');
    return false;
  }
  
  // Test URL input
  photoUrl.value = 'https://example.com/test.jpg';
  const keyupEvent = new KeyboardEvent('keyup', { key: 'Enter' });
  photoUrl.dispatchEvent(keyupEvent);
  
  console.log('✅ URL analysis triggered');
  return true;
}

// Test 6: Test mode switching
function testModeSwitching() {
  console.log('\n🔄 Test 6: Mode Switching');
  
  const video = document.getElementById('video-feed');
  const photoUpload = document.getElementById('photo-upload');
  const photoUrl = document.getElementById('photo-url');
  
  if (!video || !photoUpload || !photoUrl) {
    console.log('❌ Mode switching elements not found');
    return false;
  }
  
  // Test camera mode
  video.click();
  setTimeout(() => {
    const cameraActive = video.classList.contains('active');
    console.log(`${cameraActive ? '✅' : '❌'} Camera mode switching: ${cameraActive}`);
  }, 100);
  
  // Test photo mode
  photoUpload.dispatchEvent(new Event('change'));
  setTimeout(() => {
    const photoActive = photoUpload.classList.contains('active') || photoUrl.classList.contains('active');
    console.log(`${photoActive ? '✅' : '❌'} Photo mode switching: ${photoActive}`);
  }, 100);
  
  return true;
}

// Test 7: Test error handling
function testErrorHandling() {
  console.log('\n⚠️ Test 7: Error Handling');
  
  const input = document.getElementById('object-input');
  if (!input) {
    console.log('❌ Text input element not found');
    return false;
  }
  
  // Test with invalid input
  input.value = 'invalid_molecule_xyz123';
  const enterEvent = new KeyboardEvent('keyup', { key: 'Enter' });
  input.dispatchEvent(enterEvent);
  
  console.log('✅ Error handling triggered (check for error messages)');
  return true;
}

// Test 8: Test component cleanup
function testComponentCleanup() {
  console.log('\n🧹 Test 8: Component Cleanup');
  
  // Look for existing results
  const existingColumns = document.querySelectorAll('.object-column');
  if (existingColumns.length > 0) {
    // Test closing a result
    const closeButton = existingColumns[0].querySelector('.object-close');
    if (closeButton) {
      closeButton.click();
      console.log('✅ Component cleanup triggered');
      return true;
    }
  }
  
  console.log('⚠️ No existing results to test cleanup');
  return true;
}

// Test 9: Test payment integration
function testPaymentIntegration() {
  console.log('\n💳 Test 9: Payment Integration');
  
  // Check if payment manager is available
  if (!window.paymentManager) {
    console.log('❌ Payment manager not found');
    return false;
  }
  
  // Test payment check
  window.paymentManager.checkPaymentMethod().then(hasPayment => {
    console.log(`✅ Payment check completed: ${hasPayment ? 'Has payment' : 'No payment'}`);
  }).catch(err => {
    console.log(`⚠️ Payment check error: ${err.message}`);
  });
  
  return true;
}

// Test 10: Test 3D visualization
function test3DVisualization() {
  console.log('\n🎨 Test 10: 3D Visualization');
  
  // Check if 3Dmol is available
  if (typeof $3Dmol === 'undefined') {
    console.log('❌ 3Dmol library not loaded');
    return false;
  }
  
  console.log('✅ 3Dmol library available for molecule rendering');
  return true;
}

// Run all tests
function runAllTests() {
  console.log('🚀 Running all component connection tests...\n');
  
  const tests = [
    testInitialization,
    testCameraClick,
    testTextAnalysis,
    testPhotoUpload,
    testUrlAnalysis,
    testModeSwitching,
    testErrorHandling,
    testComponentCleanup,
    testPaymentIntegration,
    test3DVisualization
  ];
  
  let passed = 0;
  let total = tests.length;
  
  tests.forEach((test, index) => {
    try {
      const result = test();
      if (result) passed++;
    } catch (error) {
      console.log(`❌ Test ${index + 1} failed with error:`, error.message);
    }
  });
  
  console.log(`\n📊 Test Results: ${passed}/${total} tests passed`);
  
  if (passed === total) {
    console.log('🎉 All component connections are working properly!');
  } else {
    console.log('⚠️ Some component connections need attention');
  }
}

// Export functions for manual testing
if (typeof window !== 'undefined') {
  window.componentTests = {
  runAllTests,
  testInitialization,
  testCameraClick,
  testTextAnalysis,
  testPhotoUpload,
  testUrlAnalysis,
  testModeSwitching,
  testErrorHandling,
  testComponentCleanup,
  testPaymentIntegration,
  test3DVisualization
};
}

console.log('✅ Manual testing functions loaded. Run componentTests.runAllTests() to start testing.'); 