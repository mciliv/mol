// Test script to verify camera click triggers analysis
// Run this in the browser console

console.log('🧪 Testing Camera Click → Analysis Connection');

// Test 1: Check if camera elements are properly set up
function testCameraSetup() {
  console.log('\n🔍 Test 1: Camera Setup');
  
  const elements = {
    'Video Element': document.getElementById('video-feed'),
    'Camera Manager': window.cameraManager,
    'UI Manager': window.uiManager,
    'Payment Manager': window.paymentManager
  };
  
  Object.entries(elements).forEach(([name, element]) => {
    console.log(`${element ? '✅' : '❌'} ${name}: ${element ? 'Available' : 'Missing'}`);
  });
  
  return Object.values(elements).every(Boolean);
}

// Test 2: Test camera click event binding
function testCameraClickBinding() {
  console.log('\n👆 Test 2: Camera Click Event Binding');
  
  const video = document.getElementById('video-feed');
  if (!video) {
    console.log('❌ Video element not found');
    return false;
  }
  
  // Check if click event is bound
  const events = getEventListeners ? getEventListeners(video) : null;
  if (events && events.click) {
    console.log(`✅ Click event bound: ${events.click.length} listener(s)`);
  } else {
    console.log('⚠️ Click event listeners not visible (normal in production)');
  }
  
  // Test click simulation
  console.log('Simulating camera click...');
  video.click();
  
  return true;
}

// Test 3: Test camera capture functionality
function testCameraCapture() {
  console.log('\n📸 Test 3: Camera Capture');
  
  if (!window.cameraManager) {
    console.log('❌ Camera manager not available');
    return false;
  }
  
  // Check if camera is initialized
  if (window.cameraManager.video && window.cameraManager.video.srcObject) {
    console.log('✅ Camera stream active');
  } else {
    console.log('⚠️ Camera stream not active - may need permission');
  }
  
  return true;
}

// Test 4: Test analysis API endpoint
async function testAnalysisAPI() {
  console.log('\n🌐 Test 4: Analysis API');
  
  try {
    const response = await fetch('/image-molecules', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        imageBase64: 'test',
        croppedImageBase64: 'test',
        x: 100,
        y: 100,
        cropMiddleX: 50,
        cropMiddleY: 50,
        cropSize: 100
      })
    });
    
    console.log(`✅ API endpoint responds: ${response.status}`);
    return true;
  } catch (error) {
    console.log(`❌ API endpoint error: ${error.message}`);
    return false;
  }
}

// Test 5: Test complete camera click flow
function testCompleteFlow() {
  console.log('\n🔄 Test 5: Complete Camera Click Flow');
  
  const video = document.getElementById('video-feed');
  if (!video) {
    console.log('❌ Video element not found');
    return false;
  }
  
  // Listen for analysis completion event
  const analysisListener = (event) => {
    console.log('✅ Analysis completed!', event.detail);
    document.removeEventListener('imageAnalysisComplete', analysisListener);
  };
  
  document.addEventListener('imageAnalysisComplete', analysisListener);
  
  // Simulate camera click
  console.log('Simulating camera click to trigger analysis...');
  video.click();
  
  // Clean up listener after 10 seconds
  setTimeout(() => {
    document.removeEventListener('imageAnalysisComplete', analysisListener);
    console.log('⚠️ Analysis event not received within 10 seconds');
  }, 10000);
  
  return true;
}

// Run all tests
function runCameraClickTests() {
  console.log('🚀 Running Camera Click Tests...\n');
  
  const tests = [
    testCameraSetup,
    testCameraClickBinding,
    testCameraCapture,
    testAnalysisAPI,
    testCompleteFlow
  ];
  
  tests.forEach((test, index) => {
    try {
      test();
    } catch (error) {
      console.log(`❌ Test ${index + 1} failed:`, error.message);
    }
  });
  
  console.log('\n📝 Manual Test Instructions:');
  console.log('1. Make sure camera mode is activated (checkbox checked)');
  console.log('2. Grant camera permission if prompted');
  console.log('3. Click on an object in the camera view');
  console.log('4. Verify analysis starts and results appear');
  console.log('5. Check browser console for any errors');
}

// Export for manual testing
window.cameraClickTests = {
  runCameraClickTests,
  testCameraSetup,
  testCameraClickBinding,
  testCameraCapture,
  testAnalysisAPI,
  testCompleteFlow
};

console.log('✅ Camera click testing functions loaded. Run cameraClickTests.runCameraClickTests() to test.'); 