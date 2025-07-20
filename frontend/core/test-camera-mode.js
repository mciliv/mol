// Simple test to verify camera mode activation
// Run this in the browser console

console.log('🧪 Testing Camera Mode Activation');

// Test 1: Check if elements exist
function testElements() {
  console.log('\n🔍 Test 1: Check Elements');
  
  const elements = {
    'Camera Mode Checkbox': document.getElementById('camera-mode'),
    'Camera Container': document.getElementById('camera-container'),
    'Photo Mode Checkbox': document.getElementById('photo-mode'),
    'Photo Options': document.getElementById('photo-options')
  };
  
  Object.entries(elements).forEach(([name, element]) => {
    console.log(`${element ? '✅' : '❌'} ${name}: ${element ? 'Found' : 'Missing'}`);
  });
  
  return Object.values(elements).every(Boolean);
}

// Test 2: Test camera mode activation
function testCameraMode() {
  console.log('\n📷 Test 2: Camera Mode Activation');
  
  const cameraMode = document.getElementById('camera-mode');
  const cameraContainer = document.getElementById('camera-container');
  
  if (!cameraMode || !cameraContainer) {
    console.log('❌ Required elements not found');
    return false;
  }
  
  // Check initial state
  console.log(`Initial state - Checked: ${cameraMode.checked}, Active class: ${cameraContainer.classList.contains('active')}`);
  
  // Activate camera mode
  cameraMode.checked = true;
  cameraMode.dispatchEvent(new Event('change'));
  
  // Check after activation
  setTimeout(() => {
    console.log(`After activation - Checked: ${cameraMode.checked}, Active class: ${cameraContainer.classList.contains('active')}`);
    console.log(`Camera container display: ${window.getComputedStyle(cameraContainer).display}`);
    
    if (cameraContainer.classList.contains('active')) {
      console.log('✅ Camera mode activated successfully!');
    } else {
      console.log('❌ Camera mode not activated');
    }
  }, 100);
  
  return true;
}

// Test 3: Test photo mode activation
function testPhotoMode() {
  console.log('\n📸 Test 3: Photo Mode Activation');
  
  const photoMode = document.getElementById('photo-mode');
  const photoOptions = document.getElementById('photo-options');
  
  if (!photoMode || !photoOptions) {
    console.log('❌ Required elements not found');
    return false;
  }
  
  // Check initial state
  console.log(`Initial state - Checked: ${photoMode.checked}, Active class: ${photoOptions.classList.contains('active')}`);
  
  // Activate photo mode
  photoMode.checked = true;
  photoMode.dispatchEvent(new Event('change'));
  
  // Check after activation
  setTimeout(() => {
    console.log(`After activation - Checked: ${photoMode.checked}, Active class: ${photoOptions.classList.contains('active')}`);
    console.log(`Photo options display: ${window.getComputedStyle(photoOptions).display}`);
    
    if (photoOptions.classList.contains('active')) {
      console.log('✅ Photo mode activated successfully!');
    } else {
      console.log('❌ Photo mode not activated');
    }
  }, 100);
  
  return true;
}

// Test 4: Test mode switching
function testModeSwitching() {
  console.log('\n🔄 Test 4: Mode Switching');
  
  const cameraMode = document.getElementById('camera-mode');
  const photoMode = document.getElementById('photo-mode');
  const cameraContainer = document.getElementById('camera-container');
  const photoOptions = document.getElementById('photo-options');
  
  if (!cameraMode || !photoMode || !cameraContainer || !photoOptions) {
    console.log('❌ Required elements not found');
    return false;
  }
  
  // Start with camera mode
  cameraMode.checked = true;
  photoMode.checked = false;
  cameraMode.dispatchEvent(new Event('change'));
  photoMode.dispatchEvent(new Event('change'));
  
  setTimeout(() => {
    console.log('Camera mode active:', cameraContainer.classList.contains('active'));
    console.log('Photo mode active:', photoOptions.classList.contains('active'));
    
    // Switch to photo mode
    cameraMode.checked = false;
    photoMode.checked = true;
    cameraMode.dispatchEvent(new Event('change'));
    photoMode.dispatchEvent(new Event('change'));
    
    setTimeout(() => {
      console.log('After switch - Camera mode active:', cameraContainer.classList.contains('active'));
      console.log('After switch - Photo mode active:', photoOptions.classList.contains('active'));
      
      if (!cameraContainer.classList.contains('active') && photoOptions.classList.contains('active')) {
        console.log('✅ Mode switching works correctly!');
      } else {
        console.log('❌ Mode switching not working');
      }
    }, 100);
  }, 100);
  
  return true;
}

// Run all tests
function runTests() {
  console.log('🚀 Running Camera Mode Tests...\n');
  
  const tests = [
    testElements,
    testCameraMode,
    testPhotoMode,
    testModeSwitching
  ];
  
  tests.forEach((test, index) => {
    try {
      test();
    } catch (error) {
      console.log(`❌ Test ${index + 1} failed:`, error.message);
    }
  });
  
  console.log('\n📝 Manual Test Instructions:');
  console.log('1. Click the camera checkbox to activate camera mode');
  console.log('2. Verify camera view appears');
  console.log('3. Click the photo checkbox to switch to photo mode');
  console.log('4. Verify photo upload options appear');
}

// Export for manual testing
window.cameraModeTests = {
  runTests,
  testElements,
  testCameraMode,
  testPhotoMode,
  testModeSwitching
};

console.log('✅ Camera mode testing functions loaded. Run cameraModeTests.runTests() to test.'); 