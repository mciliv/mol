// Test camera click functionality and debug image analysis issues
// Run this in browser console to diagnose problems

console.log('🔍 Camera Click Analysis Debug Tool');

// Test 1: Check if camera manager is properly initialized
function testCameraManager() {
  console.log('\n📷 Test 1: Camera Manager Status');
  
  const video = document.getElementById('video-feed');
  if (!video) {
    console.log('❌ Video element not found');
    return false;
  }
  
  console.log('✅ Video element found:', video);
  console.log('📊 Video properties:', {
    videoWidth: video.videoWidth,
    videoHeight: video.videoHeight,
    clientWidth: video.clientWidth,
    clientHeight: video.clientHeight,
    srcObject: !!video.srcObject,
    paused: video.paused
  });
  
  // Check if camera manager exists
  if (window.cameraManager) {
    console.log('✅ Camera manager found');
    console.log('📊 Camera manager state:', {
      isMobile: window.cameraManager.isMobile,
      isSafari: window.cameraManager.isSafari,
      isIOS: window.cameraManager.isIOS,
      currentStream: !!window.cameraManager.currentStream
    });
  } else {
    console.log('❌ Camera manager not found');
  }
  
  return true;
}

// Test 2: Check payment manager status
function testPaymentManager() {
  console.log('\n💳 Test 2: Payment Manager Status');
  
  if (window.paymentManager) {
    console.log('✅ Payment manager found');
    
    // Check local storage
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    console.log('📊 Payment status:', {
      deviceToken: !!deviceToken,
      cardInfo: !!cardInfo,
      hasPaymentSetup: !!(deviceToken && cardInfo)
    });
    
    // Test payment check
    window.paymentManager.checkPaymentMethod().then(hasPayment => {
      console.log('✅ Payment check result:', hasPayment);
    }).catch(err => {
      console.log('❌ Payment check error:', err);
    });
    
  } else {
    console.log('❌ Payment manager not found');
  }
  
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

// Test 4: Test image analysis event
function testImageAnalysis() {
  console.log('\n🔬 Test 4: Image Analysis Event');
  
  // Listen for analysis completion event
  const analysisListener = (event) => {
    console.log('✅ Image analysis completed:', event.detail);
    document.removeEventListener('imageAnalysisComplete', analysisListener);
  };
  
  document.addEventListener('imageAnalysisComplete', analysisListener);
  
  // Simulate camera click to trigger analysis
  const video = document.getElementById('video-feed');
  if (video) {
    video.click();
    console.log('✅ Camera click triggered for analysis test');
  }
  
  // Clean up listener after 10 seconds
  setTimeout(() => {
    document.removeEventListener('imageAnalysisComplete', analysisListener);
    console.log('⚠️ Analysis event not received within 10 seconds');
  }, 10000);
  
  return true;
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

// Test 6: Check if payment is blocking analysis
function testPaymentBlocking() {
  console.log('\n🚫 Test 6: Payment Blocking Check');
  
  // Check if payment popdown is visible
  const paymentPopdown = document.getElementById('payment-popdown');
  if (paymentPopdown) {
    const isVisible = paymentPopdown.style.display !== 'none' && 
                     !paymentPopdown.classList.contains('hidden');
    console.log('📊 Payment popdown visible:', isVisible);
    
    if (isVisible) {
      console.log('🚫 Analysis blocked by payment requirement');
      console.log('💡 Solution: Complete payment setup to enable analysis');
    }
  }
  
  // Check local storage for payment setup
  const deviceToken = localStorage.getItem('molDeviceToken');
  const cardInfo = localStorage.getItem('molCardInfo');
  
  if (!deviceToken || !cardInfo) {
    console.log('🚫 Payment not set up - analysis will be blocked');
    console.log('💡 Solution: Set up payment method to enable analysis');
  } else {
    console.log('✅ Payment appears to be set up');
  }
  
  return true;
}

// Test 7: Debug camera interaction handler
function testCameraInteraction() {
  console.log('\n🎯 Test 7: Camera Interaction Handler');
  
  if (window.cameraManager && window.cameraManager.handleInteraction) {
    console.log('✅ Camera interaction handler found');
    
    // Create a mock event
    const mockEvent = {
      clientX: 100,
      clientY: 100,
      preventDefault: () => console.log('Event prevented')
    };
    
    // Test the handler (this will show payment check)
    console.log('Testing camera interaction handler...');
    window.cameraManager.handleInteraction(mockEvent).then(() => {
      console.log('✅ Camera interaction handler completed');
    }).catch(err => {
      console.log('❌ Camera interaction handler error:', err);
    });
    
  } else {
    console.log('❌ Camera interaction handler not found');
  }
  
  return true;
}

// Run all tests
function runAllTests() {
  console.log('🚀 Running all camera click tests...\n');
  
  testCameraManager();
  testPaymentManager();
  testCameraClick();
  testPaymentBlocking();
  testCameraInteraction();
  
  // Wait a bit then run analysis tests
  setTimeout(() => {
    testImageAnalysis();
    testCompleteFlow();
  }, 1000);
  
  console.log('\n📋 Test Summary:');
  console.log('- Check console for individual test results');
  console.log('- Look for ❌ errors that indicate problems');
  console.log('- Payment issues are the most common cause of blocked analysis');
}

// Export functions for manual testing
window.testCameraClick = testCameraClick;
window.testPaymentBlocking = testPaymentBlocking;
window.testCompleteFlow = testCompleteFlow;
window.runAllTests = runAllTests;

console.log('✅ Debug functions loaded. Run runAllTests() to start testing.'); 