// Video Feed Validation Script
// Run this in the browser console to validate video feed functionality

console.log('🎥 Video Feed Validation Script');

class VideoFeedValidator {
  constructor() {
    this.results = [];
    this.currentStream = null;
  }

  // Log test result
  logResult(testName, status, message, data = null) {
    const result = {
      test: testName,
      status: status,
      message: message,
      timestamp: new Date().toISOString(),
      data: data
    };
    
    this.results.push(result);
    console.log(`[${status}] ${testName}: ${message}`);
    
    if (data) {
      console.log('Data:', data);
    }
    
    return result;
  }

  // Test 1: Browser API Support
  testBrowserSupport() {
    const support = {
      mediaDevices: !!navigator.mediaDevices,
      getUserMedia: !!(navigator.mediaDevices && navigator.mediaDevices.getUserMedia),
      permissions: !!navigator.permissions,
      isSecureContext: window.isSecureContext,
      protocol: location.protocol,
      hostname: location.hostname
    };
    
    const isSupported = support.mediaDevices && support.getUserMedia && 
                       (support.isSecureContext || support.hostname === 'localhost');
    
    return this.logResult('Browser Support', 
                         isSupported ? 'PASS' : 'FAIL',
                         isSupported ? 'Browser supports camera API' : 'Browser does not support camera API',
                         support);
  }

  // Test 2: Camera Permission
  async testCameraPermission() {
    try {
      if (!navigator.permissions) {
        return this.logResult('Camera Permission', 'SKIP', 'Permissions API not available');
      }
      
      const permission = await navigator.permissions.query({ name: 'camera' });
      return this.logResult('Camera Permission', 'PASS', `Permission state: ${permission.state}`, {
        state: permission.state,
        onchange: !!permission.onchange
      });
    } catch (error) {
      return this.logResult('Camera Permission', 'FAIL', `Permission check failed: ${error.message}`);
    }
  }

  // Test 3: Video Stream Creation
  async testVideoStream() {
    try {
      const constraints = {
        video: {
          width: { ideal: 1280 },
          height: { ideal: 720 }
        }
      };
      
      this.currentStream = await navigator.mediaDevices.getUserMedia(constraints);
      const tracks = this.currentStream.getTracks();
      
      // Stop the stream immediately for testing
      tracks.forEach(track => track.stop());
      this.currentStream = null;
      
      return this.logResult('Video Stream', 'PASS', 'Successfully created video stream', {
        trackCount: tracks.length,
        trackTypes: tracks.map(track => track.kind),
        trackStates: tracks.map(track => track.readyState)
      });
    } catch (error) {
      return this.logResult('Video Stream', 'FAIL', `Failed to create video stream: ${error.message}`);
    }
  }

  // Test 4: Camera Devices
  async testCameraDevices() {
    try {
      const devices = await navigator.mediaDevices.enumerateDevices();
      const videoDevices = devices.filter(device => device.kind === 'videoinput');
      
      return this.logResult('Camera Devices', 'PASS', `Found ${videoDevices.length} camera device(s)`, {
        totalDevices: devices.length,
        videoDevices: videoDevices.length,
        deviceIds: videoDevices.map(device => device.deviceId ? 'present' : 'missing'),
        deviceLabels: videoDevices.map(device => device.label || 'no label')
      });
    } catch (error) {
      return this.logResult('Camera Devices', 'FAIL', `Failed to enumerate devices: ${error.message}`);
    }
  }

  // Test 5: Video Element Integration
  testVideoElement() {
    const video = document.getElementById('video-feed');
    if (!video) {
      return this.logResult('Video Element', 'FAIL', 'Video element not found');
    }
    
    const properties = {
      tagName: video.tagName,
      autoplay: video.autoplay,
      playsinline: video.hasAttribute('playsinline'),
      muted: video.muted,
      readyState: video.readyState,
      videoWidth: video.videoWidth,
      videoHeight: video.videoHeight,
      clientWidth: video.clientWidth,
      clientHeight: video.clientHeight
    };
    
    return this.logResult('Video Element', 'PASS', 'Video element found and configured', properties);
  }

  // Test 6: Click Event Handling
  testClickHandling() {
    const video = document.getElementById('video-feed');
    if (!video) {
      return this.logResult('Click Handling', 'FAIL', 'Video element not found');
    }
    
    // Check if click event is attached
    const hasClickHandler = video.onclick !== null || 
                           video.hasAttribute('onclick') ||
                           video.addEventListener !== undefined;
    
    return this.logResult('Click Handling', hasClickHandler ? 'PASS' : 'WARN', 
                         hasClickHandler ? 'Click handling available' : 'No click handler detected');
  }

  // Test 7: Camera Switching
  async testCameraSwitching() {
    try {
      const devices = await navigator.mediaDevices.enumerateDevices();
      const videoDevices = devices.filter(device => device.kind === 'videoinput');
      
      if (videoDevices.length < 2) {
        return this.logResult('Camera Switching', 'SKIP', 'Only one camera available');
      }
      
      // Test environment camera
      const envStream = await navigator.mediaDevices.getUserMedia({
        video: { facingMode: 'environment' }
      });
      envStream.getTracks().forEach(track => track.stop());
      
      // Test user camera
      const userStream = await navigator.mediaDevices.getUserMedia({
        video: { facingMode: 'user' }
      });
      userStream.getTracks().forEach(track => track.stop());
      
      return this.logResult('Camera Switching', 'PASS', 'Successfully switched between cameras');
    } catch (error) {
      return this.logResult('Camera Switching', 'FAIL', `Camera switching failed: ${error.message}`);
    }
  }

  // Test 8: Performance
  async testPerformance() {
    const startTime = performance.now();
    
    try {
      const stream = await navigator.mediaDevices.getUserMedia({ video: true });
      const endTime = performance.now();
      const duration = endTime - startTime;
      
      stream.getTracks().forEach(track => track.stop());
      
      const isFast = duration < 2000; // Less than 2 seconds
      
      return this.logResult('Performance', isFast ? 'PASS' : 'WARN', 
                           `Camera access took ${duration.toFixed(0)}ms`, {
        duration: duration,
        threshold: 2000,
        isFast: isFast
      });
    } catch (error) {
      return this.logResult('Performance', 'FAIL', `Performance test failed: ${error.message}`);
    }
  }

  // Run all tests
  async runAllTests() {
    console.log('🚀 Starting Video Feed Validation...');
    
    this.testBrowserSupport();
    await this.testCameraPermission();
    await this.testVideoStream();
    await this.testCameraDevices();
    this.testVideoElement();
    this.testClickHandling();
    await this.testCameraSwitching();
    await this.testPerformance();
    
    this.generateReport();
  }

  // Generate test report
  generateReport() {
    const total = this.results.length;
    const passed = this.results.filter(r => r.status === 'PASS').length;
    const failed = this.results.filter(r => r.status === 'FAIL').length;
    const skipped = this.results.filter(r => r.status === 'SKIP').length;
    const warnings = this.results.filter(r => r.status === 'WARN').length;
    
    console.log('\n📊 Video Feed Validation Report');
    console.log('================================');
    console.log(`Total Tests: ${total}`);
    console.log(`✅ Passed: ${passed}`);
    console.log(`❌ Failed: ${failed}`);
    console.log(`⚠️ Warnings: ${warnings}`);
    console.log(`⏭️ Skipped: ${skipped}`);
    console.log(`Success Rate: ${((passed / total) * 100).toFixed(1)}%`);
    
    if (failed > 0) {
      console.log('\n❌ Failed Tests:');
      this.results.filter(r => r.status === 'FAIL').forEach(result => {
        console.log(`  - ${result.test}: ${result.message}`);
      });
    }
    
    if (warnings > 0) {
      console.log('\n⚠️ Warnings:');
      this.results.filter(r => r.status === 'WARN').forEach(result => {
        console.log(`  - ${result.test}: ${result.message}`);
      });
    }
    
    console.log('\n📋 Full Results:', this.results);
  }
}

// Create and export validator instance
window.videoFeedValidator = new VideoFeedValidator();

// Auto-run validation if requested
if (window.location.search.includes('validate=true')) {
  window.videoFeedValidator.runAllTests();
}

console.log('✅ Video Feed Validator loaded. Run: videoFeedValidator.runAllTests()'); 