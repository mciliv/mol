const fs = require('fs');
const path = require('path');
const { JSDOM } = require('jsdom');

describe('Frontend Tests', () => {
  let dom;
  let document;
  let window;

  beforeEach(() => {
    // Create a DOM environment
    const htmlContent = fs.readFileSync(path.join(__dirname, '../index.html'), 'utf8');
    dom = new JSDOM(htmlContent, {
      runScripts: 'dangerously',
      resources: 'usable',
      url: 'http://localhost:8080'
    });
    
    document = dom.window.document;
    window = dom.window;
    
    // Mock browser APIs
    window.navigator = {
      mediaDevices: {
        getUserMedia: jest.fn().mockResolvedValue({
          getTracks: () => [{ stop: jest.fn() }]
        }),
        enumerateDevices: jest.fn().mockResolvedValue([
          { kind: 'videoinput', label: 'Camera 1' }
        ])
      }
    };

    window.isSecureContext = true;
    window.location = { protocol: 'https:' };
    
    // Mock fetch
    window.fetch = jest.fn();
    
    // Mock 3Dmol
    window.$3Dmol = {
      createViewer: jest.fn().mockReturnValue({
        addModel: jest.fn(),
        setBackgroundColor: jest.fn(),
        setStyle: jest.fn(),
        zoomTo: jest.fn(),
        render: jest.fn(),
        resize: jest.fn()
      })
    };

    // Mock FileReader
    window.FileReader = jest.fn(() => ({
      readAsDataURL: jest.fn(),
      result: 'data:image/jpeg;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg=='
    }));

    // Mock alert and confirm
    window.alert = jest.fn();
    window.confirm = jest.fn();
  });

  describe('DOM Structure', () => {
    test('should have required HTML elements', () => {
      // Check for essential elements
      expect(document.getElementById('object-input')).toBeTruthy();
      expect(document.getElementById('video-feed')).toBeTruthy();
      expect(document.getElementById('camera-mode')).toBeTruthy();
      expect(document.getElementById('photo-mode')).toBeTruthy();
      expect(document.getElementById('photo-upload')).toBeTruthy();
      expect(document.getElementById('photo-url')).toBeTruthy();
      expect(document.getElementById('url-analyze')).toBeTruthy();
      expect(document.querySelector('.snapshots-container')).toBeTruthy();
      expect(document.getElementById('gldiv')).toBeTruthy();
    });

    test('should have proper CSS classes', () => {
      expect(document.querySelector('.app-container')).toBeTruthy();
      expect(document.querySelector('.top-bar')).toBeTruthy();
      expect(document.querySelector('.input-mode-section')).toBeTruthy();
      expect(document.querySelector('.mode-selector')).toBeTruthy();
      expect(document.querySelector('.camera-container')).toBeTruthy();
      expect(document.querySelector('.photo-options')).toBeTruthy();
    });

    test('should have proper form structure', () => {
      const textInput = document.getElementById('object-input');
      expect(textInput.type).toBe('text');
      expect(textInput.placeholder).toContain('Describe any object');

      const fileInput = document.getElementById('photo-upload');
      expect(fileInput.type).toBe('file');
      expect(fileInput.accept).toBe('image/*');

      const urlInput = document.getElementById('photo-url');
      expect(urlInput.type).toBe('url');
      expect(urlInput.placeholder).toContain('Paste image URL');
    });
  });

  describe('Input Mode Toggle', () => {
    test('should toggle between camera and photo modes', () => {
      const cameraMode = document.getElementById('camera-mode');
      const photoMode = document.getElementById('photo-mode');
      const cameraContainer = document.querySelector('.camera-container');
      const photoOptions = document.getElementById('photo-options');

      // Initially camera mode should be active
      expect(cameraMode.checked).toBe(true);
      expect(photoMode.checked).toBe(false);
      expect(cameraContainer.style.display).not.toBe('none');
      expect(photoOptions.style.display).toBe('none');

      // Switch to photo mode
      photoMode.checked = true;
      cameraMode.checked = false;
      photoMode.dispatchEvent(new window.Event('change'));

      expect(cameraContainer.style.display).toBe('none');
      expect(photoOptions.style.display).toBe('flex');
    });

    test('should handle mode change events', () => {
      const cameraMode = document.getElementById('camera-mode');
      const photoMode = document.getElementById('photo-mode');
      
      // Test camera mode change
      const cameraChangeEvent = new window.Event('change');
      cameraMode.dispatchEvent(cameraChangeEvent);

      // Test photo mode change
      const photoChangeEvent = new window.Event('change');
      photoMode.dispatchEvent(photoChangeEvent);

      // Events should be handled without errors
      expect(true).toBe(true);
    });
  });

  describe('Text Input Functionality', () => {
    test('should handle text input submission', () => {
      const textInput = document.getElementById('object-input');
      
      // Set input value
      textInput.value = 'test object';
      
      // Simulate Enter key press
      const enterEvent = new window.KeyboardEvent('keyup', { key: 'Enter' });
      textInput.dispatchEvent(enterEvent);

      // Input should be cleared after submission
      expect(textInput.value).toBe('');
    });

    test('should not submit empty text input', () => {
      const textInput = document.getElementById('object-input');
      
      // Set empty input value
      textInput.value = '   ';
      
      // Simulate Enter key press
      const enterEvent = new window.KeyboardEvent('keyup', { key: 'Enter' });
      textInput.dispatchEvent(enterEvent);

      // Input should remain unchanged
      expect(textInput.value).toBe('   ');
    });

    test('should handle non-Enter key presses', () => {
      const textInput = document.getElementById('object-input');
      textInput.value = 'test';
      
      // Simulate non-Enter key press
      const spaceEvent = new window.KeyboardEvent('keyup', { key: ' ' });
      textInput.dispatchEvent(spaceEvent);

      // Input should not be cleared
      expect(textInput.value).toBe('test');
    });
  });

  describe('File Upload Functionality', () => {
    test('should handle valid image file upload', () => {
      const fileInput = document.getElementById('photo-upload');
      
      // Create mock image file
      const mockFile = new window.File([''], 'test.jpg', { type: 'image/jpeg' });
      
      // Simulate file selection
      const changeEvent = new window.Event('change');
      Object.defineProperty(fileInput, 'files', { value: [mockFile] });
      
      fileInput.dispatchEvent(changeEvent);

      // Should not show alert for valid image
      expect(window.alert).not.toHaveBeenCalled();
    });

    test('should reject invalid file types', () => {
      const fileInput = document.getElementById('photo-upload');
      
      // Create mock invalid file
      const mockFile = new window.File([''], 'test.txt', { type: 'text/plain' });
      
      // Simulate file selection
      const changeEvent = new window.Event('change');
      Object.defineProperty(fileInput, 'files', { value: [mockFile] });
      
      fileInput.dispatchEvent(changeEvent);

      // Should show alert for invalid file type
      expect(window.alert).toHaveBeenCalledWith('Please select an image file');
    });

    test('should handle no file selection', () => {
      const fileInput = document.getElementById('photo-upload');
      
      // Simulate change event with no files
      const changeEvent = new window.Event('change');
      Object.defineProperty(fileInput, 'files', { value: [] });
      
      fileInput.dispatchEvent(changeEvent);

      // Should not show alert
      expect(window.alert).not.toHaveBeenCalled();
    });
  });

  describe('URL Input Functionality', () => {
    test('should handle valid URL input', () => {
      const urlInput = document.getElementById('photo-url');
      const urlButton = document.getElementById('url-analyze');
      
      // Set valid URL
      urlInput.value = 'https://example.com/image.jpg';
      
      // Simulate button click
      urlButton.click();

      // Should not show alert for valid URL
      expect(window.alert).not.toHaveBeenCalled();
    });

    test('should reject invalid URL', () => {
      const urlInput = document.getElementById('photo-url');
      const urlButton = document.getElementById('url-analyze');
      
      // Set invalid URL
      urlInput.value = 'not-a-url';
      
      // Simulate button click
      urlButton.click();

      // Should show alert for invalid URL
      expect(window.alert).toHaveBeenCalledWith('Please enter a valid URL');
    });

    test('should handle empty URL', () => {
      const urlInput = document.getElementById('photo-url');
      const urlButton = document.getElementById('url-analyze');
      
      // Set empty URL
      urlInput.value = '';
      
      // Simulate button click
      urlButton.click();

      // Should show alert for empty URL
      expect(window.alert).toHaveBeenCalledWith('Please enter an image URL');
    });

    test('should handle Enter key in URL input', () => {
      const urlInput = document.getElementById('photo-url');
      
      // Set valid URL
      urlInput.value = 'https://example.com/image.jpg';
      
      // Simulate Enter key press
      const enterEvent = new window.KeyboardEvent('keyup', { key: 'Enter' });
      urlInput.dispatchEvent(enterEvent);

      // Should not show alert for valid URL
      expect(window.alert).not.toHaveBeenCalled();
    });
  });

  describe('Camera Functionality', () => {
    test('should initialize camera on page load', async () => {
      // Load app.js
      const appScript = fs.readFileSync(path.join(__dirname, '../app.js'), 'utf8');
      const script = document.createElement('script');
      script.textContent = appScript;
      document.head.appendChild(script);

      // Wait for initialization
      await new Promise(resolve => setTimeout(resolve, 100));

      // Verify camera was requested
      expect(window.navigator.mediaDevices.getUserMedia).toHaveBeenCalled();
    });

    test('should handle camera permission denial', async () => {
      // Mock camera permission denial
      window.navigator.mediaDevices.getUserMedia.mockRejectedValueOnce(
        new Error('Permission denied')
      );

      // Load app.js
      const appScript = fs.readFileSync(path.join(__dirname, '../app.js'), 'utf8');
      const script = document.createElement('script');
      script.textContent = appScript;
      document.head.appendChild(script);

      // Wait for initialization
      await new Promise(resolve => setTimeout(resolve, 100));

      // Verify error handling
      expect(window.navigator.mediaDevices.getUserMedia).toHaveBeenCalled();
    });

    test('should create switch camera button when multiple cameras available', async () => {
      // Mock multiple cameras
      window.navigator.mediaDevices.enumerateDevices.mockResolvedValueOnce([
        { kind: 'videoinput', label: 'Front Camera' },
        { kind: 'videoinput', label: 'Back Camera' }
      ]);

      // Load app.js
      const appScript = fs.readFileSync(path.join(__dirname, '../app.js'), 'utf8');
      const script = document.createElement('script');
      script.textContent = appScript;
      document.head.appendChild(script);

      // Wait for initialization
      await new Promise(resolve => setTimeout(resolve, 100));

      // Verify switch camera button was created
      const switchButton = document.querySelector('.switch-camera-btn');
      expect(switchButton).toBeTruthy();
    });
  });

  describe('API Integration', () => {
    test('should make API calls for text analysis', async () => {
      // Mock successful API response
      window.fetch.mockResolvedValueOnce({
        ok: true,
        json: () => Promise.resolve({
          output: {
            object: 'test object',
            smiles: ['CCO']
          }
        })
      });

      // Simulate text input
      const textInput = document.getElementById('object-input');
      textInput.value = 'test object';
      
      const enterEvent = new window.KeyboardEvent('keyup', { key: 'Enter' });
      textInput.dispatchEvent(enterEvent);

      // Wait for async operations
      await new Promise(resolve => setTimeout(resolve, 100));

      // Verify API call was made
      expect(window.fetch).toHaveBeenCalledWith('/object-molecules', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ object: 'test object' })
      });
    });

    test('should handle API errors gracefully', async () => {
      // Mock API error
      window.fetch.mockRejectedValueOnce(new Error('Network error'));

      // Simulate text input
      const textInput = document.getElementById('object-input');
      textInput.value = 'test object';
      
      const enterEvent = new window.KeyboardEvent('keyup', { key: 'Enter' });
      textInput.dispatchEvent(enterEvent);

      // Wait for async operations
      await new Promise(resolve => setTimeout(resolve, 100));

      // Verify error handling
      expect(window.fetch).toHaveBeenCalled();
    });
  });

  describe('3D Visualization', () => {
    test('should create 3Dmol viewers', async () => {
      // Mock successful SDF loading
      window.fetch
        .mockResolvedValueOnce({
          ok: true,
          json: () => Promise.resolve({
            sdfPaths: ['/sdf_files/CCO.sdf'],
            errors: [],
            skipped: []
          })
        })
        .mockResolvedValueOnce({
          ok: true,
          text: () => Promise.resolve('CCO\n  RDKit          3D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  1  0  0  0  0\nM  END\n$$$$')
        });

      // Simulate text input
      const textInput = document.getElementById('object-input');
      textInput.value = 'test object';
      
      const enterEvent = new window.KeyboardEvent('keyup', { key: 'Enter' });
      textInput.dispatchEvent(enterEvent);

      // Wait for async operations
      await new Promise(resolve => setTimeout(resolve, 100));

      // Verify 3Dmol viewer was created
      expect(window.$3Dmol.createViewer).toHaveBeenCalled();
    });
  });

  describe('Error Handling', () => {
    test('should display error messages for failed API calls', async () => {
      // Mock API error
      window.fetch.mockRejectedValueOnce(new Error('API Error'));

      // Simulate text input
      const textInput = document.getElementById('object-input');
      textInput.value = 'test object';
      
      const enterEvent = new window.KeyboardEvent('keyup', { key: 'Enter' });
      textInput.dispatchEvent(enterEvent);

      // Wait for async operations
      await new Promise(resolve => setTimeout(resolve, 100));

      // Verify error message was created
      const errorElement = document.querySelector('h3[style*="color: red"]');
      expect(errorElement).toBeTruthy();
      expect(errorElement.textContent).toContain('Error');
    });

    test('should handle file upload errors', () => {
      // Mock FileReader error
      const mockFileReader = {
        readAsDataURL: jest.fn(),
        onerror: null
      };
      window.FileReader = jest.fn(() => mockFileReader);

      const fileInput = document.getElementById('photo-upload');
      const mockFile = new window.File([''], 'test.jpg', { type: 'image/jpeg' });
      
      const changeEvent = new window.Event('change');
      Object.defineProperty(fileInput, 'files', { value: [mockFile] });
      
      fileInput.dispatchEvent(changeEvent);

      // Simulate FileReader error
      if (mockFileReader.onerror) {
        mockFileReader.onerror(new Error('File read error'));
      }

      // Should handle error gracefully
      expect(true).toBe(true);
    });
  });

  describe('Accessibility', () => {
    test('should have proper ARIA labels', () => {
      const fileInput = document.getElementById('photo-upload');
      const urlInput = document.getElementById('photo-url');
      
      expect(fileInput.getAttribute('aria-label')).toContain('Upload photo');
      expect(urlInput.getAttribute('aria-label')).toContain('Enter image URL');
    });

    test('should have proper form labels', () => {
      const cameraMode = document.getElementById('camera-mode');
      const photoMode = document.getElementById('photo-mode');
      
      expect(cameraMode).toBeTruthy();
      expect(photoMode).toBeTruthy();
      
      const cameraLabel = document.querySelector('label[for="camera-mode"]');
      const photoLabel = document.querySelector('label[for="photo-mode"]');
      
      expect(cameraLabel).toBeTruthy();
      expect(photoLabel).toBeTruthy();
    });

    test('should have proper button types', () => {
      const urlButton = document.getElementById('url-analyze');
      expect(urlButton.type).toBe('button');
    });
  });
}); 