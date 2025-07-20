// debug-events.js - Event debugging utilities

console.log('🐛 Debug script loading...');

class EventDebugger {
  constructor() {
    this.originalAddEventListener = EventTarget.prototype.addEventListener;
    this.eventLog = [];
    this.isLogging = false;
    console.log('🔧 EventDebugger initialized');
  }

  // Start logging all events
  startLogging() {
    this.isLogging = true;
    this.eventLog = [];
    
    // Override addEventListener to intercept all event registrations
    EventTarget.prototype.addEventListener = (type, listener, options) => {
      if (this.isLogging) {
        console.log(`🔗 Event listener added: ${type} on`, this);
      }
      return this.originalAddEventListener.call(this, type, listener, options);
    };
    
    console.log('🐛 Event debugging started. All events will be logged.');
  }

  // Stop logging
  stopLogging() {
    this.isLogging = false;
    EventTarget.prototype.addEventListener = this.originalAddEventListener;
    console.log('🐛 Event debugging stopped.');
  }

  // Log a specific event
  logEvent(event, element, handler) {
    if (!this.isLogging) return;
    
    const logEntry = {
      timestamp: Date.now(),
      type: event.type,
      target: element.tagName || element.constructor.name,
      id: element.id || 'no-id',
      className: element.className || 'no-class',
      key: event.key || null,
      code: event.code || null,
      value: element.value || null,
      handler: handler.name || 'anonymous',
      stack: new Error().stack
    };
    
    this.eventLog.push(logEntry);
    console.log('📝 Event captured:', logEntry);
  }

  // Get recent events
  getRecentEvents(count = 10) {
    return this.eventLog.slice(-count);
  }

  // Clear event log
  clearLog() {
    this.eventLog = [];
    console.log('🗑️ Event log cleared');
  }
}

// Global event debugger instance
window.eventDebugger = new EventDebugger();

// Debug specific enter key events
function debugEnterKeyEvents() {
  console.log('🎯 Setting up Enter key debugging...');
  
  // Wait for elements to be available
  const checkForElements = () => {
    // Debug object input
    const objectInput = document.getElementById('object-input');
    if (objectInput) {
      console.log('✅ Found object-input element');
      
      objectInput.addEventListener('keydown', (e) => {
        console.log(`⬇️ KEYDOWN on object-input:`, {
          key: e.key,
          code: e.code,
          target: e.target.tagName,
          value: e.target.value,
          timestamp: Date.now()
        });
      });
      
      objectInput.addEventListener('keyup', (e) => {
        console.log(`⬆️ KEYUP on object-input:`, {
          key: e.key,
          code: e.code,
          target: e.target.tagName,
          value: e.target.value,
          willTriggerHandler: e.key === 'Enter',
          timestamp: Date.now()
        });
        
        if (e.key === 'Enter') {
          console.log('🚀 Enter pressed - about to trigger text analysis');
          console.log('📊 Current app state:', {
            objectInputValue: e.target.value,
            paymentPopdown: document.getElementById('payment-popdown'),
            paymentDisplay: document.getElementById('payment-popdown')?.style.display,
            mainApp: document.getElementById('main-app-interface'),
            mainAppDisplay: document.getElementById('main-app-interface')?.style.display
          });
        }
      });
      
      objectInput.addEventListener('input', (e) => {
        console.log('📝 INPUT on object-input:', e.target.value);
      });
      
    } else {
      console.log('❌ object-input element not found');
    }

    // Debug photo URL input  
    const photoUrl = document.getElementById('photo-url');
    if (photoUrl) {
      console.log('✅ Found photo-url element');
      
      photoUrl.addEventListener('keydown', (e) => {
        console.log(`⬇️ KEYDOWN on photo-url:`, {
          key: e.key,
          code: e.code,
          value: e.target.value,
          timestamp: Date.now()
        });
      });
      
      photoUrl.addEventListener('keyup', (e) => {
        console.log(`⬆️ KEYUP on photo-url:`, {
          key: e.key,
          code: e.code,
          value: e.target.value,
          willTriggerHandler: e.key === 'Enter',
          timestamp: Date.now()
        });
        
        if (e.key === 'Enter') {
          console.log('🚀 Enter pressed on URL field - about to trigger URL analysis');
        }
      });
    } else {
      console.log('❌ photo-url element not found');
    }
  };

  // Try immediately and also wait for DOM
  checkForElements();
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', checkForElements);
  }
  
  console.log('✅ Enter key debugging setup complete');
}

// Debug form submissions
function debugFormEvents() {
  console.log('📋 Setting up form debugging...');
  
  const checkForms = () => {
    const forms = document.querySelectorAll('form');
    console.log(`Found ${forms.length} forms`);
    
    forms.forEach((form, index) => {
      form.addEventListener('submit', (e) => {
        console.log(`📝 FORM SUBMIT #${index}:`, {
          formId: form.id,
          action: form.action,
          method: form.method,
          preventDefault: e.defaultPrevented,
          timestamp: Date.now()
        });
      });
    });
  };
  
  checkForms();
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', checkForms);
  }
}

// Debug all input events
function debugAllInputs() {
  console.log('🎛️ Setting up input debugging...');
  
  const checkInputs = () => {
    const inputs = document.querySelectorAll('input, textarea');
    console.log(`Found ${inputs.length} input elements`);
    
    inputs.forEach((input, index) => {
      console.log(`Input ${index}: ${input.tagName} id="${input.id}" type="${input.type}"`);
      
      ['focus', 'blur', 'input', 'change'].forEach(eventType => {
        input.addEventListener(eventType, (e) => {
          console.log(`📝 ${eventType.toUpperCase()} on ${input.tagName}:`, {
            id: input.id,
            type: input.type,
            value: input.value,
            name: input.name,
            timestamp: Date.now()
          });
        });
      });
    });
  };
  
  checkInputs();
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', checkInputs);
  }
}

// Simple immediate debugging function
function simpleDebug() {
  console.log('🔍 Simple debug check:');
  console.log('- Document ready state:', document.readyState);
  console.log('- object-input exists:', !!document.getElementById('object-input'));
  console.log('- photo-url exists:', !!document.getElementById('photo-url'));
  console.log('- All inputs:', document.querySelectorAll('input').length);
  console.log('- Active element:', document.activeElement);
  
  // Add immediate keydown listener to document
  document.addEventListener('keydown', (e) => {
    console.log('🔑 Global keydown:', e.key, 'on', e.target.tagName, e.target.id || 'no-id');
  });
}

// Comprehensive debugging setup
function setupDebugMode() {
  console.log('🔍 Setting up comprehensive event debugging...');
  
  // Run simple debug immediately
  simpleDebug();
  
  debugEnterKeyEvents();
  debugFormEvents();
  debugAllInputs();
  
  // Add global keyboard debugger
  document.addEventListener('keydown', (e) => {
    // Log all keydowns first
    console.log('🌍 Global keydown captured:', {
      key: e.key,
      code: e.code,
      ctrl: e.ctrlKey,
      shift: e.shiftKey,
      alt: e.altKey,
      target: e.target.tagName,
      targetId: e.target.id
    });
    
    if (e.ctrlKey && e.shiftKey && e.key === 'D') {
      console.log('🎛️ Debug toggle pressed (Ctrl+Shift+D)');
      if (window.eventDebugger.isLogging) {
        window.eventDebugger.stopLogging();
      } else {
        window.eventDebugger.startLogging();
      }
    }
  });
  
  console.log('✅ Debug mode active. Press Ctrl+Shift+D to toggle detailed logging.');
  console.log('💡 Available commands:');
  console.log('  - debugEnterKeyEvents() - Track enter key presses');
  console.log('  - debugAllInputs() - Track all input events');
  console.log('  - simpleDebug() - Quick debug check');
  console.log('  - eventDebugger.startLogging() - Start comprehensive logging');
  console.log('  - eventDebugger.getRecentEvents() - Get recent events');
}

// Auto-setup with multiple triggers
console.log('🚀 Auto-setting up debug mode...');

// Try immediate setup
setupDebugMode();

// Also setup on DOMContentLoaded
if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', () => {
    console.log('📄 DOM loaded - setting up debug again');
    setupDebugMode();
  });
}

// And also on window load
window.addEventListener('load', () => {
  console.log('🪟 Window loaded - final debug setup');
  setupDebugMode();
});

// Export for global access
window.debugEnterKeyEvents = debugEnterKeyEvents;
window.debugFormEvents = debugFormEvents;
window.debugAllInputs = debugAllInputs;
window.setupDebugMode = setupDebugMode;
window.simpleDebug = simpleDebug;

console.log('🐛 Debug script loaded completely'); 