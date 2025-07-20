// Enhanced debugging utilities for Cursor/VS Code debugging
// Place breakpoints at the lines marked with "// 🔴 BREAKPOINT:" comments

// Global debugging state tracker
window.debugState = {
  lastEvent: null,
  currentInput: null,
  appState: null,
  apiCalls: [],
  errors: []
};

// Enhanced state inspector for debugging
window.inspectAppState = () => {
  const app = window.molecularApp;
  const state = {
    // Input state
    objectInputValue: document.getElementById('object-input')?.value || '',
    objectInputFocused: document.activeElement?.id === 'object-input',
    
    // UI state
    paymentPopdown: document.getElementById('payment-popdown'),
    paymentDisplay: getComputedStyle(document.getElementById('payment-popdown') || {}).display,
    mainApp: document.getElementById('main-app-interface'),
    mainAppDisplay: getComputedStyle(document.getElementById('main-app-interface') || {}).display,
    
    // App instance state
    molecularApp: app ? {
      isProcessing: app.isProcessing,
      hasPaymentSetup: app.hasPaymentSetup,
      lastAnalysis: app.lastAnalysis,
      currentAnalysisType: app.currentAnalysisType
    } : null,
    
    // API state
    recentApiCalls: window.debugState.apiCalls.slice(-3),
    recentErrors: window.debugState.errors.slice(-3),
    
    // DOM state
    allInputs: Array.from(document.querySelectorAll('input')).map(input => ({
      id: input.id,
      type: input.type,
      value: input.value,
      focused: input === document.activeElement
    })),
    
    timestamp: Date.now()
  };
  
  window.debugState.appState = state;
  return state;
};

// Enhanced API call tracker
const originalFetch = window.fetch;
window.fetch = function(...args) {
  const callInfo = {
    url: args[0],
    options: args[1] || {},
    timestamp: Date.now(),
    stack: new Error().stack
  };
  
  // 🔴 BREAKPOINT: Set breakpoint here to inspect all API calls
  console.log('🌐 API Call:', callInfo);
  window.debugState.apiCalls.push(callInfo);
  
  return originalFetch.apply(this, args)
    .then(response => {
      callInfo.response = {
        status: response.status,
        statusText: response.statusText,
        headers: Object.fromEntries(response.headers.entries())
      };
      // 🔴 BREAKPOINT: Set breakpoint here to inspect API responses
      console.log('✅ API Response:', callInfo);
      return response;
    })
    .catch(error => {
      callInfo.error = error.message;
      // 🔴 BREAKPOINT: Set breakpoint here to inspect API errors
      console.log('❌ API Error:', callInfo);
      window.debugState.errors.push(callInfo);
      throw error;
    });
};

// Enhanced error tracking
window.addEventListener('error', (event) => {
  const errorInfo = {
    message: event.message,
    filename: event.filename,
    lineno: event.lineno,
    colno: event.colno,
    error: event.error,
    timestamp: Date.now()
  };
  
  // 🔴 BREAKPOINT: Set breakpoint here to inspect JavaScript errors
  console.error('💥 JavaScript Error:', errorInfo);
  window.debugState.errors.push(errorInfo);
});

document.addEventListener('DOMContentLoaded', () => {
  console.log('🔧 Debug system loaded - breakpoint-ready');
  console.log('📍 Set breakpoints on lines marked with "🔴 BREAKPOINT:" comments');
  console.log('🎯 Key debugging functions available:');
  console.log('   • window.inspectAppState() - Get detailed app state');
  console.log('   • window.debugState - Access debugging history');
  
  const objectInput = document.getElementById('object-input');
  
  if (objectInput) {
    // Enhanced keyup handler with debugging hooks
    objectInput.addEventListener('keyup', (e) => {
      window.debugState.lastEvent = {
        type: 'keyup',
        key: e.key,
        target: e.target.id,
        value: e.target.value,
        timestamp: Date.now()
      };
      
      console.log(`⬆️ KEYUP on object-input:`, { 
        key: e.key, 
        value: e.target.value, 
        willTriggerHandler: e.key === 'Enter' 
      });
      
      if (e.key === 'Enter') {
        // 🔴 BREAKPOINT: Set breakpoint here to debug Enter key handling
        console.log('🚀 Enter pressed - about to trigger text analysis');
        const currentState = window.inspectAppState();
        console.log('📊 Current app state:', currentState);
        
        // Additional debugging info for Enter key press
        console.log('🔍 Enter Key Debug Details:', {
          inputElement: e.target,
          inputValue: e.target.value,
          inputFocused: document.activeElement === e.target,
          eventPhase: e.eventPhase,
          bubbles: e.bubbles,
          cancelable: e.cancelable,
          defaultPrevented: e.defaultPrevented,
          molecularAppExists: !!window.molecularApp,
          timestamp: Date.now()
        });
      }
    });
    
    // Enhanced input change tracking
    objectInput.addEventListener('input', (e) => {
      window.debugState.currentInput = {
        value: e.target.value,
        timestamp: Date.now(),
        length: e.target.value.length
      };
      
      // 🔴 BREAKPOINT: Set breakpoint here to debug input changes
      console.log('✏️ INPUT CHANGE:', window.debugState.currentInput);
    });
    
    // Enhanced focus/blur tracking
    objectInput.addEventListener('focus', (e) => {
      // 🔴 BREAKPOINT: Set breakpoint here to debug focus events
      console.log('🎯 FOCUS on object-input:', {
        value: e.target.value,
        timestamp: Date.now(),
        previousActiveElement: window.debugState.lastActiveElement
      });
      window.debugState.lastActiveElement = e.target;
    });
    
    objectInput.addEventListener('blur', (e) => {
      // 🔴 BREAKPOINT: Set breakpoint here to debug blur events
      const blurInfo = {
        id: e.target.id,
        type: e.target.type,
        value: e.target.value,
        name: e.target.name,
        timestamp: Date.now()
      };
      
      console.log('📝 BLUR on INPUT:', blurInfo);
      window.debugState.lastBlur = blurInfo;
    });
  }
  
  // Track all form submissions
  document.addEventListener('submit', (e) => {
    // 🔴 BREAKPOINT: Set breakpoint here to debug form submissions
    console.log('📤 FORM SUBMIT:', {
      form: e.target,
      formData: new FormData(e.target),
      timestamp: Date.now()
    });
  });
  
  // Enhanced click tracking for buttons
  document.addEventListener('click', (e) => {
    if (e.target.tagName === 'BUTTON' || e.target.type === 'button') {
      // 🔴 BREAKPOINT: Set breakpoint here to debug button clicks
      console.log('🖱️ BUTTON CLICK:', {
        button: e.target,
        id: e.target.id,
        className: e.target.className,
        textContent: e.target.textContent,
        timestamp: Date.now()
      });
    }
  });
});

// Global debugging commands
window.debugCommands = {
  // Get current state
  state: () => window.inspectAppState(),
  
  // Start detailed logging
  startDetailedLogging: () => {
    window.debugState.detailedLogging = true;
    console.log('🔍 Detailed logging started');
  },
  
  // Stop detailed logging
  stopDetailedLogging: () => {
    window.debugState.detailedLogging = false;
    console.log('🔍 Detailed logging stopped');
  },
  
  // Clear debug history
  clearHistory: () => {
    window.debugState.apiCalls = [];
    window.debugState.errors = [];
    console.log('🧹 Debug history cleared');
  },
  
  // Export debug session
  exportSession: () => {
    const session = {
      debugState: window.debugState,
      currentState: window.inspectAppState(),
      timestamp: Date.now(),
      userAgent: navigator.userAgent,
      url: window.location.href
    };
    
    console.log('📊 Debug Session Export:', session);
    return session;
  }
};

// Make debugging easier in console
window.debug = window.debugCommands;

console.log('🚀 Enhanced debugging system loaded!');
console.log('💡 Use window.debug.state() for current state');
console.log('📍 Set breakpoints on lines with "🔴 BREAKPOINT:" comments'); 