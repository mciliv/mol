// Modern Card-First Account Setup
let stripe;
let cardElement;
let currentUser = null;

// Initialize the page
document.addEventListener('DOMContentLoaded', function() {
    initializeStripe();
    checkExistingUser();
    setupEventListeners();
});

// Initialize Stripe with modern styling
async function initializeStripe() {
    try {
        // Get Stripe publishable key from server
        const response = await fetch('/stripe-config');
        const config = await response.json();
        
        stripe = Stripe(config.publishableKey);
        
        // Create card element with modern dark theme
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
            hidePostalCode: true // No postal codes needed!
        });
        
        cardElement.mount('#card-element');
        
        // Modern card error handling
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
        
    } catch (error) {
        console.error('Failed to initialize Stripe:', error);
        showError('Payment system unavailable. Please try again.');
    }
}

// Check if user is already recognized on this device
function checkExistingUser() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (deviceToken && cardInfo) {
        const card = JSON.parse(cardInfo);
        currentUser = {
            deviceToken,
            cardLast4: card.last4,
            usage: card.usage || 0,
            name: card.name
        };
        
        showExistingUser();
    } else {
        showNewUserSetup();
    }
}

// Show welcome back screen for recognized users
function showExistingUser() {
    document.getElementById('existing-user').style.display = 'block';
    document.getElementById('new-user').style.display = 'none';
    
    document.getElementById('card-last4').textContent = currentUser.cardLast4;
    document.getElementById('user-usage').textContent = currentUser.usage;
}

// Show setup screen for new users
function showNewUserSetup() {
    document.getElementById('existing-user').style.display = 'none';
    document.getElementById('new-user').style.display = 'block';
}

// Setup event listeners
function setupEventListeners() {
    // Card setup form
    document.getElementById('card-setup-form').addEventListener('submit', handleCardSetup);
    
    // Continue with existing user
    document.getElementById('continue-btn').addEventListener('click', () => {
        showSuccessAndRedirect();
    });
    
    // Use different card
    document.getElementById('use-different-card').addEventListener('click', () => {
        // Clear stored data and show new user setup
        localStorage.removeItem('molDeviceToken');
        localStorage.removeItem('molCardInfo');
        currentUser = null;
        showNewUserSetup();
    });
}

// Handle the streamlined card setup
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
        
        // Setup payment method with server (includes 3D Secure if needed)
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
        
        // Handle 3D Secure authentication if required
        if (result.requires_action) {
            const { error: confirmError } = await stripe.confirmCardSetup(
                result.client_secret
            );
            
            if (confirmError) {
                throw new Error(confirmError.message);
            }
        }
        
        // Store user info locally for future recognition
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
        
        currentUser = {
            deviceToken,
            cardLast4: cardInfo.last4,
            usage: 0,
            name: cardInfo.name
        };
        
        // Show success
        showSuccess();
        
    } catch (error) {
        console.error('Setup error:', error);
        showError(error.message);
        
        // Reset button state
        btnText.style.display = 'inline';
        btnLoading.style.display = 'none';
        setupBtn.disabled = false;
    }
}

// Generate a device-specific token for recognition
function generateDeviceToken() {
    const fingerprint = getDeviceFingerprint();
    const timestamp = Date.now();
    const random = Math.random().toString(36).substring(2);
    
    return btoa(`${fingerprint}-${timestamp}-${random}`).replace(/[+/=]/g, '');
}

// Get basic device fingerprint for recognition (privacy-conscious)
function getDeviceFingerprint() {
    const screen = `${window.screen.width}x${window.screen.height}`;
    const timezone = Intl.DateTimeFormat().resolvedOptions().timeZone;
    const language = navigator.language;
    const platform = navigator.platform;
    
    return btoa(`${screen}-${timezone}-${language}-${platform}`);
}

// Show success state
function showSuccess() {
    document.getElementById('new-user').style.display = 'none';
    document.getElementById('existing-user').style.display = 'none';
    document.getElementById('setup-success').style.display = 'block';
}

// Show success and redirect after delay
function showSuccessAndRedirect() {
    showSuccess();
    setTimeout(() => {
        window.location.href = 'index.html';
    }, 1500);
}

// Show error message
function showError(message) {
    const errorDiv = document.getElementById('card-errors');
    errorDiv.textContent = message;
    errorDiv.style.display = 'block';
    
    // Auto-hide after 5 seconds
    setTimeout(() => {
        errorDiv.style.display = 'none';
    }, 5000);
}

// Export user status for main app
window.getMolUser = function() {
    return currentUser;
};

// Check if user has valid payment setup
window.checkPaymentMethod = function() {
    return !!currentUser && !!currentUser.deviceToken;
};

// Increment usage counter
window.incrementUsage = function() {
    if (currentUser) {
        currentUser.usage++;
        const cardInfo = JSON.parse(localStorage.getItem('molCardInfo'));
        cardInfo.usage = currentUser.usage;
        localStorage.setItem('molCardInfo', JSON.stringify(cardInfo));
    }
}; 