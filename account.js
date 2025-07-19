// Account Page JavaScript
let stripe;
let cardElement;
let isAuthenticated = false;
let currentUser = null;

// Initialize the page
document.addEventListener('DOMContentLoaded', function() {
    initializeTabs();
    initializeStripe();
    checkAuthentication();
    setupEventListeners();
});

// Tab functionality
function initializeTabs() {
    const tabBtns = document.querySelectorAll('.tab-btn');
    const tabContents = document.querySelectorAll('.tab-content');

    tabBtns.forEach(btn => {
        btn.addEventListener('click', () => {
            const targetTab = btn.getAttribute('data-tab');
            
            // Update active tab button
            tabBtns.forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            
            // Update active tab content
            tabContents.forEach(content => content.classList.remove('active'));
            document.getElementById(targetTab).classList.add('active');
        });
    });
}

// Initialize Stripe
async function initializeStripe() {
    try {
        // Get Stripe publishable key from server
        const response = await fetch('/stripe-config');
        const config = await response.json();
        
        stripe = Stripe(config.publishableKey);
        
        // Create card element
        const elements = stripe.elements();
        cardElement = elements.create('card', {
            style: {
                base: {
                    fontSize: '16px',
                    color: '#ffffff',
                    '::placeholder': {
                        color: '#888888'
                    }
                }
            }
        });
        
        cardElement.mount('#card-element');
        
        // Handle card errors
        cardElement.on('change', function(event) {
            const displayError = document.getElementById('card-errors');
            if (event.error) {
                displayError.textContent = event.error.message;
            } else {
                displayError.textContent = '';
            }
        });
        
    } catch (error) {
        console.error('Failed to initialize Stripe:', error);
        showPaymentError('Failed to initialize payment system');
    }
}

// Check if user is authenticated
function checkAuthentication() {
    const savedUser = localStorage.getItem('molUser');
    if (savedUser) {
        currentUser = JSON.parse(savedUser);
        isAuthenticated = true;
        updateAccountInfo();
        showAccountTab();
    }
}

// Setup event listeners
function setupEventListeners() {
    // Sign in form
    document.getElementById('signin-form').addEventListener('submit', handleSignIn);
    
    // Register form
    document.getElementById('register-form').addEventListener('submit', handleRegister);
    
    // Payment form
    document.getElementById('payment-form').addEventListener('submit', handlePayment);
    
    // Account actions
    document.getElementById('signout-btn').addEventListener('click', handleSignOut);
    document.getElementById('delete-account-btn').addEventListener('click', handleDeleteAccount);
    
    // Password manager autofill triggers
    setupPasswordManagerTriggers();
}

// Setup password manager autofill triggers
function setupPasswordManagerTriggers() {
    const signinEmail = document.getElementById('signin-email');
    const registerEmail = document.getElementById('register-email');
    
    // Trigger password autofill when email is entered
    signinEmail.addEventListener('input', function() {
        if (this.value.includes('@')) {
            // Small delay to let password manager detect the email
            setTimeout(() => {
                const passwordField = document.getElementById('signin-password');
                passwordField.focus();
                passwordField.blur();
            }, 100);
        }
    });
    
    registerEmail.addEventListener('input', function() {
        if (this.value.includes('@')) {
            // Small delay to let password manager detect the email
            setTimeout(() => {
                const passwordField = document.getElementById('register-password');
                passwordField.focus();
                passwordField.blur();
            }, 100);
        }
    });
}

// Handle sign in
async function handleSignIn(e) {
    e.preventDefault();
    
    const email = document.getElementById('signin-email').value;
    const password = document.getElementById('signin-password').value;
    
    // Demo mode - accept any credentials
    currentUser = {
        name: email.split('@')[0],
        email: email,
        memberSince: new Date().toLocaleDateString(),
        paymentMethod: 'Not set',
        analysesUsed: 0
    };
    
    isAuthenticated = true;
    localStorage.setItem('molUser', JSON.stringify(currentUser));
    
    updateAccountInfo();
    showAccountTab();
    
    // Show success message
    showMessage('Signed in successfully!', 'success');
}

// Handle register
async function handleRegister(e) {
    e.preventDefault();
    
    const name = document.getElementById('register-name').value;
    const email = document.getElementById('register-email').value;
    const password = document.getElementById('register-password').value;
    const confirm = document.getElementById('register-confirm').value;
    
    if (password !== confirm) {
        showMessage('Passwords do not match', 'error');
        return;
    }
    
    // Demo mode - create account
    currentUser = {
        name: name,
        email: email,
        memberSince: new Date().toLocaleDateString(),
        paymentMethod: 'Not set',
        analysesUsed: 0
    };
    
    isAuthenticated = true;
    localStorage.setItem('molUser', JSON.stringify(currentUser));
    
    updateAccountInfo();
    showAccountTab();
    
    // Show success message
    showMessage('Account created successfully!', 'success');
}

// Handle payment method setup
async function handlePayment(e) {
    e.preventDefault();
    
    if (!isAuthenticated) {
        showPaymentError('Please sign in first');
        return;
    }
    
    const saveBtn = document.getElementById('save-payment-btn');
    const btnText = saveBtn.querySelector('.btn-text');
    const btnLoading = saveBtn.querySelector('.btn-loading');
    
    // Show loading state
    saveBtn.disabled = true;
    btnText.style.display = 'none';
    btnLoading.style.display = 'inline';
    
    try {
        // Create payment method
        const { paymentMethod, error } = await stripe.createPaymentMethod({
            type: 'card',
            card: cardElement,
        });
        
        if (error) {
            throw new Error(error.message);
        }
        
        // Save payment method to server
        const response = await fetch('/save-payment-method', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                paymentMethodId: paymentMethod.id,
                email: currentUser.email
            })
        });
        
        const result = await response.json();
        
        if (!response.ok) {
            throw new Error(result.error || 'Failed to save payment method');
        }
        
        // Update user info
        currentUser.paymentMethod = `Card ending in ${paymentMethod.card.last4}`;
        localStorage.setItem('molUser', JSON.stringify(currentUser));
        updateAccountInfo();
        
        // Show success
        showPaymentSuccess();
        
        // Clear form
        cardElement.clear();
        
    } catch (error) {
        console.error('Payment error:', error);
        showPaymentError(error.message);
    } finally {
        // Reset button state
        saveBtn.disabled = false;
        btnText.style.display = 'inline';
        btnLoading.style.display = 'none';
    }
}

// Handle sign out
function handleSignOut() {
    currentUser = null;
    isAuthenticated = false;
    localStorage.removeItem('molUser');
    
    // Reset forms
    document.getElementById('signin-form').reset();
    document.getElementById('register-form').reset();
    if (cardElement) {
        cardElement.clear();
    }
    
    // Show sign in tab
    showSignInTab();
    
    showMessage('Signed out successfully', 'success');
}

// Handle delete account
function handleDeleteAccount() {
    if (confirm('Are you sure you want to delete your account? This action cannot be undone.')) {
        localStorage.removeItem('molUser');
        currentUser = null;
        isAuthenticated = false;
        
        // Reset forms
        document.getElementById('signin-form').reset();
        document.getElementById('register-form').reset();
        if (cardElement) {
            cardElement.clear();
        }
        
        // Show sign in tab
        showSignInTab();
        
        showMessage('Account deleted successfully', 'success');
    }
}

// Update account information display
function updateAccountInfo() {
    if (!currentUser) return;
    
    document.getElementById('account-name').textContent = currentUser.name;
    document.getElementById('account-email').textContent = currentUser.email;
    document.getElementById('account-date').textContent = currentUser.memberSince;
    document.getElementById('account-payment').textContent = currentUser.paymentMethod;
    document.getElementById('account-usage').textContent = currentUser.analysesUsed;
    
    document.getElementById('account-info').style.display = 'block';
}

// Show account tab
function showAccountTab() {
    const accountTab = document.querySelector('[data-tab="account"]');
    accountTab.click();
}

// Show sign in tab
function showSignInTab() {
    const signinTab = document.querySelector('[data-tab="signin"]');
    signinTab.click();
}

// Show payment success
function showPaymentSuccess() {
    const status = document.getElementById('payment-status');
    const success = document.getElementById('payment-success');
    const error = document.getElementById('payment-error');
    
    status.style.display = 'block';
    success.style.display = 'block';
    error.style.display = 'none';
    
    setTimeout(() => {
        status.style.display = 'none';
    }, 5000);
}

// Show payment error
function showPaymentError(message) {
    const status = document.getElementById('payment-status');
    const success = document.getElementById('payment-success');
    const error = document.getElementById('payment-error');
    
    status.style.display = 'block';
    success.style.display = 'none';
    error.style.display = 'block';
    error.textContent = message;
    
    setTimeout(() => {
        status.style.display = 'none';
    }, 5000);
}

// Show general message
function showMessage(message, type) {
    // Create temporary message element
    const msgElement = document.createElement('div');
    msgElement.className = `message ${type}`;
    msgElement.textContent = message;
    msgElement.style.cssText = `
        position: fixed;
        top: 20px;
        right: 20px;
        padding: 15px 20px;
        border-radius: 8px;
        color: white;
        font-weight: 600;
        z-index: 1000;
        animation: slideIn 0.3s ease;
    `;
    
    if (type === 'success') {
        msgElement.style.background = '#00ff00';
        msgElement.style.color = '#000';
    } else {
        msgElement.style.background = '#ff4444';
    }
    
    document.body.appendChild(msgElement);
    
    setTimeout(() => {
        msgElement.remove();
    }, 3000);
}

// Add CSS animation
const style = document.createElement('style');
style.textContent = `
    @keyframes slideIn {
        from {
            transform: translateX(100%);
            opacity: 0;
        }
        to {
            transform: translateX(0);
            opacity: 1;
        }
    }
`;
document.head.appendChild(style); 