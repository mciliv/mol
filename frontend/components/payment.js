// payment.js - Payment and account management module

class PaymentManager {
  constructor() {
    this.stripe = null;
    this.cardElement = null;
    this.paymentRequest = null;
    this.setupInProgress = false;
  }

  // Check if we're in local development mode (DEPRECATED - use developer account instead)
  isLocalDevelopment() {
    console.warn('isLocalDevelopment() is deprecated - use developer account instead');
    return false;
  }

  // Check if current user is the developer account
  isDeveloperAccount() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      return false;
    }
    
    try {
      const card = JSON.parse(cardInfo);
      // Check for developer account indicators
      return card.brand === 'Development' || 
             card.last4 === '0000' ||
             deviceToken.includes('local_dev_token') ||
             deviceToken.includes('dev_');
    } catch (error) {
      return false;
    }
  }

  // Initialize payment setup
  async checkInitialPaymentSetup() {
    // Skip payment setup for developer account
    if (this.isDeveloperAccount()) {
      console.log('Developer account detected - skipping payment setup');
      return true;
    }

    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      this.showPaymentPopdown();
      return false;
    }
    
    try {
      const response = await fetch('/validate-payment', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ device_token: deviceToken })
      });
      
      if (!response.ok) {
        localStorage.removeItem('molDeviceToken');
        localStorage.removeItem('molCardInfo');
        this.showPaymentPopdown();
        return false;
      }
      
      const result = await response.json();
      this.updateAccountStatus(result.user);
      return true;
      
    } catch (error) {
      console.error('Initial payment check failed:', error);
      return true;
    }
  }

  // Show payment popdown
  showPaymentPopdown() {
    const popdown = document.getElementById('payment-popdown');
    const mainInterface = document.getElementById('main-app-interface');
    
    popdown.style.display = 'block';
    mainInterface.classList.add('payment-required');
    
    setTimeout(() => {
      popdown.classList.add('show');
    }, 10);
    
    this.initializePaymentSetup();
  }

  // Hide payment popdown
  hidePaymentPopdown() {
    const popdown = document.getElementById('payment-popdown');
    const mainInterface = document.getElementById('main-app-interface');
    
    popdown.classList.remove('show');
    mainInterface.classList.remove('payment-required');
    
    setTimeout(() => {
      popdown.style.display = 'none';
    }, 400);
  }

  // Update account status display
  updateAccountStatus(user) {
    const accountStatus = document.getElementById('account-status');
    const accountName = document.getElementById('account-name');
    
    if (user && user.name) {
      accountName.textContent = user.name;
    } else {
      accountName.textContent = 'Account';
    }
    
    accountStatus.style.display = 'flex';
  }

  // Check if payment is required
  isPaymentRequired() {
    // Skip payment requirement for developer account
    if (this.isDeveloperAccount()) {
      return false;
    }
    
    const paymentPopdown = document.getElementById('payment-popdown');
    return paymentPopdown && paymentPopdown.classList.contains('show');
  }

  // Check payment method validity
  async checkPaymentMethod() {
    // Skip payment check for developer account
    if (this.isDeveloperAccount()) {
      return true;
    }
    
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      this.showPaymentPopdown();
      return false;
    }
    
    try {
      const response = await fetch('/validate-payment', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ device_token: deviceToken })
      });
      
      if (!response.ok) {
        localStorage.removeItem('molDeviceToken');
        localStorage.removeItem('molCardInfo');
        this.showPaymentPopdown();
        return false;
      }
      
      const result = await response.json();
      const localCardInfo = JSON.parse(cardInfo);
      localCardInfo.usage = result.user.usage;
      localStorage.setItem('molCardInfo', JSON.stringify(localCardInfo));
      
      return true;
      
    } catch (error) {
      console.error('Payment validation error:', error);
      return true;
    }
  }

  // Increment usage counter
  async incrementUsage() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    if (!deviceToken) return;
    
    try {
      const response = await fetch('/increment-usage', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ device_token: deviceToken })
      });
      
      if (response.ok) {
        const result = await response.json();
        const cardInfo = JSON.parse(localStorage.getItem('molCardInfo') || '{}');
        cardInfo.usage = result.usage;
        localStorage.setItem('molCardInfo', JSON.stringify(cardInfo));
        
        console.log(`📊 Analysis complete - Total usage: ${result.usage}`);
      }
    } catch (error) {
      console.error('Usage increment error:', error);
      const cardInfo = JSON.parse(localStorage.getItem('molCardInfo') || '{}');
      cardInfo.usage = (cardInfo.usage || 0) + 1;
      localStorage.setItem('molCardInfo', JSON.stringify(cardInfo));
    }
  }

  // Reset local development user (DEPRECATED)
  resetLocalDevUser() {
    console.warn('resetLocalDevUser() is deprecated - use developer account instead');
    this.setupDeveloperAccount();
  }

  // Initialize payment setup UI
  initializePaymentSetup() {
    // Payment setup logic would go here
    console.log('Payment setup initialized');
  }

  // Setup local development user (DEPRECATED)
  setupLocalDevUser() {
    console.warn('setupLocalDevUser() is deprecated - use developer account instead');
    this.setupDeveloperAccount();
  }

  // Get local development user info (DEPRECATED)
  getLocalDevUser() {
    console.warn('getLocalDevUser() is deprecated - use developer account instead');
    return this.getDeveloperAccount();
  }

  // Setup developer account (your account)
  setupDeveloperAccount() {
    const developerUser = {
      name: 'Developer',
      email: 'developer@mol.com',
      usage: 0,
      device_token: 'dev_token_' + Date.now(),
      card_info: {
        last4: '0000',
        brand: 'Development',
        usage: 0
      }
    };

    // Store developer account data
    localStorage.setItem('molDeviceToken', developerUser.device_token);
    localStorage.setItem('molCardInfo', JSON.stringify(developerUser.card_info));
    localStorage.setItem('molDeveloperUser', JSON.stringify(developerUser));

    // Update UI to show developer account
    this.updateAccountStatus(developerUser);
    
    console.log('Developer account setup complete');
    return developerUser;
  }

  // Get developer account info
  getDeveloperAccount() {
    const developerUser = localStorage.getItem('molDeveloperUser');
    return developerUser ? JSON.parse(developerUser) : null;
  }
}

// Export singleton instance
export const paymentManager = new PaymentManager(); 