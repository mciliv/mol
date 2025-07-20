// payment.js - Payment and account management module

class PaymentManager {
  constructor() {
    this.stripe = null;
    this.cardElement = null;
    this.paymentRequest = null;
    this.setupInProgress = false;
  }

  // Initialize payment setup
  async checkInitialPaymentSetup() {
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
    const paymentPopdown = document.getElementById('payment-popdown');
    return paymentPopdown && paymentPopdown.classList.contains('show');
  }

  // Check payment method validity
  async checkPaymentMethod() {
    // TEMPORARY: Bypass payment check for testing
    console.log('🔓 Payment check bypassed for testing');
    return true;
    
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

  // Initialize payment setup UI
  initializePaymentSetup() {
    // Payment setup logic would go here
    console.log('Payment setup initialized');
  }
}

// Export singleton instance
export const paymentManager = new PaymentManager(); 