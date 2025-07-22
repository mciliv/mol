class PaymentManager {
  constructor() {
    this.stripe = null;
    this.cardElement = null;
    this.paymentRequest = null;
    this.setupInProgress = false;
  }

  isLocalDevelopment() {
    console.warn('isLocalDevelopment() is deprecated - use developer account instead');
    return false;
  }

  isDeveloperAccount() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      return false;
    }
    
    try {
      const card = JSON.parse(cardInfo);
      return card.brand === 'Development' || 
             card.last4 === '0000' ||
             deviceToken.includes('local_dev_token') ||
             deviceToken.includes('dev_');
    } catch (error) {
      return false;
    }
  }

  async checkInitialPaymentSetup() {
    // Check if payment is enabled via toggle
    const paymentEnabled = localStorage.getItem('molPaymentEnabled') === 'true';
    
    if (!paymentEnabled) {
      console.log('Payment disabled via toggle - no payment required');
      return true;
    }

    // Payment is enabled - check for actual payment method
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      console.log('Payment enabled but no method found - showing payment popdown');
      this.showPaymentPopdown();
      return false;
    }
    
    console.log('Payment method found - setup complete');
    return true;
  }

  showPaymentPopdown() {
    console.log('🔄 Showing payment popdown...');
    const popdown = document.getElementById('payment-popdown');
    const mainInterface = document.getElementById('main-app-interface');
    
    console.log('📊 Elements found:', { popdown: !!popdown, mainInterface: !!mainInterface });
    
    if (popdown) {
      popdown.style.display = 'block';
      popdown.classList.remove('hidden');
      
      setTimeout(() => {
        popdown.classList.add('show');
        console.log('✅ Payment popdown should now be visible');
      }, 10);
    } else {
      console.error('❌ Payment popdown element not found!');
    }
    
    if (mainInterface) {
      mainInterface.classList.add('payment-required');
    }
    
    this.initializePaymentSetup();
  }

  hidePaymentPopdown() {
    console.log('🔄 Hiding payment popdown...');
    const popdown = document.getElementById('payment-popdown');
    const mainInterface = document.getElementById('main-app-interface');
    
    if (popdown) {
      popdown.classList.remove('show');
      setTimeout(() => {
        popdown.style.display = 'none';
        popdown.classList.add('hidden');
      }, 300);
    }
    
    if (mainInterface) {
      mainInterface.classList.remove('payment-required');
    }
    
    console.log('✅ Payment popdown hidden');
  }



  updateAccountStatus(user) {
    const accountStatus = document.getElementById('account-status');
    const accountName = document.getElementById('account-name');
    
    if (accountStatus) {
      // Always show the account status
      accountStatus.classList.remove('hidden');
      accountStatus.style.display = 'flex';
      
      if (user && user.name) {
        accountName.textContent = user.name;
        accountStatus.style.color = '#00d4ff'; // Blue when set up
      } else {
        accountName.textContent = 'Add Card';
        accountStatus.style.color = '#ffa500'; // Orange when not set up
      }
      
      // Add click handler to show payment popdown when no payment set up
      const paymentEnabled = localStorage.getItem('molPaymentEnabled') === 'true';
      if (paymentEnabled) {
        const deviceToken = localStorage.getItem('molDeviceToken');
        const cardInfo = localStorage.getItem('molCardInfo');
        
        if (!deviceToken || !cardInfo) {
          // No payment set up - make it clickable to show popdown
          accountStatus.style.cursor = 'pointer';
          accountStatus.onclick = () => {
            this.showPaymentPopdown();
          };
        } else {
          // Payment set up - remove click handler
          accountStatus.style.cursor = 'default';
          accountStatus.onclick = null;
        }
      }
    }
  }

  isPaymentRequired() {
    // Check if payment is enabled via toggle
    const paymentEnabled = localStorage.getItem('molPaymentEnabled') === 'true';
    
    if (!paymentEnabled) {
      return false;
    }
    
    const paymentPopdown = document.getElementById('payment-popdown');
    return paymentPopdown && paymentPopdown.classList.contains('show');
  }

  async checkPaymentMethod() {
    // Check if payment is enabled via toggle
    const paymentEnabled = localStorage.getItem('molPaymentEnabled') === 'true';
    
    if (!paymentEnabled) {
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

  resetLocalDevUser() {
    console.warn('resetLocalDevUser() is deprecated - use developer account instead');
    this.setupDeveloperAccount();
  }

  initializePaymentSetup() {
    console.log('Payment setup initialized');
  }

  setupLocalDevUser() {
    console.warn('setupLocalDevUser() is deprecated - use developer account instead');
    this.setupDeveloperAccount();
  }

  getLocalDevUser() {
    console.warn('getLocalDevUser() is deprecated - use developer account instead');
    return this.getDeveloperAccount();
  }

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

    localStorage.setItem('molDeviceToken', developerUser.device_token);
    localStorage.setItem('molCardInfo', JSON.stringify(developerUser.card_info));
    localStorage.setItem('molDeveloperUser', JSON.stringify(developerUser));

    this.updateAccountStatus(developerUser);
    
    console.log('Developer account setup complete');
    return developerUser;
  }

  getDeveloperAccount() {
    const developerUser = localStorage.getItem('molDeveloperUser');
    return developerUser ? JSON.parse(developerUser) : null;
  }
  
  // Function to clear payment setup for testing
  clearPaymentSetup() {
    localStorage.removeItem('molDeviceToken');
    localStorage.removeItem('molCardInfo');
    localStorage.removeItem('molDeveloperUser');
    console.log('🧹 Payment setup cleared for testing');
  }
}

// Export singleton instance
export const paymentManager = new PaymentManager(); 