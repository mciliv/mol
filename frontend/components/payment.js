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
    // Check for actual payment method regardless of toggle
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      console.log('No payment method found - showing payment modal for new user');
      this.showPaymentModal();
      return false;
    }
    
    console.log('Payment method found - setup complete');
    return true;
  }

  showPaymentModal() {
    console.log('🔄 Showing payment modal...');
    const modal = document.getElementById('payment-modal');
    const mainInterface = document.getElementById('main-app-interface');
    
    console.log('📊 Elements found:', { modal: !!modal, mainInterface: !!mainInterface });
    
    if (modal) {
      modal.style.display = 'block';
      modal.classList.remove('hidden');
      
      setTimeout(() => {
        modal.classList.add('show');
        console.log('✅ Payment modal should now be visible');
      }, 10);
    } else {
      console.error('❌ Payment modal element not found!');
    }
    
    if (mainInterface) {
      mainInterface.classList.add('payment-required');
      mainInterface.classList.add('modal-showing');
    }
    
    this.initializePaymentSetup();
  }

  hidePaymentModal() {
    console.log('🔄 Hiding payment modal...');
    const modal = document.getElementById('payment-modal');
    const mainInterface = document.getElementById('main-app-interface');
    
    if (modal) {
      modal.classList.remove('show');
      setTimeout(() => {
        modal.style.display = 'none';
        modal.classList.add('hidden');
      }, 300);
    }
    
    if (mainInterface) {
      mainInterface.classList.remove('payment-required');
      mainInterface.classList.remove('modal-showing');
    }
    
    console.log('✅ Payment modal hidden');
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
      
      // Add click handler to show payment modal when no payment set up
      const deviceToken = localStorage.getItem('molDeviceToken');
      const cardInfo = localStorage.getItem('molCardInfo');
      
      if (!deviceToken || !cardInfo) {
        // No payment set up - make it clickable to show modal
        accountStatus.style.cursor = 'pointer';
        accountStatus.onclick = () => {
          this.showPaymentModal();
        };
      } else {
        // Payment set up - make it clickable to show card management
        accountStatus.style.cursor = 'pointer';
        accountStatus.onclick = () => {
          this.showCardManagementModal();
        };
      }
    }
  }

  // Card Management Modal Methods
  showCardManagementModal() {
    console.log('🔄 Showing card management modal...');
    const modal = document.getElementById('card-management-modal');
    const mainInterface = document.getElementById('main-app-interface');
    
    if (modal) {
      modal.style.display = 'block';
      modal.classList.remove('hidden');
      
      setTimeout(() => {
        modal.classList.add('show');
        console.log('✅ Card management modal should now be visible');
      }, 10);
    } else {
      console.error('❌ Card management modal element not found!');
    }
    
    if (mainInterface) {
      mainInterface.classList.add('modal-showing');
    }
    
    this.loadUserCards();
  }

  hideCardManagementModal() {
    console.log('🔄 Hiding card management modal...');
    const modal = document.getElementById('card-management-modal');
    const mainInterface = document.getElementById('main-app-interface');
    
    if (modal) {
      modal.classList.remove('show');
      setTimeout(() => {
        modal.style.display = 'none';
        modal.classList.add('hidden');
      }, 300);
    }
    
    if (mainInterface) {
      mainInterface.classList.remove('modal-showing');
    }
    
    console.log('✅ Card management modal hidden');
  }

  async loadUserCards() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    if (!deviceToken) {
      console.error('No device token found');
      return;
    }

    try {
      const response = await fetch(`/get-payment-methods?device_token=${encodeURIComponent(deviceToken)}`, {
        method: 'GET',
        headers: { 'Content-Type': 'application/json' }
      });

      if (!response.ok) {
        throw new Error(`Failed to load cards: ${response.status}`);
      }

      const result = await response.json();
      this.displayCards(result.payment_methods, result.default_method);
      
    } catch (error) {
      console.error('Failed to load user cards:', error);
      this.showCardError('Failed to load payment methods');
    }
  }

  displayCards(cards, defaultMethodId) {
    const cardsList = document.getElementById('cards-list');
    if (!cardsList) return;

    cardsList.innerHTML = '';

    if (!cards || cards.length === 0) {
      cardsList.innerHTML = `
        <div class="empty-state">
          <p>No payment methods</p>
          <button type="button" class="text-button primary" onclick="paymentManager.showAddCardForm()">Add your first payment method</button>
        </div>
      `;
      return;
    }

    cards.forEach(card => {
      const methodElement = document.createElement('div');
      methodElement.className = `payment-method ${card.is_default ? 'default' : ''}`;
      methodElement.innerHTML = `
        <div class="method-info">
          <div class="method-primary">
            <span class="method-brand">${this.formatCardBrand(card.brand)}</span>
            <span class="method-number">•••• ${card.last4}</span>
            ${card.is_default ? '<span class="default-label">Default</span>' : ''}
          </div>
          <div class="method-secondary">
            Expires ${card.exp_month}/${card.exp_year}
          </div>
        </div>
        <div class="method-actions">
          <button type="button" class="action-link" onclick="paymentManager.editCard('${card.id}')">Edit</button>
          ${!card.is_default ? `<button type="button" class="action-link" onclick="paymentManager.setDefaultCard('${card.id}')">Set as default</button>` : ''}
          <button type="button" class="action-link delete" onclick="paymentManager.deleteCard('${card.id}')">Delete</button>
        </div>
      `;
      cardsList.appendChild(methodElement);
    });
  }

  formatCardBrand(brand) {
    const brands = {
      'visa': 'Visa',
      'mastercard': 'Mastercard',
      'amex': 'American Express',
      'discover': 'Discover',
      'diners': 'Diners Club',
      'jcb': 'JCB',
      'unionpay': 'UnionPay'
    };
    return brands[brand] || brand.charAt(0).toUpperCase() + brand.slice(1);
  }

  showAddCardForm() {
    const currentCards = document.getElementById('current-cards');
    const editForm = document.getElementById('edit-card-form');
    const editTitle = document.getElementById('edit-card-title');
    
    if (currentCards) currentCards.classList.add('hidden');
    if (editForm) editForm.classList.remove('hidden');
    if (editTitle) editTitle.textContent = 'Add New Payment Method';
    
    this.currentEditingCardId = null;
    this.setupCardForm('edit-card-element');
  }

  editCard(cardId) {
    const currentCards = document.getElementById('current-cards');
    const editForm = document.getElementById('edit-card-form');
    const editTitle = document.getElementById('edit-card-title');
    
    if (currentCards) currentCards.classList.add('hidden');
    if (editForm) editForm.classList.remove('hidden');
    if (editTitle) editTitle.textContent = 'Edit Payment Method';
    
    this.currentEditingCardId = cardId;
    this.setupCardForm('edit-card-element');
  }

  async setDefaultCard(cardId) {
    const deviceToken = localStorage.getItem('molDeviceToken');
    if (!deviceToken) return;

    try {
      const response = await fetch('/set-default-payment-method', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ 
          device_token: deviceToken,
          payment_method_id: cardId 
        })
      });

      if (!response.ok) {
        throw new Error(`Failed to set default card: ${response.status}`);
      }

      console.log('✅ Default card updated');
      this.loadUserCards(); // Refresh the cards list
      
    } catch (error) {
      console.error('Failed to set default card:', error);
      this.showCardError('Failed to update default card');
    }
  }

  async deleteCard(cardId) {
    if (!confirm('Are you sure you want to delete this payment method?')) {
      return;
    }

    const deviceToken = localStorage.getItem('molDeviceToken');
    if (!deviceToken) return;

    try {
      const response = await fetch('/delete-payment-method', {
        method: 'DELETE',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ 
          device_token: deviceToken,
          payment_method_id: cardId 
        })
      });

      if (!response.ok) {
        throw new Error(`Failed to delete card: ${response.status}`);
      }

      console.log('✅ Card deleted');
      this.loadUserCards(); // Refresh the cards list
      
      // If this was the only card, clear local storage and show payment setup
      const result = await response.json();
      if (result.message && result.message.includes('last card')) {
        localStorage.removeItem('molDeviceToken');
        localStorage.removeItem('molCardInfo');
        this.hideCardManagementModal();
        this.showPaymentModal();
      }
      
    } catch (error) {
      console.error('Failed to delete card:', error);
      this.showCardError('Failed to delete payment method');
    }
  }

  cancelCardEdit() {
    const currentCards = document.getElementById('current-cards');
    const editForm = document.getElementById('edit-card-form');
    
    if (currentCards) currentCards.classList.remove('hidden');
    if (editForm) editForm.classList.add('hidden');
    
    this.currentEditingCardId = null;
  }

  setupCardForm(elementId) {
    // This would integrate with Stripe Elements
    // For demo purposes, we'll create a simple placeholder
    const cardElement = document.getElementById(elementId);
    if (cardElement) {
      cardElement.innerHTML = `
        <div style="padding: 12px; background: rgba(255,255,255,0.1); border-radius: 6px; color: #fff;">
          <input type="text" placeholder="Card number" style="background: transparent; border: none; color: #fff; width: 100%; margin-bottom: 8px;">
          <div style="display: flex; gap: 8px;">
            <input type="text" placeholder="MM/YY" style="background: transparent; border: none; color: #fff; flex: 1;">
            <input type="text" placeholder="CVC" style="background: transparent; border: none; color: #fff; flex: 1;">
          </div>
        </div>
      `;
    }
  }

  async saveCardChanges() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    const userName = document.getElementById('edit-user-name').value;
    
    if (!deviceToken) return;

    try {
      const saveBtn = document.getElementById('save-card-btn');
      const btnText = saveBtn.querySelector('.btn-text');
      const btnLoading = saveBtn.querySelector('.btn-loading');
      
      btnText.classList.add('hidden');
      btnLoading.classList.remove('hidden');
      saveBtn.disabled = true;

      const endpoint = this.currentEditingCardId ? '/update-payment-method' : '/setup-payment-method';
      const payload = {
        device_token: deviceToken,
        payment_method: 'pm_demo_' + Date.now(), // Demo payment method
        name: userName
      };

      const response = await fetch(endpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload)
      });

      if (!response.ok) {
        throw new Error(`Failed to save card: ${response.status}`);
      }

      console.log('✅ Card saved successfully');
      this.cancelCardEdit();
      this.loadUserCards();
      
    } catch (error) {
      console.error('Failed to save card:', error);
      this.showCardError('Failed to save payment method');
    } finally {
      const saveBtn = document.getElementById('save-card-btn');
      if (saveBtn) {
        const btnText = saveBtn.querySelector('.btn-text');
        const btnLoading = saveBtn.querySelector('.btn-loading');
        
        btnText.classList.remove('hidden');
        btnLoading.classList.add('hidden');
        saveBtn.disabled = false;
      }
    }
  }

  showCardError(message) {
    // Simple error display - could be enhanced with a toast notification
    alert(message);
  }

  isPaymentRequired() {
    const paymentModal = document.getElementById('payment-modal');
    return paymentModal && paymentModal.classList.contains('show');
  }

  async checkPaymentMethod() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      this.showPaymentModal();
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
        this.showPaymentModal();
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