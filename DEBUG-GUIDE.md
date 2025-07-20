# 🐛 Cursor Debugging Guide for Molecular Analysis App

## 🎯 **Quick Start - Set These Breakpoints First**

### **Primary Debugging Points** (🔴 BREAKPOINT markers in code)

1. **Enter Key Handler** - `app.js:91`
   ```js
   // 🔴 BREAKPOINT: Set breakpoint here to debug text analysis trigger
   console.log('🚀 Text analysis triggered from Enter key');
   ```

2. **Main Analysis Flow** - `app.js:104`
   ```js
   // 🔴 BREAKPOINT: Set breakpoint here to debug main analysis flow
   console.log('🔬 Starting handleTextAnalysis');
   ```

3. **Payment Check** - `app.js:119`
   ```js
   // 🔴 BREAKPOINT: Set breakpoint here to debug payment check
   if (!this.hasPaymentSetup) {
   ```

4. **API Call** - `app.js:133`
   ```js
   // 🔴 BREAKPOINT: Set breakpoint here to debug API call preparation
   console.log('🌐 Preparing API call for text analysis');
   ```

5. **API Response** - `app.js:140`
   ```js
   // 🔴 BREAKPOINT: Set breakpoint here to debug API response
   console.log('📡 API response received:', {
   ```

---

## 🚀 **How to Debug in Cursor**

### **Step 1: Start the App with Debugger**
```bash
./run dev
```
The app automatically starts with Node.js debugging enabled.

### **Step 2: Open Cursor Debugger**
1. **Press `Ctrl+Shift+D`** (or `Cmd+Shift+D` on Mac)
2. Click **"Run and Debug"** panel
3. VS Code auto-attaches to the running Node.js process

### **Step 3: Set Breakpoints**
1. **Open `app.js`** in Cursor
2. **Find lines marked with** `🔴 BREAKPOINT:` comments
3. **Click in the gutter** (left of line numbers) to set red dot breakpoints
4. **Set breakpoints on these key lines:**
   - Line 91: Enter key trigger
   - Line 104: Analysis start
   - Line 119: Payment check
   - Line 133: API preparation
   - Line 140: API response

### **Step 4: Set Browser Breakpoints**
1. **Open browser** → `http://localhost:8080`
2. **Open DevTools** (`F12`)
3. **Go to Sources tab**
4. **Find `debug-events.js`**
5. **Set breakpoints on lines marked** `🔴 BREAKPOINT:`

---

## 🔍 **Debugging Workflow**

### **When You Type and Press Enter:**

1. **Browser Event** (debug-events.js)
   - Keyup event fires
   - Console logs event details
   - Triggers app analysis

2. **Server-Side Flow** (app.js) 
   - Enter key handler activates
   - Payment check occurs
   - API call prepared
   - Response processed

### **What You'll See at Each Breakpoint:**

**🎯 Breakpoint 1 (Enter Key):**
```js
Variables panel shows:
- e.key = "Enter"
- this.objectInput.value = "your input"
- this.isProcessing = false/true
- this.hasPaymentSetup = false/true
```

**🎯 Breakpoint 2 (Analysis Start):**
```js
Variables panel shows:
- inputValue = trimmed input
- this.isProcessing = current state
- this.currentAnalysisType = 'text'
```

**🎯 Breakpoint 3 (Payment Check):**
```js
Variables panel shows:
- this.hasPaymentSetup = boolean
- this.paymentPopdown = DOM element
```

**🎯 Breakpoint 4 (API Call):**
```js
Variables panel shows:
- inputValue = what's being sent
- fetch request details
- request body JSON
```

**🎯 Breakpoint 5 (API Response):**
```js
Variables panel shows:
- response.status = 200/400/500
- response.ok = true/false
- result = parsed JSON data
```

---

## 🛠️ **Enhanced Debugging Commands**

### **In Browser Console:**
```js
// Get current app state
window.debug.state()

// Export full debug session
window.debug.exportSession()

// Check if app is processing
window.inspectAppState()

// Clear debug history
window.debug.clearHistory()
```

### **In Cursor Debug Console:**
```js
// Inspect the main app instance
molecularApp

// Check current processing state
molecularApp.isProcessing

// View payment setup
molecularApp.hasPaymentSetup

// Check last analysis
molecularApp.lastAnalysis
```

---

## 📊 **What Information to Collect**

### **When Sharing Debug Results with AI:**

**🔴 Essential Info:**
1. **Breakpoint location** (file:line)
2. **Variable values** from Variables panel
3. **Console output** at that moment
4. **Call stack** from Call Stack panel
5. **Any error messages**

**📋 Template for Reporting:**
```
BREAKPOINT: app.js:91 (Enter key handler)

VARIABLES:
- e.key: "Enter"
- this.objectInput.value: "caffeine"
- this.isProcessing: false
- this.hasPaymentSetup: true

CONSOLE OUTPUT:
🚀 Text analysis triggered from Enter key
📊 App state before analysis: {...}

CALL STACK:
setupTextAnalysis (app.js:91)
EventTarget.addEventListener (<anonymous>)

ISSUE: Payment check failing even though hasPaymentSetup is true
```

---

## 🚨 **Common Debugging Scenarios**

### **❌ Enter Key Not Working:**
**Breakpoints to check:**
1. `debug-events.js:64` - Is keyup firing?
2. `app.js:91` - Is Enter key handler reached?
3. Check if `this.objectInput` exists

### **💳 Payment Issues:**
**Breakpoints to check:**
1. `app.js:119` - Payment check logic
2. Check `this.hasPaymentSetup` value
3. Inspect payment popdown display state

### **🌐 API Call Problems:**
**Breakpoints to check:**
1. `app.js:133` - Request preparation
2. `app.js:140` - Response handling
3. Check network tab for actual HTTP requests

### **🔄 Processing State Issues:**
**Breakpoints to check:**
1. `app.js:125` - Processing start
2. `app.js:161` - Processing cleanup
3. Check `this.isProcessing` value changes

---

## 💡 **Pro Tips**

1. **Use Conditional Breakpoints:** Right-click breakpoint → Add condition
   ```js
   this.objectInput.value === "specific-test-case"
   ```

2. **Watch Expressions:** Add variables to Watch panel
   ```js
   this.isProcessing
   this.hasPaymentSetup
   this.objectInput.value
   ```

3. **Call Stack Navigation:** Click on any function in Call Stack to jump there

4. **Console Evaluation:** Use Debug Console to run code in current context
   ```js
   this.objectInput.value = "test"
   this.handleTextAnalysis()
   ```

5. **Step Through Code:**
   - `F10` - Step Over
   - `F11` - Step Into
   - `Shift+F11` - Step Out
   - `F5` - Continue

---

## 🎯 **Next Steps After Setting Breakpoints**

1. **Test the flow:** Type something and press Enter
2. **Step through each breakpoint** with F10/F11
3. **Inspect variables** at each step
4. **Note any unexpected values**
5. **Export debug session** if you need to share details

This setup gives you complete visibility into what happens when you press Enter! 🕵️‍♂️ 