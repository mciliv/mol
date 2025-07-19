# Password Manager Compatibility Rules

## Universal Password Manager Best Practices

### 1. Form Structure
- Use semantic HTML5 form elements
- Include proper `action` and `method` attributes
- Add `autocomplete="on"` to forms
- Use descriptive `id` and `name` attributes

### 2. Email/Username Fields
```html
<!-- Standard login email field -->
<input type="email" 
       id="email" 
       name="email" 
       autocomplete="username email" 
       required>
```

### 3. Password Fields
```html
<!-- Login password -->
<input type="password" 
       id="password" 
       name="password" 
       autocomplete="current-password" 
       required>

<!-- Registration password -->
<input type="password" 
       id="password" 
       name="password" 
       autocomplete="new-password" 
       required>
```

### 4. Platform-Specific Attributes

#### Apple Keychain (Safari, iOS)
- Prefers `autocomplete="username"` for email fields
- Recognizes `name="login"` hidden fields
- Works best with `action` URLs matching saved domains
- Supports `autocomplete="email"` as secondary attribute

#### Chrome Password Manager
- Uses `autocomplete="username"` for login fields
- Recognizes `type="email"` and `type="password"`
- Works with `name` attributes matching common patterns
- Supports `autocomplete="email"` as fallback

#### Firefox Password Manager
- Prefers `autocomplete="username"` for login
- Recognizes `type="email"` and `type="password"`
- Works with semantic form structure
- Supports multiple autocomplete attributes

#### Edge Password Manager
- Uses `autocomplete="username"` for login fields
- Recognizes `type="email"` and `type="password"`
- Works with `name` attributes
- Supports `autocomplete="email"` as secondary

#### 1Password, LastPass, Bitwarden
- Use standard `autocomplete` attributes
- Recognize `type="email"` and `type="password"`
- Work with semantic form structure
- Support multiple autocomplete values

### 5. Universal Implementation
```html
<form action="https://domain.com/login" method="post" autocomplete="on">
  <input type="hidden" name="login" value="1">
  <input type="hidden" name="remember" value="1">
  
  <label for="email">Email</label>
  <input type="email" 
         id="email" 
         name="email" 
         autocomplete="username email" 
         required>
  
  <label for="password">Password</label>
  <input type="password" 
         id="password" 
         name="password" 
         autocomplete="current-password" 
         required>
  
  <button type="submit">Sign In</button>
</form>
```

### 6. Registration Forms
```html
<form action="https://domain.com/register" method="post" autocomplete="on">
  <label for="name">Name</label>
  <input type="text" 
         id="name" 
         name="name" 
         autocomplete="name" 
         required>
  
  <label for="email">Email</label>
  <input type="email" 
         id="email" 
         name="email" 
         autocomplete="username email" 
         required>
  
  <label for="password">Password</label>
  <input type="password" 
         id="password" 
         name="password" 
         autocomplete="new-password" 
         required>
  
  <label for="confirm">Confirm Password</label>
  <input type="password" 
         id="confirm" 
         name="confirm" 
         autocomplete="new-password" 
         required>
  
  <button type="submit">Create Account</button>
</form>
```

### 7. Testing Checklist
- [ ] Test in Safari (Apple Keychain)
- [ ] Test in Chrome (Chrome Password Manager)
- [ ] Test in Firefox (Firefox Password Manager)
- [ ] Test in Edge (Edge Password Manager)
- [ ] Test with 1Password extension
- [ ] Test with LastPass extension
- [ ] Test with Bitwarden extension
- [ ] Test on iOS Safari
- [ ] Test on Android Chrome

### 8. Common Issues & Solutions
- **Apple Keychain not working**: Add `name="login"` hidden field
- **Chrome not suggesting**: Use `autocomplete="username"` for email
- **Firefox not recognizing**: Ensure proper form structure
- **Extensions not working**: Use standard autocomplete attributes
- **Mobile not working**: Test with mobile browsers

### 9. Priority Order
1. `autocomplete="username"` (primary for all platforms)
2. `autocomplete="email"` (secondary for Apple/Chrome)
3. `type="email"` and `type="password"`
4. Proper `name` attributes
5. Semantic form structure
6. Domain-specific `action` URLs 